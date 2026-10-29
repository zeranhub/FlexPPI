"""All-chain overlap, dependence components, and fixed external-family purge."""
from pathlib import Path
import json
import hashlib
from collections import defaultdict,Counter
import pandas as pd
from prepare_corrected_labels import sha
from run_corrected_static import corrected_records
HERE=Path(__file__).resolve().parent
OUT=HERE/"sequence_audit"

class Components:
    def __init__(self,items): self.parent={x:x for x in items}
    def find(self,x):
        while self.parent[x]!=x:
            self.parent[x]=self.parent[self.parent[x]]
            x=self.parent[x]
        return x
    def join(self,a,b): self.parent[self.find(a)]=self.find(b)

def main():
    chains=json.loads((OUT/"chain_inventory.json").read_text(encoding="utf-8"))
    items={r["sample_id"] for r in chains}
    graph=Components(items)
    sequence_records=defaultdict(set)
    sequence_neighbors=defaultdict(set)
    biological_records=defaultdict(set)
    record_chains=defaultdict(list)
    for r in chains:
        sequence_records[r["sequence_id"]].add(r["sample_id"])
        record_chains[r["sample_id"]].append(r)
        biological_records[(r["dataset"],r["biological_system"])].add(r["sample_id"])
    def join_group(group):
        ordered=sorted(group)
        for item in ordered[1:]:graph.join(ordered[0],item)
    for group in sequence_records.values():join_group(group)
    for group in biological_records.values():join_group(group)
    accepted=0
    pairs=set()
    # MMseqs prints fractional identity/coverage to three decimals. Enforce the
    # underlying integer aligned span for the coverage threshold as well.
    with (OUT/"all_chain_similarity.tsv").open(encoding="utf-8") as f:
        for line in f:
            fields=line.rstrip().split("\t")
            assert len(fields)==14
            a,b=fields[:2]
            identity=float(fields[2])
            qspan=(int(fields[5])-int(fields[4])+1)/int(fields[6])
            tspan=(int(fields[8])-int(fields[7])+1)/int(fields[9])
            if identity<.3 or min(qspan,tspan)<.7:continue
            sequence_neighbors[a].add(b)
            sequence_neighbors[b].add(a)
            if a==b:continue
            pair=tuple(sorted([a,b]))
            if pair in pairs:continue
            pairs.add(pair);accepted+=1
            join_group(sequence_records[a] | sequence_records[b])
    for key in sequence_records:sequence_neighbors[key].add(key)
    components=defaultdict(set)
    for sample in items:components[graph.find(sample)].add(sample)
    component_ids={root:"component_"+hashlib.sha256("\n".join(sorted(members)).encode()).hexdigest()[:16]
                   for root,members in components.items()}
    external={r["sample_id"] for r in chains if r["role"]=="external"}
    training=items-external
    external_components={graph.find(x) for x in external}
    removed={x for x in training if graph.find(x) in external_components}
    retained=training-removed
    audit_rows=[]
    for sample,observations in sorted(record_chains.items()):
        first=observations[0]
        if first["role"]=="training":
            allowed={r["sample_id"] for r in chains if r["role"]=="training" and r["fold"]!=first["fold"]}
        else:allowed=training
        seen_long=0
        long=0
        exact_short=0
        chain_rows=[]
        for r in observations:
            seq=r["sequence_id"]
            hits=set().union(*(sequence_records[k] for k in sequence_neighbors[seq])) & allowed
            if r["length"]>=20:
                long+=1;seen_long+=bool(hits)
            else:
                hits=sequence_records[seq]&allowed
                exact_short+=bool(hits)
            chain_rows.append(dict(**r,training_overlap_records=len(hits),exact_training_records=len(sequence_records[seq]&allowed)))
        novelty="no_long_chains" if not long else "all_long_chains_unseen" if not seen_long else "all_long_chains_seen" if seen_long==long else "partly_seen_long_chains"
        audit_rows.append(dict(sample_id=sample,role=first["role"],dataset=first["dataset"],fold=first["fold"],
            long_chains=long,seen_long_chains=seen_long,exact_short_seen_chains=exact_short,novelty=novelty,
            component_id=component_ids[graph.find(sample)],component_size=len(components[graph.find(sample)]),
            purged_training_member=sample in removed))
        # Every per-chain decision is preserved alongside the aggregate.
        for r in chain_rows: r["component_id"]=component_ids[graph.find(sample)]
        if "all_chain_rows" not in locals():all_chain_rows=[]
        all_chain_rows.extend(chain_rows)
    df=pd.DataFrame(audit_rows)
    df.to_csv(OUT/"record_overlap_audit.csv",index=False)
    pd.DataFrame(all_chain_rows).to_csv(OUT/"chain_overlap_audit.csv",index=False)
    df.groupby(["role","dataset","novelty"]).size().rename("records").to_csv(OUT/"novelty_counts.csv")
    original=corrected_records()
    rows=[r for r in original if r["sample_id"] in retained]
    multi=sum(len(record_chains[r["sample_id"]])>1 for r in rows)
    report={"training_records":len(training),"external_records":len(external),"accepted_nonself_unique_sequence_edges":accepted,
        "all_unique_long_queries":3770,"sequence_search_sha256":sha(OUT/"all_chain_similarity.tsv"),
        "search_input_sha256":sha(OUT/"long_unique.fasta"),"mmseqs_version":"eec9c354be4276d2373996af2e50808b1390d527",
        "component_rule":"Any qualifying long-chain sequence edge, any exact sequence including short peptides, or same biological_system within source.",
        "components":len(components),"largest_component_records":max(map(len,components.values())),
        "external_components":len(external_components),"purged_training_records":len(removed),"retained_training_records":len(rows),
        "retained_by_source":dict(Counter(r["dataset"] for r in rows)),"retained_multichain_records":multi,
        "strict_fit_feasible":len(rows)>=500 and multi>=20,"prespecified_plan_sha256":sha(HERE/"strict_external_prespecified_plan.json"),
        "limitations":"All chains searched at sensitivity7.5, not a mathematical guarantee of detecting every evolutionary homolog. Components may be conservatively joined by shared partners. External source-out alone is not unseen-family evidence."}
    (OUT/"purged_training_local.jsonl").write_text("\n".join(json.dumps(r,ensure_ascii=False) for r in rows)+"\n",encoding="utf-8")
    (OUT/"sequence_audit_summary.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
