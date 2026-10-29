from pathlib import Path
import json,subprocess,urllib.request,hashlib,datetime,re

ROOT=Path(__file__).parent
RAW=ROOT/'raw';RAW.mkdir(exist_ok=True)
REPO='StatPhysBio/tcr_antigen_design'
PINNED_COMMIT='cfd85f3fefa652f0f00d3fd1e58cb2b4256f5ce2'
def gh(path):
    p=subprocess.run(['gh','api',path],capture_output=True,check=True)
    return p.stdout
def download(url,name):
    p=RAW/name;p.parent.mkdir(parents=True,exist_ok=True)
    req=urllib.request.Request(url,headers={'User-Agent':'FlexPPI-source-audit/1.0'})
    data=urllib.request.urlopen(req,timeout=50).read();p.write_bytes(data)
    return p

if __name__=='__main__':
    commit=json.loads(gh(f'repos/{REPO}/commits/{PINNED_COMMIT}'))
    (RAW/'github_commit.json').write_text(json.dumps(commit,indent=2),encoding='utf8')
    sha=commit['sha']
    tree=json.loads(gh(f'repos/{REPO}/git/trees/{sha}?recursive=1'))
    (RAW/'github_tree.json').write_text(json.dumps(tree,indent=2),encoding='utf8')
    candidates=[]
    for x in tree['tree']:
        s=x['path']
        if x['type']=='blob' and s.startswith('mutation_effects/') and '/results/' not in s and '/plots/' not in s and any(q in s.lower() for q in ['1g4','nyeso','ny_eso','tax','a6','chen','ding','cole']):
            candidates.append({'path':s,'size':x.get('size')})
    (ROOT/'source_candidates.json').write_text(json.dumps({'repository':REPO,'commit':sha,'files':candidates},indent=2),encoding='utf8')
    download(f'https://raw.githubusercontent.com/{REPO}/{sha}/README.md','github_README.md')
    download('https://pmc.ncbi.nlm.nih.gov/articles/PMC12557503/','pnas_PMC12557503.html')
    print(json.dumps({'commit':sha,'tree_files':len(tree['tree']),'candidate_source_files':len(candidates)},indent=2))
