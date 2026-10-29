# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Published affinity dots -> BarComparison/log-scale behavior -> param inherit
# (b) Per-residue curves -> LineTrend/plot_trend.py -> param inherit
# (c,d) Signed residue maps -> heatmap/plot_comparison.py -> param inherit
# No asset has the four-state, sequence-aligned multi-region data structure.
# Companion PNGs inspected; preserve clean spines, aligned cells, zero-centred colors.
# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})
# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"
# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

from pathlib import Path
from datetime import datetime, timezone
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle, Patch

ROOT = Path(__file__).resolve().parent
IN = ROOT/'inventory'
OUT = ROOT/'case_figure'
OUT.mkdir(exist_ok=True)
STATES = ['3hg1', '4jff', '4jfd', '4jfe']
NAMES = ['MEL5 · ELA', 'α24β17 · ELA', 'α24β17 · G4A', 'α24β17 · I7A']
COLORS = [CATEGORICAL[3], CATEGORICAL[4], CATEGORICAL[0], CATEGORICAL[2]]
STYLES = ['-', '-', '--', ':']
PAIRS = [('3hg1', '4jff'), ('4jff', '4jfd'), ('4jff', '4jfe')]
PAIR_NAMES = ['TCR maturation', 'Peptide G4A', 'Peptide I7A']
SHA = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()

def label(ax, letter, title):
    ax.set_title(title, loc='left', pad=10)
    ax.text(-.055, 1.13, letter, transform=ax.transAxes, fontsize=11, fontweight='bold', ha='right')

def main():
    source_names = ['MEL5_four_states_complete_peptide.csv', 'MEL5_three_contrasts_all_six_CDR.csv', 'MEL5_three_contrasts_four_state_union_pMHC_interface.csv']
    pep, cdr, iface = [pd.read_csv(IN/name) for name in source_names]
    assert len(pep) == 40 and len(cdr) == 135 and len(iface) == 90
    fig = plt.figure(figsize=(183/25.4, 185/25.4))
    # One dominant full-peptide panel, two complete regional maps below.
    a = fig.add_axes([.19, .685, .255, .235])
    b = fig.add_axes([.56, .685, .415, .235])
    c = fig.add_axes([.19, .405, .785, .13])
    d = fig.add_axes([.19, .145, .785, .13])
    kd = np.array([18e-6, .6e-9, 36e-6, 31e-6])
    for i, (value, color) in enumerate(zip(kd, COLORS)):
        a.scatter(value, i, s=24, color=color, zorder=3)
        a.text(value*1.6, i, ['18 μM', '0.6 nM', '36 μM', '31 μM'][i], va='center', fontsize=7)
    a.set_xscale('log')
    a.set_xlim(1e-10, 1e-3)
    a.set_ylim(3.6, -.6)
    a.set_yticks(range(4), NAMES, fontsize=7)
    a.set_xticks([1e-9, 1e-7, 1e-5, 1e-3])
    a.set_xlabel('Published SPR $K_D$ (M)')
    a.tick_params(axis='y', length=0)
    label(a, 'a', 'Binding affinity at 25 °C')
    a.text(0, -.39, 'Lower $K_D$ = stronger binding', transform=a.transAxes, fontsize=6.5)
    for state, name, color, style in zip(STATES, NAMES, COLORS, STYLES):
        q = pep[pep.pdb.eq(state)].sort_values('id')
        assert len(q) == 10
        b.plot(q.id, q.rmsf_full, color=color, linestyle=style, marker='o', markersize=2.7, linewidth=1.1, label=name)
    b.set_xlim(.7, 10.3)
    b.set_xticks(range(1, 11))
    b.set_ylim(.57, 1.02)
    b.set_ylabel('Predicted RMSF (Å)')
    b.set_xlabel('Peptide position')
    label(b, 'b', 'Complete peptide dynamics')
    b.axvline(4, color=GREY, lw=.6, alpha=.6, zorder=0)
    b.axvline(10, color=GREY, lw=.6, alpha=.6, zorder=0)
    b.legend(loc='upper center', bbox_to_anchor=(.5, -.32), ncol=2, fontsize=6.5, columnspacing=.9, handlelength=2)
    b.text(.02, .98, 'P4: recognition region\nP10: C-terminal HLA anchor', va='top', transform=b.transAxes, fontsize=6.3)
    source_rows = [dict(panel='a', state=s, state_name=n, KD_M=float(v), source='Madura et al. 2013 JBC Table2/Fig1, doi10.1074/jbc.M113.464560') for s,n,v in zip(STATES,NAMES,kd)]
    for r in pep.to_dict('records'):
        source_rows.append(dict(panel='b', state=r['pdb'], peptide_position=r['id'], RMSF_A=r['rmsf_full']))
    cmap = LinearSegmentedColormap.from_list('signed_delta', DIVERGING).copy()
    cmap.set_bad('#D9D9D9')
    norm = TwoSlopeNorm(vmin=-.45, vcenter=0, vmax=.45)
    ref = cdr[(cdr.pdb_before.eq('3hg1')) & (cdr.pdb_after.eq('4jff'))].copy().reset_index(drop=True)
    ref_keys = [(r.chain_before, str(r.id_before)) for r in ref.itertuples()]
    engineered_to_ref = {(r.chain_after,str(r.id_after)):key for key,r in zip(ref_keys,ref.itertuples())}
    # IDs may be parsed as int or float by pandas; normalize explicitly.
    def rid(x): return str(int(x))
    ref_keys = [(r.chain_before, rid(r.id_before)) for r in ref.itertuples()]
    engineered_to_ref = {(r.chain_after,rid(r.id_after)):key for key,r in zip(ref_keys,ref.itertuples())}
    cmat = []
    for pair_index, (before, after) in enumerate(PAIRS):
        q = cdr[cdr.pdb_before.eq(before) & cdr.pdb_after.eq(after)]
        values = {}
        for r in q.itertuples():
            key = (r.chain_before, rid(r.id_before))
            if before != '3hg1': key = engineered_to_ref[key]
            values[key] = r.rmsf_change_after_minus_before
            source_rows.append(dict(panel='c', contrast=PAIR_NAMES[pair_index], reference_chain=key[0], reference_id=key[1], before_id=r.id_before, after_id=r.id_after, delta_RMSF_A=r.rmsf_change_after_minus_before))
        assert set(values) == set(ref_keys)
        cmat.append([values[key] for key in ref_keys])
    image_c = c.imshow(cmat, cmap=cmap, norm=norm, aspect='auto', interpolation='nearest')
    c.set_yticks(range(3), PAIR_NAMES, fontsize=7)
    c.set_xticks(range(45), [key[1] for key in ref_keys], rotation=90, fontsize=5.5)
    c.tick_params(length=0)
    label(c, 'c', 'All six recognition loops · predicted ΔRMSF')
    starts = [0, 6, 12, 21, 27, 32]
    ends = [6, 12, 21, 27, 32, 45]
    for start, end, name in zip(starts, ends, ['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']):
        c.text(((start+end)/2)/45, 1.055, name, transform=c.transAxes, ha='center', fontsize=6.5)
        if start: c.axvline(start-.5, color='white', lw=1.5)
    for sp in c.spines.values(): sp.set_visible(False)
    c.text(0, -.42, 'MEL5 reference residue numbering; all 45 aligned positions retained', transform=c.transAxes, fontsize=6.5)
    ic = iface[(iface.pdb_before.eq('3hg1')) & (iface.pdb_after.eq('4jff'))].reset_index(drop=True)
    ikeys = [(r.chain_before, rid(r.id_before)) for r in ic.itertuples()]
    imat = []
    for pair_index,(before,after) in enumerate(PAIRS):
        q = iface[iface.pdb_before.eq(before) & iface.pdb_after.eq(after)]
        values = {}
        for r in q.itertuples():
            key = (r.chain_before, rid(r.id_before))
            eligible = bool(r.hchi_valid and r.same_aa)
            value = float(r.hchi_change_after_minus_before) if eligible else np.nan
            values[key] = value
            source_rows.append(dict(panel='d', contrast=PAIR_NAMES[pair_index], chain=key[0], author_id=key[1], same_AA_valid_Hchi=eligible, delta_Hchi_nats=value))
        assert set(values) == set(ikeys)
        imat.append([values[key] for key in ikeys])
    assert np.isfinite(imat).sum(axis=1).tolist() == [23,23,22]
    d.imshow(imat, cmap=cmap, norm=norm, aspect='auto', interpolation='nearest')
    d.set_yticks(range(3), PAIR_NAMES, fontsize=7)
    d.set_xticks(range(30), [f'{r.aa_before}{int(r.id_before)}' for r in ic.itertuples()], rotation=90, fontsize=5.8)
    d.tick_params(length=0)
    label(d, 'd', 'Entire pMHC interface · predicted ΔHχ')
    d.axvline(21.5, color='white', lw=2)
    d.text(11/30, 1.055, 'HLA', transform=d.transAxes, ha='center', fontsize=6.5)
    d.text(26/30, 1.055, 'Peptide', transform=d.transAxes, ha='center', fontsize=6.5)
    for sp in d.spines.values(): sp.set_visible(False)
    # Q72 is a pre-existing fixed anchor; preserve all states including its adverse pattern.
    q72 = ikeys.index(('A','72'))
    d.add_patch(Rectangle((q72-.5,-.5),1,3,fill=False,edgecolor=BLACK,lw=.65))
    cax1 = fig.add_axes([.79,.335,.185,.012])
    cb1 = fig.colorbar(image_c,cax=cax1,orientation='horizontal',ticks=[-.45,0,.45])
    cb1.ax.tick_params(labelsize=6,length=2)
    cb1.ax.set_xlabel('ΔRMSF (Å)', fontsize=6.5, labelpad=1)
    cax2 = fig.add_axes([.79,.065,.185,.012])
    cb2 = fig.colorbar(image_c,cax=cax2,orientation='horizontal',ticks=[-.45,0,.45])
    cb2.ax.tick_params(labelsize=6,length=2)
    cb2.ax.set_xlabel('ΔHχ (nats)', fontsize=6.5, labelpad=1)
    d.text(0, -.48, 'Four-state interface union: 30 positions; same-AA Hχ valid at 23 / 23 / 22.', transform=d.transAxes, fontsize=6.2)
    fig.text(.19,.017,'Δ = after − reference: blue, decrease; red, increase. Grey: Hχ not comparable.\nPublished affinity and model dynamics are separate measurements; no affinity classifier is fitted.',fontsize=6.3)
    save_cns_figure(fig, OUT/'MEL5_opposite_affinity_case_series')
    fig.savefig(OUT/'MEL5_opposite_affinity_case_series.svg',bbox_inches='tight')
    plt.close(fig)
    pd.DataFrame(source_rows).to_csv(OUT/'figure_exact_source.csv',index=False,encoding='utf-8-sig')
    caption = ('MEL5 lineage four-state analysis. (a) Published SPR dissociation constants at 25 C: MEL5-ELA 18 uM; alpha24beta17-ELA 0.6 nM; peptide G4A 36 uM; I7A 31 uM (Madura et al., JBC 2013). These are reported point estimates; no unreported errors were added. (b) All ten model-predicted peptide RMSF values in each observed complex. (c) Three fixed contrasts across all 45 aligned CDR residues. The six CDR definitions are projected from MEL5 by sequence alignment. (d) Same-amino-acid chemically defined Hchi changes across the union of all four pMHC interfaces, with nearest TCR C-alpha distance <=8 A. All 30 positions are retained; grey cells denote ineligible comparisons. Q72 is outlined as a previously named anchor. All heatmaps use after-minus-reference differences. Full observed-structure predictions use the same 2430-record fold0 joint model, without ITrace training. The panels show spatial dynamic predictions and published binding outcomes, not a validated affinity-direction classifier. Existing MD has mixed support, including an incorrect G4A P10 RMSF decrease and failure of the Q72 anchor in weak states; see the complete report and unfiltered source tables.')
    (OUT/'caption_en.txt').write_text(caption,encoding='utf-8')
    manifest = dict(status='generated_pending_visual_QA', created_utc=datetime.now(timezone.utc).isoformat(),core_claim='Complete spatial predictions distinguish regional responses hidden by a maximum statistic; this does not establish opposite-affinity prediction.',archetype='asymmetric_mixed',source_files=[dict(path=str(IN/n),sha256=SHA(IN/n)) for n in source_names], n=dict(observed_states=4,fixed_contrasts=3,peptide_values=40,CDR_comparisons=135,pMHC_interface_comparisons=90,Hchi_comparable_by_pair=[23,23,22]),all_data_preserved=True,primary_model_checkpoint_sha256='ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5',no_new_inference=True,manuscript_edited=False,statistics='Published affinity point estimates, one deterministic prediction per state; existing MD summary mean of 3 trajectories is archived separately. Residues are not independent biological replicates. No significance tests or fitted affinity model.',MD_temperature_K=310,SPR_temperature_C=25,caption=caption,asset_route='Parameter inheritance after inspecting production scripts and companion PNGs; incompatible data and scientific axis meanings prevent a native unmodified asset execution.',files=[dict(path=str(p),sha256=SHA(p)) for p in sorted(OUT.glob('*')) if p.is_file()])
    (OUT/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Figure generated with all40 peptide values,135 CDR and90 pMHC comparisons.')

if __name__ == '__main__': main()
