from pathlib import Path
from hashlib import sha256
import shutil,json
import pandas as pd

root=Path(__file__).resolve().parent
source=root.parent/'revision_v110_compact_five_panels'
previous=root.parent/'paperbanana_from_scratch_20261027'
ref=root/'data/PaperBananaBench/diagram';ref.mkdir(parents=True,exist_ok=True)
shutil.copy2(previous/'data/PaperBananaBench/diagram/ref.json',ref/'ref.json')
style=root/'style_guides';style.mkdir(exist_ok=True)
fig2=root.parents[2]/'figure2_dense_20261024/fig2_dense.png'
assert fig2.exists()
shutil.copy2(fig2,root/'manuscript_figure2_style_reference.png')
inventory=pd.read_csv(source/'figure3/main_original_panels_a_to_e_values.csv',float_precision='round_trip')
assert len(inventory)==767 and set(inventory.panel)=={'b','c','e'}
inventory.to_csv(root/'frozen_numeric_inventory_all767.csv',index=False)
profiles=pd.read_csv(source/'case_profiles_all97_with_masks.csv',float_precision='round_trip')
sasa=pd.read_csv(source/'case_SASA_all97_starting_complex.csv',float_precision='round_trip')
stats=json.loads((source/'direct_relationship_statistics.json').read_text(encoding='utf-8'))

method='''TASK: FROM-SCRATCH PAGE LAYOUT for a quantitative scientific RESULTS figure in a protein flexibility manuscript. Generate a genuinely NEW coherent composition, rather than polishing or reproducing an old figure. No previous Figure3 canvas, panel coordinates or complete Figure3 screenshot is supplied. The supplied Figure2 is a STYLE REFERENCE ONLY: follow its fonts, colors, white space discipline and quantitative plotting language, not its scientific contents. This is NOT a machine-learning architecture, flowchart or conceptual mechanism diagram.

SCIENTIFIC READING PATH
Two complementary structure-resolved examples show whether the model predicts backbone RMSF and joint side-chain entropy Hχ correctly, with starting-complex solvent exposure as environmental context. Then show the limited global RMSF–Hχ correspondence and the side-chain prediction's relation to explicit angular occupancy. Do not invent activation mechanisms, causal solvent-exposure effects or experiments.

COMPOSITIONAL GOAL
Replan ALL positions, widths, heights, ordering, shared legends and white space. Make the two cases easy to compare: their sequence traces and SASA subtracks should have equal usable widths/heights and matching axes, avoiding a zigzag reading order. The two real molecular contexts must remain adjacent in the same structural region and be large enough to read; place them to associate naturally with their respective profiles. Density information is a compact supporting plot, not a wide bottom banner. Four χ scatterplots must stay IN A SINGLE HORIZONTAL ROW spanning the figure width; never 2x2. Use a compact data-led hierarchy, balanced fill, coherent alignments and sensible type sizes. Do not force five equal boxes; you may form coherent subregions but retain five logical groups. No mechanism cartoons, tables, local structural zooms, decorations, tinted cards, thick frames, arrows connecting panels, overall figure title, caption prose inside the image, gray explanatory microtext or newly added evidence. Make an intelligent new page, not a wireframe. Do not reintroduce the deleted six-Gln scatter/SASA panel, Intra curves or 83-record comparisons.

FIVE LOGICAL EVIDENCE GROUPS (labels may be reordered, these are NOT fixed positions):

STRUCTURES: TWO adjacent COMPLETE genuine ribbon contexts: CaM–CaMKI 1MXE and ubiquitin–UBA 2QHO. Analyzed chains purple #8B72A5, partners light gray #BAC1C7, key positions teal #3A8D8C. CaMKI W303/M316 and ubiquitin Q41/Q62 are labeled. No local zooms, no isolated single residue sticks. At concept stage molecular depictions are composition placeholders; native implementation will retain the true PyMOL renders. Do not invent molecular chemistry as scientific evidence.

CASE1: CaMKI regulatory peptide native residues294–318, 25 total positions. One dual-axis line plot: blue #2166AC RMSF left axis0–2.4Å, ochre #B87936 Hχ right axis0–4.5nats. COLOR means METRIC, LINE STYLE means SOURCE: Full solid, MD dashed. Both MD and Full shown for both metrics. Hχ has21 valid positions, Ala299/306/309/311 remain gaps, not zero, not interpolation; tiny markers preserve isolated valid positions. No Intra. Directly underneath this trace put a compact, same-width, position-aligned TEAL/GRAY #3A8D8C side-chain SASA bar track, axis0–160Å², all25 positions, only bottom track repeats residue numbers. SASA is from exact STARTING COMPLEX, computed PyMOL3.0/probe1.4Å/dot_density4, NOT model output or MD average. Do not show obsolete CaM shielding%. RMSF Full MAE0.115Å; Hχ Full MAE0.254nats (21 valid). W303/M316 positional indications allowed; statistics can be in caption rather than cluttering plot.

CASE2: Ubiquitin chainA residues1–72. SAME dual-axis blue RMSF/ochre Hχ, SAME0–2.4Å/0–4.5nats and SAME Full solid/MD dashed styles as case1. RMSF72 valid, Hχ66 valid. G10/A28/G35/A46/G47/G53 are Hχ gaps. Valid Pro19/37/38 MD Hχ zero values MUST remain, with Full nearzero values. Directly underneath put same-height, same-width side-chain SASA0–160Å² bar track aligned to all72 positions, same color as CASE1 bars. Gly side-chain SASA0 is due to no sidechain heavy atoms, not burial. A28/A46 still have real SASA bars despite undefined Hχ. Q41/Q62 labels/reference positions allowed. RMSF Full MAE0.040Å; Hχ Full MAE0.080nats.

GLOBAL RELATIONSHIP: ONE compact overlaid joint density plot with small upper/right marginal histograms, common x=RMSF0–10Å and y=Hχ0–4.5nats. Both MD and Full use SAME824,988 jointly valid residue positions across3,063 existing heldout records. MD gray #62676B and Full blue #2166AC, SOFT TRANSPARENT LAYERS: nonzero alpha0.10–0.35, overlap subtly darker blue-gray, zero cells white. ALL nonzero gray and blue hexbin cells visible, including sparse MD right tail. NO contour lines, regression, confidence ellipses, thresholds or 90% envelopes. MD Spearmanρ0.209, Fullρ0.262, computed on all824988; main view extends0–10Å (MD4234 and Full906 higher-RMSF values are in full-range companion). Upper/right histograms gray dashed MD, blue solid Full. This is global distribution shape, NOT residual scatter, not predicted-vs-reference parity, not difference plot. It must be modest width, roughly one case panel wide or smaller, not a full-width band.

ANGLE OCCUPANCY: Four scatter plots chi1/chi2/chi3/chi4 IN ONE ROW, equal widths/heights/common x/y scales. X=Predicted JOINT Hχ deviation after amino-acid-type centering (nats), Y=MD SINGLE-ANGLE fine-bin occupancy entropy deviation after same centering (nats). Xlim[-.7,.5], ylim[-1,.75]. χ1 n66ρ.50, χ2 n52ρ.45, χ3 n23ρ.52, χ4 n10ρ.13. Preserve weak chi4, no forced trendline. Shared x/y label once. Other points blue circles, Q41 teal diamond and Q2 teal square. Q41/Q2 here differs from the structure's Q41/Q62 environment labels. All151 single-angle observations retained. No four separate prediction heads; one joint Hχ prediction reused for defined angles. Subtle dotted zero guides acceptable.

MANUSCRIPT FIGURE2 STYLE TO MATCH
183mm double-column width, compact HEIGHT appropriate to evidence (try about160–190mm without shrinking text). White background; Arial9pt axes/ticks/legend,10pt short titles,11pt bold lowercase panel letters; final manuscript scaledwidth167.6mm must still read. Thin0.6pt axes, outticks2.5pt, no upper spine and only essential right spine for dual-axis cases. No pale gray captions or technical production notes in image. Molecular colors purple/gray/teal, measured metric colors blue RMSF/ochre Hχ, Full solid/MD dashed. The density separately compares sources with gray MD/blue Full and must have its own small local legend. Shared legends near relevant plots. Do not erase either semantic distinction. Reference Figure2 supplies this exact plotting language.

NUMERICAL REFERENCE
Layout is the generative tool's task; true data-native curves and structures will be reimplemented afterwards. You may render a design concept rather than force exact individual pixel placement, but definitions, axes, sample sizes, weak findings and evidence inventory must remain correct. Do not print any 'concept', 'placeholder', model notes, numbers table or production remarks in the figure. No old layout is supplied. The following full current sequence values and scatter coordinates define the scientific content; they do NOT dictate positions.
'''
def records_no_nan(df):
    return json.loads(df.to_json(orient='records',double_precision=15))
inventory_json={'profiles_all97':records_no_nan(profiles),'static_sidechain_SASA_all97':records_no_nan(sasa),
 'scatter_coordinate_records':records_no_nan(inventory[inventory.panel.eq('e')]),
 'global_relationship':{'n':824988,'records':3063,'MD_rho':stats['panels']['MD']['rho'],'Full_rho':stats['panels']['Full']['rho']}}
method+='\n'+json.dumps(inventory_json,ensure_ascii=False,separators=(',',':'))
(root/'scientific_content.txt').write_text(method,encoding='utf-8')
intent='''From-scratch layout proposal for a compact scientific results figure. Retain five current evidence groups: adjacent CaM–CaMKI/ubiquitin–UBA ribbon structures; TWO matched dual-axis RMSF/Hχ sequence plots with aligned starting-complex side-chain SASA subtracks; ONE small overlaid MD/Full joint density with marginal histograms and NO contours; FOUR χ occupancy scatters across one horizontal row. Invent page-wide spatial logic and clean visual hierarchy, improve case comparability, reduce zigzag reading and whitespace. Coordinate exactly with supplied manuscript Figure2 style: Arial9/10/11,183mm width, white canvas, blue/ochre metric colors, solid/dashed source styles, purple-gray molecular components. No six-Gln panel, tables, local zooms, Intra curves, mechanism diagrams, decorations or gray microtext. Do not reproduce any previous Figure3. Show a complete new design rather than a wireframe. Numeric source is for preservation, not fixed layout.'''
(root/'visual_intent.txt').write_text(intent,encoding='utf-8')
(style/'neurips2025_diagram_style_guide.md').write_text('''# FlexPPI manuscript Figure2 style binding
This filename is the author pipeline convention. The target is quantitative protein-science RESULTS, not AI architecture. The supplied Figure2 image is style only. White, Arial9pt axes/ticks/legend,10pt short titles,11pt bold panel letters,183mm width, about160–190mm high. Fine0.6pt axes, outticks2.5pt, zero extra decorations. Analyzed chain purple8B72A5/partnergrayBAC1C7/keyteal3A8D8C. In both case traces blue2166AC=RMSF left axis0–2.4Å; ochreB87936=Hχ right axis0–4.5nats; Fullsolid/MDdashed. Same-height teal sidechainSASA bars0–160 below each trace. In local GLOBAL DENSITY legend gray62676B=MD/blue2166AC=Full, transparent alpha.10–.35, no contours. All sequence gaps and legitimate zeros remain. Four single-angle scatters must be one full-width horizontal row. TWO structures stay adjacent and substantial; case profiles equal size; density modest width; thoughtful page-wide reading flow. No overall title, tables, paragraphs, local structural zooms, mechanism/cartoon, colored cards, shadows or gray microtext. Do not resurrect deleted sixGln/Intra/83record panels. PaperBanana proposes layout; final source-native drawing will use exact real input data and molecular renders.
''',encoding='utf-8')

runner=(previous/'run_from_scratch.py').read_text(encoding='utf-8')
old_intent='Independent composition A: invent your strongest compact unified evidence composition, with asymmetric data-first hierarchy and no fixed inherited grid.'
new_intent='Independent composition A: invent an elegant new data-first composition organized for side-by-side comparison of the two complete cases. The two structures stay adjacent, the two trace+SASA blocks match sizes, the four chi plots stay one full row, and the density stays compact. No fixed inherited grid.'
assert old_intent in runner;runner=runner.replace(old_intent,new_intent)
old_intent='Independent composition B: invent a distinct spatial organization with an especially strong coherent side-chain evidence region and a complementary backbone case. Keep both structures adjacent; do not let any single sparse plot dominate the page.'
new_intent='Independent composition B: explore a genuinely different page-wide hierarchy from A, optimizing a coherent integrated case-evidence region, adjacent substantial structures, matched trace+SASA blocks, a compact global density, and one full row of four chi plots. Eliminate zigzag reading and inherited panel positions.'
assert old_intent in runner;runner=runner.replace(old_intent,new_intent)
# Supply the actual manuscript Figure2 as a style-only few-shot input to the
# original Planner, alongside the author's automatically selected examples.
# No author module or system prompt is modified; no old Figure3 image is supplied.
needle='    install_visualizer_checkpoints(processor.visualizer_agent, args.output)'
injection='''    original_planner_process = processor.planner_agent.process
    async def planner_with_manuscript_style(data):
        refs = json.loads((Path(__file__).resolve().parent / "data/PaperBananaBench/diagram/ref.json").read_text(encoding="utf-8"))
        by_id = {item["id"]: item for item in refs}
        examples = [by_id[item] for item in data.get("top10_references", []) if item in by_id]
        style_example = {"id":"manuscript_Figure2_style_only", "content":"STYLE ONLY: Existing Figure2 from the same manuscript. Match its compact quantitative visual language, Arial typography, blue/ochre data colors, delicate axes and panel labels; its scientific content and layout are not the target.", "visual_intent":"Use this visual grammar for a NEW protein case-results composition. Never reproduce this figure's plots or duplicate its layout.", "path_to_gt_image":str(Path(__file__).resolve().parent / "manuscript_figure2_style_reference.png")}
        data["retrieved_examples"] = [style_example] + examples
        return await original_planner_process(data)
    processor.planner_agent.process = planner_with_manuscript_style
'''+needle
assert needle in runner;runner=runner.replace(needle,injection)
runner=runner.replace('old_figure_supplied=False, local_style_override=True, source_commit=TOOL_COMMIT, entrypoint=',
 'old_figure_supplied=False, local_style_override=True, manuscript_figure2_style_image=True, author_agents_unmodified=True, source_commit=TOOL_COMMIT, entrypoint=')
(root/'run_from_scratch.py').write_text(runner,encoding='utf-8')
files=[root/'scientific_content.txt',root/'visual_intent.txt',root/'run_from_scratch.py',root/'manuscript_figure2_style_reference.png',style/'neurips2025_diagram_style_guide.md',ref/'ref.json',root/'frozen_numeric_inventory_all767.csv']
manifest={'workflow':'PaperBanana author demo_full from scratch with auto retrieval, Planner style-only Figure2 injection,2 independent layout intents',
 'old_Figure3_supplied':False,'Figure2_style_supplied':True,'native_input_values':767,'source':str(source),
 'files':[{'path':str(p),'sha256':sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in files]}
(root/'input_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'prepared':str(root),'content_characters':len(method),'numeric_values':767,'old_figure3_supplied':False,'Figure2_style_image':True},ensure_ascii=False))
