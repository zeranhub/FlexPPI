from pathlib import Path
import csv
import hashlib
import json
import shutil

root = Path(__file__).resolve().parent
tools = Path(r'C:\Users\DELL\.cache\flexppi-figure-tools')
repo = tools / 'PaperBanana'
source = root.parent / 'revision_v98_structure_fill' / 'figure3' / 'source_data'
root.mkdir(parents=True, exist_ok=True)

method = '''Create an entirely NEW multi-panel scientific RESULTS figure from scratch for a protein flexibility prediction paper, suitable for Nature Communications or Genome Biology. No old complete figure, screenshot, canvas, panel coordinates, or previous layout is supplied. Plan the visual composition afresh from the following scientific evidence. This is a results figure with molecular structures and quantitative plots, NOT a neural-network architecture, flow chart, or infographic about the PaperBanana agents.

SCIENTIFIC ARGUMENT
Complex context improves backbone flexibility prediction in a regulatory peptide; a complementary ubiquitin example demonstrates different conformational diversity for the same side-chain type in different local environments. The predicted side-chain entropy is then related to explicit angular occupancy descriptors, first within the example and then across 83 records. Lead the reader from prediction quality to physical interpretation and broader occupancy relevance. No activation or binding-thermodynamic mechanism is proven.

DESIGN FREEDOM
Invent a compact, information-rich composition with a clear reading path, strong spatial relationships, generous usable plot areas, and balanced molecular and quantitative evidence. Let evidence density determine panel sizes. Do not make a rigid array of six equal boxes. You may combine related blocks, vary sizes and orientations, and choose panel lettering. Keep the two molecular structures together, complete, adjacent, and visually substantial. Avoid huge empty molecular margins and whitespace gaps. The complete 2x2 angular scatter family needs substantial readable area. Use visual proximity rather than decorative arrows. No tables, mechanistic cartoons, local structural zooms, decorative cards, shadows, colored boxes, or gray explanatory microtext. Concise titles, necessary axes, shared legends, residue labels, and selected statistics are sufficient. Do not put a large overall title, explanatory paragraphs, or this production note inside the image.

EVIDENCE INVENTORY (these are scientific roles, not fixed panel positions)
1. Molecular context: complete ribbon structures of CaM-CaMKI (PDB 1MXE, CaMKI W303 and M316) and ubiquitin-UBA (PDB 2QHO, ubiquitin Q41 and Q62). Both analyzed chains purple, partners light gray, target residues teal. Molecular shapes at this generative concept stage are composition references; final publication production will use genuine molecular renders. No local magnification insets or extra isolated sticks.

2. CaMKI backbone: absolute RMSF in Angstrom across all 25 residues 294-318, vertical baseline zero. Three distinguishable curves: MD reference, Full model, Intra model. Full follows the low central RMSF and more mobile ends. Full MAE=0.115 Angstrom; Intra MAE=2.383 Angstrom. W303: MD0.606, Full0.469, Intra3.023 Angstrom. M316: MD0.624, Full0.689, Intra3.883 Angstrom. Keep the 25-residue trace readable. Use W303/M316 to connect structure and profile. Full and Intra are two model-context settings on the SAME bound input, not bound/free experimental states. No CaMKI Hchi plot.

3. Ubiquitin side chains: absolute joint Hchi in nats over all 72 sequence positions, MD and Full curves, baseline zero. There are 66 valid sites. Missing positions remain gaps; valid proline zero values remain zero. Full MAE=0.080 nats. Highlight Q41 and Q62. Show residue-specific prediction quality across the complete sequence. No ubiquitin RMSF plot. Do not substitute delta values for the absolute profiles.

4. Same-type side-chain environments: all six ubiquitin Gln residues, side-chain SASA in square Angstrom versus MD and Full Hchi in nats, concise paired symbols. Data (residue/SASA/MD/Full): Q2/71.24/3.009/2.990; Q31/80.56/2.439/2.462; Q40/70.57/2.186/2.554; Q41/0.00/1.890/2.000; Q49/69.23/3.004/2.945; Q62/100.67/2.631/2.684. Emphasize Q41/Q62 without hiding other Gln. These data do NOT establish a monotonic relationship. No regression line or causal arrow. Neither Q41 nor Q62 directly contacts UBA within 5 Angstrom. The interpretation is local environment-dependent diversity, not direct partner contact at these sites.

5. Ubiquitin angular occupancy association: a complete legible 2x2 scatter family for chi1, chi2, chi3, chi4. X: predicted JOINT Hchi deviation from amino-acid-type mean, nats. Y: MD SINGLE-angle fine-bin occupancy entropy deviation from amino-acid-type mean, nats. Same joint prediction reused in four comparisons; no four angle-specific prediction heads. chi1 n=66 Spearman rho=0.50; chi2 n=52 rho=0.45; chi3 n=23 rho=0.52; chi4 n=10 rho=0.13. Both positive and negative deviations. Weak chi4 remains visibly weak. Emphasize Q41 as a teal diamond and Q2 as a teal square, other sites blue circles. Q41/Q2 are a distinct low/high prediction example; do not confuse them with Q41/Q62 environment comparison. Subtle zero guides allowed; no forced regression fit or claim of uniformly strong agreement.

6. Broader occupancy association across 83 RECORDS: show all record-level correlations and mean with bootstrap 95% interval for each of four groups. Full Hchi and Full RMSF are compared against two MD angular occupancy descriptors. Fine-bin entropy: Hchi mean rho0.510, CI0.491-0.527; RMSF mean0.321, CI0.305-0.336. Low-F* angular width: Hchi mean0.400, CI0.382-0.416; RMSF mean0.239, CI0.221-0.255. There are83 values in EACH group. The correlations control amino-acid type within each record. Keep visible record dispersion and uncertainty; a sparse bar chart of four means is insufficient. Hchi is more directly related to angular occupancy than backbone RMSF. F*=-ln[p/pmax] is a dimensionless relative occupancy potential. Low-F* width is NOT binding free energy. Do not relabel83 records as83 proteins, complexes, independent trajectories or wet experiments.

STYLE
White background, black readable text, restrained thin axes, no grids except subtle zero guides, clear font hierarchy, disciplined alignment. Coordinated palette: analyzed molecular chain #8B72A5; partner #BAC1C7; key residues #3A8D8C; Full prediction #2166AC; Intra prediction #B87936; MD #222222. Do not reuse model-method colors to encode different molecular partners. In the cross-record comparison use marker form (solid circle Hchi versus open triangle RMSF) plus labels to distinguish two Full predictors. Proper Greek chi and rho (rho must not turn into p). Target roughly180mm double-column width, about200mm height, typography readable at164mm width, minimum8pt equivalent. Shared legends concise and near relevant evidence. Prefer space for the data rather than tiny points inside large margins.

PRODUCTION STATUS
Generate a convincing complete DESIGN CONCEPT from scratch, not an empty wireframe and not a polish of a supplied whole figure. Exact numerical inventory follows for scientific reference; it does not dictate layout. This generative concept need not replicate every individual point exactly, but it must retain all evidence groups, definitions, sample sizes, statistics, weak findings and plot directions. Genuine structures and exact curves/scatter will later be implemented from frozen sources using native plotting. Do not print this production-status paragraph in the image. Never invent new experiments or a causal mechanism.
'''

with (source / 'all_plotted_values.csv').open(encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 853
roles = {'b':'backbone_profile','c':'sidechain_profile','d':'same_type_environment','e':'angular_occupancy','f':'cross_record_correlations'}
inventory = {}
for row in rows:
    key = '|'.join([roles[row['panel']], row['kind'], row['series'], row['chi'], row['metric']]).rstrip('|')
    item = [row['identity'], float(row['value'])]
    if row['SASA_A2']:
        item.append(float(row['SASA_A2']))
    inventory.setdefault(key, []).append(item)
(root / 'numeric_inventory.json').write_text(json.dumps(inventory, ensure_ascii=False, indent=2), encoding='utf-8')
method += '\nFROZEN NUMERICAL REFERENCE (grouped lists are [identity,value], or [identity,value,SASA] for same-type environments):\n' + json.dumps(inventory, ensure_ascii=False, separators=(',', ':'))
(root / 'scientific_content.txt').write_text(method, encoding='utf-8')

caption = '''Figure 3. Context-sensitive protein flexibility prediction connects backbone constraints, side-chain environments, and angular occupancy. Create a complete new, compact, coherent scientific results figure from the supplied content. The neighboring CaM-CaMKI and ubiquitin-UBA molecular contexts introduce complementary absolute RMSF and Hchi profiles. Include all six Gln environments, the complete chi1-chi4 scatter family, and record-level occupancy correlations across83 records. Use the scientific palette, white background, readable black text, genuine statistics, and no tables, gray microtext, local zooms, mechanism cartoons or decorative architecture. No old figure has been supplied. Design a new composition rather than recreating a familiar layout. Do not print an overall figure title or caption. This is a design concept for later data-faithful native production.'''
(root / 'visual_intent.txt').write_text(caption, encoding='utf-8')

# Preserve author agents/prompts; use a separate experiment work directory for
# the user-specific journal style. Reference images remain their exact originals.
ref_source = repo / 'data' / 'PaperBananaBench' / 'diagram'
references = json.loads((ref_source / 'ref.json').read_text(encoding='utf-8'))
assert len(references) == 298
for record in references:
    path = (ref_source / record['path_to_gt_image']).resolve()
    assert path.is_file()
    record['path_to_gt_image'] = path.as_posix()
ref_dest = root / 'data' / 'PaperBananaBench' / 'diagram'
ref_dest.mkdir(parents=True, exist_ok=True)
(ref_dest / 'ref.json').write_text(json.dumps(references, ensure_ascii=False), encoding='utf-8')
style_dir = root / 'style_guides'
style_dir.mkdir(exist_ok=True)
(style_dir / 'neurips2025_diagram_style_guide.md').write_text('''# Journal-specific style for this experiment
The filename is the author pipeline convention. The actual target is Nature Communications / Genome Biology scientific RESULTS figures, not AI architecture diagrams.
Use a pure white canvas and black sans-serif text. No cards, colored boxes, icon metaphors, decorative arrows, tables, mechanism cartoons or gray microtext. Arrange evidence coherently with asymmetric areas determined by information density, fine aligned axes, compact shared legends, and clear panel lettering. Keep two complete molecular ribbon contexts adjacent, visually substantial, and connected by spatial proximity to quantitative evidence. Do not make empty plot slots. Both sequence profiles use absolute values from zero. Keep the full four-angle scatter family readable, including weak chi4; all six Gln environments and83 record-level distributions remain present. Use analyzed-chain purple #8B72A5, partnergray #BAC1C7, keysite teal #3A8D8C, Fullblue #2166AC, Intrabrown #B87936, MDblack #222222. Methods versus molecular components have distinct color semantics. Brief titles, axes, necessary legend symbols and selected exact statistics only. No overall title inside the image. Target180mm width, about200mm height, text equivalent to at least8pt at164mm manuscript width. Do not add scientific claims beyond the input. Improve evidence relationships and visual hierarchy, not ornamental finish.
''', encoding='utf-8')

wrapper = (tools / 'run_paperbanana.py').read_text(encoding='utf-8')
assert 'work_dir=args.repo,' in wrapper and 'max_concurrent=1' in wrapper
wrapper = wrapper.replace('work_dir=args.repo,', 'work_dir=Path(__file__).resolve().parent,')
wrapper = wrapper.replace('reference_file = args.repo / "data"', 'reference_file = Path(__file__).resolve().parent / "data"')
wrapper = wrapper.replace('max_concurrent=1', 'max_concurrent=2')
wrapper = wrapper.replace('caption=caption, content=method, visual_intent=caption,', '''caption=caption, content=method, visual_intent=caption + ("\\nIndependent composition A: invent your strongest compact unified evidence composition, with asymmetric data-first hierarchy and no fixed inherited grid." if index == 0 else "\\nIndependent composition B: invent a distinct spatial organization with an especially strong coherent side-chain evidence region and a complementary backbone case. Keep both structures adjacent; do not let any single sparse plot dominate the page."),''')
wrapper = wrapper.replace('source_commit=TOOL_COMMIT, entrypoint=', 'old_figure_supplied=False, local_style_override=True, source_commit=TOOL_COMMIT, entrypoint=')
(root / 'run_from_scratch.py').write_text(wrapper, encoding='utf-8')
files = [source / 'all_plotted_values.csv', tools / 'run_paperbanana.py', root / 'scientific_content.txt', root / 'visual_intent.txt', root / 'numeric_inventory.json', ref_source / 'ref.json', style_dir / 'neurips2025_diagram_style_guide.md', root / 'run_from_scratch.py']
manifest = {'workflow':'author demo_full, auto retrieval, independent concepts, no old complete figure',
    'source_images_supplied':False, 'plot_numbers_in_inventory':853, 'reference_examples':len(references),
    'runner_adaptations':['experiment work_dir and journal style','absolute paths to unchanged reference images','two concurrent independent composition intents','no author agents/prompts modified'],
    'files':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in files]}
(root / 'input_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'prepared':str(root),'content_characters':len(method),'reference_count':len(references),'old_figure_supplied':False},ensure_ascii=False))
