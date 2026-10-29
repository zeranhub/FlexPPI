"""Stage the reviewed manifest and existing local changes; never pushes."""
from pathlib import Path
import csv, json, subprocess

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
rows=list(csv.DictReader((HERE/'selected_files.csv').open(encoding='utf-8-sig')))
paths={row['path'] for row in rows}
required={
 'src/flexppi/model/network.py',
 'experiments/core_corrected_pipeline/training_payload/dataset/structures/A00LZ.pdb',
 'experiments/_sources/bcr_case_choice_20261028/6A77_all_positions_predictions_and_MD.csv',
 'experiments/_sources/bcr_case_choice_20261028/hcd3_functional_analysis_20261029/md_audit/all_h3_partner_minimum_distances.npz',
 'experiments/_sources/bcr_case_choice_20261028/figure4_extension_20261029/revision_e_v3/plot_e_v3.py',
 'paper/revision_v14_20261010/paper_v129_20261029_中文修订_抗原识别与抗体动态.docx',
 'paper/revision_v14_20261010/figure4_20261029_TCR与抗体合并预览_v3.pdf',
 'paper/revision_v14_20261010/figure4e_20261029_局部结构重绘_v3.png',
 'paper/revision_v14_20261010/section33_expansion_20261025/model_first_story_search_20261025/revision_v124_figure4_compact_palette/structure/MEL5_affinity_maturation.pse',
}
assert required <= paths, sorted(required-paths)
assert all(int(row['bytes'])<100*1024*1024 for row in rows)
assert not any(path.endswith('reference/robo1_patent.txt') for path in paths)
assert all(row['sha256'] for row in rows)
# This script only modifies the local index. Repository privacy must be checked
# separately immediately before any network push.
commands=[
 ['git','-c','gc.auto=0','--literal-pathspecs','add','-f','--pathspec-from-file=docs/backups/2026-10-29/selected_paths.nul','--pathspec-file-nul'],
 ['git','-c','gc.auto=0','add','-u'],
 ['git','-c','gc.auto=0','add','--','docs/backups/2026-10-29/selected_files.csv','docs/backups/2026-10-29/selection_summary.json','docs/backups/2026-10-29/README.md','docs/backups/2026-10-29/stage_reviewed_snapshot.py'],
]
for command in commands:
    result=subprocess.run(command,cwd=ROOT,capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode('utf-8','replace')[-3000:])
    print(json.dumps({'step':command[command.index('add')+1:][:1],'status':'staged','diagnostic_lines':len(result.stderr.splitlines())}),flush=True)
print(subprocess.check_output(['git','diff','--cached','--shortstat'],cwd=ROOT).decode())
