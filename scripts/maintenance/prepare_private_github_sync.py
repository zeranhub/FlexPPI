"""Prepare a source-and-evidence backup manifest; does not upload or delete files.

Run from any directory. Generated manifests list original repository-relative
paths, sizes and SHA256 hashes. Runtime environments and bulk simulations stay
local. Git staging/pushing is a separate, explicitly reviewed operation.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
CODE = {'.py', '.pml', '.ps1', '.sh', '.bash', '.bat', '.cmd', '.r', '.R',
        '.js', '.ts', '.mjs', '.ipynb', '.yml', '.yaml', '.toml', '.ini',
        '.cfg', '.md', '.rst', '.bib', '.tex', '.json', '.csv', '.tsv', '.txt',
        '.drawio', '.cif', '.pdb', '.fasta', '.fa', '.npy', '.npz', '.lock'}
PRUNE = {'.git', '__pycache__', 'node_modules', 'site-packages', '.pytest_cache',
         '.ruff_cache', '.mypy_cache', '.cache', 'dependencies', 'vendor',
         'third_party', 'third-party', '.venv', 'venv', 'env', 'Lib', 'Include',
         'include', 'bin', 'share', 'dist', 'build', '_trajectory_staging',
         'plot_runtime', 'runtime', 'runtime_py312', '待删除', 'artifacts',
         'output', 'tmp', 'checkpoints', 'ideaspark_run', '.codex-ppt-inspect'}
PRUNE.update({'tool_runtime', 'deps', 'upstream_archive', 'upstream_prody'})
PRUNE_PATHS = (
    'experiments/_sources/local_model_benchmark_20261023/fastprotflex/records',
    'experiments/_sources/local_model_benchmark_20261023/backflip/inputs',
    'experiments/_sources/local_model_benchmark_20261023/backflip/results',
    'experiments/_sources/energy_landscape_v21_20261012/trajectory_analysis/cache',
    'experiments/_sources/baselines/cohort/chains',
    'experiments/_sources/baselines/cohort/structures',
    'experiments/_sources/mel5_tcr_binding_20261028/ergo/download',
)
QA = {'render', 'renders', 'qa_render', 'rendered', 'render_pages'}
BINARY = {'.png', '.pdf', '.svg', '.docx', '.pse', '.pptx'}
CORE_TEXT = CODE - {'.pdb', '.cif', '.npy', '.npz', '.csv', '.tsv', '.json'}
FOCAL_CASES = ('experiments/_sources/bcr_case_choice_20261028/',
               'experiments/_sources/mel5_tcr_binding_20261028/',
               'experiments/_sources/mel5_absolute_affinity_20261028/')
EXCLUDED_FILES = {
    # Third-party HTML cache contains a website API identifier; citation metadata
    # and our analysis notes are retained, without republishing that raw page.
    'experiments/_sources/bcr_case_choice_20261028/reference/robo1_patent.txt',
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def selected_export(rel: str) -> bool:
    name = Path(rel).name
    # Current manuscript and reviewed Figure 4 sources/exports.
    if name == 'paper_v129_20261029_中文修订_抗原识别与抗体动态.docx':
        return True
    if name.startswith(('figure4_v124_20261027_紧凑蓝橙灰.',
                        'figure4_20261029_TCR与抗体合并预览_v3.',
                        'figure4e_20261029_局部结构重绘_v3.',
                        'figure4e_g_20261029_抗体动态预览_v3.',
                        'figure4a_e_20261029_风格对照_v3.')):
        return True
    if 'bcr_case_choice_20261028/figure4_extension_20261029/revision_e_v3/' in rel:
        return True
    if '/revision_v124_figure4_compact_palette/structure/' in rel:
        return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='docs/backups/2026-10-29')
    parser.add_argument('--hash', action='store_true')
    args = parser.parse_args()
    target = (ROOT / args.output).resolve()
    target.relative_to(ROOT)
    target.mkdir(parents=True, exist_ok=True)
    tracked_raw = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT)
    tracked = set(tracked_raw.decode('utf-8').split('\0')) - {''}
    selected, excluded = [], collections.Counter()
    skipped_dirs = []
    for current, directories, files in os.walk(ROOT):
        parent = Path(current)
        for directory in directories[:]:
            child = parent / directory
            relative = child.relative_to(ROOT).as_posix()
            skip = (relative in PRUNE_PATHS or directory in PRUNE or directory.startswith(('.venv-', 'runtime_', 'mmseqs_',
                    '.render_check', '_qa_render')) or child.is_symlink() or child.is_junction()
                    or (child / 'pyvenv.cfg').exists()
                    or ('paper' in child.relative_to(ROOT).parts and directory in QA))
            if skip:
                directories.remove(directory)
                skipped_dirs.append(relative)
        for name in files:
            path = parent / name
            rel = path.relative_to(ROOT).as_posix()
            if rel in EXCLUDED_FILES:
                excluded['third_party_page_with_embedded_key'] += 1
                continue
            if path.is_symlink() or target == parent or target in parent.parents:
                continue
            lower = name.lower()
            if (lower.startswith('.env') and lower != '.env.example') or lower.endswith(('.pem', '.p12', '.pfx', '.key')):
                excluded['credential_filename'] += 1
                continue
            size = path.stat().st_size
            suffix = path.suffix
            special_evidence = rel == 'experiments/_sources/bcr_case_choice_20261028/hcd3_functional_analysis_20261029/md_audit/all_h3_partner_minimum_distances.npz'
            if size >= (50 if special_evidence else 25) * 1024 * 1024:
                excluded['file_at_least_25_MiB'] += 1
                continue
            parts = Path(rel).parts
            # Bulk training structures are recoverable externally; retain the
            # exact canonical structure required by the current 6A77 case.
            if '/training_payload/dataset/' in '/' + rel and name != 'A00LZ.pdb':
                excluded['bulk_training_data'] += 1
                continue
            if parts[0] in {'experiments', 'paper'}:
                focal = rel.startswith(FOCAL_CASES)
                if suffix in CORE_TEXT:
                    keep = True
                elif suffix in {'.csv', '.tsv', '.json'}:
                    batch = any(part in {'predictions', 'table_parts', 'cache',
                                'public_candidates', 'prediction_records'} for part in parts)
                    keep = focal or (size < 2 * 1024 * 1024 and not batch)
                elif suffix in {'.pdb', '.cif', '.npy', '.npz'}:
                    keep = (focal or name == 'A00LZ.pdb'
                            or '/revision_v124_figure4_compact_palette/structure/' in rel
                            or rel in tracked)
                else:
                    keep = suffix in BINARY and selected_export(rel)
                if suffix == '.html':
                    keep = selected_export(rel)
            else:
                keep = suffix in CODE or rel in tracked or parent == ROOT
            if name in {'LICENSE', 'COPYING', 'NOTICE'}:
                keep = True
            if not keep:
                excluded['unselected_generated_or_bulk_format'] += 1
                continue
            item = {'path': rel, 'bytes': size, 'sha256': digest(path) if args.hash else ''}
            selected.append(item)
    selected.sort(key=lambda item: item['path'])
    with (target / 'selected_files.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['path', 'bytes', 'sha256'])
        writer.writeheader(); writer.writerows(selected)
    (target / 'selected_paths.nul').write_bytes(b''.join(item['path'].encode('utf-8') + b'\0' for item in selected))
    groups = collections.Counter()
    totals = collections.Counter()
    for item in selected:
        group = item['path'].split('/')[0]
        groups[group] += 1; totals[group] += item['bytes']
    summary = {'date': '2026-10-29', 'repository': 'zeranhub/FlexPPI',
               'required_visibility': 'PRIVATE', 'files': len(selected),
               'bytes': sum(item['bytes'] for item in selected),
               'sha256_computed': args.hash, 'groups': dict(groups),
               'group_bytes': dict(totals), 'excluded_file_counts': dict(excluded),
               'excluded_directories': sorted(skipped_dirs),
               'limit': 'Source and evidence backup, not a complete MD/checkpoint archive.'}
    (target / 'selection_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({key: value for key, value in summary.items() if key != 'excluded_directories'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
