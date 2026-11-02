"""Verify current-paper publication contents without running scientific calculations."""
from pathlib import Path
from zipfile import ZipFile
from urllib.parse import unquote
import ast
import hashlib
import json
import re
import subprocess


ROOT = Path(__file__).resolve().parents[2]


def verify():
    try:
        repo = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()).resolve()
    except (subprocess.CalledProcessError, FileNotFoundError):
        repo = None
    if repo == ROOT:
        raw = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT)
        tracked = set(raw.decode('utf-8').split('\0')) - {''}
    else:
        # Also verify downloaded source archives and isolated export directories.
        tracked = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                   if p.is_file() and not any(x in p.parts for x in ['.git', '__pycache__', '.pytest_cache'])}
    for rel in tracked:
        assert (ROOT / rel).is_file(), 'Missing tracked file: ' + rel
        assert '.egg-info/' not in rel, rel
        if rel.startswith('paper/'):
            assert rel == 'paper/README.md' or rel.startswith('paper/current/'), rel
        assert not re.search(r'/(?:kix[^/]*|bcr_[^/]*|kras_[^/]*|protac[^/]*|pdz3[^/]*)/', rel, re.I), rel
        assert '/oof_analysis/' not in rel, rel
        if rel.endswith('.py'):
            ast.parse((ROOT / rel).read_text(encoding='utf-8-sig'), filename=rel)

    paper = ROOT / 'paper/current'
    manifest = json.loads((paper / 'manifest.json').read_text(encoding='utf-8'))
    assert manifest['manuscript_version'] == 'v136'
    for record in manifest['files']:
        path = paper / record['path']
        assert path.relative_to(ROOT).as_posix() in tracked, str(path)
        data = path.read_bytes()
        assert len(data) == record['bytes'], str(path)
        assert hashlib.sha256(data).hexdigest() == record['sha256'], str(path)

    with ZipFile(paper / 'manuscript.docx') as doc:
        images = [hashlib.sha256(doc.read(p)).hexdigest() for p in doc.namelist() if p.startswith('word/media/')]
    names = ['Figure1', 'Figure2', 'Figure3', 'Figure4', 'FigureS1', 'FigureS2', 'FigureS3', 'FigureS4', 'FigureS5']
    hashes = [hashlib.sha256((paper / 'figures' / name / 'figure.png').read_bytes()).hexdigest() for name in names]
    assert images == hashes and len(images) == 9

    registry = json.loads((ROOT / 'experiments/steps.json').read_text(encoding='utf-8'))
    assert set(registry) == {'performance', 'environment', 'labels', 'training', 'family', 'occupancy', 'mel5', 'baselines'}
    checked_steps = 0
    for group in registry.values():
        for record in group['steps'].values():
            for rel in [record['script'], *record.get('inputs', [])]:
                path = (ROOT / 'experiments' / rel).resolve()
                assert path.relative_to(ROOT).as_posix() in tracked, str(path)
            checked_steps += 1

    docs = ['README.md', 'BACKUP_MANIFEST.md', 'experiments/README.md', 'docs/PROJECT_LAYOUT.md', 'docs/REPRODUCIBILITY.md', 'docs/REPOSITORY_SCOPE.md', 'paper/README.md', 'paper/current/README.md', 'paper/current/figures/README.md']
    links = 0
    for rel in docs:
        path = ROOT / rel
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text(encoding='utf-8')):
            if re.match(r'[a-z]+://|#', target):
                continue
            resolved = (path.parent / unquote(target.split('#')[0])).resolve()
            key = resolved.relative_to(ROOT).as_posix()
            assert key in tracked or any(f.startswith(key.rstrip('/') + '/') for f in tracked), (rel, target)
            links += 1
    return dict(status='passed', tracked_files=len(tracked), paper_figures=4, supplementary_figures=5,
                manuscript_images_identical=9, manifest_files_checked=len(manifest['files']),
                registered_steps_checked=checked_steps, entry_links_checked=links)


if __name__ == '__main__':
    print(json.dumps(verify(), ensure_ascii=False, indent=2))
