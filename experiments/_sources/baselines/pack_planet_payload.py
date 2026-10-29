"""Package verified local source/weights and minimum offline ESM3 snapshot."""
import hashlib
import json
from pathlib import Path
import tarfile

OUTPUT = Path(__file__).resolve().parent
MODELS = Path(r'E:\Project\protein_flexibility_models')
ESM_ID = '47f0545b2b6daf26a93439a3cd610f4f7f3d5478'
ROCKET_ID = '7831d318182d4f4fc1acbdafdbf23862d8bc1657'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    included = []
    source = MODELS/'PLANET-MD'
    for path in sorted((source/'planet_md').rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts:
            included.append((path, 'PLANET-MD/'+path.relative_to(source).as_posix()))
    for filename in ['README.md','LICENSE','pyproject.toml']:
        included.append((source/filename, 'PLANET-MD/'+filename))
    esm_repo = MODELS/'hf_cache/hub/models--EvolutionaryScale--esm3-sm-open-v1'
    for path in sorted((esm_repo/'snapshots'/ESM_ID).rglob('*')):
        if not path.is_file():
            continue
        if path.stat().st_size > 100000000 and path.name != 'esm3_sm_open_v1.pth':
            continue
        included.append((path,'hf_cache/hub/models--EvolutionaryScale--esm3-sm-open-v1/'+path.relative_to(esm_repo).as_posix()))
    included.append((esm_repo/'refs/main','hf_cache/hub/models--EvolutionaryScale--esm3-sm-open-v1/refs/main'))
    rocket = MODELS/'hf_cache/hub/models--samsl--rocketshp/snapshots'/ROCKET_ID/'checkpoints'
    for filename in ['rocketshp_pretrained_20250505_v1.ckpt','rocketshp_pretrained_20250507_v1_seq.ckpt']:
        included.append((rocket/filename,'weights/'+filename))
    inventory = []
    with tarfile.open(OUTPUT/'planet_payload.tar', 'w', dereference=True) as archive:
        for number, (path, name) in enumerate(included):
            inventory.append(dict(path=name, bytes=path.stat().st_size, sha256=sha(path)))
            archive.add(path, arcname=name, recursive=False)
            if path.stat().st_size > 100000000:
                print('Packaged', name, path.stat().st_size, flush=True)
    report = dict(status='frozen_payload', upstream_planet_commit='917999482343c6400bef51180c817a8af01ee871',
        esm_snapshot=ESM_ID, rocket_snapshot=ROCKET_ID, files=inventory,
        bytes=sum(item['bytes'] for item in inventory), payload_sha256=sha(OUTPUT/'planet_payload.tar'),
        policy='Minimal offline snapshot contains sequence model and structure encoder plus all small metadata; unused function/structure decoder weights omitted. Set HF_HUB_OFFLINE=1 and HF_HUB_CACHE to extracted hf_cache/hub. No API token required for already-cached weights.')
    (OUTPUT/'planet_payload_freeze.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Complete',report['bytes'],report['payload_sha256'],flush=True)


if __name__ == '__main__':
    main()
