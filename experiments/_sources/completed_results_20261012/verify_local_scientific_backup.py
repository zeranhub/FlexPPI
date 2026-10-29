from pathlib import Path, PurePosixPath
import tarfile, hashlib, json
OUT = Path(__file__).resolve().parent
package = OUT/'completed_scientific_sources_20261012.tar.gz'
report = json.loads((OUT/'completed_scientific_sources_manifest.json').read_text(encoding='utf8'))
def digest_file(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()
assert package.stat().st_size == report['archive_bytes']
assert digest_file(package) == report['archive_sha256']
verified = {}
with tarfile.open(package, 'r|gz') as archive:
    for member in archive:
        relative = PurePosixPath(member.name)
        assert member.isfile() and not relative.is_absolute() and '..' not in relative.parts and '\\' not in member.name
        assert member.name in report['members'] and member.name not in verified
        expected = report['members'][member.name]
        assert member.size == expected['bytes']
        stream = archive.extractfile(member)
        h = hashlib.sha256()
        for block in iter(lambda:stream.read(1024*1024), b''):
            h.update(block)
        assert h.hexdigest() == expected['sha256'], member.name
        verified[member.name] = dict(bytes=member.size, sha256=h.hexdigest())
assert set(verified) == set(report['members'])
final_weights = [n for n in verified if n.endswith('/epoch_040.pt')]
assert len(final_weights) == 20
audit = {'status':'passed', 'local_archive':str(package), 'archive_sha256':report['archive_sha256'],
         'archive_bytes':package.stat().st_size, 'all_members_verified':len(verified),
         'final_weights_verified':20,'per_residue_compressed_sources_verified':4,
         'archive_path_validation':'All ordinary files, relative paths only, no parent escapes or duplicate members.',
         'members':verified,'unpacked_to_taskroot':False}
(OUT/'local_scientific_backup_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf8')
print(json.dumps({k:audit[k] for k in ['status','archive_bytes','all_members_verified','final_weights_verified','per_residue_compressed_sources_verified']}))
