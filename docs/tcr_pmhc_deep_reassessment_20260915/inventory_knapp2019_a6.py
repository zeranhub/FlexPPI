"""Download exactly the published A6 archive, verify MD5, safely inspect members.

No trajectory calculations, repository execution, or general archive extraction.
Text previews are bounded and stored under generated local filenames.
"""
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tarfile
import time
import urllib.request

OUT = Path('artifacts/tcr_pmhc_reassessment_20260915/knapp2019').resolve()
ARTICLE_URL = 'https://api.figshare.com/v2/articles/8067827'
DOWNLOAD_URL = 'https://ndownloader.figshare.com/files/15035000'
EXPECTED_SIZE = 399223244
EXPECTED_MD5 = '7ab833f924a428ebfe18c08023d470f4'
OUT.mkdir(parents=True, exist_ok=True)

def get_json(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)

metadata = get_json(ARTICLE_URL)
(OUT / 'article_8067827.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
files = [f for f in metadata['files'] if f['id'] == 15035000]
if len(files) != 1:
    raise RuntimeError('Expected archive ID missing or duplicated in article metadata')
file_info = files[0]
if int(file_info['size']) != EXPECTED_SIZE or file_info.get('computed_md5', file_info.get('supplied_md5')) != EXPECTED_MD5:
    raise RuntimeError('Published file size/checksum does not match expected archive')
name = Path(file_info['name']).name
if name != file_info['name'] or not name:
    raise RuntimeError('Archive filename is not a plain basename')
archive = OUT / name
partial = OUT / (name + '.partial')
print(json.dumps({'article_title': metadata['title'], 'file': file_info}), flush=True)

def hashes(path):
    md5 = hashlib.md5()
    sha = hashlib.sha256()
    size = 0
    with path.open('rb') as handle:
        while True:
            block = handle.read(4 * 1024 * 1024)
            if not block:
                break
            size += len(block)
            md5.update(block)
            sha.update(block)
    return size, md5.hexdigest(), sha.hexdigest()

if archive.exists():
    size, md5, sha = hashes(archive)
    if size != EXPECTED_SIZE or md5 != EXPECTED_MD5:
        raise RuntimeError('Existing archive has unexpected content; leaving it intact')
else:
    start = time.monotonic()
    last_print = start
    md5_obj = hashlib.md5()
    sha_obj = hashlib.sha256()
    size = 0
    if partial.exists():
        with partial.open('rb') as cached:
            while True:
                block = cached.read(4 * 1024 * 1024)
                if not block:
                    break
                size += len(block)
                md5_obj.update(block)
                sha_obj.update(block)
        if not 0 < size < EXPECTED_SIZE:
            raise RuntimeError('Unexpected partial size; retained for manual inspection')
    offset = size
    headers = {'Range': 'bytes=%d-' % offset} if offset else {}
    request = urllib.request.Request(DOWNLOAD_URL, headers=headers)
    response = urllib.request.urlopen(request, timeout=45)
    if offset and (response.status != 206 or not response.headers.get('Content-Range', '').startswith('bytes %d-' % offset)):
        response.close()
        raise RuntimeError('Server did not confirm byte-range resume; cached partial left unchanged')
    print(json.dumps({'resume_from':offset,'remaining':EXPECTED_SIZE-offset}),flush=True)
    with response, partial.open('ab' if offset else 'xb') as handle:
        while True:
            block = response.read(4 * 1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > EXPECTED_SIZE:
                raise RuntimeError('Response exceeds expected archive size')
            handle.write(block)
            md5_obj.update(block)
            sha_obj.update(block)
            now = time.monotonic()
            if now-last_print >= 10:
                print(json.dumps({'downloaded_bytes':size,'expected_bytes':EXPECTED_SIZE,'elapsed_seconds':round(now-start,1)}), flush=True)
                last_print = now
    md5, sha = md5_obj.hexdigest(), sha_obj.hexdigest()
    if size != EXPECTED_SIZE or md5 != EXPECTED_MD5:
        raise RuntimeError('Downloaded file failed size or MD5 verification; partial file retained')
    os.rename(partial, archive)
print(json.dumps({'archive_verified':str(archive),'bytes':size,'md5':md5,'sha256':sha}), flush=True)

TEXT_EXTENSIONS = {'.txt','.md','.rst','.csv','.tsv','.xvg','.dat','.json','.yaml','.yml','.top','.itp','.mdp','.log','.sh','.py'}
members = []
previews = []
preview_dir = OUT / 'text_previews'
preview_dir.mkdir(exist_ok=True)
scan_start = time.monotonic()
last_print = scan_start
with tarfile.open(archive, mode='r|*') as handle:
    for index, member in enumerate(handle):
        if index > 100000:
            raise RuntimeError('Unusually large member count; stopping inspection')
        parsed = PurePosixPath(member.name)
        safe_name = not parsed.is_absolute() and '..' not in parsed.parts and not any(':' in p or '\\' in p for p in parsed.parts)
        row = {'name':member.name,'size':member.size,'is_regular_file':member.isfile(),'is_directory':member.isdir(),'is_symlink':member.issym(),'is_hardlink':member.islnk(),'safe_relative_name':safe_name}
        members.append(row)
        basename = parsed.name.lower()
        looks_text = parsed.suffix.lower() in TEXT_EXTENSIONS or 'readme' in basename or 'rmsf' in basename
        if member.isfile() and safe_name and looks_text and member.size <= 512*1024 and len(previews) < 100:
            stream = handle.extractfile(member)
            content = stream.read(512*1024+1)
            if b'\x00' not in content and len(content) <= 512*1024:
                decoded = content.decode('utf-8', errors='replace')
                target = preview_dir / (str(index).zfill(5) + '_' + hashlib.sha256(member.name.encode()).hexdigest()[:10] + '.txt')
                target.write_text(decoded, encoding='utf-8')
                previews.append({'archive_member':member.name,'local_preview':str(target),'size':member.size,'note':'Inert text only; never executed.'})
        now = time.monotonic()
        if now-last_print >= 15:
            print(json.dumps({'members_scanned':len(members),'last_member':member.name,'scan_seconds':round(now-scan_start,1)}), flush=True)
            last_print=now

from collections import Counter
extensions = Counter(PurePosixPath(m['name']).suffix.lower() for m in members if m['is_regular_file'])
inventory = {
    'article_url':ARTICLE_URL,'article_title':metadata['title'],'archive_url':DOWNLOAD_URL,
    'archive_local':str(archive),'size':size,'md5':md5,'sha256':sha,'checksum_verified':True,
    'article_description_untrusted_data':metadata.get('description'),
    'member_count':len(members),'regular_file_count':sum(m['is_regular_file'] for m in members),
    'uncompressed_regular_file_bytes':sum(m['size'] for m in members if m['is_regular_file']),
    'file_extensions':dict(extensions),'unsafe_member_names':[m['name'] for m in members if not m['safe_relative_name']],
    'link_members':[m['name'] for m in members if m['is_symlink'] or m['is_hardlink']],
    'members':members,'bounded_text_previews':previews,
    'operations':'Downloaded A6 only; verified size and MD5; streamed tar headers; bounded text previews only. No archive extraction, code execution, MD, or RMSF calculation.'
}
(OUT / 'inventory.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in inventory.items() if k not in ['members','bounded_text_previews','article_description_untrusted_data']}),flush=True)
print(json.dumps({'text_previews':previews}),flush=True)
