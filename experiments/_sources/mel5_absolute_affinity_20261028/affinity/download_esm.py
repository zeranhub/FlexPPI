"""Fetch the public frozen ESM checkpoint with verified HTTP byte ranges."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib
import json
import re
import time

import requests

ROOT = Path(__file__).resolve().parent
DEST = ROOT / "vendor" / "esm2_t33_650M_UR50D"
URL = "https://huggingface.co/facebook/esm2_t33_650M_UR50D/resolve/08e4846e537177426273712802403f7ba8261b6c/model.safetensors"
EXPECTED_SHA = "a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0"
PARTS = DEST / "download_parts"
PARTS.mkdir(exist_ok=True)
response = requests.get(URL, headers={"Range": "bytes=0-0"}, timeout=(10, 40))
response.raise_for_status()
size = int(response.headers["Content-Range"].split("/")[-1])
chunk_size = 64 * 1024 * 1024
chunks = [(start, min(start + chunk_size, size) - 1) for start in range(0, size, chunk_size)]


def fetch(pair):
    start, stop = pair
    path = PARTS / f"{start:012d}.part"
    expected_size = stop - start + 1
    if path.exists() and path.stat().st_size == expected_size:
        return path
    for attempt in range(3):
        try:
            result = requests.get(URL, headers={"Range": f"bytes={start}-{stop}"}, timeout=(15, 60))
            assert result.status_code == 206, result.status_code
            assert result.headers["Content-Range"] == f"bytes {start}-{stop}/{size}"
            assert len(result.content) == expected_size
            path.write_bytes(result.content)
            return path
        except Exception:
            if attempt == 2:
                raise
    raise RuntimeError("unreachable")


started = time.monotonic()
with ThreadPoolExecutor(max_workers=4) as pool:
    futures = [pool.submit(fetch, chunk) for chunk in chunks]
    for i, future in enumerate(as_completed(futures), 1):
        future.result()
        print(f"Downloaded {i}/{len(chunks)} chunks in {time.monotonic()-started:.1f}s", flush=True)
target = DEST / "model.safetensors"
digest = hashlib.sha256()
with target.open("wb") as output:
    for start, _ in chunks:
        part = PARTS / f"{start:012d}.part"
        with part.open("rb") as handle:
            while buffer := handle.read(4 * 1024 * 1024):
                digest.update(buffer)
                output.write(buffer)
assert target.stat().st_size == size
assert digest.hexdigest() == EXPECTED_SHA, digest.hexdigest()
for start, _ in chunks:
    (PARTS / f"{start:012d}.part").unlink()
(ROOT / "esm_download_audit.json").write_text(json.dumps({
    "url": URL, "bytes": size, "sha256": digest.hexdigest(),
    "revision": "08e4846e537177426273712802403f7ba8261b6c",
    "elapsed_seconds": time.monotonic() - started,
}, indent=2), encoding="utf-8")
print("Verified ESM checkpoint", digest.hexdigest(), flush=True)
