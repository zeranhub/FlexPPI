"""Download KIX pilot MD trajectories from the remote server with resume and checks.

The SSH password is requested interactively and never stored.  Files are written
to ``<local-root>/<state>/production/segment_XXXX.dcd`` plus a topology PDB per
state, skipping any local file whose size already matches the remote file.
When the locally saved analysis summaries contain SHA256 values, every downloaded
file is verified against them.

Example
-------
python scripts/kix_md/fetch_trajectories.py --host 172.16.196.191 --user ranze ^
    --remote-root /ranze/runs/flexppi_kix_md_20260914 ^
    --local-root artifacts/external/kix_md_20260914
"""

from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import sys
from pathlib import Path

import paramiko

ROOT = Path(__file__).resolve().parents[2]
STATES = ("apo", "myb", "mll", "both")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_hashes(state: str) -> dict:
    summary = ROOT / "artifacts/kix_md_20260914/pilot_results_20260916" / state / "summary.json"
    if not summary.exists():
        return {}
    data = json.loads(summary.read_text(encoding="utf-8"))
    output = {str(row["filename"]): str(row["sha256"]) for row in data.get("segments", [])}
    if data.get("protein_pdb_sha256"):
        output["protein.pdb"] = str(data["protein_pdb_sha256"])
    return output


def remote_size(sftp: paramiko.SFTPClient, path: str) -> int:
    try:
        return int(sftp.stat(path).st_size)
    except FileNotFoundError:
        return -1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", required=True)
    parser.add_argument("--port", type=int, default=22)
    parser.add_argument("--user", required=True)
    parser.add_argument("--remote-root", default="/ranze/runs/flexppi_kix_md_20260914")
    parser.add_argument("--local-root", default=str(ROOT / "artifacts/external/kix_md_20260914"))
    parser.add_argument("--states", nargs="+", default=list(STATES))
    parser.add_argument("--replica", default="rep01")
    parser.add_argument("--segments", type=int, default=100)
    parser.add_argument("--verify-sha", action="store_true", default=True)
    parser.add_argument("--no-verify-sha", dest="verify_sha", action="store_false")
    parser.add_argument("--password-stdin", action="store_true",
                        help="read the password from one stdin line instead of prompting")
    parser.add_argument("--list", action="store_true",
                        help="only report remote file counts and sizes, then exit")
    args = parser.parse_args()

    local_root = Path(args.local_root)
    if args.password_stdin:
        password = sys.stdin.readline().rstrip("\r\n")
    else:
        password = getpass.getpass("SSH password (not saved): ")
    client = paramiko.SSHClient()
    client.load_host_keys(str(Path.home() / ".ssh" / "known_hosts"))
    client.set_missing_host_key_policy(paramiko.RejectPolicy())
    client.connect(args.host, port=args.port, username=args.user, password=password,
                   look_for_keys=False, allow_agent=False, timeout=20, auth_timeout=30)
    password = None
    client.get_transport().set_keepalive(30)
    sftp = client.open_sftp()

    failures = []
    verified = 0
    downloaded = 0
    skipped = 0
    if args.list:
        total_bytes = 0
        missing = []
        for state in args.states:
            for index in range(args.segments):
                remote_path = "{}/runs/{}/{}/production/segment_{:04d}.dcd".format(
                    args.remote_root, state, args.replica, index)
                size = remote_size(sftp, remote_path)
                if size < 0:
                    missing.append(remote_path)
                else:
                    total_bytes += size
            topology = "{}/prepared/{}/protein.pdb".format(args.remote_root, state)
            size = remote_size(sftp, topology)
            if size < 0:
                missing.append(topology)
            else:
                total_bytes += size
        sftp.close()
        client.close()
        print(json.dumps({"states": args.states, "segments_per_state": args.segments,
                          "total_gib": round(total_bytes / 1024 ** 3, 3),
                          "missing": missing}, indent=2))
        return
    for state in args.states:
        hashes = expected_hashes(state) if args.verify_sha else {}
        jobs = []
        for index in range(args.segments):
            name = "segment_{:04d}.dcd".format(index)
            jobs.append((state,
                         "{}/runs/{}/{}/production/{}".format(args.remote_root, state, args.replica, name),
                         local_root / state / "production" / name,
                         hashes.get(name)))
        jobs.append((state,
                     "{}/prepared/{}/protein.pdb".format(args.remote_root, state),
                     local_root / state / "protein.pdb",
                     hashes.get("protein.pdb")))
        for _, remote_path, local_path, expected in jobs:
            local_path.parent.mkdir(parents=True, exist_ok=True)
            size = remote_size(sftp, remote_path)
            if size < 0:
                failures.append({"remote": remote_path, "error": "missing on server"})
                continue
            if local_path.exists() and local_path.stat().st_size == size:
                if expected and sha256(local_path) != expected:
                    failures.append({"remote": remote_path, "error": "local hash mismatch"})
                    continue
                skipped += 1
                continue
            part = local_path.with_suffix(local_path.suffix + ".part")
            state_label = "{}".format(local_path.relative_to(local_root))
            last = {"percent": -10}

            def callback(transferred: int, total: int) -> None:
                percent = int(100 * transferred / max(total, 1))
                if percent >= last["percent"] + 10:
                    last["percent"] = percent
                    print("  {} {}%".format(state_label, percent), flush=True)

            sftp.get(remote_path, str(part), callback=callback)
            part.replace(local_path)
            actual = sha256(local_path) if expected else None
            if expected and actual != expected:
                failures.append({"remote": remote_path, "error": "sha256 mismatch",
                                 "expected": expected, "actual": actual})
                continue
            verified += 1 if expected else 0
            downloaded += 1
    sftp.close()
    client.close()
    print(json.dumps({"downloaded": downloaded, "skipped_existing": skipped,
                      "sha256_verified": verified, "failures": failures}, indent=2))
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
