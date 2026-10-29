"""Small bounded HTTP probes and streaming audit; never download full weights."""
import ctypes
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import zipfile

ROOT = Path(__file__).resolve().parent.parent
DIRECTORY = ROOT / "weights/flexpert"
OFFICIAL_SOURCE = ROOT / "vendor/Flexpert-main/download_flexpert_weights.sh"
SOURCE_COMMIT = "66c5adb348b53183a1e75e7f1ad5fdb395b6ba00"
BASE = "https://data.ciirc.cvut.cz/public/projects/2025Flexpert/flexpert-weights/"


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        blocks = 0
        for data in iter(lambda: handle.read(1024**2), b""):
            result.update(data)
            blocks += 1
            if blocks % 64 == 0:
                time.sleep(0.005)
    return result.hexdigest()


def headers(path):
    data = path.read_text(encoding="utf-8-sig")
    result = {}
    for line in data.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            result[key.strip().lower()] = value.strip()
        elif line.startswith("HTTP/"):
            result["http_status_line"] = line
    return result


def file_record(path):
    before = path.stat()
    result = {"path": str(path), "logical_size_bytes": before.st_size,
              "mtime_ns": before.st_mtime_ns, "sha256": digest(path)}
    after = path.stat()
    result["unchanged_during_hash"] = before.st_size == after.st_size and before.st_mtime_ns == after.st_mtime_ns
    assert result["unchanged_during_hash"], str(path)
    return result


def bounded_probe(variant, remote):
    probe = DIRECTORY / f"probe_{variant}_range8MiB_20261024.bin"
    response = DIRECTORY / f"probe_{variant}_range8MiB_20261024.headers.txt"
    metrics = DIRECTORY / f"probe_{variant}_range8MiB_20261024.metrics.json"
    maximum = 8 * 1024**2
    command = ["curl.exe", "--location", "--silent", "--show-error", "--fail",
               "--connect-timeout", "8", "--max-time", "60", "--range", f"0-{maximum-1}",
               "--max-filesize", str(maximum), "--header", "If-Range: " + remote["etag"],
               "--dump-header", str(response), "--output", str(probe), "--write-out", "%{json}",
               BASE + f"flexpert_{variant}_weights.bin"]
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    process = subprocess.run(command, capture_output=True, text=True, timeout=75)
    parsed = json.loads(process.stdout) if process.stdout.strip().startswith("{") else {}
    keep_metrics = {"http_code", "url_effective", "remote_ip", "size_download", "speed_download",
                    "time_namelookup", "time_connect", "time_appconnect", "time_starttransfer",
                    "time_total", "response_code", "exitcode", "errormsg", "num_redirects", "http_version"}
    parsed = {key: value for key, value in parsed.items() if key in keep_metrics}
    saved_headers = headers(response) if response.exists() else {}
    size = probe.stat().st_size if probe.exists() else 0
    valid = (process.returncode == 0 and parsed.get("http_code") == 206
             and saved_headers.get("content-range") == f"bytes 0-{maximum-1}/{remote['content-length']}"
             and saved_headers.get("etag") == remote["etag"] and size == maximum)
    record = {"started_utc": started_utc, "wall_seconds": time.perf_counter() - started,
              "exit_code": process.returncode, "stderr": process.stderr,
              "request": {"URL": command[-1], "range_start": 0, "range_end": maximum-1,
                          "maximum_transfer_bytes": maximum, "maximum_transfer_seconds": 60,
                          "if_range_ETag": remote["etag"]},
              "curl_metrics": parsed, "response_headers": saved_headers,
              "actual_probe_bytes": size, "complete_exact_range_and_ETag_verified": valid,
              "probe_sha256": digest(probe) if probe.exists() else None}
    metrics.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def main():
    # Low process priority; do not contend with ongoing local scientific jobs.
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
    created = datetime.now(timezone.utc).isoformat()
    disk_before = shutil.disk_usage(DIRECTORY)
    local = []
    for name in ("flexpert_3d_weights.bin", "flexpert_3d_weights.bin.aria2", "flexpert_seq_weights.bin",
                 "download_probe.bin", "range_probe_64MB.bin", "download_3d.log"):
        local.append(file_record(DIRECTORY / name))
    remote = {variant: headers(DIRECTORY / f"probe_{variant}_head_20261024.headers.txt") for variant in ("3d", "seq")}
    small_range_checks = []
    for variant in ("3d", "seq"):
        local_path = DIRECTORY / f"flexpert_{variant}_weights.bin"
        probe = DIRECTORY / f"probe_{variant}_range64k_20261024.bin"
        with local_path.open("rb") as handle:
            prefix = handle.read(65536)
        small_range_checks.append({"variant": variant, "range": "0-65535",
                                   "response_headers": headers(probe.with_name(probe.stem + ".headers.txt")),
                                   "probe_file": file_record(probe), "saved_prefix_equal_remote_probe": prefix == probe.read_bytes()})
    seq_path = DIRECTORY / "flexpert_seq_weights.bin"
    tail_probe = DIRECTORY / "probe_seq_prefix_tail64k_20261024.bin"
    with seq_path.open("rb") as handle:
        handle.seek(seq_path.stat().st_size - 65536)
        tail = handle.read(65536)
    small_range_checks.append({"variant": "seq", "range": "11948032-12013567", "probe_file": file_record(tail_probe),
                               "response_headers": headers(tail_probe.with_name(tail_probe.stem + ".headers.txt")),
                               "saved_prefix_tail_equal_remote_probe": tail == tail_probe.read_bytes()})
    probes = {variant: bounded_probe(variant, info) for variant, info in remote.items()}
    partial = {}
    for variant in ("3d", "seq"):
        path = DIRECTORY / f"flexpert_{variant}_weights.bin"
        partial[variant] = {
            "status": "partial_not_loadable_as_complete_weights", "path": str(path),
            "local_logical_size_bytes": path.stat().st_size, "remote_total_bytes": int(remote[variant]["content-length"]),
            "aria2_control_file_present": path.with_name(path.name + ".aria2").exists(),
            "zip_central_directory_complete": zipfile.is_zipfile(path),
            "warning": "Logical file size and local partial SHA are not a completed-download percentage or an official full-weight digest. Aria2 can leave noncontiguous pieces and zero-filled holes.",
            "official_remote_SHA256_available": False,
            "source_commit": SOURCE_COMMIT, "remote_ETag": remote[variant]["etag"],
        }
        (DIRECTORY / f"flexpert_{variant}_weights.bin.partial.json").write_text(json.dumps(partial[variant], indent=2), encoding="utf-8")
    audit = {
        "created_utc": created, "source_script": file_record(OFFICIAL_SOURCE), "source_repository_commit": SOURCE_COMMIT,
        "scope": "Read-only local weight inspection and bounded official HTTP probes; no full download, model load, inference, or ongoing-job changes.",
        "local_files": local, "official_HEAD_responses": remote, "small_range_checks": small_range_checks,
        "bounded_8MiB_probes": probes, "partial_markers": partial,
        "E_drive_bytes_before_probes": {"total": disk_before.total, "used": disk_before.used, "free": disk_before.free},
        "E_drive_bytes_after_probes": dict(zip(("total", "used", "free"), shutil.disk_usage(DIRECTORY))),
        "revision_limit": "Official CIIRC static weight URLs provide ETag and Last-Modified, not an immutable HF revision, x-repo-commit, X-Linked-ETag, or LFS SHA256. Preserve this absence; do not invent a revision/hash. Conditional If-Range binds each probe to the exact observed ETag.",
        "recovery_plan": {
            "3d": "Retain original .bin and .bin.aria2 together. Recheck original URL ETag/length before resuming with aria2's control file. Never curl -C - based on 4,242,281,755 logical bytes; that would skip noncontiguous missing pieces.",
            "seq": "Current 12,013,568-byte file has no aria2 map. Its beginning and contiguous-prefix tail match remote bounded probes. Prefer restart into a distinct aria2-managed .part file for maximal auditability; curl -C - is only appropriate after validating the entire local prefix to the same entity, not merely sampling.",
            "completion": "Require exit 0, exact remote length and unchanged ETag, no aria2 remaining-control marker, full SHA256, and ZIP CRC/inventory validation before model load. Remote ETag is not a SHA256; no official cryptographic whole-file digest is provided by this endpoint.",
        },
    }
    output = ROOT / "references/flexpert_weight_recovery_audit_20261024.json"
    output.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({"audit": str(output), "local_files": local[:3], "probes": {key: {
        "complete": value["complete_exact_range_and_ETag_verified"], "bytes": value["actual_probe_bytes"],
        "wall_seconds": value["wall_seconds"], "speed_download_bytes_per_second": value["curl_metrics"].get("speed_download"),
        "curl_seconds": value["curl_metrics"].get("time_total")} for key, value in probes.items()},
        "E_drive_free_bytes": shutil.disk_usage(DIRECTORY).free}, indent=2), flush=True)


if __name__ == "__main__":
    main()
