"""Read-only remote collection of COMPLETE frozen neural exports, never inference.

--preflight saves a new local readiness report and exact queue snapshots only.
Without it, both successful queues and all three complete exports are mandatory.
Remote Python reads files and streams a gzip tar to stdout; it writes no files.
CRC, safe paths, exact names/sizes, file SHA256, and an unchanged second remote
snapshot are required before an atomic new local mirror is published. Existing
different mirrors, links, special files, and partial exports are never accepted.
The separate audit_corrected_neural_exports.py must still verify model results.
"""
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
import argparse
import base64
import gzip
import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import tarfile
import uuid
import zipfile

BASE = Path(__file__).resolve().parents[1]
PIPE = BASE
HOST = "hit-gpu"
REMOTE = "/root/flexppi_revision_20261002/corrected_training_v1"
DIRECTORIES = ("evaluation", "external_evaluation", "strict_external_evaluation")
STATES = ("queue_state.json", "strict_queue_state.json")
FREEZES = {
    "training_freeze.json": (PIPE / "training_payload/training_freeze.json", "5b994bb59684e3a33c97f4a1801ee9a606e5dd856db9cfe23bf18515697707d3"),
    "strict_external_freeze.json": (PIPE / "strict_external_payload/strict_external_freeze.json", "16be2a49c6cd2b295a8f150de5e0354d3f6e20028de05f530ae8908c26dbacea"),
    "external/external_payload_freeze.json": (PIPE / "external_payload/external/external_payload_freeze.json", "4dc8c3d560ac7e11eb2bef1be5bab431c58be32281f2712e4c8019715fe5781b"),
}
READERS = {
    "evaluate_corrected_models.py": (PIPE / "training_payload/evaluate_corrected_models.py", "59be9322393daae296610d5e684cc6d62ea9d815b2c88a3f59db260cd4698fa7"),
    "evaluate_external_models.py": (PIPE / "external_payload/evaluate_external_models.py", "d809d798b2b7339d3aa429208da6f2a6f34748c370d9dfa4e19ee43fd467ded5"),
    "evaluate_strict_external_models.py": (PIPE / "strict_external_payload/evaluate_strict_external_models.py", "87623f0466024e92673a2f0b2d622c6763b6b6079d8196d6888c2985367a1c64"),
    "run_training_queue.py": (PIPE / "training_payload/run_training_queue.py", "9c6784ae117195825c191b973a506fef742a457666809cd0ded66fce5ad7a827"),
    "run_strict_external_queue.py": (PIPE / "strict_external_payload/run_strict_external_queue.py", "36eaacb76f2b895f5596e22d0e9e23e030ef268a6139430755917b86504b86a9"),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def require_unlinked(path):
    """Reject Windows junction/reparse points as well as ordinary symlinks."""
    for candidate in [Path(path)] + list(Path(path).parents):
        if candidate.exists() or candidate.is_symlink():
            info = candidate.lstat()
            require(not candidate.is_symlink() and
                    not getattr(info, "st_file_attributes", 0) & getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400),
                    "Linked/reparse local path refused: " + str(candidate))


def safe_name(name):
    """Canonical POSIX relative names also safe on the Windows destination."""
    require(isinstance(name, str) and bool(name), "Empty/non-string archive path")
    require("\\" not in name and ":" not in name and "\x00" not in name,
            "Nonportable archive path")
    parts = name.split("/")
    require(all(p and p not in (".", "..") and not p.endswith((" ", ".")) for p in parts),
            "Noncanonical archive path")
    forbidden = {"CON", "PRN", "AUX", "NUL", "CLOCK$"}
    forbidden.update("COM" + str(i) for i in range(1, 10))
    forbidden.update("LPT" + str(i) for i in range(1, 10))
    require(all(p.split(".", 1)[0].upper() not in forbidden and
                not any(ord(c) < 32 or c in '<>"|?*' for c in p) for p in parts),
            "Windows-unsafe archive component")
    require(str(PurePosixPath(name)) == name and not PurePosixPath(name).is_absolute(),
            "Archive path is not canonical and relative")
    return name


def expected_sources():
    """Bind the remote readers/source/configs to local frozen evidence."""
    require_unlinked(PIPE)
    sources = {}
    for relative, (path, expected) in {**FREEZES, **READERS}.items():
        require(sha(path) == expected, "Local frozen source changed: " + str(path))
        sources[relative] = expected
    for relative, (path, _) in FREEZES.items():
        freeze = read_json(path)
        local_root = path.parent
        for name, expected in freeze.get("execution_hashes", {}).items():
            name = safe_name(name.replace("\\", "/"))
            require(sha(local_root.joinpath(*name.split("/"))) == expected,
                    "Local frozen execution changed: " + name)
            remote_name = name if not relative.startswith("external/") else "external/" + name
            require(remote_name not in sources or sources[remote_name] == expected,
                    "Conflicting remote source bindings")
            sources[remote_name] = expected
        # Input hashes are read from the already hash-bound freeze, rather than
        # inferred from live remote files. Pending preflight skips reading these
        # large immutable inputs; a COMPLETE snapshot verifies every one.
        for key in ("input_hashes", "data_sha256"):
            hashes = freeze.get(key, {})
            if key == "input_hashes" and relative == "training_freeze.json":
                require(len(hashes) == 6126, "Frozen training input hash count differs")
            if key == "data_sha256" and relative.startswith("external/"):
                require(len(hashes) == 64, "Frozen external input hash count differs")
            for name, expected in hashes.items():
                name = safe_name(name.replace("\\", "/"))
                remote_name = "external/" + name if relative.startswith("external/") else name
                require(re.fullmatch("[0-9a-f]{64}", expected) is not None and
                        (remote_name.startswith("dataset/") or
                         remote_name.startswith(("external/structures/", "external/labels/"))),
                        "Malformed frozen input path/hash")
                require(remote_name not in sources or sources[remote_name] == expected,
                        "Conflicting frozen input bindings")
                sources[remote_name] = expected
    extra = {
        "splits/all_ready.jsonl": PIPE / "training_payload/splits/all_ready.jsonl",
        "external/ready_test.jsonl": PIPE / "external_payload/external/ready_test.jsonl",
        "external/extraction_summary.json": PIPE / "external_payload/external/extraction_summary.json",
        "external/prespecified_protocol.json": PIPE / "external_payload/external/prespecified_protocol.json",
    }
    for name, path in extra.items():
        actual = sha(path)
        require(name not in sources or sources[name] == actual, "Local manifest binding differs")
        sources[name] = actual
    return sources


# This remote program only uses Python's standard library. No project import,
# model load, GPU query, training, subprocess, temp file, or remote mutation.
REMOTE_PROGRAM = r'''
import base64, hashlib, json, math, os, stat, sys, tarfile
from pathlib import Path, PurePosixPath
sys.dont_write_bytecode = True
C = json.loads(base64.b64decode(CONFIG_B64))
ROOT = Path(C["root"])
DIRS = ("evaluation", "external_evaluation", "strict_external_evaluation")
STATES = ("queue_state.json", "strict_queue_state.json")

def demand(condition, message):
    if not condition: raise ValueError(message)

def safe(name):
    demand(isinstance(name, str) and name and "\\" not in name and ":" not in name and "\x00" not in name, "Unsafe relative name")
    parts = name.split("/")
    bad = {"CON", "PRN", "AUX", "NUL", "CLOCK$"} | {"COM"+str(i) for i in range(1,10)} | {"LPT"+str(i) for i in range(1,10)}
    demand(all(p and p not in (".", "..") and not p.endswith((" ", ".")) and p.split(".",1)[0].upper() not in bad and not any(ord(c)<32 or c in '<>"|?*' for c in p) for p in parts), "Nonportable relative name")
    demand(str(PurePosixPath(name)) == name and not PurePosixPath(name).is_absolute(), "Noncanonical relative name")
    return name

def path_for(name):
    safe(name)
    demand(ROOT.is_dir() and not ROOT.is_symlink() and ROOT.resolve() == ROOT, "Remote execution root is not the frozen real directory")
    p = ROOT
    for part in name.split("/"):
        p = p / part
        demand(not p.is_symlink(), "Symlink refused: " + name)
    resolved = p.resolve()
    # System python3 on hit-gpu predates pathlib.is_relative_to (Python 3.9).
    demand(resolved == ROOT or ROOT in resolved.parents, "Source escaped root")
    return p

def signature(st):
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_mode)

def file_info(name, include_text=False):
    p = path_for(name)
    before = p.lstat()
    demand(stat.S_ISREG(before.st_mode), "Nonregular source refused: " + name)
    digest = hashlib.sha256()
    chunks = [] if include_text else None
    with p.open("rb") as handle:
        opened = os.fstat(handle.fileno())
        demand(signature(opened) == signature(before), "Source changed while opening: " + name)
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
            if include_text: chunks.append(block)
        demand(signature(os.fstat(handle.fileno())) == signature(before), "Source changed while reading: " + name)
    demand(signature(path_for(name).lstat()) == signature(before), "Source changed after reading: " + name)
    result = dict(bytes=before.st_size, sha256=digest.hexdigest())
    if include_text: result["text"] = b"".join(chunks).decode("utf-8")
    return result

def source_files(ready):
    result = {}
    for name, expected in sorted(C["sources"].items()):
        if not ready and name.startswith(("dataset/", "external/structures/", "external/labels/")):
            continue
        info = file_info(name)
        demand(info["sha256"] == expected, "Remote frozen source mismatch: " + name)
        result[name] = info
    if ready:
        for arm in ("full", "intrachain"):
            members = [(arm,k) for k in range(5)] + [(arm+"_seed"+str(s),0) for s in (11,23)] + [("strict_"+arm,0)]
            for member, fold in members:
                name = "runs/%s/fold_%d/epoch_040.pt" % (member,fold)
                result[name] = file_info(name)
    return result

def captured_text(name, table):
    demand(name in table, "Text source missing from captured snapshot: " + name)
    info = file_info(name, True)
    demand({k:info[k] for k in ("bytes","sha256")} == table[name],
           "Text source changed since captured snapshot: " + name)
    return info["text"]

def json_file(name, table):
    return json.loads(captured_text(name, table))

def lines(name, table):
    return [json.loads(x) for x in captured_text(name, table).splitlines() if x.strip()]

def walk_error(error):
    raise error

def queue_ready(states):
    p, s = states.get("queue_state.json",{}), states.get("strict_queue_state.json",{})
    if p.get("status") != "complete" or s.get("status") != "complete": return False
    demand(p.get("evaluation_exit_code") == p.get("external_exit_code") == 0 and s.get("evaluation_exit_code") == 0, "Export queue completed with a failed exit")
    expected = {(v,k) for v in ("full","intrachain") for k in range(5)} | {(v+"_seed"+str(seed),0) for v in ("full","intrachain") for seed in (11,23)}
    jobs = p.get("jobs",[])
    demand(len(jobs) == 14 and {(j["variant"],j["fold"]) for j in jobs} == expected, "Incomplete/duplicate primary jobs")
    for j in jobs:
        demand(j["status"] in ("complete","already_complete"), "Primary fit not complete")
        if j["status"] == "complete": demand(j.get("exit_code") == 0, "Primary training exit failed")
    jobs = s.get("jobs",[])
    demand(len(jobs)==2 and {j["variant"] for j in jobs}=={"full","intrachain"} and all(j["status"]=="complete" for j in jobs), "Strict fits incomplete/duplicated")
    demand(all(j.get("exit_code",0)==0 for j in jobs), "Strict training exit failed")
    return True

def export_names():
    names = []
    for directory in DIRS:
        start = path_for(directory)
        demand(start.is_dir(), "Complete export directory missing")
        for here, dirs, files in os.walk(start, followlinks=False, onerror=walk_error):
            for child in dirs:
                p = Path(here)/child
                demand(not p.is_symlink() and p.is_dir(), "Linked or special export directory")
                safe(p.relative_to(ROOT).as_posix())
            for child in files:
                name = (Path(here)/child).relative_to(ROOT).as_posix()
                demand(stat.S_ISREG(path_for(name).lstat().st_mode), "Special export file")
                names.append(safe(name))
    names += list(STATES)
    demand(len(names)==len(set(n.casefold() for n in names)), "Duplicate/case-colliding export filenames")
    return sorted(names)

def export_contract(files, sources):
    cohort = lines("splits/all_ready.jsonl", sources)
    records = {r["sample_id"]:r for r in cohort}
    external = lines("external/ready_test.jsonl", sources)
    ext = {r["sample_id"]:r for r in external}
    demand(len(records)==len(cohort)==3063 and len(ext)==len(external)==32, "Frozen record membership differs")
    expected_npz = set()
    for scope, strict in (("evaluation",False),("external_evaluation",False),("strict_external_evaluation",True)):
        audit = json_file(scope+"/evaluation_audit.json", files)
        code = "evaluate_corrected_models.py" if scope=="evaluation" else "evaluate_strict_external_models.py" if strict else "evaluate_external_models.py"
        freeze = "strict_external_freeze.json" if strict else "training_freeze.json"
        demand(audit.get("status")=="complete" and audit.get("epoch")==40 and audit.get("export_code_sha256")==sources[code]["sha256"], "Incomplete or wrong-reader export audit")
        key = "input_freeze_sha256" if scope=="evaluation" else "training_freeze_sha256"
        demand(audit[key]==sources[freeze]["sha256"], "Export training freeze differs")
        if scope=="evaluation":
            demand(audit["oof_trained_fits"]==10 and audit["additional_seed_fits"]==4 and audit["seed"]==7 and audit["records_per_variant"]==3063, "OOF fit/support audit differs")
            expected_keys = {(v,int(r["metadata"]["fold"]),7,k) for v in ("full","intrachain") for k,r in records.items()}
            required = ("per_record_metrics.csv","oof_per_residue.csv.gz","source_fold_macro_metrics.csv","five_fold_mean_sample_sd.csv")
        else:
            demand(audit["scored_records"]==audit["prespecified_records"]==32 and audit["checkpoints"]==(2 if strict else 14) and not audit["extraction_failures"], "External count/failure audit differs")
            demand(audit["external_payload_sha256"]==sources["external/external_payload_freeze.json"]["sha256"] and audit["prespecified_protocol_sha256"]==sources["external/prespecified_protocol.json"]["sha256"], "External input freeze differs")
            members = [(v,0,7) for v in ("full","intrachain")] if strict else [(v,k,7) for v in ("full","intrachain") for k in range(5)] + [(v,0,s) for v in ("full","intrachain") for s in (11,23)]
            expected_keys = {(v,k,s,x) for v,k,s in members for x in ext}
            required = ("per_record_metrics.csv","macro_metrics.csv","per_residue.csv.gz")
            for arm in ("full","intrachain"):
                expected_npz.update(scope+"/"+arm+"/ensemble/"+sample+".npz" for sample in ext)
        for filename in required+("prediction_manifest.jsonl","evaluation_audit.json"):
            demand(scope+"/"+filename in files, "Complete output table missing")
        observed = set()
        for row in lines(scope+"/prediction_manifest.jsonl", files):
            is_oof = scope=="evaluation"
            needed = {"sample_id","variant","fold","checkpoint_sha256","input_sha256","label_sha256","cache_path"} | ({"record","dataset"} if is_oof else {"seed","forward_seconds","residue_count"})
            demand(needed.issubset(row), "Incomplete prediction manifest row")
            key = (row["variant"],int(row["fold"]),7 if is_oof else int(row["seed"]),row["sample_id"])
            demand(key in expected_keys and key not in observed, "Unexpected/duplicate prediction key")
            observed.add(key)
            arm,fold,seed,sample = key
            member = "strict_"+arm if strict else arm if seed==7 else arm+"_seed"+str(seed)
            relative = scope+"/"+arm+"/"+sample+".npz" if is_oof else scope+"/"+member+"/fold_"+str(fold)+"/"+sample+".npz"
            demand(row["cache_path"]==str(ROOT)+"/"+relative and relative in files, "Manifest cache path/membership differs")
            demand(row["checkpoint_sha256"]==sources["runs/%s/fold_%d/epoch_040.pt"%(member,fold)]["sha256"], "Manifest checkpoint source hash differs")
            r = records[sample] if is_oof else ext[sample]
            if is_oof: demand(row["record"]==r and row["dataset"]==r["dataset"], "OOF source identity differs")
            else: demand(math.isfinite(row["forward_seconds"]) and row["forward_seconds"]>0 and row["residue_count"]==r["metadata"]["residue_count"], "Invalid forward metadata")
            for field in ("structure_path","label_path"):
                demand(r[field].startswith(str(ROOT)+"/"), "Frozen input outside root")
                relative_input = safe(r[field][len(str(ROOT))+1:])
                demand(relative_input in sources and relative_input in C["sources"] and
                       sources[relative_input]["sha256"] == C["sources"][relative_input] ==
                       row["input_sha256" if field=="structure_path" else "label_sha256"],
                       "Prediction input/label differs from frozen input SHA")
            expected_npz.add(relative)
        demand(observed==expected_keys, "Incomplete prediction manifest")
    demand({x for x in files if x.endswith(".npz")}==expected_npz and len(expected_npz)==6766, "Missing/extra prediction NPZ")

def snapshot():
    state_files, states, pending = {}, {}, []
    for name in STATES:
        p=path_for(name)
        if not p.exists():
            pending.append(name+" missing"); continue
        info=file_info(name,True)
        state_files[name]=info
        try: states[name]=json.loads(info["text"])
        except json.JSONDecodeError: pending.append(name+" not yet a complete JSON state")
    ready = not pending and queue_ready(states)
    sources=source_files(ready)
    result=dict(ready=bool(ready), states=states, state_files=state_files, source_files=sources,
                pending_reasons=pending or ([] if ready else ["Both primary and strict queue/export statuses must be complete"]), files={})
    if ready:
        files={name:file_info(name) for name in export_names()}
        export_contract(files,sources)
        for name,info in state_files.items():
            demand(files[name]["sha256"]==info["sha256"], "Queue changed during complete snapshot")
        result["files"]=files
    return result

if C["action"]=="snapshot":
    print(json.dumps(snapshot(),ensure_ascii=True))
elif C["action"]=="archive":
    before=snapshot()
    demand(before==C["snapshot"], "Source changed before archive streaming")
    with tarfile.open(fileobj=sys.stdout.buffer,mode="w|gz",format=tarfile.USTAR_FORMAT) as archive:
        for name,expected in sorted(before["files"].items()):
            p=path_for(name)
            demand(file_info(name)==expected,"Source changed before archive member")
            with p.open("rb") as handle:
                info=archive.gettarinfo(str(p),arcname=name)
                demand(info.isreg() and info.size==expected["bytes"],"Nonregular/changing archive source")
                archive.addfile(info,handle)
            demand(file_info(name)==expected,"Source changed after archive member")
else: raise ValueError("Unknown collector action")
'''


def remote_script(sources, action="snapshot", snapshot=None):
    config = dict(root=REMOTE, sources=sources, action=action)
    if snapshot is not None:
        config["snapshot"] = snapshot
    encoded = base64.b64encode(json.dumps(config, sort_keys=True).encode()).decode("ascii")
    return ("CONFIG_B64 = " + repr(encoded) + "\n" + REMOTE_PROGRAM).encode("utf-8")


def remote_snapshot(sources):
    result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", HOST, "python3 -"],
                            input=remote_script(sources), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=300, check=False)
    require(result.returncode == 0, "Read-only remote snapshot failed: " + result.stderr.decode("utf-8", errors="replace"))
    snapshot = json.loads(result.stdout.decode("utf-8"))
    require(isinstance(snapshot.get("ready"), bool), "Invalid remote readiness response")
    return snapshot


def verify_file_table(table):
    require(len(table) == len({name.casefold() for name in table}), "Case-colliding file names")
    for name, info in table.items():
        safe_name(name)
        require(name in STATES or name.split("/")[0] in DIRECTORIES, "File outside allowed export scope")
        require(isinstance(info["bytes"], int) and info["bytes"] >= 0 and
                re.fullmatch("[0-9a-f]{64}", info["sha256"]) is not None, "Malformed remote file size/SHA")


def read_crc(path):
    # Read to EOF independently of tarfile, which may stop at tar end markers
    # before consuming the gzip trailer. GzipFile verifies its CRC/size there.
    with gzip.open(path, "rb") as stream:
        for _ in iter(lambda: stream.read(1024 * 1024), b""):
            pass


def verify_tar_stream(path, table):
    """Validate every raw header and terminal padding, including gzip CRC."""
    seen = set()
    with gzip.open(path, "rb") as stream:
        while True:
            header = stream.read(512)
            require(len(header) == 512, "Truncated tar header/end markers")
            if header == bytes(512):
                require(stream.read(512) == bytes(512), "Missing second tar end marker")
                # TarFile may silently stop on zero blocks or a bad trailing
                # header; only zero padding is allowed after the explicit EOF.
                for padding in iter(lambda: stream.read(1024 * 1024), b""):
                    require(not any(padding), "Unexpected data after tar end markers")
                break
            member = tarfile.TarInfo.frombuf(header, encoding="utf-8", errors="strict")
            name = safe_name(member.name)
            require(member.type in (tarfile.REGTYPE, tarfile.AREGTYPE) and not member.linkname,
                    "Archive header is not an ordinary regular file")
            require(name in table and name not in seen and member.size == table[name]["bytes"],
                    "Raw archive header names/sizes differ")
            seen.add(name)
            remaining = ((member.size + 511) // 512) * 512
            while remaining:
                block = stream.read(min(remaining, 1024 * 1024))
                require(bool(block), "Truncated tar member content")
                remaining -= len(block)
    require(seen == set(table), "Raw tar stream omitted expected files")


def inspect_and_extract(archive_path, table, staging):
    verify_file_table(table)
    require_unlinked(staging)
    verify_tar_stream(archive_path, table)
    seen = set()
    # No extractall: each validated regular member is streamed into a fresh
    # exclusive file beneath the owned staging directory.
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            name = safe_name(member.name)
            require(member.type in (tarfile.REGTYPE, tarfile.AREGTYPE) and member.sparse is None and
                    not member.linkname and not member.pax_headers,
                    "Links/special/sparse/PAX archive members are refused")
            require(name in table and name not in seen and member.size == table[name]["bytes"],
                    "Unexpected/duplicate/wrong-size archive member")
            seen.add(name)
            target = staging.joinpath(*name.split("/"))
            require(target.resolve().is_relative_to(staging.resolve()), "Archive destination escaped staging")
            target.parent.mkdir(parents=True, exist_ok=True)
            require_unlinked(target.parent)
            digest, count = hashlib.sha256(), 0
            source = archive.extractfile(member)
            require(source is not None, "Missing regular archive member stream")
            with source, target.open("xb") as output:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    count += len(block)
                    digest.update(block)
                    output.write(block)
            require(count == member.size and digest.hexdigest() == table[name]["sha256"],
                    "Transferred archive member differs from remote SHA")
            if name.endswith(".npz"):
                with zipfile.ZipFile(target) as zipped:
                    require(zipped.testzip() is None, "NPZ ZIP CRC mismatch")
            elif name.endswith(".gz"):
                read_crc(target)
    require(seen == set(table), "Missing archive members")
    return dict(archive_gzip_crc=True, tar_header_checksums=True, exact_members=len(seen),
                npz_zip_crc=True, compressed_csv_crc=True, source_sha256=True)


def mirror_table(root):
    require_unlinked(root)
    require(root.is_dir() and not root.is_symlink(), "Existing mirror is not a regular directory")
    result = {}
    def walk_error(error):
        raise error
    for here, dirs, files in os.walk(root, followlinks=False, onerror=walk_error):
        for child in dirs:
            p = Path(here) / child
            require_unlinked(p)
            require(not p.is_symlink() and p.is_dir(), "Existing mirror contains a directory link")
        for child in files:
            p = Path(here) / child
            require_unlinked(p)
            require(not p.is_symlink() and stat.S_ISREG(p.lstat().st_mode), "Existing mirror contains linked/special file")
            name = safe_name(p.relative_to(root).as_posix())
            result[name] = dict(bytes=p.stat().st_size, sha256=sha(p))
    require(len(result) == len({name.casefold() for name in result}), "Mirror has case aliases")
    return result


def save_snapshot(directory, snapshot, prefix):
    (directory / (prefix + ".json")).write_text(json.dumps(snapshot, indent=2, ensure_ascii=True), encoding="utf-8")
    for name, info in snapshot["state_files"].items():
        path = directory / (prefix + "_" + name + ".snapshot")
        path.write_bytes(info["text"].encode("utf-8"))
        require(sha(path) == info["sha256"], "Local exact queue snapshot changed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true", help="Read readiness and save QA only; never transfer outputs")
    parser.add_argument("--mirror-path", type=Path, default=PIPE / "neural_exports")
    args = parser.parse_args()
    tag = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    report_dir = BASE / "qa/neural_export_collection" / tag
    report_dir.mkdir(parents=True, exist_ok=False)
    report = dict(status="started; no collection success", checked_utc=datetime.now(timezone.utc).isoformat(),
                  collector_sha256=sha(__file__), host=HOST, remote_root=REMOTE,
                  preflight_only=args.preflight, new_fits=0, new_forward_calls=0,
                  remote_files_modified=0, predictions_transferred=False,
                  final_prediction_audit_complete=False, accuracy_evaluation_complete=False)
    exit_code = 0
    try:
        sources = expected_sources()
        before = remote_snapshot(sources)
        save_snapshot(report_dir, before, "before")
        report.update(source_exports_complete=before["ready"], queue_statuses={n: s.get("status") for n, s in before["states"].items()},
                      local_frozen_source_bindings=sources, remote_source_sha256=before["source_files"],
                      pending_reasons=before["pending_reasons"])
        if not before["ready"]:
            report["status"] = "pending; incomplete queues/exports; no predictions collected"
            exit_code = 0 if args.preflight else 2
        elif args.preflight:
            report["status"] = "ready; preflight only; outputs not transferred or scientifically audited"
        else:
            # The no-replacement rename guarantee used below is native Windows
            # semantics; POSIX rename may replace a concurrently created empty
            # directory. Refuse formal publication outside this supported OS.
            require(os.name == "nt", "Formal atomic publication currently requires native Windows")
            table = before["files"]
            verify_file_table(table)
            target = args.mirror_path.absolute()
            require_unlinked(target)
            require(target.parent.resolve() == PIPE.resolve() and
                    (target.name == "neural_exports" or target.name.startswith("neural_exports_")),
                    "Mirror must be a dedicated neural_exports directory beneath the pipeline")
            safe_name(target.name)
            require(not target.is_symlink() and target.parent.resolve() == target.parent,
                    "Linked mirror or mirror parent is refused")
            if target.exists():
                require(mirror_table(target) == table, "Existing mirror differs; overwrite refused")
                after = remote_snapshot(sources)
                save_snapshot(report_dir, after, "after")
                require(after == before, "Remote queue/source/export changed before mirror reuse")
                report.update(status="existing identical mirror verified; no overwrite or transfer", mirror=str(target),
                              predictions_transferred=False, collection_verified=True)
            else:
                archive_path = target.parent / (target.name + "." + uuid.uuid4().hex + ".tar.gz.partial")
                staging = target.parent / (target.name + "." + uuid.uuid4().hex + ".partial")
                staging.mkdir(exist_ok=False)
                report.update(archive_partial=str(archive_path), staging_partial=str(staging))
                # Binary stdout directly to an exclusive local file: never
                # PowerShell redirection or text decoding of the tar stream.
                with archive_path.open("xb") as output:
                    result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", HOST, "python3 -"],
                                            input=remote_script(sources, "archive", before), stdout=output,
                                            stderr=subprocess.PIPE, timeout=7200, check=False)
                require(result.returncode == 0, "Remote archive transfer failed; partial retained: " + result.stderr.decode("utf-8", errors="replace"))
                report["archive_sha256"] = sha(archive_path)
                report["transfer_checks"] = inspect_and_extract(archive_path, table, staging)
                require(mirror_table(staging) == table, "Extracted mirror differs from source manifest")
                after = remote_snapshot(sources)
                save_snapshot(report_dir, after, "after")
                require(after == before, "Remote queues/sources/exports changed during transfer; partial retained")
                if target.exists():
                    require(mirror_table(target) == table, "Concurrent destination differs; refused overwrite")
                    report["identical_concurrent_destination_reused"] = True
                else:
                    # Same-parent directory rename is atomic. No replace,
                    # deletion, or merge of any prior mirror is permitted.
                    staging.rename(target)
                require(mirror_table(target) == table, "Published mirror no longer matches verified source")
                report.update(status="collected_verified; separate scientific output audit still required",
                              mirror=str(target), predictions_transferred=True, collection_verified=True,
                              collected_files_sha256=table)
        report["completed_utc"] = datetime.now(timezone.utc).isoformat()
    except (Exception, KeyboardInterrupt) as error:
        report.update(status="failed or interrupted; no collection success", collection_verified=False,
                      final_prediction_audit_complete=False, error_type=type(error).__name__, error=str(error),
                      completed_utc=datetime.now(timezone.utc).isoformat(), partial_policy="All owned partial files retained; never treated as final outputs")
        exit_code = 1
    finally:
        path = report_dir / "collection_audit.json"
        path.write_text(json.dumps(report, indent=2, ensure_ascii=True), encoding="utf-8")
        print(report["status"], flush=True)
        print("Read-only collection QA:", path, flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
