"""Launch the two authorized independent PLANET CPU queues under nohup."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

ROOT = Path('/root/flexppi_completion_20261010')
BASE = ROOT / 'baselines'
PYTHON = ROOT / 'planet_env' / 'bin' / 'python'
WRAPPER = BASE / 'run_planet_cpu_queue.py'
EXPECTED_WRAPPER_SHA256 = '034a7fdfc5b6bc39164bb43f23518d09586863af186ece495275965c3910e602'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    assert sha(WRAPPER) == EXPECTED_WRAPPER_SHA256
    entries = []
    for variant in ['full', 'seq']:
        preflight = BASE / ('planet_' + variant + '_cpu_queue_preflight.json')
        audit = json.loads(preflight.read_text())
        assert audit['status'] == 'passed' and audit['chains_verified'] == 4973
        assert audit['wrapper_sha256'] == EXPECTED_WRAPPER_SHA256
        log = BASE / ('planet_' + variant + '_cpu_queue.log')
        command = ['nohup', 'env', 'CUDA_VISIBLE_DEVICES=', 'OMP_NUM_THREADS=2',
                   'MKL_NUM_THREADS=2', 'OPENBLAS_NUM_THREADS=2', 'PYTHONDONTWRITEBYTECODE=1',
                   str(PYTHON), '-u', str(WRAPPER), '--variant', variant]
        with log.open('ab', buffering=0) as output:
            process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=output,
                                       stderr=subprocess.STDOUT, start_new_session=True, close_fds=True)
        entries.append(dict(variant=variant, queue_parent_pid=process.pid, launcher_pid=os.getpid(),
            command=command, command_shell=shlex.join(command), log=str(log),
            state=str(BASE / ('planet_' + variant + '_cpu_queue_state.json')),
            preflight_sha256=sha(preflight), preflight_checked_utc=audit['checked_utc']))
    record = dict(status='dispatched', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  wrapper_sha256=EXPECTED_WRAPPER_SHA256, launcher_sha256=sha(__file__), entries=entries)
    target = BASE / 'planet_cpu_queue_launch_20261010.json'
    target.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
