"""Read-only SSH inventory; password is prompted and never saved or printed."""
import getpass
import json
import shlex
from datetime import datetime, timezone
from pathlib import Path

import paramiko

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'output/tcr_external_validation_search_20260917/remote_capacity.json'


def main():
    password = getpass.getpass('SSH password (not saved): ')
    client = paramiko.SSHClient()
    client.load_host_keys(str(Path.home() / '.ssh/known_hosts'))
    client.set_missing_host_key_policy(paramiko.RejectPolicy())
    client.connect('hpc.hrbmu.edu.cn', port=30352, username='root', password=password,
                   allow_agent=False, look_for_keys=False, timeout=15,
                   auth_timeout=20, banner_timeout=20)
    del password
    checks = {
        'host': 'hostname; uname -a',
        'gpu': 'nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu --format=csv',
        'gpu_processes': 'nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory --format=csv',
        'cpu_memory': 'lscpu | head -22; free -h',
        'disk': 'df -h / /root /tmp; ls -ld /data /workspace /root/autodl-tmp /root/FlexPPI /root/miniconda3 /opt/conda 2>/dev/null',
        'software': 'for x in gmx gmx_mpi sbatch sinfo squeue conda python python3 nvcc; do command -v "$x" || true; done',
        'gromacs': 'if command -v gmx >/dev/null 2>&1; then gmx --version; elif command -v gmx_mpi >/dev/null 2>&1; then gmx_mpi --version; fi',
        'scheduler': 'if command -v sinfo >/dev/null 2>&1; then sinfo -h -o "%P %a %l %D %G"; fi',
        'software_locations': 'ls -d /usr/local/gromacs* /opt/gromacs* /root/*/bin/gmx* /root/miniconda3/envs/*/bin/gmx* /opt/conda/envs/*/bin/gmx* 2>/dev/null',
    }
    report = {'retrieved_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Read-only inventory, no installations or simulation jobs launched', 'checks': {}}
    try:
        for name, command in checks.items():
            _, stdout, stderr = client.exec_command('bash -lc ' + shlex.quote(command), timeout=35)
            output, error = stdout.read().decode('utf-8', 'replace'), stderr.read().decode('utf-8', 'replace')
            report['checks'][name] = {'command': command, 'stdout': output, 'stderr': error,
                                     'exit_status': stdout.channel.recv_exit_status()}
            print(name + ':\n' + output + error, flush=True)
    finally:
        client.close()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print('Saved read-only inventory: ' + str(OUT), flush=True)


if __name__ == '__main__':
    main()
