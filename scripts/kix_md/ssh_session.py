"""Small interactive SSH/SFTP client. Credentials are read without echo, never saved.

Host keys must already be verified in the user's OpenSSH known_hosts file.
After authentication, accept one JSON request per line: exec, put, get, or close.
"""
import argparse
import getpass
import json
import sys
from pathlib import Path

import paramiko


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', required=True)
    parser.add_argument('--port', type=int, default=22)
    parser.add_argument('--user', required=True)
    args = parser.parse_args()
    client = paramiko.SSHClient()
    client.load_host_keys(str(Path.home() / '.ssh' / 'known_hosts'))
    client.set_missing_host_key_policy(paramiko.RejectPolicy())
    password = getpass.getpass('SSH password (not saved): ')
    client.connect(args.host, port=args.port, username=args.user,
                   password=password, look_for_keys=False, allow_agent=False,
                   timeout=20, auth_timeout=30)
    password = None
    client.get_transport().set_keepalive(30)
    sftp = client.open_sftp()
    print(json.dumps({'status': 'connected', 'host': args.host}), flush=True)
    for line in sys.stdin:
        try:
            request = json.loads(line)
            op = request['op']
            if op == 'close':
                break
            if op == 'exec':
                _, stdout, stderr = client.exec_command(request['command'], timeout=request.get('timeout', 120))
                output = stdout.read().decode('utf-8', errors='replace')
                error = stderr.read().decode('utf-8', errors='replace')
                result = {'exit_code': stdout.channel.recv_exit_status(), 'stdout': output, 'stderr': error}
            elif op == 'put':
                sftp.put(request['local'], request['remote'])
                result = {'uploaded': request['remote'], 'bytes': sftp.stat(request['remote']).st_size}
            elif op == 'get':
                sftp.get(request['remote'], request['local'])
                result = {'downloaded': request['local']}
            else:
                raise ValueError('Unsupported operation')
            print(json.dumps(result, ensure_ascii=False), flush=True)
        except Exception as exc:
            print(json.dumps({'error_type': type(exc).__name__, 'error': str(exc)}), flush=True)
    sftp.close()
    client.close()


if __name__ == '__main__':
    main()
