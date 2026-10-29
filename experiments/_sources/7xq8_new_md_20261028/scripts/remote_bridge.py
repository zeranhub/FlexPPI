"""Session-only SSH bridge. The first stdin line is a password, never logged/saved.

Subsequent lines are JSON operations. Host keys must already be trusted by OpenSSH.
"""
import json
import sys
import getpass
from pathlib import Path
import paramiko

sys.stdout.reconfigure(encoding="utf-8")
client = paramiko.SSHClient()
client.load_host_keys(str(Path.home() / ".ssh" / "known_hosts"))
client.set_missing_host_key_policy(paramiko.RejectPolicy())
print("AUTH_INPUT_READY", flush=True)
password = getpass.getpass("SSH password (session memory only): ")
try:
    client.connect("hpc.hrbmu.edu.cn", port=30352, username="root",
                   password=password, allow_agent=False, look_for_keys=False,
                   timeout=15, auth_timeout=20, banner_timeout=20)
finally:
    password = None
print("SSH_AUTHENTICATED", flush=True)
client.get_transport().set_keepalive(30)
try:
    for line in sys.stdin:
        req = json.loads(line)
        try:
            if req["op"] == "exec":
                stdin, stdout, stderr = client.exec_command(req["cmd"], timeout=req.get("timeout", 45))
                stdin.close()
                result = {"stdout": stdout.read().decode("utf-8", "replace"),
                          "stderr": stderr.read().decode("utf-8", "replace"),
                          "exit_code": stdout.channel.recv_exit_status()}
            elif req["op"] in ("upload", "download"):
                with client.open_sftp() as sftp:
                    if req["op"] == "upload":
                        sftp.put(req["local"], req["remote"])
                    else:
                        sftp.get(req["remote"], req["local"])
                result = {"transfer": "complete"}
            elif req["op"] == "close":
                break
            else:
                raise ValueError("Unknown operation")
            print(json.dumps(result, ensure_ascii=False), flush=True)
            if "save" in req:
                Path(req["save"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        except Exception as exc:
            print(json.dumps({"error_type": type(exc).__name__, "message": str(exc)}), flush=True)
finally:
    client.close()
