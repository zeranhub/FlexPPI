"""Dispatch the manuscript's existing programs without importing scientific packages."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

EXPERIMENTS = Path(__file__).resolve().parent
PROJECT = EXPERIMENTS.parent


def main(section: str) -> int:
    registry = json.loads((EXPERIMENTS / "steps.json").read_text(encoding="utf-8"))
    record = registry[section]
    steps = record["steps"]
    parser = argparse.ArgumentParser(
        description=section.replace("_", ".") + " " + record["title"],
        epilog="先用 --list 或 --check 查看。执行必须显式 --step；原程序参数放在 -- 之后，例如 --step NAME -- --help。",
    )
    parser.add_argument("--list", action="store_true", help="列出实际实现、输入和用途，不运行实验")
    parser.add_argument("--check", action="store_true", help="只检查脚本及登记的本地输入是否存在")
    parser.add_argument("--step", choices=list(steps), help="明确选择一个原程序运行")
    parser.add_argument("--python", default=sys.executable, help="运行原程序的Python解释器，默认当前解释器")
    arguments = sys.argv[1:]
    if "--" in arguments:
        boundary = arguments.index("--")
        own, forwarded = arguments[:boundary], arguments[boundary + 1:]
        options = parser.parse_args(own)
    else:
        options, forwarded = parser.parse_known_args(arguments)
    if (options.list or options.check) and forwarded:
        parser.error("--list/--check 不接受原程序参数；查看原帮助请使用 --step NAME -- --help")
    if options.check:
        failed = False
        for name, item in steps.items():
            checks = [item["script"], *item.get("inputs", [])]
            missing = [path for path in checks if not (EXPERIMENTS / path).exists()]
            failed |= bool(missing)
            print(("MISSING " if missing else "OK ") + name)
            for path in missing:
                print("  " + path)
        print("只检查登记路径；不验证外部轨迹、科学环境、冻结哈希或GPU可用性。")
        return 1 if failed else 0
    if options.list or not options.step:
        if forwarded:
            parser.error("运行原程序参数前请指定 --step")
        print(parser.description)
        for name, item in steps.items():
            print("\n" + name + ": " + item["description"])
            print("  " + str((EXPERIMENTS / item["script"]).resolve()))
            if item.get("notes"):
                print("  " + item["notes"])
        print("\n未执行实验。使用 --step NAME -- 原程序参数 执行一个步骤。")
        return 0
    item = steps[options.step]
    script = (EXPERIMENTS / item["script"]).resolve()
    if not script.is_file():
        parser.error("实现文件不存在：" + str(script))
    environment = os.environ.copy()
    imports = [str(EXPERIMENTS / "core_corrected_pipeline" / "training_payload" / "src"),
               str(PROJECT / "scripts"), str(PROJECT / "src")]
    if environment.get("PYTHONPATH"):
        imports.append(environment["PYTHONPATH"])
    environment["PYTHONPATH"] = os.pathsep.join(imports)
    # Source scripts retain their own __file__ and their original CLI semantics.
    # No scientific module is imported and no result is changed by this dispatcher.
    cwd = script.parent if item.get("cwd") == "script" else PROJECT
    return subprocess.run([options.python, str(script), *item.get("arguments", []), *forwarded], cwd=str(cwd),
                          env=environment).returncode
