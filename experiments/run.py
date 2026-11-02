"""Explicitly dispatch workflows retained in manuscript v136."""
import json
from pathlib import Path
import sys
from _dispatch import main

if __name__ == '__main__':
    registry = json.loads((Path(__file__).parent / 'steps.json').read_text(encoding='utf-8'))
    if len(sys.argv) < 2 or sys.argv[1] in ['--help', '-h']:
        print('Usage: python experiments/run.py WORKFLOW [--list | --check | --step NAME -- ARGUMENTS]')
        for name, record in registry.items():
            print(name + ': ' + record['title'])
        raise SystemExit(0)
    workflow = sys.argv.pop(1)
    if workflow not in registry:
        raise SystemExit('Unknown workflow: ' + workflow)
    raise SystemExit(main(workflow))
