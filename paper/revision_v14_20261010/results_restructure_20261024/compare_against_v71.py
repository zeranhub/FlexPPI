"""Compare a final v72 canonical export with fully inspected v71 r3."""
from pathlib import Path
import subprocess
import sys

base = Path(__file__).resolve().parent
reference = base.parent / "table1_to_supplement_20261024" / "render_v71_final_r3"
current = Path(sys.argv[1]) if len(sys.argv) > 1 else base / "render_v72_final"
subprocess.run([sys.executable, str(base / "compare_render_runs.py"), str(reference), str(current)], check=True)
