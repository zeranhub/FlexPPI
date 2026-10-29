# Open with PyMOL: File > Run Script, or @absolute_path/rebuild_relative.pml.
# Exact complete-protein scenes are rebuilt without repeating ray tracing.
# To regenerate PNGs, run render_from_bundle.py without --build-only.
python
from pathlib import Path
import runpy
_bundle = Path(__script__).resolve().parent
_render_ns = runpy.run_path(str(_bundle / "render_from_bundle.py"))
_render_ns["render"](_bundle, build_only=True, command=cmd)
python end
