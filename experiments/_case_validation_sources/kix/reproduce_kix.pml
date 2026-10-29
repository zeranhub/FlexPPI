# Open with PyMOL: File > Run Script, or @absolute_path/reproduce_kix.pml
python
from pathlib import Path
kix_script_dir = Path(__script__).resolve().parent
kix_helper = kix_script_dir / "reproduce_kix.py"
kix_namespace = {"__file__": str(kix_helper), "__name__": "kix_reproduction", "_KIX_BUILD_ONLY": True}
exec(compile(kix_helper.read_text(encoding="utf-8"), str(kix_helper), "exec"), kix_namespace)
kix_namespace["build_kix"](cmd)
python end
