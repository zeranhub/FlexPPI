"""Run the three portable PML entrypoints from an external working directory."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import numpy as np
from pymol import cmd

BUNDLE = Path(__file__).resolve().parent


def verify_entry(scene, views):
    assert set(cmd.get_scene_list()) == set(views)
    ca = {chain: cmd.count_atoms(f"complex and chain {chain} and name CA and rep cartoon")
          for chain in ["A", "B"]}
    assert ca == {"A": 72, "B": 47}
    difference = float(np.max(np.abs(np.asarray(cmd.get_view())-views[scene])))
    assert difference < 1e-4, (scene, difference)
    return dict(expected_default_scene=scene, saved_scene_count=len(views),
        complete_complex_CA_cartoon=ca, camera_max_difference=difference,
        enabled_objects=cmd.get_names("objects", enabled_only=1))


def main():
    assert Path.cwd().resolve() != BUNDLE.resolve()
    # Use PyMOL's own sequential PML parser. Enqueuing a Python block with
    # cmd.do from a separate Python main thread can race the following checks.
    folder = BUNDLE.as_posix()
    wrapper = f'''@{folder}/rebuild_relative.pml
python
from pathlib import Path
import json
import runpy
_vbundle=Path(r"{folder}")
_smoke=runpy.run_path(str(_vbundle/"smoke_test_portable_entrypoints.py"))
_views=json.loads((_vbundle/"views.json").read_text(encoding="utf-8"))
_report=dict(status="pass",executed_from=str(Path.cwd()),tested_bundle=str(_vbundle),entrypoints={{}})
_report["entrypoints"]["rebuild_relative.pml"]=_smoke["verify_entry"]("01_complex",_views)
python end
@{folder}/f_Q41_local.pml
python
_report["entrypoints"]["f_Q41_local.pml"]=_smoke["verify_entry"]("02_Q41_MD_observed",_views)
python end
@{folder}/f_Q62_local.pml
python
_report["entrypoints"]["f_Q62_local.pml"]=_smoke["verify_entry"]("02_Q62_MD_observed",_views)
(_vbundle/"portable_entrypoint_smoke_test.json").write_text(json.dumps(_report,indent=2)+"\\n",encoding="utf-8")
python end
quit
'''
    report_path = BUNDLE / "portable_entrypoint_smoke_test.json"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pml", encoding="utf-8",
                                     delete=False) as handle:
        handle.write(wrapper)
        wrapper_path = Path(handle.name)
    try:
        report_path.unlink(missing_ok=True)
        subprocess.run([sys.executable, "-m", "pymol", "-cq", str(wrapper_path)],
                       cwd=Path.cwd(), check=True, timeout=60)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        assert report["status"] == "pass"
    finally:
        wrapper_path.unlink(missing_ok=True)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
