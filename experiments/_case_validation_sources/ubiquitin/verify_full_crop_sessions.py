"""Reopen the delivered PyMOL sessions and verify the full-protein crop chain."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np
from PIL import Image
import pymol2

BUNDLE = Path(__file__).resolve().parent


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    manifest = json.loads((BUNDLE / "crop_manifest.json").read_text(encoding="utf-8"))
    views = json.loads((BUNDLE / "views.json").read_text(encoding="utf-8"))
    assert len(views) == 5
    report = dict(status="pass", cases={}, sessions={})
    for name, plan in manifest["cases"].items():
        full_path = BUNDLE / "renders" / plan["full_image"]
        crop_path = BUNDLE / "renders" / plan["crop_image"]
        with Image.open(full_path) as original, Image.open(crop_path) as crop:
            assert list(original.size) == plan["full_size"]
            assert list(crop.size) == plan["crop_size"]
            assert original.crop(tuple(plan["crop_box"])).tobytes() == crop.tobytes()
        assert sha256(full_path) == plan["full_image_sha256"]
        assert sha256(crop_path) == plan["crop_image_sha256"]
        assert plan["full_view"][15:17] == plan["focus_view"][15:17]
        assert plan["full_view"][11] == plan["focus_view"][11]
        report["cases"][name] = dict(full_image=full_path.name, crop_image=crop_path.name,
            crop_box=plan["crop_box"], exact_pixel_crop=True,
            full_and_focus_same_distance_and_depth=True)
    for relative, expected in manifest["input_PDB_sha256_before"].items():
        assert sha256(BUNDLE / relative) == expected
    with pymol2.PyMOL() as pm:
        cmd = pm.cmd
        for session_name, initial_scene in [
                ("ubiquitin_2QHO.pse", "01_complex"),
                ("f_Q41_local.pse", "02_Q41_MD_observed"),
                ("f_Q62_local.pse", "02_Q62_MD_observed")]:
            cmd.reinitialize()
            cmd.load(str(BUNDLE / session_name))
            assert set(cmd.get_scene_list()) == set(views)
            assert np.max(np.abs(np.asarray(cmd.get_view())-views[initial_scene])) < 1e-4
            complex_coordinates = np.asarray(cmd.get_coords("complex")).copy()
            assert complex_coordinates.shape == (935, 3)
            verified_scenes = {}
            for scene, expected_view in views.items():
                cmd.scene(scene, "recall", animate=0)
                actual_view = list(cmd.get_view())
                difference = float(np.max(np.abs(np.asarray(actual_view)-expected_view)))
                assert difference < 1e-4, (session_name, scene, difference)
                ca_count = {chain: cmd.count_atoms(
                    f"complex and chain {chain} and name CA and rep cartoon")
                            for chain in ["A", "B"]}
                assert ca_count == {"A": 72, "B": 47}, (scene, ca_count)
                assert "complex" in cmd.get_names("objects", enabled_only=1)
                assert np.array_equal(np.asarray(cmd.get_coords("complex")), complex_coordinates)
                if scene.startswith("03_"):
                    site = "Q41" if "Q41" in scene else "Q62"
                    plan = manifest["cases"][site]
                    full_view = actual_view
                    rotation = np.asarray(full_view[:9]).reshape(3, 3)
                    rotated = (complex_coordinates-np.asarray(full_view[12:15])) @ rotation
                    w, h = plan["full_size"]
                    field_h = 2*abs(full_view[11])*math.tan(math.radians(full_view[17]/2))
                    field_w = field_h*w/h
                    x_pixels = w/2+(rotated[:,0]+full_view[9])*w/field_w
                    y_pixels = h/2-(rotated[:,1]+full_view[10])*h/field_h
                    assert x_pixels.min() > 0 and x_pixels.max() < w
                    assert y_pixels.min() > 0 and y_pixels.max() < h
                    depth = -full_view[11]-rotated[:,2]
                    assert depth.min() > full_view[15] and depth.max() < full_view[16]
                verified_scenes[scene] = dict(camera_max_difference=difference,
                    CA_cartoon_by_chain=ca_count, CA_cartoon_total=sum(ca_count.values()))
            cmd.load(str(BUNDLE / "inputs/A003R_input_structure.pdb"), "verification_source", zoom=0)
            source_coordinates = np.asarray(cmd.get_coords("verification_source and not hydro"))
            assert source_coordinates.shape == complex_coordinates.shape
            source_difference = float(np.max(np.abs(source_coordinates-complex_coordinates)))
            assert source_difference < 1e-5
            report["sessions"][session_name] = dict(default_scene=initial_scene,
                source_coordinate_max_difference_A=source_difference,
                scenes=verified_scenes, saved_scene_count=len(views))
    (BUNDLE / "full_crop_session_verification.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    verify()
