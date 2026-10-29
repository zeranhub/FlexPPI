"""Render a complete protein first, then obtain a local image by exact pixel crop.

The saved focus camera keeps the complete protein and its full depth interval.
Only its field of view and image-plane translation change. Coordinates, objects,
and cartoon representations are never removed to make a local panel.
"""
from pathlib import Path
import hashlib
import math
from PIL import Image, ImageDraw


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plan_crop(c, full_selection, focus_view, full_size=(6000, 4800),
              crop_aspect=(5, 4), margin_A=4.0):
    view = list(focus_view)
    assert view[17] > 0, "An orthoscopic focus camera is required"
    points = c.get_coords(full_selection)
    assert points is not None and len(points)
    # get_view rotation is column-major: world coordinate i contributes to
    # screen coordinate j through view[3*i+j].
    rotated = [[sum((float(p[i])-view[12+i])*view[3*i+j]
                    for i in range(3)) for j in range(3)] for p in points]
    lo = [min(p[j] for p in rotated) for j in range(3)]
    hi = [max(p[j] for p in rotated) for j in range(3)]
    desired_h = 2*abs(view[11])*math.tan(math.radians(view[17]/2))
    desired_w = desired_h*crop_aspect[0]/crop_aspect[1]
    # Include the entire desired crop as well as all atoms in the full frame.
    target_x, target_y = -view[9], -view[10]
    left = min(lo[0]-margin_A, target_x-desired_w/2-margin_A)
    right = max(hi[0]+margin_A, target_x+desired_w/2+margin_A)
    bottom = min(lo[1]-margin_A, target_y-desired_h/2-margin_A)
    top = max(hi[1]+margin_A, target_y+desired_h/2+margin_A)
    W, H = map(int, full_size)
    full_h = max(top-bottom, (right-left)*H/W)
    # Keep the complete object in front of the camera, including its depth.
    full_h = max(full_h, 2*(hi[2]+margin_A+1)*math.tan(math.radians(view[17]/2)))
    full_w = full_h*W/H
    center_x, center_y = (left+right)/2, (bottom+top)/2
    distance = full_h/(2*math.tan(math.radians(view[17]/2)))
    full = view.copy()
    full[9:12] = [-center_x, -center_y, -distance]
    full[15:17] = [distance-hi[2]-margin_A, distance-lo[2]+margin_A]
    c.viewport(W, H)
    c.set_view(full)
    full = list(c.get_view())
    full_h = 2*abs(full[11])*math.tan(math.radians(full[17]/2))
    full_w = full_h*W/H
    crop_h = max(crop_aspect[1], round(H*desired_h/full_h/crop_aspect[1])*crop_aspect[1])
    crop_w = crop_h*crop_aspect[0]//crop_aspect[1]
    cx = W/2+(target_x+full[9])*W/full_w
    cy = H/2-(target_y+full[10])*H/full_h
    x0, y0 = round(cx-crop_w/2), round(cy-crop_h/2)
    box = [x0, y0, x0+crop_w, y0+crop_h]
    assert 0 <= box[0] < box[2] <= W and 0 <= box[1] < box[3] <= H
    # A direct PyMOL focus view matches the pixel crop without narrowing the
    # near/far clipping planes or moving the camera into the protein.
    focus = full.copy()
    focus[9] -= ((box[0]+box[2])/2-W/2)*full_w/W
    focus[10] -= (H/2-(box[1]+box[3])/2)*full_h/H
    focus[17] = math.degrees(2*math.atan(math.tan(math.radians(full[17]/2))*crop_h/H))
    c.viewport(crop_w, crop_h)
    c.set_view(focus)
    focus = list(c.get_view())
    pixels = [[W/2+(p[0]+full[9])*W/full_w,
               H/2-(p[1]+full[10])*H/full_h] for p in rotated]
    bounds = [min(p[0] for p in pixels), min(p[1] for p in pixels),
              max(p[0] for p in pixels), max(p[1] for p in pixels)]
    assert bounds[0] > 0 and bounds[1] > 0 and bounds[2] < W and bounds[3] < H
    ca = c.count_atoms(full_selection + " and name CA")
    ca_cartoon = c.count_atoms(full_selection + " and name CA and rep cartoon")
    assert ca == ca_cartoon, (full_selection, ca, ca_cartoon)
    return dict(full_selection=full_selection, full_size=[W, H],
                crop_box=box, crop_size=[crop_w, crop_h], full_view=full,
                focus_view=focus, full_atom_projected_bounds_px=bounds,
                all_atoms_inside_full_frame=True, all_atoms_inside_depth_interval=True,
                full_CA_count=ca, full_CA_cartoon_count=ca_cartoon,
                depth_planes_identical_in_full_and_focus=True,
                operation="PyMOL full-protein ray image, then PIL.Image.crop; no resampling")


def render_crop(c, plan, full_path, crop_path, dpi=600, preview_path=None):
    full_path, crop_path = Path(full_path), Path(crop_path)
    c.viewport(*plan["full_size"])
    c.set_view(plan["full_view"])
    c.png(str(full_path), width=plan["full_size"][0], height=plan["full_size"][1],
          dpi=dpi, ray=1, quiet=1)
    with Image.open(full_path) as original:
        assert list(original.size) == plan["full_size"]
        cropped = original.crop(tuple(plan["crop_box"]))
        cropped.save(crop_path, dpi=(dpi, dpi))
        if preview_path is not None:
            preview = original.convert("RGB")
            draw = ImageDraw.Draw(preview)
            draw.rectangle(plan["crop_box"], outline="#B2182B", width=max(3, original.width//500))
            preview.thumbnail((1500, 1500))
            preview.save(preview_path)
    with Image.open(full_path) as original, Image.open(crop_path) as cropped:
        assert original.crop(tuple(plan["crop_box"])).tobytes() == cropped.tobytes()
    plan.update(full_image=full_path.name, crop_image=crop_path.name,
                full_image_sha256=sha256(full_path), crop_image_sha256=sha256(crop_path),
                exact_pixel_crop_verified=True, dpi=dpi)
    c.viewport(*plan["crop_size"])
    c.set_view(plan["focus_view"])
    return plan
