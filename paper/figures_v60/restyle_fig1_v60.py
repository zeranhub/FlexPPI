"""Apply a restrained native-Draw.io tint change, preserving every other token."""
from pathlib import Path
import colorsys
import hashlib
import json
import re
import shutil
import subprocess
from lxml import etree as ET

PROJECT = Path(__file__).resolve().parents[1]
BASELINE = Path('E:/Project/FlexPPI/paper/figures_v59/01_模型框架')
OUTPUT = PROJECT / 'figures_v60' / '01_模型框架'
APP = Path('C:/Program Files/draw.io/draw.io.exe')

# The lavender/pink identities stay visible. No new teal blocks or dark-blue fills.
EXPLICIT = {
    '#FCEEEF': '#FBEFEA',
    '#FFC8C4': '#F6CCBE', '#FFC4C0': '#F6C8BA',
    '#FFC5C2': '#F6C9BC', '#FFC3C0': '#F6C7BA',
    '#CF7772': '#C78370', '#D27670': '#CA826E',
    '#A99AB9': '#9FA3B9', '#A394B2': '#979DB5',
    '#E1E5F1': '#DEE6F1', '#DEE4F0': '#DCE5F0',
    '#DDDFF0': '#DAE1F0',
}

def tint(token):
    if token.upper() in EXPLICIT:
        return EXPLICIT[token.upper()]
    rgb = [int(token[i:i+2], 16) / 255 for i in [1,3,5]]
    hue, lightness, saturation = colorsys.rgb_to_hls(*rgb)
    degrees = hue * 360
    if 235 <= degrees <= 310 and saturation > .05:
        # A small 18% hue interpolation toward blue; retain the source lightness.
        hue = (degrees + .18 * (225 - degrees)) / 360
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation * .98)
    elif (degrees >= 340 or degrees <= 20) and lightness > .65 and saturation > .1:
        degrees = degrees - 360 if degrees >= 340 else degrees
        rgb = colorsys.hls_to_rgb((degrees + .18*(28-degrees))/360, lightness, saturation*.97)
    else:
        return token
    return '#' + ''.join(f'{round(v*255):02X}' for v in rgb)

def main():
    OUTPUT.mkdir(exist_ok=True, parents=True)
    old_path = BASELINE / 'fig1_framework.drawio'
    original = old_path.read_text(encoding='utf-8')
    before = ET.fromstring(original.encode('utf-8'))
    changes = []
    color_pattern = re.compile(r'(?P<key>fillColor|strokeColor|gradientColor)=(?P<hex>#[0-9a-fA-F]{6})(?=;)')
    # XML is not reserialised: only native colour properties in style attributes change.
    def style_edit(match):
        style = match.group(1)
        def color_edit(color):
            source = color.group('hex')
            target = tint(source)
            if source != target:
                changes.append({'property':color.group('key'),'from':source,'to':target})
            return color.group('key') + '=' + target
        return 'style="' + color_pattern.sub(color_edit, style) + '"'
    updated = re.sub(r'style="([^"]*)"', style_edit, original)
    after = ET.fromstring(updated.encode('utf-8'))
    left, right = before.findall('.//mxCell'), after.findall('.//mxCell')
    assert len(left) == len(right)
    changed_cells = []
    for a,b in zip(left,right):
        assert a.get('id') == b.get('id')
        aa, bb = dict(a.attrib), dict(b.attrib)
        sa,sb=aa.pop('style',''),bb.pop('style','')
        assert aa == bb, a.get('id')
        assert ET.tostring(a.find('mxGeometry')) == ET.tostring(b.find('mxGeometry')) if a.find('mxGeometry') is not None else b.find('mxGeometry') is None
        assert color_pattern.sub('COLOR',sa) == color_pattern.sub('COLOR',sb), a.get('id')
        if 'image=' in sa:
            assert sa == sb, 'Embedded image must stay byte-identical'
        if sa != sb:
            changed_cells.append(a.get('id'))
    destination = OUTPUT / old_path.name
    destination.write_text(updated,encoding='utf-8')
    shutil.copy2(BASELINE/'caption_zh.md',OUTPUT/'caption_zh.md')
    for extension in ['png','pdf','svg']:
        command=[str(APP),'--export','--format',extension,'--border','0',
                 '--output',str(OUTPUT/('fig1_framework.'+extension)),str(destination)]
        if extension=='png':
            command[1:1]=['--scale','2']
        completed=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',creationflags=subprocess.CREATE_NO_WINDOW,timeout=90)
        (OUTPUT/('export_'+extension+'.log')).write_text(completed.stdout+completed.stderr,encoding='utf-8')
        if completed.returncode:
            raise RuntimeError(completed.stdout+completed.stderr)
        assert (OUTPUT/('fig1_framework.'+extension)).is_file(), extension
    audit={'status':'source_verified_visual_pending','baseline':str(BASELINE),
        'cell_count':len(left),'changed_style_cells':changed_cells,'color_edits':changes,
        'geometry_values_images_arrows_and_fonts_unchanged':True,
        'embedded_image_styles_byte_identical':True,
        'only_native_colour_properties_changed':True,
        'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUTPUT.glob('fig1_framework.*')}}
    (OUTPUT/'tint_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'cells':len(left),'edited_cells':len(changed_cells),'colour_properties':len(changes),'output':str(OUTPUT)},ensure_ascii=False))

if __name__=='__main__':
    main()
