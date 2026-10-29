"""Translate labels only; inherit original plotting functions and exact inputs."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import shutil
import sys
import zipfile

BASE = Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010')
OUT = Path(__file__).resolve().parent
RUNTIME = BASE / 'experiment_completion/completed_results_20261012/plot_runtime'
sys.path.insert(0, str(RUNTIME))
import numpy as np
import pandas as pd
from PIL import Image, ImageChops


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replace_exact(source, mapping):
    for old, new in mapping.items():
        assert old in source, old
        source = source.replace(old, new)
    return source


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def inspect_visible_english(fig):
    from matplotlib.text import Text
    labels = [t.get_text() for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    chinese = [x for x in labels if re.search(r'[\u3400-\u9fff]', x)]
    assert not chinese, chinese
    return labels


def install_export_audit(mod):
    old = mod.export

    def checked_export(fig, *args, **kwargs):
        visible = inspect_visible_english(fig)
        old(fig, *args, **kwargs)
        stem = args[0]
        (mod.OUT / (stem + '_visible_labels.json')).write_text(
            json.dumps(visible, ensure_ascii=False, indent=2), encoding='utf-8')

    mod.export = checked_export


def compare_csv(new, old):
    a = pd.read_csv(new)
    b = pd.read_csv(old)
    assert a.shape == b.shape and list(a.columns) == list(b.columns), (new, a.shape, b.shape)
    differing_text = []
    for c in a:
        if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
            assert np.array_equal(a[c].to_numpy(), b[c].to_numpy(), equal_nan=True), (new, c)
        elif not a[c].fillna('').equals(b[c].fillna('')):
            pairs = list(zip(a.loc[a[c].fillna('') != b[c].fillna(''), c], b.loc[a[c].fillna('') != b[c].fillna(''), c]))
            assert c == 'label' and all(x == 'Unclassified' and y == '未分类' for x, y in pairs), (new, c, pairs[:3])
            differing_text.append(c)
    return {'file': new.name, 'rows': len(a), 'columns': len(a.columns), 'all_numeric_cells_unchanged': True,
            'translated_text_columns': differing_text, 'new_sha256': sha(new), 'original_sha256': sha(old)}


def main():
    for name in ['fig2', 'fig3', 'fig4', 'figS1']:
        (OUT / name).mkdir(parents=True, exist_ok=True)
    original = BASE / 'figures/redraw_quantitative.py'
    source = original.read_text(encoding='utf-8')
    source = replace_exact(source, {
        "font_manager.findfont('Microsoft YaHei', fallback_to_default=False)": "font_manager.findfont('Arial', fallback_to_default=False)",
        "'font.sans-serif':['Microsoft YaHei', 'Arial']": "'font.sans-serif':['Arial', 'Helvetica', 'Liberation Sans']",
        "'累计记录比例'": "'Fraction of records'",
        "'MAE (Å；近零线性，远端对数)'": "'MAE (Å)'",
        "'：整体误差分布'": "': error distribution'",
        "'记录内 Spearman ρ'": "'Within-record Spearman ρ'",
        "'：残基排序相关'": "': residue ranking'",
        "'实际输入的蛋白链数'": "'Number of input chains'",
        "'：按输入链数分层'": "': effect of chain count'",
        "'记录内参考 '+TN[task]+' 均值 ('": "'Mean reference '+TN[task]+' ('",
        "'Full MAE / 参考均值 (%)'": "'Full MAE / reference mean (%)'",
        "'记录内参考 RMSF 均值 (Å；对数)'": "'Mean reference RMSF (Å)'",
        "ax.set_ylabel(ax.get_ylabel()+'；对数')": "ax.set_ylabel(ax.get_ylabel())",
        "'：动态幅度与绝对误差'": "': amplitude and error'",
        "'：相对于参考幅度的误差'": "': relative error'",
        "'氨基酸类型（按 χ 角数排列）'": "'Amino acid (ordered by χ count)'",
        "'Full：同类残基内的预测误差'": "'Full: within-type error'",
        "'相对对照的 MAE 降幅 (%)'": "'MAE reduction (%)'",
        "'MAE 降幅 (%)\\n近零线性，远端对数'": "'MAE reduction (%)'",
        "'同类残基：Full 相对对照的收益'": "'Full: within-type gain'",
    })
    # Keep every inherited function, but call only the requested figures.
    source = source.replace("if __name__=='__main__':main()", "if __name__=='__main__' and False:main()")
    script = OUT / 'fig2/redraw_quantitative_english.py'
    script.write_text(source, encoding='utf-8')
    shutil.copyfile(original, OUT / 'fig2/redraw_quantitative_original.py')
    shutil.copyfile(script, OUT / 'fig3/redraw_quantitative_english.py')
    mod = load_module(script, 'english_quant')
    install_export_audit(mod)
    per, chains, labels, aa, pairs, external = mod.prepare()
    mod.fig2(per)
    mod.OUT = OUT / 'fig3'
    mod.fig3(per, labels, aa, pairs)
    for name in ['fig2', 'fig3']:
        folder = OUT / name
        items = {k: v for k, v in mod.MANIFEST['figures'].items() if k.startswith(name + '_')}
        (folder / 'figure_manifest.json').write_text(json.dumps(items, ensure_ascii=False, indent=2), encoding='utf-8')

    f4_original = BASE / 'experiment_completion/energy_landscape_v21_20261012/figures/plot_joint_complementarity.py'
    source = f4_original.read_text(encoding='utf-8')
    # Only the active new Figure 4 function is rendered; inherited functions remain complete.
    source = source.replace("'font.sans-serif':['Microsoft YaHei', 'Arial']", "'font.sans-serif':['Arial', 'Helvetica', 'Liberation Sans']")
    source = replace_exact(source, {
        "data=OUT.parent/'oof_analysis'": "data=Path(r'E:\\Project\\FlexPPI\\experiments\\_sources\\energy_landscape_v21_20261012\\oof_analysis')",
        "label='轨迹参考'": "label='MD reference'",
        "label='折外预测'": "label='Prediction'",
        "label='逐记录'": "label='Per record'",
        "label='整体均值'": "label='Mean'",
        "'两读数的 Spearman ρ'": "'RMSF–Hχ Spearman ρ'",
        "'累计记录比例'": "'Fraction of records'",
        "'均值 0.391'": "'Mean 0.391'",
        "'均值 0.517'": "'Mean 0.517'",
        "'3,063 条记录'": "'3,063 records'",
        "'两指标的记录内相关'": "'Within-record correlation'",
        "label='未分类'": "label='Unclassified'",
        "'记录等权占比 (%)'": "'Mean record fraction (%)'",
        "f'未分类：参考 {means_ref.unclassified:.1f}%\\n预测 {means_pred.unclassified:.1f}%'": "f'Unclassified: reference {means_ref.unclassified:.1f}%\\nprediction {means_pred.unclassified:.1f}%'",
        "'同类型位点的组合占比'": "'Dynamic combinations'",
        "'参考组合／均值／记录数'": "'Reference class / mean / records'",
        "'组合召回率 (%)'": "'Combination recall (%)'",
        "'四类组合的预测召回'": "'Combination recall'",
    })
    source = source.replace("if __name__=='__main__':fig4_joint_complementarity()", "if __name__=='__main__' and False:fig4_joint_complementarity()")
    script = OUT / 'fig4/plot_joint_complementarity_english.py'
    script.write_text(source, encoding='utf-8')
    shutil.copyfile(f4_original, OUT / 'fig4/plot_joint_complementarity_original.py')
    mod4 = load_module(script, 'english_joint')
    install_export_audit(mod4)
    mod4.fig4_joint_complementarity()

    # S1 in the active document already uses English and needs no redraw.
    oldroot = Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002')
    s1_vector = oldroot / 'figures/writing_revision_20261004/figS_corrected_observed_label_variation'
    for ext in ['png', 'pdf', 'svg']:
        shutil.copyfile(str(s1_vector) + '.' + ext, OUT / ('figS1/figS1_observed_label_variation.' + ext))
    s1source = oldroot / 'figures/corrected_label_variation'
    for name in ['plot_label_variation.py', 'per_record_labels.csv', 'by_source_labels.csv',
                 'by_source_amino_acid_labels.csv', 'label_analysis_audit.json',
                 'plot_ecdf_all_record_points.csv', 'plot_hchi_normalized_mean_matrix.csv',
                 'plot_zero_variance_count_evidence.csv']:
        shutil.copyfile(s1source / name, OUT / 'figS1' / name)
    shutil.copyfile(oldroot / 'qa/writing_revision_20261004/polish_figure_text.py', OUT / 'figS1/polish_figure_text.py')
    with zipfile.ZipFile(BASE / 'paper_v38_20261016_中文修订_图7标注说明与图例调整_三线表.docx') as z:
        active = z.read('word/media/image7.png')
    png = OUT / 'figS1/active_v38_image7.png'
    png.write_bytes(active)
    a = Image.open(png).convert('RGBA')
    b = Image.open(OUT / 'figS1/figS1_observed_label_variation.png').convert('RGBA')
    assert a.size == b.size and ImageChops.difference(a, b).getbbox() is None
    import xml.etree.ElementTree as ET
    tree = ET.parse(OUT / 'figS1/figS1_observed_label_variation.svg')
    texts = [''.join(t.itertext()) for t in tree.findall('.//{http://www.w3.org/2000/svg}text')]
    assert not any(re.search(r'[\u3400-\u9fff]', x) for x in texts)
    s1audit = {'status': 'passed', 'already_English_in_v38': True, 'active_document_image': 'word/media/image7.png',
               'pixel_identical_to_vector_companion_png': True, 'pixel_size': list(a.size),
               'translation_required': False, 'all_data_and_layout_unchanged': True,
               'visible_labels': texts, 'png_sha256': sha(OUT / 'figS1/figS1_observed_label_variation.png')}
    (OUT / 'figS1/data_and_language_audit.json').write_text(json.dumps(s1audit, ensure_ascii=False, indent=2), encoding='utf-8')

    results = []
    for name in ['fig2', 'fig3', 'fig4']:
        comparisons = []
        oldfolder = BASE / 'figures' if name in ['fig2', 'fig3'] else BASE / 'experiment_completion/energy_landscape_v22_maintext_20261012/figures'
        for file in (OUT / name).glob(name + '_*.csv'):
            old = oldfolder / file.name
            if old.exists():
                comparisons.append(compare_csv(file, old))
        assert comparisons
        report = {'status': 'passed', 'edit_scope': 'English labels and font only; original plotted values and transforms preserved',
                  'data_comparisons': comparisons, 'random_sampling': False, 'new_inference': False,
                  'source_script_sha256': sha(original if name in ['fig2', 'fig3'] else f4_original),
                  'all_numeric_cells_unchanged': True, 'all_visible_labels_English': True,
                  'visual_QA': 'pending full image review'}
        (OUT / name / 'data_and_language_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
        results.append({'figure': name, 'compared_csvs': len(comparisons), 'rows_compared': sum(c['rows'] for c in comparisons)})
    print(json.dumps({'status': 'rendered; pending visual review', 'figures': results, 'S1': 'already English; retained pixel-identical'}, indent=2))


if __name__ == '__main__':
    main()
