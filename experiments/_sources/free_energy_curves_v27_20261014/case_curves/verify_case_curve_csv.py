"""Independent CSV reconstruction: standard-library histogram/log arithmetic.

The extraction script is not imported. File hashes, prediction ranks, all bins,
per-origin widths and window entropy are checked from saved numeric sources.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import csv
import math
import hashlib
import json
from collections import Counter, defaultdict
import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(name):
    with (HERE / name).open(encoding='utf8', newline='') as handle:
        return list(csv.DictReader(handle))


def close(actual, expected, tolerance=1e-12):
    assert abs(float(actual) - float(expected)) <= tolerance, (actual, expected)


def count_angles(rows, chi, origin=0.):
    counter = Counter(math.floor((float(r[chi + '_degrees']) + 180 - origin) % 360 / 5) for r in rows)
    assert all(0 <= b < 72 for b in counter)
    assert sum(counter.values()) == len(rows)
    return counter


plan = json.loads((HERE / 'analysis_plan.json').read_text(encoding='utf8'))
audit = json.loads((HERE / 'extraction_audit.json').read_text(encoding='utf8'))
assert audit['status'] == 'complete' and audit['plan_sha256'] == sha(HERE / 'analysis_plan.json')
for path, digest in plan['source_hash_files'].items(): assert sha(path) == digest
for name, item in audit['outputs'].items():
    assert sha(HERE / name) == item['sha256']
    assert len(read(name)) == item['rows']

inventory = read('selection_inventory_prediction_only.csv')
summary = {r['residue']: r for r in read('selection_inventory_with_MD_outcomes.csv')}
cache_path = next(Path(p) for p in plan['source_hash_files'] if 'neural_exports' in p)
cache = np.load(cache_path, allow_pickle=False)
for row in inventory:
    index = int(row['residue_index'])
    assert str(cache['chain'][index]) == 'A' and str(cache['amino_acid'][index]) == 'Q' and cache['valid'][index, 1]
    assert row['residue'] == 'Q' + str(cache['author_residue_id'][index])
    close(row['predicted_joint_Hchi_nats'], cache['mean'][index, 1])
ordered = sorted(inventory, key=lambda r: (float(r['predicted_joint_Hchi_nats']), int(r['residue_index'])))
assert ordered[0]['residue'] == plan['low_predicted_residue'] == 'Q41'
assert ordered[-1]['residue'] == plan['high_predicted_residue'] == 'Q2'
assert [r['residue'] for r in inventory] == [r['residue'] for r in plan['prediction_inventory']]

angles = {}
for grid, name in [('primary_501frames', 'angles_501frames_all6GLN.csv'), ('label_source_5000frames', 'angles_5000frames_all6GLN.csv')]:
    grouped = defaultdict(list)
    for row in read(name): grouped[row['residue']].append(row)
    angles[grid] = grouped
    assert sorted(grouped) == sorted(summary)
    for site, rows in grouped.items():
        if grid == 'primary_501frames':
            assert len(rows) == 501
            assert [int(r['frame_zero_based']) for r in rows] == list(range(0, 1001, 2))
            for i, row in enumerate(rows): close(row['time_ns'], .2 * i)
        else:
            assert len(rows) == 5000
            assert [int(r['frame_zero_based']) for r in rows] == list(range(5000))
            for i, row in enumerate(rows): close(row['time_ns'], .1 * i)
for site, short in angles['primary_501frames'].items():
    long = angles['label_source_5000frames'][site]
    for row in short:
        previous = long[int(row['frame_zero_based'])]
        for key in row: assert row[key] == previous[key]

checked_bins = 0
max_Fstar_error = 0.
for name in ('histograms_501frames.csv', 'histograms_5000frames.csv'):
    grouped = defaultdict(list)
    for row in read(name): grouped[(row['grid'], row['residue'], row['window'], row['chi'])].append(row)
    for (grid, site, window, chi), rows in grouped.items():
        raw = angles[grid][site]
        if window not in ('pooled', '0'):
            w = int(window) - 1
            raw = raw[w * 1000:(w + 1) * 1000]
        counter = count_angles(raw, chi)
        maximum = max(counter.values())
        assert len(rows) == 72
        for i, row in enumerate(rows):
            assert int(row['count']) == counter.get(i, 0)
            assert int(row['n_frames']) == len(raw)
            close(row['angle_lower_degrees'], -180 + i * 5)
            close(row['angle_upper_degrees'], -175 + i * 5)
            close(row['angle_center_degrees'], -177.5 + i * 5)
            count = counter.get(i, 0)
            close(row['probability'], count / len(raw))
            close(row['density_per_degree'], count / len(raw) / 5)
            if count:
                expected = math.log(maximum / count)
                close(row['F_star'], expected)
                max_Fstar_error = max(max_Fstar_error, abs(float(row['F_star']) - expected))
                assert row['low_Fstar_le1'] == str(expected <= 1.)
            else:
                assert row['F_star'] == '' and row['low_Fstar_le1'] == 'False'
            checked_bins += 1

widths = defaultdict(list)
per_chi = defaultdict(list)
for row in read('width_per_angle_origin.csv'):
    grid, site, chi, origin = row['grid'], row['residue'], row['chi'], float(row['origin_degrees'])
    counter = count_angles(angles[grid][site], chi, origin)
    maximum = max(counter.values())
    count_low = sum(math.log(maximum / count) <= 1. for count in counter.values())
    width = count_low * 5
    assert width == float(row['width_degrees'])
    assert maximum == int(row['maximum_bin_count'])
    assert len(counter) == int(row['occupied_bins']) and count_low == int(row['qualifying_bins'])
    widths[(grid, site)].append(width)
    per_chi[(grid, site, chi)].append(width)
for (grid, site), items in widths.items():
    assert len(items) == 12
    column = 'MD_width_501frames_deg' if grid == 'primary_501frames' else 'MD_width_5000frames_deg'
    close(summary[site][column], math.fsum(items) / len(items))
for site in summary:
    for q in range(1, 4):
        items = per_chi[('primary_501frames', site, f'chi{q}')]
        assert len(items) == 4
        close(summary[site][f'MD_width_501_chi{q}_deg'], math.fsum(items) / len(items))

joint_windows = {(r['residue'], int(r['window'])): r for r in read('joint_entropy_original_five_windows.csv')}
max_entropy_error = 0.
for site, rows in angles['label_source_5000frames'].items():
    values = []
    for w in range(5):
        raw = rows[w * 1000:(w + 1) * 1000]
        # Saved float32 degrees are restored exactly, preserving original endpoint arithmetic.
        points = np.array([[float(r[f'chi{q}_degrees']) for q in range(1, 4)] for r in raw], np.float32)
        bins = np.clip(np.floor((points + np.float32(180)) / np.float32(120)).astype(int), 0, 2)
        counter = Counter(tuple(row) for row in bins.tolist())
        entropy = math.fsum(-(count / 1000) * math.log(count / 1000) for count in counter.values())
        close(joint_windows[(site, w + 1)]['Hchi_nats'], entropy)
        max_entropy_error = max(max_entropy_error, abs(float(joint_windows[(site, w + 1)]['Hchi_nats']) - entropy))
        values.append(entropy)
    mean = math.fsum(values) / 5
    variance = math.fsum((x - mean) ** 2 for x in values) / 4
    close(summary[site]['MD_joint_Hchi_mean5windows_nats'], mean)
    close(summary[site]['original_joint_label_variance'], variance)
    close(summary[site]['frozen_MD_joint_Hchi_target_nats'], mean, 2e-6)
    close(summary[site]['prediction_absolute_error_nats'], abs(float(summary[site]['predicted_joint_Hchi_nats']) - mean))

result = dict(status='passed', verification_code_sha256=sha(__file__), extraction_audit_sha256=sha(HERE / 'extraction_audit.json'),
    source_sha256_all_match=True, model_selection_and_all6_inventory_verified=True,
    exact_primary_grid_and_original_source_subset_verified=True, histogram_bins_checked=checked_bins,
    maximum_Fstar_numeric_error=max_Fstar_error, original_windows_checked=len(joint_windows),
    maximum_joint_entropy_numeric_error=max_entropy_error, per_origin_width_rows_checked=len(read('width_per_angle_origin.csv')),
    all6_original_joint_label_mean_and_sample_variance_verified=True,
    files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in HERE.iterdir() if p.is_file() and p.name != 'independent_csv_verification.json'})
(HERE / 'independent_csv_verification.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
print(json.dumps({k: v for k, v in result.items() if k != 'files'}, ensure_ascii=False, indent=2))
