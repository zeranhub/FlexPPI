"""Exact 47-feature definitions and balancing weights from the historical static probe.
Original source SHA256: 1c06cddde24af07824e553a654397bb031210a1b52db90bcd014dafbd7aa9dcd
No label, dataset or simulation-condition features are used.
"""
import numpy as np
from scipy.spatial.distance import cdist

def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-8)
    v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
    return np.arctan2((np.cross(b1, v) * w).sum(-1), (v * w).sum(-1))

def features(s):
    n = s.num_residues
    names, values = [], []

    def add(name, v):
        names.append(name)
        values.append(np.asarray(v, np.float32))

    # Sequence features see the same parsed sequence as the frozen model.
    for delta in range(-4, 5):
        v = np.full(n, 21, np.float32)  # 21 is a chain-boundary token, not a residue.
        for chain in np.unique(s.chain_index):
            ids = np.flatnonzero(s.chain_index == chain)
            positions = np.arange(len(ids)) + delta
            ok = (positions >= 0) & (positions < len(ids))
            v[ids[ok]] = s.aatype[ids[positions[ok]]]
        add("aa_offset_%+d" % delta, v)
    seq_cat = names.copy()
    length = np.zeros(n)
    fraction = np.zeros(n)
    end_distance = np.zeros(n)
    composition = np.zeros((n, 20))
    for chain in np.unique(s.chain_index):
        ids = np.flatnonzero(s.chain_index == chain)
        length[ids] = len(ids)
        fraction[ids] = np.arange(len(ids)) / max(len(ids) - 1, 1)
        end_distance[ids] = np.minimum(np.arange(len(ids)), np.arange(len(ids))[::-1])
        composition[ids] = np.bincount(s.aatype[ids], minlength=21)[:20] / len(ids)
    add("log_chain_length", np.log1p(length))
    add("fraction_chain_position", fraction)
    add("log_distance_to_chain_end", np.log1p(end_distance))
    positional = names[len(seq_cat):].copy()
    for a in range(20):
        add("chain_aa_fraction_%02d" % a, composition[:, a])
    seq = names.copy()

    # All geometry comes from N/CA/C, never side-chain atoms or trajectories.
    ca = s.backbone[:, 1].astype(np.float64)
    dist = cdist(ca, ca)
    same = s.chain_index[:, None] == s.chain_index[None, :]
    sep = np.abs(s.residue_index[:, None] - s.residue_index[None, :])
    valid_pair = s.frame_mask[:, None] & s.frame_mask[None, :]
    before = len(names)
    for radius in (6, 8, 10, 12):
        add("intra_contacts_%dA" % radius, np.log1p(((dist < radius) & same & (sep > 2) & valid_pair).sum(1)))
    for radius in (8, 12):
        add("intra_longrange_contacts_%dA" % radius, np.log1p(((dist < radius) & same & (sep >= 12) & valid_pair).sum(1)))
    phi, psi = np.full(n, np.nan), np.full(n, np.nan)
    for chain in np.unique(s.chain_index):
        ids = np.flatnonzero(s.chain_index == chain)
        prev, curr = ids[:-1], ids[1:]
        # Do not connect residues across a missing segment in the input structure.
        bonded = (np.linalg.norm(s.backbone[prev, 2] - s.backbone[curr, 0], axis=1) < 2.0)
        bonded &= s.frame_mask[prev] & s.frame_mask[curr]
        prev, curr = prev[bonded], curr[bonded]
        phi[curr] = dihedral(s.backbone[prev, 2], s.backbone[curr, 0], s.backbone[curr, 1], s.backbone[curr, 2])
        psi[prev] = dihedral(s.backbone[prev, 0], s.backbone[prev, 1], s.backbone[prev, 2], s.backbone[curr, 0])
    for name, angles in (("phi", phi), ("psi", psi)):
        add(name + "_sin", np.sin(angles))
        add(name + "_cos", np.cos(angles))
    v, w = s.backbone[:, 0] - ca, s.backbone[:, 2] - ca
    cosine = (v * w).sum(1) / np.maximum(np.linalg.norm(v, axis=1) * np.linalg.norm(w, axis=1), 1e-8)
    cosine[~s.frame_mask] = np.nan
    add("N_CA_C_cosine", cosine)
    intra = names[before:].copy()
    before = len(names)
    for radius in (6, 8, 10, 12):
        add("partner_contacts_%dA" % radius, np.log1p(((dist < radius) & ~same & valid_pair).sum(1)))
    partner = names[before:].copy()
    sets = {
        "sequence": seq,
        "geometry": positional + intra + partner,
        "sequence_intra": seq + intra,
        "sequence_complex": seq + intra + partner,
    }
    return np.stack(values, axis=1), names, sets, seq_cat

def weights(records, record_ids, train_mask):
    # Same expected source/system/record balance as the neural run; each metric
    # has uniform residues within a record. Unlike neural NLL, the loss is MSE.
    from collections import Counter
    counts = np.bincount(record_ids[train_mask], minlength=len(records))
    active = [r for r in records if counts[r["record_index"]] > 0]
    system_records = Counter((r["dataset"], r["biological_system"]) for r in active)
    source_systems = Counter(source for source, system in system_records)
    per_record = np.zeros(len(records))
    for r in active:
        idx = r["record_index"]
        per_record[idx] = 1 / (len(source_systems) * source_systems[r["dataset"]]
                               * system_records[(r["dataset"], r["biological_system"])] * counts[idx])
    w = per_record[record_ids[train_mask]]
    return w / w.mean()
