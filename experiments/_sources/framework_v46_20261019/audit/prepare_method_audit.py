"""Record independently inspected source evidence and exact prose clarifications."""
from pathlib import Path
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
PAPER = ROOT.parents[2]
CODE = Path(r'E:\Project\FlexPPI')
WORKTREE = Path(r'C:\Users\DELL\.codex\worktrees\7471\FlexPPI')
DOCX = PAPER / 'paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx'
W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
M = '{http://schemas.openxmlformats.org/officeDocument/2006/math}'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

with zipfile.ZipFile(DOCX) as z:
    body = ET.fromstring(z.read('word/document.xml')).find(W+'body')
paragraphs = body.findall(W+'p')
text = lambda p: ''.join(x.text or '' for x in p.iter(W+'t'))
replacements = []
changes = {
    8: ('同一编码器处理输入中的全部有效链，两个读出头分别生成RMSF与Hχ的均值和残差尺度。',
        '同一编码器处理输入中的全部有效链，联合读出层为RMSF与Hχ分别生成均值和残差尺度。'),
    12: ('两个读出头分别生成目标均值及用于异方差回归的残差尺度。',
         '联合读出层为两项目标分别生成均值及用于异方差回归的残差尺度。'),
    29: ('默认联合模型共享几何编码器，并为RMSF和Hχ保留各自的均值与对数方差读出参数。',
         '默认联合模型采用共享几何编码器与联合读出层，为RMSF和Hχ分别输出均值与对数方差。'),
}
for index, (oldpart, newpart) in changes.items():
    p = paragraphs[index]
    assert not list(p.iter(M+'oMath'))
    old = text(p)
    assert oldpart in old
    new = old.replace(oldpart, newpart, 1)
    if index == 29:
        assert '方差头提供式(3)中的残差尺度。' in new
        new = new.replace('方差头提供式(3)中的残差尺度。', '预测方差提供式(3)中的残差尺度。')
    replacements.append({'prefix': f'P{index}', 'old_exact': old, 'new_exact': new})
(ROOT/'prose_clarifications.json').write_text(json.dumps({
    'source_document': str(DOCX), 'source_sha256': sha(DOCX),
    'rationale': 'The implementation has one StateHead MLP with shared intermediate mapping and target-specific final output channels, rather than two independent MLP heads.',
    'replacements': replacements,
}, ensure_ascii=False, indent=2), encoding='utf-8')

files = ['src/flexppi/model/network.py', 'src/flexppi/model/heads.py', 'src/flexppi/model/geometry.py',
         'src/flexppi/model/ipa.py', 'src/flexppi/model/sequence.py', 'src/flexppi/losses.py',
         'src/flexppi/data/labels.py', 'src/flexppi/data/structures.py']
file_checks = []
for relative in files:
    a, b = CODE/relative, WORKTREE/relative
    file_checks.append({'path': str(a), 'sha256': sha(a), 'worktree_byte_identical': a.read_bytes() == b.read_bytes(),
                        'worktree_text_identical': a.read_text(encoding='utf-8') == b.read_text(encoding='utf-8')})
assert all(x['worktree_text_identical'] for x in file_checks)
configs = []
for fold in range(5):
    path = CODE / f'artifacts/performance_first_v2/fold_{fold}/resolved_config.json'
    config = json.loads(path.read_text(encoding='utf-8'))
    c = config['model']
    assert c['sequence_backend'] == 'learned'
    assert (c['ipa_blocks'],c['ipa_heads'],c['single_dim'],c['pair_dim']) == (4,8,128,48)
    assert not c['use_condition_features'] and not c['enable_response_head']
    assert c['metrics'] == ['rmsf','chi_entropy']
    configs.append({'fold':fold,'path':str(path),'sha256':sha(path),'model':c})

def evidence(relative, first, last):
    return {'path':str(CODE/relative),'line_start':first,'line_end':last}

audit = {
    'status':'passed', 'scope':'Architecture/content audit for Figure 1 redraw; no training, inference, or MD performed.',
    'source_document':str(DOCX), 'source_document_sha256':sha(DOCX),
    'code_files':file_checks, 'actual_five_fold_resolved_configs':configs,
    'input': {
        'truth':'Amino-acid identity, chain membership, continuous within-chain residue order, N/CA/C coordinates. Learned AA embeddings only; no ESM, MD input, side-chain coordinate input, temperature or dataset conditioning in the evaluated model.',
        'single_features':'128-dimensional learned AA embedding concatenated with an independently learned 32-dimensional type embedding, linearly projected to 128 dimensions.',
        'pair_features':'29 raw features: 24 Cα-distance RBF channels (centres 0–24 Å); same-chain flag; abs within-chain position difference clipped to 32 then divided by 32; unit CA_j−CA_i direction in query-i local frame. Two-layer MLP projects to 48 dimensions.',
        'cross_chain_order':'The ordinal difference is computed across chains too but disambiguated by same-chain flag; do not label it a meaningful cross-chain sequence separation.',
        'evidence':[evidence('src/flexppi/model/network.py',27,45),evidence('src/flexppi/model/network.py',68,90),evidence('src/flexppi/model/sequence.py',15,28),evidence('src/flexppi/model/geometry.py',59,99),evidence('src/flexppi/data/structures.py',115,135)],
        'compact_labels':['Amino-acid identity','N / Cα / C','Residue features, 128 d','Pair geometry, 48 d','Distance · chain relation · local direction'],
    },
    'frames':{
        'truth':'t_i = CA_i; x_i = normalize(C_i−CA_i); y_i = normalize((N_i−CA_i)−dot(N_i−CA_i,x_i)x_i); z_i = normalize(x_i×y_i); R_i = [x_i y_i z_i] (basis vectors as columns).',
        'transforms':'p_global = R_i q_local + t_i; q_local = R_i^T (p_global−t_i).',
        'fixed':'Input coordinates, R_i/t_i and pair features are computed once and held fixed through all four blocks; only residue representations update.',
        'global_origin':'The model uses the supplied coordinate reference; it does not require CA-centroid centring. A drawing may use an arbitrary clearly specified global origin. The old caption O=all-CA mean is a drawing convention and is not an encoder operation.',
        'evidence':[evidence('src/flexppi/model/geometry.py',13,55),evidence('src/flexppi/model/network.py',79,90)],
        'compact_labels':['Local frame at Cα','x: Cα → C','y: orthogonalized Cα → N','z = x × y','Local → global','Global → query frame','Fixed input frames'],
    },
    'ipa':{
        'truth':'Each of four pre-normalized residual blocks has eight attention heads. Each head uses 16 scalar channels, 4 query/key latent points and 8 value latent points.',
        'logits':'Scaled scalar q_i·k_j plus learned pair bias plus negative positively-weighted squared distances between query-i and key-j latent points after transforming each from its own fixed frame to global coordinates. Masked softmax normalizes across key residues j.',
        'aggregation':'The same per-head attention weights aggregate scalar values v_j, pair features z_ij and global value points. Weighted global points are transformed back to the query-i local frame; local x/y/z coordinates and point norms are concatenated with scalar and pair outputs. A linear map projects this merged result to the residue feature dimension.',
        'update':'Pre-LayerNorm → IPA → dropout + residual; then pre-LayerNorm → Linear(128,512) → SiLU → Linear(512,128) → dropout + residual. Final LayerNorm after four blocks.',
        'latent_point_boundary':'Q/K/V points are learned latent geometric probes; they are not physical atoms or generated protein coordinates. Invariance is to common rigid transformations, not to changing relative chain arrangements.',
        'evidence':[evidence('src/flexppi/model/ipa.py',35,46),evidence('src/flexppi/model/ipa.py',57,108),evidence('src/flexppi/model/ipa.py',123,148),evidence('src/flexppi/model/network.py',46,61)],
        'compact_labels':['Scalar similarity','Pair bias','Point distance','Attention weights','Weighted sums','To query frame','Merge + project','Residual update','IPA block × 4','8 attention heads'],
    },
    'readout':{
        'truth':'One StateHead MLP (LayerNorm → Linear(128,128) → SiLU → Linear(128,4)) serves both metrics. Shared intermediate mapping, target-specific final output channels. Outputs for each target are softplus mean and clipped predicted log variance [−8,6].',
        'avoid':'Do not draw two independent multi-layer prediction heads. A single Joint readout box may feed two labelled output ports RMSF (μ, log σ²) and Hχ (μ, log σ²).',
        'evidence':[evidence('src/flexppi/model/network.py',60,61),evidence('src/flexppi/model/network.py',89,95),evidence('src/flexppi/model/heads.py',9,26)],
        'compact_labels':['Joint readout','RMSF: μ, log σ²','Hχ: μ, log σ²'],
    },
    'supervision':{
        'truth':'MD windows or independent repeats produce a target mean y and observed sample variance v_obs per residue/metric. Both enter heteroscedastic Gaussian NLL. MD observed variance is not directly regressed by the model variance output.',
        'formula':'v_total = exp(log σ²_pred) + max(v_obs,0); NLL = 0.5[(μ_pred−y)^2/v_total + ln(v_total)] (v_total numerically clamped to [1e−6,1e6]).',
        'labels':'RMSF: chain-wise fitted Cα displacements, Å. Hχ: entropy of joint fixed-bin χ states, natural log, nats; not sum of marginal entropies and not a χ angle prediction.',
        'reduction':'Mean over valid sites per record/metric, then equal mean over available metrics, then mean over records. Target mask limits supervision, not all-chain context encoding.',
        'evidence':[evidence('src/flexppi/data/labels.py',174,189),evidence('src/flexppi/data/labels.py',291,313),evidence('src/flexppi/losses.py',18,49),evidence('src/flexppi/losses.py',75,84)],
        'compact_labels':['MD observations','Target mean y','Observed variance v_obs','Noise-aware Gaussian NLL'],
    },
    'old_figure_risks':[
        'Old main prose says two readout heads, which can imply independent MLPs; clarify one joint layer with target-specific output channels.',
        'Old panel b spends most space on global axes, translation and three local frames, while it does not show how N/CA/C define a frame; use one N–CA–C inset and consistent true-coordinate axes.',
        'Global origin at the CA centroid is a display convention, not preprocessing required by the model; do not add a centring operation to the architecture.',
        'Pair-feature grid must not imply pair updates or triangle operations; pair representation is fixed per input through blocks.',
        'Point branch must explicitly return to query local coordinates before merge, and latent points must not appear to be predicted atoms.',
        'Old MD supervision mean+variance → NLL path is too compressed; show observed variance and predicted variance combine, not an arrow that implies direct variance regression.',
        'Old example timestamps 0/100/200 ns could imply one common observation duration; source datasets have distinct observation designs. Prefer unnumbered MD snapshots or a correctly labelled specific illustration.',
    ],
    'protected_data':'Existing 4g9f 840-residue predictions and per-chain 11-residue smoothing must be preserved, including all residue positions and actual chain boundaries. Final independent check is pending redesigned assets.',
    'prose_clarification_file':str(ROOT/'prose_clarifications.json'),
}
(ROOT/'method_content.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','method_audit':str(ROOT/'method_content.json'),'prose_clarifications':str(ROOT/'prose_clarifications.json'),'replacement_count':len(replacements),'five_actual_fold_configs_checked':len(configs)},ensure_ascii=False))
