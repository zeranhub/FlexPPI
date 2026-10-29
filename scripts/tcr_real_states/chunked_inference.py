"""Memory-bounded exact all-key inference for the frozen structure model.

Only query rows are split. Every query attends to every valid residue; there is
no spatial cutoff, sequence truncation, sampling, or parameter change. Geometry
features are recomputed per query block to avoid keeping an N by N pair tensor.
The independent validation compares against the unchanged native forward_state.
"""
import math
import torch
from flexppi.model.geometry import frames_from_backbone, local_to_global, global_to_local


def pair_rows(encoder, batch, rotation, mask, start, end):
    ca = batch['backbone'][..., 1, :]
    difference = ca[:, None, :, :] - ca[:, start:end, None, :]
    distance = torch.sqrt(torch.sum(difference*difference, dim=-1) + 1e-8)
    rbf = torch.exp(-((distance[..., None] - encoder.centers)/encoder.width)**2)
    chains, indices = batch['chain_index'], batch['residue_index']
    same_chain = (chains[:, start:end, None] == chains[:, None, :]).to(ca.dtype)
    separation = torch.abs(indices[:, start:end, None] - indices[:, None, :]).to(ca.dtype)
    separation = torch.clamp(separation, max=32.) / 32.
    inverse = rotation[:, start:end].transpose(-1, -2)[:, :, None, :, :]
    local = torch.matmul(inverse, difference.unsqueeze(-1)).squeeze(-1)
    direction = local / torch.sqrt(torch.sum(local*local, dim=-1, keepdim=True)+1e-8)
    raw = torch.cat((rbf, same_chain[..., None], separation[..., None], direction), dim=-1)
    pair = encoder.projection(raw)
    valid = mask[:, start:end, None] & mask[:, None, :]
    return pair * valid[..., None].to(pair.dtype)


@torch.inference_mode()
def forward_chunked(model, batch, query_chunk=48):
    if model.training:
        raise ValueError('Frozen eval mode is required')
    if query_chunk < 1:
        raise ValueError('query_chunk must be positive')
    padding = batch['padding_mask'].bool()
    mask = batch['geometry_mask'].bool() & padding
    sequence = model.sequence_encoder(batch['aatype'], batch['chain_index'], padding)
    single = model.single_input(torch.cat((sequence, model.aa_embedding(batch['aatype'])), dim=-1))
    if model.condition_mlp is not None:
        single = single + model.condition_mlp(batch['condition'])[:, None, :]
    rotation, translation, frame_mask = frames_from_backbone(batch['backbone'], batch['backbone_mask'])
    mask = mask & frame_mask
    b, n, _ = single.shape
    for block in model.blocks:
        ipa = block.ipa
        normalized = block.input_norm(single)
        q = ipa.linear_q(normalized).view(b, n, ipa.heads, ipa.scalar_dim)
        kv = ipa.linear_kv(normalized).view(b, n, ipa.heads, ipa.scalar_dim*2)
        k, v = torch.split(kv, ipa.scalar_dim, dim=-1)
        qlocal = ipa.linear_q_points(normalized).view(b, n, ipa.heads, ipa.qk_points, 3)
        kvlocal = ipa.linear_kv_points(normalized).view(b, n, ipa.heads, ipa.qk_points+ipa.v_points, 3)
        klocal, vlocal = torch.split(kvlocal, [ipa.qk_points, ipa.v_points], dim=-2)
        qpoints = local_to_global(qlocal, rotation, translation)
        kpoints = local_to_global(klocal, rotation, translation)
        vpoints = local_to_global(vlocal, rotation, translation)
        scale = torch.nn.functional.softplus(ipa.point_weights)
        scale = scale * math.sqrt(1. / (3.*(ipa.qk_points*9./2.)))
        output = torch.empty_like(single)
        for start in range(0, n, query_chunk):
            end = min(start+query_chunk, n)
            size = end-start
            pair = pair_rows(model.pair_encoder, batch, rotation, mask, start, end)
            scalar_logits = torch.einsum('bihc,bjhc->bhij', q[:, start:end], k)
            scalar_logits = scalar_logits * math.sqrt(1./(3.*ipa.scalar_dim))
            pair_logits = ipa.linear_pair_bias(pair).permute(0, 3, 1, 2) * math.sqrt(1./3.)
            difference = qpoints[:, start:end, None] - kpoints[:, None, :]
            d2 = torch.sum(difference*difference, dim=-1)
            point_logits = -.5*torch.sum(d2 * scale[None, None, None, :, None], dim=-1)
            point_logits = point_logits.permute(0, 3, 1, 2)
            logits = scalar_logits + pair_logits + point_logits
            valid = mask[:, None, start:end, None] & mask[:, None, None, :]
            attention = torch.softmax(logits.masked_fill(~valid, -1e5), dim=-1)
            scalar_output = torch.einsum('bhij,bjhc->bihc', attention, v)
            point_global = torch.einsum('bhij,bjhpc->bihpc', attention, vpoints)
            point_local = global_to_local(point_global, rotation[:, start:end], translation[:, start:end])
            point_norm = torch.sqrt(torch.sum(point_local*point_local, dim=-1)+ipa.eps)
            pair_output = torch.einsum('bhij,bijc->bihc', attention, pair)
            combined = torch.cat((scalar_output.reshape(b, size, -1),
                                  point_local.reshape(b, size, -1),
                                  point_norm.reshape(b, size, -1),
                                  pair_output.reshape(b, size, -1)), dim=-1)
            output[:, start:end] = ipa.linear_out(combined) * mask[:, start:end, None].to(single.dtype)
        single = single + block.ipa_dropout(output)
        single = single + block.transition_dropout(block.transition(block.transition_norm(single)))
        single = single * mask[..., None].to(single.dtype)
    representation = model.final_norm(single)
    mean, log_variance = model.state_head(representation)
    return {'representation': representation, 'mean': mean, 'log_variance': log_variance, 'mask': mask}
