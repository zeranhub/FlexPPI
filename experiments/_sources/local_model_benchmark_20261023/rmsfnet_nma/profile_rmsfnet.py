"""Low-memory CPU layout microbenchmark; original model/weights remain unmodified."""
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key] = '2'
sys.dont_write_bytecode = True
import numpy as np
from scipy.stats import spearmanr
from cohort_common import ROOT, COHORT, OLD, sha, load_records, local_input, record_key, write_json


def load_adapter():
    sys.path.insert(0, 'E:/Project/FlexPPI/tmp/third_party/RMSF-net')
    spec = importlib.util.spec_from_file_location('frozen_adapter', OLD/'rmsf_net_pdb_adapter.py')
    adapter = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = adapter
    spec.loader.exec_module(adapter)
    adapter.torch.set_num_threads(2)
    adapter.torch.set_num_interop_threads(1)
    return adapter


def main():
    adapter = load_adapter()
    torch = adapter.torch
    rows, _ = load_records(True)
    output = ROOT/'rmsfnet/profile'
    output.mkdir(parents=True, exist_ok=True)
    weights = Path('E:/Project/FlexPPI/tmp/third_party/RMSF-net/model_cpu/sim_model.pth')
    model = adapter.load_model(weights)
    atoms, residues = adapter.parse_pdb(str(local_input(rows[0])))
    origin, shape_xyz = adapter.grid_geometry(atoms)
    density = adapter.simulated_density(atoms, origin, shape_xyz)
    indices = adapter.atom_indices_zyx(atoms, origin, density.shape)
    cells = adapter.selected_cells(indices)
    counts = np.ceil(np.asarray(density.shape)/adapter.CORE_SIZE).astype(int)
    back = counts*adapter.CORE_SIZE-np.asarray(density.shape)+adapter.BOX_HALO
    padded = np.pad(density, tuple((adapter.BOX_HALO,int(v)) for v in back))
    # Four patches are enough to inspect batch four without a large complete assembly allocation.
    patches = np.stack([adapter.extract_patch(padded, cell) for cell in cells[:4]])[:,None]
    results = dict(model_weights_sha256=sha(weights), dtype='float32', threads=2,
        mkldnn_enabled=torch.backends.mkldnn.enabled, torch_version=torch.__version__,
        script_sha256=sha(__file__), input_patch_sample=rows[0]['sample_id'], microbenchmarks=[])
    reference = None
    with torch.inference_mode():
        for layout in ['contiguous','channels_last_3d']:
            memory_format = torch.contiguous_format if layout == 'contiguous' else torch.channels_last_3d
            model.to(memory_format=memory_format)
            for batch_size in [1,2,4]:
                tensor = torch.from_numpy(patches[:batch_size]).contiguous(memory_format=memory_format)
                tick = time.perf_counter()
                value = model(tensor).cpu().numpy().copy()
                cold = time.perf_counter()-tick
                timings=[]
                for repeat in range(2):
                    tick = time.perf_counter()
                    value = model(tensor).cpu().numpy().copy()
                    timings.append(time.perf_counter()-tick)
                if reference is None:
                    reference = value[0].copy()
                item = dict(layout=layout, batch_size=batch_size, cold_seconds=cold,
                    repeat_seconds=timings, seconds_per_patch=float(np.median(timings)/batch_size),
                    first_patch_log10_maxabs=float(np.max(np.abs(value[0]-reference))))
                results['microbenchmarks'].append(item)
                write_json(output/'microbenchmarks.json', results)
                print(json.dumps(item), flush=True)


if __name__ == '__main__':
    main()
