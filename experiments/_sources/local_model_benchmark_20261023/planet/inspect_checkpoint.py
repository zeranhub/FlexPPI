"""Read released torch-checkpoint hyperparameters without allocating tensors.

This constrained unpickler reads only metadata in the torch ZIP archive. It
replaces tensor rebuild operations with shape descriptors and rejects all
unrecognised globals. It does not import torch or read tensor-storage payloads.
"""
import argparse
from collections import OrderedDict
import io
import json
import math
from pathlib import Path
import pickle
import struct
import zipfile


class TensorSpec:
    def __init__(self, storage, offset, shape, stride, *args):
        self.storage = storage
        self.offset = offset
        self.shape = tuple(shape)
        self.stride = tuple(stride)


class MetadataUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) == ("collections", "OrderedDict"):
            return OrderedDict
        if module == "torch" and name.endswith("Storage"):
            return name
        if module == "torch._utils" and name in {"_rebuild_tensor", "_rebuild_tensor_v2", "_rebuild_tensor_v3"}:
            return TensorSpec
        if module == "torch._utils" and name == "_rebuild_parameter":
            return lambda data, *args: data
        raise pickle.UnpicklingError(f"Unsupported metadata global {module}.{name}")

    def persistent_load(self, value):
        assert value[0] == "storage"
        return {"type": value[1], "key": value[2], "location": value[3], "elements": value[4]}


def inspect(path):
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        metadata_name = next(name for name in archive.namelist() if name.endswith("/data.pkl"))
        checkpoint = MetadataUnpickler(io.BytesIO(archive.read(metadata_name))).load()
        default_temp = checkpoint["state_dict"]["child_model.default_temp"]
        storage_name = metadata_name[:-len("data.pkl")] + "data/" + default_temp.storage["key"]
        raw_temp = archive.read(storage_name)
        assert default_temp.storage["type"] == "FloatStorage"
        default_temp_value = struct.unpack_from("<f", raw_temp, 4 * default_temp.offset)[0]
    state = checkpoint["state_dict"]
    tensors = {key: value for key, value in state.items() if isinstance(value, TensorSpec)}
    counts = {}
    for value in tensors.values():
        dtype = value.storage["type"]
        counts[dtype] = counts.get(dtype, 0) + math.prod(value.shape)
    return {
        "path": str(path), "file_bytes": path.stat().st_size,
        "hyper_parameters": checkpoint["hyper_parameters"],
        "state_tensor_count": len(tensors), "state_elements_by_storage_dtype": counts,
        "default_temp_shape": list(tensors["child_model.default_temp"].shape),
        "default_temp_K": default_temp_value,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--output")
    args = parser.parse_args()
    result = [inspect(path) for path in args.paths]
    text = json.dumps(result, indent=2)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    print(text)
