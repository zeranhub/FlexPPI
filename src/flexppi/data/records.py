"""JSON-lines manifest records.

Manifests deliberately store absolute paths. The two source datasets live on
different drives in the user's environment, and silently resolving a relative
path against the current working directory would be error-prone.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional


@dataclass
class StateRecord:
    sample_id: str
    dataset: str
    structure_path: str
    label_path: Optional[str] = None
    trajectory_paths: List[str] = field(default_factory=list)
    temperature: float = 300.0
    window_ns: float = 100.0
    target_chains: List[str] = field(default_factory=list)
    source_group: str = "unknown"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> "StateRecord":
        return cls(**value)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PairRecord:
    pair_id: str
    apo: StateRecord
    holo: StateRecord
    apo_indices: List[int]
    holo_indices: List[int]
    sequence_identity: float
    coverage: float
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, value: Dict[str, Any]) -> "PairRecord":
        copied = dict(value)
        copied["apo"] = StateRecord.from_dict(copied["apo"])
        copied["holo"] = StateRecord.from_dict(copied["holo"])
        return cls(**copied)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def read_jsonl(path: str) -> Iterator[Dict[str, Any]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                yield json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise ValueError("Invalid JSON at {}:{}".format(path, line_number)) from exc


def write_jsonl(path: str, records: Iterable[Dict[str, Any]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True))
            handle.write("\n")


def remap_path(path: str, mappings: Optional[Mapping[str, str]] = None) -> str:
    """Resolve manifests copied between machines without rewriting source files.

    A manifest produced on Linux can contain paths such as ``/research/data/...``
    while evaluation happens on Windows.  Prefix matching is performed on a
    slash-normalised representation, then the unmatched suffix is appended to
    the local target prefix with :class:`pathlib.Path`.
    """

    if not path or not mappings or Path(path).exists():
        return path
    normalised = path.replace("\\", "/").rstrip("/")
    candidates = sorted(mappings.items(), key=lambda item: len(item[0]), reverse=True)
    for source, target in candidates:
        source_normalised = str(source).replace("\\", "/").rstrip("/")
        if normalised == source_normalised:
            return str(Path(target))
        prefix = source_normalised + "/"
        if normalised.startswith(prefix):
            relative = normalised[len(prefix):]
            return str(Path(target).joinpath(*relative.split("/")))
    return path


def remap_state_record(
    record: StateRecord,
    mappings: Optional[Mapping[str, str]] = None,
) -> StateRecord:
    copied = StateRecord.from_dict(record.to_dict())
    copied.structure_path = remap_path(copied.structure_path, mappings)
    if copied.label_path:
        copied.label_path = remap_path(copied.label_path, mappings)
    copied.trajectory_paths = [remap_path(path, mappings) for path in copied.trajectory_paths]
    return copied


def remap_pair_record(
    record: PairRecord,
    mappings: Optional[Mapping[str, str]] = None,
) -> PairRecord:
    copied = PairRecord.from_dict(record.to_dict())
    copied.apo = remap_state_record(copied.apo, mappings)
    copied.holo = remap_state_record(copied.holo, mappings)
    return copied


def load_state_records(path: str) -> List[StateRecord]:
    return [StateRecord.from_dict(value) for value in read_jsonl(path)]


def load_pair_records(path: str) -> List[PairRecord]:
    return [PairRecord.from_dict(value) for value in read_jsonl(path)]
