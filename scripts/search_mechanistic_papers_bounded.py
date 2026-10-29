"""Run all six paper-search skill backends with individual wall-time limits.

The prior aggregate session was lost without retrievable output. Isolating the
backends preserves successful results if a different API stalls. No retries.
"""
import concurrent.futures
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts/apo_holo_mechanistic_literature_v2_20260913"
SKILL = Path("C:/Users/DELL/.codex/skills/paper_search/scripts")
SOURCES = ["semantic_scholar", "open_alex", "arxiv", "openreview", "crossref", "dblp"]
QUERY = "Hidden dynamic allostery"


def run_source(source):
    code = ("import sys,json; sys.path.insert(0," + repr(str(SKILL)) + "); "
            "from search_papers import search_papers; "
            "r=search_papers(query=" + repr(QUERY) + ",start_year=2000,end_year=2026,max_results=2,sources=[" + repr(source) + "]); "
            "print('RESULT_JSON='+json.dumps(r,ensure_ascii=True))")
    try:
        process = subprocess.run([sys.executable, "-u", "-c", code], capture_output=True, text=True,
                                 encoding="utf-8", errors="replace", timeout=55)
        return dict(source=source, returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
    except subprocess.TimeoutExpired as exc:
        def decode(value):
            return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else (value or "")
        return dict(source=source, returncode=None, stdout=decode(exc.stdout), stderr=decode(exc.stderr),
                    error=f"[{source}] outer wall-time limit exceeded (55 s); no complete result recovered")


if __name__ == "__main__":
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(run_source, SOURCES))
    OUT.mkdir(parents=True, exist_ok=True)
    record = dict(query=QUERY, years=[2000, 2026], max_per_source=2, results=results)
    (OUT / "multisource_search.json").write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(record, ensure_ascii=True, indent=2))
