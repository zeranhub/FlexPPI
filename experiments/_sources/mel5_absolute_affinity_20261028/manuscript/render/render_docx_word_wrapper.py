"""Use Word's read-only PDF export with the packaged DOCX PNG renderer.

The managed renderer remains unmodified. Its PDF-to-PNG, page numbering,
optional PDF persistence, and DPI calculation are used as supplied.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.stderr.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
RENDERER = Path(r"C:\Users\DELL\.codex\plugins\cache\openai-primary-runtime\documents\26.1025.11240\skills\documents\render_docx.py")
POPPLER = Path(r"C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin")
os.environ["PATH"] = str(POPPLER) + os.pathsep + os.environ.get("PATH", "")

spec = importlib.util.spec_from_file_location("codex_packaged_render_docx", RENDERER)
renderer = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(renderer)


def word_convert_to_pdf(doc_path, user_profile, convert_tmp_dir, stem, verbose):
    pdf_path = str(Path(convert_tmp_dir) / (stem + ".pdf"))
    command = [
        "powershell.exe", "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass",
        "-File", str(HERE / "export_word_readonly.ps1"),
        "-DocumentPath", str(doc_path), "-PdfPath", pdf_path,
    ]
    result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=600)
    log = f"Word COM exit: {result.returncode}\n{result.stdout}\n{result.stderr}"
    if verbose:
        print(log)
    if result.returncode != 0 or not Path(pdf_path).is_file() or Path(pdf_path).stat().st_size == 0:
        return "", log
    return pdf_path, log


renderer.convert_to_pdf = word_convert_to_pdf
if __name__ == "__main__":
    renderer.main()
