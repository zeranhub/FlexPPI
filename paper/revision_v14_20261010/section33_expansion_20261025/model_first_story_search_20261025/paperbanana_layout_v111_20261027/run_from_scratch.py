"""Run the author's actual PaperBanana agents on custom method/caption files.

Credentials are read only into the process environment. They are never saved
into the output provenance, author result, or command-line arguments.
"""
from __future__ import annotations

import argparse
import asyncio
import base64
from datetime import datetime, timezone
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import re
import sys


TOOL_COMMIT = "f07f5f8a04d3025709d9ee3ec41d3d325a2e8529"


class GenerationFailure(RuntimeError):
    """Safe API diagnostic; deliberately excludes request/response payloads."""

    def __init__(self, *, model, code, reason, attempts):
        self.model = model
        self.code = code
        self.reason = reason
        self.attempts = attempts
        self.message = f"model={model}; reason={reason}; attempts={attempts}"
        super().__init__(self.message)


def failure_summary(error):
    return {
        "type": type(error).__name__,
        "status": getattr(error, "status_code", getattr(error, "code", "unavailable")),
        "reason": getattr(error, "reason", "pipeline_failed"),
        "model": getattr(error, "model", None),
        "attempts": getattr(error, "attempts", None),
    }


def is_zero_quota(error):
    """Read only quota fields, never headers, from a Gemini SDK exception."""
    message = str(getattr(error, "message", "") or "")
    if re.search(r"\blimit\s*[:=]?\s*0(?:\D|$)", message, re.IGNORECASE):
        return True

    def inspect(item):
        if isinstance(item, dict):
            for key, value in item.items():
                quota_key = re.sub(r"[^a-z]", "", str(key).lower())
                if quota_key in {"quotalimitvalue", "limit", "quotalimit", "quotavalue"}:
                    if value == 0 or value == "0":
                        return True
                if inspect(value):
                    return True
        elif isinstance(item, (list, tuple)):
            return any(inspect(value) for value in item)
        return False

    return inspect(getattr(error, "details", None))


def install_gemini_failfast(generation_utils, *, pause=asyncio.sleep):
    """Guard real requests instead of making an extra paid preflight request.

    This changes only the runner's retry/error behavior. Author agents, their
    prompts, the SDK request conversion, and the image pipeline remain intact.
    """
    from google.genai import types

    async def guarded_call(model_name, contents, config, max_attempts=5,
                           retry_delay=5, error_context=""):
        if generation_utils.gemini_client is None:
            raise GenerationFailure(model=model_name, code=401,
                                    reason="missing_google_api_key", attempts=0)
        target_count = int(config.candidate_count or 1)
        request_config = config.model_copy(deep=True)
        request_config.candidate_count = min(target_count, 8)
        request_config.http_options = (request_config.http_options.model_copy(deep=True)
                                       if request_config.http_options else types.HttpOptions())
        # One SDK transport attempt per wrapper attempt; avoid nested retries.
        request_config.http_options.retry_options = types.HttpRetryOptions(attempts=1)
        converted = generation_utils._convert_to_gemini_parts(contents)
        for attempt in range(1, 3):
            try:
                response = await generation_utils.gemini_client.aio.models.generate_content(
                    model=model_name, contents=converted, config=request_config
                )
            except Exception as error:
                code = getattr(error, "status_code", getattr(error, "code", None))
                try:
                    code = int(code)
                except (ValueError, TypeError):
                    code = "unavailable"
                zero_quota = code == 429 and is_zero_quota(error)
                if zero_quota:
                    reason = "quota_limit_zero_enable_billing_or_use_eligible_model"
                elif code in (401, 403):
                    reason = "authentication_or_permission_denied"
                elif code == 404:
                    reason = "model_not_found_or_unavailable"
                elif code == 429:
                    reason = "rate_or_quota_limit_exhausted"
                else:
                    reason = "api_request_failed"
                if code == 429 and not zero_quota and attempt == 1:
                    print(f"API retry: model={model_name}; status=429; reason=rate_limit; retry=1/1", flush=True)
                    await pause(min(max(float(retry_delay), 0), 5))
                    continue
                print(f"API stopped: model={model_name}; status={code}; reason={reason}; attempts={attempt}", flush=True)
                raise GenerationFailure(model=model_name, code=code,
                                        reason=reason, attempts=attempt) from None

            parts = [part for candidate in (response.candidates or [])
                     for part in (getattr(candidate.content, "parts", None) or [])
                     if candidate.content is not None]
            image_mode = "image" in model_name or "nanoviz" in model_name
            if image_mode:
                values = [base64.b64encode(part.inline_data.data).decode("ascii")
                          for part in parts
                          if getattr(part, "inline_data", None) is not None
                          and part.inline_data.data]
            else:
                values = [part.text for part in parts
                          if isinstance(getattr(part, "text", None), str) and part.text.strip()]
            if len(values) < target_count or any(value.strip() == "Error" for value in values):
                raise GenerationFailure(model=model_name, code="empty_response",
                                        reason="no_valid_image" if image_mode else "no_valid_text",
                                        attempts=attempt)
            return values[:target_count]

    generation_utils.call_gemini_with_retry_async = guarded_call


class RedactedStream:
    def __init__(self, stream, secret):
        self.stream = stream
        self.secret = secret
        self.pending = ""

    def write(self, message):
        original_length = len(message)
        self.pending += message.replace(self.secret, "[REDACTED]") if self.secret else message
        while "\n" in self.pending:
            line, self.pending = self.pending.split("\n", 1)
            # Author modules sometimes print raw model/error responses. Only emit
            # our stage/output summaries, never prompts, base64 or HTTP headers.
            allowed = ("Stage:", "Validated author pipeline:", "Dry run passed.",
                       "Reference benchmark is absent;", "Saved candidate ",
                       "Author PaperBanana pipeline finished.", "Run failed:",
                       "API retry:", "API stopped:", "Saved diagnostics:", "Saved checkpoint:")
            if line.startswith(allowed):
                self.stream.write(line[:500] + "\n")
        return original_length

    def flush(self):
        return self.stream.flush()

    def __getattr__(self, name):
        return getattr(self.stream, name)


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", required=True, type=Path)
    parser.add_argument("--caption", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repo", "--paperbanana-root", dest="repo", type=Path, default=Path(__file__).parent / "PaperBanana")
    parser.add_argument("--key-file", type=Path)
    parser.add_argument("--main-model", default="gemini-3.1-pro-preview")
    parser.add_argument("--image-model", default="gemini-3.1-flash-image-preview")
    parser.add_argument("--aspect-ratio", default="16:9")
    parser.add_argument("--resolution", choices=["1K", "2K", "4K"], default="2K")
    parser.add_argument("--critic-rounds", type=int, default=2)
    parser.add_argument("--candidates", type=int, default=1)
    parser.add_argument("--retrieval", choices=["auto", "none"], default="auto")
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and imports without API calls.")
    return parser.parse_args()


def save_result(result, output, candidate_number, *, error=None, require_final=True):
    from PIL import Image

    metadata = {}
    image_files = {}
    for key, value in result.items():
        if "base64" in key and isinstance(value, str) and value:
            raw = value.split(",", 1)[1] if value.startswith("data:") else value
            try:
                image = Image.open(BytesIO(base64.b64decode(raw)))
                image.load()
            except Exception:
                metadata[key] = "[Undecodable image omitted]"
                continue
            safe_key = re.sub(r"[^A-Za-z0-9_.-]", "_", key)
            path = output / f"candidate_{candidate_number}_{safe_key}.png"
            image.save(path)
            image_files[key] = path
            metadata[key] = {"file": path.name, "size": list(image.size)}
        else:
            metadata[key] = value
    final_key = result.get("eval_image_field", "")
    final_path = image_files.get(final_key)
    if error is not None:
        metadata["pipeline_error"] = failure_summary(error)
    if final_path is None:
        metadata.setdefault("pipeline_error", {
            "type": "MissingFinalImage", "reason": "author_pipeline_produced_no_final_image",
            "requested_field": final_key,
        })
    # Write diagnostics before checking final output, so failed runs retain the
    # useful descriptions/stages without writing any image base64 or credentials.
    metadata_path = output / f"candidate_{candidate_number}_pipeline.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    if final_path is None:
        if not require_final:
            return
        raise RuntimeError(f"Author pipeline produced no final image; requested field: {final_key!r}")
    final_image = Image.open(final_path)
    final_image.save(output / f"candidate_{candidate_number}_final.png")
    if candidate_number == 1:
        final_image.save(output / "framework_final.png")
    print(f"Saved candidate {candidate_number}: {final_path} ({final_image.width}x{final_image.height})", flush=True)


def install_visualizer_checkpoints(visualizer, output):
    """Save returned intermediate images before the next remote review starts."""
    from PIL import Image

    original_process = visualizer.process

    def without_base64(value):
        if isinstance(value, dict):
            return {key: without_base64(item) for key, item in value.items()
                    if "base64" not in str(key).lower()}
        if isinstance(value, (list, tuple)):
            return [without_base64(item) for item in value]
        return value

    async def checkpointed_process(data):
        result = await original_process(data)
        root = output / "checkpoints"
        root.mkdir(parents=True, exist_ok=True)
        candidate = re.sub(r"[^A-Za-z0-9_.-]", "_", str(result.get("candidate_id", 0)))
        pass_number = 1
        while True:
            folder = root / f"candidate_{candidate}_visualizer_{pass_number:03d}"
            try:
                folder.mkdir()
                break
            except FileExistsError:
                pass_number += 1
        images = {}
        for key, value in result.items():
            if "base64" not in key.lower() or not isinstance(value, str) or not value:
                continue
            raw = value.split(",", 1)[1] if value.startswith("data:") else value
            try:
                image = Image.open(BytesIO(base64.b64decode(raw)))
                image.load()
            except Exception:
                continue
            filename = re.sub(r"[^A-Za-z0-9_.-]", "_", key) + ".png"
            image.save(folder / filename)
            images[key] = {"file": filename, "size": list(image.size)}
        metadata = {
            "checkpoint_status": "visualizer_complete_pending_review",
            "images": images,
            "pipeline": without_base64(result),
        }
        (folder / "pipeline.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
        )
        print(f"Saved checkpoint: {folder}; images={len(images)}; pending review", flush=True)
        return result

    visualizer.process = checkpointed_process


async def run(args):
    # Imports resolve to the unmodified author source snapshot.
    from agents.critic_agent import CriticAgent
    from agents.planner_agent import PlannerAgent
    from agents.polish_agent import PolishAgent
    from agents.retriever_agent import RetrieverAgent
    from agents.stylist_agent import StylistAgent
    from agents.vanilla_agent import VanillaAgent
    from agents.visualizer_agent import VisualizerAgent
    from utils.config import ExpConfig
    from utils.legacy_generation_options import generation_additional_info
    from utils.paperviz_processor import PaperVizProcessor
    from utils import generation_utils

    install_gemini_failfast(generation_utils)

    method = args.method.read_text(encoding="utf-8-sig").strip()
    caption = args.caption.read_text(encoding="utf-8-sig").strip()
    if not method or not caption:
        raise ValueError("Both method and caption files must contain text.")
    if args.candidates < 1 or args.critic_rounds < 0:
        raise ValueError("Candidates must be positive and critic rounds cannot be negative.")
    reference_file = Path(__file__).resolve().parent / "data" / "PaperBananaBench" / "diagram" / "ref.json"
    retrieval = args.retrieval
    if retrieval == "auto" and not reference_file.is_file():
        retrieval = "none"
        print("Reference benchmark is absent; using the author's supported retrieval='none' mode.", flush=True)

    config = ExpConfig(
        dataset_name="PaperBananaBench", task_name="diagram", split_name="demo",
        exp_mode="demo_full", retrieval_setting=retrieval,
        max_critic_rounds=args.critic_rounds,
        main_model_name=args.main_model, image_gen_model_name=args.image_model,
        work_dir=Path(__file__).resolve().parent,
    )
    processor = PaperVizProcessor(
        exp_config=config,
        vanilla_agent=VanillaAgent(exp_config=config),
        planner_agent=PlannerAgent(exp_config=config),
        visualizer_agent=VisualizerAgent(exp_config=config),
        stylist_agent=StylistAgent(exp_config=config),
        critic_agent=CriticAgent(exp_config=config),
        retriever_agent=RetrieverAgent(exp_config=config),
        polish_agent=PolishAgent(exp_config=config),
    )
    original_planner_process = processor.planner_agent.process
    async def planner_with_manuscript_style(data):
        refs = json.loads((Path(__file__).resolve().parent / "data/PaperBananaBench/diagram/ref.json").read_text(encoding="utf-8"))
        by_id = {item["id"]: item for item in refs}
        examples = [by_id[item] for item in data.get("top10_references", []) if item in by_id]
        style_example = {"id":"manuscript_Figure2_style_only", "content":"STYLE ONLY: Existing Figure2 from the same manuscript. Match its compact quantitative visual language, Arial typography, blue/ochre data colors, delicate axes and panel labels; its scientific content and layout are not the target.", "visual_intent":"Use this visual grammar for a NEW protein case-results composition. Never reproduce this figure's plots or duplicate its layout.", "path_to_gt_image":str(Path(__file__).resolve().parent / "manuscript_figure2_style_reference.png")}
        data["retrieved_examples"] = [style_example] + examples
        return await original_planner_process(data)
    processor.planner_agent.process = planner_with_manuscript_style
    install_visualizer_checkpoints(processor.visualizer_agent, args.output)
    records = [dict(
        filename=f"flexppi_framework_candidate_{index}", candidate_id=index,
        caption=caption, content=method, visual_intent=caption + ("\nIndependent composition A: invent an elegant new data-first composition organized for side-by-side comparison of the two complete cases. The two structures stay adjacent, the two trace+SASA blocks match sizes, the four chi plots stay one full row, and the density stays compact. No fixed inherited grid." if index == 0 else "\nIndependent composition B: explore a genuinely different page-wide hierarchy from A, optimizing a coherent integrated case-evidence region, adjacent substantial structures, matched trace+SASA blocks, a compact global density, and one full row of four chi plots. Eliminate zigzag reading and inherited panel positions."),
        additional_info=generation_additional_info(args.aspect_ratio, args.resolution),
        max_critic_rounds=args.critic_rounds, task_name="diagram",
    ) for index in range(args.candidates)]
    print(f"Validated author pipeline: demo_full; retrieval={retrieval}; candidates={args.candidates}; critic_rounds={args.critic_rounds}; resolution={args.resolution}", flush=True)
    if args.dry_run:
        print("Dry run passed. No generation or paid API request was made.", flush=True)
        return

    args.output.mkdir(parents=True, exist_ok=True)
    provenance = dict(
        source="https://github.com/dwzhu-pku/PaperBanana",
        old_figure_supplied=False, local_style_override=True, manuscript_figure2_style_image=True, author_agents_unmodified=True, source_commit=TOOL_COMMIT, entrypoint="utils.paperviz_processor.PaperVizProcessor.process_queries_batch",
        acquisition="GitHub official source archive (git clone transport reset)",
        started_at=datetime.now(timezone.utc).isoformat(),
        method_sha256=sha256(method.encode()).hexdigest(), caption_sha256=sha256(caption.encode()).hexdigest(),
        main_model=args.main_model, image_model=args.image_model,
        aspect_ratio=args.aspect_ratio, resolution=args.resolution,
        exp_mode="demo_full", retrieval=retrieval, candidates=args.candidates,
        critic_rounds=args.critic_rounds,
        request_guard="No extra preflight request; quota-zero/auth/model errors stop immediately; other 429 retries once; SDK retries disabled.",
    )
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    previous_stage = None

    def progress(status):
        nonlocal previous_stage
        stage = status.get("stage") or status.get("active_stage", "")
        if stage and stage != previous_stage:
            print(f"Stage: {stage}", flush=True)
            previous_stage = stage

    completed = 0
    try:
        async for result in processor.process_queries_batch(records, max_concurrent=2, do_eval=False, progress_callback=progress):
            completed += 1
            save_result(result, args.output, completed)
    except Exception as error:
        diagnostic = failure_summary(error)
        diagnostic["stage"] = previous_stage
        (args.output / "failure.json").write_text(json.dumps(diagnostic, indent=2), encoding="utf-8")
        for index, partial_record in enumerate(records, 1):
            save_result(partial_record, args.output, index, error=error, require_final=False)
        print(f"Saved diagnostics: {args.output / 'failure.json'}", flush=True)
        raise
    if completed != args.candidates:
        raise RuntimeError(f"Only {completed} of {args.candidates} candidates completed.")
    print(f"Author PaperBanana pipeline finished. Final figure: {args.output / 'framework_final.png'}", flush=True)


def main():
    args = arguments()
    args.repo = args.repo.resolve()
    args.method = args.method.resolve()
    args.caption = args.caption.resolve()
    args.output = args.output.resolve()
    if not (args.repo / "utils" / "paperviz_processor.py").is_file():
        raise SystemExit("Author source repository is absent.")
    secret = ""
    if not args.dry_run:
        secret = args.key_file.read_text(encoding="utf-8-sig").strip() if args.key_file else os.environ.get("GOOGLE_API_KEY", "").strip()
        if not secret:
            raise SystemExit("Set GOOGLE_API_KEY or provide --key-file; no API request was made.")
        os.environ["GOOGLE_API_KEY"] = secret
    sys.stdout = RedactedStream(sys.stdout, secret)
    sys.stderr = RedactedStream(sys.stderr, secret)
    sys.path.insert(0, str(args.repo))
    os.chdir(args.repo)
    asyncio.run(run(args))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        code = getattr(error, "status_code", getattr(error, "code", "unavailable"))
        fallback = str(error) if type(error).__module__ == "builtins" else ""
        reason = str(getattr(error, "message", fallback))[:240].replace("\n", " ")
        print(f"Run failed: {type(error).__name__}; status={code}; reason={reason}", flush=True)
        raise SystemExit(1)
