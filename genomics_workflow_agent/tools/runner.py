from __future__ import annotations

import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

STATUS_PLANNED = "planned"
STATUS_SKIPPED = "skipped"
STATUS_SUCCEEDED = "succeeded"
STATUS_FAILED = "failed"
STATUS_ERROR = "error"


def validate_outputs(expected: list[str]) -> dict[str, Any]:
    present = [p for p in expected if Path(p).exists()]
    missing = [p for p in expected if not Path(p).exists()]
    return {
        "expected": expected,
        "present": present,
        "missing": missing,
        "all_present": len(missing) == 0,
    }


def run_command(
    cmd: list[str],
    *,
    cwd: str | Path | None = None,
    dry_run: bool = True,
    timeout: int = 3600,
    provenance_dir: Path | None = None,
    label: str = "",
    expected_outputs: list[str] | None = None,
    capture_stdout_path: str | Path | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "label": label,
        "command": cmd,
        "command_str": " ".join(str(c) for c in cmd),
        "cwd": str(cwd) if cwd else str(Path.cwd()),
        "dry_run": dry_run,
        "started_at": None,
        "ended_at": None,
        "return_code": None,
        "stdout_snippet": None,
        "stderr_snippet": None,
        "error": None,
        "executed": False,
        "status": STATUS_PLANNED,
        "output_validation": None,
    }

    if dry_run:
        record["started_at"] = datetime.now(timezone.utc).isoformat()
        record["ended_at"] = record["started_at"]
        record["status"] = STATUS_PLANNED
        if expected_outputs:
            record["output_validation"] = {"expected": expected_outputs, "note": "not validated in dry-run"}
        _write_provenance(record, provenance_dir, label)
        return record

    record["started_at"] = datetime.now(timezone.utc).isoformat()
    t0 = time.monotonic()
    try:
        result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        elapsed = round(time.monotonic() - t0, 2)
        record["ended_at"] = datetime.now(timezone.utc).isoformat()
        record["return_code"] = result.returncode
        record["stdout_snippet"] = (result.stdout or "")[-2000:] or None
        record["stderr_snippet"] = (result.stderr or "")[-2000:] or None
        record["runtime_s"] = elapsed
        record["executed"] = True

        if capture_stdout_path and result.stdout:
            cap_path = Path(capture_stdout_path)
            cap_path.parent.mkdir(parents=True, exist_ok=True)
            cap_path.write_text(result.stdout, encoding="utf-8")
            record["stdout_file"] = str(cap_path)

        if result.returncode != 0:
            record["error"] = f"Non-zero exit code: {result.returncode}"
            record["status"] = STATUS_FAILED
        else:
            if expected_outputs:
                validation = validate_outputs(expected_outputs)
                record["output_validation"] = validation
                record["status"] = STATUS_SUCCEEDED if validation["all_present"] else STATUS_FAILED
                if not validation["all_present"]:
                    record["error"] = f"Expected outputs missing: {validation['missing']}"
            else:
                record["status"] = STATUS_SUCCEEDED

    except subprocess.TimeoutExpired:
        record["ended_at"] = datetime.now(timezone.utc).isoformat()
        record["error"] = f"Command timed out after {timeout}s"
        record["executed"] = True
        record["status"] = STATUS_ERROR
    except FileNotFoundError:
        record["ended_at"] = datetime.now(timezone.utc).isoformat()
        record["error"] = f"Executable not found: {cmd[0]}"
        record["status"] = STATUS_ERROR
    except Exception as e:
        record["ended_at"] = datetime.now(timezone.utc).isoformat()
        record["error"] = str(e)
        record["status"] = STATUS_ERROR

    _write_provenance(record, provenance_dir, label)
    return record


def _write_provenance(record: dict, provenance_dir: Path | None, label: str) -> None:
    if provenance_dir is None:
        return
    provenance_dir = Path(provenance_dir)
    provenance_dir.mkdir(parents=True, exist_ok=True)
    safe_label = (label or "cmd").replace(" ", "_").replace("/", "_")[:60]
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
    out_path = provenance_dir / f"provenance_{safe_label}_{ts}.json"
    out_path.write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")


def assert_no_silent_failure(record: dict) -> None:
    if record.get("executed") and record.get("return_code") not in (None, 0):
        raise RuntimeError(
            f"Command failed (exit {record['return_code']}): {record['command_str']}\n"
            f"stderr: {record.get('stderr_snippet', '')}"
        )


def execution_summary(step_results: list[dict]) -> dict[str, Any]:
    counts = {STATUS_PLANNED: 0, STATUS_SKIPPED: 0, STATUS_SUCCEEDED: 0,
              STATUS_FAILED: 0, STATUS_ERROR: 0}
    for r in step_results:
        status = r.get("status", STATUS_PLANNED)
        counts[status] = counts.get(status, 0) + 1
    failed = [r for r in step_results if r.get("status") in (STATUS_FAILED, STATUS_ERROR)]
    return {
        "total": len(step_results),
        "counts": counts,
        "failed_steps": [r.get("label") for r in failed],
        "overall_status": "failed" if failed else ("planned" if counts[STATUS_PLANNED] > 0 else "succeeded"),
    }
