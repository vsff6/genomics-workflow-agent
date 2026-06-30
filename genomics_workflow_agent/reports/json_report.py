"""JSON report writer with provenance."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def write_json_report(
    data: dict[str, Any],
    output_path: str | Path,
    *,
    tool: str = "genomics_workflow_agent",
    version: str = "0.3.0",
) -> Path:
    """Write structured JSON report, injecting metadata."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    report = {
        "meta": {
            "tool": tool,
            "version": version,
            "generated_at": datetime.now(timezone.utc).isoformat(),
        },
        **data,
    }

    output_path.write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    return output_path
