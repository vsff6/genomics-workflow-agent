"""Markdown report generator for workflow plans and inspection results."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def write_markdown_report(
    plan: dict[str, Any],
    output_path: str | Path,
    *,
    title: str | None = None,
) -> Path:
    """Generate a Markdown report from a workflow plan dict."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    workflow = plan.get("workflow", "unknown")
    report_title = title or f"Genomics Workflow Report: {workflow}"
    lines = [
        f"# {report_title}",
        f"\n**Generated**: {datetime.now(timezone.utc).isoformat()}",
        f"\n**Workflow**: `{workflow}`",
        f"\n**Input**: `{plan.get('input_dir', 'N/A')}`",
        f"\n**Output**: `{plan.get('output_dir', 'N/A')}`",
        f"\n**Mode**: {'DRY RUN (no commands executed)' if plan.get('dry_run') else 'EXECUTE'}",
        "\n---\n",
    ]

    # Input summary
    if "samples_detected" in plan:
        lines.append(f"## Input Summary\n\n- Samples detected: **{plan['samples_detected']}**")
    if "bam_files_found" in plan:
        lines.append(f"- BAM files: **{len(plan['bam_files_found'])}**")
    if "vcf_files_found" in plan:
        lines.append(f"- VCF files: **{len(plan['vcf_files_found'])}**")

    # Tool status
    if "tool_status" in plan:
        lines.append("\n## Tool Availability\n")
        for tool, info in plan["tool_status"].items():
            icon = "✓" if info["available"] else "✗"
            ver = info.get("version") or "not found"
            lines.append(f"- {icon} **{tool}**: {ver}")

    # Blockers
    blockers = plan.get("blockers", [])
    if blockers:
        lines.append("\n## Blockers\n")
        lines.append("> Resolve these before executing the workflow.\n")
        for b in blockers:
            lines.append(f"- **[BLOCKER]** {b}")

    # Warnings
    warnings = plan.get("warnings", [])
    if warnings:
        lines.append("\n## Warnings\n")
        for w in warnings:
            lines.append(f"- [WARN] {w}")

    # Samplesheet
    ss = plan.get("samplesheet")
    if ss:
        lines.append("\n## Samplesheet\n")
        if ss.get("created"):
            lines.append(f"- **Path**: `{ss['path']}`")
            lines.append(f"- **Rows**: {ss['rows']}")
            for w in ss.get("warnings", []):
                lines.append(f"- [WARN] {w}")
        else:
            lines.append(f"- Not created: {ss.get('reason', 'unknown')}")

    # Command
    cmd_str = plan.get("command_str") or (
        " ".join(plan["command"]) if isinstance(plan.get("command"), list) else None
    )
    if cmd_str:
        lines.append("\n## Planned Command\n")
        lines.append("```bash")
        lines.append(cmd_str)
        lines.append("```\n")
        if plan.get("dry_run"):
            lines.append("> **This command was NOT executed** (dry-run mode).\n")

    # Steps
    steps = plan.get("steps", [])
    if steps:
        lines.append(f"\n## Workflow Steps ({len(steps)} planned)\n")
        for i, step in enumerate(steps, 1):
            name = step.get("name", f"step_{i}")
            desc = step.get("description", "")
            lines.append(f"### Step {i}: {name}\n")
            if desc:
                lines.append(f"{desc}\n")
            cmd = step.get("command")
            if cmd:
                lines.append("```bash")
                lines.append(" ".join(str(c) for c in cmd))
                lines.append("```")
            if step.get("note"):
                lines.append(f"\n> {step['note']}")
            if step.get("skipped"):
                lines.append(f"\n> **SKIPPED**: {step.get('reason', '')}")
            lines.append("")

    # Skipped steps
    skipped = plan.get("skipped_steps", [])
    if skipped:
        lines.append("\n## Skipped Steps\n")
        for s in skipped:
            lines.append(f"- **{s.get('step', '?')}**: {s.get('reason', '')} — Install: `{s.get('install', 'see docs')}`")

    # Biological caveats
    caveats = plan.get("biological_caveats", [])
    if caveats:
        lines.append("\n## Biological Caveats\n")
        lines.append(
            "> These apply regardless of whether the workflow completes successfully. "
            "Successful execution is not biological or clinical validation.\n"
        )
        for c in caveats:
            lines.append(f"- {c}")

    # Clinical disclaimer
    disclaimer = plan.get("clinical_disclaimer")
    if disclaimer:
        lines.append(f"\n## Clinical Disclaimer\n\n> {disclaimer}\n")

    # Next actions
    next_actions = plan.get("next_actions", [])
    if next_actions:
        lines.append("\n## Recommended Next Actions\n")
        for a in next_actions:
            lines.append(f"1. {a}")

    # Limitations
    limitations = plan.get("limitations", [])
    if limitations:
        lines.append("\n## Limitations\n")
        for lim in limitations:
            lines.append(f"- {lim}")

    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path


def write_inspection_report(inspection: dict[str, Any], output_path: str | Path) -> Path:
    """Generate a Markdown report from a directory inspection result."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# File Inspection Report",
        f"\n**Path**: `{inspection.get('input_path', 'N/A')}`",
        f"\n**Inspected at**: {inspection.get('inspected_at', 'N/A')}",
        f"\n**Total files**: {inspection.get('total_files', 0)}",
        f"\n**Total size**: {inspection.get('total_size_mb', 0)} MB",
        f"\n**Workflow guess**: `{inspection.get('workflow_guess', 'unknown')}`",
        "\n---\n",
        "## File Type Summary\n",
    ]

    for ftype, count in (inspection.get("file_type_counts") or {}).items():
        lines.append(f"- **{ftype}**: {count}")

    files = inspection.get("files", [])
    if files:
        lines.append("\n## Files\n")
        lines.append("| Name | Type | Size (MB) |")
        lines.append("|------|------|-----------|")
        for f in files[:50]:
            lines.append(f"| `{f['name']}` | {f['type']} | {f['size_mb']} |")
        if len(files) > 50:
            lines.append(f"\n*... and {len(files) - 50} more files*")

    warnings = inspection.get("warnings", [])
    if warnings:
        lines.append("\n## Warnings\n")
        for w in warnings:
            lines.append(f"- {w}")

    output_path.write_text("\n".join(lines), encoding="utf-8")
    return output_path
