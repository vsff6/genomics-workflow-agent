"""Render an InterpretationResult dict to Markdown."""
from __future__ import annotations

from typing import Any


def render_interpretation_md(interp: dict[str, Any]) -> str:
    """Render a JSON interpretation result dict to a Markdown string."""
    lines: list[str] = [
        "## Biological interpretation and hypotheses",
        "",
        f"> **Interpretation version**: {interp.get('interpretation_version', '?')} "
        f"| **Workflow**: `{interp.get('workflow', '?')}`",
        "",
        f"**Scope**: {interp.get('scope', '')}",
        "",
        "> **Note**: This section contains deterministic scaffolds for human review. "
        "No LLM was called to produce these interpretations. "
        "Findings and hypotheses are rule-based and require expert review before any action is taken. "
        "**Clinical claims are not made here.**",
        "",
    ]

    limitations = interp.get("limitations", [])
    if limitations:
        lines.append("### Limitations")
        for lim in limitations:
            lines.append(f"- {lim}")
        lines.append("")

    findings = interp.get("findings", [])
    if findings:
        lines.append(f"### Findings ({len(findings)})")
        lines.append("")
        for finding in findings:
            lines.append(f"#### {finding['finding_id']}: {finding['sample']}")
            lines.append("")
            lines.append(f"**Observation**: {finding['observation']}")
            lines.append(f"**Evidence source**: `{finding['evidence_source']}`")
            lines.append("")

            tech = finding.get("technical_explanations", [])
            if tech:
                lines.append("**Technical explanations**:")
                for t in tech:
                    lines.append(f"- {t}")
                lines.append("")

            bio = finding.get("plausible_biological_explanations", [])
            if bio:
                lines.append("**Plausible biological explanations**:")
                for b in bio:
                    lines.append(f"- {b}")
                lines.append("")

            meta = finding.get("metadata_needed", [])
            if meta:
                lines.append("**Metadata needed to distinguish artifact from biology**:")
                for m in meta:
                    lines.append(f"- {m}")
                lines.append("")

            val = finding.get("recommended_validation", [])
            if val:
                lines.append("**Recommended validation**:")
                for v in val:
                    lines.append(f"- {v}")
                lines.append("")

            lines.append(f"**Recommended action**: {finding.get('recommended_action', '')}")
            lines.append(f"**Confidence**: {finding.get('confidence', '?')}")
            lines.append(
                f"**Should filter**: {'Yes' if finding.get('should_filter') else 'No'} | "
                f"**Preserve until review**: {'Yes' if finding.get('should_preserve_until_review') else 'No'}"
            )
            lines.append("")
            lines.append("> **Safety note**: No clinical claims are made. "
                         "These are research-grade hypotheses for expert review.")
            lines.append("")

    hypotheses = interp.get("hypotheses", [])
    if hypotheses:
        lines.append(f"### Hypotheses ({len(hypotheses)})")
        lines.append("")
        for hyp in hypotheses:
            lines.append(f"#### {hyp['hypothesis_id']}")
            lines.append("")
            lines.append(f"**Statement**: {hyp['statement']}")
            lines.append(f"**Interpretation type**: {hyp.get('interpretation_type', '?')} "
                         f"| **Confidence**: {hyp.get('confidence', '?')}")
            lines.append(f"**Clinical claim**: {'Yes' if hyp.get('clinical_claim') else 'No (research only)'}")
            lines.append("")

            sup = hyp.get("supporting_observations", [])
            if sup:
                lines.append("**Supporting observations**:")
                for s in sup:
                    lines.append(f"- {s}")
                lines.append("")

            alt = hyp.get("alternative_explanations", [])
            if alt:
                lines.append("**Alternative explanations**:")
                for a in alt:
                    lines.append(f"- {a}")
                lines.append("")

            vsteps = hyp.get("validation_steps", [])
            if vsteps:
                lines.append("**Validation steps**:")
                for v in vsteps:
                    lines.append(f"- {v}")
                lines.append("")

    val_recs = interp.get("validation_recommendations", [])
    if val_recs:
        lines.append("### Overall validation recommendations")
        for rec in val_recs:
            lines.append(f"- {rec}")
        lines.append("")

    safety_flags = interp.get("safety_flags", [])
    if safety_flags:
        lines.append("### Safety flags")
        for flag in safety_flags:
            lines.append(f"- {flag}")
        lines.append("")

    lines.append(
        f"> `clinical_claims_allowed`: **{interp.get('clinical_claims_allowed', False)}**"
    )
    lines.append("")

    return "\n".join(lines)
