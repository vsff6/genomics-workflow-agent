"""Biological interpretation scaffold for genomics QC outputs."""
from __future__ import annotations

from genomics_workflow_agent.interpretation.hypothesis_generator import generate_interpretation
from genomics_workflow_agent.interpretation.render import render_interpretation_md
from genomics_workflow_agent.interpretation.safety import scan_result_dict, scan_text

__all__ = [
    "generate_interpretation",
    "render_interpretation_md",
    "scan_result_dict",
    "scan_text",
]
