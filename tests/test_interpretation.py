"""Tests for the biological interpretation scaffold."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _obs(category: str, status: str, sample: str = "sample_A", source: str = "decision_engine") -> dict:
    return {
        "source": source,
        "sample": sample,
        "category": category,
        "status": status,
        "severity": "warning",
        "message": f"{sample}: {category} {status.upper()} test observation.",
        "evidence": {"module": category, "fastqc_status": status},
        "suggested_action": "",
    }


# ---------------------------------------------------------------------------
# Import / public API
# ---------------------------------------------------------------------------

def test_public_api_import():
    from genomics_workflow_agent import generate_interpretation
    assert callable(generate_interpretation)


def test_interpretation_module_import():
    from genomics_workflow_agent.interpretation import (
        generate_interpretation,
        render_interpretation_md,
        scan_result_dict,
        scan_text,
    )
    assert callable(generate_interpretation)
    assert callable(render_interpretation_md)


# ---------------------------------------------------------------------------
# JSON serializability
# ---------------------------------------------------------------------------

def test_fastq_interpretation_is_json_serializable():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("gc_content", "fail")],
        decisions=[],
    )
    serialized = json.dumps(result)
    assert isinstance(serialized, str)
    parsed = json.loads(serialized)
    assert parsed["clinical_claims_allowed"] is False


def test_variant_interpretation_is_json_serializable():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="variant-qc",
        observations=[_obs("alignment", "warn", source="variant_decision_engine")],
        decisions=[],
    )
    serialized = json.dumps(result)
    assert isinstance(serialized, str)
    parsed = json.loads(serialized)
    assert parsed["clinical_claims_allowed"] is False


# ---------------------------------------------------------------------------
# FASTQ QC rules
# ---------------------------------------------------------------------------

def test_gc_fail_produces_biological_review_finding():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("gc_content", "fail")],
        decisions=[],
    )
    assert len(result["findings"]) >= 1
    gc_finding = result["findings"][0]
    assert "GC" in gc_finding["observation"] or "gc_content" in gc_finding["evidence_source"].lower()
    # Must NOT recommend automatic filtering
    assert gc_finding["should_filter"] is False
    assert gc_finding["should_preserve_until_review"] is True
    # Must have biological explanations
    assert len(gc_finding["plausible_biological_explanations"]) > 0


def test_gc_fail_no_auto_filter():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("gc_content", "fail")],
        decisions=[],
    )
    for finding in result["findings"]:
        assert finding["should_filter"] is False, (
            f"Finding {finding['finding_id']} incorrectly recommends filtering for GC content"
        )


def test_gc_fail_has_hypothesis():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("gc_content", "fail")],
        decisions=[],
    )
    assert len(result["hypotheses"]) >= 1
    hyp = result["hypotheses"][0]
    assert hyp["clinical_claim"] is False
    assert hyp["interpretation_type"] == "ambiguous"


def test_adapter_warn_produces_technical_finding():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("adapter_content", "warn")],
        decisions=[],
    )
    assert len(result["findings"]) >= 1
    finding = result["findings"][0]
    assert "adapter" in finding["observation"].lower() or "Adapter" in finding["observation"]
    # Adapter is a technical artifact - should_preserve_until_review may be False
    assert len(finding["technical_explanations"]) > 0
    assert finding["should_filter"] is False


def test_adapter_fail_produces_trimming_recommendation():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("adapter_content", "fail")],
        decisions=[],
    )
    finding = result["findings"][0]
    action = finding["recommended_action"].lower()
    assert "trim" in action


def test_all_pass_no_biological_warning():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("overall_qc", "pass")],
        decisions=[],
    )
    # Should produce an all-pass finding, no strong biological hypothesis
    assert len(result["findings"]) >= 1
    assert len(result["hypotheses"]) == 0


def test_empty_observations_fastq():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(workflow="fastq-qc", observations=[], decisions=[])
    # Should produce an all-pass finding
    assert len(result["findings"]) >= 1


# ---------------------------------------------------------------------------
# Variant QC rules
# ---------------------------------------------------------------------------

def test_zero_vcf_not_evidence_of_absence():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="variant-qc",
        observations=[_obs("vcf_content", "fail", source="variant_decision_engine")],
        decisions=[],
    )
    assert len(result["findings"]) >= 1
    finding = result["findings"][0]
    text = json.dumps(finding)
    # Must explicitly state this is not evidence of absence
    assert "absence" in text.lower() or "zero records" in text.lower()
    assert finding["should_preserve_until_review"] is True


def test_low_coverage_produces_cautious_finding():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="variant-qc",
        observations=[_obs("coverage", "warn", source="variant_decision_engine")],
        decisions=[],
    )
    assert len(result["findings"]) >= 1
    finding = result["findings"][0]
    assert finding["should_filter"] is False
    # Must not infer genotype or make clinical claims
    text = json.dumps(finding).lower()
    assert "genotype" not in text or "infer genotype" in text  # cautionary language OK


def test_low_mapping_has_hypothesis():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="variant-qc",
        observations=[_obs("alignment", "warn", source="variant_decision_engine")],
        decisions=[],
    )
    assert len(result["hypotheses"]) >= 1
    hyp = result["hypotheses"][0]
    assert hyp["clinical_claim"] is False


def test_variant_no_issues_produces_acceptance_finding():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="variant-qc",
        observations=[_obs("alignment", "pass", source="variant_decision_engine")],
        decisions=[],
    )
    assert len(result["findings"]) >= 1
    finding = result["findings"][0]
    assert finding["should_filter"] is False


# ---------------------------------------------------------------------------
# Safety: no forbidden clinical terms
# ---------------------------------------------------------------------------

CLINICAL_TERMS = [
    "pathogenic",
    "benign",
    "likely pathogenic",
    "likely benign",
    "disease-causing",
    "therapy recommendation",
    "clinical action",
    "medical action",
]


@pytest.mark.parametrize("workflow,category,status", [
    ("fastq-qc", "gc_content", "fail"),
    ("fastq-qc", "adapter_content", "warn"),
    ("fastq-qc", "per_base_quality", "fail"),
    ("fastq-qc", "overrepresented_sequences", "warn"),
    ("variant-qc", "alignment", "warn"),
    ("variant-qc", "vcf_content", "fail"),
    ("variant-qc", "coverage", "warn"),
    ("variant-qc", "contig_coverage", "warn"),
])
def test_no_forbidden_clinical_terms(workflow, category, status):
    from genomics_workflow_agent import generate_interpretation
    from genomics_workflow_agent.interpretation.safety import scan_result_dict

    source = "variant_decision_engine" if workflow == "variant-qc" else "decision_engine"
    result = generate_interpretation(
        workflow=workflow,
        observations=[_obs(category, status, source=source)],
        decisions=[],
    )
    found = scan_result_dict(result)
    assert found == [], f"Forbidden clinical terms found in {workflow}/{category}: {found}"


def test_clinical_claims_allowed_always_false():
    from genomics_workflow_agent import generate_interpretation

    for workflow in ("fastq-qc", "variant-qc"):
        result = generate_interpretation(workflow=workflow, observations=[], decisions=[])
        assert result["clinical_claims_allowed"] is False


def test_hypothesis_clinical_claim_false():
    from genomics_workflow_agent import generate_interpretation

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("gc_content", "fail")],
        decisions=[],
    )
    for hyp in result["hypotheses"]:
        assert hyp["clinical_claim"] is False


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------

def test_render_produces_markdown():
    from genomics_workflow_agent.interpretation import generate_interpretation, render_interpretation_md

    result = generate_interpretation(
        workflow="fastq-qc",
        observations=[_obs("gc_content", "fail")],
        decisions=[],
    )
    md = render_interpretation_md(result)
    assert "## Biological interpretation and hypotheses" in md
    assert "clinical_claims_allowed" in md
    assert "false" in md.lower() or "False" in md


# ---------------------------------------------------------------------------
# CLI: interpret subcommand
# ---------------------------------------------------------------------------

def test_cli_interpret_writes_json_and_md(tmp_path):
    # Create a minimal agent_report.json
    report = {
        "workflow": "fastq-qc",
        "observations": [_obs("gc_content", "fail")],
        "decisions": [],
        "recommended_actions": [],
        "warnings": [],
        "limitations": [],
    }
    report_path = tmp_path / "agent_report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")

    out_dir = tmp_path / "interpretation"
    result = subprocess.run(
        [
            sys.executable, "-m", "genomics_workflow_agent",
            "interpret",
            "--input", str(report_path),
            "--out", str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    assert (out_dir / "interpretation_report.json").exists()
    assert (out_dir / "interpretation_report.md").exists()

    parsed = json.loads((out_dir / "interpretation_report.json").read_text(encoding="utf-8"))
    assert parsed["clinical_claims_allowed"] is False
    assert len(parsed["findings"]) >= 1


def test_cli_interpret_variant(tmp_path):
    report = {
        "workflow": "variant-qc",
        "observations": [_obs("vcf_content", "fail", source="variant_decision_engine")],
        "decisions": [],
        "recommended_actions": [],
        "warnings": [],
        "limitations": [],
    }
    report_path = tmp_path / "variant_agent_report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")

    out_dir = tmp_path / "interpretation"
    result = subprocess.run(
        [
            sys.executable, "-m", "genomics_workflow_agent",
            "interpret",
            "--input", str(report_path),
            "--out", str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    assert (out_dir / "interpretation_report.json").exists()
    parsed = json.loads((out_dir / "interpretation_report.json").read_text(encoding="utf-8"))
    assert parsed["workflow"] == "variant-qc"


def test_cli_interpret_workflow_override(tmp_path):
    # Report without workflow field - must use --workflow flag
    report = {
        "observations": [_obs("coverage", "warn", source="variant_decision_engine")],
        "decisions": [],
    }
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")

    out_dir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable, "-m", "genomics_workflow_agent",
            "interpret",
            "--input", str(report_path),
            "--workflow", "variant-qc",
            "--out", str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"
    parsed = json.loads((out_dir / "interpretation_report.json").read_text(encoding="utf-8"))
    assert parsed["workflow"] == "variant-qc"


# ---------------------------------------------------------------------------
# Safety scanner unit tests
# ---------------------------------------------------------------------------

def test_safety_scanner_catches_forbidden_terms():
    from genomics_workflow_agent.interpretation.safety import scan_text

    found = scan_text("This variant is pathogenic and requires treatment.")
    assert "pathogenic" in found
    assert "treatment" in found


def test_safety_scanner_clean_text():
    from genomics_workflow_agent.interpretation.safety import scan_text

    found = scan_text("Review GC content distribution with biological context.")
    assert found == []


# ---------------------------------------------------------------------------
# Schema structure
# ---------------------------------------------------------------------------

def test_schema_file_exists():
    schema_path = Path(__file__).parent.parent / \
        "genomics_workflow_agent" / "schemas" / "interpretation_result.schema.json"
    assert schema_path.exists()
    data = json.loads(schema_path.read_text(encoding="utf-8"))
    assert data["title"] == "InterpretationResult"
    assert "findings" in data["properties"]
    assert "hypotheses" in data["properties"]


# ---------------------------------------------------------------------------
# Integration: agent report JSON includes biological_interpretation
# ---------------------------------------------------------------------------

def test_fastq_agent_report_includes_interpretation(tmp_path):
    from genomics_workflow_agent.agent.state import AgentState, Observation
    from genomics_workflow_agent.agent.fastq_agent import write_agent_report_json

    state = AgentState(input_path=str(tmp_path), workflow="fastq-qc")
    state.observations.append(Observation(
        source="decision_engine",
        sample="sample_A",
        category="gc_content",
        status="fail",
        severity="warning",
        message="sample_A: GC content FAIL.",
        evidence={"module": "Per sequence GC content", "fastqc_status": "fail"},
        suggested_action="Review.",
    ))

    json_path = write_agent_report_json(state, tmp_path / "agent_report.json")
    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert "biological_interpretation" in data
    interp = data["biological_interpretation"]
    assert interp["clinical_claims_allowed"] is False
    assert len(interp["findings"]) >= 1


def test_variant_agent_report_includes_interpretation(tmp_path):
    from genomics_workflow_agent.agent.state import AgentState, Observation
    from genomics_workflow_agent.agent.variant_agent import write_variant_agent_report_json

    state = AgentState(input_path=str(tmp_path), workflow="variant-qc")
    state.observations.append(Observation(
        source="variant_decision_engine",
        sample="sample_B",
        category="vcf_content",
        status="fail",
        severity="critical",
        message="sample_B: VCF contains zero variant records.",
        evidence={"n_records": 0},
        suggested_action="Investigate.",
    ))

    json_path = write_variant_agent_report_json(state, tmp_path / "variant_agent_report.json")
    data = json.loads(json_path.read_text(encoding="utf-8"))
    assert "biological_interpretation" in data
    interp = data["biological_interpretation"]
    assert interp["clinical_claims_allowed"] is False
    assert len(interp["findings"]) >= 1
