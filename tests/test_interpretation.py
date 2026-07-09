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


# ---------------------------------------------------------------------------
# RNA-seq QC interpretation tests
# ---------------------------------------------------------------------------

def _rnaseq_obs(
    category: str,
    status: str,
    sample: str = "sample_R",
    evidence: dict | None = None,
) -> dict:
    return {
        "source": "rnaseq_qc",
        "sample": sample,
        "category": category,
        "status": status,
        "severity": "warning",
        "message": f"{sample}: {category} {status.upper()} test observation.",
        "evidence": evidence or {},
        "suggested_action": "",
    }


class TestRnaseqInterpretation:
    """RNA-seq QC interpretation rules."""

    # --- rRNA fraction ---

    def test_rrna_fail_produces_finding(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        assert result["workflow"] == "rna-seq-qc"
        assert len(result["findings"]) >= 1
        f = result["findings"][0]
        assert "rRNA" in f["observation"] or "rrna" in f["observation"].lower()

    def test_rrna_fail_no_auto_filter(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        for finding in result["findings"]:
            assert finding["should_filter"] is False, (
                f"Finding {finding['finding_id']} incorrectly sets should_filter=True for rRNA"
            )

    def test_rrna_fail_preserve_until_review(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        assert result["findings"][0]["should_preserve_until_review"] is True

    def test_rrna_warn_triggers_finding(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "warn")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1

    def test_rrna_finding_has_technical_and_biological_explanations(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        f = result["findings"][0]
        assert len(f["technical_explanations"]) > 0
        assert len(f["plausible_biological_explanations"]) > 0

    def test_rrna_numeric_threshold_triggers(self):
        """Numeric evidence above threshold triggers interpretation even without explicit status."""
        from genomics_workflow_agent import generate_interpretation

        obs = _rnaseq_obs(
            "rrna_fraction", "pass",
            evidence={"rrna_fraction": 0.25},
        )
        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[obs],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        assert result["findings"][0]["finding_id"].startswith("RNASEQ_RRNA")

    def test_rrna_hypothesis_not_clinical(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        for hyp in result["hypotheses"]:
            assert hyp["clinical_claim"] is False

    # --- Gene-body / 3' bias ---

    def test_gene_body_bias_fail_produces_finding(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("gene_body_coverage", "fail")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        f = result["findings"][0]
        assert "coverage" in f["observation"].lower() or "bias" in f["observation"].lower()

    def test_three_prime_bias_triggers_same_rule(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("three_prime_bias", "warn")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        f = result["findings"][0]
        assert f["finding_id"].startswith("RNASEQ_BIAS")

    def test_gene_body_bias_mentions_rin(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("gene_body_coverage", "fail")],
            decisions=[],
        )
        text = json.dumps(result)
        assert "RIN" in text or "RNA integrity" in text.lower()

    def test_gene_body_bias_confidence_high(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("gene_body_coverage", "fail")],
            decisions=[],
        )
        assert result["findings"][0]["confidence"] == "high"

    # --- Intronic / intergenic mapping ---

    def test_intronic_mapping_fail_produces_finding(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("intronic_mapping", "warn")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        f = result["findings"][0]
        assert f["finding_id"].startswith("RNASEQ_INTRONIC")

    def test_intronic_mapping_has_gdna_technical_explanation(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("intronic_mapping", "warn")],
            decisions=[],
        )
        tech = " ".join(result["findings"][0]["technical_explanations"]).lower()
        assert "gdna" in tech or "dna" in tech or "dnase" in tech

    def test_intronic_mapping_has_intron_retention_biological_explanation(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("intronic_mapping", "warn")],
            decisions=[],
        )
        bio = " ".join(result["findings"][0]["plausible_biological_explanations"]).lower()
        assert "intron retention" in bio or "intron" in bio

    def test_intergenic_mapping_triggers_same_rule(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("intergenic_mapping", "warn")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        assert result["findings"][0]["finding_id"].startswith("RNASEQ_INTRONIC")

    def test_intronic_no_auto_filter(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("intronic_mapping", "fail")],
            decisions=[],
        )
        for f in result["findings"]:
            assert f["should_filter"] is False

    def test_intronic_numeric_threshold_triggers(self):
        from genomics_workflow_agent import generate_interpretation

        obs = _rnaseq_obs(
            "intronic_mapping", "pass",
            evidence={"intronic_mapping_rate": 0.22},
        )
        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[obs],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        assert result["findings"][0]["finding_id"].startswith("RNASEQ_INTRONIC")

    # --- Low mapping / multimapping ---

    def test_low_mapping_produces_finding(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("mapping_rate", "warn")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        assert result["findings"][0]["finding_id"].startswith("RNASEQ_MAP")

    def test_multimapping_triggers_low_mapping_rule(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("multimapping_rate", "fail")],
            decisions=[],
        )
        assert len(result["findings"]) >= 1
        f = result["findings"][0]
        assert f["finding_id"].startswith("RNASEQ_MAP")

    def test_low_mapping_mentions_kraken(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("mapping_rate", "fail")],
            decisions=[],
        )
        text = json.dumps(result)
        assert "Kraken" in text or "taxonomic" in text.lower()

    def test_low_mapping_no_auto_filter(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("mapping_rate", "fail")],
            decisions=[],
        )
        for f in result["findings"]:
            assert f["should_filter"] is False

    # --- All pass ---

    def test_rnaseq_all_pass_no_flags(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[],
            decisions=[],
        )
        assert len(result["findings"]) == 1
        f = result["findings"][0]
        assert f["should_filter"] is False
        assert f["should_preserve_until_review"] is False
        assert len(result["hypotheses"]) == 0

    # --- combined findings ---

    def test_multiple_rnaseq_flags_produce_multiple_findings(self):
        from genomics_workflow_agent import generate_interpretation

        obs = [
            _rnaseq_obs("rrna_fraction", "fail", sample="S1"),
            _rnaseq_obs("gene_body_coverage", "warn", sample="S1"),
            _rnaseq_obs("intronic_mapping", "warn", sample="S1"),
            _rnaseq_obs("mapping_rate", "fail", sample="S1"),
        ]
        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=obs,
            decisions=[],
        )
        assert len(result["findings"]) >= 4
        assert len(result["hypotheses"]) >= 4

    # --- clinical safety ---

    def test_rnaseq_no_forbidden_clinical_terms(self):
        from genomics_workflow_agent import generate_interpretation
        from genomics_workflow_agent.interpretation.safety import scan_result_dict

        categories = [
            ("rrna_fraction", "fail"),
            ("gene_body_coverage", "warn"),
            ("intronic_mapping", "warn"),
            ("mapping_rate", "fail"),
            ("multimapping_rate", "fail"),
        ]
        for category, status in categories:
            result = generate_interpretation(
                workflow="rna-seq-qc",
                observations=[_rnaseq_obs(category, status)],
                decisions=[],
            )
            found = scan_result_dict(result)
            assert found == [], (
                f"Forbidden clinical terms in rna-seq-qc/{category}: {found}"
            )

    def test_rnaseq_clinical_claims_allowed_false(self):
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        assert result["clinical_claims_allowed"] is False

    # --- JSON schema validation ---

    def test_rnaseq_output_validates_against_schema(self):
        """Output must conform to interpretation_result.schema.json."""
        from genomics_workflow_agent import generate_interpretation

        schema_path = (
            Path(__file__).parent.parent
            / "genomics_workflow_agent"
            / "schemas"
            / "interpretation_result.schema.json"
        )
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[
                _rnaseq_obs("rrna_fraction", "fail"),
                _rnaseq_obs("gene_body_coverage", "warn"),
            ],
            decisions=[],
        )

        # Validate required top-level fields
        for required_field in schema.get("required", []):
            assert required_field in result, f"Missing required field: {required_field}"

        # Validate finding structure
        finding_props = schema["properties"]["findings"]["items"].get("required", [])
        for finding in result["findings"]:
            for req in finding_props:
                assert req in finding, f"Finding missing required field: {req}"

        # Validate hypothesis structure
        hyp_props = schema["properties"]["hypotheses"]["items"].get("required", [])
        for hyp in result["hypotheses"]:
            for req in hyp_props:
                assert req in hyp, f"Hypothesis missing required field: {req}"

        # clinical_claims_allowed must be false per schema
        assert result["clinical_claims_allowed"] is False

    # --- Safety scanner: injection test ---

    def test_safety_scanner_catches_injected_forbidden_term(self):
        """Prove the scanner catches forbidden terms in any string, including mock inputs."""
        from genomics_workflow_agent.interpretation.safety import scan_text, scan_result_dict

        # Injected forbidden term in a mock "observation message"
        dirty_text = "This RNA-seq result is diagnostic for a pathogenic mutation."
        found = scan_text(dirty_text)
        assert "diagnostic" in found
        assert "pathogenic" in found

        # Also verify that a dict containing the dirty text is caught
        dirty_dict = {"message": dirty_text, "workflow": "rna-seq-qc"}
        found_dict = scan_result_dict(dirty_dict)
        assert "diagnostic" in found_dict or "pathogenic" in found_dict

    def test_clean_rnaseq_output_passes_safety_scanner(self):
        from genomics_workflow_agent import generate_interpretation
        from genomics_workflow_agent.interpretation.safety import scan_result_dict

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[
                _rnaseq_obs("rrna_fraction", "fail"),
                _rnaseq_obs("intronic_mapping", "warn"),
                _rnaseq_obs("mapping_rate", "fail"),
                _rnaseq_obs("gene_body_coverage", "warn"),
            ],
            decisions=[],
        )
        found = scan_result_dict(result)
        assert found == [], f"Forbidden clinical terms in RNA-seq output: {found}"

    # --- CLI: rna-seq-qc interpret ---

    def test_cli_interpret_rnaseq(self, tmp_path):
        report = {
            "workflow": "rna-seq-qc",
            "observations": [_rnaseq_obs("rrna_fraction", "fail")],
            "decisions": [],
        }
        report_path = tmp_path / "rnaseq_report.json"
        report_path.write_text(json.dumps(report), encoding="utf-8")

        out_dir = tmp_path / "interpretation"
        proc = subprocess.run(
            [
                sys.executable, "-m", "genomics_workflow_agent",
                "interpret",
                "--input", str(report_path),
                "--out", str(out_dir),
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0, f"CLI failed: {proc.stderr}"
        assert (out_dir / "interpretation_report.json").exists()
        assert (out_dir / "interpretation_report.md").exists()

        parsed = json.loads((out_dir / "interpretation_report.json").read_text())
        assert parsed["workflow"] == "rna-seq-qc"
        assert parsed["clinical_claims_allowed"] is False

    # --- Render ---

    def test_rnaseq_render_markdown_has_required_sections(self):
        from genomics_workflow_agent import generate_interpretation
        from genomics_workflow_agent.interpretation import render_interpretation_md

        result = generate_interpretation(
            workflow="rna-seq-qc",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        md = render_interpretation_md(result)
        assert "## Biological interpretation and hypotheses" in md
        assert "rna-seq-qc" in md
        assert "clinical_claims_allowed" in md
        # Technical and biological sections present
        assert "Technical explanations" in md
        assert "biological explanations" in md.lower()

    def test_rnaseq_workflow_alias_rnaseq(self):
        """workflow='rnaseq' (no hyphen) is also accepted."""
        from genomics_workflow_agent import generate_interpretation

        result = generate_interpretation(
            workflow="rnaseq",
            observations=[_rnaseq_obs("rrna_fraction", "fail")],
            decisions=[],
        )
        assert result["workflow"] == "rna-seq-qc"
        assert result["clinical_claims_allowed"] is False
