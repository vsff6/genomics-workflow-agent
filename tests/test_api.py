"""Tests for the public Python API (genomics_workflow_agent.api)."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

REQUIRED_KEYS = {
    "status", "workflow", "input_path", "outdir",
    "summary", "warnings", "errors", "paths", "provenance_paths",
}


def _is_json_serializable(obj) -> bool:
    try:
        json.dumps(obj)
        return True
    except (TypeError, ValueError):
        return False


class TestApiImports:
    def test_top_level_imports(self):
        from genomics_workflow_agent import (
            inspect_inputs,
            plan_workflow,
            run_workflow,
            run_fastq_qc_agent,
            run_variant_qc_agent,
            write_report,
        )
        for fn in [inspect_inputs, plan_workflow, run_workflow,
                   run_fastq_qc_agent, run_variant_qc_agent, write_report]:
            assert callable(fn)

    def test_api_module_imports(self):
        from genomics_workflow_agent.api import (
            inspect_inputs,
            plan_workflow,
            run_workflow,
            run_fastq_qc_agent,
            run_variant_qc_agent,
            write_report,
        )
        assert True  # if we get here, imports are fine


class TestInspectInputs:
    def test_returns_dict_with_required_keys(self, tmp_path):
        from genomics_workflow_agent import inspect_inputs
        result = inspect_inputs(tmp_path)
        assert isinstance(result, dict)
        assert REQUIRED_KEYS.issubset(result.keys())

    def test_result_is_json_serializable(self, tmp_path):
        from genomics_workflow_agent import inspect_inputs
        result = inspect_inputs(tmp_path)
        assert _is_json_serializable(result)

    def test_missing_input_returns_failure(self):
        from genomics_workflow_agent import inspect_inputs
        result = inspect_inputs("/nonexistent/path/that/does/not/exist")
        assert result["status"] == "failed"
        assert len(result["errors"]) > 0
        assert _is_json_serializable(result)

    def test_success_status_for_existing_dir(self, tmp_path):
        from genomics_workflow_agent import inspect_inputs
        result = inspect_inputs(tmp_path)
        assert result["status"] == "success"


class TestPlanWorkflow:
    def test_returns_dict_with_required_keys(self, tmp_path):
        from genomics_workflow_agent import plan_workflow
        result = plan_workflow(tmp_path, workflow="fastq-qc", outdir=tmp_path / "plan_out")
        assert isinstance(result, dict)
        assert REQUIRED_KEYS.issubset(result.keys())

    def test_result_is_json_serializable(self, tmp_path):
        from genomics_workflow_agent import plan_workflow
        result = plan_workflow(tmp_path, outdir=tmp_path / "plan_out")
        assert _is_json_serializable(result)

    def test_missing_input_returns_failure(self, tmp_path):
        from genomics_workflow_agent import plan_workflow
        result = plan_workflow("/no/such/path", outdir=tmp_path / "out")
        assert result["status"] == "failed"
        assert len(result["errors"]) > 0

    def test_unsupported_workflow_returns_failure(self, tmp_path):
        from genomics_workflow_agent import plan_workflow
        result = plan_workflow(tmp_path, workflow="not-a-workflow", outdir=tmp_path / "out")
        assert result["status"] == "failed"
        assert any("not-a-workflow" in e or "Unsupported" in e for e in result["errors"])
        assert _is_json_serializable(result)

    def test_status_is_dry_run(self, tmp_path):
        from genomics_workflow_agent import plan_workflow
        result = plan_workflow(tmp_path, workflow="fastq-qc", outdir=tmp_path / "out")
        assert result["status"] == "dry_run"

    def test_plan_json_is_written(self, tmp_path):
        from genomics_workflow_agent import plan_workflow
        out = tmp_path / "plan_out"
        plan_workflow(tmp_path, workflow="fastq-qc", outdir=out)
        assert (out / "plan.json").exists()


class TestRunWorkflow:
    def test_returns_dict_with_required_keys(self, tmp_path):
        from genomics_workflow_agent import run_workflow
        result = run_workflow(tmp_path, workflow="fastq-qc", outdir=tmp_path / "out")
        assert isinstance(result, dict)
        assert REQUIRED_KEYS.issubset(result.keys())

    def test_result_is_json_serializable(self, tmp_path):
        from genomics_workflow_agent import run_workflow
        result = run_workflow(tmp_path, outdir=tmp_path / "out")
        assert _is_json_serializable(result)

    def test_missing_input_returns_failure(self, tmp_path):
        from genomics_workflow_agent import run_workflow
        result = run_workflow("/no/such/path", outdir=tmp_path / "out")
        assert result["status"] == "failed"

    def test_unsupported_workflow_returns_failure(self, tmp_path):
        from genomics_workflow_agent import run_workflow
        result = run_workflow(tmp_path, workflow="bad-workflow", outdir=tmp_path / "out")
        assert result["status"] == "failed"
        assert _is_json_serializable(result)

    def test_dry_run_does_not_call_subprocess(self, tmp_path):
        from genomics_workflow_agent import run_workflow
        fake_tools = {t: {"available": False, "version": None}
                      for t in ["fastqc", "multiqc", "fastp", "cutadapt",
                                "nextflow", "docker", "singularity", "conda"]}
        with patch("genomics_workflow_agent.workflows.fastq_qc.check_tools",
                   return_value=fake_tools):
            with patch("subprocess.run") as mock_run:
                run_workflow(tmp_path, execute=False, outdir=tmp_path / "out")
                mock_run.assert_not_called()

    def test_status_is_dry_run_when_not_execute(self, tmp_path):
        from genomics_workflow_agent import run_workflow
        result = run_workflow(tmp_path, execute=False, outdir=tmp_path / "out")
        assert result["status"] == "dry_run"


class TestRunFastqQcAgent:
    def test_returns_dict_with_required_keys(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        result = run_fastq_qc_agent(tmp_path, outdir=tmp_path / "out")
        assert isinstance(result, dict)
        assert REQUIRED_KEYS.issubset(result.keys())

    def test_result_is_json_serializable(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        result = run_fastq_qc_agent(tmp_path, outdir=tmp_path / "out")
        assert _is_json_serializable(result)

    def test_missing_input_returns_failure(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        result = run_fastq_qc_agent("/no/such/path", outdir=tmp_path / "out")
        assert result["status"] == "failed"
        assert _is_json_serializable(result)

    def test_auto_trim_without_execute_returns_failure(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        result = run_fastq_qc_agent(
            tmp_path, outdir=tmp_path / "out", auto_trim=True, execute=False
        )
        assert result["status"] == "failed"
        assert any("auto_trim" in e or "execute" in e.lower() for e in result["errors"])

    def test_dry_run_status(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        result = run_fastq_qc_agent(tmp_path, outdir=tmp_path / "out", execute=False)
        assert result["status"] == "dry_run"

    def test_dry_run_does_not_call_subprocess(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        fake_tools = {t: {"available": False, "version": None}
                      for t in ["fastqc", "multiqc", "fastp", "cutadapt",
                                "nextflow", "docker", "singularity", "conda"]}
        with patch("genomics_workflow_agent.workflows.fastq_qc.check_tools",
                   return_value=fake_tools):
            with patch("subprocess.run") as mock_run:
                run_fastq_qc_agent(tmp_path, outdir=tmp_path / "out", execute=False)
                mock_run.assert_not_called()

    def test_observations_and_decisions_in_result(self, tmp_path):
        from genomics_workflow_agent import run_fastq_qc_agent
        result = run_fastq_qc_agent(tmp_path, outdir=tmp_path / "out")
        assert "observations" in result
        assert "decisions" in result
        assert "recommended_actions" in result
        assert isinstance(result["observations"], list)


class TestRunVariantQcAgent:
    def test_returns_dict_with_required_keys(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        result = run_variant_qc_agent(tmp_path, outdir=tmp_path / "out")
        assert isinstance(result, dict)
        assert REQUIRED_KEYS.issubset(result.keys())

    def test_result_is_json_serializable(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        result = run_variant_qc_agent(tmp_path, outdir=tmp_path / "out")
        assert _is_json_serializable(result)

    def test_missing_input_returns_failure(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        result = run_variant_qc_agent("/no/such/path", outdir=tmp_path / "out")
        assert result["status"] == "failed"
        assert _is_json_serializable(result)

    def test_dry_run_status(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        result = run_variant_qc_agent(tmp_path, outdir=tmp_path / "out", execute=False)
        assert result["status"] == "dry_run"

    def test_dry_run_does_not_execute_bioinformatics_tools(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        fake_tools = {t: {"available": False, "version": None}
                      for t in ["samtools", "bcftools", "mosdepth", "nextflow",
                                "docker", "singularity", "conda"]}
        with patch("genomics_workflow_agent.workflows.variant_qc.check_tools",
                   return_value=fake_tools):
            with patch("subprocess.run") as mock_run:
                run_variant_qc_agent(tmp_path, outdir=tmp_path / "out", execute=False)
                mock_run.assert_not_called()

    def test_no_clinical_claims_in_observations(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        result = run_variant_qc_agent(tmp_path, outdir=tmp_path / "out")
        for obs in result.get("observations", []):
            msg = obs.get("message", "").lower()
            for forbidden in ["pathogenic", "benign", "diagnosis", "diagnoses"]:
                assert forbidden not in msg, \
                    f"Clinical term {forbidden!r} found in observation: {msg}"

    def test_report_json_is_written(self, tmp_path):
        from genomics_workflow_agent import run_variant_qc_agent
        out = tmp_path / "out"
        run_variant_qc_agent(tmp_path, outdir=out)
        assert (out / "variant_agent_report.json").exists()


class TestWriteReport:
    def test_returns_dict_with_required_keys(self, tmp_path):
        from genomics_workflow_agent import write_report
        result = write_report(tmp_path)
        assert isinstance(result, dict)
        assert REQUIRED_KEYS.issubset(result.keys())

    def test_result_is_json_serializable(self, tmp_path):
        from genomics_workflow_agent import write_report
        result = write_report(tmp_path)
        assert _is_json_serializable(result)

    def test_missing_results_dir_returns_failure(self, tmp_path):
        from genomics_workflow_agent import write_report
        result = write_report("/no/such/path")
        assert result["status"] == "failed"
        assert _is_json_serializable(result)

    def test_success_writes_final_report(self, tmp_path):
        from genomics_workflow_agent import write_report
        # put a dummy json file in tmp_path so there is something to aggregate
        (tmp_path / "run_report.json").write_text(
            json.dumps({"workflow": "test"}), encoding="utf-8"
        )
        result = write_report(tmp_path)
        assert result["status"] == "success"
        assert Path(result["paths"]["final_report_json"]).exists()
