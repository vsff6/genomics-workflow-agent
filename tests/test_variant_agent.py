"""Tests for variant-qc parsers, decision engine, agent, and CLI."""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
SAMTOOLS_FIXTURES = FIXTURES / "samtools"
BCFTOOLS_FIXTURES = FIXTURES / "bcftools"
MOSDEPTH_FIXTURES = FIXTURES / "mosdepth"


class TestSamtoolsFlagstatParser:
    def test_normal_flagstat_parse_ok(self):
        from genomics_workflow_agent.parsers.samtools import parse_flagstat
        text = (SAMTOOLS_FIXTURES / "sample1_flagstat.txt").read_text(encoding="utf-8")
        result = parse_flagstat(text, sample="sample1")
        assert result["parse_ok"] is True
        assert result["total_reads"] == 100000
        assert result["mapped_reads"] == 98500
        assert result["mapped_pct"] == pytest.approx(98.50)

    def test_low_mapping_flagstat(self):
        from genomics_workflow_agent.parsers.samtools import parse_flagstat
        text = (SAMTOOLS_FIXTURES / "lowmap_flagstat.txt").read_text(encoding="utf-8")
        result = parse_flagstat(text, sample="lowmap")
        assert result["parse_ok"] is True
        assert result["mapped_pct"] == pytest.approx(40.0)

    def test_malformed_flagstat_does_not_crash(self):
        from genomics_workflow_agent.parsers.samtools import parse_flagstat
        text = (SAMTOOLS_FIXTURES / "malformed_flagstat.txt").read_text(encoding="utf-8")
        result = parse_flagstat(text, sample="bad")
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_empty_flagstat(self):
        from genomics_workflow_agent.parsers.samtools import parse_flagstat
        result = parse_flagstat("", sample="empty")
        assert result["parse_ok"] is False
        assert "Empty" in result["errors"][0]

    def test_flagstat_file_function(self):
        from genomics_workflow_agent.parsers.samtools import parse_flagstat_file
        path = SAMTOOLS_FIXTURES / "sample1_flagstat.txt"
        result = parse_flagstat_file(path)
        assert result["parse_ok"] is True
        assert result["sample"] == "sample1"

    def test_properly_paired_extracted(self):
        from genomics_workflow_agent.parsers.samtools import parse_flagstat
        text = (SAMTOOLS_FIXTURES / "sample1_flagstat.txt").read_text(encoding="utf-8")
        result = parse_flagstat(text)
        assert result["properly_paired"] == 95000
        assert result["properly_paired_pct"] == pytest.approx(95.0)


class TestSamtoolsIdxstatsParser:
    def test_normal_idxstats(self):
        from genomics_workflow_agent.parsers.samtools import parse_idxstats
        text = (SAMTOOLS_FIXTURES / "sample1_idxstats.txt").read_text(encoding="utf-8")
        result = parse_idxstats(text, sample="sample1")
        assert result["parse_ok"] is True
        assert len(result["contigs"]) == 6
        assert result["total_mapped"] > 0
        # '*' contig should not appear in zero_read_contigs
        assert "*" not in result["zero_read_contigs"]

    def test_empty_idxstats(self):
        from genomics_workflow_agent.parsers.samtools import parse_idxstats
        result = parse_idxstats("", sample="empty")
        assert result["parse_ok"] is False

    def test_idxstats_file_function(self):
        from genomics_workflow_agent.parsers.samtools import parse_idxstats_file
        result = parse_idxstats_file(SAMTOOLS_FIXTURES / "sample1_idxstats.txt")
        assert result["parse_ok"] is True
        assert result["sample"] == "sample1"


class TestSamtoolsStatsParser:
    def test_empty_stats_does_not_crash(self):
        from genomics_workflow_agent.parsers.samtools import parse_stats
        result = parse_stats("")
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_stats_with_sn_lines(self):
        from genomics_workflow_agent.parsers.samtools import parse_stats
        text = (
            "SN\traw total sequences:\t100000\n"
            "SN\treads mapped:\t98500\n"
            "SN\taverage length:\t150\n"
            "SN\terror rate:\t1.5e-03\n"
        )
        result = parse_stats(text, sample="s1")
        assert result["parse_ok"] is True
        assert result["raw_sn"]["raw total sequences"] == 100000
        assert result["summary"]["mapped_reads"] == 98500


class TestBcftoolsStatsParser:
    def test_normal_bcftools_stats(self):
        from genomics_workflow_agent.parsers.bcftools import parse_bcftools_stats
        text = (BCFTOOLS_FIXTURES / "sample1_bcftools_stats.txt").read_text(encoding="utf-8")
        result = parse_bcftools_stats(text, sample="sample1")
        assert result["parse_ok"] is True
        assert result["n_records"] == 45678
        assert result["n_snps"] == 38000
        assert result["n_indels"] == 6000
        assert result["ts_tv"] == pytest.approx(2.15)

    def test_zero_records_bcftools_stats(self):
        from genomics_workflow_agent.parsers.bcftools import parse_bcftools_stats
        text = (BCFTOOLS_FIXTURES / "zero_records_bcftools_stats.txt").read_text(encoding="utf-8")
        result = parse_bcftools_stats(text, sample="empty_vcf")
        assert result["parse_ok"] is True
        assert result["n_records"] == 0
        assert result["n_snps"] == 0

    def test_malformed_bcftools_stats_does_not_crash(self):
        from genomics_workflow_agent.parsers.bcftools import parse_bcftools_stats
        text = (BCFTOOLS_FIXTURES / "malformed_bcftools_stats.txt").read_text(encoding="utf-8")
        result = parse_bcftools_stats(text, sample="bad")
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_empty_bcftools_stats(self):
        from genomics_workflow_agent.parsers.bcftools import parse_bcftools_stats
        result = parse_bcftools_stats("", sample="empty")
        assert result["parse_ok"] is False

    def test_bcftools_stats_file_function(self):
        from genomics_workflow_agent.parsers.bcftools import parse_bcftools_stats_file
        result = parse_bcftools_stats_file(BCFTOOLS_FIXTURES / "sample1_bcftools_stats.txt")
        assert result["parse_ok"] is True
        assert result["sample"] == "sample1"


class TestMosdepthParser:
    def test_normal_mosdepth_summary(self):
        from genomics_workflow_agent.parsers.mosdepth import parse_mosdepth_summary
        text = (MOSDEPTH_FIXTURES / "sample1.mosdepth.summary.txt").read_text(encoding="utf-8")
        result = parse_mosdepth_summary(text, sample="sample1")
        assert result["parse_ok"] is True
        assert result["mean_coverage"] == pytest.approx(30.0)
        assert len(result["regions"]) > 0

    def test_low_coverage_mosdepth(self):
        from genomics_workflow_agent.parsers.mosdepth import parse_mosdepth_summary
        text = (MOSDEPTH_FIXTURES / "lowcov.mosdepth.summary.txt").read_text(encoding="utf-8")
        result = parse_mosdepth_summary(text, sample="lowcov")
        assert result["parse_ok"] is True
        assert result["mean_coverage"] == pytest.approx(3.0)

    def test_empty_mosdepth_does_not_crash(self):
        from genomics_workflow_agent.parsers.mosdepth import parse_mosdepth_summary
        result = parse_mosdepth_summary("", sample="empty")
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_mosdepth_file_function(self):
        from genomics_workflow_agent.parsers.mosdepth import parse_mosdepth_summary_file
        result = parse_mosdepth_summary_file(MOSDEPTH_FIXTURES / "sample1.mosdepth.summary.txt")
        assert result["parse_ok"] is True


class TestVariantDecisionEngine:
    def _make_flagstat(self, sample, mapped_pct, total=100000):
        return {
            "sample": sample,
            "parse_ok": True,
            "errors": [],
            "total_reads": total,
            "mapped_reads": int(total * mapped_pct / 100),
            "mapped_pct": float(mapped_pct),
        }

    def _make_bcftools(self, sample, n_records, n_snps=None):
        return {
            "sample": sample,
            "parse_ok": True,
            "errors": [],
            "n_records": n_records,
            "n_snps": n_snps or n_records,
            "n_indels": 0,
            "n_multiallelic": 0,
            "ts_tv": 2.1,
        }

    def test_low_mapping_triggers_review_decision(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        fs = self._make_flagstat("sA", mapped_pct=40.0)
        result = evaluate_variant_qc_results(
            {"flagstat": [fs], "idxstats": [], "stats": []}, [], []
        )
        obs = result["observations"]
        assert any(o.category == "alignment" and o.status == "warn" for o in obs)
        assert any(d.decision_type == "review" for d in result["decisions"])

    def test_zero_vcf_records_triggers_review_decision(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        bc = self._make_bcftools("sA", n_records=0)
        result = evaluate_variant_qc_results(
            {"flagstat": [], "idxstats": [], "stats": []}, [bc], []
        )
        obs = result["observations"]
        assert any(o.category == "vcf_content" and o.status == "fail" for o in obs)
        assert any(d.decision_type == "review" for d in result["decisions"])

    def test_good_data_produces_accept_decision(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        fs = self._make_flagstat("sA", mapped_pct=98.5)
        bc = self._make_bcftools("sA", n_records=45000, n_snps=38000)
        result = evaluate_variant_qc_results(
            {"flagstat": [fs], "idxstats": [], "stats": []}, [bc], []
        )
        assert any(d.decision_type == "accept" for d in result["decisions"])

    def test_no_data_produces_missing_observation(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        result = evaluate_variant_qc_results(
            {"flagstat": [], "idxstats": [], "stats": []}, [], []
        )
        assert any(o.status == "missing" for o in result["observations"])

    def test_clinical_disclaimer_present(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        result = evaluate_variant_qc_results(
            {"flagstat": [], "idxstats": [], "stats": []}, [], []
        )
        assert "clinical_disclaimer" in result
        disc = result["clinical_disclaimer"].lower()
        assert "not" in disc
        assert "clinical" in disc

    def test_no_assertive_clinical_language_in_disclaimer(self):
        from genomics_workflow_agent.agent.variant_decision_engine import CLINICAL_DISCLAIMER
        lower = CLINICAL_DISCLAIMER.lower()
        # Ensure the disclaimer doesn't make positive clinical assertions
        for forbidden in ["is pathogenic", "is benign", "causes disease", "is diagnostic"]:
            assert forbidden not in lower, f"Assertive clinical term found: {forbidden!r}"

    def test_low_coverage_triggers_warning_observation(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        md = {
            "sample": "sA",
            "parse_ok": True,
            "errors": [],
            "mean_coverage": 3.0,
            "regions": [],
        }
        result = evaluate_variant_qc_results(
            {"flagstat": [], "idxstats": [], "stats": []}, [], [md]
        )
        assert any(o.category == "coverage" and o.status == "warn" for o in result["observations"])

    def test_malformed_flagstat_does_not_crash_engine(self):
        from genomics_workflow_agent.agent.variant_decision_engine import evaluate_variant_qc_results
        bad_fs = {
            "sample": "bad",
            "parse_ok": False,
            "errors": ["No recognizable flagstat lines found"],
        }
        result = evaluate_variant_qc_results(
            {"flagstat": [bad_fs], "idxstats": [], "stats": []}, [], []
        )
        assert isinstance(result["observations"], list)


class TestVariantAgentDryRun:
    def test_dry_run_returns_agent_state(self, tmp_path):
        from genomics_workflow_agent.agent.variant_agent import run_variant_agent
        state = run_variant_agent(tmp_path, tmp_path / "out", execute=False)
        assert state.workflow == "variant-qc"
        assert any("Dry-run" in w for w in state.warnings)
        assert len(state.observations) > 0

    def test_dry_run_writes_reports(self, tmp_path):
        from genomics_workflow_agent.agent.variant_agent import (
            run_variant_agent,
            write_variant_agent_report_json,
            write_variant_agent_report_md,
        )
        out = tmp_path / "out"
        state = run_variant_agent(tmp_path, out, execute=False)
        json_path = write_variant_agent_report_json(state, out / "variant_agent_report.json")
        md_path = write_variant_agent_report_md(state, out / "variant_agent_report.md")
        assert json_path.exists()
        assert md_path.exists()
        data = json.loads(json_path.read_text(encoding="utf-8"))
        assert "clinical_disclaimer" in data

    def test_dry_run_does_not_execute_bioinformatics_tools(self, tmp_path):
        from genomics_workflow_agent.agent.variant_agent import run_variant_agent
        fake_tools = {t: {"available": False, "version": None}
                      for t in ["samtools", "bcftools", "mosdepth", "nextflow",
                                "docker", "singularity", "conda"]}
        with patch("genomics_workflow_agent.workflows.variant_qc.check_tools",
                   return_value=fake_tools):
            with patch("subprocess.run") as mock_run:
                run_variant_agent(tmp_path, tmp_path / "out", execute=False)
                mock_run.assert_not_called()

    def test_report_contains_no_clinical_claims(self, tmp_path):
        from genomics_workflow_agent.agent.variant_agent import (
            run_variant_agent,
            write_variant_agent_report_md,
        )
        out = tmp_path / "out"
        state = run_variant_agent(tmp_path, out, execute=False)
        md_path = write_variant_agent_report_md(state, out / "variant_agent_report.md")
        text = md_path.read_text(encoding="utf-8").lower()
        # Ensure no positive clinical assertions appear (saying something "is diagnostic" etc.)
        for forbidden in ["is pathogenic", "is benign", "causes disease", "is diagnostic"]:
            assert forbidden not in text, f"Found clinical assertion: {forbidden!r}"


class TestVariantAgentCLI:
    def test_cli_variant_agent_dry_run(self, tmp_path):
        from genomics_workflow_agent.cli import main
        argv = [
            "agent",
            "--input", str(tmp_path),
            "--workflow", "variant-qc",
            "--out", str(tmp_path / "out"),
        ]
        with patch.object(sys, "argv", ["genomics_workflow_agent"] + argv):
            rc = main()
        assert rc == 0
        assert (tmp_path / "out" / "variant_agent_report.json").exists()

    def test_cli_unsupported_workflow_returns_error(self, tmp_path):
        from genomics_workflow_agent.cli import cmd_agent
        import argparse
        args = argparse.Namespace(
            workflow="rnaseq",
            input=str(tmp_path),
            out=str(tmp_path / "out"),
            execute=False,
        )
        rc = cmd_agent(args)
        assert rc == 1
