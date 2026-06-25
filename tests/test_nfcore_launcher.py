"""
Tests for tools/nfcore_launcher.py.

All tests run without Nextflow, Docker, Singularity, or any network access.
Pure-function tests import directly; CLI tests invoke via subprocess.
"""

import json
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))

from nfcore_launcher import (
    BIOLOGICAL_CAVEATS,
    build_atacseq_samplesheet,
    build_nfcore_command,
    build_rnaseq_samplesheet,
    build_sarek_samplesheet,
    detect_executors,
    find_fastq_files,
    parse_multiqc_output,
)

REPO_ROOT = Path(__file__).parent.parent
TOOLS_DIR = REPO_ROOT / "tools"
EXAMPLES_DIR = REPO_ROOT / "examples"
LOG = logging.getLogger("test_nfcore_launcher")


def run_launcher(*args, **kwargs):
    cmd = [sys.executable, str(TOOLS_DIR / "nfcore_launcher.py")] + list(args)
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60, **kwargs)


def synthetic_fastqs(tmp_path, pairs=("sample1",)):
    """Write minimal synthetic paired-end FASTQ files and return the directory."""
    for name in pairs:
        (tmp_path / f"{name}_R1.fastq").write_text("@r1\nACGT\n+\nIIII\n")
        (tmp_path / f"{name}_R2.fastq").write_text("@r2\nACGT\n+\nIIII\n")
    return tmp_path


# ──────────────────────────────────────────────────────────
# CLI: help and dry-run basics
# ──────────────────────────────────────────────────────────

class TestCLIHelp:
    def test_help_exits_zero(self):
        r = run_launcher("--help")
        assert r.returncode == 0

    def test_help_mentions_supported_workflows(self):
        r = run_launcher("--help")
        output = r.stdout + r.stderr
        assert "rnaseq" in output
        assert "sarek" in output
        assert "atacseq" in output


class TestDryRunOutputs:
    def test_dry_run_exits_zero_without_nextflow(self, tmp_path):
        r = run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                         "--output-dir", str(tmp_path), "--dry-run")
        assert r.returncode == 0

    def test_json_written(self, tmp_path):
        run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        json_path = tmp_path / "nfcore_plan.json"
        assert json_path.exists(), "nfcore_plan.json not created"
        data = json.loads(json_path.read_text())
        assert data["workflow"] == "nf-core/rnaseq"
        assert data["mode"] == "dry_run"

    def test_markdown_written(self, tmp_path):
        run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        md_path = tmp_path / "nfcore_plan.md"
        assert md_path.exists(), "nfcore_plan.md not created"
        assert "nf-core/rnaseq" in md_path.read_text()

    def test_commands_sh_written(self, tmp_path):
        run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        sh_path = tmp_path / "commands.sh"
        assert sh_path.exists(), "commands.sh not created"
        assert "nextflow run nf-core/rnaseq" in sh_path.read_text()

    def test_missing_nextflow_is_blocker_not_crash(self, tmp_path):
        """Missing nextflow must appear in blockers in JSON, not cause a crash."""
        r = run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                         "--output-dir", str(tmp_path), "--dry-run")
        assert r.returncode == 0
        if shutil.which("nextflow") is None:
            data = json.loads((tmp_path / "nfcore_plan.json").read_text())
            assert len(data["blockers"]) > 0
            assert any("nextflow" in b.lower() for b in data["blockers"])


# ──────────────────────────────────────────────────────────
# CLI: --run refuses with blockers
# ──────────────────────────────────────────────────────────

class TestRunRefusesWithBlockers:
    def test_run_exits_nonzero_when_fasta_missing(self, tmp_path):
        """--run with a nonexistent --fasta is always a blocker."""
        r = run_launcher(
            "--workflow", "rnaseq",
            "--fasta", "/nonexistent/genome.fa",
            "--output-dir", str(tmp_path),
            "--run",
        )
        assert r.returncode != 0

    def test_run_nonzero_still_writes_json(self, tmp_path):
        run_launcher(
            "--workflow", "rnaseq",
            "--fasta", "/nonexistent/genome.fa",
            "--output-dir", str(tmp_path),
            "--run",
        )
        data = json.loads((tmp_path / "nfcore_plan.json").read_text())
        assert len(data["blockers"]) > 0
        assert data["run_result"]["executed"] is False


# ──────────────────────────────────────────────────────────
# Pure function: find_fastq_files
# ──────────────────────────────────────────────────────────

class TestFindFastqFiles:
    def test_finds_fastq_in_examples(self):
        files = find_fastq_files(EXAMPLES_DIR)
        assert any(f.suffix == ".fastq" or str(f).endswith(".fastq.gz") for f in files)

    def test_skips_non_fastq(self, tmp_path):
        (tmp_path / "data.csv").write_text("a,b")
        (tmp_path / "reads.fastq").write_text("@r1\nA\n+\nI\n")
        files = find_fastq_files(tmp_path)
        assert len(files) == 1
        assert files[0].name == "reads.fastq"

    def test_nonexistent_dir_returns_empty(self, tmp_path):
        files = find_fastq_files(tmp_path / "missing")
        assert files == []


# ──────────────────────────────────────────────────────────
# Pure function: build_rnaseq_samplesheet
# ──────────────────────────────────────────────────────────

class TestRnaseqSamplesheet:
    def test_pair_detection(self, tmp_path):
        synthetic_fastqs(tmp_path, ["mysample"])
        ss_path = tmp_path / "ss.csv"
        result = build_rnaseq_samplesheet(tmp_path, ss_path, LOG)
        assert result["created"] is True
        assert result["rows"] == 1
        content = ss_path.read_text()
        assert "mysample" in content

    def test_header_columns(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        build_rnaseq_samplesheet(tmp_path, ss_path, LOG)
        assert "sample,fastq_1,fastq_2,strandedness" in ss_path.read_text()

    def test_strandedness_default_is_auto(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        build_rnaseq_samplesheet(tmp_path, ss_path, LOG)
        assert "auto" in ss_path.read_text()

    def test_warns_about_strandedness(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        result = build_rnaseq_samplesheet(tmp_path, ss_path, LOG)
        assert result["created"] is True
        warnings_text = " ".join(result.get("warnings", [])).lower()
        assert "strandedness" in warnings_text

    def test_no_fastq_returns_not_created(self, tmp_path):
        ss_path = tmp_path / "ss.csv"
        result = build_rnaseq_samplesheet(tmp_path, ss_path, LOG)
        assert result["created"] is False

    def test_from_examples_dir(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            ss_path = Path(td) / "ss.csv"
            result = build_rnaseq_samplesheet(EXAMPLES_DIR, ss_path, LOG)
            if result["created"]:
                content = ss_path.read_text()
                assert "sample,fastq_1,fastq_2,strandedness" in content


# ──────────────────────────────────────────────────────────
# Pure function: build_sarek_samplesheet
# ──────────────────────────────────────────────────────────

class TestSarekSamplesheet:
    def test_created(self, tmp_path):
        synthetic_fastqs(tmp_path, ["patient_sample"])
        ss_path = tmp_path / "ss.csv"
        result = build_sarek_samplesheet(tmp_path, ss_path, LOG)
        assert result["created"] is True

    def test_header_columns(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        build_sarek_samplesheet(tmp_path, ss_path, LOG)
        assert "patient,sex,status,sample,lane,fastq_1,fastq_2" in ss_path.read_text()

    def test_patient_id_placeholder(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        build_sarek_samplesheet(tmp_path, ss_path, LOG)
        assert "PATIENT_ID" in ss_path.read_text()

    def test_warns_tumor_normal(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        result = build_sarek_samplesheet(tmp_path, ss_path, LOG)
        warnings_text = " ".join(result.get("warnings", [])).lower()
        assert "tumor" in warnings_text or "normal" in warnings_text

    def test_warns_manual_review(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        result = build_sarek_samplesheet(tmp_path, ss_path, LOG)
        warnings_text = " ".join(result.get("warnings", [])).lower()
        assert "review" in warnings_text or "manual" in warnings_text


# ──────────────────────────────────────────────────────────
# Pure function: build_atacseq_samplesheet
# ──────────────────────────────────────────────────────────

class TestAtacseqSamplesheet:
    def test_created(self, tmp_path):
        synthetic_fastqs(tmp_path, ["atac_rep1"])
        ss_path = tmp_path / "ss.csv"
        result = build_atacseq_samplesheet(tmp_path, ss_path, LOG)
        assert result["created"] is True

    def test_header_columns(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        build_atacseq_samplesheet(tmp_path, ss_path, LOG)
        assert "sample,fastq_1,fastq_2,replicate" in ss_path.read_text()

    def test_warns_replicate(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        result = build_atacseq_samplesheet(tmp_path, ss_path, LOG)
        warnings_text = " ".join(result.get("warnings", [])).lower()
        assert "replicate" in warnings_text

    def test_warns_blacklist(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        result = build_atacseq_samplesheet(tmp_path, ss_path, LOG)
        warnings_text = " ".join(result.get("warnings", [])).lower()
        assert "blacklist" in warnings_text

    def test_warns_control(self, tmp_path):
        synthetic_fastqs(tmp_path, ["s1"])
        ss_path = tmp_path / "ss.csv"
        result = build_atacseq_samplesheet(tmp_path, ss_path, LOG)
        warnings_text = " ".join(result.get("warnings", [])).lower()
        assert "control" in warnings_text or "replicate" in warnings_text


# ──────────────────────────────────────────────────────────
# Pure function: parse_multiqc_output
# ──────────────────────────────────────────────────────────

class TestMultiQCParser:
    def test_missing_dir_returns_skipped(self, tmp_path):
        result = parse_multiqc_output(tmp_path / "nonexistent", LOG)
        assert result["multiqc_report_html"] is None
        assert len(result["skipped"]) > 0

    def test_empty_dir_no_html(self, tmp_path):
        result = parse_multiqc_output(tmp_path, LOG)
        assert result["multiqc_report_html"] is None

    def test_finds_report_html(self, tmp_path):
        (tmp_path / "multiqc_report.html").write_text("<html></html>")
        result = parse_multiqc_output(tmp_path, LOG)
        assert result["multiqc_report_html"] is not None

    def test_parses_general_stats(self, tmp_path):
        data_dir = tmp_path / "multiqc_data"
        data_dir.mkdir()
        (data_dir / "multiqc_general_stats.txt").write_text(
            "Sample\tmetric_a\tmetric_b\nSAMPLE1\t1.0\t2.0\n"
        )
        result = parse_multiqc_output(tmp_path, LOG)
        assert result["general_stats_present"] is True
        assert result["summary"]["general_stats_rows"] == 1

    def test_missing_html_appears_in_skipped(self, tmp_path):
        result = parse_multiqc_output(tmp_path, LOG)
        assert any("multiqc_report.html" in s for s in result["skipped"])


# ──────────────────────────────────────────────────────────
# Pure function: build_nfcore_command
# ──────────────────────────────────────────────────────────

class TestBuildNfcoreCommand:
    def _args(self, **kwargs):
        class Args:
            pass
        a = Args()
        a.genome = kwargs.get("genome", "GRCh38")
        a.fasta = kwargs.get("fasta", None)
        a.gtf = kwargs.get("gtf", None)
        a.bed = kwargs.get("bed", None)
        a.profile = kwargs.get("profile", "docker")
        a.output_dir = kwargs.get("output_dir", "/tmp/out")
        a.samplesheet = kwargs.get("samplesheet", None)
        a.max_cpus = kwargs.get("max_cpus", None)
        a.max_memory = kwargs.get("max_memory", None)
        a.max_time = kwargs.get("max_time", None)
        a.resume = kwargs.get("resume", False)
        a.extra_args = kwargs.get("extra_args", None)
        return a

    def test_starts_with_nextflow_run(self):
        cmd = build_nfcore_command("rnaseq", self._args(), "/tmp/ss.csv")
        assert cmd[0] == "nextflow"
        assert cmd[1] == "run"
        assert cmd[2] == "nf-core/rnaseq"

    def test_includes_genome(self):
        cmd = build_nfcore_command("rnaseq", self._args(genome="GRCh38"), None)
        assert "--genome" in cmd
        assert "GRCh38" in cmd

    def test_includes_profile(self):
        cmd = build_nfcore_command("rnaseq", self._args(profile="singularity"), None)
        assert "-profile" in cmd
        assert "singularity" in cmd

    def test_includes_samplesheet(self):
        cmd = build_nfcore_command("rnaseq", self._args(), "/path/to/ss.csv")
        assert "--input" in cmd
        idx = cmd.index("--input")
        assert cmd[idx + 1] == "/path/to/ss.csv"


# ──────────────────────────────────────────────────────────
# Biological caveats and no-clinical-claims
# ──────────────────────────────────────────────────────────

class TestBiologicalCaveats:
    def test_rnaseq_caveats_section_in_markdown(self, tmp_path):
        run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        md = (tmp_path / "nfcore_plan.md").read_text()
        assert "Biological and Experimental Caveats" in md

    def test_rnaseq_strandedness_caveat(self, tmp_path):
        run_launcher("--workflow", "rnaseq", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        md = (tmp_path / "nfcore_plan.md").read_text().lower()
        assert "strandedness" in md

    def test_sarek_caveats_in_json(self, tmp_path):
        run_launcher("--workflow", "sarek", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        data = json.loads((tmp_path / "nfcore_plan.json").read_text())
        assert len(data["biological_caveats"]) > 0

    def test_sarek_no_clinical_claims(self, tmp_path):
        run_launcher("--workflow", "sarek", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        data = json.loads((tmp_path / "nfcore_plan.json").read_text())
        caveats = " ".join(data["biological_caveats"]).lower()
        # Must warn about clinical limitations
        assert "clinical" in caveats
        # Must not make diagnostic claims
        assert "diagnosis" not in caveats
        assert "diagnose" not in caveats
        assert "medical decision" not in caveats or "must not" in caveats

    def test_atacseq_peak_caveat_in_json(self, tmp_path):
        run_launcher("--workflow", "atacseq", "--genome", "GRCh38",
                     "--output-dir", str(tmp_path), "--dry-run")
        data = json.loads((tmp_path / "nfcore_plan.json").read_text())
        caveats = " ".join(data["biological_caveats"]).lower()
        assert "peak" in caveats or "frip" in caveats

    def test_biological_caveats_constant_all_workflows(self):
        for workflow in ["rnaseq", "sarek", "atacseq"]:
            assert workflow in BIOLOGICAL_CAVEATS
            assert len(BIOLOGICAL_CAVEATS[workflow]) >= 3


# ──────────────────────────────────────────────────────────
# Detect executors (pure path check — no actual subprocesses required)
# ──────────────────────────────────────────────────────────

class TestDetectExecutors:
    def test_returns_all_expected_tools(self):
        result = detect_executors()
        for tool in ["nextflow", "docker", "singularity", "apptainer", "conda"]:
            assert tool in result

    def test_available_field_is_bool(self):
        result = detect_executors()
        for info in result.values():
            assert isinstance(info["available"], bool)
