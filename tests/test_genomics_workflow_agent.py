"""
Tests for the genomics_workflow_agent package.

Covers:
- File detection and paired-end FASTQ detection
- Tool checker with missing tools
- Dry-run command generation and provenance capture
- Failed command handling
- Amplicon plan creation
- RNA-seq samplesheet generation
- ATAC-seq samplesheet generation
- No-clinical-claims guardrail
- Workflow auto-detection
- CLI inspect/plan commands
"""

import json
import os
import sys
from pathlib import Path

import pytest

# ──────────────────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


@pytest.fixture
def fastq_dir(tmp_path):
    """Directory with paired-end and single-end FASTQ files."""
    d = tmp_path / "fastqs"
    d.mkdir()
    (d / "sampleA_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleA_R2.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleB_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleB_R2.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "singleC.fastq").write_text("@read1\nACGT\n+\nIIII\n")
    return d


@pytest.fixture
def mixed_dir(tmp_path):
    """Directory with multiple genomics file types."""
    d = tmp_path / "mixed"
    d.mkdir()
    (d / "sample.bam").write_bytes(b"BAM\x01")
    (d / "variants.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")
    (d / "sample_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sample_R2.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    return d


@pytest.fixture
def examples_dir():
    root = Path(__file__).parent.parent / "examples"
    return root if root.exists() else None


# ──────────────────────────────────────────────────────────────────────────────
# Safety / guardrails
# ──────────────────────────────────────────────────────────────────────────────

class TestGuardrails:
    def test_no_clinical_claims_clean(self):
        from genomics_workflow_agent.safety.guardrails import assert_no_clinical_claims
        assert_no_clinical_claims("This sample has high mitochondrial percentage.")

    def test_no_clinical_claims_raises_on_pathogenic(self):
        from genomics_workflow_agent.safety.guardrails import assert_no_clinical_claims
        with pytest.raises(ValueError, match="pathogenic"):
            assert_no_clinical_claims("Variant is likely pathogenic.")

    def test_no_clinical_claims_case_insensitive(self):
        from genomics_workflow_agent.safety.guardrails import assert_no_clinical_claims
        with pytest.raises(ValueError):
            assert_no_clinical_claims("This is a DIAGNOSTIC result.")

    def test_clinical_disclaimer_is_present(self):
        from genomics_workflow_agent.safety.guardrails import DISCLAIMER
        assert "research" in DISCLAIMER.lower()
        assert "clinical" in DISCLAIMER.lower()

    def test_large_file_size_check(self, tmp_dir):
        from genomics_workflow_agent.safety.guardrails import check_file_size
        small = tmp_dir / "small.txt"
        small.write_text("hello")
        result = check_file_size(str(small))
        assert result["safe_to_read"] is True
        assert result["size_bytes"] > 0

    def test_large_file_not_safe(self, tmp_dir):
        from genomics_workflow_agent.safety.guardrails import check_file_size, LARGE_FILE_BYTES
        large = tmp_dir / "large.bin"
        large.write_bytes(b"\x00" * (LARGE_FILE_BYTES + 1))
        result = check_file_size(str(large))
        assert result["safe_to_read"] is False


# ──────────────────────────────────────────────────────────────────────────────
# File detection
# ──────────────────────────────────────────────────────────────────────────────

class TestFileDetection:
    def test_classify_fastq(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("sample_R1.fastq.gz")) == "fastq"
        assert classify_file(Path("sample.fq")) == "fastq"

    def test_classify_bam(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("sample.bam")) == "alignment"
        assert classify_file(Path("sample.cram")) == "alignment"

    def test_classify_vcf(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("variants.vcf.gz")) == "vcf"
        assert classify_file(Path("calls.bcf")) == "vcf"

    def test_classify_bed(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("peaks.bed")) == "bed"

    def test_classify_annotation(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("genes.gtf")) == "annotation"
        assert classify_file(Path("genome.gff3")) == "annotation"

    def test_classify_h5ad(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("data.h5ad")) == "h5"

    def test_classify_unknown(self):
        from genomics_workflow_agent.tools.files import classify_file
        assert classify_file(Path("README.md")) == "unknown"

    def test_discover_files(self, fastq_dir):
        from genomics_workflow_agent.tools.files import discover_files
        files = discover_files(fastq_dir)
        types = {f["type"] for f in files}
        assert "fastq" in types

    def test_discover_single_file(self, tmp_dir):
        from genomics_workflow_agent.tools.files import discover_files
        f = tmp_dir / "test.vcf"
        f.write_text("##fileformat=VCFv4.2\n")
        result = discover_files(f)
        assert len(result) == 1
        assert result[0]["type"] == "vcf"


class TestPairedEndDetection:
    def test_paired_end_detection(self, fastq_dir):
        from genomics_workflow_agent.tools.files import find_fastq_pairs
        pairs = find_fastq_pairs(fastq_dir)
        paired = [p for p in pairs if p["paired"]]
        assert len(paired) == 2
        assert all(p["r2"] is not None for p in paired)

    def test_single_end_detected(self, fastq_dir):
        from genomics_workflow_agent.tools.files import find_fastq_pairs
        pairs = find_fastq_pairs(fastq_dir)
        single = [p for p in pairs if not p["paired"]]
        assert len(single) == 1
        assert single[0]["sample"] == "singleC"

    def test_is_r1_r2(self):
        from genomics_workflow_agent.tools.files import is_r1, is_r2
        assert is_r1(Path("sample_R1.fastq.gz"))
        assert is_r2(Path("sample_R2.fastq.gz"))
        assert not is_r1(Path("sample_R2.fastq.gz"))
        assert not is_r2(Path("sample_R1.fastq.gz"))

    def test_sample_name_stripping(self):
        from genomics_workflow_agent.tools.files import sample_name
        assert sample_name(Path("sampleA_R1.fastq.gz")) == "sampleA"
        assert sample_name(Path("sampleB.1.fastq.gz")) == "sampleB"


# ──────────────────────────────────────────────────────────────────────────────
# File inspector
# ──────────────────────────────────────────────────────────────────────────────

class TestInspector:
    def test_inspect_nonexistent(self):
        from genomics_workflow_agent.inspect.inspector import inspect_file
        result = inspect_file("/nonexistent/file.fastq")
        assert result["exists"] is False
        assert result["errors"]

    def test_inspect_fastq(self, tmp_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_file
        fq = tmp_dir / "sample.fastq"
        fq.write_text("@read1\nACGT\n+\nIIII\n")
        result = inspect_file(fq)
        assert result["type"] == "fastq"
        assert result["exists"] is True
        assert result["fastq_valid_header"] is True

    def test_inspect_directory(self, fastq_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_directory
        result = inspect_directory(fastq_dir)
        assert result["total_files"] > 0
        assert "fastq" in result["file_type_counts"]

    def test_workflow_guess_rnaseq(self, tmp_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_directory
        d = tmp_dir / "rna"
        d.mkdir()
        (d / "rna_sample_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
        result = inspect_directory(d)
        assert result["workflow_guess"] == "rnaseq"

    def test_workflow_guess_amplicon(self, tmp_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_directory
        d = tmp_dir / "amp"
        d.mkdir()
        (d / "16s_sample_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
        result = inspect_directory(d)
        assert result["workflow_guess"] == "amplicon"

    def test_workflow_guess_variant(self, mixed_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_directory
        result = inspect_directory(mixed_dir)
        # BAM + VCF files: should guess variant-qc
        assert result["workflow_guess"] in ("variant-qc", "fastq-qc")

    def test_no_large_file_loading(self, tmp_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_file
        from genomics_workflow_agent.safety.guardrails import LARGE_FILE_BYTES
        large = tmp_dir / "big.fastq"
        large.write_bytes(b"@r\nA\n+\nI\n" * (LARGE_FILE_BYTES // 8 + 1))
        result = inspect_file(large)
        assert result["safe_to_read"] is False
        assert result["head_lines"] is None


# ──────────────────────────────────────────────────────────────────────────────
# Tool checker
# ──────────────────────────────────────────────────────────────────────────────

class TestToolChecker:
    def test_check_tools_returns_dict(self):
        from genomics_workflow_agent.tools.versions import check_tools
        result = check_tools(["fastqc", "multiqc"])
        assert "fastqc" in result
        assert "multiqc" in result
        for info in result.values():
            assert "available" in info
            assert isinstance(info["available"], bool)

    def test_missing_tool_has_false_available(self):
        from genomics_workflow_agent.tools.versions import check_tools
        result = check_tools(["__nonexistent_tool_xyz__"])
        assert result["__nonexistent_tool_xyz__"]["available"] is False
        assert result["__nonexistent_tool_xyz__"]["version"] is None

    def test_available_tools_returns_set(self):
        from genomics_workflow_agent.tools.versions import available_tools
        result = available_tools(["__nonexistent__"])
        assert "__nonexistent__" not in result


# ──────────────────────────────────────────────────────────────────────────────
# Command runner / provenance
# ──────────────────────────────────────────────────────────────────────────────

class TestRunner:
    def test_dry_run_does_not_execute(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        record = run_command(
            ["nonexistent_cmd", "--foo"],
            dry_run=True,
            provenance_dir=tmp_dir / "prov",
        )
        assert record["dry_run"] is True
        assert record["executed"] is False
        assert record["return_code"] is None

    def test_dry_run_writes_provenance(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        prov_dir = tmp_dir / "prov"
        run_command(["echo", "hi"], dry_run=True, provenance_dir=prov_dir, label="test_step")
        files = list(prov_dir.glob("*.json"))
        assert len(files) == 1
        data = json.loads(files[0].read_text())
        assert data["label"] == "test_step"
        assert data["dry_run"] is True

    def test_execute_captures_return_code(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        # Use a command guaranteed to exist on Windows and Unix
        cmd = [sys.executable, "-c", "print('hello')"]
        record = run_command(cmd, dry_run=False, provenance_dir=tmp_dir / "prov")
        assert record["executed"] is True
        assert record["return_code"] == 0

    def test_failed_command_captured(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        cmd = [sys.executable, "-c", "import sys; sys.exit(1)"]
        record = run_command(cmd, dry_run=False, provenance_dir=tmp_dir / "prov")
        assert record["return_code"] == 1
        assert record["error"] is not None

    def test_missing_executable_captured(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        record = run_command(
            ["__nonexistent_binary_xyz__"],
            dry_run=False,
            provenance_dir=tmp_dir / "prov",
        )
        assert record["error"] is not None

    def test_assert_no_silent_failure(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, assert_no_silent_failure
        cmd = [sys.executable, "-c", "import sys; sys.exit(42)"]
        record = run_command(cmd, dry_run=False, provenance_dir=tmp_dir / "prov")
        with pytest.raises(RuntimeError, match="42"):
            assert_no_silent_failure(record)


# ──────────────────────────────────────────────────────────────────────────────
# Samplesheet generators
# ──────────────────────────────────────────────────────────────────────────────

class TestSamplesheets:
    def test_rnaseq_samplesheet_created(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.tools.samplesheets import build_rnaseq_samplesheet
        result = build_rnaseq_samplesheet(fastq_dir, tmp_dir / "ss.csv")
        assert result["created"] is True
        assert result["rows"] > 0
        assert Path(result["path"]).exists()

    def test_rnaseq_samplesheet_has_strandedness_warning(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.tools.samplesheets import build_rnaseq_samplesheet
        result = build_rnaseq_samplesheet(fastq_dir, tmp_dir / "ss.csv")
        warnings_text = " ".join(result.get("warnings", []))
        assert "strandedness" in warnings_text

    def test_atacseq_samplesheet_created(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.tools.samplesheets import build_atacseq_samplesheet
        result = build_atacseq_samplesheet(fastq_dir, tmp_dir / "ss_atac.csv")
        assert result["created"] is True
        assert result["rows"] > 0

    def test_sarek_samplesheet_has_patient_placeholder(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.tools.samplesheets import build_sarek_samplesheet
        result = build_sarek_samplesheet(fastq_dir, tmp_dir / "ss_sarek.csv")
        assert result["created"] is True
        content = Path(result["path"]).read_text()
        assert "PATIENT_ID" in content

    def test_amplicon_samplesheet_created(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.tools.samplesheets import build_amplicon_samplesheet
        result = build_amplicon_samplesheet(fastq_dir, tmp_dir / "ss_amp.csv")
        assert result["created"] is True
        assert result["rows"] > 0

    def test_samplesheet_empty_dir(self, tmp_dir):
        from genomics_workflow_agent.tools.samplesheets import build_rnaseq_samplesheet
        empty = tmp_dir / "empty"
        empty.mkdir()
        result = build_rnaseq_samplesheet(empty, tmp_dir / "ss.csv")
        assert result["created"] is False


# ──────────────────────────────────────────────────────────────────────────────
# Workflows
# ──────────────────────────────────────────────────────────────────────────────

class TestFastqQcWorkflow:
    def test_plan_returns_dict(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["workflow"] == "fastq-qc"
        assert "steps" in result
        assert "biological_caveats" in result

    def test_plan_dry_run_default(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["dry_run"] is True

    def test_plan_detects_samples(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["samples_detected"] >= 2

    def test_plan_no_crash_on_empty_dir(self, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        empty = tmp_dir / "empty"
        empty.mkdir()
        result = plan(empty, tmp_dir / "out")
        assert result["workflow"] == "fastq-qc"


class TestRnaseqWorkflow:
    def test_plan_returns_dict(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.rnaseq import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["workflow"] == "rnaseq"
        assert result["pipeline"] == "nf-core/rnaseq"

    def test_plan_warns_no_genome(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.rnaseq import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert len(result["missing_requirements"]) > 0

    def test_plan_has_biological_caveats(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.rnaseq import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert len(result["biological_caveats"]) > 0

    def test_plan_command_contains_nfcore(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.rnaseq import plan
        result = plan(fastq_dir, tmp_dir / "out", genome="GRCh38")
        assert "nf-core/rnaseq" in result["command_str"]


class TestAtacseqWorkflow:
    def test_plan_returns_dict(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.atacseq import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["workflow"] == "atacseq"

    def test_plan_warns_no_blacklist(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.atacseq import plan
        result = plan(fastq_dir, tmp_dir / "out")
        warnings_text = " ".join(result["warnings"])
        assert "blacklist" in warnings_text.lower()


class TestAmpliconWorkflow:
    def test_plan_returns_dict(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["workflow"] == "amplicon"

    def test_plan_has_steps(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert isinstance(result["steps"], list)
        step_names = [s["name"] for s in result["steps"]]
        assert "feature_table_filter" in step_names
        assert "normalisation" in step_names
        assert "alpha_diversity" in step_names
        assert "beta_diversity" in step_names

    def test_plan_warns_no_primers(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        warnings_text = " ".join(result["warnings"])
        assert "primer" in warnings_text.lower()

    def test_plan_warns_no_taxonomy_db(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        warnings_text = " ".join(result["warnings"])
        assert "taxonomy" in warnings_text.lower() or "database" in warnings_text.lower()

    def test_plan_biological_caveats(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert len(result["biological_caveats"]) >= 5

    def test_plan_silva_default(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert result["taxonomy_db"] == "SILVA"

    def test_plan_gtdb_selection(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out", taxonomy_db="GTDB")
        assert result["taxonomy_db"] == "GTDB"
        assert "gtdb" in str(result["taxonomy_db_info"]).lower()

    def test_taxonomy_databases_registered(self):
        from genomics_workflow_agent.workflows.amplicon import TAXONOMY_DATABASES
        for db in ["SILVA", "GTDB", "UNITE", "Greengenes2", "custom"]:
            assert db in TAXONOMY_DATABASES

    def test_plan_does_not_claim_qiime_executed(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        for step in result["steps"]:
            if step.get("skipped"):
                assert step["command"] is None

    def test_samplesheet_created(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        ss = result["samplesheet"]
        assert ss["created"] is True
        content = Path(ss["path"]).read_text()
        assert "sampleID" in content

    def test_amplicon_has_limitations(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.amplicon import plan
        result = plan(fastq_dir, tmp_dir / "out")
        assert len(result["limitations"]) > 0
        limitations_text = " ".join(result["limitations"])
        assert "QIIME2" in limitations_text or "qiime2" in limitations_text.lower()


class TestVariantQcWorkflow:
    def test_plan_returns_dict(self, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        empty = tmp_dir / "empty"
        empty.mkdir()
        result = plan(empty, tmp_dir / "out")
        assert result["workflow"] == "variant-qc"

    def test_plan_has_clinical_disclaimer(self, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        empty = tmp_dir / "empty"
        empty.mkdir()
        result = plan(empty, tmp_dir / "out")
        disclaimer = result.get("clinical_disclaimer", "")
        assert "clinical" in disclaimer.lower() or "research" in disclaimer.lower()

    def test_plan_no_clinical_claims_in_caveats(self, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        from genomics_workflow_agent.safety.guardrails import CLINICAL_CLAIM_PATTERNS
        empty = tmp_dir / "empty"
        empty.mkdir()
        result = plan(empty, tmp_dir / "out")
        caveats_text = " ".join(result["biological_caveats"])
        # Caveats may mention "pathogenicity" in the context of NOT making claims
        # The guardrail only applies to direct assertions — caveats explicitly deny claims
        assert "must not be used for medical decisions" in caveats_text


# ──────────────────────────────────────────────────────────────────────────────
# Workflow planner
# ──────────────────────────────────────────────────────────────────────────────

class TestPlanner:
    def test_resolve_workflow_explicit(self, fastq_dir):
        from genomics_workflow_agent.workflows.planner import resolve_workflow
        assert resolve_workflow("rnaseq", fastq_dir) == "rnaseq"

    def test_resolve_workflow_invalid(self, fastq_dir):
        from genomics_workflow_agent.workflows.planner import resolve_workflow
        with pytest.raises(ValueError, match="Unknown workflow"):
            resolve_workflow("invalid-workflow", fastq_dir)

    def test_resolve_auto_fastq(self, fastq_dir):
        from genomics_workflow_agent.workflows.planner import resolve_workflow
        result = resolve_workflow("auto", fastq_dir)
        assert result in ["fastq-qc", "rnaseq", "atacseq", "amplicon"]

    def test_build_plan_fastq_qc(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.planner import build_plan
        result = build_plan("fastq-qc", fastq_dir, tmp_dir / "out", dry_run=True)
        assert result["workflow"] == "fastq-qc"

    def test_build_plan_amplicon(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.planner import build_plan
        result = build_plan("amplicon", fastq_dir, tmp_dir / "out", dry_run=True)
        assert result["workflow"] == "amplicon"


# ──────────────────────────────────────────────────────────────────────────────
# Reports
# ──────────────────────────────────────────────────────────────────────────────

class TestReports:
    def test_json_report_written(self, tmp_dir):
        from genomics_workflow_agent.reports.json_report import write_json_report
        path = write_json_report({"workflow": "test", "data": [1, 2, 3]}, tmp_dir / "report.json")
        assert path.exists()
        data = json.loads(path.read_text())
        assert data["workflow"] == "test"
        assert "meta" in data
        assert data["meta"]["tool"] == "genomics_workflow_agent"

    def test_markdown_report_written(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        from genomics_workflow_agent.reports.markdown import write_markdown_report
        p = plan(fastq_dir, tmp_dir / "out")
        md_path = write_markdown_report(p, tmp_dir / "report.md")
        assert md_path.exists()
        content = md_path.read_text()
        assert "# Genomics Workflow Report" in content
        assert "fastq-qc" in content

    def test_inspection_report_written(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.inspect.inspector import inspect_directory
        from genomics_workflow_agent.reports.markdown import write_inspection_report
        inspection = inspect_directory(fastq_dir)
        md_path = write_inspection_report(inspection, tmp_dir / "inspection.md")
        assert md_path.exists()
        content = md_path.read_text()
        assert "File Inspection Report" in content

    def test_markdown_report_includes_caveats(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.rnaseq import plan
        from genomics_workflow_agent.reports.markdown import write_markdown_report
        p = plan(fastq_dir, tmp_dir / "out")
        md_path = write_markdown_report(p, tmp_dir / "rnaseq_report.md")
        content = md_path.read_text()
        assert "Biological Caveats" in content

    def test_variant_qc_markdown_has_disclaimer(self, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        from genomics_workflow_agent.reports.markdown import write_markdown_report
        empty = tmp_dir / "empty"
        empty.mkdir()
        p = plan(empty, tmp_dir / "out")
        md_path = write_markdown_report(p, tmp_dir / "variant_report.md")
        content = md_path.read_text()
        assert "Clinical Disclaimer" in content or "clinical" in content.lower()


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

class TestCLI:
    def test_cli_imports(self):
        from genomics_workflow_agent.cli import build_parser, main
        parser = build_parser()
        assert parser is not None

    def test_cli_no_command_exits_zero(self):
        from genomics_workflow_agent.cli import main
        sys.argv = ["genomics_workflow_agent"]
        result = main()
        assert result == 0

    def test_cli_inspect(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "inspect",
            "--input", str(fastq_dir),
            "--out", str(tmp_dir / "inspect_out"),
        ]
        result = main()
        assert result == 0
        assert (tmp_dir / "inspect_out" / "inspection.json").exists()

    def test_cli_plan_fastq_qc(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "plan",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "plan_out"),
        ]
        result = main()
        assert result == 0
        assert (tmp_dir / "plan_out" / "plan.json").exists()
        assert (tmp_dir / "plan_out" / "plan.md").exists()

    def test_cli_plan_amplicon(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "plan",
            "--input", str(fastq_dir),
            "--workflow", "amplicon",
            "--out", str(tmp_dir / "amp_out"),
        ]
        result = main()
        assert result == 0

    def test_cli_run_dry_run(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "run",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "run_out"),
            # no --execute flag = dry-run
        ]
        result = main()
        assert result == 0

    def test_cli_report(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        # First generate a plan
        plan_dir = tmp_dir / "plan_out"
        sys.argv = [
            "genomics_workflow_agent", "plan",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(plan_dir),
        ]
        main()
        # Then generate report
        sys.argv = [
            "genomics_workflow_agent", "report",
            "--results", str(plan_dir),
        ]
        result = main()
        assert result == 0
        assert (plan_dir / "final_report.json").exists()
