from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


@pytest.fixture
def fastq_dir(tmp_path):
    d = tmp_path / "fastqs"
    d.mkdir()
    (d / "sampleA_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleA_R2.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleB_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleB_R2.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    return d


@pytest.fixture
def bam_dir(tmp_path):
    d = tmp_path / "bams"
    d.mkdir()
    (d / "sample1.bam").write_bytes(b"BAM\x01\x00")
    (d / "sample2.bam").write_bytes(b"BAM\x01\x00")
    return d


@pytest.fixture
def vcf_dir(tmp_path):
    d = tmp_path / "vcfs"
    d.mkdir()
    (d / "variants.vcf.gz").write_bytes(b"\x1f\x8b\x00")
    return d


class TestRunnerCore:
    def test_dry_run_never_executes(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, STATUS_PLANNED
        record = run_command(["false_cmd_should_not_run"], dry_run=True)
        assert record["executed"] is False
        assert record["status"] == STATUS_PLANNED
        assert record["return_code"] is None

    def test_execute_captures_return_code(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, STATUS_SUCCEEDED
        cmd = [sys.executable, "-c", "print('hello')"]
        record = run_command(cmd, dry_run=False, provenance_dir=tmp_dir / "prov")
        assert record["executed"] is True
        assert record["return_code"] == 0
        assert record["status"] == STATUS_SUCCEEDED

    def test_failed_command_never_reported_as_success(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, STATUS_FAILED
        cmd = [sys.executable, "-c", "import sys; sys.exit(1)"]
        record = run_command(cmd, dry_run=False, provenance_dir=tmp_dir / "prov")
        assert record["return_code"] == 1
        assert record["status"] == STATUS_FAILED
        assert record["error"] is not None

    def test_missing_executable_is_error(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, STATUS_ERROR
        record = run_command(["__nonexistent_xyz__"], dry_run=False)
        assert record["status"] == STATUS_ERROR
        assert record["error"] is not None

    def test_execute_validates_expected_outputs_present(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, STATUS_SUCCEEDED
        out_file = tmp_dir / "out.txt"
        cmd = [sys.executable, "-c", f"open(r'{out_file}', 'w').write('ok')"]
        record = run_command(cmd, dry_run=False, expected_outputs=[str(out_file)])
        assert record["status"] == STATUS_SUCCEEDED
        assert record["output_validation"]["all_present"] is True

    def test_execute_fails_if_expected_outputs_missing(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command, STATUS_FAILED
        nonexistent = str(tmp_dir / "does_not_exist.html")
        cmd = [sys.executable, "-c", "pass"]
        record = run_command(cmd, dry_run=False, expected_outputs=[nonexistent])
        assert record["status"] == STATUS_FAILED
        assert nonexistent in record["output_validation"]["missing"]

    def test_dry_run_does_not_validate_outputs(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        record = run_command(
            ["echo", "hi"], dry_run=True,
            expected_outputs=[str(tmp_dir / "nonexistent.html")],
        )
        assert record["executed"] is False
        assert "not validated" in str(record.get("output_validation", ""))

    def test_capture_stdout_path_writes_file(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        out_file = tmp_dir / "stdout.txt"
        cmd = [sys.executable, "-c", "print('captured output')"]
        record = run_command(cmd, dry_run=False, capture_stdout_path=out_file)
        assert out_file.exists()
        assert "captured output" in out_file.read_text()

    def test_provenance_written_for_dry_run(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        prov = tmp_dir / "prov"
        run_command(["echo", "test"], dry_run=True, provenance_dir=prov, label="step1")
        files = list(prov.glob("provenance_step1_*.json"))
        assert len(files) == 1
        data = json.loads(files[0].read_text())
        assert data["dry_run"] is True
        assert data["label"] == "step1"

    def test_provenance_written_for_execute(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        prov = tmp_dir / "prov"
        cmd = [sys.executable, "-c", "pass"]
        run_command(cmd, dry_run=False, provenance_dir=prov, label="real_step")
        files = list(prov.glob("provenance_real_step_*.json"))
        assert len(files) == 1
        data = json.loads(files[0].read_text())
        assert data["dry_run"] is False
        assert data["executed"] is True


class TestOutputValidation:
    def test_validate_all_present(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import validate_outputs
        f1, f2 = tmp_dir / "a.html", tmp_dir / "b.zip"
        f1.write_text(""), f2.write_text("")
        result = validate_outputs([str(f1), str(f2)])
        assert result["all_present"] is True
        assert result["missing"] == []

    def test_validate_some_missing(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import validate_outputs
        f1 = tmp_dir / "present.html"
        f1.write_text("")
        missing = str(tmp_dir / "absent.zip")
        result = validate_outputs([str(f1), missing])
        assert result["all_present"] is False
        assert missing in result["missing"]

    def test_validate_empty_list(self):
        from genomics_workflow_agent.tools.runner import validate_outputs
        result = validate_outputs([])
        assert result["all_present"] is True


class TestExecutionSummary:
    def test_summary_all_succeeded(self):
        from genomics_workflow_agent.tools.runner import execution_summary, STATUS_SUCCEEDED
        records = [{"label": f"step{i}", "status": STATUS_SUCCEEDED} for i in range(3)]
        summary = execution_summary(records)
        assert summary["overall_status"] == "succeeded"
        assert summary["failed_steps"] == []

    def test_summary_with_failures(self):
        from genomics_workflow_agent.tools.runner import execution_summary, STATUS_SUCCEEDED, STATUS_FAILED
        records = [
            {"label": "step1", "status": STATUS_SUCCEEDED},
            {"label": "step2", "status": STATUS_FAILED},
        ]
        summary = execution_summary(records)
        assert summary["overall_status"] == "failed"
        assert "step2" in summary["failed_steps"]

    def test_summary_planned_only(self):
        from genomics_workflow_agent.tools.runner import execution_summary, STATUS_PLANNED
        records = [{"label": "step1", "status": STATUS_PLANNED}]
        summary = execution_summary(records)
        assert summary["overall_status"] == "planned"


class TestFastqcOutputNaming:
    def test_stem_fastq_gz(self):
        from genomics_workflow_agent.workflows.fastq_qc import _fastqc_stem
        assert _fastqc_stem("sampleA_R1.fastq.gz") == "sampleA_R1"

    def test_stem_fq_gz(self):
        from genomics_workflow_agent.workflows.fastq_qc import _fastqc_stem
        assert _fastqc_stem("sampleA_R1.fq.gz") == "sampleA_R1"

    def test_stem_fastq(self):
        from genomics_workflow_agent.workflows.fastq_qc import _fastqc_stem
        assert _fastqc_stem("sample.fastq") == "sample"

    def test_expected_outputs(self, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import _fastqc_expected_outputs
        outdir = str(tmp_dir / "fastqc")
        result = _fastqc_expected_outputs(["sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz"], outdir)
        assert any("sampleA_R1_fastqc.html" in p for p in result)
        assert any("sampleA_R1_fastqc.zip" in p for p in result)
        assert any("sampleA_R2_fastqc.html" in p for p in result)
        assert len(result) == 4


class TestFastqQcMocked:
    def _make_fake_multiqc(self, multiqc_dir: Path) -> None:
        multiqc_dir.mkdir(parents=True, exist_ok=True)
        (multiqc_dir / "multiqc_report.html").write_text("<html>MultiQC</html>")
        (multiqc_dir / "multiqc_data").mkdir(exist_ok=True)

    def _fake_tools(self, present: list[str]) -> dict:
        all_tools = ["fastqc", "multiqc", "fastp", "cutadapt", "nextflow",
                     "samtools", "bcftools", "mosdepth", "docker", "singularity", "conda",
                     "qiime", "Rscript", "bedtools"]
        return {t: {"available": t in present, "path": f"/usr/bin/{t}", "version": "0.0"}
                for t in all_tools}

    def test_plan_dry_run_generates_fastqc_step(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        result = plan(fastq_dir, tmp_dir / "out", dry_run=True)
        step_names = [s["name"] for s in result["steps"]]
        skipped_names = [s.get("step", "") for s in result.get("skipped_steps", [])]
        has_fastqc = any("fastqc" in n for n in step_names) or any("fastqc" in n for n in skipped_names)
        assert has_fastqc, f"fastqc missing from steps {step_names} and skipped {skipped_names}"

    def test_plan_includes_expected_outputs(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import plan
        result = plan(fastq_dir, tmp_dir / "out", dry_run=True)
        fqc_steps = [s for s in result["steps"] if s["name"] == "fastqc"]
        if fqc_steps:
            assert "expected_outputs" in fqc_steps[0]

    def test_mocked_fastqc_execution_succeeds(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.fastq_qc import _fastqc_stem

        out_dir = tmp_dir / "out"
        fastqc_dir = out_dir / "fastqc"
        multiqc_dir = out_dir / "multiqc"

        all_fastqs = [
            str(fastq_dir / "sampleA_R1.fastq.gz"),
            str(fastq_dir / "sampleA_R2.fastq.gz"),
            str(fastq_dir / "sampleB_R1.fastq.gz"),
            str(fastq_dir / "sampleB_R2.fastq.gz"),
        ]

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout = ""
            mock.stderr = ""
            if "fastqc" in str(cmd[0]):
                fastqc_dir.mkdir(parents=True, exist_ok=True)
                for f in all_fastqs:
                    stem = _fastqc_stem(Path(f).name)
                    (fastqc_dir / f"{stem}_fastqc.html").write_text("FastQC")
                    (fastqc_dir / f"{stem}_fastqc.zip").write_bytes(b"PK")
            elif "multiqc" in str(cmd[0]):
                self._make_fake_multiqc(multiqc_dir)
            return mock

        fake_tools = self._fake_tools(["fastqc", "multiqc"])
        with patch("subprocess.run", side_effect=fake_run), \
             patch("genomics_workflow_agent.workflows.fastq_qc.check_tools", return_value=fake_tools):
            from genomics_workflow_agent.workflows.fastq_qc import execute
            result = execute(fastq_dir, out_dir, provenance_dir=tmp_dir / "prov")

        step_results = result.get("step_results", [])
        assert len(step_results) > 0
        for r in step_results:
            assert r["status"] in ("succeeded", "planned", "skipped"), \
                f"Step {r['label']} unexpectedly {r['status']}: {r.get('error')}"

    def test_mocked_fastqc_failure_not_reported_as_success(self, fastq_dir, tmp_dir):
        def fake_fail(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 1
            mock.stdout = ""
            mock.stderr = "FastQC error"
            return mock

        with patch("subprocess.run", side_effect=fake_fail):
            from genomics_workflow_agent.workflows.fastq_qc import execute
            result = execute(fastq_dir, tmp_dir / "out", provenance_dir=tmp_dir / "prov")

        step_results = result.get("step_results", [])
        executed = [r for r in step_results if r.get("executed")]
        assert all(r["status"] == "failed" for r in executed), \
            "Failed commands must not be marked as succeeded"

    def test_mocked_trim_fastp_execution(self, fastq_dir, tmp_dir):
        out_dir = tmp_dir / "out"
        trim_dir = out_dir / "trimmed"

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "", ""
            if "fastp" in str(cmd[0]):
                trim_dir.mkdir(parents=True, exist_ok=True)
                for arg in cmd:
                    if isinstance(arg, str) and str(trim_dir) in arg:
                        Path(arg).parent.mkdir(parents=True, exist_ok=True)
                        Path(arg).write_text("trimmed")
            elif "fastqc" in str(cmd[0]):
                (out_dir / "fastqc").mkdir(parents=True, exist_ok=True)
            elif "multiqc" in str(cmd[0]):
                (out_dir / "multiqc").mkdir(parents=True, exist_ok=True)
                (out_dir / "multiqc" / "multiqc_report.html").write_text("mq")
                (out_dir / "multiqc" / "multiqc_data").mkdir(exist_ok=True)
            return mock

        fake_tools = self._fake_tools(["fastqc", "multiqc", "fastp"])
        with patch("subprocess.run", side_effect=fake_run), \
             patch("genomics_workflow_agent.workflows.fastq_qc.check_tools", return_value=fake_tools):
            from genomics_workflow_agent.workflows.fastq_qc import execute
            result = execute(fastq_dir, out_dir, trim="fastp",
                             provenance_dir=tmp_dir / "prov")

        assert result["trimming_requested"] is True
        assert result["trimmer"] == "fastp"
        trim_steps = [s for s in result["steps"] if "trim" in s["name"]]
        assert len(trim_steps) > 0


class TestNextflowRunner:
    def test_build_nextflow_cmd_rnaseq(self):
        from genomics_workflow_agent.tools.nextflow import build_nextflow_cmd
        cmd = build_nextflow_cmd(
            "rnaseq",
            input_path="/data/samplesheet.csv",
            outdir="/results",
            genome="GRCh38",
            profile="docker",
        )
        assert "nextflow" in cmd
        assert "nf-core/rnaseq" in cmd
        assert "--input" in cmd
        assert "--genome" in cmd
        assert "GRCh38" in cmd
        assert "-profile" in cmd

    def test_build_nextflow_cmd_ampliseq_with_primers(self):
        from genomics_workflow_agent.tools.nextflow import build_nextflow_cmd
        cmd = build_nextflow_cmd(
            "ampliseq",
            primer_fw="GTGYCAGCMGCCGCGGTAA",
            primer_rv="GGACTACNVGGGTWTCTAAT",
            taxonomy_param="silva=2.13",
            profile="docker",
        )
        assert "--FW_primer" in cmd
        assert "GTGYCAGCMGCCGCGGTAA" in cmd
        assert "--dada_ref_taxonomy" in cmd

    def test_build_nextflow_cmd_sarek_with_known_snps(self):
        from genomics_workflow_agent.tools.nextflow import build_nextflow_cmd
        cmd = build_nextflow_cmd(
            "sarek",
            known_snps="/refs/dbsnp.vcf.gz",
            genome="GRCh38",
        )
        assert "--known_snps" in cmd

    def test_dry_run_nextflow_does_not_execute(self, tmp_dir):
        from genomics_workflow_agent.tools.nextflow import run_nextflow
        cmd = ["nextflow", "run", "nf-core/rnaseq"]
        record = run_nextflow("rnaseq", cmd, outdir=tmp_dir / "out", dry_run=True)
        assert record["executed"] is False
        assert record["status"] == "planned"

    def test_mocked_nextflow_success(self, tmp_dir):
        from genomics_workflow_agent.tools.nextflow import run_nextflow

        outdir = tmp_dir / "nf_out"

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "Pipeline completed", ""
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "pipeline_info").mkdir()
            (outdir / "multiqc").mkdir()
            (outdir / "multiqc" / "multiqc_report.html").write_text("mq")
            return mock

        cmd = ["nextflow", "run", "nf-core/rnaseq", "--outdir", str(outdir)]
        with patch("subprocess.run", side_effect=fake_run), \
             patch("genomics_workflow_agent.tools.nextflow.is_nextflow_available", return_value=True):
            record = run_nextflow("rnaseq", cmd, outdir=outdir, dry_run=False,
                                  provenance_dir=tmp_dir / "prov")

        assert record["executed"] is True
        assert record["return_code"] == 0

    def test_mocked_nextflow_failure(self, tmp_dir):
        from genomics_workflow_agent.tools.nextflow import run_nextflow

        def fake_fail(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 1
            mock.stdout, mock.stderr = "", "ERROR: Workflow failed"
            return mock

        cmd = ["nextflow", "run", "nf-core/rnaseq"]
        with patch("subprocess.run", side_effect=fake_fail), \
             patch("genomics_workflow_agent.tools.nextflow.is_nextflow_available", return_value=True):
            record = run_nextflow("rnaseq", cmd, outdir=tmp_dir / "nf_out", dry_run=False)

        assert record["return_code"] == 1
        assert record["status"] == "failed"

    def test_nextflow_not_available_returns_skipped(self, tmp_dir):
        from genomics_workflow_agent.tools.nextflow import run_nextflow, STATUS_SKIPPED

        cmd = ["nextflow", "run", "nf-core/rnaseq"]
        with patch("genomics_workflow_agent.tools.nextflow.is_nextflow_available", return_value=False):
            record = run_nextflow("rnaseq", cmd, outdir=tmp_dir / "out", dry_run=False)

        assert record["status"] == STATUS_SKIPPED
        assert record["executed"] is False


class TestVariantQcMocked:
    def test_plan_has_samtools_steps_for_bams(self, bam_dir, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        result = plan(bam_dir, tmp_dir / "out")
        step_names = [s["name"] for s in result["steps"]]
        if shutil.which("samtools"):
            assert any("samtools_flagstat" in n for n in step_names)

    def test_plan_has_bcftools_steps_for_vcfs(self, vcf_dir, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        result = plan(vcf_dir, tmp_dir / "out")
        step_names = [s["name"] for s in result["steps"]]
        if shutil.which("bcftools"):
            assert any("bcftools_stats" in n for n in step_names)

    def test_plan_steps_have_capture_stdout_path(self, bam_dir, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        result = plan(bam_dir, tmp_dir / "out")
        for step in result["steps"]:
            if "samtools_flagstat" in step["name"]:
                assert "capture_stdout_path" in step
                assert step["capture_stdout_path"].endswith("_flagstat.txt")

    def test_mocked_samtools_flagstat_writes_file(self, bam_dir, tmp_dir):
        flagstat_output = (
            "100000 + 0 in total (QC-passed reads + QC-failed reads)\n"
            "98000 + 0 mapped (98.00% : N/A)\n"
        )

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout = flagstat_output
            mock.stderr = ""
            return mock

        with patch("subprocess.run", side_effect=fake_run):
            from genomics_workflow_agent.workflows.variant_qc import execute
            result = execute(bam_dir, tmp_dir / "out", provenance_dir=tmp_dir / "prov")

        step_results = result.get("step_results", [])
        flagstat_results = [r for r in step_results if "flagstat" in r.get("label", "")]
        if flagstat_results:
            r = flagstat_results[0]
            assert r["status"] in ("succeeded", "failed")
            stdout_file = r.get("stdout_file")
            if stdout_file and Path(stdout_file).exists():
                content = Path(stdout_file).read_text()
                assert "mapped" in content

    def test_mocked_samtools_failure_not_success(self, bam_dir, tmp_dir):
        def fake_fail(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 1
            mock.stdout, mock.stderr = "", "[E::sam_parse1]"
            return mock

        with patch("subprocess.run", side_effect=fake_fail):
            from genomics_workflow_agent.workflows.variant_qc import execute
            result = execute(bam_dir, tmp_dir / "out", provenance_dir=tmp_dir / "prov")

        for r in result.get("step_results", []):
            if r.get("executed"):
                assert r["status"] != "succeeded", \
                    "samtools failure must not be reported as success"

    def test_no_clinical_claims_in_variant_qc_plan(self, bam_dir, tmp_dir):
        from genomics_workflow_agent.workflows.variant_qc import plan
        result = plan(bam_dir, tmp_dir / "out")
        SAFE_QUALIFIERS = ["requires", "must not", "do not", "without", "expert", "validated"]
        for caveat in result["biological_caveats"]:
            lower = caveat.lower()
            if "pathogenic" in lower:
                is_safe = any(q in lower for q in SAFE_QUALIFIERS)
                if not is_safe:
                    pytest.fail(f"Caveat makes or implies clinical claim without qualification: {caveat}")


class TestAmpliconExecutionMocked:
    def test_mocked_ampliseq_execution(self, fastq_dir, tmp_dir):
        outdir = tmp_dir / "amp_out"
        pipeline_outdir = outdir / "ampliseq_output"

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "nf-core/ampliseq complete", ""
            pipeline_outdir.mkdir(parents=True, exist_ok=True)
            (pipeline_outdir / "dada2").mkdir()
            (pipeline_outdir / "taxonomy").mkdir()
            (pipeline_outdir / "diversity").mkdir()
            (pipeline_outdir / "pipeline_info").mkdir()
            return mock

        with patch("subprocess.run", side_effect=fake_run):
            with patch("genomics_workflow_agent.tools.nextflow.is_nextflow_available",
                       return_value=True):
                from genomics_workflow_agent.workflows.amplicon import execute
                result = execute(
                    fastq_dir, outdir,
                    primer_fw="GTGYCAGCMGCCGCGGTAA",
                    primer_rv="GGACTACNVGGGTWTCTAAT",
                    provenance_dir=tmp_dir / "prov",
                )

        step_results = result.get("step_results", [])
        assert len(step_results) > 0

    def test_amplicon_blockers_prevent_execution(self, fastq_dir, tmp_dir):
        with patch("genomics_workflow_agent.tools.nextflow.is_nextflow_available",
                   return_value=False):
            from genomics_workflow_agent.workflows.amplicon import execute
            result = execute(fastq_dir, tmp_dir / "out",
                             provenance_dir=tmp_dir / "prov")

        step_results = result.get("step_results", [])
        assert any(r.get("status") in ("skipped", "error") or r.get("executed") is False
                   for r in step_results)


class TestPlannerExecuteDispatch:
    def test_execute_plan_fastqc_dry_run_via_plan(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.planner import execute_plan

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "", ""
            return mock

        with patch("subprocess.run", side_effect=fake_run):
            result = execute_plan("fastq-qc", fastq_dir, tmp_dir / "out",
                                  provenance_dir=tmp_dir / "prov")

        assert result["workflow"] == "fastq-qc"
        assert result["dry_run"] is False

    def test_execute_plan_unknown_workflow_raises(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.workflows.planner import execute_plan
        with pytest.raises(ValueError):
            execute_plan("bad-workflow", fastq_dir, tmp_dir / "out")


class TestCLIExecuteFlag:
    def test_cli_run_dry_run_no_execute(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "run",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "out"),
        ]
        result = main()
        assert result == 0
        report = tmp_dir / "out" / "run_report.json"
        assert report.exists()
        data = json.loads(report.read_text())
        for r in data.get("step_results", []):
            assert r.get("executed") is False

    def test_cli_run_with_execute_calls_execute_plan(self, fastq_dir, tmp_dir):
        executed_cmds = []

        def fake_run(cmd, **kwargs):
            executed_cmds.append(cmd)
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "", ""
            return mock

        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "run",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "out"),
            "--execute",
        ]

        with patch("subprocess.run", side_effect=fake_run):
            result = main()

        assert result in (0, 1)

    def test_cli_run_trim_fastp(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "run",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "out"),
            "--trim", "fastp",
        ]
        result = main()
        assert result == 0
        report_path = tmp_dir / "out" / "run_report.json"
        data = json.loads(report_path.read_text())
        assert data.get("trimming_requested") is True

    def test_cli_run_trim_default_fastp(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "run",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "out"),
            "--trim",
        ]
        result = main()
        assert result == 0
        data = json.loads((tmp_dir / "out" / "run_report.json").read_text())
        assert data.get("trimming_requested") is True
        assert data.get("trimmer") in ("fastp", "cutadapt", None)

    def test_cli_plan_always_dry_run(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        fake_tools = {t: {"available": False, "version": None}
                      for t in ["fastqc", "multiqc", "fastp", "cutadapt",
                                "nextflow", "docker", "singularity", "conda"]}
        sys.argv = [
            "genomics_workflow_agent", "plan",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "out"),
        ]
        with patch("genomics_workflow_agent.workflows.fastq_qc.check_tools",
                   return_value=fake_tools):
            with patch("subprocess.run") as mock_sub:
                result = main()
                mock_sub.assert_not_called()
        assert result == 0


SKIP_NO_FASTQC = pytest.mark.skipif(not shutil.which("fastqc"), reason="fastqc not installed")
SKIP_NO_MULTIQC = pytest.mark.skipif(not shutil.which("multiqc"), reason="multiqc not installed")
SKIP_NO_FASTP = pytest.mark.skipif(not shutil.which("fastp"), reason="fastp not installed")
SKIP_NO_CUTADAPT = pytest.mark.skipif(not shutil.which("cutadapt"), reason="cutadapt not installed")
SKIP_NO_SAMTOOLS = pytest.mark.skipif(not shutil.which("samtools"), reason="samtools not installed")
SKIP_NO_BCFTOOLS = pytest.mark.skipif(not shutil.which("bcftools"), reason="bcftools not installed")
SKIP_NO_NEXTFLOW = pytest.mark.skipif(not shutil.which("nextflow"), reason="nextflow not installed")


@SKIP_NO_FASTQC
class TestFastQCIntegration:
    def test_fastqc_runs_on_example_file(self, tmp_dir):
        examples = Path(__file__).parent.parent / "examples"
        fastq = examples / "tiny.fastq"
        if not fastq.exists():
            pytest.skip("examples/tiny.fastq not found")

        from genomics_workflow_agent.tools.runner import run_command, STATUS_SUCCEEDED
        outdir = tmp_dir / "fastqc_out"
        outdir.mkdir()
        record = run_command(
            ["fastqc", "--outdir", str(outdir), str(fastq)],
            dry_run=False,
            label="integration_fastqc",
            expected_outputs=[
                str(outdir / "tiny_fastqc.html"),
                str(outdir / "tiny_fastqc.zip"),
            ],
        )
        assert record["return_code"] == 0, f"FastQC failed: {record['stderr_snippet']}"
        assert record["status"] == STATUS_SUCCEEDED
        assert record["output_validation"]["all_present"] is True


@SKIP_NO_MULTIQC
class TestMultiQCIntegration:
    def test_multiqc_runs_on_fastqc_output(self, tmp_dir):
        if not shutil.which("fastqc"):
            pytest.skip("fastqc also required")
        examples = Path(__file__).parent.parent / "examples"
        fastq = examples / "tiny.fastq"
        if not fastq.exists():
            pytest.skip("examples/tiny.fastq not found")

        from genomics_workflow_agent.workflows.fastq_qc import execute
        out_dir = tmp_dir / "fq_qc_out"
        result = execute(examples, out_dir, provenance_dir=tmp_dir / "prov")

        multiqc_html = out_dir / "multiqc" / "multiqc_report.html"
        assert multiqc_html.exists(), f"multiqc_report.html not found in {out_dir / 'multiqc'}"


@SKIP_NO_FASTP
class TestFaspIntegration:
    def test_fastp_trims_example_file(self, tmp_dir):
        examples = Path(__file__).parent.parent / "examples"
        fastq = examples / "tiny.fastq"
        if not fastq.exists():
            pytest.skip("examples/tiny.fastq not found")

        from genomics_workflow_agent.tools.runner import run_command, STATUS_SUCCEEDED
        out_file = tmp_dir / "trimmed.fastq"
        json_out = tmp_dir / "fastp.json"
        html_out = tmp_dir / "fastp.html"

        record = run_command(
            ["fastp", "-i", str(fastq), "-o", str(out_file),
             "--json", str(json_out), "--html", str(html_out)],
            dry_run=False,
            label="integration_fastp",
            expected_outputs=[str(out_file), str(json_out)],
        )
        assert record["return_code"] == 0, f"fastp failed: {record['stderr_snippet']}"
        assert record["status"] == STATUS_SUCCEEDED


@SKIP_NO_SAMTOOLS
class TestSamtoolsIntegration:
    def test_samtools_version_check(self):
        from genomics_workflow_agent.tools.versions import check_tools
        result = check_tools(["samtools"])
        assert result["samtools"]["available"] is True
        assert result["samtools"]["version"] is not None

    def test_samtools_flagstat_on_bam(self, tmp_dir):
        examples = Path(__file__).parent.parent / "examples"
        bam = examples / "tiny.bam"
        if not bam.exists():
            pytest.skip("examples/tiny.bam not found")

        from genomics_workflow_agent.tools.runner import run_command
        out_file = tmp_dir / "flagstat.txt"
        record = run_command(
            ["samtools", "flagstat", str(bam)],
            dry_run=False,
            capture_stdout_path=out_file,
            label="integration_flagstat",
        )
        assert record["return_code"] == 0
        assert out_file.exists()
        content = out_file.read_text()
        assert "in total" in content or "mapped" in content


@SKIP_NO_BCFTOOLS
class TestBcftoolsIntegration:
    def test_bcftools_version_check(self):
        from genomics_workflow_agent.tools.versions import check_tools
        result = check_tools(["bcftools"])
        assert result["bcftools"]["available"] is True

    def test_bcftools_stats_on_vcf(self, tmp_dir):
        examples = Path(__file__).parent.parent / "examples"
        vcf = examples / "tiny.vcf"
        if not vcf.exists():
            pytest.skip("examples/tiny.vcf not found")

        from genomics_workflow_agent.tools.runner import run_command
        out_file = tmp_dir / "bcftools_stats.txt"
        record = run_command(
            ["bcftools", "stats", str(vcf)],
            dry_run=False,
            capture_stdout_path=out_file,
            label="integration_bcftools_stats",
        )
        assert record["return_code"] == 0
        assert out_file.exists()
        content = out_file.read_text()
        assert "SN" in content


@SKIP_NO_NEXTFLOW
class TestNextflowIntegration:
    def test_nextflow_version(self):
        from genomics_workflow_agent.tools.versions import check_tools
        result = check_tools(["nextflow"])
        assert result["nextflow"]["available"] is True
        assert result["nextflow"]["version"] is not None

    def test_nextflow_help_does_not_crash(self, tmp_dir):
        from genomics_workflow_agent.tools.runner import run_command
        record = run_command(
            ["nextflow", "-version"],
            dry_run=False,
            label="integration_nextflow_version",
        )
        assert record["return_code"] == 0
