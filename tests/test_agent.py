from __future__ import annotations

import io
import json
import sys
import zipfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

FIXTURES = Path(__file__).parent / "fixtures" / "fastqc"


def _load_fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def _make_fastqc_zip(tmp_path: Path, sample: str, txt_content: str) -> Path:
    """Create a minimal FastQC zip at tmp_path/<sample>_fastqc.zip."""
    zip_path = tmp_path / f"{sample}_fastqc.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr(f"{sample}_fastqc/fastqc_data.txt", txt_content)
    return zip_path


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path


@pytest.fixture
def all_pass_txt() -> str:
    return _load_fixture("all_pass_fastqc_data.txt")


@pytest.fixture
def adapter_warn_txt() -> str:
    return _load_fixture("adapter_warn_fastqc_data.txt")


@pytest.fixture
def qual_fail_txt() -> str:
    return _load_fixture("qual_fail_fastqc_data.txt")


@pytest.fixture
def malformed_txt() -> str:
    return _load_fixture("malformed_fastqc_data.txt")


@pytest.fixture
def fastq_dir(tmp_path):
    d = tmp_path / "fastqs"
    d.mkdir()
    (d / "sampleA_R1.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    (d / "sampleA_R2.fastq.gz").write_bytes(b"\x1f\x8b\x00")
    return d


class TestFastqcParser:
    def test_parse_all_pass(self, all_pass_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        result = parse_fastqc_txt(all_pass_txt, sample="sampleA_R1")
        assert result["parse_ok"] is True
        assert "Per base sequence quality" in result["modules"]
        assert result["modules"]["Per base sequence quality"]["status"] == "pass"
        assert result["modules"]["Adapter Content"]["status"] == "pass"
        assert "Per base sequence quality" in result["summary"]["pass"]

    def test_parse_adapter_warn(self, adapter_warn_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        result = parse_fastqc_txt(adapter_warn_txt, sample="sampleB_R1")
        assert result["modules"]["Adapter Content"]["status"] == "warn"
        assert "Adapter Content" in result["summary"]["warn"]
        assert result["parse_ok"] is True

    def test_parse_qual_fail(self, qual_fail_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        result = parse_fastqc_txt(qual_fail_txt, sample="sampleC_R1")
        assert result["modules"]["Per base sequence quality"]["status"] == "fail"
        assert "Per base sequence quality" in result["summary"]["fail"]
        assert result["modules"]["Overrepresented sequences"]["status"] == "warn"

    def test_parse_malformed_does_not_crash(self, malformed_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        result = parse_fastqc_txt(malformed_txt, sample="broken")
        assert isinstance(result, dict)
        assert "modules" in result
        assert "errors" in result

    def test_missing_module_returns_none(self, all_pass_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt, module_status
        result = parse_fastqc_txt(all_pass_txt, sample="s")
        assert module_status(result, "Nonexistent Module") is None

    def test_parse_zip(self, tmp_dir, all_pass_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_zip
        zip_path = _make_fastqc_zip(tmp_dir, "sampleA_R1", all_pass_txt)
        result = parse_fastqc_zip(zip_path)
        assert result["parse_ok"] is True
        assert result["sample"] == "sampleA_R1"
        assert "Per base sequence quality" in result["modules"]

    def test_parse_zip_missing_data_file(self, tmp_dir):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_zip
        zip_path = tmp_dir / "empty_fastqc.zip"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("something_else.txt", "not fastqc data")
        result = parse_fastqc_zip(zip_path)
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_parse_zip_not_found(self, tmp_dir):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_zip
        result = parse_fastqc_zip(tmp_dir / "does_not_exist_fastqc.zip")
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_parse_zip_bad_zip(self, tmp_dir):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_zip
        bad_zip = tmp_dir / "corrupt_fastqc.zip"
        bad_zip.write_bytes(b"this is not a zip file at all")
        result = parse_fastqc_zip(bad_zip)
        assert result["parse_ok"] is False
        assert len(result["errors"]) > 0

    def test_parse_fastqc_dir(self, tmp_dir, all_pass_txt, adapter_warn_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_dir
        _make_fastqc_zip(tmp_dir, "sampleA_R1", all_pass_txt)
        _make_fastqc_zip(tmp_dir, "sampleB_R1", adapter_warn_txt)
        results = parse_fastqc_dir(tmp_dir)
        assert len(results) == 2
        assert all(r["parse_ok"] for r in results)

    def test_parse_fastqc_dir_empty(self, tmp_dir):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_dir
        results = parse_fastqc_dir(tmp_dir)
        assert results == []

    def test_parse_fastqc_dir_missing(self, tmp_dir):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_dir
        results = parse_fastqc_dir(tmp_dir / "does_not_exist")
        assert results == []


class TestMultiqcParser:
    def test_parse_missing_directory(self, tmp_dir):
        from genomics_workflow_agent.parsers.multiqc import parse_multiqc_output
        result = parse_multiqc_output(tmp_dir / "no_multiqc_here")
        assert isinstance(result, dict)
        assert "errors" in result
        assert result["parse_ok"] is False

    def test_parse_multiqc_fastqc_tsv(self, tmp_dir):
        from genomics_workflow_agent.parsers.multiqc import parse_multiqc_fastqc_tsv
        data_dir = tmp_dir / "multiqc_data"
        data_dir.mkdir()
        (data_dir / "multiqc_fastqc.txt").write_text(
            "Sample\tper_base_sequence_quality\tadapter_content\n"
            "sampleA\tpass\tpass\n"
            "sampleB\tpass\twarn\n",
            encoding="utf-8",
        )
        result = parse_multiqc_fastqc_tsv(data_dir)
        assert result["parse_ok"] is True
        assert "sampleA" in result["samples"]
        assert result["samples"]["sampleB"]["adapter_content"] == "warn"

    def test_parse_multiqc_fastqc_tsv_missing(self, tmp_dir):
        from genomics_workflow_agent.parsers.multiqc import parse_multiqc_fastqc_tsv
        result = parse_multiqc_fastqc_tsv(tmp_dir)
        assert result["parse_ok"] is False

    def test_parse_multiqc_output_with_html(self, tmp_dir):
        from genomics_workflow_agent.parsers.multiqc import parse_multiqc_output
        (tmp_dir / "multiqc_report.html").write_text("<html>MultiQC</html>")
        (tmp_dir / "multiqc_data").mkdir()
        result = parse_multiqc_output(tmp_dir)
        assert result["html_report"] is not None
        assert result["parse_ok"] is True


class TestDecisionEngine:
    def test_no_results_recommends_run_fastqc(self):
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results
        result = evaluate_fastqc_results([])
        obs_statuses = [o.status for o in result["observations"]]
        assert "missing" in obs_statuses
        actions = [a.action for a in result["recommended_actions"]]
        assert any("FastQC" in a for a in actions)

    def test_all_pass_recommends_no_trimming(self, all_pass_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results
        parsed = [parse_fastqc_txt(all_pass_txt, sample="sampleA_R1")]
        result = evaluate_fastqc_results(parsed)
        decision_types = [d.decision_type for d in result["decisions"]]
        assert "trim" not in decision_types
        assert "accept" in decision_types

    def test_adapter_warn_recommends_trimming(self, adapter_warn_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results
        parsed = [parse_fastqc_txt(adapter_warn_txt, sample="sampleB_R1")]
        result = evaluate_fastqc_results(parsed)
        decision_types = [d.decision_type for d in result["decisions"]]
        assert "trim" in decision_types

    def test_qual_fail_recommends_trim_and_review(self, qual_fail_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results
        parsed = [parse_fastqc_txt(qual_fail_txt, sample="sampleC_R1")]
        result = evaluate_fastqc_results(parsed)
        decision_types = [d.decision_type for d in result["decisions"]]
        assert "trim" in decision_types
        assert "review" in decision_types

    def test_no_clinical_or_biological_conclusions(self, qual_fail_txt, adapter_warn_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results

        parsed = [
            parse_fastqc_txt(qual_fail_txt, sample="sampleC_R1"),
            parse_fastqc_txt(adapter_warn_txt, sample="sampleB_R1"),
        ]
        result = evaluate_fastqc_results(parsed)

        # Must not make clinical claims
        FORBIDDEN = ["diagnosis", "pathogenic", "benign", "disease", "clinical finding"]
        for obs in result["observations"]:
            lower = obs.message.lower()
            for word in FORBIDDEN:
                assert word not in lower, (
                    f"Observation contains forbidden clinical term '{word}': {obs.message}"
                )
        for dec in result["decisions"]:
            lower = dec.reason.lower()
            for word in FORBIDDEN:
                assert word not in lower, (
                    f"Decision contains forbidden clinical term '{word}': {dec.reason}"
                )

    def test_gc_fail_does_not_recommend_auto_filter(self):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results

        txt = (FIXTURES / "all_pass_fastqc_data.txt").read_text()
        txt = txt.replace(">>Per sequence GC content\tpass", ">>Per sequence GC content\tfail")
        parsed = [parse_fastqc_txt(txt, sample="gc_fail_sample")]
        result = evaluate_fastqc_results(parsed)

        # GC fail must not trigger a trim decision
        decision_types = [d.decision_type for d in result["decisions"]]
        assert "trim" not in decision_types, "GC content failure must not trigger automatic trimming"
        # Should warn, not silently pass
        gc_obs = [o for o in result["observations"] if o.category == "gc_content"]
        assert len(gc_obs) > 0

    def test_trimming_decision_has_safety_notes(self, adapter_warn_txt):
        from genomics_workflow_agent.parsers.fastqc import parse_fastqc_txt
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results
        parsed = [parse_fastqc_txt(adapter_warn_txt, sample="s")]
        result = evaluate_fastqc_results(parsed)
        trim_decisions = [d for d in result["decisions"] if d.decision_type == "trim"]
        assert trim_decisions
        assert trim_decisions[0].safety_notes

    def test_parse_error_result_becomes_observation(self):
        from genomics_workflow_agent.agent.decision_engine import evaluate_fastqc_results
        bad_result = {
            "sample": "broken",
            "source": "broken_fastqc.zip",
            "modules": {},
            "summary": {"pass": [], "warn": [], "fail": []},
            "errors": ["Could not open zip"],
            "parse_ok": False,
        }
        result = evaluate_fastqc_results([bad_result])
        obs_categories = [o.category for o in result["observations"]]
        assert "parse_error" in obs_categories


class TestAgentState:
    def test_state_serializes_to_dict(self):
        from genomics_workflow_agent.agent.state import AgentState, Observation, Decision, RecommendedAction
        state = AgentState(input_path="/data", workflow="fastq-qc")
        state.observations.append(Observation(
            source="test", sample="s1", category="adapter_content",
            status="warn", severity="warning", message="Adapter warn",
        ))
        state.decisions.append(Decision(
            action="trim_reads", decision_type="trim", reason="Adapter warn",
        ))
        state.recommended_actions.append(RecommendedAction(
            action="Run fastp", priority="high", reason="Adapters found",
        ))
        d = state.to_dict()
        assert d["input_path"] == "/data"
        assert len(d["observations"]) == 1
        assert d["observations"][0]["category"] == "adapter_content"
        assert len(d["decisions"]) == 1
        assert len(d["recommended_actions"]) == 1


class TestFastqAgent:
    def test_dry_run_does_not_call_subprocess(self, fastq_dir, tmp_dir):
        with patch("subprocess.run") as mock_sub:
            from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent
            state = run_fastq_agent(fastq_dir, tmp_dir / "out", execute=False)
            mock_sub.assert_not_called()

        assert state.workflow == "fastq-qc"
        assert any("dry-run" in w.lower() or "dry run" in w.lower() for w in state.warnings)

    def test_dry_run_has_dry_run_observation(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent
        state = run_fastq_agent(fastq_dir, tmp_dir / "out", execute=False)
        obs_categories = [o.category for o in state.observations]
        assert "dry_run" in obs_categories

    def test_auto_trim_without_execute_raises(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent
        with pytest.raises(ValueError, match="--auto-trim requires --execute"):
            run_fastq_agent(fastq_dir, tmp_dir / "out", execute=False, auto_trim=True)

    def test_mocked_execute_parses_outputs(self, fastq_dir, tmp_dir, all_pass_txt):
        """Mock fastqc to create real-parseable output files, then check agent parses them."""
        out_dir = tmp_dir / "out"
        fastqc_dir = out_dir / "fastqc"
        multiqc_dir = out_dir / "multiqc"

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "", ""
            if "fastqc" in str(cmd[0]):
                fastqc_dir.mkdir(parents=True, exist_ok=True)
                for fname in ["sampleA_R1", "sampleA_R2"]:
                    _make_fastqc_zip(fastqc_dir, fname, all_pass_txt)
                    (fastqc_dir / f"{fname}_fastqc.html").write_text("<html>FastQC</html>")
            elif "multiqc" in str(cmd[0]):
                multiqc_dir.mkdir(parents=True, exist_ok=True)
                (multiqc_dir / "multiqc_report.html").write_text("<html>MultiQC</html>")
                (multiqc_dir / "multiqc_data").mkdir(exist_ok=True)
            return mock

        fake_tools = {t: {"available": t in ["fastqc", "multiqc"], "path": f"/usr/bin/{t}", "version": "0"}
                      for t in ["fastqc", "multiqc", "fastp", "cutadapt", "samtools",
                                "bcftools", "mosdepth", "nextflow", "docker", "singularity", "conda",
                                "qiime", "Rscript", "bedtools"]}

        with patch("subprocess.run", side_effect=fake_run), \
             patch("genomics_workflow_agent.workflows.fastq_qc.check_tools", return_value=fake_tools):
            from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent
            state = run_fastq_agent(fastq_dir, out_dir, execute=True)

        assert len(state.observations) > 0
        assert len(state.decisions) > 0

    def test_mocked_execute_adapter_warn_triggers_trim_decision(self, fastq_dir, tmp_dir, adapter_warn_txt):
        out_dir = tmp_dir / "out"
        fastqc_dir = out_dir / "fastqc"
        multiqc_dir = out_dir / "multiqc"

        def fake_run(cmd, **kwargs):
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "", ""
            if "fastqc" in str(cmd[0]):
                fastqc_dir.mkdir(parents=True, exist_ok=True)
                for fname in ["sampleA_R1", "sampleA_R2"]:
                    _make_fastqc_zip(fastqc_dir, fname, adapter_warn_txt)
                    (fastqc_dir / f"{fname}_fastqc.html").write_text("<html>FastQC</html>")
            elif "multiqc" in str(cmd[0]):
                multiqc_dir.mkdir(parents=True, exist_ok=True)
                (multiqc_dir / "multiqc_report.html").write_text("<html>MultiQC</html>")
                (multiqc_dir / "multiqc_data").mkdir()
            return mock

        fake_tools = {t: {"available": t in ["fastqc", "multiqc"], "path": f"/usr/bin/{t}", "version": "0"}
                      for t in ["fastqc", "multiqc", "fastp", "cutadapt", "samtools",
                                "bcftools", "mosdepth", "nextflow", "docker", "singularity", "conda",
                                "qiime", "Rscript", "bedtools"]}

        with patch("subprocess.run", side_effect=fake_run), \
             patch("genomics_workflow_agent.workflows.fastq_qc.check_tools", return_value=fake_tools):
            from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent
            state = run_fastq_agent(fastq_dir, out_dir, execute=True, auto_trim=False)

        trim_decisions = [d for d in state.decisions if d.decision_type == "trim"]
        assert len(trim_decisions) > 0
        assert not trim_decisions[0].executed

    def test_auto_trim_requires_trim_decision_to_execute(self, fastq_dir, tmp_dir, all_pass_txt):
        """auto-trim must NOT run if decision engine says no trimming needed."""
        out_dir = tmp_dir / "out"
        fastqc_dir = out_dir / "fastqc"
        multiqc_dir = out_dir / "multiqc"

        trim_was_called = []

        def fake_run(cmd, **kwargs):
            if any(t in str(cmd[0]) for t in ["fastp", "cutadapt"]):
                trim_was_called.append(cmd)
            mock = MagicMock()
            mock.returncode = 0
            mock.stdout, mock.stderr = "", ""
            if "fastqc" in str(cmd[0]):
                fastqc_dir.mkdir(parents=True, exist_ok=True)
                for fname in ["sampleA_R1", "sampleA_R2"]:
                    _make_fastqc_zip(fastqc_dir, fname, all_pass_txt)
                    (fastqc_dir / f"{fname}_fastqc.html").write_text("<html>FastQC</html>")
            elif "multiqc" in str(cmd[0]):
                multiqc_dir.mkdir(parents=True, exist_ok=True)
                (multiqc_dir / "multiqc_report.html").write_text("<html>MultiQC</html>")
                (multiqc_dir / "multiqc_data").mkdir()
            return mock

        fake_tools = {t: {"available": t in ["fastqc", "multiqc", "fastp"], "path": f"/usr/bin/{t}", "version": "0"}
                      for t in ["fastqc", "multiqc", "fastp", "cutadapt", "samtools",
                                "bcftools", "mosdepth", "nextflow", "docker", "singularity", "conda",
                                "qiime", "Rscript", "bedtools"]}

        with patch("subprocess.run", side_effect=fake_run), \
             patch("genomics_workflow_agent.workflows.fastq_qc.check_tools", return_value=fake_tools):
            from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent
            state = run_fastq_agent(fastq_dir, out_dir, execute=True, auto_trim=True)

        assert not trim_was_called, "Trimmer was called despite no trimming decision"

    def test_agent_report_json_written(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent, write_agent_report_json
        state = run_fastq_agent(fastq_dir, tmp_dir / "out", execute=False)
        json_path = write_agent_report_json(state, tmp_dir / "out" / "agent_report.json")
        assert json_path.exists()
        data = json.loads(json_path.read_text())
        assert "observations" in data
        assert "decisions" in data
        assert "clinical_disclaimer" in data

    def test_agent_report_md_written(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.agent.fastq_agent import run_fastq_agent, write_agent_report_md
        state = run_fastq_agent(fastq_dir, tmp_dir / "out", execute=False)
        md_path = write_agent_report_md(state, tmp_dir / "out" / "agent_report.md")
        assert md_path.exists()
        content = md_path.read_text()
        assert "FASTQ QC" in content
        assert "Disclaimer" in content


class TestCLIAgentCommand:
    def test_cli_agent_dry_run(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "agent",
            "--input", str(fastq_dir),
            "--workflow", "fastq-qc",
            "--out", str(tmp_dir / "out"),
        ]
        result = main()
        assert result == 0
        assert (tmp_dir / "out" / "agent_report.json").exists()
        assert (tmp_dir / "out" / "agent_report.md").exists()

    def test_cli_agent_auto_trim_without_execute_exits_1(self, fastq_dir, tmp_dir):
        from genomics_workflow_agent.cli import main
        sys.argv = [
            "genomics_workflow_agent", "agent",
            "--input", str(fastq_dir),
            "--out", str(tmp_dir / "out"),
            "--auto-trim",
        ]
        result = main()
        assert result == 1

    def test_cli_agent_help(self):
        from genomics_workflow_agent.cli import build_parser
        parser = build_parser()
        # Parsing --help would sys.exit(0); just confirm agent is registered
        sub_parsers_action = next(
            a for a in parser._actions if hasattr(a, "_name_parser_map")
        )
        assert "agent" in sub_parsers_action._name_parser_map

    def test_cli_agent_dry_run_no_subprocess(self, fastq_dir, tmp_dir):
        with patch("subprocess.run") as mock_sub:
            from genomics_workflow_agent.cli import main
            sys.argv = [
                "genomics_workflow_agent", "agent",
                "--input", str(fastq_dir),
                "--out", str(tmp_dir / "out"),
            ]
            result = main()
            mock_sub.assert_not_called()
        assert result == 0
