"""Tests for tools/check_environment.py."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
TOOL = REPO_ROOT / "tools" / "check_environment.py"
DEMO_SCRIPT = REPO_ROOT / "examples" / "run_tiny_demo.sh"


def run_tool(*args):
    return subprocess.run(
        [sys.executable, str(TOOL)] + list(args),
        capture_output=True,
        text=True,
    )


class TestHelp:
    def test_help_exits_zero(self):
        r = run_tool("--help")
        assert r.returncode == 0

    def test_help_mentions_output_dir(self):
        r = run_tool("--help")
        assert "output-dir" in r.stdout or "output_dir" in r.stdout

    def test_help_mentions_verbose(self):
        r = run_tool("--help")
        assert "verbose" in r.stdout


class TestJsonOutput:
    def test_json_file_created(self, tmp_path):
        r = run_tool("--output-dir", str(tmp_path))
        assert r.returncode in (0, 1), f"Unexpected exit code {r.returncode}\n{r.stderr}"
        assert (tmp_path / "environment_check.json").exists(), "JSON file not created"

    def test_json_is_valid(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        assert isinstance(data, dict)

    def test_json_has_required_top_level_keys(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        for key in ("tool", "version", "generated", "python", "required_packages",
                    "optional_packages", "external_tools", "output_directory",
                    "official_skills", "summary"):
            assert key in data, f"Missing top-level key: {key}"

    def test_json_summary_has_ready_field(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        assert "ready" in data["summary"], "summary.ready field missing"
        assert isinstance(data["summary"]["ready"], bool)

    def test_json_python_version_present(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        assert "version" in data["python"]
        assert data["python"]["version"].startswith(
            f"{sys.version_info.major}.{sys.version_info.minor}"
        )

    def test_json_required_packages_list(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        pkgs = data["required_packages"]
        assert isinstance(pkgs, list)
        assert len(pkgs) > 0
        for p in pkgs:
            assert "name" in p
            assert "status" in p
            assert p["status"] in ("ok", "MISSING", "warning")

    def test_json_optional_packages_are_warnings_not_failures(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        for p in data["optional_packages"]:
            if p["status"] != "ok":
                assert p["status"] == "warning", (
                    f"Optional package {p['name']} has status {p['status']} - "
                    "missing optional packages must be 'warning', not 'MISSING' or 'FAIL'"
                )

    def test_json_external_tools_are_warnings_not_failures(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        for t in data["external_tools"]:
            if t["status"] != "ok":
                assert t["status"] == "warning", (
                    f"CLI tool {t['name']} has status {t['status']} - "
                    "missing CLI tools must be 'warning', not 'MISSING' or 'FAIL'"
                )

    def test_json_official_skills_has_manual_status(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        skills = data["official_skills"]
        assert skills["status"] == "manual", "Official skills status should be 'manual' (cannot verify from script)"
        assert "single-cell-rna-qc@life-sciences" in str(skills)

    def test_missing_optional_tools_do_not_cause_failure(self, tmp_path):
        """Exit code should be 0 when all required packages are present, even if optional tools are missing."""
        run_tool("--output-dir", str(tmp_path))
        with open(tmp_path / "environment_check.json", encoding="utf-8") as f:
            data = json.load(f)
        n_required_missing = sum(
            1 for p in data["required_packages"] if p["status"] not in ("ok",)
        )
        if n_required_missing == 0:
            # If required packages are all OK, the tool should exit 0
            r = run_tool("--output-dir", str(tmp_path))
            assert r.returncode == 0, (
                "Tool should exit 0 when all required packages are present "
                f"(exit code was {r.returncode})"
            )


class TestMarkdownOutput:
    def test_markdown_file_created(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        assert (tmp_path / "environment_check.md").exists(), "Markdown file not created"

    def test_markdown_has_headings(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        content = (tmp_path / "environment_check.md").read_text(encoding="utf-8")
        assert "# Environment Check Report" in content
        assert "## Required Python Packages" in content
        assert "## Optional Python Packages" in content
        assert "## External CLI Tools" in content
        assert "## Summary" in content

    def test_markdown_notes_optional_tools_are_warnings(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        content = (tmp_path / "environment_check.md").read_text(encoding="utf-8")
        assert "optional" in content.lower() or "warning" in content.lower()

    def test_markdown_mentions_official_skills(self, tmp_path):
        run_tool("--output-dir", str(tmp_path))
        content = (tmp_path / "environment_check.md").read_text(encoding="utf-8")
        assert "single-cell-rna-qc@life-sciences" in content


class TestDemoScript:
    def test_demo_script_exists(self):
        assert DEMO_SCRIPT.exists(), f"run_tiny_demo.sh not found at {DEMO_SCRIPT}"

    def test_demo_script_is_executable_text(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert content.startswith("#!/usr/bin/env bash"), "Script missing bash shebang"

    def test_demo_script_references_check_environment(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "check_environment.py" in content

    def test_demo_script_references_inspect_file(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "inspect_file.py" in content

    def test_demo_script_references_scrna_qc(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "scrna_qc_local.py" in content

    def test_demo_script_references_atac_qc(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "atac_qc_local.py" in content

    def test_demo_script_references_wgs_qc(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "wgs_vcf_qc_local.py" in content

    def test_demo_script_references_report_builder(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "report_builder.py" in content

    def test_demo_script_targets_final_report(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "final_report.md" in content

    def test_demo_script_mentions_official_skill_preference(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "single-cell-rna-qc@life-sciences" in content, (
            "Demo script must note that local scRNA QC is a fallback "
            "and single-cell-rna-qc@life-sciences is preferred"
        )

    def test_demo_script_uses_set_e(self):
        content = DEMO_SCRIPT.read_text(encoding="utf-8")
        assert "set -e" in content, "Script should use set -e for early exit on failure"
