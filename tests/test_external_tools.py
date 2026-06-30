"""
test_external_tools.py - Parser-level and integration tests for external tool support.

Tests do NOT require samtools, bcftools, or bedtools to be installed.
They test:
  - Pure parser functions directly (parse_flagstat_output, parse_idxstats_output,
    parse_bcftools_stats)
  - Graceful degradation when external tools are absent
  - Commands recorded in JSON when tools are used
  - Biological caveats present in both JSON and Markdown outputs
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
FIXTURES = REPO_ROOT / "tests" / "fixtures"
WGS_TOOL = REPO_ROOT / "tools" / "wgs_vcf_qc_local.py"
ATAC_TOOL = REPO_ROOT / "tools" / "atac_qc_local.py"
EXAMPLES = REPO_ROOT / "examples"

# Import pure parser functions directly (no subprocess, no external tools)
sys.path.insert(0, str(REPO_ROOT / "tools"))
from wgs_vcf_qc_local import (
    parse_flagstat_output,
    parse_idxstats_output,
    parse_samtools_stats_summary,
    parse_bcftools_stats,
    compare_vcf_stats,
)


def run_wgs(*args):
    return subprocess.run(
        [sys.executable, str(WGS_TOOL)] + list(args),
        capture_output=True, text=True,
    )


def run_atac(*args):
    return subprocess.run(
        [sys.executable, str(ATAC_TOOL)] + list(args),
        capture_output=True, text=True,
    )


# ──────────────────────────────────────────────────────────
# Parse fixture: samtools flagstat
# ──────────────────────────────────────────────────────────

class TestParseFlagstat:
    def test_total_reads(self):
        text = (FIXTURES / "samtools_flagstat.txt").read_text()
        parsed = parse_flagstat_output(text)
        assert parsed["in_total"]["qc_passed"] == 500000

    def test_mapped_reads(self):
        text = (FIXTURES / "samtools_flagstat.txt").read_text()
        parsed = parse_flagstat_output(text)
        assert parsed["mapped"]["qc_passed"] == 490000

    def test_duplicate_count(self):
        text = (FIXTURES / "samtools_flagstat.txt").read_text()
        parsed = parse_flagstat_output(text)
        assert parsed["duplicates"]["qc_passed"] == 12500

    def test_pct_mapped_derived(self):
        text = (FIXTURES / "samtools_flagstat.txt").read_text()
        parsed = parse_flagstat_output(text)
        assert "_pct_mapped" in parsed
        assert parsed["_pct_mapped"] == pytest.approx(98.0, abs=0.1)

    def test_pct_duplicate_derived(self):
        text = (FIXTURES / "samtools_flagstat.txt").read_text()
        parsed = parse_flagstat_output(text)
        assert "_pct_duplicate" in parsed
        assert parsed["_pct_duplicate"] == pytest.approx(2.5, abs=0.1)

    def test_empty_string_returns_dict(self):
        parsed = parse_flagstat_output("")
        assert isinstance(parsed, dict)

    def test_malformed_lines_skipped(self):
        text = "not a valid line\n500000 + 0 in total (QC-passed reads + QC-failed reads)\n"
        parsed = parse_flagstat_output(text)
        assert parsed["in_total"]["qc_passed"] == 500000


# ──────────────────────────────────────────────────────────
# Parse fixture: samtools idxstats
# ──────────────────────────────────────────────────────────

class TestParseIdxstats:
    def test_record_count(self):
        text = (FIXTURES / "samtools_idxstats.txt").read_text()
        records = parse_idxstats_output(text)
        assert len(records) == 6  # 5 chroms + *

    def test_first_chrom(self):
        text = (FIXTURES / "samtools_idxstats.txt").read_text()
        records = parse_idxstats_output(text)
        assert records[0]["chrom"] == "chr1"
        assert records[0]["mapped"] == 150000
        assert records[0]["length"] == 248956422

    def test_unmapped_star_row(self):
        text = (FIXTURES / "samtools_idxstats.txt").read_text()
        records = parse_idxstats_output(text)
        star = [r for r in records if r["chrom"] == "*"]
        assert len(star) == 1
        assert star[0]["unmapped"] == 3000

    def test_empty_string_returns_list(self):
        records = parse_idxstats_output("")
        assert isinstance(records, list)
        assert len(records) == 0


# ──────────────────────────────────────────────────────────
# Parse fixture: bcftools stats
# ──────────────────────────────────────────────────────────

class TestParseBcftoolsStats:
    def test_snp_count(self):
        text = (FIXTURES / "bcftools_stats.txt").read_text()
        parsed = parse_bcftools_stats(text)
        assert parsed.get("number of SNPs") == 6

    def test_indel_count(self):
        text = (FIXTURES / "bcftools_stats.txt").read_text()
        parsed = parse_bcftools_stats(text)
        assert parsed.get("number of indels") == 1

    def test_titv_present(self):
        text = (FIXTURES / "bcftools_stats.txt").read_text()
        parsed = parse_bcftools_stats(text)
        assert "ts/tv" in parsed
        assert float(parsed["ts/tv"]) == pytest.approx(2.0, abs=0.01)

    def test_sample_count(self):
        text = (FIXTURES / "bcftools_stats.txt").read_text()
        parsed = parse_bcftools_stats(text)
        assert parsed.get("number of samples") == 1

    def test_empty_string_returns_dict(self):
        parsed = parse_bcftools_stats("")
        assert isinstance(parsed, dict)

    def test_comment_lines_ignored(self):
        text = "# comment\n# another comment\nSN\t0\tnumber of SNPs:\t42\n"
        parsed = parse_bcftools_stats(text)
        assert parsed.get("number of SNPs") == 42


# ──────────────────────────────────────────────────────────
# compare_vcf_stats - discrepancy detection
# ──────────────────────────────────────────────────────────

class TestCompareVcfStats:
    def test_no_discrepancy_when_equal(self):
        local = {"n_snps": 6, "n_indels": 1, "ti_tv_ratio": 2.0}
        bc = {"number of SNPs": 6, "number of indels": 1, "ts/tv": 2.0}
        warnings = compare_vcf_stats(local, bc)
        assert warnings == []

    def test_snp_discrepancy_detected(self):
        local = {"n_snps": 5, "n_indels": 1, "ti_tv_ratio": 2.0}
        bc = {"number of SNPs": 6, "number of indels": 1, "ts/tv": 2.0}
        warnings = compare_vcf_stats(local, bc)
        assert any(w["metric"] == "n_snps" for w in warnings)

    def test_titv_discrepancy_detected(self):
        local = {"n_snps": 6, "n_indels": 1, "ti_tv_ratio": 1.5}
        bc = {"number of SNPs": 6, "number of indels": 1, "ts/tv": 2.0}
        warnings = compare_vcf_stats(local, bc)
        assert any(w["metric"] == "ti_tv_ratio" for w in warnings)

    def test_discrepancy_includes_both_values(self):
        local = {"n_snps": 5, "n_indels": 1, "ti_tv_ratio": 2.0}
        bc = {"number of SNPs": 6, "number of indels": 1, "ts/tv": 2.0}
        warnings = compare_vcf_stats(local, bc)
        snp_w = next(w for w in warnings if w["metric"] == "n_snps")
        assert snp_w["local_value"] == 5
        assert snp_w["bcftools_value"] == 6


# ──────────────────────────────────────────────────────────
# Graceful degradation: samtools missing
# ──────────────────────────────────────────────────────────

class TestSamtoolsMissingDegradation:
    """When no BAM is provided, BAM/CRAM QC must be skipped with biological caveat."""

    def test_runs_without_bam(self, tmp_path):
        r = run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        assert r.returncode == 0

    def test_bam_skipped_has_biological_caveat(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        bam_skipped = [
            s for s in data["skipped_metrics"]
            if "BAM" in s["metric"] or "bam" in s["metric"].lower()
        ]
        assert len(bam_skipped) > 0
        for s in bam_skipped:
            assert "missing_biological_conclusion" in s
            assert len(s["missing_biological_conclusion"]) > 20

    def test_samtools_availability_recorded_in_json(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        assert "samtools_available" in data


# ──────────────────────────────────────────────────────────
# Graceful degradation: bcftools missing
# ──────────────────────────────────────────────────────────

class TestBcftoolsMissingDegradation:
    """Local VCF parser must work and report correctly when bcftools is absent."""

    def test_vcf_parsed_successfully(self, tmp_path):
        r = run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        assert r.returncode == 0
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        assert data["vcf_metrics"]["n_snps"] >= 5

    def test_bcftools_availability_recorded(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        assert "bcftools_available" in data

    def test_bcftools_metrics_field_present(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        # Field must exist (either results or skipped-with-caveat)
        assert "bcftools_metrics" in data

    def test_bcftools_skipped_has_biological_caveat(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        bc = data["bcftools_metrics"]
        # If bcftools not installed, must have skipped=True with caveat
        if bc.get("skipped"):
            assert "missing_biological_conclusion" in bc
            assert len(bc["missing_biological_conclusion"]) > 20


# ──────────────────────────────────────────────────────────
# Graceful degradation: bedtools missing (ATAC)
# ──────────────────────────────────────────────────────────

class TestBedtoolsMissingDegradation:
    """ATAC QC must run without bedtools, local FRiP fallback intact."""

    def test_atac_runs_without_bedtools(self, tmp_path):
        r = run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--peaks", str(EXAMPLES / "tiny_peaks.bed"),
            "--output-dir", str(tmp_path),
        )
        assert r.returncode == 0

    def test_local_frip_still_computed(self, tmp_path):
        run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--peaks", str(EXAMPLES / "tiny_peaks.bed"),
            "--output-dir", str(tmp_path),
        )
        with open(tmp_path / "atac_qc_summary.json") as f:
            data = json.load(f)
        frip = data["metrics"]["frip"].get("frip")
        assert frip is not None
        assert 0 <= frip <= 1

    def test_bedtools_availability_recorded(self, tmp_path):
        run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--output-dir", str(tmp_path),
        )
        with open(tmp_path / "atac_qc_summary.json") as f:
            data = json.load(f)
        assert "bedtools_available" in data

    def test_blacklist_skipped_has_biological_caveat(self, tmp_path):
        r = run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--blacklist", str(EXAMPLES / "tiny_peaks.bed"),  # reuse as mock blacklist
            "--output-dir", str(tmp_path),
        )
        assert r.returncode == 0
        with open(tmp_path / "atac_qc_summary.json") as f:
            data = json.load(f)
        # Either bedtools ran (blacklist_overlap has results) or was skipped with caveat
        skipped = data.get("skipped_metrics", [])
        bl_skipped = [s for s in skipped if "blacklist" in s["metric"].lower() or "Blacklist" in s["metric"]]
        bl_overlap = data.get("blacklist_overlap")
        # At least one of: computed result OR skipped metric with caveat
        if bl_skipped:
            for s in bl_skipped:
                assert "missing_biological_conclusion" in s
        else:
            # bedtools ran and produced a result
            assert bl_overlap is not None


# ──────────────────────────────────────────────────────────
# Commands recorded when external tools are used
# ──────────────────────────────────────────────────────────

class TestCommandsRecorded:
    def test_tool_metadata_in_json(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        assert data["tool"] == "wgs_vcf_qc_local.py"
        assert "version" in data
        assert "generated" in data

    def test_samtools_qc_field_present(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        # samtools_qc key always present (may be empty dict if no BAM)
        assert "samtools_qc" in data

    def test_atac_tool_metadata_in_json(self, tmp_path):
        run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--output-dir", str(tmp_path),
        )
        with open(tmp_path / "atac_qc_summary.json") as f:
            data = json.load(f)
        assert data["tool"] == "atac_qc_local.py"
        assert "version" in data


# ──────────────────────────────────────────────────────────
# Biological caveats in JSON and Markdown
# ──────────────────────────────────────────────────────────

class TestBiologicalCaveats:
    def test_all_wgs_skipped_metrics_have_biological_conclusion(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        for s in data["skipped_metrics"]:
            assert "missing_biological_conclusion" in s, (
                f"Skipped metric '{s['metric']}' lacks missing_biological_conclusion"
            )

    def test_wgs_biological_caveats_in_markdown(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        md = (tmp_path / "wgs_vcf_qc_report.md").read_text(encoding="utf-8")
        assert "Missing biological conclusion" in md

    def test_all_atac_skipped_metrics_have_biological_conclusion(self, tmp_path):
        run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--output-dir", str(tmp_path),
        )
        with open(tmp_path / "atac_qc_summary.json") as f:
            data = json.load(f)
        for s in data["skipped_metrics"]:
            assert "missing_biological_conclusion" in s, (
                f"Skipped metric '{s['metric']}' lacks missing_biological_conclusion"
            )

    def test_atac_biological_caveats_in_markdown(self, tmp_path):
        run_atac(
            "--fragments", str(EXAMPLES / "tiny_fragments.tsv"),
            "--output-dir", str(tmp_path),
        )
        md = (tmp_path / "atac_qc_report.md").read_text(encoding="utf-8")
        assert "Missing biological conclusion" in md or "biological" in md.lower()

    def test_wgs_clinical_disclaimer_in_json(self, tmp_path):
        run_wgs("--vcf", str(EXAMPLES / "tiny.vcf"), "--output-dir", str(tmp_path))
        with open(tmp_path / "wgs_vcf_qc_summary.json") as f:
            data = json.load(f)
        assert "clinical_disclaimer" in data
        assert "No clinical claims" in data["clinical_disclaimer"]
