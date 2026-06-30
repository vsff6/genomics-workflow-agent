"""
Variant/WGS QC workflow — samtools, bcftools, mosdepth, and nf-core/sarek planning.

Execution model:
  plan()    → always dry-run; returns steps with expected_outputs
  execute() → runs samtools/bcftools steps; sarek via Nextflow (dry-run by default unless --execute)

External tools required:
  samtools (flagstat, idxstats, stats)
  bcftools (stats, index)
  mosdepth (coverage, optional)
  nextflow (for sarek, optional)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomics_workflow_agent.safety.guardrails import DISCLAIMER
from genomics_workflow_agent.tools.nextflow import build_nextflow_cmd, run_nextflow
from genomics_workflow_agent.tools.runner import (
    STATUS_SKIPPED,
    run_command,
)
from genomics_workflow_agent.tools.samplesheets import build_sarek_samplesheet
from genomics_workflow_agent.tools.versions import check_tools

OPTIONAL_TOOLS = ["samtools", "bcftools", "nextflow", "docker", "singularity", "conda", "mosdepth"]

BIOLOGICAL_CAVEATS = [
    "No clinical interpretation is provided. Variant calls must not be used for medical decisions.",
    "Tumor/normal pairing, pedigree structure, and germline/somatic context must be confirmed by the analyst.",
    "Reference build and known-sites VCF compatibility must be verified for BQSR.",
    "Ti/Tv ratio expected values vary by capture kit, read length, variant caller, and genome region.",
    "Het/hom ratio is affected by ancestry, ploidy, copy number variation, and filtering stringency.",
    "Coverage uniformity depends on library preparation, target capture, and GC content — not biology alone.",
    "Pathogenicity requires expert review with validated annotation sources (ClinVar, ACMG, expert panels).",
]

_BAM_EXTS = {".bam", ".cram"}
_VCF_EXTS = {".vcf", ".vcf.gz", ".bcf"}
_FASTQ_EXTS = {".fastq.gz", ".fq.gz", ".fastq", ".fq"}


def _find_files(directory: Path, exts: set) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in exts
                  or any(str(p.name).endswith(e) for e in exts))


def plan(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    genome: str | None = None,
    fasta: str | None = None,
    known_sites: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    dry_run: bool = True,
) -> dict[str, Any]:
    """Return a WGS/variant QC execution plan."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    qc_outdir = output_dir / "variant_qc"

    tools = check_tools(OPTIONAL_TOOLS)
    warnings: list[str] = []
    blockers: list[str] = []
    skipped: list[dict] = []
    steps: list[dict] = []

    bam_files = _find_files(input_dir, _BAM_EXTS)
    vcf_files = _find_files(input_dir, _VCF_EXTS)
    fastq_files = _find_files(input_dir, _FASTQ_EXTS)

    if bam_files:
        if tools["samtools"]["available"]:
            for bam in bam_files[:10]:
                stem = bam.stem
                steps += [
                    {
                        "name": f"samtools_flagstat_{stem}",
                        "description": f"Alignment statistics: {bam.name}",
                        "command": ["samtools", "flagstat", str(bam)],
                        "output_dir": str(qc_outdir),
                        "capture_stdout_path": str(qc_outdir / f"{stem}_flagstat.txt"),
                        "expected_outputs": [str(qc_outdir / f"{stem}_flagstat.txt")],
                        "dry_run": dry_run,
                        "required_tools": ["samtools"],
                    },
                    {
                        "name": f"samtools_idxstats_{stem}",
                        "description": f"Index statistics: {bam.name}",
                        "command": ["samtools", "idxstats", str(bam)],
                        "output_dir": str(qc_outdir),
                        "capture_stdout_path": str(qc_outdir / f"{stem}_idxstats.txt"),
                        "expected_outputs": [str(qc_outdir / f"{stem}_idxstats.txt")],
                        "dry_run": dry_run,
                        "required_tools": ["samtools"],
                    },
                    {
                        "name": f"samtools_stats_{stem}",
                        "description": f"Detailed alignment stats: {bam.name}",
                        "command": ["samtools", "stats", str(bam)],
                        "output_dir": str(qc_outdir),
                        "capture_stdout_path": str(qc_outdir / f"{stem}_stats.txt"),
                        "expected_outputs": [str(qc_outdir / f"{stem}_stats.txt")],
                        "dry_run": dry_run,
                        "required_tools": ["samtools"],
                    },
                ]
            if len(bam_files) > 10:
                warnings.append(f"Only planning samtools commands for first 10 of {len(bam_files)} BAM files")
        else:
            skipped.append({
                "step": "samtools_flagstat/idxstats/stats",
                "reason": "samtools not found",
                "status": STATUS_SKIPPED,
                "install": "conda install -c bioconda samtools",
            })

        # mosdepth
        if tools["mosdepth"]["available"]:
            for bam in bam_files[:10]:
                stem = bam.stem
                mosdepth_prefix = str(qc_outdir / f"{stem}_mosdepth")
                steps.append({
                    "name": f"mosdepth_{stem}",
                    "description": f"Coverage statistics: {bam.name}",
                    "command": ["mosdepth", "--threads", "4", mosdepth_prefix, str(bam)],
                    "output_dir": str(qc_outdir),
                    "expected_outputs": [
                        mosdepth_prefix + ".mosdepth.summary.txt",
                        mosdepth_prefix + ".mosdepth.global.dist.txt",
                    ],
                    "dry_run": dry_run,
                    "required_tools": ["mosdepth"],
                })
        else:
            warnings.append("mosdepth not found — coverage statistics skipped. Install: conda install -c bioconda mosdepth")

    if vcf_files:
        if tools["bcftools"]["available"]:
            for vcf in vcf_files[:10]:
                stem = vcf.name.replace(".vcf.gz", "").replace(".vcf", "").replace(".bcf", "")
                steps.append({
                    "name": f"bcftools_stats_{stem}",
                    "description": f"VCF statistics: {vcf.name}",
                    "command": ["bcftools", "stats", str(vcf)],
                    "output_dir": str(qc_outdir),
                    "capture_stdout_path": str(qc_outdir / f"{stem}_bcftools_stats.txt"),
                    "expected_outputs": [str(qc_outdir / f"{stem}_bcftools_stats.txt")],
                    "dry_run": dry_run,
                    "required_tools": ["bcftools"],
                })
        else:
            skipped.append({
                "step": "bcftools_stats",
                "reason": "bcftools not found",
                "status": STATUS_SKIPPED,
                "install": "conda install -c bioconda bcftools",
            })

    sarek_plan: dict = {"created": False}
    sarek_cmd: list[str] = []

    if fastq_files:
        ss_path = output_dir / "samplesheet_sarek.csv"
        samplesheet_result = build_sarek_samplesheet(input_dir, ss_path)
        pipeline_outdir = str(output_dir / "sarek_output")

        if not genome and not fasta:
            warnings.append("No reference genome specified for sarek — command will be incomplete")

        sarek_cmd = build_nextflow_cmd(
            "sarek",
            input_path=str(ss_path) if samplesheet_result.get("created") else None,
            outdir=pipeline_outdir,
            genome=genome,
            fasta=fasta,
            known_snps=known_sites,
            profile=profile,
            resume=resume,
            extra_args=extra_args,
        )

        sarek_plan = {
            "samplesheet": samplesheet_result,
            "command": sarek_cmd,
            "command_str": " ".join(sarek_cmd),
            "pipeline_outdir": pipeline_outdir,
            "note": "Sarek samplesheet requires manual review — PATIENT_ID, sex, status must be set manually",
        }

        if not tools["nextflow"]["available"]:
            blockers.append("nextflow not found — sarek cannot be executed")

    if not bam_files and not vcf_files and not fastq_files:
        warnings.append("No BAM, VCF, or FASTQ files detected in input directory")

    return {
        "workflow": "variant-qc",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "qc_outdir": str(qc_outdir),
        "dry_run": dry_run,
        "bam_files_found": [str(f) for f in bam_files],
        "vcf_files_found": [str(f) for f in vcf_files],
        "fastq_files_found": len(fastq_files),
        "tool_status": tools,
        "steps": steps,
        "skipped_steps": skipped,
        "sarek_plan": sarek_plan,
        "blockers": blockers,
        "warnings": warnings,
        "biological_caveats": BIOLOGICAL_CAVEATS,
        "clinical_disclaimer": DISCLAIMER,
        "next_actions": [
            "Review sarek samplesheet — PATIENT_ID, sex, and status must be set manually",
            "Confirm genome build and known-sites VCF for BQSR",
            "Ensure BAM files are sorted and indexed before running samtools",
            "After QC: review MultiQC output and invoke biology-interpretation-reviewer",
            "Do not interpret variants clinically without validated annotation and expert review",
        ],
        "limitations": [
            "No clinical claims are made from any output.",
            "Sarek samplesheet requires manual review before execution.",
            "Coverage metrics from mosdepth are not parsed — review text outputs directly.",
        ],
    }


def execute(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    genome: str | None = None,
    fasta: str | None = None,
    known_sites: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    provenance_dir: Path | None = None,
    timeout: int = 3600,
) -> dict[str, Any]:
    """
    Execute variant QC steps (samtools/bcftools/mosdepth).

    Sarek is planned but NOT executed here — use nextflow execute for that.
    Returns annotated step records with status and output validation.
    """
    output_dir = Path(output_dir)
    qc_outdir = output_dir / "variant_qc"

    wf_plan = plan(
        input_dir, output_dir,
        genome=genome, fasta=fasta, known_sites=known_sites,
        profile=profile, resume=resume, extra_args=extra_args, dry_run=False,
    )

    # Create output directory
    qc_outdir.mkdir(parents=True, exist_ok=True)

    step_results = []
    for step in wf_plan["steps"]:
        cmd = step.get("command")
        if not cmd:
            continue

        record = run_command(
            cmd,
            dry_run=False,
            timeout=timeout,
            provenance_dir=provenance_dir,
            label=step["name"],
            expected_outputs=step.get("expected_outputs"),
            capture_stdout_path=step.get("capture_stdout_path"),
        )
        step["execution"] = record
        step_results.append(record)

    wf_plan["step_results"] = step_results
    wf_plan["dry_run"] = False
    return wf_plan
