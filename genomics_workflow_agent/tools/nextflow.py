from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from genomics_workflow_agent.tools.runner import (
    STATUS_FAILED,
    STATUS_PLANNED,
    STATUS_SKIPPED,
    STATUS_SUCCEEDED,
    run_command,
    validate_outputs,
)

SUPPORTED_PIPELINES = ["rnaseq", "sarek", "atacseq", "ampliseq"]


def is_nextflow_available() -> bool:
    return shutil.which("nextflow") is not None


def build_nextflow_cmd(
    pipeline: str,
    *,
    input_path: str | None = None,
    outdir: str | None = None,
    genome: str | None = None,
    fasta: str | None = None,
    gtf: str | None = None,
    blacklist: str | None = None,
    known_snps: str | None = None,
    primer_fw: str | None = None,
    primer_rv: str | None = None,
    taxonomy_param: str | None = None,
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    cmd = ["nextflow", "run", f"nf-core/{pipeline}"]

    if input_path:
        cmd += ["--input", input_path]
    if outdir:
        cmd += ["--outdir", outdir]
    if genome:
        cmd += ["--genome", genome]
    if fasta:
        cmd += ["--fasta", fasta]
    if gtf:
        cmd += ["--gtf", gtf]
    if blacklist:
        cmd += ["--blacklist", blacklist]
    if known_snps:
        cmd += ["--known_snps", known_snps]
    if primer_fw:
        cmd += ["--FW_primer", primer_fw]
    if primer_rv:
        cmd += ["--RV_primer", primer_rv]
    if taxonomy_param:
        cmd += ["--dada_ref_taxonomy", taxonomy_param]

    cmd += ["-profile", profile]
    if resume:
        cmd += ["-resume"]
    if extra_args:
        cmd += extra_args

    return cmd


def run_nextflow(
    pipeline: str,
    cmd: list[str],
    *,
    outdir: str | Path,
    dry_run: bool = True,
    provenance_dir: Path | None = None,
    timeout: int = 86400,
) -> dict[str, Any]:
    outdir = Path(outdir)
    label = f"nextflow_{pipeline}"

    if not is_nextflow_available() and not dry_run:
        return {
            "label": label,
            "command": cmd,
            "command_str": " ".join(cmd),
            "dry_run": dry_run,
            "executed": False,
            "status": STATUS_SKIPPED,
            "error": "nextflow not found in PATH. Install: https://www.nextflow.io",
            "return_code": None,
        }

    record = run_command(
        cmd,
        dry_run=dry_run,
        timeout=timeout,
        provenance_dir=provenance_dir,
        label=label,
        expected_outputs=[str(outdir)] if not dry_run else None,
    )

    if record.get("executed") and record.get("return_code") == 0:
        pipeline_validation = _validate_pipeline_outputs(pipeline, outdir)
        record["pipeline_output_validation"] = pipeline_validation
        if not pipeline_validation["has_multiqc"] and not pipeline_validation["has_pipeline_info"]:
            record["status"] = STATUS_FAILED
            record["error"] = (
                record.get("error") or
                f"Pipeline outdir exists but expected outputs not found in {outdir}"
            )

    return record


def _validate_pipeline_outputs(pipeline: str, outdir: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "outdir": str(outdir),
        "outdir_exists": outdir.exists(),
        "has_multiqc": False,
        "has_pipeline_info": False,
        "has_multiqc_html": False,
        "pipeline_specific": {},
    }

    if not outdir.exists():
        return result

    result["has_multiqc"] = (outdir / "multiqc").exists()
    result["has_pipeline_info"] = (outdir / "pipeline_info").exists()

    html_reports = list(outdir.rglob("multiqc_report.html"))
    result["has_multiqc_html"] = len(html_reports) > 0
    if html_reports:
        result["multiqc_html"] = str(html_reports[0])

    if pipeline == "rnaseq":
        result["pipeline_specific"]["has_star_salmon"] = (outdir / "star_salmon").exists()
        result["pipeline_specific"]["count_files"] = [str(p) for p in outdir.rglob("*.sf")][:5]
    elif pipeline == "atacseq":
        result["pipeline_specific"]["has_peaks"] = (outdir / "bwa" / "mergedLibrary" / "macs2").exists()
    elif pipeline == "sarek":
        result["pipeline_specific"]["has_variant_calling"] = (outdir / "variant_calling").exists()
    elif pipeline == "ampliseq":
        result["pipeline_specific"]["has_dada2"] = (outdir / "dada2").exists()
        result["pipeline_specific"]["has_taxonomy"] = (outdir / "taxonomy").exists()
        result["pipeline_specific"]["has_diversity"] = (outdir / "diversity").exists()

    return result
