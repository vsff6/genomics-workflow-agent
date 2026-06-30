from __future__ import annotations

from pathlib import Path
from typing import Any

from genomics_workflow_agent.tools.files import find_fastq_pairs
from genomics_workflow_agent.tools.runner import STATUS_SKIPPED, run_command
from genomics_workflow_agent.tools.versions import check_tools

REQUIRED_TOOLS = ["fastqc", "multiqc"]
OPTIONAL_TOOLS = ["fastp", "cutadapt"]

BIOLOGICAL_CAVEATS = [
    "Per-base quality scores reflect library chemistry and instrument, not biology directly.",
    "Low-quality tails are common for longer reads — trimming thresholds must match downstream tool requirements.",
    "Adapter content depends on insert size — very short inserts are expected in some protocols (e.g. ATAC-seq).",
    "High duplication rates may reflect library complexity, PCR amplification, or genuine biological signal.",
    "GC content outliers may indicate contamination, AT-rich organisms, or biologically meaningful regions.",
]


def _fastqc_stem(filename: str) -> str:
    # FastQC strips .gz then .fastq/.fq to produce the output stem
    name = Path(filename).name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _fastqc_expected_outputs(fastq_files: list[str], fastqc_outdir: str) -> list[str]:
    outdir = Path(fastqc_outdir)
    expected = []
    for f in fastq_files:
        stem = _fastqc_stem(f)
        expected.append(str(outdir / f"{stem}_fastqc.html"))
        expected.append(str(outdir / f"{stem}_fastqc.zip"))
    return expected


def plan(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    trim: str | bool | None = None,
    trimmer: str = "fastp",
    dry_run: bool = True,
) -> dict[str, Any]:
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    if isinstance(trim, str) and trim not in (False, "false", "no", "none", ""):
        trimmer = trim
        do_trim = True
    elif trim is True:
        do_trim = True
    else:
        do_trim = False

    tools = check_tools(REQUIRED_TOOLS + OPTIONAL_TOOLS)
    pairs = find_fastq_pairs(input_dir) if input_dir.is_dir() else []
    all_fastqs = sorted({p["r1"] for p in pairs} | {p["r2"] for p in pairs if p["r2"]})

    steps: list[dict] = []
    warnings: list[str] = []
    skipped: list[dict] = []
    fastqc_outdir = str(output_dir / "fastqc")

    if tools["fastqc"]["available"]:
        if all_fastqs:
            cmd = ["fastqc", "--outdir", fastqc_outdir, "--threads", "4"] + all_fastqs
            expected = _fastqc_expected_outputs(all_fastqs, fastqc_outdir)
        else:
            cmd = ["fastqc", "--outdir", fastqc_outdir, "--threads", "4",
                   str(input_dir / "*.fastq.gz")]
            expected = []
            warnings.append("No FASTQ files detected — FastQC command uses a glob placeholder")
        steps.append({
            "name": "fastqc",
            "description": "Per-file FASTQ quality control",
            "command": cmd,
            "output_dir": fastqc_outdir,
            "expected_outputs": expected,
            "dry_run": dry_run,
            "required_tools": ["fastqc"],
        })
    else:
        skipped.append({
            "step": "fastqc",
            "reason": "fastqc not found in PATH",
            "status": STATUS_SKIPPED,
            "install": "conda install -c bioconda fastqc",
        })

    multiqc_outdir = str(output_dir / "multiqc")
    if tools["multiqc"]["available"]:
        multiqc_input = fastqc_outdir if tools["fastqc"]["available"] else str(output_dir)
        steps.append({
            "name": "multiqc",
            "description": "Aggregate QC report across all samples",
            "command": ["multiqc", multiqc_input, "--outdir", multiqc_outdir, "--force"],
            "output_dir": multiqc_outdir,
            "expected_outputs": [
                str(Path(multiqc_outdir) / "multiqc_report.html"),
                str(Path(multiqc_outdir) / "multiqc_data"),
            ],
            "dry_run": dry_run,
            "required_tools": ["multiqc"],
        })
    else:
        skipped.append({
            "step": "multiqc",
            "reason": "multiqc not found in PATH",
            "status": STATUS_SKIPPED,
            "install": "pip install multiqc",
        })

    if do_trim:
        trim_steps, trim_skipped, trim_warnings = _build_trim_steps(
            trimmer, tools, pairs, output_dir, dry_run
        )
        steps.extend(trim_steps)
        skipped.extend(trim_skipped)
        warnings.extend(trim_warnings)

    return {
        "workflow": "fastq-qc",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "dry_run": dry_run,
        "trimming_requested": do_trim,
        "trimmer": trimmer if do_trim else None,
        "samples_detected": len(pairs),
        "fastq_files": all_fastqs,
        "tool_status": tools,
        "steps": steps,
        "skipped_steps": skipped,
        "warnings": warnings,
        "biological_caveats": BIOLOGICAL_CAVEATS,
        "limitations": [
            "FastQC does not detect all adapter types — check MultiQC adapter contamination plot.",
            "Trimming is not applied automatically unless --trim is passed.",
            "Trimming parameters depend on library type and downstream requirements.",
        ],
    }


def _build_trim_steps(
    trimmer: str,
    tools: dict,
    pairs: list[dict],
    output_dir: Path,
    dry_run: bool,
) -> tuple[list[dict], list[dict], list[str]]:
    steps: list[dict] = []
    skipped: list[dict] = []
    warnings: list[str] = []

    chosen = trimmer if tools.get(trimmer, {}).get("available") else None
    if chosen is None:
        alt = "cutadapt" if trimmer == "fastp" else "fastp"
        if tools.get(alt, {}).get("available"):
            chosen = alt
            warnings.append(f"{trimmer} not found; using {alt} instead")

    if chosen is None:
        skipped.append({
            "step": "trimming",
            "reason": f"Neither {trimmer} nor alternative trimmer found in PATH",
            "status": STATUS_SKIPPED,
            "install": "conda install -c bioconda fastp  OR  pip install cutadapt",
        })
        return steps, skipped, warnings

    trim_dir = output_dir / "trimmed"

    for pair in pairs:
        sample = pair["sample"]
        r1 = pair["r1"]
        r2 = pair.get("r2")
        paired = pair.get("paired", False)

        r1_out = str(trim_dir / Path(r1).name)
        r2_out = str(trim_dir / Path(r2).name) if r2 else None
        expected = [r1_out]
        if r2_out:
            expected.append(r2_out)

        if chosen == "fastp":
            cmd = ["fastp", "-i", r1, "-o", r1_out, "--thread", "4",
                   "--json", str(trim_dir / f"{sample}_fastp.json"),
                   "--html", str(trim_dir / f"{sample}_fastp.html")]
            if paired and r2:
                cmd += ["-I", r2, "-O", r2_out]
            expected += [
                str(trim_dir / f"{sample}_fastp.json"),
                str(trim_dir / f"{sample}_fastp.html"),
            ]
        else:
            if paired and r2:
                cmd = ["cutadapt", "-o", r1_out, "-p", r2_out, r1, r2]
            else:
                cmd = ["cutadapt", "-o", r1_out, r1]

        steps.append({
            "name": f"trim_{chosen}_{sample}",
            "description": f"Trim adapters/{chosen} — {sample}",
            "command": cmd,
            "output_dir": str(trim_dir),
            "expected_outputs": expected,
            "dry_run": dry_run,
            "required_tools": [chosen],
        })

    return steps, skipped, warnings


def execute(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    trim: str | bool | None = None,
    trimmer: str = "fastp",
    provenance_dir: Path | None = None,
    timeout: int = 3600,
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    wf_plan = plan(input_dir, output_dir, trim=trim, trimmer=trimmer, dry_run=False)

    for step in wf_plan["steps"]:
        step_outdir = step.get("output_dir")
        if step_outdir:
            Path(step_outdir).mkdir(parents=True, exist_ok=True)

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
        )
        step["execution"] = record
        step_results.append(record)

    wf_plan["step_results"] = step_results
    wf_plan["dry_run"] = False
    return wf_plan
