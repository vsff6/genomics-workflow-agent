"""
Amplicon/microbiome workflow — planning and execution via nf-core/ampliseq.

Execution priority:
  1. nf-core/ampliseq (preferred — fully validated pipeline)
  2. QIIME2 direct (if available and ampliseq not requested)
  3. R/DADA2 (if available and QIIME2 not available)

Execution requires nextflow in PATH for nf-core/ampliseq.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from genomics_workflow_agent.tools.nextflow import build_nextflow_cmd, run_nextflow
from genomics_workflow_agent.tools.runner import STATUS_SKIPPED
from genomics_workflow_agent.tools.samplesheets import build_amplicon_samplesheet
from genomics_workflow_agent.tools.versions import check_tools

OPTIONAL_TOOLS = ["nextflow", "qiime", "Rscript", "fastqc", "multiqc", "fastp", "cutadapt",
                  "docker", "singularity", "conda"]

TAXONOMY_DATABASES = {
    "SILVA": {
        "description": "SILVA SSU/LSU rRNA database (16S, 18S, 23S/28S)",
        "url": "https://www.arb-silva.de/",
        "versions": ["138.1", "138", "132"],
        "qiime2_classifier": "silva-138-99-seqs-classifier.qza",
        "ampliseq_param": "silva=2.13",
    },
    "GTDB": {
        "description": "Genome Taxonomy Database (prokaryotic phylogeny)",
        "url": "https://gtdb.ecogenomic.org/",
        "versions": ["r220", "r214", "r207"],
        "qiime2_classifier": "gtdb-ssu-r220-classifier.qza",
        "ampliseq_param": "gtdb=r220",
    },
    "UNITE": {
        "description": "Fungal ITS database (ITS1, ITS2, full ITS)",
        "url": "https://unite.ut.ee/",
        "versions": ["9.0", "8.3"],
        "qiime2_classifier": "unite-ver9-seqs-classifier.qza",
        "ampliseq_param": "unite-fungi=9.0",
    },
    "Greengenes2": {
        "description": "Greengenes2 (16S, 16S+ITS)",
        "url": "http://greengenes2.ucsd.edu/",
        "versions": ["2022.10", "2024.09"],
        "qiime2_classifier": "gg-2022.10.backbone.full-length.nb.sklearn-1.4.2.qza",
        "ampliseq_param": "greengenes2=2022.10",
    },
    "custom": {
        "description": "User-provided custom database",
        "url": None,
        "versions": ["user-defined"],
        "qiime2_classifier": "<path-to-classifier.qza>",
        "ampliseq_param": "<path-to-fasta>",
    },
}

BIOLOGICAL_CAVEATS = [
    "16S rRNA V-region primer choice determines taxonomic resolution and community biases.",
    "DADA2 error models are run-specific — denoise each sequencing run separately before merging.",
    "Rarefaction discards reads and introduces variance — evaluate necessity from depth distribution.",
    "Alpha diversity metrics (Shannon, Faith PD) have different sensitivity to rare taxa — report multiple.",
    "Beta diversity ordination (Bray-Curtis, UniFrac) is sensitive to data transformation and normalisation.",
    "Taxonomy confidence thresholds affect rare taxon detection — do not use generic defaults without validation.",
    "Phylogenetic diversity metrics (Faith PD, UniFrac) require a phylogenetic tree.",
    "Amplicon sequencing is not metagenomics — functional inference from taxonomy is indirect.",
    "Database version and taxonomic nomenclature affect reproducibility — always record version.",
]


def plan(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    primer_fw: str | None = None,
    primer_rv: str | None = None,
    taxonomy_db: str = "SILVA",
    taxonomy_db_path: str | None = None,
    denoiser: str = "dada2",
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    dry_run: bool = True,
) -> dict[str, Any]:
    """Return an amplicon workflow execution plan (defaults to nf-core/ampliseq)."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    pipeline_outdir = str(output_dir / "ampliseq_output")

    tools = check_tools(OPTIONAL_TOOLS)
    warnings: list[str] = []
    blockers: list[str] = []

    ss_path = output_dir / "samplesheet_amplicon.csv"
    samplesheet_result = build_amplicon_samplesheet(input_dir, ss_path)

    db_info = TAXONOMY_DATABASES.get(taxonomy_db, TAXONOMY_DATABASES["custom"])
    taxonomy_param = taxonomy_db_path or db_info.get("ampliseq_param", "silva=2.13")

    if not primer_fw or not primer_rv:
        warnings.append(
            "Primer sequences not provided. nf-core/ampliseq requires --FW_primer and --RV_primer. "
            "Common 16S V4: 515F=GTGYCAGCMGCCGCGGTAA / 806R=GGACTACNVGGGTWTCTAAT"
        )

    if not taxonomy_db_path and taxonomy_db != "custom":
        warnings.append(
            f"Taxonomy database path not provided. Using built-in ampliseq parameter '{taxonomy_param}'. "
            f"Download {taxonomy_db}: {db_info.get('url', 'see documentation')}"
        )

    if not tools["nextflow"]["available"]:
        blockers.append("nextflow not found in PATH. Install: https://www.nextflow.io")

    container_ok = tools["docker"]["available"] or tools["singularity"]["available"]
    if not container_ok and not tools["conda"]["available"]:
        warnings.append("No container runtime or conda found for Nextflow")

    cmd = build_nextflow_cmd(
        "ampliseq",
        input_path=str(ss_path) if samplesheet_result.get("created") else None,
        outdir=pipeline_outdir,
        primer_fw=primer_fw,
        primer_rv=primer_rv,
        taxonomy_param=taxonomy_param,
        profile=profile,
        resume=resume,
        extra_args=extra_args,
    )

    nfcore_step = {
        "name": "nextflow_ampliseq",
        "description": "Run nf-core/ampliseq via Nextflow (DADA2 + taxonomy + diversity)",
        "command": cmd,
        "output_dir": pipeline_outdir,
        "expected_outputs": [pipeline_outdir],
        "dry_run": dry_run,
        "required_tools": ["nextflow"],
    }

    # Planning steps for reference (always shown, even in execute mode)
    planning_steps = [
        {
            "name": "feature_table_filter",
            "description": "Filter feature table: remove low-frequency ASVs, contaminants, mitochondria",
            "note": "After ampliseq: filter in QIIME2 or phyloseq. Min frequency 10, min samples 2.",
            "command": None,
            "dry_run": True,
        },
        {
            "name": "normalisation",
            "description": "Normalisation strategy",
            "note": "Options: rarefaction, CSS (metagenomeSeq), DESeq2 VST, CLR. Choose based on analysis goals.",
            "options": {
                "rarefaction": "qiime diversity alpha-rarefaction",
                "css": "phyloseq::normalizeSampleCounts",
                "clr": "vegan::decostand(method='clr')",
            },
            "command": None,
            "dry_run": True,
        },
        {
            "name": "alpha_diversity",
            "description": "Alpha diversity (within-sample richness and evenness)",
            "note": "nf-core/ampliseq computes Shannon, Observed features, Faith PD, Pielou evenness automatically.",
            "metrics": ["shannon", "observed_features", "faith_pd", "pielou_e"],
            "command": None,
            "dry_run": True,
        },
        {
            "name": "beta_diversity",
            "description": "Beta diversity (between-sample community composition)",
            "note": "nf-core/ampliseq computes Bray-Curtis, weighted/unweighted UniFrac, Jaccard automatically.",
            "metrics": ["bray_curtis", "weighted_unifrac", "unweighted_unifrac", "jaccard"],
            "command": None,
            "dry_run": True,
        },
    ]

    return {
        "workflow": "amplicon",
        "pipeline": "nf-core/ampliseq",
        "denoiser": denoiser,
        "taxonomy_db": taxonomy_db,
        "taxonomy_db_info": db_info,
        "taxonomy_db_path": taxonomy_db_path,
        "taxonomy_param": taxonomy_param,
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "pipeline_outdir": pipeline_outdir,
        "dry_run": dry_run,
        "primers": {"forward": primer_fw, "reverse": primer_rv},
        "samplesheet": samplesheet_result,
        "steps": [nfcore_step] + planning_steps,
        "nfcore_ampliseq_command": cmd,
        "nfcore_ampliseq_command_str": " ".join(cmd),
        "command": cmd,
        "command_str": " ".join(cmd),
        "tool_status": tools,
        "blockers": blockers,
        "warnings": warnings,
        "biological_caveats": BIOLOGICAL_CAVEATS,
        "next_actions": [
            "Confirm primer sequences from experimental records",
            "Review samplesheet before execution",
            "Run FastQC on raw reads before submitting to pipeline",
            "After pipeline: review DADA2 stats, taxonomy bar plots, diversity outputs",
            "Invoke biology-interpretation-reviewer on diversity outputs",
        ],
        "limitations": [
            "nf-core/ampliseq is the primary execution path. Direct QIIME2 execution is planned only.",
            "Taxonomy databases must be downloaded separately if providing a local path.",
            "Phylogenetic diversity requires a tree — ampliseq produces one via QIIME2 internally.",
            "This framework does not replace QIIME2, DADA2, or nf-core/ampliseq.",
        ],
    }


def execute(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    primer_fw: str | None = None,
    primer_rv: str | None = None,
    taxonomy_db: str = "SILVA",
    taxonomy_db_path: str | None = None,
    denoiser: str = "dada2",
    profile: str = "docker",
    resume: bool = False,
    extra_args: list[str] | None = None,
    provenance_dir: Path | None = None,
    timeout: int = 86400,
) -> dict[str, Any]:
    """Execute nf-core/ampliseq. Requires nextflow in PATH."""
    output_dir = Path(output_dir)
    pipeline_outdir = output_dir / "ampliseq_output"

    wf_plan = plan(
        input_dir, output_dir,
        primer_fw=primer_fw, primer_rv=primer_rv,
        taxonomy_db=taxonomy_db, taxonomy_db_path=taxonomy_db_path,
        denoiser=denoiser, profile=profile, resume=resume,
        extra_args=extra_args, dry_run=False,
    )

    if wf_plan["blockers"]:
        wf_plan["step_results"] = [{
            "label": "nextflow_ampliseq", "status": STATUS_SKIPPED,
            "error": f"Blocked: {wf_plan['blockers']}", "executed": False,
        }]
        return wf_plan

    record = run_nextflow(
        "ampliseq", wf_plan["command"],
        outdir=pipeline_outdir, dry_run=False,
        provenance_dir=provenance_dir, timeout=timeout,
    )

    wf_plan["step_results"] = [record]
    wf_plan["dry_run"] = False
    if record.get("pipeline_output_validation"):
        wf_plan["pipeline_output_validation"] = record["pipeline_output_validation"]
    return wf_plan
