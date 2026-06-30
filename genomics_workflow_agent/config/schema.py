"""Configuration schema and defaults for workflow runs."""

from __future__ import annotations

WORKFLOW_DEFAULTS: dict[str, dict] = {
    "fastq-qc": {
        "trim": False,
        "trimmer": "fastp",
    },
    "rnaseq": {
        "profile": "docker",
        "genome": None,
        "fasta": None,
        "gtf": None,
    },
    "atacseq": {
        "profile": "docker",
        "genome": None,
        "fasta": None,
        "gtf": None,
        "blacklist": None,
    },
    "amplicon": {
        "denoiser": "dada2",
        "taxonomy_db": "SILVA",
        "taxonomy_db_path": None,
        "primer_fw": None,
        "primer_rv": None,
        "profile": "docker",
    },
    "variant-qc": {
        "profile": "docker",
        "genome": None,
        "fasta": None,
        "known_sites": None,
    },
}

SUPPORTED_WORKFLOWS = list(WORKFLOW_DEFAULTS.keys())
