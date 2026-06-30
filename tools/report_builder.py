"""
report_builder.py - Combine Markdown report sections into a final analysis report.

Assembles outputs from file inspection, QC tools, nf-core summaries, and
biological interpretation into a single reproducible Markdown report.
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

VERSION = "1.0.0"
TOOL_NAME = "report_builder.py"
TEMPLATE_PATH = Path(__file__).parent.parent / "templates" / "report_template.md"


def load_json_safe(path: Path) -> dict:
    try:
        with open(path) as f:
            return json.load(f)
    except Exception as e:
        return {"error": str(e), "path": str(path)}


def load_md_section(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except Exception as e:
        return f"[Could not load {path}: {e}]"


def collect_sections(args) -> dict:
    sections = {}

    if args.inspect_dir:
        d = Path(args.inspect_dir)
        md_file = d / "file_inventory.md"
        json_file = d / "file_inventory.json"
        sections["file_inventory"] = load_md_section(md_file) if md_file.exists() else None
        sections["file_inventory_json"] = load_json_safe(json_file) if json_file.exists() else {}

    if args.scrna_qc_dir:
        d = Path(args.scrna_qc_dir)
        md_file = d / "qc_report.md"
        json_file = d / "summary.json"
        sections["scrna_qc"] = load_md_section(md_file) if md_file.exists() else None
        sections["scrna_qc_json"] = load_json_safe(json_file) if json_file.exists() else {}

    if args.atac_qc_dir:
        d = Path(args.atac_qc_dir)
        md_file = d / "atac_qc_report.md"
        json_file = d / "atac_qc_summary.json"
        sections["atac_qc"] = load_md_section(md_file) if md_file.exists() else None
        sections["atac_qc_json"] = load_json_safe(json_file) if json_file.exists() else {}

    if args.wgs_qc_dir:
        d = Path(args.wgs_qc_dir)
        md_file = d / "wgs_vcf_qc_report.md"
        json_file = d / "wgs_vcf_qc_summary.json"
        sections["wgs_qc"] = load_md_section(md_file) if md_file.exists() else None
        sections["wgs_qc_json"] = load_json_safe(json_file) if json_file.exists() else {}

    if args.reference_check_dir:
        d = Path(args.reference_check_dir)
        md_file = d / "reference_validation.md"
        sections["reference_check"] = load_md_section(md_file) if md_file.exists() else None

    if args.bio_interpretation_file:
        sections["bio_interpretation"] = load_md_section(Path(args.bio_interpretation_file))

    return sections


def build_report(args, sections: dict) -> str:
    title = args.title or "Genomics Analysis Report"
    timestamp = datetime.now().isoformat()

    lines = [
        f"# {title}",
        f"\nGenerated: {timestamp}",
        f"\n**Tool**: `{TOOL_NAME}` v{VERSION}",
        "\n---",
        "\n## Dataset Summary\n",
        f"- **Analysis date**: {timestamp}",
        f"- **Species**: {args.species or 'not specified'}",
        f"- **Tissue**: {args.tissue or 'not specified'}",
        f"- **Disease/condition**: {args.disease or 'not specified'}",
        f"- **Genome build**: {args.genome_build or 'not specified'}",
        f"- **Protocol**: {args.protocol or 'not specified'}",
        f"- **Notes**: {args.notes or 'none'}",
    ]

    lines.append("\n## File Provenance\n")
    if sections.get("file_inventory"):
        lines.append(sections["file_inventory"])
    else:
        lines.append("_File inventory not provided. Run `python tools/inspect_file.py` first._")

    lines.append("\n## Metadata and Assumptions\n")
    lines.append("Review assumptions documented in each QC tool section below.")
    lines.append("If metadata was missing, tools proceed with conservative labeled assumptions.")

    lines.append("\n## Official Claude Life Sciences Tools Used\n")
    if args.official_tools_used:
        for t in args.official_tools_used.split(","):
            lines.append(f"- {t.strip()}")
    else:
        lines.append("- None reported (local fallback tools used).")
    lines.append("\n> If `single-cell-rna-qc@life-sciences`, `scvi-tools@life-sciences`, or `nextflow-development@life-sciences` were used, add them here.")

    lines.append("\n## Local Tools Used\n")
    tools_used = []
    if sections.get("scrna_qc"):
        tools_used.append("`tools/scrna_qc_local.py`")
    if sections.get("atac_qc"):
        tools_used.append("`tools/atac_qc_local.py`")
    if sections.get("wgs_qc"):
        tools_used.append("`tools/wgs_vcf_qc_local.py`")
    if sections.get("reference_check"):
        tools_used.append("`tools/reference_validator.py`")
    if not tools_used:
        tools_used = ["_No local QC tool outputs provided._"]
    for t in tools_used:
        lines.append(f"- {t}")

    lines.append("\n## Genome Build and Annotation Version\n")
    lines.append(f"- **Genome build**: {args.genome_build or '[NOT SPECIFIED]'}")
    lines.append(f"- **Annotation version**: {args.annotation_version or '[NOT SPECIFIED]'}")
    lines.append("\n> Never assume genome build. Always confirm from file headers, alignment parameters, or provider metadata.")

    lines.append("\n## Reference Files Used\n")
    if sections.get("reference_check"):
        lines.append(sections["reference_check"])
    else:
        lines.append("_No reference validation report provided._")

    lines.append("\n## QC Metrics\n")
    if sections.get("scrna_qc"):
        lines.append("### scRNA-seq QC\n")
        lines.append(sections["scrna_qc"])
    if sections.get("atac_qc"):
        lines.append("\n### ATAC-seq QC\n")
        lines.append(sections["atac_qc"])
    if sections.get("wgs_qc"):
        lines.append("\n### WGS/VCF QC\n")
        lines.append(sections["wgs_qc"])
    if not any([sections.get("scrna_qc"), sections.get("atac_qc"), sections.get("wgs_qc")]):
        lines.append("_No QC metrics sections provided._")

    lines.append("\n## Plots Generated\n")
    for key in ["scrna_qc_json", "atac_qc_json", "wgs_qc_json"]:
        data = sections.get(key, {})
        plots = data.get("plots", [])
        for p in plots:
            if p:
                lines.append(f"- `{p}`")

    lines.append("\n## Recommended Filtering Parameters\n")
    lines.append("See individual QC sections above for filter recommendations.")
    lines.append("All filter recommendations are advisory only. Apply only after biological review.")

    lines.append("\n## Biological Justification\n")
    lines.append("For every proposed filter, the QC tools document:")
    lines.append("1. What was observed?")
    lines.append("2. What technical artifact could explain it?")
    lines.append("3. What biological state could also explain it?")
    lines.append("4. What metadata would help distinguish artifact from biology?")
    lines.append("5. What validation should be done?")
    lines.append("6. Should data be filtered, flagged, stratified, or preserved?")

    lines.append("\n## Technical Artifact vs. Plausible Biology\n")
    if sections.get("bio_interpretation"):
        lines.append(sections["bio_interpretation"])
    else:
        lines.append("_Biological interpretation not yet provided._")
        lines.append("\n> Run the `biology-interpretation-reviewer` agent or `biological-interpretation-report` skill to generate this section.")
        lines.append("\n| Observation | Possible Technical Explanation | Possible Biological Explanation | Evidence Supporting Artifact | Evidence Supporting Biology | Recommended Follow-up | Confidence |")
        lines.append("|-------------|-------------------------------|--------------------------------|------------------------------|-----------------------------|-----------------------|------------|")
        lines.append("| _pending_ | | | | | | |")

    lines.append("\n## Skipped Metrics and Missing Inputs\n")
    for key in ["scrna_qc_json", "atac_qc_json", "wgs_qc_json"]:
        data = sections.get(key, {})
        skipped = data.get("skipped_metrics", [])
        if skipped:
            lines.append(f"\n**From {key.replace('_json','')}:**")
            for s in skipped:
                lines.append(f"- {s.get('metric', '?')}: {s.get('reason', '?')}")

    lines.append("\n## Limitations\n")
    lines.append("- Results are from local QC tools and may not replace full nf-core pipeline outputs.")
    lines.append("- Official Claude Life Sciences skills should be used when available and compatible.")
    lines.append("- No clinical or pathogenicity interpretation is provided.")
    lines.append("- Doublet detection and ambient RNA estimation not included in local scRNA QC.")
    lines.append("- Coverage calculations require external tools (mosdepth, samtools depth).")
    lines.append("- Filtering recommendations require biological validation before application.")

    lines.append("\n## Suggested Next Analyses\n")
    lines.append("_To be populated by the `biology-interpretation-reviewer` agent._")
    lines.append("- [ ] Doublet detection (scrublet, DoubletFinder)")
    lines.append("- [ ] Ambient RNA estimation (SoupX, DecontX)")
    lines.append("- [ ] Marker gene validation for proposed clusters")
    lines.append("- [ ] Reference mapping or label transfer (if reference available)")
    lines.append("- [ ] Sensitivity analysis for filter thresholds")

    lines.append("\n## Commands Run\n")
    lines.append("```bash")
    lines.append(f"# Report generated: {timestamp}")
    lines.append(f"python tools/report_builder.py {' '.join(sys.argv[1:])}")
    lines.append("```")

    lines.append("\n## Software Versions\n")
    lines.append(f"- Python: {sys.version}")
    lines.append(f"- `{TOOL_NAME}`: v{VERSION}")

    for key, label in [("scrna_qc_json", "scrna_qc_local.py"), ("atac_qc_json", "atac_qc_local.py"),
                        ("wgs_qc_json", "wgs_vcf_qc_local.py")]:
        data = sections.get(key, {})
        if data.get("version"):
            lines.append(f"- `{label}`: v{data['version']}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Combine QC tool outputs into a final Markdown report.",
    )
    parser.add_argument("--title", default="Genomics QC Report", help="Report title")
    parser.add_argument("--species", help="Species")
    parser.add_argument("--tissue", help="Tissue")
    parser.add_argument("--disease", help="Disease/condition")
    parser.add_argument("--genome-build", help="Genome build")
    parser.add_argument("--annotation-version", help="Annotation version")
    parser.add_argument("--protocol", help="Protocol")
    parser.add_argument("--notes", help="Additional notes")
    parser.add_argument("--inspect-dir", help="Output dir from inspect_file.py")
    parser.add_argument("--scrna-qc-dir", help="Output dir from scrna_qc_local.py")
    parser.add_argument("--atac-qc-dir", help="Output dir from atac_qc_local.py")
    parser.add_argument("--wgs-qc-dir", help="Output dir from wgs_vcf_qc_local.py")
    parser.add_argument("--reference-check-dir", help="Output dir from reference_validator.py")
    parser.add_argument("--bio-interpretation-file", help="Markdown file with biological interpretation section")
    parser.add_argument("--official-tools-used", help="Comma-separated list of official skills used")
    parser.add_argument("--output-dir", default="reports", help="Output directory")
    parser.add_argument("--output-name", default="final_report.md", help="Output filename")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sections = collect_sections(args)
    report = build_report(args, sections)

    report_path = out_dir / args.output_name
    with open(report_path, "w", encoding="utf-8", errors="replace") as f:
        f.write(report)

    print(f"Report written: {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
