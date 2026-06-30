"""
inspect_file.py - Low-overhead genomic file inspector.

Detects file type, compression, size, dimensions (where safe), and likely assay.
Never reads entire large files into memory. Uses sampling and index-aware methods.
"""

import argparse
import gzip
import hashlib
import json
import os
import struct
import sys
from datetime import datetime
from pathlib import Path


try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

try:
    import anndata
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False


VERSION = "1.0.0"
TOOL_NAME = "inspect_file.py"

LARGE_FILE_BYTES = 50 * 1024 * 1024  
SAMPLE_LINES = 5


def file_size_mb(path: Path) -> float:
    return path.stat().st_size / (1024 * 1024)


def is_gzipped(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def is_bgzipped(path: Path) -> bool:
    """BGZF has extra field with BC tag."""
    try:
        with open(path, "rb") as f:
            header = f.read(18)
        if len(header) < 18:
            return False
        if header[:2] != b"\x1f\x8b":
            return False
        # Check BGZF extra field
        if header[3] & 4:
            xlen = struct.unpack("<H", header[10:12])[0]
            return xlen >= 6
        return False
    except Exception:
        return False


def safe_head(path: Path, n: int = SAMPLE_LINES, compressed: bool = False) -> list:
    """Read first n lines of a text file safely."""
    lines = []
    try:
        opener = gzip.open if compressed else open
        mode = "rt" if compressed else "r"
        with opener(path, mode, errors="replace") as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                lines.append(line.rstrip("\n"))
    except Exception as e:
        lines = [f"[ERROR reading file: {e}]"]
    return lines


def detect_extension(path: Path) -> tuple:
    """Return (base_ext, compressed) tuple."""
    name = path.name.lower()
    compressed = False
    if name.endswith(".gz") or name.endswith(".bz2") or name.endswith(".bgz"):
        compressed = True
        name = name.rsplit(".", 1)[0]
    ext = "." + name.rsplit(".", 1)[-1] if "." in name else ""
    return ext, compressed


def estimate_csv_dims(path: Path, compressed: bool, sep: str = ",") -> dict:
    """Estimate rows and columns by counting header cols and sampling row count."""
    try:
        opener = gzip.open if compressed else open
        mode = "rt" if compressed else "r"
        ncols = 0
        nrows = 0
        with opener(path, mode, errors="replace") as f:
            for i, line in enumerate(f):
                if i == 0:
                    ncols = len(line.split(sep))
                nrows += 1
                if nrows > 100_000:  
                    nrows = f">{nrows}"
                    break
        return {"rows": nrows, "cols": ncols}
    except Exception as e:
        return {"rows": "unknown", "cols": "unknown", "error": str(e)}


def inspect_fasta(path: Path, compressed: bool) -> dict:
    lines = safe_head(path, 10, compressed)
    n_headers = sum(1 for l in lines if l.startswith(">"))
    return {
        "type": "FASTA",
        "header_sample": lines[:4],
        "headers_in_sample": n_headers,
        "notes": "Use samtools faidx for indexing. Never load full genome into context.",
    }


def inspect_fastq(path: Path, compressed: bool) -> dict:
    lines = safe_head(path, 8, compressed)
    return {
        "type": "FASTQ",
        "header_sample": lines[:4],
        "notes": "Use FastQC or fastp for QC. Never load into context.",
    }


def inspect_vcf(path: Path, compressed: bool) -> dict:
    lines = safe_head(path, 20, compressed)
    meta = [l for l in lines if l.startswith("##")]
    header = next((l for l in lines if l.startswith("#CHROM")), None)
    samples = []
    if header:
        cols = header.split("\t")
        samples = cols[9:] if len(cols) > 9 else []
    return {
        "type": "VCF",
        "meta_lines_in_sample": len(meta),
        "header": header,
        "samples": samples,
        "n_samples": len(samples),
        "notes": "Use bcftools stats for full stats. Sample header only - never load full VCF.",
    }


def inspect_bed(path: Path, compressed: bool) -> dict:
    lines = safe_head(path, 5, compressed)
    first_data = next((l for l in lines if not l.startswith("#") and l.strip()), None)
    ncols = len(first_data.split("\t")) if first_data else "unknown"
    return {
        "type": "BED",
        "first_data_line": first_data,
        "n_cols": ncols,
        "notes": "BED3 has chrom/start/end. More columns may be name, score, strand, etc.",
    }


def inspect_gtf_gff(path: Path, compressed: bool, ext: str) -> dict:
    lines = safe_head(path, 10, compressed)
    meta = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#")]
    ftype = "GTF" if ext in (".gtf",) else "GFF"
    return {
        "type": ftype,
        "meta_lines_in_sample": len(meta),
        "first_data_lines": data[:3],
        "notes": "Required for TSS enrichment and gene annotation. Never load fully.",
    }


def inspect_mtx(path: Path) -> dict:
    lines = safe_head(path, 3, False)
    header = lines[0] if lines else ""
    dims_line = lines[1] if len(lines) > 1 else ""
    dims = dims_line.split() if dims_line and not dims_line.startswith("%") else []
    return {
        "type": "MTX (Market Exchange Format)",
        "header": header,
        "dimensions_line": dims_line,
        "dims_parsed": {"rows": dims[0], "cols": dims[1], "nnz": dims[2]} if len(dims) >= 3 else {},
        "notes": "Typically paired with barcodes.tsv and features.tsv in 10x directory.",
    }


def inspect_10x_dir(path: Path) -> dict:
    expected = {
        "matrix.mtx.gz": False,
        "barcodes.tsv.gz": False,
        "features.tsv.gz": False,
        "genes.tsv.gz": False,  # older format
        "matrix.mtx": False,
        "barcodes.tsv": False,
        "features.tsv": False,
    }
    found = []
    for fname in expected:
        fpath = path / fname
        if fpath.exists():
            expected[fname] = True
            found.append(fname)

    n_barcodes = "unknown"
    bc_file = path / "barcodes.tsv.gz"
    if not bc_file.exists():
        bc_file = path / "barcodes.tsv"
    if bc_file.exists():
        try:
            opener = gzip.open if str(bc_file).endswith(".gz") else open
            with opener(bc_file, "rt") as f:
                n_barcodes = sum(1 for _ in f)
        except Exception:
            pass

    return {
        "type": "10x Genomics directory",
        "files_found": found,
        "n_barcodes": n_barcodes,
        "is_complete": expected.get("matrix.mtx.gz") or expected.get("matrix.mtx"),
        "notes": "Use single-cell-rna-qc@life-sciences or scrna_qc_local.py for QC.",
    }


def inspect_h5ad(path: Path) -> dict:
    result = {"type": "h5ad (AnnData)", "notes": "scRNA-seq or other single-cell data."}

    if HAS_ANNDATA:
        try:
            import anndata as ad
            adata = ad.read_h5ad(path, backed="r")
            result["shape"] = list(adata.shape)
            result["n_obs"] = adata.n_obs
            result["n_vars"] = adata.n_vars
            result["obs_keys"] = list(adata.obs.columns)[:20]
            result["var_keys"] = list(adata.var.columns)[:20]
            result["obsm_keys"] = list(adata.obsm.keys())
            result["uns_keys"] = list(adata.uns.keys())[:10]
            adata.file.close()
        except Exception as e:
            result["error"] = str(e)
            result["notes"] += f" Could not read with anndata: {e}"
    elif HAS_H5PY:
        try:
            with h5py.File(path, "r") as f:
                result["h5_keys"] = list(f.keys())
                if "X" in f:
                    result["X_shape"] = list(f["X"].shape) if hasattr(f["X"], "shape") else "sparse"
                if "obs" in f:
                    result["obs_keys"] = list(f["obs"].keys())[:20]
                if "var" in f:
                    result["var_keys"] = list(f["var"].keys())[:20]
        except Exception as e:
            result["error"] = str(e)
    else:
        result["notes"] += " Install anndata or h5py for detailed inspection."

    return result


def inspect_bam_cram(path: Path, file_type: str) -> dict:
    result = {"type": file_type}

    if HAS_PYSAM:
        try:
            with pysam.AlignmentFile(str(path), "rb" if file_type == "BAM" else "rc") as f:
                header = f.header.to_dict()
                sq = header.get("SQ", [])
                result["n_references"] = len(sq)
                result["references_sample"] = [s["SN"] for s in sq[:5]]
                rg = header.get("RG", [])
                result["n_read_groups"] = len(rg)
                result["read_groups_sample"] = rg[:3]
                pg = header.get("PG", [])
                result["programs"] = [p.get("ID", "?") for p in pg]
                result["notes"] = f"{file_type} opened with pysam. Full stats: samtools flagstat."
        except Exception as e:
            result["error"] = str(e)
            result["notes"] = f"pysam failed: {e}. Try: samtools view -H {path}"
    else:
        result["notes"] = (
            f"pysam not available. Install pysam or run: samtools view -H {path}"
        )

    return result


def inspect_fragments(path: Path, compressed: bool) -> dict:
    lines = safe_head(path, 5, compressed)
    meta = [l for l in lines if l.startswith("#")]
    data = [l for l in lines if not l.startswith("#") and l.strip()]
    first = data[0].split("\t") if data else []
    return {
        "type": "Fragments file (scATAC)",
        "meta_lines_in_sample": len(meta),
        "first_data_cols": first,
        "expected_cols": "chrom, start, end, barcode, duplicates",
        "notes": "Use atac_qc_local.py for QC. Index with tabix if not already indexed.",
    }


def guess_assay(file_type: str, path: Path) -> str:
    name = path.name.lower()
    mapping = {
        "FASTQ": "Sequencing reads (assay unknown)",
        "FASTA": "Reference genome or sequences",
        "h5ad (AnnData)": "Single-cell (likely scRNA-seq, scATAC, or multiome)",
        "10x Genomics directory": "Single-cell (10x Genomics)",
        "MTX (Market Exchange Format)": "Single-cell count matrix",
        "VCF": "Variant calls (WGS/WES/targeted)",
        "BED": "Genomic intervals (peaks, regions, blacklist, etc.)",
        "GTF": "Gene annotation",
        "GFF": "Gene annotation",
        "BAM": "Aligned reads",
        "CRAM": "Aligned reads (compressed)",
        "Fragments file (scATAC)": "Single-cell ATAC-seq",
    }
    guess = mapping.get(file_type, "Unknown")
    if "atac" in name or "fragment" in name or "peak" in name:
        guess += " - likely ATAC-seq"
    if "rna" in name or "count" in name or "matrix" in name:
        guess += " - likely RNA-seq/scRNA-seq"
    if "wgs" in name or "wes" in name or "variant" in name or "snp" in name:
        guess += " - likely WGS/WES"
    return guess


def inspect_path(path: Path) -> dict:
    result = {
        "path": str(path),
        "exists": path.exists(),
        "size_mb": None,
        "compressed": False,
        "type": "unknown",
        "assay_guess": "unknown",
        "details": {},
        "warnings": [],
    }

    if not path.exists():
        result["warnings"].append("File does not exist.")
        return result

    if path.is_dir():
        result["type"] = "directory"
        result["size_mb"] = sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / (1024 * 1024)
        # Check if it's a 10x directory
        if (path / "matrix.mtx.gz").exists() or (path / "matrix.mtx").exists():
            result["details"] = inspect_10x_dir(path)
            result["type"] = "10x Genomics directory"
            result["assay_guess"] = guess_assay(result["type"], path)
        return result

    result["size_mb"] = round(file_size_mb(path), 2)
    if result["size_mb"] > LARGE_FILE_BYTES / (1024 * 1024):
        result["warnings"].append(f"Large file ({result['size_mb']:.1f} MB). Sampling only.")

    ext, compressed = detect_extension(path)
    result["compressed"] = compressed
    result["extension"] = ext

    bgz = is_bgzipped(path)
    gz = is_gzipped(path) and not bgz
    result["bgzipped"] = bgz
    result["gzipped"] = gz

    # Dispatch by extension
    if ext in (".fa", ".fna", ".fasta"):
        result["details"] = inspect_fasta(path, compressed)
        result["type"] = "FASTA"
    elif ext in (".fq", ".fastq"):
        result["details"] = inspect_fastq(path, compressed)
        result["type"] = "FASTQ"
    elif ext == ".vcf":
        result["details"] = inspect_vcf(path, compressed)
        result["type"] = "VCF"
    elif ext == ".bed":
        result["details"] = inspect_bed(path, compressed)
        result["type"] = "BED"
    elif ext == ".gtf":
        result["details"] = inspect_gtf_gff(path, compressed, ext)
        result["type"] = "GTF"
    elif ext in (".gff", ".gff3"):
        result["details"] = inspect_gtf_gff(path, compressed, ext)
        result["type"] = "GFF"
    elif ext == ".mtx":
        result["details"] = inspect_mtx(path)
        result["type"] = "MTX (Market Exchange Format)"
    elif ext == ".h5ad":
        result["details"] = inspect_h5ad(path)
        result["type"] = "h5ad (AnnData)"
    elif ext in (".bam",):
        result["details"] = inspect_bam_cram(path, "BAM")
        result["type"] = "BAM"
    elif ext in (".cram",):
        result["details"] = inspect_bam_cram(path, "CRAM")
        result["type"] = "CRAM"
    elif ext in (".csv",):
        dims = estimate_csv_dims(path, compressed, ",")
        result["details"] = {"dims": dims, "delimiter": ",", "header_sample": safe_head(path, 3, compressed)}
        result["type"] = "CSV"
    elif ext in (".tsv", ".txt"):
        # Check if it's a fragments file
        name_lower = path.name.lower()
        if "fragment" in name_lower:
            result["details"] = inspect_fragments(path, compressed)
            result["type"] = "Fragments file (scATAC)"
        else:
            dims = estimate_csv_dims(path, compressed, "\t")
            result["details"] = {"dims": dims, "delimiter": "\\t", "header_sample": safe_head(path, 3, compressed)}
            result["type"] = "TSV"
    else:
        result["details"] = {"header_sample": safe_head(path, 3, compressed)}
        result["warnings"].append(f"Unknown extension '{ext}'. Showing header sample only.")

    result["assay_guess"] = guess_assay(result["type"], path)
    return result


def build_markdown(results: list, args) -> str:
    lines = ["# File Inspection Report", f"\nGenerated: {datetime.now().isoformat()}\n"]
    lines.append("## File Inventory\n")
    lines.append("| File | Type | Size (MB) | Compressed | Assay Guess |")
    lines.append("|------|------|-----------|------------|-------------|")
    for r in results:
        fname = Path(r["path"]).name
        size = r.get("size_mb", "?")
        compressed = "yes" if r.get("compressed") or r.get("bgzipped") or r.get("gzipped") else "no"
        lines.append(f"| {fname} | {r['type']} | {size} | {compressed} | {r['assay_guess']} |")

    lines.append("\n## Details\n")
    for r in results:
        lines.append(f"### {Path(r['path']).name}\n")
        lines.append(f"- **Path**: `{r['path']}`")
        lines.append(f"- **Type**: {r['type']}")
        lines.append(f"- **Size**: {r.get('size_mb', '?')} MB")
        if r.get("warnings"):
            for w in r["warnings"]:
                lines.append(f"- **WARNING**: {w}")
        details = r.get("details", {})
        for k, v in details.items():
            if k not in ("header_sample", "first_data_lines", "first_data_cols"):
                lines.append(f"- **{k}**: {v}")
        if "header_sample" in details:
            lines.append("\n**Header sample:**\n```")
            for line in details["header_sample"]:
                lines.append(line)
            lines.append("```")
        lines.append("")

    lines.append("## Missing Metadata Checklist\n")
    lines.append("Review the following for each file/dataset:\n")
    for item in ["Species", "Genome build", "Tissue", "Disease/condition", "Batch labels",
                 "Sample identifiers", "Protocol/chemistry", "Cell barcode format"]:
        lines.append(f"- [ ] {item}")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Low-overhead genomic file inspector.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, nargs="+", help="File(s) or directory to inspect.")
    parser.add_argument("--output-dir", default=".", help="Directory for output files.")
    parser.add_argument("--json", action="store_true", help="Write JSON summary.")
    parser.add_argument("--markdown", action="store_true", help="Write Markdown summary.")
    parser.add_argument("--verbose", action="store_true", help="Print full details to stdout.")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for inp in args.input:
        path = Path(inp)
        result = inspect_path(path)
        results.append(result)

    summary = {
        "tool": TOOL_NAME,
        "version": VERSION,
        "python_version": sys.version,
        "pysam_available": HAS_PYSAM,
        "h5py_available": HAS_H5PY,
        "anndata_available": HAS_ANNDATA,
        "generated": datetime.now().isoformat(),
        "files": results,
    }

    if args.verbose or not (args.json or args.markdown):
        for r in results:
            print(f"\n{'='*60}")
            print(f"File: {r['path']}")
            print(f"Type: {r['type']}")
            print(f"Size: {r.get('size_mb', '?')} MB")
            print(f"Assay guess: {r['assay_guess']}")
            if r.get("warnings"):
                for w in r["warnings"]:
                    print(f"WARNING: {w}")
            if r.get("details"):
                import pprint
                pprint.pprint(r["details"])

    if args.json or True:  # always write JSON
        json_path = out_dir / "file_inventory.json"
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"JSON written: {json_path}")

    if args.markdown or True:  # always write markdown
        md = build_markdown(results, args)
        md_path = out_dir / "file_inventory.md"
        with open(md_path, "w") as f:
            f.write(md)
        print(f"Markdown written: {md_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
