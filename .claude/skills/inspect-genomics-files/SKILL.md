# Skill: inspect-genomics-files

## Purpose
Low-overhead inspection of genomic files. Produces a file inventory, identifies types and likely assays, and recommends the next workflow. No official Claude Life Sciences skill replaces this; it is always built locally.

## When to use
- At the start of any genomics session when file types are unknown.
- Before delegating to any QC specialist or pipeline.
- When the user says "what files do I have?" or "inspect this directory."

## Official tool to prefer
None. This skill has no official Life Sciences equivalent.

## Required inputs
- File path or directory path

## Optional inputs
- `--output-dir` for saving the inventory
- `--verbose` for additional header/schema detail

## Workflow steps

1. **Enumerate files**: Use Glob to list all files in the target path.
2. **For each file**:
   - Check extension and compression
   - Get file size
   - Read first 5–10 lines only for text files
   - Use `tools/inspect_file.py` for structured inspection
3. **Build inventory table**: path, type, size, compression, dimensions, assay guess, notes
4. **Check metadata signals**: chromosome names, gene ID style, barcode structure, sample labels
5. **Generate missing metadata checklist**
6. **Recommend next workflow** based on detected file types

## Command
```bash
python tools/inspect_file.py \
  --input <path> \
  --output-dir <output_dir> \
  --json \
  --markdown
```

## Expected outputs
- `file_inventory.json` - machine-readable inventory
- `file_inventory.md` - human-readable inventory table
- Missing metadata checklist
- Recommended next workflow

## Failure modes
- Binary files (BAM, CRAM, h5ad) cannot be sampled with text tools; always use `inspect_file.py`
- Compressed files need decompression-aware reading
- MTX files require paired barcodes/features files

## Reproducibility requirements
- Record file paths, sizes, and modification timestamps
- Record tool version and Python version
- Save JSON output to output directory

## Biological reasoning checklist
- [ ] Species identified from chromosome names or gene IDs?
- [ ] Genome build clues (chromosome sizes, contig names)?
- [ ] Assay type confirmed from header/content?
- [ ] Sample identifiers present?
- [ ] Condition/batch labels present?
- [ ] Barcode format consistent with expected protocol?
- [ ] Feature annotation style (Ensembl vs gene symbol)?
