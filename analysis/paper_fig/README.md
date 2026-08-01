# Paper Figure Scripts

These scripts reproduce publication figures from `simitall` output files. They
call the package API directly and export the source data used by every panel.
Generated results belong under `analysis/results/`, which is ignored by Git.

## GWAS-linked RNA-seq and eQTL figure

From the repository root:

```bash
conda activate simitall
Rscript analysis/paper_fig/fig6_rnaseq_eqtl_simulation.R
```

Optional arguments:

```bash
Rscript analysis/paper_fig/fig6_rnaseq_eqtl_simulation.R \
  --out_dir analysis/results/rnaseq_eqtl \
  --seed 42
```

The script writes:

- Simulated GWAS, RNA-seq, ASE, and eQTL result tables
- A six-panel PDF and PNG figure
- A one-row result summary TSV
- One CSV source-data table for every figure panel

The default runner uses the lightweight genome bundled with the package. For a
paper-scale analysis, replace `demo_genome` and the simulation dimensions in
the script with the selected reference, annotation, cohort size, and replicate
design.

## Donor-aware single-cell RNA-seq and cell-type eQTL figure

Run the complete GWAS donor to single-cell workflow with:

```bash
conda activate simitall
Rscript analysis/paper_fig/fig7_scrnaseq_celltype_eqtl.R
```

Optional arguments:

```bash
Rscript analysis/paper_fig/fig7_scrnaseq_celltype_eqtl.R \
  --out_dir analysis/results/scrnaseq_celltype_eqtl \
  --seed 52 \
  --backend native
```

Use `--backend auto` to select Splatter when installed and otherwise use the
native simulator. The runner writes sparse counts, cell/gene/truth metadata,
donor-cell-type pseudobulk matrices, cell-type eQTL results, a six-panel PDF
and PNG, and six panel-level source-data CSV files.

## Annotation-aware ChIP-seq figure

Run the TF binding, matched-input, and differential-binding example with:

```bash
conda activate simitall
Rscript analysis/paper_fig/fig8_chipseq_simulation.R --backend native
```

Optional arguments:

```bash
Rscript analysis/paper_fig/fig8_chipseq_simulation.R \
  --out_dir analysis/results/chipseq \
  --seed 62 \
  --backend auto
```

The runner writes regulatory target annotations, truth peaks, ChIP and input
FASTQ files, read positions, peak counts, differential-binding truth, library
QC, a six-panel PDF and PNG, and six panel-level source-data CSV files. Auto
uses ChIPsim for the bundled small reference when installed; native is the
memory-efficient option for chromosome-scale references.
