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

## Maize NAM chr10 drought GxE, RNA-seq, GWAS, and genomic selection

This portfolio-scale runner starts from a synthetic eight-founder NAM panel on
one chromosome, retains breeding ancestry and breakpoints, introduces a focal
GxE locus under irrigated and drought environments, and uses the same NAM
individuals for bulk RNA-seq/eQTL truth, drought GWAS, and one GBLUP selection
cycle.

```bash
conda activate simitall
Rscript analysis/paper_fig/fig12_maize_nam_chr10_gxe_multiomics.R \
  --out_dir analysis/results/maize_nam_chr10_gxe_multiomics --seed 1201
```

The default uses a bundled synthetic chr10 panel and creates an explicit
`SIM_DROUGHT_RESPONSE_CHR10` feature plus background genes. To anchor the
focal feature to a real, compatible chr10 annotation, provide a local GFF3:

```bash
Rscript analysis/paper_fig/fig12_maize_nam_chr10_gxe_multiomics.R \
  --gff3 /path/to/maize_chr10.gff3 \
  --out_dir analysis/results/maize_nam_real_annotation --seed 1201
```

Outputs include founder/breakpoint truth, drought and irrigated phenotypes,
trait-QTL truth, kinship-aware GWAS and recovery metrics, bulk RNA-seq counts,
eQTL truth, selected parents/crosses, and a `figures/` directory containing
individual phenotype, Manhattan, QQ, RNA-seq PCA, library-size, volcano, and
selection plots plus a combined overview. This is an integrated synthetic demonstration, not evidence that the
default focal feature is a real drought gene or that its programmed effects
match any named maize population.

## Human relatedness and GWAS-calibration benchmark

This runner creates a synthetic chr10 cohort with two ancestry groups and
sibling-pair relatedness, deliberately programs phenotype confounding without
any causal marker, then compares an unadjusted scan with PC-plus-family/kinship-proxy adjustment.

```bash
conda activate simitall
Rscript analysis/paper_fig/fig13_human_relatedness_gwas_calibration.R \
  --out_dir analysis/results/human_relatedness_gwas_calibration --seed 1301
```

The validation summary records genomic inflation (`lambda`) for both scans and
checks that adjustment reduces inflation for this fixed synthetic design. The
outputs and figure are calibration teaching material only: they are not a
human disease simulation, clinical model, or a substitute for a validated
mixed-model association workflow on a real cohort.
