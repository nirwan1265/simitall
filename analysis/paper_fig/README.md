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
