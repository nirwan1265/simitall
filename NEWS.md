# simitall 0.1.0

* Added one-command, profile-based bootstrap installers for Apple Silicon and
  Intel macOS, Linux/WSL2, and native Windows PowerShell. Profiles install
  minimal, population, omics, sequencing, or full dependency stacks.
* Added `check_simitall_dependencies()` and
  `install_simitall_dependencies()` for dependency audits and explicit
  profile installation from R. Native Windows sequencing/full installations
  delegate to WSL2 for Bioconda-only tools.
* Added annotation-aware TF, histone-mark, and nucleosome ChIP-seq simulation
  with GFF3/BED targets, matched input controls, biological and technical
  replicates, differential binding, PCR duplicates, sequencing errors, FASTQ,
  BED, peak-count, QC, and truth outputs.
* Added native and optional ChIPsim read-density backends plus a six-panel
  ChIP-seq publication figure and reproducible example runner.
* Added standalone bulk and single-cell RNA-seq experiment simulators that do
  not require GWAS genotypes, with explicit technical, biological, mixed, and
  custom replicate designs.
* Added separate biological- and library-level effects, standalone condition
  and batch truth, donor-level and library-level pseudobulk, and reusable
  single-cell publication plotting for one-sample experiments.
* Added donor-aware single-cell RNA-seq simulation from GWAS genotypes with
  shared and cell-type-specific eQTLs, marker programs, pseudotime, condition,
  batch, donor, latent, dropout, ambient-RNA, and doublet effects.
* Added native and optional Splatter technical backends, sparse 10x-style
  output, optional `SingleCellExperiment` export, donor-cell-type pseudobulk,
  cell-type eQTL benchmarking, and six-panel publication figures.

- Established the project as an installable R package.
- Added R-first genome, annotation, read, assembly, GWAS, phenotype, breeding,
  and SimuPOP interfaces.
- Added GWAS-linked RNA-seq simulation with cis/trans eQTLs, condition and
  batch effects, genotype-by-condition interactions, latent confounders,
  allele-specific expression, optional Rsubread FASTQ generation, and eQTL
  recovery benchmarks.
- Added RNA-seq/eQTL result summaries, six-panel publication figures, exported
  panel source data, and a reproducible end-to-end paper figure runner.
- Added lightweight package examples and on-demand reference download helpers.
