# simitall (development)

- Factor levels follow the spec's `factors.levels` (reference first), then
  natural order for levels it does not list, as in Bioana. Before, levels were
  always in natural order, so `trt` became the reference before `untrt` and the
  planted effects had the opposite sign of what the spec meant. The outcome,
  genetics, assays, truth files and the outcome figure all use the new order.
- `simitall_power_sweep()` and `simitall_from_question()` take
  `allow_fallback` (default `FALSE`, as in Bioana) and warn when Bioana
  refuses an analysis because its citable package is not installed.
- All analysis moved to Bioana (github.com/nirwan1265/bioana): GWAS, genomic
  prediction and selection, peak calling and differential peaks, single-cell
  clustering and tests, eQTL benchmarks, RNA-seq results summary, GO
  annotation (with the gene-coordinate catalogues), assembly evaluation and
  `bioana_fit_biodes()`. Simitall only simulates; `simitall_power_sweep()`
  and `simitall_from_question()` call `bioana::bioana_analyze()`, and
  `biopaint_biodes()` reads Bioana's output.
- `data_manifest.json`: every run lists its data sets (data type, files, truth
  in Bioana's format under `truth/`) for Bioana (github.com/nirwan1265/bioana).
- Biodes 0.5.1: routing reads `assay_kinds` (DNase/FAIRE -> ATAC, DAP-seq/SELEX
  -> ChIP; amplicon genotyping, repertoires, CRISPR screens, reporter assays
  planned), breeding notation sets backcross and selfing generations, and
  survival outcomes simulate (log-normal time with censoring, `event` column).
- Backcross / NIL lines (rare donor allele) no longer stop the genetics layer:
  QTN are drawn among the commonest polymorphic markers, with a note.
- `simulate_merip()`: m6A MeRIP-seq as IP and input counts per window (peaks
  and background), with known peaks, changed peaks and expression changes that
  move IP and input together. Biodes designs with `m6a_seq` / `merip` simulate it.
- "How many per group?": when a question has no replicate count,
  `simitall_from_question()` chooses `design.reps` with `simitall_power_sweep()`
  (reads `power.reps_range` and `power.solve_for`). The sweep now also reports,
  for count assays, the share of genes with a true 2-fold change that are found.
- Bioana's count test uses empirical-Bayes moderated t-tests (as limma eBayes),
  so designs with few replicates are tested fairly.
- A design whose genetics assay has no genotype factor (BSA bulks, populations)
  no longer stops: the rest is simulated and genetics is listed as
  `needs_genotype_factor` in `sim_truth.json$not_simulated`.
- `simulate_methylation()`: bisulfite-sequencing methylation (WGBS, RRBS,
  EM-seq) as methylated and covering read counts per cytosine, with known DMRs
  (`dmr_truth.tsv`) and per-site truth. Three engines: a real template of the
  same kingdom (mouse ENCODE WGBS), a native engine with CG/CHG/CHH chromatin
  states for plants (`inst/suite/methylation_presets.json`), and beta values
  for 450k/EPIC arrays (TCGA template). Biodes designs with methylation assays
  now simulate it; `simulator_parameters.tsv` has an `ask` column whose
  questions `simitall_routes()` prints.
- Real-data templates for count simulation: `simitall_templates()`,
  `simitall_template()`, `simulate_counts_from_template()` (binomial thinning
  of real samples, or NB learned from them). Small RNA-seq (human and mouse
  miRNA atlas) and 16S/ITS microbiome (plant root, human HMP) now simulate
  from them for Biodes designs.
- Simitall now creates data for a Biodes design from the spec's fields:
  `inst/suite/simulators.tsv` routes `primary_outcome.type` (continuous,
  binary, count, ordinal), genetics (`assays` gwas/qtl_mapping/eqtl,
  `population.mating`) and omics assays (RNA-seq, eQTL, single-cell, ChIP-seq,
  ATAC-seq) to simulators, with parameters in
  `inst/suite/simulator_parameters.tsv`. New: `simitall_routes()`,
  `simitall_simulator_parameters()`, `params` in `simitall_from_biodes()`.
- Removed the prompt readers (`simitall_ask()`, `launch_simitall_agent()`,
  `search_simitall_knowledge()`, `simitall_agent_tools()`, the keyword spec,
  lexicon, prompt model and capability catalog): Biodes reads the user's
  question. They and the Shiny app are kept in `old_model/` (not built or
  loaded); see `old_model/README.md` to run them.

# simitall 0.1.0

- Added a retrieval-grounded, read-only agent interface through
  `search_simitall_knowledge()`, `simitall_ask()`, and
  `simitall_agent_tools()`. Added an optional local Shiny chat interface that
  uses a locally stored OpenAI API key and never executes simulations.

- Added mixed-model GWAS analysis, Manhattan/QQ plotting, and causal-truth
  benchmarking through `rrBLUP`.
- Added GBLUP, RR-BLUP, Bayesian, and random-forest genomic prediction,
  cross-validation, genomic kinship, diversity-aware parent selection,
  crossing-plan design, and one-round genomic-selection orchestration.
- GWAS cohort simulation now writes causal-variant truth and per-sample true
  breeding values.

* Added package-backed breeding phenotype wrappers: `simstudy` provides
  ordinal, Poisson, negative-binomial, and shared family effects, while
  `pedtricks` provides direct, maternal, and paternal pedigree effects with
  component-level truth output.
* Added multi-chromosome breeding with `founder|chromosome` FASTA panels,
  chromosome-specific recombination maps, independent assortment,
  chromosome-aware VCF output, founder-ancestry tracts, and breakpoint truth.
  Added a bundled three-chromosome founder panel and generator controls for
  chromosome counts and lengths.
* Expanded the breeding tutorial with runnable F1, F2, F2:3-style,
  F2-derived S3, backcross, RIL-SSD, RIL-sib, DH, NIL, NAM, MAGIC, and SimuPOP
  examples. Corrected F2 and MAGIC generation expansion and made empty VCF
  truth sets safe for monomorphic final populations.
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
