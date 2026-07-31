# simitall

`simitall` (SIMulate IT ALL) is an R-first toolkit for building coordinated,
truth-aware genomics simulations. It connects genome and annotation generation,
DNA read simulation, hybrid assembly evaluation, GWAS cohorts, phenotype
architectures, breeding designs, and SimuPOP mating schemes through one R API.

The goal is not to replace mature simulators. `simitall` provides a reproducible
R layer that makes those tools work together, keeps parameters in one analysis,
and emits compatible truth files for benchmarking.

## What it can simulate

- Random genomes or reference-derived genomes with tandem and motif repeats
- Uniform, Poisson, or fixed repeat spacing with optional GC matching
- Haploid through polyploid genome copies with SNP and indel divergence
- Genes, CDS features, operons, promoters, TSSs, terminators, rRNA/tRNA
  clusters, riboswitches, CRISPR arrays, plasmids, and regulatory elements
- Illumina reads with ART
- PacBio CLR and HiFi/CCS reads with PBSIM/PBSIM3
- Oxford Nanopore reads with Badread
- Hybrid assembly grids with Unicycler and evaluation with QUAST
- GWAS cohorts with LD blocks, recombination maps, subpopulations, and
  quantitative or binary phenotypes
- Advanced phenotype architectures through `simplePHENOTYPES`
- F1, F2, backcross, selfing, RIL, NIL, doubled-haploid, NAM, and MAGIC
  populations
- SimuPOP mating schemes through `reticulate`, with VCF-compatible genotype
  output and sample metadata

Truth outputs include FASTA, VCF, GFF3, BED, TSV, and JSON files. Functions
write results to user-selected output paths; generated analyses are not stored
inside the package source tree.

## Installation on Apple Silicon

The R package is small. Large reference genomes and command-line tools are
installed or downloaded separately.

```bash
cd /Users/nirwantandukar/Documents/Github/simitall

conda env create -f environment-macos.yml
conda activate simitall

R -q -e 'install.packages("remotes")'
R -q -e 'remotes::install_local(".", dependencies = TRUE)'
```

QUAST and SimuPOP are installed with `pip` in the supplied environment because
native Conda builds are not consistently available for `osx-arm64`.

For package development:

```r
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::document()
devtools::test()
devtools::check()
```

## Quick start

Load the package and locate its bundled, lightweight example files:

```r
library(simitall)

demo_genome <- system.file(
  "extdata", "examples", "demo_genome.fa",
  package = "simitall"
)
demo_panel <- system.file(
  "extdata", "panels", "demo_panel.fa",
  package = "simitall"
)
```

Generate a synthetic genome with controlled repeats:

```r
dir.create("results", showWarnings = FALSE)

simulate_genome(
  out_fa = "results/synthetic_genome.fa",
  random_genome = TRUE,
  random_length = 100000,
  random_gc = 0.51,
  mode = "both",
  n_events = 10,
  seg_len = 500,
  copies = 3,
  spacing_distribution = "poisson",
  spacing_mean = 5000,
  ploidy = 2,
  snp_rate = 0.001,
  indel_rate = 0.0001,
  seed = 1
)
```

The genome simulator can also modify an existing reference:

```r
simulate_genome(
  in_fa = demo_genome,
  out_fa = "results/reference_with_repeats.fa",
  n_events = 3,
  seg_len = 100,
  copies = 4,
  spacing_distribution = "fixed",
  fixed_spacing = 1000,
  gc_preserve = TRUE
)
```

Generate annotations:

```r
simulate_annotations(
  in_fa = "results/synthetic_genome.fa",
  out_prefix = "results/synthetic_genome",
  n_genes = 80,
  rrna_clusters = 1,
  trna_per_cluster = 8,
  riboswitch_count = 4,
  crispr_count = 1,
  plasmid_count = 1,
  regulatory_default_human = TRUE,
  seed = 1
)
```

## Sequencing and assembly

All external programs are called from R:

```r
sim_illumina_art(
  ref_fa = "results/synthetic_genome.fa",
  outprefix = "results/reads/illumina_",
  cov = 30
)

sim_pacbio(
  ref_fa = "results/synthetic_genome.fa",
  outdir = "results/reads/pacbio_hifi",
  cov = 20,
  type = "HIFI"
)

sim_nanopore(
  ref_fa = "results/synthetic_genome.fa",
  outdir = "results/reads/nanopore",
  cov = 20
)
```

Run a small coverage grid, hybrid assemblies, and QUAST:

```r
run_grid_both(
  simref_fa = "results/synthetic_genome.fa",
  tag = "demo",
  ill_covs = c(20, 40),
  pb_covs = c(10, 20),
  output_dir = "results/reads"
)

run_unicycler_grid(
  tag = "demo",
  reads_dir = "results/reads",
  output_dir = "results/assemblies",
  ill_covs = c(20, 40),
  pb_covs = c(10, 20)
)

run_quast_grid(
  tag = "demo",
  ref_fa = "results/synthetic_genome.fa",
  assemblies_dir = "results/assemblies",
  output_dir = "results/evaluation",
  ill_covs = c(20, 40),
  pb_covs = c(10, 20)
)

summarize_quast(
  quast_root = "results/evaluation/demo",
  out_csv = "results/quast_summary.csv"
)
```

## GWAS and phenotypes

```r
simulate_gwas_cohort(
  genome_fa = demo_genome,
  out_prefix = "results/gwas/demo",
  n_samples = 200,
  snp_rate = 0.01,
  n_subpops = 2,
  ld_blocks = TRUE,
  phenotype = "quantitative",
  n_causal = 10,
  seed = 1
)
```

Advanced phenotype simulation is optional and requires
`simplePHENOTYPES`:

```r
simulate_phenotypes(
  geno_file = "results/gwas/demo.vcf",
  out_prefix = "results/phenotypes/demo",
  h2 = 0.5,
  n_add_qtn = 20,
  n_traits = 2
)
```

## Breeding and SimuPOP

Generate a founder panel or use the bundled demo:

```r
generate_random_haplotype_panel(
  out_fa = "results/founders.fa",
  n_haplotypes = 8,
  length = 50000,
  snp_rate = 0.005,
  seed = 1
)

simulate_breeding(
  haplotype_fa = demo_panel,
  out_prefix = "results/breeding/ril",
  scheme = "RIL",
  n_offspring = 100,
  generations = 6,
  interference_shape = 2,
  genotyping_error = 0.005,
  missing_rate = 0.01,
  seed = 1
)
```

Run a SimuPOP preset:

```r
simupop_api(
  out_prefix = "results/simupop/f2",
  config = list(
    population = list(size = 100, ploidy = 2, loci = 200),
    preset = "F2",
    generations = 2,
    export_vcf = TRUE
  )
)
```

`simupop_api()` also supports configured random, monogamous, polygamous,
selfing, hermaphroditic, clonal, conditional, heterogeneous, pedigree, and
controlled mating schemes. A `python_hook` field can define specialized
SimuPOP parent choosers or frequency trajectories before a run.

## Example data

Small synthetic files are bundled under `inst/extdata`:

- `examples/demo_genome.fa`: tiny genome for smoke tests
- `panels/demo_panel.fa`: aligned synthetic founder haplotypes
- `markers/markers_curated.tsv`: antibiotic-marker accession manifest
- `regulatory/regulatory_panel_human.tsv`: synthetic human-inspired
  enhancer/silencer placement panel
- `config/simupop_example.json`: example SimuPOP configuration

Download larger references only when needed:

```r
fetch_ref_files("references")
fetch_marker_panel("references/antibiotic_markers.fa")
```

`fetch_ref_files()` retrieves E. coli K-12 MG1655, human GRCh38 chromosome 1,
maize B73 chromosome 1, and Arabidopsis chromosome 1 from NCBI. Users can pass
their own FASTA, GFF3, recombination maps, regulatory panels, marker panels,
and aligned founder haplotypes to the corresponding simulators.

## Package layout

```text
simitall/
├── R/                  package functions and internal simulation engines
├── man/                generated R help pages
├── inst/extdata/       lightweight manifests and demo files
├── tests/testthat/     automated tests
├── DESCRIPTION         package metadata and dependencies
├── NAMESPACE           generated exports
└── README.md           overview and examples
```

## Reproducibility and attribution

Set `seed` in each simulator and record tool versions for publication. The
package orchestrates external software; publications should cite both
`simitall` and the underlying tools used in an analysis, including ART, PBSIM
or PBSIM3, Badread, Unicycler, QUAST, SimuPOP, and simplePHENOTYPES as
applicable.

## Status

`simitall` is under active development. Simulated annotations and regulatory
elements are benchmarking truth, not biological claims. Large-genome and
high-coverage runs can require substantial memory, storage, and compute time.

## License

MIT
