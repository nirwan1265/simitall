# SimItAll

`SimItAll` (SIMulate IT ALL) is an R framework for building
truth-aware genomics simulations, analyzing their outputs, and benchmarking
results against the truth used to generate them. It connects genome and
annotation generation, DNA read simulation and assembly, breeding populations,
GWAS and genomic selection, bulk and single-cell RNA-seq/eQTLs, ChIP-seq, and
reproducible reports through one R API. In short, simulate all seq analysis.

SimItAll keeps inputs, parameters, random seeds, provenance, and truth outputs
together so that an entire study can be rerun and checked.

SimItAll is the **Simitall** head of the Biodes suite. Biodes
([BioDesignGPT](https://github.com/nirwan1265/BioDesignGPT)) reads the
experiment in plain language and writes a design folder (`study_spec.json`,
`field_layout.csv`); SimItAll creates data for every unit of that layout:

```r
library(simitall)
simitall_routes("~/biodes_runs/my_design")        # what will be simulated, with which values
simitall_from_biodes("~/biodes_runs/my_design",   # outcome + genetics + omics, with sim_truth.json
  params = list(suite_genetics = list(heritability = 0.6, phenotype_model = "gxe")))
```

What runs is looked up from the spec's fields in
[`inst/suite/simulators.tsv`](inst/suite/simulators.tsv) (outcome type,
genetics from `assays`/`population.mating`, RNA-seq, eQTL, single-cell,
ChIP-seq, ATAC-seq); parameters and their defaults are in
[`inst/suite/simulator_parameters.tsv`](inst/suite/simulator_parameters.tsv).
See [docs/BIODES_BRIDGE.md](docs/BIODES_BRIDGE.md).

**Start here:** [Install](#install-simitall) | [Detailed tutorial](docs/TUTORIAL.md) | [Validation and benchmarks](docs/VALIDATION.md)

## Install simitall

SimItAll keeps R dependencies, Python modules, and external scientific
software outside this repository. Select a profile that matches the work you
need; `full` installs the supported union for the operating system.

```bash
git clone https://github.com/nirwan1265/simitall.git
cd simitall

# macOS or Linux
bash install_simitall_macos.sh population

# Windows PowerShell
powershell -ExecutionPolicy Bypass -File install_simitall_windows.ps1 -Profile population
```

| Profile | Includes |
| --- | --- |
| `minimal` | R package and basic simulation dependencies |
| `population` | SimuPOP, breeding, phenotype simulation, GWAS/genomic-selection R packages, TASSEL, and GEMMA where the platform supports it |
| `omics` | Bulk RNA-seq, single-cell, ChIP-seq, and supporting R/Bioconductor packages |
| `sequencing` | ART, PBSIM/PBSIM3, Badread, Unicycler, QUAST, and helper tools |
| `full` | Every supported profile component for the current operating system |

TASSEL is installed from Bioconda in the `population` profile. GEMMA is added
on Linux and Intel macOS; native Windows and Apple Silicon users retain the
portable R mixed-model workflow and receive a clear GEMMA fallback message.
Use `check_simitall_dependencies("population")` after activation to inspect
what is available.

## What SimItAll Does

| `SIMULATE` | `ANALYZE` |
| --- | --- |
| Create genomes, annotations, reads, founder panels, pedigrees, populations, phenotypes, and multi-omics truth, directly or for a Biodes design. | Run and benchmark assembly, GWAS, eQTL, and genomic-selection analyses against known simulated truth. |

## Simulate, Analyze, and Validate

SimItAll supports the complete study loop: create data with known truth, run
an analysis, and check the result against the truth and design requirements.
It can also analyze compatible user-supplied inputs, but real-data conclusions
remain the user's scientific responsibility.

| Study area | Simulations | Analyses and outputs |
| --- | --- | --- |
| **Genomes and populations** | Genomes, annotations, founder panels, pedigrees, designed crosses, NAM/MAGIC/RIL/NIL/DH populations, ancestry, recombination, and structural truth | Segregation, ancestry, breakpoint, LD, relatedness, and breeding-composition summaries |
| **Traits and association** | Quantitative, binary, count, ordinal, family-effect, and GxE traits with causal truth | GWAS, QTL/causal-truth recovery, calibration, kinship/PC adjustment, and genomic prediction/selection |
| **Bulk and single-cell transcriptomics** | Bulk RNA-seq, eQTLs, ASE, batches, cell types, pseudotime, dropout, ambient RNA, doublets, and cell-type eQTLs | Differential expression, donor-aware pseudobulk, eQTL benchmarking, clustering, PCA, volcano, and QC figures |
| **Chromatin** | ChIP-seq and ATAC-seq reads, peaks, matched controls, conditions, and differential-accessibility truth | Peak calling, differential peaks, enrichment/QC, and precision/recall/F1 against retained truth |
| **Gene and functional annotation** | Compact assembly-aware gene-coordinate catalogues and annotation-compatible synthetic features | Interval-to-gene overlap, local GFF3/GTF parsing, Gene Ontology GAF import, evidence-preserving gene-to-GO joins, and GO-term summaries |
| **Sequencing and assembly** | Illumina short reads, PacBio HiFi/CLR reads, Nanopore reads, hybrid read sets, coverage grids, and synthetic references | Assembly with Unicycler, QUAST comparison, coverage/read-length QC, and truth-aware assembly metrics |

### Gene Coordinates and GO Annotations

SimItAll includes compact genome-wide, gene-only coordinate catalogues for
Arabidopsis, barley, rice, wheat, maize, and fly. They retain chromosome,
start/end, strand, stable gene ID, source ID, gene name, and biotype while
avoiding the size of full GFF3/GTF files. Small GFF3 and Arabidopsis GO examples
ship with the package for reproducible demonstrations.

For functional annotation, load a dated local Gene Ontology Annotation File
(GAF) with `read_go_gaf()`, join it to coordinates with `join_gene_go()`, and
retain GAF evidence codes in the result. The repository includes the supplied
full GAF downloads under `resources/go_annotations/`, but they are deliberately
not installed with the R package.

For an organism not represented in the repository, download its current
organism-specific GAF and matching ontology from the
[Gene Ontology annotation downloads](https://geneontology.org/docs/download-go-annotations/downloads/), then supply a reviewed gene-ID mapping when its GFF3/GTF identifiers differ from the GAF identifiers. SimItAll does not infer GO assignments or silently guess identifier conversions.

```r
genes <- read_gene_coordinates("Arabidopsis_thaliana")
go <- read_go_gaf("resources/go_annotations/ARATH-mod.gaf.gz")
gene_go <- join_gene_go(genes, go)
```

### Real-Data Import Catalogue

SimItAll also maintains a machine-readable catalogue of public, real-data
examples for bulk RNA-seq, ChIP-seq, GWAS, and single-cell RNA-seq. It records
download sources, expected inputs, companion gene-annotation requirements, and
import/QC checks without bundling large public datasets in the package.

Use `inst/extdata/real_data_catalog.json` to select a documented
example and keep downloaded files outside Git. Each local download folder
should include a completed `manifest.json` based on
`data-raw/real_data_manifests/manifest.template.json`; it records provenance,
reference assembly, annotation release, gene-ID namespace, and sample-ID
mappings.

A JSON manifest points to annotation resources; it is not a substitute for a
GFF3/GTF. Before reporting gene, peak, or locus annotations, SimItAll must
confirm compatible assembly, chromosome labels, coordinates, and gene IDs. See
[real-data import examples](docs/REAL_DATA.md) for the initial public datasets
and download instructions.

### Sequencing Platforms

Use `simitall_sequencing_profiles()` to see the validated platform labels:

- **Illumina:** `MiSeq`, `MiSeq i100`, `NextSeq 1000/2000`, and `NovaSeq X`.
- **PacBio:** `Sequel IIe HiFi`, `Revio HiFi`, and `Vega HiFi`.
- **Oxford Nanopore:** Nanopore simulation through Badread, with configurable
  error/length models for MinION-, GridION-, or PromethION-style study designs.

Platform labels select transparent simulator settings and defaults. They do not
claim to reproduce every current instrument chemistry; empirical error models
should be supplied when instrument-specific accuracy matters.

### Validation and Benchmarking

Validation lives in **[docs/VALIDATION.md](docs/VALIDATION.md)**, not on this
project overview page. It contains the reproducible figures, fixed-seed
benchmarks, source-data expectations, and interpretation boundaries.

Each supported workflow is also described by a four-part validation contract:

1. Simulation truth retained by the run.
2. Required result files.
3. Required figures.
4. Pass criteria for the programmed mechanism.

See [the validation-contract README](inst/agent/knowledge/validation_contracts/README.md)
for the contract rules and [docs/VALIDATION.md](docs/VALIDATION.md) for the
human-readable demonstrations.

## Capability Details

- Random genomes or reference-derived genomes with tandem and motif repeats
- Uniform, Poisson, or fixed repeat spacing with optional GC matching
- Haploid through polyploid genome copies with SNP and indel divergence
- Genes, CDS features, operons, promoters, TSSs, terminators, rRNA/tRNA
  clusters, riboswitches, CRISPR arrays, plasmids, and regulatory elements
- GWAS cohorts with LD blocks, recombination maps, subpopulations, and
  quantitative or binary phenotypes
- Mixed-model GWAS with kinship, principal components, fixed effects,
  Manhattan/QQ plots, and causal-variant recovery benchmarks
- GBLUP, RR-BLUP, Bayesian, and random-forest genomic prediction with
  leakage-aware cross-validation, parent selection, and crossing plans
- Advanced quantitative and binary phenotype architectures through
  `simplePHENOTYPES`
- Ordinal and count traits plus shared breeding-family effects through
  `simstudy`
- Direct, maternal, and paternal pedigree effects through `pedtricks`
- F1, F2, backcross, selfing, RIL, NIL, doubled-haploid, NAM, and MAGIC
  populations
- SimuPOP mating schemes through `reticulate`, with VCF-compatible genotype
  output and sample metadata
- Illumina reads with ART
- PacBio CLR and HiFi/CCS reads with PBSIM/PBSIM3
- Oxford Nanopore reads with Badread
- Hybrid assembly grids with Unicycler and evaluation with QUAST
- Standalone bulk RNA-seq experiments with technical, biological, or mixed
  replicates, even when no genotype data are available
- RNA-seq from the same GWAS individuals, with cis/trans eQTLs, condition and
  batch effects, genotype-by-condition interactions, latent confounders, and
  allele-specific expression
- eQTL association testing with causal-truth precision and recall summaries
- Standalone single-cell RNA-seq experiments with technical, biological, or
  mixed replicates
- Single-cell RNA-seq from GWAS donors with cell types, marker programs,
  pseudotime, cell-type-specific eQTLs, dropout, ambient RNA, and doublets
- Donor-by-cell-type pseudobulk and cell-type eQTL benchmarking without
  treating cells as independent biological replicates
- Annotation-aware TF, histone-mark, and nucleosome ChIP-seq with matched
  inputs, differential binding, FASTQ reads, peak truth, and publication plots

Truth outputs include FASTA, VCF, GFF3, BED, TSV, and JSON files. Functions
write results to user-selected output paths; generated analyses are not stored
inside the package source tree.

## Detailed tutorial

The full, step-by-step simulation guide is intentionally kept off this front
page so the project overview stays easy to scan. It includes runnable R
examples for genome and annotation generation, GWAS and traits, multi-chromosome
breeding, sequencing and assembly, bulk RNA-seq/eQTLs, single-cell RNA-seq, and
ChIP-seq.

**[Open the detailed simulation tutorial](docs/TUTORIAL.md)**

For generated figures, fixed-seed demonstrations, validation summaries, and
reproduction commands, see **[Results and Validation](docs/VALIDATION.md)**.

## Package layout

```text
simitall/
├── analysis/           paper-figure runners and generated output folders
├── docs/               detailed tutorial and validation documentation
├── R/                  package functions and internal simulation engines
├── man/                generated R help pages
├── inst/extdata/       lightweight manifests and demo files
├── tests/testthat/     automated tests
├── DESCRIPTION         package metadata and dependencies
├── NAMESPACE           generated exports
├── README.md           project overview and installation
├── docs/TUTORIAL.md    detailed simulation workflows
└── docs/VALIDATION.md  reproducible figures and validation notes
```

## Reproducibility and attribution

Set `seed` in each simulator and record tool versions for publication. The
package orchestrates external software; publications should cite both
`simitall` and the underlying tools used in an analysis, including ART, PBSIM
or PBSIM3, Badread, Unicycler, QUAST, SimuPOP, and simplePHENOTYPES as
applicable. Analyses that generate RNA-seq FASTQ files should also cite
`Rsubread`/`simReads`. Single-cell analyses using `backend = "splatter"`
should cite Splatter in addition to `simitall`. ChIP-seq analyses using
`backend = "chipsim"` should cite ChIPsim.

## Future work

Planned, roughly in order (gains are from the 2000-question check in
`docs/SIMULATION_COVERAGE.md`):

**Better templates (real data to learn from)**

- Plant WGBS template: plants now use the native CG/CHG/CHH engine. Two to four
  Arabidopsis WGBS samples with all three contexts (e.g. GEO GSE43857, per-site
  methylated/total read counts) would let plant methylation learn from real
  data, as mouse does. A maize or sorghum methylome would be better still.
- Human WGBS template: human questions now borrow the mouse ENCODE template.
- Rice, potato and soybean miRNA libraries with sample labels (PmiREN) to add
  to the plant small RNA-seq templates.

**New simulators**

- MeDIP (18 questions): enrichment counts per window, like ChIP IP/input.
- m6A / MeRIP-seq (+40): IP and input counts per peak.
- Ribo-seq (+26): RNA and footprint counts per gene with translational-efficiency truth.
- Hi-C (+21): contact matrices with planted loop and TAD changes.
- Degradome (+20), scATAC-seq (+16), MNase-seq (+10), spatial (+9),
  metabolomics/lipidomics (+10), RIP/CLIP (+5).
- Genetics with no lines (24 of the second 1000 questions): BSA / pool-seq
  bulks (allele frequencies by phenotype bulk) and population resequencing
  (populations with drift and selection, evolve-and-resequence). Today the rest
  of such a design is simulated and genetics is listed as `needs_genotype_factor`.
- Survival and compositional outcomes.

**Design questions**

- "How many per group?" questions have no replicate count; send them to
  `simitall_power_sweep()` over a range of replicates instead of assuming 3.
- Questions whose reviewed spec has no factors or no level count stop before
  simulation; Biodes should ask the user for these.

**Read level**

- Read-level engines (FASTQ/BAM) for the count-level simulators that do not
  have one yet (methylation, small RNA, microbiome).

## Status

`simitall` is under active development. Simulated annotations and regulatory
elements are benchmarking truth, not biological claims. Large-genome and
high-coverage runs can require substantial memory, storage, and compute time.

## License

MIT
