# Simitall

<p align="center">
  <img src="docs/bio_suite.svg" alt="BioCore at the center, shared by four heads: Biodes designs the experiment, Simitall simulates the data, Bioana analyzes the results, Biopaint makes the figures" width="860">
</p>

**Simulate the experiment before you run it, with the answer key.** Simitall is the simulation head of a small
suite of free, offline R tools for plant and biological research. It reads a design written by
[Biodes](https://github.com/nirwan1265/BioDesignGPT) and simulates realistic data for every plot, plant, seed or
sample in that design. The effects in the data are planted, so their true sizes are known. Everything planted,
and every assumption made, is written to `sim_truth.json`. [Bioana](https://github.com/nirwan1265/Bioana) then
analyses the data and is scored against that truth.

What this is good for: checking that an analysis finds what is really there before any real data exist, seeing
how often a design detects an effect of a given size (power by simulation), and teaching or testing statistical
methods on data whose answer is known.

## The suite

| Part | Job |
|---|---|
| **[BioCore](https://github.com/nirwan1265/BioCore)** | Shared rules and StudySpec; Shiny app that runs design → simulate → analyze from one text box; worked examples |
| **[Biodes](https://github.com/nirwan1265/BioDesignGPT)** | Design: unit audit, power, randomized layout, `study_spec.json` |
| **Simitall** (this repo) | Simulate data under that design with planted, known truth |
| **[Bioana](https://github.com/nirwan1265/Bioana)** | Analyse with the model that matches the design and score it against Simitall's truth |
| **Biopaint** | Figures (in progress) |

Worked examples of the whole flow live in [BioCore](https://github.com/nirwan1265/BioCore).

## What Simitall reads

Only the three files Biodes `write_design_package()` writes:

| File | Used for |
|---|---|
| `study_spec.json` | factors and levels, the analysis model (`analysis.model`, `analysis.terms`), goals, outcome type and scale, power values, survival, covariate, within-subject and genetics settings |
| `field_layout.csv` | the randomized units: block, plot, subplot, treatment levels, positions |
| `sample_metadata_template.tsv` | one row per measured unit; Simitall fills its outcome columns |

Nothing is assumed silently. Where the spec leaves a value open (effect size, residual SD, heritability,
number of genes), Simitall uses the default in [`inst/suite/simulator_parameters.tsv`](inst/suite/simulator_parameters.tsv)
and records it under `assumptions` in `sim_truth.json`.

## What Simitall writes

Into the design folder (or `outdir`):

| File | Contents |
|---|---|
| `sample_metadata_template.tsv` | the same sample sheet with `primary_outcome` (and secondary traits, raters, time points, covariates) filled in |
| `counts.tsv` | feature × sample count matrix for bulk RNA-seq; other assays write `<assay>/unit_counts.tsv` |
| `genotypes.vcf` | genotypes of the design's lines, when the design has a genotype factor |
| `genetic_truth/` | `qtn_truth.tsv` (planted QTNs and their effects), `genetic_values.tsv` (true genetic value per line) and, for panels, `population_structure.tsv` |
| `truth/` | per-feature truth for omics goals (which genes carry the planted effect) |
| `sim_truth.json` | the answer key (below) |
| `data_manifest.json` | which file holds which data, and where its truth is; read by Bioana |

### `sim_truth.json`

| Key | What it holds |
|---|---|
| `columns`, `model`, `design`, `outcome_type`, `seed` | how the spec's model terms map to sample-sheet columns, the model simulated, the design, outcome type, seed |
| `goal_truth` | the true answer to the study's question: goal type (treatment effect, interaction, trend, regression, survival, ...), the terms, `true_effect`, and per-level values such as `difference_of_differences` |
| `treatment_effects` | planted effect of every level of every treatment factor |
| `interaction_effects` | planted interaction effects |
| `adjustment_effects` | planted effects of adjustment covariates (ANCOVA, sex, baseline) |
| `intercept`, `residual_sd`, `random_sd`, `block_effects` | the error structure: residual SD, random-effect SDs (e.g. whole plots in a split plot), realised block effects |
| `outcome_model` | the generating model for non-normal outcomes; for survival, per-group `event_free_survival_at`, `event_free_rmst`, `cif_event_at`, `median`, conditional and `marginal` versions, and `measure_notes` |
| `outcome_scale` | unit, baseline, residual SD and CV of the outcome (standardised when the spec gives none) |
| `regression_truth` | planted slopes for continuous predictors |
| `within_subject` | correlation structure for repeated measures (e.g. AR(1) and its rho) |
| `secondary_traits`, `agreement` | correlated secondary traits and planted rater agreement |
| `genetics` | heritability, phenotype model and paths into `genetic_truth/` |
| `assumptions` | every value Simitall chose because the spec left it open, in plain words |

These key names are a contract with Bioana and BioCore; they are not renamed without notice.

## What it can simulate

**Outcomes:** continuous, binary (yes/no), count, ordinal (scores), nominal (classes), and time-to-event with
censoring, competing events and frailty. Repeated measures with correlated time points (AR(1)), covariates
(regression, ANCOVA), secondary traits and multiple raters.

**Designs:** anything Biodes writes, including CRD, RCBD, alpha lattice, split
plots, factorials, multi-location and multi-season trials, repeated measures.

**Genetics**, when the design has a genotype factor:

| Path | When | How |
|---|---|---|
| Diversity panel | GWAS, GBS, WGS, capture, or a breeding study with no cross given | structured panel with subpopulations (built in, pure R) |
| Cross | RIL, DH, NAM, MAGIC, backcross/NIL, F2, selfing, sib mating; QTL mapping, BSA-seq, pool-seq | `simulate_breeding()` with the scheme in `inst/suite/mating_schemes.tsv` (pure R) |
| simuPOP | forward population simulations with your own mating schemes | `simupop_api()`, through Python |

Phenotype models: additive, additive + dominance, epistatic, genotype × environment, or all together, scaled to a
chosen heritability. The panel and cross paths need no Python. The simuPOP path needs `reticulate` and a Python
environment with simuPOP; `reticulate` is optional (Suggests), so Simitall installs and runs without Python,
including in webR.

**Omics assays:** bulk RNA-seq, eQTL, single-cell RNA-seq, small RNA, ChIP-seq, ATAC-seq, DNA methylation, m6A
and microbiome amplicons, each with planted features.

## Install

```r
# install.packages("remotes")
remotes::install_github("nirwan1265/simitall")
remotes::install_github("nirwan1265/BioDesignGPT")   # Biodes, to write designs
remotes::install_github("nirwan1265/Bioana")         # optional: analyse and score
```

Optional, for the simuPOP path only: `install.packages("reticulate")` and a conda environment named `simitall`
with simuPOP (`conda create -n simitall -c conda-forge simupop`). Install profiles for the external
simulators and sequencing tools are in the [detailed README](docs/README_DETAILED.md).

## Use

```r
library(simitall)

# a folder written by Biodes write_design_package()
simitall_from_biodes("my_design", seed = 2026)

truth <- jsonlite::read_json("my_design/sim_truth.json")
truth$goal_truth                                    # the true answer to the study's question

simitall_power_check("my_design", n_sim = 100)      # share of simulated experiments where Bioana finds it
```

Change what is planted with arguments (`effect`, `residual_sd`, `interaction_effect`, `n_genes`) or with
`params`, a list of overrides named by simulator, e.g. `params = list(suite_genetics = list(heritability = 0.3))`.

## More

- [Detailed README](docs/README_DETAILED.md): install profiles, every simulator, the R API, validation, benchmarks
- [How Simitall reads a Biodes design](docs/BIODES_BRIDGE.md) | [Coverage over published study designs](docs/SIMULATION_COVERAGE.md) | [Validation](docs/VALIDATION.md)

Simitall is part of the Biodes suite by Nirwan Tandukar (NC State University). MIT licence.
