# SimItAll: AI Agent Handoff

**Read this before making changes.** This repository is actively edited by the
project owner with more than one coding assistant. The working tree can be
intentionally dirty. Preserve existing work, inspect `git status` first, and
do not revert, reset, rename, or reformat unrelated changes.

Last handoff update: 2026-11-22.

## Project Purpose

SimItAll is the **Simitall head** of the Biodes suite: it creates data, with
known truth, for an experiment that Biodes has designed.

| Head | Repo | Job |
| --- | --- | --- |
| Biodes | github.com/nirwan1265/BioDesignGPT (package `studyforge`) | Reads the user's question, audits the unit, picks the design, writes `study_spec.json` + `field_layout.csv` |
| **Simitall** | this repo | Simulates data for every unit of that layout, writes `sim_truth.json` |
| Bioana | this repo (`bioana_fit_biodes()`) | Fits the design's model, checks recovery against the truth |
| Biopaint | this repo (`biopaint_biodes()`) | Figures |

**Simitall does not read prompts.** Biodes reads the user's words; Simitall
reads Biodes' StudySpec fields. The old prompt readers (keyword spec, legacy
router, lexicon reader, BERT tagger, capability catalog) and `simitall_ask()`
were removed on 2026-11-22 and saved in `old_model/` (not built, loaded or
tested; `old_model/README.md` says how to run them). Do not add a new prompt reader here: a missing field is
added in Biodes (BioDesignGPT `docs/SUITE_INTEGRATION.md`).

## Collaboration Rules

- Treat current working-tree changes as another collaborator's work unless the
  user explicitly asks to change them.
- Never use `git reset --hard`, `git checkout --`, or destructive cleanup.
- Do not invent biological results, package availability, input files, or GO
  assignments.
- No hard-coding: vocabulary, routes, parameters and defaults live in the
  TSV/JSON rule files under `inst/suite/`; R code holds the logic.
- Deterministic code computes every number. Every default used is written to
  `sim_truth.json$assumptions`.
- A request Simitall cannot simulate is reported (`planned` / `not_supported`
  in `simulators.tsv`, and `sim_truth.json$not_simulated`), never replaced by
  a different simulation.
- Commit as nirwan1265 (Nirwan Tandukar), no AI attribution lines.

## Main Architecture

```text
Biodes design folder (study_spec.json, field_layout.csv, sample_metadata_template.tsv)
  -> simitall_routes(dir)          which simulators run, with which parameter values
  -> simitall_from_biodes(dir)     outcome + genetics + assays for every unit
       outcome   primary_outcome.type -> continuous / binary / count / ordinal
       genetics  assays (gwas, qtl_mapping, eqtl, ...), population.mating, study_type
                 -> genotypes of the genotype factor's levels + genetic values
       assays    assays -> RNA-seq, eQTL, single-cell, ChIP-seq, ATAC-seq
  -> sim_truth.json, genotypes.vcf, genetic_truth/, counts.tsv, <assay>/
```

| Area | Location | Notes |
| --- | --- | --- |
| Routes | `inst/suite/simulators.tsv`, `R/suite_routes.R` | One row per route: spec field, values (Biodes' canonical ids, normalised to lower case with `_`), simulator, status. Several genetics rows: the cross wins over the panel. eQTL turns on the panel. |
| Parameters | `inst/suite/simulator_parameters.tsv` | Argument, default, type, spec field that overrides it. Order: `params` > spec field > table default > the function's own default. `simitall_simulator_parameters(sim)` lists every argument. |
| Outcome | `R/suite_simulate.R` | Effects from `analysis.model` (treatment levels spread to the minimum detectable difference, block/random variance, subsamples, attrition), plus main effects for `analysis.terms.uncovered_factors` (e.g. tissue) and a shared effect for the rows of one unit, then the outcome family. |
| Genetics | `R/suite_genetics.R`, `inst/suite/genetics.json`, `mating_schemes.tsv`, `phenotype_models.tsv` | Genotype term = treatment term whose `analysis.terms` meaning matches `genotype_pattern`. Panel: Balding-Nichols, inbred or outbred by `system`. Cross: `generate_random_haplotype_panel()` + `simulate_breeding()`. Genetic values (additive / dominance / epistasis / GxE) scaled to plot-level h2 = Vg / (Vg + residual SD^2); they replace the even-spread genotype effects. |
| Assays | `R/suite_assays.R` | One library per unit, named by `sample_id`. Condition = first non-genotype treatment term, batch = first block term. One count assay owns `counts.tsv` (eQTL > RNA-seq > single-cell); the others write `<assay>/unit_counts.tsv`. |
| Layout | Biodes `studyforge::layout_from_spec()` | Simitall builds no layouts. Biodes keeps every factor as its own column (factorial/split crossing, within-unit factors such as tissue repeat each unit's row) and stops with `missing: ...` when the spec cannot support a layout; Simitall reports those as questions for the user. Simitall's old builder is in `old_model/`. |
| From a question | `R/suite_question.R` | `simitall_from_question()` asks Biodes (`studyforge_ask`), builds the layout, simulates, fits. |
| Fit, power, figures | `R/suite_analyze.R`, `R/suite_power.R`, `R/suite_figures.R` | Bioana, power sweep, Biopaint. |
| Simulators | `R/*_api.R`, `R/internal_*.R`, `R/peak_studies.R`, `R/scrna_studies.R`, `R/breeding_*.R` | The package simulators the routes call. |
| Old reader + Shiny app | `old_model/` | **Paused.** Saved copy of the removed prompt readers and the Shiny agent; excluded from the package. |

### Adding a simulator

1. Add a row to `simulators.tsv` keyed on a Biodes field value (ask Biodes for
   the field if it does not exist).
2. Add its arguments to `simulator_parameters.tsv`.
3. Add a runner in `R/suite_assays.R` (`.suite_assay_runner()`), writing one
   library per unit named by `sample_id` and returning its truth.
4. Add a test in `tests/testthat/test-suite-simulators.R` that builds a design
   with `make_design_folder()` and checks unit alignment and truth.

Known gaps: methylation and microbiome (`planned`), survival and compositional
outcomes (`planned`); see `docs/SIMULATION_COVERAGE.md`; genotype-factor level order follows natural sort, so the
reference level is the alphabetically first one.

## Data-Source Policy

For real chromosome, organism, gene, or cohort claims, ask for or use a
compatible real panel/FASTA/VCF/GFF3/GTF/recombination map as appropriate. A
synthetic chromosome or panel must be labeled synthetic and never presented as
the real organism's genome.

## Validation Contract

Every major workflow should define these four pieces:

1. **Truth**: known causal loci, ancestry, peaks, expression program, or other
   programmed mechanism.
2. **Required results**: the tables/files that the simulation or analysis must
   create.
3. **Required figures**: diagnostic plots appropriate to the requested method.
4. **Pass criteria**: file existence plus meaningful design/mechanism checks.

Examples:

| Workflow | Required figures | Core validation |
| --- | --- | --- |
| GWAS | Manhattan, QQ, phenotype distribution | sample count, valid allele frequencies, causal-truth recovery/calibration where truth exists |
| Bulk RNA-seq | PCA/sample QC, volcano/MA, DEG summary | count dimensions, condition/replicate integrity, DE truth recovery if simulated |
| scRNA-seq | cell QC, UMAP/cluster view, marker or pseudobulk summaries | donor/cell counts, cell-type truth recovery, no pseudoreplication |
| ChIP/ATAC-seq | peak/QC, differential peaks, enrichment | control/replicate integrity, peak truth precision/recall when simulated |
| Breeding | ancestry/LD/recombination or trait/selection plots | valid pedigree, expected segregation, diversity/relatedness and gain checks |
| Human family risk | probability comparison, trait distribution, homozygosity | relationship order and probability mechanism agree with the programmed model |

Validation specifications live in
`inst/agent/knowledge/validation_contracts/`. Reproducible figures and
fixed-seed demonstrations belong in `docs/VALIDATION.md`, not as a large result
figure on the root README.

## Figures and Outputs

- Put generated diagnostic figures in an `figures/` subdirectory of the run
  output directory.
- Produce individual figures as well as a combined overview when useful.
- Write tabular source data/truth/metrics next to figures whenever practical.
- Keep generated outputs out of the package source tree unless they are a
  small, intentional fixed demo.


## Reference knowledge

`inst/agent/knowledge/` (papers, species, tools, techniques, workflows,
validation contracts) and `inst/agent/paper_records/method_records.tsv` are
kept as reference material; no code reads them now. Only rows with
`review_status == "reviewed"` may be used as evidence. They are candidates to
move to Biodes/BioCore.

## Gene Annotation and GO

### Lightweight package data

- `inst/extdata/annotations/`: tiny GFF3 gene-only examples for demos.
- `inst/extdata/gene_coordinates/`: compact full gene-coordinate catalogues for
  Arabidopsis, barley, rice, wheat, maize, and fly.
- `inst/extdata/go/arabidopsis_tair10_go_demo.gaf.gz`: a small matched GAF demo.
- `inst/extdata/annotations/go_sources.tsv`: GO source and identifier guidance.

The compact coordinate tables retain `seqname`, `start_bp`, `end_bp`, `strand`,
`gene_id`, `source_gene_id`, `gene_name`, and `gene_biotype`. They are suitable
for interval overlap and ID joins, not transcript/exon/CDS analysis.

### GO API

```r
coordinates <- read_gene_coordinates("Arabidopsis_thaliana")
annotations <- read_go_gaf("/path/to/ARATH-mod.gaf.gz")
gene_go <- join_gene_go(coordinates, annotations)
term_summary <- summarize_go_terms(gene_go, by = "go_id_aspect")
```

`join_gene_go()` favors stable IDs. In `auto` mode it only falls back to an
exact symbol for unmapped genes; it does not use loose synonyms. Synonym joins
need a reviewed organism-specific ID mapping. In particular, the included
maize and wheat GAF downloads use UniProt-oriented IDs, so a compatible
Ensembl/NAM-5.0 or IWGSC-to-UniProt mapping is required.

### Full GO resources

`resources/go_annotations/` is a source-controlled copy of the 21 local GAF
downloads supplied on 2026-10-07. It includes checksums and header dates in
`manifest.tsv`. It is intentionally excluded from the built R package through
`.Rbuildignore`; keep it in GitHub, not `inst/extdata/`.


## Testing and Checks

```r
devtools::document()                          # after changing exported R API
devtools::test(filter = "suite")              # Biodes bridge, routes, simulators
devtools::test(filter = "suite-simulators")   # outcome / genetics / assay routes
devtools::test(reporter = "summary")          # everything, before a release
```

Optional packages such as `simstudy`, `pedtricks`, `rrBLUP` and `edgeR` may
cause optional tests to skip. State this plainly; do not report skipped tests
as passing coverage.
