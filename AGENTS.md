# SimItAll: AI Agent Handoff

**Read this before making changes.** This repository is actively edited by the
project owner with more than one coding assistant. The working tree can be
intentionally dirty. Preserve existing work, inspect `git status` first, and
do not revert, reset, rename, or reformat unrelated changes.

Last handoff update: 2026-10-07.

## Project Purpose

SimItAll is an R package and local Shiny application for two connected jobs:

1. **Simulate** genomic, population, phenotype, sequencing, and multi-omics
   data while retaining truth.
2. **Analyze and validate** synthetic or user-supplied data with appropriate
   figures, result tables, and explicit pass criteria.

The Shiny interface is called the **SimItAll agent**. It accepts a plain
language study question, determines the supported workflow, identifies missing
inputs, shows reviewable R code, and runs the displayed code only when the user
presses **RUN**.

## Collaboration Rules

- Treat current working-tree changes as another collaborator's work unless the
  user explicitly asks to change them.
- Never use `git reset --hard`, `git checkout --`, or destructive cleanup.
- Do not invent biological results, package availability, input files, or GO
  assignments.
- Use an explicit synthetic input only when the user asks to simulate/generate
  data or selects the synthetic data source. Otherwise request compatible real
  inputs.
- `Use package demo only` means use packaged demo resources; it does not mean
  use stale code from a previous prompt. The user must press **Send** after
  changing prompt or data source so the recipe is regenerated, then press
  **RUN**.
- A user request for GWAS must not silently run RNA-seq, selection, ChIP-seq,
  or any other unrequested workflow. Parse the requested analyses first.

## Main Architecture

| Area | Main locations | Notes |
| --- | --- | --- |
| Public R API | `R/` | Roxygen comments generate `NAMESPACE` and `man/`; run `devtools::document()` after exported API changes. |
| Agent entry point | `R/agent_api.R` | `simitall_ask()` returns plans, preflight information, and reviewable code. |
| Structured prompt compiler | `R/agent_spec.R` | Preferred route for supported prompts: extract -> validate -> resolve inputs -> compile plan and recipe from one spec. |
| Legacy/generic workflow routing | `R/agent_workflow_compiler.R` | Use for broader workflow compilation where supported. |
| Shiny application | `inst/shiny/simitall-agent/app.R` | Supports local, Ollama, and OpenAI modes, prompt tabs, tuning controls, local run status, uploads, and output downloads. |
| Fixed demo/validation runners | `analysis/paper_fig/` | Reproducible demonstrations such as human IRF6 and NAM GWAS. They are examples/validation runners, not the only way the agent should answer a new question. |
| Agent knowledge | `inst/agent/knowledge/` | Local Rmd knowledge records used for retrieval and planning. |
| Prompt aliases | `inst/agent/aliases.tsv` | Audited wording/inflection aliases, for example `inbred -> inbreeding` and `epistatic -> epistasis`. |
| Tests | `tests/testthat/` | Add a focused test for new user-visible behavior. |

## Prompt and Recipe Contract

The agent is not meant to merely select a hard-coded script from a question.
For a supported prompt, it should build a validated specification and compile
the plan and R code from that same specification.

The desired flow is:

```text
User wording
  -> spelling/word-form normalization
  -> audited aliases
  -> structured concept extraction
  -> input resolution and validation
  -> recipe + plan from the same specification
  -> explicit RUN by the user
  -> outputs, figures, truth, and pass criteria
```

Use fixed scripts only as versioned demos or validation references. Do not make
new scientific capabilities depend on adding one script per possible prompt.

### Current Supported Structured Study Families

- NAM and diversity-panel GWAS, quantitative traits, and genomic selection.
- Human family-risk and relatedness/IRF6 demonstration workflows.
- Bulk RNA-seq and eQTL workflows.
- Single-cell RNA-seq, pseudobulk, cell-type eQTL, and relevant benchmarks.
- ATAC-seq and ChIP-seq peak workflows.
- Breeding populations, selfing, backcrossing, pedigree effects, and phenotype
  distributions.
- Genome/read simulation, assembly, and validation.
- Gene-coordinate and Gene Ontology annotation joins.

If a requested prompt is not supported, say so clearly and provide a plan or
ask for the missing data. Do not fabricate runnable code.

## Data-Source Policy

The Shiny data-source choices have distinct meanings:

| Choice | Meaning |
| --- | --- |
| Ask for required data | No simulation code requiring real inputs should run until compatible inputs are supplied. |
| Generate synthetic inputs | Generate toy/reference-like data only when the user explicitly requested simulation. |
| I supplied local files | Use the stated local input paths after checking compatibility. |
| Use package demo only | Use the small packaged demo panels/maps/annotations. |

For real chromosome, organism, gene, or cohort claims, ask for or use a
compatible real panel/FASTA/VCF/GFF3/GTF/recombination map as appropriate. A
synthetic chromosome must be labeled synthetic and never presented as the real
organism's genome.

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
- Make figures downloadable through the Shiny run-output archive; do not force
  the app to display every figure inline.
- Write tabular source data/truth/metrics next to figures whenever practical.
- Keep generated outputs out of the package source tree unless they are a
  small, intentional fixed demo.

## Knowledge Base

The local knowledge base is intentionally more than a paper list. It includes
papers, species notes, tools/tutorials, techniques, workflows, evaluations,
terminology, and validation contracts.

Important directories:

```text
inst/agent/knowledge/papers/
inst/agent/knowledge/tools/
inst/agent/knowledge/techniques/
inst/agent/knowledge/workflows/
inst/agent/knowledge/validation_contracts/
inst/agent/paper_records/method_records.tsv
```

Paper-record evidence labels matter:

- `metadata_curated_summary`: routing-level title/metadata curation; not a
  full-text paper extraction and not a source of paper-specific numeric claims.
- `source_reviewed_from_user_pdf` and `citation_and_user_manual_reviewed`:
  may support only the explicitly documented source-specific details.
- `user_vignette_reviewed`: tutorial/vignette-grounded operational notes.

When adding a term users may phrase differently, add an auditable alias in
`inst/agent/aliases.tsv`; do not rely on an LLM alone to infer every inflection.

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

## Documentation Layout

- `README.md`: concise public project overview; no large example figure or
  giant hard-coded recipe.
- `docs/VALIDATION.md`: reproducible validation figures, expected outcomes, and
  pass criteria.
- `docs/TUTORIAL.md`: detailed runnable tutorial.
- `inst/extdata/README.md`: package-demo data boundaries.
- This file: operational handoff for coding agents.

## Testing and Checks

Run focused tests while developing:

```r
devtools::document()                      # after changing exported R API
devtools::test(filter = "agent-spec")    # prompt/spec routing
devtools::test(filter = "extdata")       # packaged annotations and GO demo
devtools::test(filter = "scrnaseq")      # single-cell behavior
devtools::test(filter = "sequencing-profiles")
```

Then run the appropriate full suite before a release:

```r
devtools::test(reporter = "summary")
```

Optional packages such as `simstudy` and `pedtricks` may cause optional tests
to skip. State this plainly; do not report skipped tests as passing coverage.

## Current Product Decisions

- Local mode is the default and does not need an API key.
- Ollama and OpenAI modes may improve explanation, but never bypass input
  validation or code safety checks.
- The user wants simple, keyword-faithful outputs. For example, “binomial
  phenotype GWAS in NAM” means a binary-trait GWAS, not a generic F1/F2/F3
  demo with unrelated extra analyses.
- Recipes should run from a source checkout or installed package without
  writing a mirrored absolute directory tree into Downloads. Prefer an
  explicit, user-selected output directory and relative paths inside it.
- Real human outputs are illustrative synthetic research demonstrations, never
  clinical prediction or genetic counseling.
