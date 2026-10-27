# SimItAll: AI Agent Handoff

**Read this before making changes.** This repository is actively edited by the
project owner with more than one coding assistant. The working tree can be
intentionally dirty. Preserve existing work, inspect `git status` first, and
do not revert, reset, rename, or reformat unrelated changes.

Last handoff update: 2026-10-26.

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

### Next Architecture: StudySpec + Capability Registry

The agent must generalize from a prompt by compiling a **runtime StudySpec**,
not by picking one hard-coded recipe from isolated keywords. A StudySpec is a
plain JSON-compatible R list with fields such as `primary_analysis`,
`organism`, `assembly`, `region`, `data_source`, `inputs`, `design`, `modules`,
`parameters`, `outputs`, `figures`, `truth`, `pass_criteria`, `assumptions`,
and `unresolved_questions`.

Use separate, auditable static registries for different jobs:

| Registry | Suggested location | Purpose |
| --- | --- | --- |
| Alias lexicon | `inst/agent/aliases.tsv` | Map wording and inflections to canonical concepts only. |
| StudySpec schema | `inst/agent/capabilities/study_spec.schema.json` | Define valid fields and status for a partial or final StudySpec. |
| Analysis catalog | `inst/agent/capabilities/analysis_catalog.json` | Map aliases to analysis families, subanalyses/goals, candidate models, typed parameters, inputs, results, figures, QC, tools, and support status. |
| Capability registry | `inst/agent/capabilities/workflows.json` | Map an analysis/design to required inputs, compatible modules, safe defaults, outputs, figures, and validators. |
| Input-type registry | `inst/agent/capabilities/input_types.json` | Describe VCF/BCF, GTF/GFF3, FASTA, map, count matrix, BED, and metadata compatibility checks. |
| Sample naming contract | `inst/agent/capabilities/sample_naming.json` | Give new experimental data a reversible sample-ID convention and authoritative metadata columns. |
| Default policy | `inst/agent/capabilities/defaults.json` | Declare defaults and whether each may be applied automatically or requires clarification. |
| Tutorial map | `inst/agent/capabilities/tutorial_map.json` | Map concepts to local tutorials, papers, workflows, validation contracts, and tool cards for retrieval. |
| Reviewed examples | `inst/agent/capabilities/study_spec_examples.jsonl` | Prompt-to-expected-StudySpec cases for regression evaluation and later fine-tuning. |

Registries contain identifiers and metadata, never executable R expressions.
R code resolves approved module identifiers to package functions/templates.
This keeps configuration reviewable and prevents a registry record from
becoming arbitrary code execution.

#### Agreed planning model

Treat organism/reference, population or experimental design, and analysis as
separate dimensions of one request, rather than baking a species into every
workflow. Identify the goal (simulate, analyze, or both), organism and genome
assembly/region where applicable, input source, population or experiment design
where relevant, and explicitly requested modules. Analysis definitions should
eventually provide typed parameters, allowable values, conditional requirements,
units, defaults and provenance, required inputs, compatibility checks,
outputs, and implemented-versus-planned status. Species rules apply where
biologically relevant (e.g. ploidy, assembly, mating feasibility); common
RNA-seq steps can be shared without assuming references are interchangeable.
Only add prerequisite modules when necessary and disclose them; never run
unrequested independent analyses. Preserve supplied values and ask targeted
questions when a consequential scientific choice is missing.

For example, a synthetic human chr10 first-cousin inbreeding request might
specify `generations = 10` and request pedigree-based inbreeding coefficients
plus marker homozygosity. Population size and exact pedigree construction must
be resolved explicitly. Generations are not interchangeable with the current
multi-generation `n_rounds`; marker homozygosity is not itself the pedigree
coefficient. A synthetic marker panel is not real chromosome data. The current
catalog now supplies a planning-only generic resolver for unsupported
structured prompts: it identifies a direct catalogued analysis, preserves only
explicitly requested secondary analyses, exposes parameters/required choices,
and does not emit code. It does not yet implement typed final-spec validation
or general-purpose executable composition.

### Current implementation plan: generic composition, not prompt functions

Do **not** add one R function, script, or parser branch for each new natural-
language prompt. The next implementation work is to make the registry-driven
compiler generic: extract a partial StudySpec, select reusable primitive modules
(e.g. population initialization, pedigree construction, meiosis/recombination,
marker/ROH summaries, and plotting), validate their typed parameters, then
compose an approved recipe. A prompt-specific example belongs in tests, not as
an API function.

`pedigree_inbreeding` in `workflows.json` is therefore a planned generic module,
not a human-chr10 or first-cousin implementation. It declares organism,
assembly/region when needed, data source, mating scheme, generations, population
size, marker count, and seed. The compiler must distinguish pedigree F, marker
homozygosity, and ROH-based measures and must ask about a missing mating scheme
rather than choosing one. First-cousin, sibling, and selfing examples should all
exercise the same generic primitive modules after they exist.

Build in this order:

1. Keep the machine-valid StudySpec schema and typed input/default registries
   synchronized with catalog records.
2. Expand the planning-only catalog resolver into a final-spec validator that
   handles conditional inputs and clarification answers.
3. Implement generic pedigree and inheritance primitives with explicit inputs,
   rather than human- or chromosome-named workflow functions.
4. Select only compatible primitives and report
   missing required choices.
5. Compile validated specs into recipes and use prompt examples only as
   regression tests.

### Analysis-catalog rules

`analysis_catalog.json` is the large general knowledge layer. It currently
records broad analysis families and subanalyses across pedigree/inbreeding,
breeding, phenotype simulation, GWAS, linkage QTL, bulk RNA-seq,
time-series RNA-seq, alternative splicing, ASE, eQTL/sQTL, single-cell RNA-seq,
ChIP-seq/ATAC-seq, genomic selection, population genetics, sequencing,
assembly, annotation/GO, multi-omics, and metagenomics. It is deliberately
extensible rather than a claim that every conceivable biological analysis is
already represented.

The expected composition is:

```text
aliases -> canonical concept -> analysis family -> subanalysis/goal
        -> candidate model -> compatible input contracts + typed parameters
        -> outputs + figures + QC -> execution-status decision
```

For example, `RNA-seq` alone must not select a model. The resolver must decide
whether the user requests bulk count simulation, differential expression,
time-series trajectories, alternative splicing, ASE, eQTL, or power analysis;
then it must require the corresponding experimental design and compatible
counts/reads/annotation/metadata. `inbreeding` likewise identifies a family but
cannot choose selfing, sibling mating, cousin mating, bottleneck, or a custom
pedigree. The registry’s `execution_status` is decisive: `partial` can only
offer an implemented subset after validation; `knowledge_catalog_only` can
explain and scope a study but cannot produce runnable code; `planned` must
report a gap.

For user-supplied experimental data, recommend the reversible naming pattern
`S001_R01_TrtControl_Ti00`: `S` is a unique sample/library, `R` is a
biological replicate, `Trt` is treatment/condition, and `Ti` is time. Never
use a bare `T`, which is ambiguous. Omit segments that do not apply, and keep
a metadata table as the authority with `sample_id`, `biological_replicate`, and
applicable `treatment`, `time`, `batch`, `subject_id`, and
`technical_replicate` columns. Existing cohorts need not be renamed: map their
original IDs in that table. Count-matrix columns, alignment sample names, and
metadata `sample_id` must match exactly.

Local tutorials and reviewed papers are retrieval evidence, not training by
themselves. `study_spec_examples.jsonl` supplies reviewed prompt-to-JSON
targets for regression evaluation; build a sufficiently reviewed dataset before
using it for an Ollama fine-tune. If a local LLM is used first, it may propose a
partial StudySpec only. Deterministic alias, catalog, input-contract, and
support-status validation remains the authority.

#### Tutorial-first catalog maintenance

Before adding or revising an analysis record, read the most relevant local
tutorial, workflow, and validation contract under `inst/agent/knowledge/`.
Capture the documented analysis unit, biological replication unit, required
design variables, candidate models, expected results, figures, QC checks, and
scope limits in the catalog. A tool name or a paper title alone is not enough
evidence for a runnable route.

For example, the local RNA-seq material establishes that bulk count simulation
needs an explicit experimental design and biological replicates; time-course
analysis needs time, condition/group, a declared time representation, and
subject-aware handling when samples repeat; bulk eQTL needs donor-aligned
genotype, expression, and metadata; and single-cell eQTL must keep the donor,
not the cell, as the biological genetic replicate. These details belong in the
typed StudySpec and input contracts, not in a prompt-specific script.

Use external primary documentation only to fill a gap in the local knowledge
base. Record it in `inst/agent/capabilities/external_sources.json` with its
scope, and never let it imply that a tool is installed or runnable. The local
tutorial map must continue to point at the exact supporting files.

#### Catalog acceptance gate

A catalog entry is not a working agent capability until all three checks pass:

1. Retrieval returns the relevant local tutorial/validation evidence without
   pulling in an unrelated workflow.
2. The generic extractor and resolver produce a partial StudySpec with the
   right primary analysis, only the explicitly requested secondary analyses,
   and an explicit list of unresolved scientific choices or missing files.
3. A route is marked runnable only after generic primitives compile it,
   validation figures/truth/pass criteria are emitted, and a regression test
   passes. A `knowledge_catalog_only` or `planned` route must stop at a plan.

Do not treat a broadly similar legacy fallback as success. For example, a
plain bulk RNA-seq count-simulation request must not silently add eQTL, GWAS,
selection, or a maize NAM workflow merely because those terms occur in a
retrieved document. Add the observed prompt and expected partial StudySpec to
`study_spec_examples.jsonl` when a mismatch is found; use it as a regression
case after review.

#### Partial StudySpec behavior

`simitall_ask()` now attaches `catalog_study_spec` for recognized catalog
requests. It is a schema-compatible, `do_not_run = TRUE` record that preserves
the question, recognized analysis, input contracts, supplied paths, resolved
parameters, expected outputs/figures, assumptions, and unanswered choices.
For a broad RNA-seq request, the catalog's `clarification_mode = "ask_goal"`
returns a targeted goal choice in Shiny instead of guessing a count, DE,
time-course, splicing, ASE, eQTL, or power workflow. Selecting a goal sends a
new explicit request back through the same generic resolver.

#### Compiler stages and precedence

```text
prompt + selected data source + uploaded paths
  -> normalize spelling and aliases
  -> identify primary analysis and secondary analyses
  -> extract a partial StudySpec
  -> enrich it from capabilities and input metadata
  -> validate compatibility and resolve defaults
  -> ask targeted questions or report required files
  -> compile plan, recipe, figures, truth, and pass criteria from the final spec
```

Explicit user details and actual supplied file types outrank broad keyword
matches. For example, `phased VCF + chr14 map + GTF + 10 generations` must
route to chromosome/pedigree/annotation analysis; it must not be replaced by a
generic human family-risk demonstration merely because the wording contains
`first cousin` or `probability`.

Aliases identify concepts; they do **not** choose a scientific design by
themselves. `inbred -> inbreeding` should select the `inbreeding` concept, then
the compiler must distinguish selfing, sibling mating, cousin mating,
population bottlenecks, and user-specified pedigrees from the rest of the
StudySpec. Do not equate `inbreeding` with one fixed mating scheme.

#### Defaults versus clarification

The current agent can state required inputs in its plan/preflight, but it does
not yet maintain a complete spec-driven clarification dialogue. Build this
next. Apply only transparent technical defaults automatically, and record them
under `assumptions` in the displayed plan:

- Safe automatic defaults: seed, output directory, figure directory, file
  naming, and standard QC/validation figures for an explicitly named analysis.
- Defaults that require clear disclosure: number of replicates, marker density,
  sample size, and synthetic population parameters.
- Values that require a user choice or compatible input: real-vs-synthetic
  source, genome assembly, mating/pedigree design, effect/trait model, causal
  locus assumptions, clinical interpretation, and the definition of an
  ambiguous outcome such as “most affected gene.”

For an unresolved scientific choice, return a short targeted question with a
safe provisional interpretation; do not create runnable code until it is
resolved. Missing required files must produce an explicit input checklist, not
a silent synthetic substitution.

#### Required regression example

The human chr14 prompt is a required agent regression case. Given a phased
chr14 founder VCF, matched recombination map, and GRCh38 GTF, the compiled spec
must select pedigree inheritance, recombination, ROH/homozygosity, gene
annotation, and gene-burden ranking. It must request absent files, define or
ask for the gene-ranking metric, label results as synthetic research outputs,
and never emit a generic disease-risk-locus recipe.

**Observed baseline gap (2026-10-16):** the structured compiler currently
marks this chr14 prompt unsupported, but the legacy local-plan fallback can
still emit an input-dependent generic human **chr10** inbreeding template. This
is incorrect: it does not use the requested chromosome, phased haplotypes,
GTF, ROH, or gene ranking. Until the new capability is implemented, this input
shape must stop with an explicit unsupported/required-input response and no
recipe. Add a regression test for both the compiler result and the full
`simitall_ask()`/Shiny Send response.

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

### Capability registry groundwork

`inst/agent/capabilities/workflows.json` is a curated, versioned inventory of
workflow concepts, candidate aliases, available package functions, conceptual
inputs, questions, recommended outputs/figures/checks, guardrails, and links to
existing validation contracts and paper summaries. See the adjacent README for
field meanings and `.simitall_agent_capability_registry()` for the internal
loader. `tests/testthat/test-agent-capability-registry.R` checks function and
resource references and key scientific distinctions.

`analysis_catalog.json` is connected to a generic **planning-only** resolver
used by `simitall_ask()` when the older structured compiler does not support a
request. It can select an explicit primary analysis, preserve explicit
secondary analyses, identify required choices, and prevent a legacy fallback
from adding unrelated work. It does not execute a catalog entry or promise
code. `workflows.json` remains a static inventory and is not a general recipe
compiler. Entries marked `partial` do not guarantee a complete recipe;
`planned` entries are not executable routes.
Recommended figures and results are not promises. In particular, the human
chr14 multigeneration ROH/gene-ranking request is still unsupported, and the
existing legacy fallback must not be mistaken for a correct chr14 workflow.
The older generic chr10 inbreeding fallback uses three rounds and marker
homozygosity without pedigree F. Do not present it as satisfying a coefficient request.
Paper summaries provide background, not evidence a tool is installed. Keep
real, synthetic, and bundled-demo inputs distinct, and do not equate
inbreeding with a particular mating scheme. Future work should connect the
registry to structured extraction, conditional input checks, clarification,
safe routing, and tests for missing inputs and unsupported requests.

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
devtools::test(filter = "agent-capability-registry") # registry and references
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
