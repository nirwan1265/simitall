# simitall Agent Knowledge Base

This directory is a retrieval-ready knowledge base for the `simitall` agent.
It supports explanation and workflow planning; it does not replace the package
documentation or silently run analysis. Each workflow must remain reproducible
as direct R code with recorded parameters, random seeds, inputs, outputs, and
truth sets.

## Directory Contract

| Directory | Contents | Evidence standard |
|---|---|---|
| `papers/` | One faithful, source-traceable summary per paper | Cite page/section evidence and preserve uncertainty |
| `techniques/` | Cross-paper method guides | Separate accepted practice, package behavior, and inference |
| `tools/` | Package/software cards plus tutorial map | Record tutorial/vignette provenance and executable boundary |
| `species/` | Organism-aware modeling notes | Do not transfer biological defaults across organisms without saying so |
| `workflows/` | Reproducible recipes using `simitall` | Include inputs, outputs, code, QC, truth benchmarks, and limitations |
| `validation_contracts/` | Four-part workflow checks | Define simulation truth, required results, figures, and pass criteria |
| `terminology/` | Definitions and aliases | Prefer unambiguous names and state context-dependent meanings |
| `evaluations/` | Agent behavior tests | Specify the expected safe, evidence-aware response |

## Source Priority

1. The user's input data, parameters, and experimental design.
2. `simitall` function documentation and installed software help.
3. Primary papers, official manuals, and package vignettes.
4. Review papers for context, not for unverified numeric defaults.

## Agent Rules

- Never claim that a simulated association is a biological discovery.
- Distinguish reference-aware truth benchmarking from real-data validation.
- Ask for a species, ploidy, population design, and intended endpoint when
  those change the scientific interpretation.
- Do not apply human, plant, bacterial, or single-cell defaults universally.
- Treat a named panel as a versioned data resource, not as a synonym for a
  species. Preserve its sample manifest, genome build, marker filters, and
  population/family labels.
- Distinguish a pedigree graph, a genomic-kinship dendrogram, and a
  phylogenetic tree. They require different inputs and answer different
  questions.
- For admixed human studies, distinguish ordinary structure-adjusted GWAS,
  admixture mapping, and local-ancestry-aware association. Never claim that
  local ancestry has been inferred unless phased cohort data, reference panels,
  and an explicit inference method are provided.
- Surface software limitations before proposing a workflow. For example,
  Unicycler is aimed at bacterial assemblies and PBSIM3 multi-pass output is
  not final HiFi consensus until processed by PacBio CCS software.

## Document Metadata

New `.Rmd` documents should begin with YAML fields for `id`, `document_type`,
`status`, `domains`, `related_functions`, `source_priority`, and
`last_reviewed`. The `papers/` documents already use a richer paper-specific
schema.
