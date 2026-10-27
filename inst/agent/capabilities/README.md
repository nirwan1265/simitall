# Study-Planning Knowledge Catalog

These files define the knowledge layer for general prompt interpretation. They
do **not** add one R function per user prompt. `analysis_catalog.json` is
currently connected to a generic planning-only resolver used by
`simitall_ask()` when its older structured compiler does not support a request.
That resolver identifies a direct analysis match, preserves only explicit
secondary analyses, and lists required choices without producing a recipe.
Adding a JSON entry never enables code execution.

| File | Role |
| --- | --- |
| `analysis_catalog.json` | Analysis families, subanalyses, aliases, goals, candidate models, typed parameters, inputs, outputs, figures, QC, tool candidates, references, and execution status. |
| `input_types.json` | VCF/BCF, FASTA, GTF/GFF3, maps, count matrices, metadata, alignments, pedigree, and other input contracts. |
| `defaults.json` | Transparent defaults, questions that must be asked, and assumptions that must never be inferred. |
| `sample_naming.json` | Reversible naming for new sample IDs and the authoritative metadata-column contract for count, alignment, and phenotype data. |
| `study_spec.schema.json` | JSON Schema for a partial/final prompt-derived StudySpec. |
| `tutorial_map.json` | Retrieval map from canonical concepts to local tutorials, validation contracts, tools, workflows, and paper summaries. |
| `external_sources.json` | Method documentation consulted where a concept was not adequately represented locally. It never declares a tool installed. |
| `workflows.json` | Existing inventory of package-level routes and their current support status. |

The internal `.simitall_agent_analysis_catalog()`,
`.simitall_agent_input_type_registry()`, and `.simitall_agent_default_policy()`
load the primary static catalogs. They are deliberately internal until the
StudySpec resolver has a stable public API.

`analysis_catalog.json` is organized as:

```text
user terms -> canonical analysis family -> subanalysis/goal -> candidate model
           -> required inputs + typed parameters -> results + figures + QC
```

For example, “RNA-seq” routes to the `bulk_rnaseq` family only provisionally;
the resolver must distinguish count simulation, differential expression,
time-series analysis, alternative splicing, ASE, eQTL, or power analysis.
The user’s goal and compatible input types decide the route. Similarly,
“inbreeding” identifies a concept but cannot choose selfing, sibling mating,
cousin mating, or a custom pedigree on the user's behalf.

`partial` means listed package functions exist but does not certify a complete
recipe. `knowledge_catalog_only` means a method can be explained and scoped,
but not run. `planned` means no execution support is claimed.

Paths under `knowledge_sources` and `paper_sources` are relative to
`inst/agent/knowledge/`. Paper summaries provide scientific context, not proof
that a particular backend is installed or implemented. The validation
contracts are recommendations, not guarantees of produced files or statistical
success. Existing `inst/agent/aliases.tsv` remains separate; catalog aliases
are candidate workflow cues, not unconditional synonym replacement.

## Sample Metadata Convention

For a new experiment, use the recommended ID pattern
`S001_R01_TrtControl_Ti00`. `S` is a unique sample/library, `R` is the
biological replicate, `Trt` is treatment/condition, and `Ti` is time. Do not
use a bare `T`, because it cannot distinguish treatment from time. `B01` may
denote batch and `Sub012` may denote a repeated-measures subject/donor.

Omit non-applicable segments and retain a metadata table as the source of
truth. Its required columns are `sample_id` and `biological_replicate`; add
`treatment`, `time`, `batch`, `subject_id`, and `technical_replicate` when
applicable. Existing studies may preserve their original IDs if this table maps
them cleanly. Matrix columns and metadata `sample_id` values must match exactly.

Next integration should extend the current partial resolver into a final
StudySpec validator: normalize aliases, retrieve records, validate conditional
requirements and supported combinations, ask about consequential missing
choices, disclose defaults, and refuse unsupported routes rather than
substituting a different analysis or source. Only then may generic primitives
compile a validated route into a recipe. Tutorials and papers teach the
LLM/retriever terminology and rationale; reviewed prompt-to-StudySpec examples
are needed before any fine-tuning. An LLM may propose JSON, but the deterministic
resolver and contracts must remain the execution authority.

The current partial resolver returns a schema-compatible `catalog_study_spec`
with `do_not_run = true`. For families marked `clarification_mode: "ask_goal"`,
such as broad bulk RNA-seq, it asks the user to choose a goal before selecting a
subanalysis. This is a planning interaction only; no catalog record can run
until generic primitives and final-spec validation exist.

## Curation and Acceptance

Curate a record from the local knowledge base first. Read the related tutorial,
workflow, and validation contract before describing inputs, replication units,
models, outputs, figures, or pass criteria in `analysis_catalog.json`. Use an
external primary method document only when local material lacks that detail;
record its scope in `external_sources.json`. Neither kind of evidence proves
that a backend is installed or implemented.

For every reviewed prompt in `study_spec_examples.jsonl`, test three layers:

1. Retrieval returns relevant evidence.
2. Generic extraction/resolution emits the expected partial StudySpec and
   preserves only the analyses requested by the user.
3. Only an entry with implemented generic primitives may compile and run.

A failed or unsupported route is useful evaluation data, not a reason to add a
one-off function or silently substitute a related workflow. Update the reviewed
example with the intended partial spec, then implement reusable primitives and
the resolver rule that makes it pass.
