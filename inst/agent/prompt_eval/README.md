# Prompt -> StudySpec evaluation set

This folder holds the examples used to measure how well the SimItAll agent
reads a plain-language request, and later to train a small prompt-reading
model. It follows the AGENTS.md rule: build a sufficiently *reviewed* dataset
before any fine-tuning; deterministic aliases, catalog and validation remain
the authority.

| File | Purpose |
|---|---|
| `study_spec_examples_draft.jsonl` | Draft examples (`review_status: draft`). Reviewed rows can later move to `capabilities/study_spec_examples.jsonl`. |
| `coverage_map.json` | Target number of examples per catalog analysis, combinations, traps and not-yet-supported techniques, and how many exist. |
| `review_batch1.md`, `review_batch2.md` | Human-readable review sheets (batch 1: genetics/breeding; batch 2: omics, GWAS and inbreeding across species, combinations, traps, unsupported). |
| `score_examples.R` | Runs the local agent on every example and scores it (`Rscript inst/agent/prompt_eval/score_examples.R [dev|test]`). |
| `results/` | Latest score tables (regenerated; safe to delete). |

Generators live in `data-raw/agent_eval/` (one script per batch). `score_lexicon.R [dev|test]` compares organism/technique recognition with the current agent; `score_reader.R [dev|test]` scores the full lexicon reader, the agent and both merged with the same scorer; `score_model.R preds.json [dev|test]` scores lexicon vs lexicon + trained model. `review_sheet.py <batch>` rebuilds a review sheet and the `have` counts in `coverage_map.json`; run `score_examples.R` first to include what the current agent picks.

## Row format

```json
{"id": "...", "batch": 1, "split": "dev|test", "review_status": "draft|reviewed",
 "workflow": "<catalog subanalysis | combination | trap | unsupported>",
 "situation": "simulate|provided|ask|demo|trap|unsupported",
 "style": "expert|beginner|terse|chatty|nonnative",
 "source": "tutorial section, workflow Rmd or AGENTS.md rule it is grounded in",
 "prompt": "...",
 "expected_spec": { ...only the StudySpec fields the prompt determines... },
 "must_ask": ["dotted fields the agent must not fill silently"],
 "must_not": ["analyses that must not appear in the plan"],
 "notes": "..."}
```

Conventions: `organism` is a lowercase common name (`maize`, `human`);
crossing structure goes under `design` (`scheme`, `founder_count`,
`selfing_generations`, `backcross_generations`, `population_size`); typed
catalog parameters go under `parameters`. Only explicitly requested secondary
analyses are listed.

## Scoring

The agent exposes two partial views (legacy keyword spec and catalog
StudySpec); `.simitall_eval_view()` merges them into StudySpec vocabulary.

| Metric | Meaning |
|---|---|
| `reading_recall` | Share of expected facts (organism, design, parameters, region) recovered |
| `goal_and_source` | goal (simulate/analyze/both) and data source inferred from the wording |
| `primary_analysis` | Correct main analysis |
| `secondary_exact` | Exactly the requested extra analyses, nothing more |
| `must_not_violations` | Unrequested analyses added (e.g. RNA-seq on a GWAS request) |
| `stops_when_it_should` | Unsupported / unsafe requests stop at a plan |
| `examples_with_invented_values` | A must-ask field was filled without disclosure |
| `examples_defaulted_instead_of_asking` | A must-ask field was filled from a disclosed default |

`split = test` rows are frozen: do not tune rules or models while looking at them.
