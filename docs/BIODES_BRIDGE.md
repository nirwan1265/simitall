# Biodes -> Simitall -> Bioana: what Simitall reads from Biodes

Repos: Biodes = github.com/nirwan1265/BioDesignGPT (R package `studyforge`);
Simitall + Bioana = github.com/nirwan1265/simitall. Hand-off contract: Biodes
`docs/SUITE_INTEGRATION.md` (handoff_version 0.4).

## Two ways Simitall gets a design from Biodes

### 1. From a design folder: `simitall_from_biodes(dir)`

Biodes writes a folder with `write_design_package(book, spec, dir)`:

| file | used for |
|---|---|
| `study_spec.json` | `handoff_version`, `design.chosen`, `design.randomization_seed`, `analysis.model`, `analysis.terms`, `power.*`, `primary_outcome`, `subsamples_per_unit`, `assays` |
| `field_layout.csv` | one row per experimental unit; columns named after the model terms |
| `sample_metadata_template.tsv` | same rows plus `sample_id` |

Simitall simulates the outcome from `analysis.model` and writes
`sample_metadata_template.tsv` (filled), `sim_truth.json` and, for RNA-seq,
`counts.tsv`. `bioana_fit_biodes(dir)` then fits the same model and writes
`results.tsv`, `fit_summary.json`, `recovery.tsv` (`de_results.tsv` for counts).

```r
simitall_from_biodes("~/biodes_runs/my_design")
bioana_fit_biodes("~/biodes_runs/my_design")
```

### 2. From a question: `simitall_from_question(question, answers = list(...))`

Calls `studyforge::studyforge_ask(question, fix = answers)` and uses only:

- `$spec`: the study details (factors, reps, locations, outcome, power values);
- `$candidates[1]`: the design id (or `design =` to choose another candidate);
- `$questions`, `$validation$missing`: shown as open questions.

Simitall then **builds the layout itself from `$spec`**
(`simitall_layout_from_spec(spec, design)`), writes the folder with Biodes'
`write_design_package()` (so `analysis.model` comes from Biodes' catalog), and
runs Simitall and Bioana. Biodes' `$recipe` (template layout code) is not used.

```r
r <- simitall_from_question("Greenhouse drought on 6 sorghum genotypes ...",
  answers = list(factors = "genotype:6, water:2", design.reps = 4,
                 power.minimum_detectable_difference = 5, power.residual_sd = 6))
ps <- simitall_power_sweep(r)
```

When details are missing, `simitall_from_question()` prints suggestions for
each one: what to measure and possible experimental units from Biodes'
registries (`traits`, `experimental_units`, each cited), Biodes' setup
hazards, and Simitall's `inst/suite/suggestions.tsv` (what to compare, omics
outcomes, effect size / SD), matched to the question, system, organism and
assays. `use_suggestions = TRUE` fills every missing item with its first
suggestion, prints what it filled, and runs; your own `answers` always win.

Figures: `biopaint_biodes(dir)` (run automatically by `simitall_from_question()`)
writes `figures/layout.png`, `outcome.png`, `recovery.png`, `power.png` (after a
power sweep) and `volcano.png` (count assays), using whatever the folder holds.

`simitall_power_sweep(dir_or_run)` rebuilds the layout for several values of
`design.reps`, simulates and fits each many times, and writes
`power_sweep.tsv` / `power_sweep.json` (power per model term).

## Spec fields Simitall reads (Biodes should fill these)

| field | meaning | if missing |
|---|---|---|
| `factors[].name`, `factors[].n_levels` | study factors and level counts | asks; level counts are never invented |
| `factors[].levels` | level names | `<name>_1`, `<name>_2`, ... |
| `factors[].applied_to` | `main_plot`, `whole_plot`, `distance_from_line`, `strip` put the factor on main plots | second factor is the main plot |
| factor named time / day / week ... | measured within each unit | see `inst/suite/factor_roles.tsv` |
| factor named check / checks | check entries (augmented, diagonal checks) | asks, for designs that need checks |
| factor named location / site / env | environments | `design.locations` |
| `design.reps` | replicates, blocks or benches | default per design, listed as an assumption |
| `design.blocks` | blocks for designs with unreplicated entries + checks | `design.reps` if > 1, else default |
| `design.locations` | environments for multi-location designs | asks |
| `design.block_size` | incomplete-block size (alpha, lattices) | divisor of the entry count nearest its square root |
| `design.n_cols` | field grid width (row-column, spatial) | square-ish grid |
| `power.minimum_detectable_difference` | treatment effect simulated | `simulation_defaults.json` |
| `power.residual_sd` | residual SD | `simulation_defaults.json` |
| `power.interaction_effect` | interaction between treatment factors (G x W, treatment x time) | ratio x main effect, `simulation_defaults.json` |
| `power.attrition` | share of units lost | 0 |
| `primary_outcome.name`, `.baseline_mean` | outcome name and intercept | intercept 100 |
| `subsamples_per_unit` | leaves per plant etc. (not replicates) | 1 |
| `assays` | `rnaseq` and similar give `counts.tsv` | outcome only |

`analysis.model` must use Biodes' generic terms (`treatment treatment2 block
rep row col unit subject time env period sequence check`). Simitall's layouts
name their columns with these terms and Biodes' `write_design_package()` lists
them in `analysis.terms.layout_columns`. An explicit `analysis.terms.columns`
map (`{"treatment": "main"}`) is also read.

## Rule files in Simitall (edit these, not the code)

| file | what it says |
|---|---|
| `inst/suite/layouts.tsv` | one row per Biodes catalog design (74): treatment structure, blocking, within-unit structure, per environment, units/subjects, observations per unit, default reps, note |
| `inst/suite/factor_roles.tsv` | which spec factor name is time, check, environment, period |
| `inst/suite/handoff.json` | supported handoff versions, count assays, unsupported model tokens, check names, main-plot `applied_to` values |
| `inst/suite/model_terms.json` | model term -> role and candidate columns (old names A/B, main/sub, genotype still work) |
| `inst/suite/simulation_defaults.json` | effect, residual SD, interaction ratio, block/random/subsample SDs, count settings |
| `inst/suite/power_defaults.json` | reps tried, simulations, alpha, target power, time limit |
| `inst/suite/suggestions.tsv` | suggestions for missing factors, unit, outcome, effect/SD: field, context regex, suggestion, answer, source |

When Biodes adds a design to `design_catalog.json`, add a row to
`layouts.tsv` with the same id. A design without a row (or with
`treatments = none`) is reported as having no layout rule, never replaced by
another design.

## What stops a run (and what Biodes can do about it)

- **Asks for details**: a factor's level count, check entries, environments, or
  at least 2 reps for a blocked design. Biodes: read these from the question or
  ask them.
- **Layout does not match the plan**: the study has more factors than the
  chosen design's model has treatment terms (e.g. genotype x water x time fitted
  as `y ~ treatment + block`). Biodes: choose candidates whose model carries
  every study factor, or say why one is left out.
- **No layout rule**: `mapping_population` (needs lines per family, families)
  and `split_split_plot` (catalog model has two treatment terms; needs three).
