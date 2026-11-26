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
`counts.tsv` and the other data, and lists them in `data_manifest.json`.
`bioana::bioana_analyze(dir)` (github.com/nirwan1265/bioana) then analyses
every data set with the same model and writes `bioana/<id>/results.tsv`,
`bioana/outcome/recovery.tsv` and `bioana/bioana_summary.json`.

```r
simitall_from_biodes("~/biodes_runs/my_design")
bioana::bioana_analyze("~/biodes_runs/my_design")
```

### 2. From a question: `simitall_from_question(question, answers = list(...))`

Calls `studyforge::studyforge_ask(question, fix = answers)` and uses only:

- `$spec`: the study details (factors, reps, locations, outcome, power values);
- `$candidates[1]`: the design id (or `design =` to choose another candidate);
- `$questions`, `$validation$missing`: shown as open questions.

The layout comes from **Biodes `studyforge::layout_from_spec(spec)`** (every
factor kept as its own column); Simitall writes the folder with Biodes'
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

## What Simitall simulates for a design (routes)

`simitall_routes(dir)` lists, without simulating, what `simitall_from_biodes(dir)`
will run and the value and source of every parameter. Routes are rows of
`inst/suite/simulators.tsv`, matched on spec fields (Biodes' canonical ids,
lower case, `-` read as `_`):

| layer | spec field | values | simulator |
|---|---|---|---|
| outcome | `primary_outcome.type` | continuous (or missing), binary, count, ordinal | built in: Gaussian from `analysis.model`, liability threshold, negative binomial, ordinal thresholds |
| genetics | `assays` | gwas, gbs, wgs, amplicon, targeted_capture | diversity panel of the genotype factor's levels |
| genetics | `population.mating` | ril, dh, nam, magic, backcross, selfing, sib_mating, biparental, f2 | `simulate_breeding()` cross (`mating_schemes.tsv`) |
| genetics | `assays` | qtl_mapping, bsa_seq, pool_seq | cross (scheme from `population.mating`, else F2) |
| assay | `assays` | rnaseq (+ synonyms) | `simulate_rnaseq_experiment()` |
| assay | `assays` | eqtl | `simulate_rnaseq_from_gwas()`, each unit carries its line's genotypes |
| assay | `assays` | single_cell | `simulate_scrnaseq_experiment()`, pseudobulk per unit |
| assay | `assays` | chip_seq, atac_seq | `simulate_chipseq()`, `simulate_atacseq()` on a synthetic genome |
| assay | `assays` | small_rnaseq | `simulate_counts_from_template()` with a same-species miRNA template (human, mouse; Arabidopsis, maize, sorghum, barley, tomato from PmiREN), else native NB |
| assay | `assays` | amplicon (16S/ITS), metagenomics | `simulate_counts_from_template()` with a host-kingdom 16S template (plant root, human HMP) |
| assay | `assays` | methylation, wgbs, rrbs, bisulfite_seq, em_seq | `simulate_methylation()`: methylated + covering reads per cytosine with known DMRs; same-kingdom WGBS template (mouse ENCODE for animals), else native CG/CHG/CHH regions (plants) |
| assay | `assays` | methylation_array | `simulate_methylation(engine = "array")`: beta values learned from TCGA 450k, known differentially methylated positions |
| assay | `assays` | m6a_seq, merip_seq, merip | `simulate_merip()`: IP and input counts per window (m6A peaks + background) with changed peaks and expression changes; expression from a same-species RNA-seq template (human airway) else native |
| - | `assays` | medip | planned |
| treatment | `assays` | crispr, overexpression, complementation, ... | a way to make the lines: a treatment, not data |

Factors the model does not hold (`analysis.terms.uncovered_factors`, e.g. tissue measured within each plant) get main effects, and the rows of one unit share a unit effect (`sim_truth.json$within_unit_sd`). For assays, the condition is the combination of the first non-genotype treatment and these factors (e.g. `water_1.tissue_2`).

Genetics needs a genotype term: a treatment term whose `analysis.terms` meaning
(or factor name) matches `genotype_pattern` in `inst/suite/genetics.json`
(genotype, line, accession, variety, ril, ...). Its levels become the lines;
their genetic values (phenotype model from `phenotype_models.tsv`: additive,
additive_dominance, epistatic, gxe, full) replace the even-spread genotype
effects, scaled to plot-level heritability h2 = Vg / (Vg + `power.residual_sd`^2).

Spec fields read by genetics: `population.mating`, `population.generations`
(`F6` = 5 selfing generations), `population.parents` (a list of names gives the
founder count), `system` (inbred or outbred lines, `genetics.json`). Fields
Biodes does not have yet (heritability, phenotype model, QTN counts) come from
`simulator_parameters.tsv` and are listed in `sim_truth.json$assumptions`;
`params = list(suite_genetics = list(heritability = 0.7))` sets them.

Outputs: `genotypes.vcf` (lines), `genetic_truth/` (`genetic_values.tsv`,
`qtn_truth.tsv`, `population_structure.tsv` or `founder_ancestry.tsv`),
`counts.tsv` (eQTL > RNA-seq > single-cell; other count assays write
`<assay>/unit_counts.tsv`), `eqtl/`, `single_cell/`, `chipseq/`, `atacseq/`,
and `sim_truth.json` with `routes`, `genetics`, `assays`, `not_simulated`.

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
| `assays` | routes in `inst/suite/simulators.tsv` (see above) | outcome only |
| `analysis.primary_goal`, `goal_terms`, `goal_effect[goal]` (`size`, `shape`, `peak`, `half_time`, `scope`) | which truth to plant (`inst/suite/goals.tsv`), written to `sim_truth.json` `goal_truth` | `treatment_effect`; size = the treatment effect |
| `factors[].values`, `.unit`, `.scale`, `.log_base`, `.log_offset` | sampling times, dose values (trend on log_base(x + log_offset) for `scale = log`) | numbers in the level labels; offset = smallest positive level when a level is 0 |
| `primary_outcome.unit`, `.baseline_mean`, `power.residual_sd`, `power.cv` | the outcome's scale (`sim_truth.json` `outcome_scale`); SD from CV % x baseline mean | a unit with no SD: standardized (SD 1, mean 0, effects in SDs); no unit: SD 10, mean 100 |
| `analysis.terms.repeated_factors`, `.time` | the occasion a unit is measured at; correlation over occasions from `suite_outcome` `within_subject` (none / ar1 / random_slope), in `sim_truth.json` `within_subject` | none: shared unit effect only |
| `primary_outcome.type = survival`, `.event_codes`; `power.target_hazard_ratio`, `.loss_fraction`; `analysis.survival.landmarks`, `.rmst_tau`, `.life_table` | survival times and event codes (`suite_outcome` `survival_model`), truth per group in `outcome_model` | Weibull PH, HR 0.6, median 30, 20% administratively censored |
| `primary_outcome.type = nominal`, `.categories`; `analysis.expected_ratio` | classes from a multinomial logit, or the segregation ratio (`segregation_distortion`) | 3 equal classes |
| `secondary_outcomes`, `analysis.correlation` | extra trait columns with a within-design residual correlation (`trait_residual_cor`) | none |
| `analysis.agreement` (`raters`, `scale`, `levels`) | one column per rater: true value + bias + noise, or a misclassification matrix | none |
| `analysis.covariates` [{term, name, column, unit, transform, mean, sd}], terms `covariate`... | covariate columns per unit (`suite_covariates`), slopes planted; truth `goal_truth.true_slopes` / `true_log_odds_ratios`, `regression_truth` | mean 50, SD 10, correlation 0.3; slope = treatment effect per SD |
| `goal_effect[goal].scale` | identity: the size is on the outcome; log2: it sizes the goal's count features | outcome |
| `goal_effect.time_course.group`, `.group_factor` | with `scope = within_group`, the one level that changes over time | the non-reference levels |

`goal_truth.true_effect` is on the outcome scale and follows Biodes' registry
definition: interaction = difference of differences between the last and
reference levels of both goal terms; time_course = largest change between
time points (`time_pattern` per level); trend = change per unit of the factor
(per log10 unit of x + log_offset on a log scale; with no stated size, the
treatment effect over the factor's range). `model_terms` are the model terms that test the goal; the
power sweep chooses replicates for them. Count assays carry the goal in 10% of
features (`goal_truth.features`, `truth/<assay>_goal.tsv`: `feature_id, term,
effect, change`, log2).

`analysis.model` must use Biodes' generic terms (`treatment treatment2 block
rep row col unit subject time env period sequence check`). Simitall's layouts
name their columns with these terms and Biodes' `write_design_package()` lists
them in `analysis.terms.layout_columns`. An explicit `analysis.terms.columns`
map (`{"treatment": "main"}`) is also read.

## Rule files in Simitall (edit these, not the code)

| file | what it says |
|---|---|
| `inst/suite/factor_roles.tsv` | which spec factor name is time, check, environment, period |
| `inst/suite/handoff.json` | supported handoff versions, unsupported model tokens, check names, main-plot `applied_to` values |
| `inst/suite/model_terms.json` | model term -> role and candidate columns (old names A/B, main/sub, genotype still work) |
| `inst/suite/simulation_defaults.json` | effect, residual SD, interaction ratio, block/random/subsample SDs |
| `inst/suite/simulators.tsv` | spec field value -> simulator, status |
| `inst/suite/simulator_parameters.tsv` | each simulator's arguments, defaults, spec field |
| `inst/suite/genetics.json`, `mating_schemes.tsv`, `phenotype_models.tsv` | genotype term, inbred systems, cross schemes, genetic architectures |
| `inst/suite/power_defaults.json` | reps tried, simulations, alpha, target power, time limit |
| `inst/suite/suggestions.tsv` | suggestions for a missing unit, outcome, effect/SD (never factors: what to compare is the user's question) |

Layouts are Biodes' job (`layout_from_spec()`); when it stops with
`missing: ...`, Simitall shows that as a question for the user.

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

## "How many per group?"

When the spec has no `design.reps` (or `power.solve_for` is `"reps"`),
`simitall_from_question()` runs `simitall_power_sweep()` over `power.reps_range`
(default `power_defaults.json`: 2, 3, 4, 5, 6, 8) and lays the design out with
the smallest count where every treatment main effect (and, for count assays, the
genes with a true 2-fold change) reaches power 0.8. The sweep is kept in
`power_sweep.tsv`, the choice in `power_choice.json`.

## Hand-off to Bioana: `data_manifest.json`

After simulating, `simitall_from_biodes()` writes `data_manifest.json`: one
entry per data set with its `data_type` (a row of Bioana's
`inst/rules/analyses.tsv`), its files and its truth in Bioana's format under
`truth/` (features: `feature_id, condition, effect`; regions: `region_id, chr,
start, end, condition, effect`). `bioana::bioana_analyze(dir)` then runs the
matching analysis for each and scores it against the truth
(github.com/nirwan1265/bioana). Code: `R/suite_manifest.R`.

| Simitall route | data_type | files |
|---|---|---|
| outcome (any type) | outcome | sample_metadata_template.tsv (`primary_outcome`, `event` for survival) |
| assay_rnaseq, assay_eqtl, assay_single_cell, assay_small_rnaseq, assay_microbiome | feature_counts | counts.tsv or `<assay>/unit_counts.tsv` |
| assay_chip_seq | peak_counts_long | chipseq/chip.peak_counts.tsv + chip.sample_metadata.tsv |
| assay_atac_seq | read_positions | atacseq/atac.read_positions.tsv + library_to_unit.tsv |
| assay_methylation (sequencing / array) | bisulfite_counts / methylation_beta | methylation/*.tsv.gz |
| assay_m6a | ip_input_counts | m6a/m6a.ip.tsv.gz, m6a.input.tsv.gz |
| genetics | genotypes | genotypes.vcf (+ `genotype_term`, the layout column of the lines) |

A new simulator adds its manifest entry there; Bioana adds the analysis.
