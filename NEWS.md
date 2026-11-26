# simitall (development)

- Biodes 0.6.7: a template already expanded to one row per subsample
  (`subsample` column) keeps its rows and sample ids; units come from the
  layout. Survival with subsamples (seeds in a tray cell,
  `analysis.survival.cluster`) gives each unit a shared log-normal frailty
  (`cluster_log_hazard_sd`, 0.3). A yes/no trait with
  `analysis.expected_ratio` draws class labels from the ratio and is written
  as a two-class nominal outcome. Raters in long format
  (`analysis.agreement` `rater_column`, `score_column`, `unit_column`): each
  unit's true score is read once per rater into the score column.
- Categorical outcomes, traits and raters (Bioana
  docs/contracts/categorical_correlation_agreement.md):
  - `primary_outcome.type = "nominal"` (new route `outcome_nominal`):
    multinomial logit on the design, classes from
    `primary_outcome.categories`; truth `outcome_model.group_probability`.
    With `analysis.expected_ratio` (e.g. 3:1) classes follow the ratio;
    `segregation_distortion` shifts it (truth `true_ratio`, `distorted`).
  - Binary outcomes on the logit scale with the unit's effect when units are
    measured more than once (`binary_model = "auto"`, or `"logit"`); truth
    the conditional odds ratios (`odds_ratio`, default 3) and the unit
    effect SD on the logit scale. Single measurements keep the liability
    threshold.
  - `secondary_outcomes` [{name, column, unit, baseline_mean, residual_sd}]:
    extra trait columns sharing `secondary_effect_share` of the design
    effects, residuals correlated `trait_residual_cor`; truth
    `secondary_traits.within_design_correlation`.
  - `analysis.agreement` {raters, scale}: continuous ratings = true value +
    rater bias + noise (truth ICC2, ICC3); nominal/ordinal ratings through a
    misclassification matrix (`rater_accuracy`; truth expected kappa and, for
    ordinal, quadratic weighted kappa). Long format is not simulated yet.
- Survival strata that are not model terms (Biodes 01d41ef: stated blocks
  under 10 go to `analysis.survival.strata` with no `+ block`) get their own
  baseline hazard: a log-hazard shift per stratum with SD
  `strata_log_hazard_sd` (0.3), recorded in `outcome_model.strata`.
- Survival with a stated follow-up and no stated median: the reference
  group's median is set so that `1 - power.event_probability` (else
  `censored_fraction`) of it is still without the event at the end of
  follow-up (`outcome_model.median_time_source`), instead of a fixed 30 that
  could exceed the follow-up. Interval scoring: an event is seen at the
  next check; one after the last check before the follow-up ends is censored
  there; a lost unit is censored at the check before it was lost. Group
  medians beyond the follow-up are still reported.
- Biodes 0.6.6 survival fields: `analysis.survival.logrank_weight`
  `fh_0_1` (stated delayed effect) simulates `delayed_effect`, `fh_1_0`
  the new `early_effect` (HR until `change_time`, then 1), unless `params`
  name a model; the follow-up length is read from `primary_outcome.censoring`
  ("administrative at 21 days"); `analysis.survival.strata` are recorded with
  each stratum's log-hazard shift. Survival groups in the truth are in level
  order (reference first); the RMST-difference goal names the two groups.
- Survival outcomes (Bioana docs/contracts/survival.md): `suite_outcome`
  `survival_model` = `weibull_ph` (new default; proportional hazards, shape
  `weibull_shape` 1.5), `exponential`, `lognormal_aft` (the old model),
  `delayed_effect` (HR 1 until `change_time`) or `crossing` (HR, then 1/HR).
  The full treatment effect has hazard ratio `hazard_ratio` (0.6;
  `power.target_hazard_ratio`, or the survival goal's size on
  `hazard_ratio`/`time_ratio`). Censoring is administrative at `follow_up`
  plus random loss (`loss_fraction`, `power.loss_fraction`), counted apart.
  `primary_outcome.event_codes` 2+ add competing events (event column 0/1/2,
  Weibull cause-specific hazards with no treatment effect). Blocks and random
  terms act as a log-normal frailty (variance in truth).
  `analysis.survival.life_table` rounds times up to the scoring grid.
  `sim_truth.json` `outcome_model` gives per group the median, S(t) at
  `analysis.survival.landmarks`, RMST at `rmst_tau`, CIF at the landmarks
  (competing), the hazard and time ratios per level, `proportional_hazards`,
  and the censoring counts; `goal_truth` for the survival goal reports the
  HR (and its log), time ratio or RMST difference.
- Repeated measures can break sphericity: `suite_outcome` `within_subject`
  = `ar1` (a unit's residuals correlate `within_rho`^lag over the occasions,
  default rho 0.7) or `random_slope` (each unit its own slope over time, SD
  `slope_sd_ratio` x residual SD over the full range). Default `none` keeps
  the shared unit effect only (compound symmetry), with the same random draws
  as before. `sim_truth.json` `within_subject` records the structure, rho or
  slope SD, the occasion column and the unit columns. Layout columns unique on
  every row (a `sample_id`) are ignored when finding the units; when `ar1` or
  `random_slope` cannot apply (no repeated factor, no unit seen twice) it
  warns and records `{"structure": "independent", "requested", "reason"}`.
- Hidden factors in bulk RNA-seq are off by default in the suite
  (`simulate_rnaseq_experiment` `n_latent = 0`): with random scores per unit,
  two hidden factors in a 16-sample design could line up with the treatment
  and shift many genes with no planted effect (Bioana's `called_not_planted`).
  `n_latent` and `latent_effect_sd` are now listed in
  `simulator_parameters.tsv` for bulk and single-cell; when used, scores per
  unit and loadings per gene are written (`rnaseq/latent_scores.tsv`,
  `rnaseq/latent_loadings.tsv`) and recorded in `sim_truth.json` under the
  assay's `counts.hidden_factors`.
- Biodes 0.6.4 outcome scale: `primary_outcome.unit`, `baseline_mean`,
  `power.residual_sd` and `power.cv` (%, times the baseline mean when no SD)
  set the outcome's scale. An outcome with a unit but no SD (or one Biodes
  asked for) is simulated standardized (residual SD 1, mean 0, effects in
  SDs) instead of SD 10 / mean 100; `sim_truth.json` `outcome_scale` says
  which, with the stated effect kept as text. Count outcomes keep their
  defaults.
- Biodes 0.6.3: `goal_effect.trend.per_unit` (e.g. "kg N/ha", "log10 uM")
  is written to `goal_truth.per_unit` and its effect scale; a within-group
  time course changes only the `group` level Biodes names. The power sweep
  reads Bioana's goal test on the outcome (`goal_recovery$detected`) as the
  goal row, and `goal_recovery$tpr` for count features.
- Biodes 0.6.2: trend on log-scaled factors uses log_base(x + log_offset),
  the predictor Bioana fits (0/1/10/100 uM -> log10 of 1/2/11/101);
  `goal_effect.scale` decides what the size sizes (identity: the outcome;
  log2: every goal feature in the count assays, a trend size per log10 unit);
  no stated size = the treatment effect over the factor's range. Time course
  `within_group` changes only `group` (a level of `group_factor`); a stated
  peak with no shape is transient, a half-time alone saturating. Assay-only
  designs report the goal per feature and keep the outcome-scale values under
  `goal_truth.outcome_truth`. The design reader reads factor `values`,
  `log_base` and `log_offset` and turns level labels into values for numeric
  model terms (same as Bioana's `read_design.R`).
- Analysis goals (Biodes schema 0.6.0, `analysis.primary_goal`): the truth
  planted follows the primary goal and is written to `sim_truth.json` ->
  `goal_truth` (`goal`, `terms`, `columns`, `test`, `definition`, `size`,
  `true_effect`, `effect_scale`, `model_terms`, `status`). Interaction: the
  difference of differences of the extreme levels = the goal size.
  Time course: shape (monotone, saturating, transient, delayed; peak and
  half-time converted to the time factor's unit) on the real sampling times,
  largest change = the goal size; `within_group` changes only the
  non-reference groups. Trend: linear in the factor's values (log10 for
  `scale = log`, a zero dose one decade below the lowest). Count assays
  (RNA-seq, small RNA, microbiome, single-cell pseudobulk, eQTL) get the
  same pattern in 10% of features by binomial thinning, truth in
  `truth/<assay>_goal.tsv` and the manifest's `truth.goal`. Genetics,
  survival, ranking, variance components and microbiome goals point to the
  truth already written; registry goals with `simulated: false` say
  `not_simulated`. No goals = `treatment_effect`, as before. What each goal
  plants: `inst/suite/goals.tsv`.
- The power sweep marks the goal's terms (`goal`, `goal_term` columns) and
  the replicate choice in `simitall_from_question()` uses them; it also reads
  Bioana's `goal_recovery$tpr` for goal features when Bioana reports it.
- Reads `handoff_version` 0.6.
- Paired designs (Biodes `analysis.terms.paired`: the repeated occasion is the
  only factor and each block is one animal) get no extra unit effect; the
  block already is the animal.
- Biodes 0.5.2 `repeated_factors` (the same units at every time / stage /
  age): a repeated factor the model does not hold is simulated like tissue
  (an effect per level, time spaced by sampling time) and the rows of one unit
  share a unit effect; a repeated factor the model holds also gets the shared
  unit effect when the model has no random term for the unit. `[]`
  (destructive sampling) keeps the occasion factor randomized to units.
- Time factors follow the real sampling times: when the levels carry numbers
  (`0 h`, `24 h`, `day 10`) or the spec gives `factors[].values`, the time
  effect and its interactions grow with time (0, 6, 24, 48 h are not evenly
  spaced). Generic labels (`time_1`, ...) stay evenly spaced. `stage` and
  other named states are categorical.
- Factor levels follow the spec's `factors.levels` (reference first), then
  natural order for levels it does not list, as in Bioana. Before, levels were
  always in natural order, so `trt` became the reference before `untrt` and the
  planted effects had the opposite sign of what the spec meant. The outcome,
  genetics, assays, truth files and the outcome figure all use the new order.
- `simitall_power_sweep()` and `simitall_from_question()` take
  `allow_fallback` (default `FALSE`, as in Bioana) and warn when Bioana
  refuses an analysis because its citable package is not installed.
- All analysis moved to Bioana (github.com/nirwan1265/bioana): GWAS, genomic
  prediction and selection, peak calling and differential peaks, single-cell
  clustering and tests, eQTL benchmarks, RNA-seq results summary, GO
  annotation (with the gene-coordinate catalogues), assembly evaluation and
  `bioana_fit_biodes()`. Simitall only simulates; `simitall_power_sweep()`
  and `simitall_from_question()` call `bioana::bioana_analyze()`, and
  `biopaint_biodes()` reads Bioana's output.
- `data_manifest.json`: every run lists its data sets (data type, files, truth
  in Bioana's format under `truth/`) for Bioana (github.com/nirwan1265/bioana).
- Biodes 0.5.1: routing reads `assay_kinds` (DNase/FAIRE -> ATAC, DAP-seq/SELEX
  -> ChIP; amplicon genotyping, repertoires, CRISPR screens, reporter assays
  planned), breeding notation sets backcross and selfing generations, and
  survival outcomes simulate (log-normal time with censoring, `event` column).
- Backcross / NIL lines (rare donor allele) no longer stop the genetics layer:
  QTN are drawn among the commonest polymorphic markers, with a note.
- `simulate_merip()`: m6A MeRIP-seq as IP and input counts per window (peaks
  and background), with known peaks, changed peaks and expression changes that
  move IP and input together. Biodes designs with `m6a_seq` / `merip` simulate it.
- "How many per group?": when a question has no replicate count,
  `simitall_from_question()` chooses `design.reps` with `simitall_power_sweep()`
  (reads `power.reps_range` and `power.solve_for`). The sweep now also reports,
  for count assays, the share of genes with a true 2-fold change that are found.
- Bioana's count test uses empirical-Bayes moderated t-tests (as limma eBayes),
  so designs with few replicates are tested fairly.
- A design whose genetics assay has no genotype factor (BSA bulks, populations)
  no longer stops: the rest is simulated and genetics is listed as
  `needs_genotype_factor` in `sim_truth.json$not_simulated`.
- `simulate_methylation()`: bisulfite-sequencing methylation (WGBS, RRBS,
  EM-seq) as methylated and covering read counts per cytosine, with known DMRs
  (`dmr_truth.tsv`) and per-site truth. Three engines: a real template of the
  same kingdom (mouse ENCODE WGBS), a native engine with CG/CHG/CHH chromatin
  states for plants (`inst/suite/methylation_presets.json`), and beta values
  for 450k/EPIC arrays (TCGA template). Biodes designs with methylation assays
  now simulate it; `simulator_parameters.tsv` has an `ask` column whose
  questions `simitall_routes()` prints.
- Real-data templates for count simulation: `simitall_templates()`,
  `simitall_template()`, `simulate_counts_from_template()` (binomial thinning
  of real samples, or NB learned from them). Small RNA-seq (human and mouse
  miRNA atlas) and 16S/ITS microbiome (plant root, human HMP) now simulate
  from them for Biodes designs.
- Simitall now creates data for a Biodes design from the spec's fields:
  `inst/suite/simulators.tsv` routes `primary_outcome.type` (continuous,
  binary, count, ordinal), genetics (`assays` gwas/qtl_mapping/eqtl,
  `population.mating`) and omics assays (RNA-seq, eQTL, single-cell, ChIP-seq,
  ATAC-seq) to simulators, with parameters in
  `inst/suite/simulator_parameters.tsv`. New: `simitall_routes()`,
  `simitall_simulator_parameters()`, `params` in `simitall_from_biodes()`.
- Removed the prompt readers (`simitall_ask()`, `launch_simitall_agent()`,
  `search_simitall_knowledge()`, `simitall_agent_tools()`, the keyword spec,
  lexicon, prompt model and capability catalog): Biodes reads the user's
  question. They and the Shiny app are kept in `old_model/` (not built or
  loaded); see `old_model/README.md` to run them.

# simitall 0.1.0

- Added a retrieval-grounded, read-only agent interface through
  `search_simitall_knowledge()`, `simitall_ask()`, and
  `simitall_agent_tools()`. Added an optional local Shiny chat interface that
  uses a locally stored OpenAI API key and never executes simulations.

- Added mixed-model GWAS analysis, Manhattan/QQ plotting, and causal-truth
  benchmarking through `rrBLUP`.
- Added GBLUP, RR-BLUP, Bayesian, and random-forest genomic prediction,
  cross-validation, genomic kinship, diversity-aware parent selection,
  crossing-plan design, and one-round genomic-selection orchestration.
- GWAS cohort simulation now writes causal-variant truth and per-sample true
  breeding values.

* Added package-backed breeding phenotype wrappers: `simstudy` provides
  ordinal, Poisson, negative-binomial, and shared family effects, while
  `pedtricks` provides direct, maternal, and paternal pedigree effects with
  component-level truth output.
* Added multi-chromosome breeding with `founder|chromosome` FASTA panels,
  chromosome-specific recombination maps, independent assortment,
  chromosome-aware VCF output, founder-ancestry tracts, and breakpoint truth.
  Added a bundled three-chromosome founder panel and generator controls for
  chromosome counts and lengths.
* Expanded the breeding tutorial with runnable F1, F2, F2:3-style,
  F2-derived S3, backcross, RIL-SSD, RIL-sib, DH, NIL, NAM, MAGIC, and SimuPOP
  examples. Corrected F2 and MAGIC generation expansion and made empty VCF
  truth sets safe for monomorphic final populations.
* Added one-command, profile-based bootstrap installers for Apple Silicon and
  Intel macOS, Linux/WSL2, and native Windows PowerShell. Profiles install
  minimal, population, omics, sequencing, or full dependency stacks.
* Added `check_simitall_dependencies()` and
  `install_simitall_dependencies()` for dependency audits and explicit
  profile installation from R. Native Windows sequencing/full installations
  delegate to WSL2 for Bioconda-only tools.
* Added annotation-aware TF, histone-mark, and nucleosome ChIP-seq simulation
  with GFF3/BED targets, matched input controls, biological and technical
  replicates, differential binding, PCR duplicates, sequencing errors, FASTQ,
  BED, peak-count, QC, and truth outputs.
* Added native and optional ChIPsim read-density backends plus a six-panel
  ChIP-seq publication figure and reproducible example runner.
* Added standalone bulk and single-cell RNA-seq experiment simulators that do
  not require GWAS genotypes, with explicit technical, biological, mixed, and
  custom replicate designs.
* Added separate biological- and library-level effects, standalone condition
  and batch truth, donor-level and library-level pseudobulk, and reusable
  single-cell publication plotting for one-sample experiments.
* Added donor-aware single-cell RNA-seq simulation from GWAS genotypes with
  shared and cell-type-specific eQTLs, marker programs, pseudotime, condition,
  batch, donor, latent, dropout, ambient-RNA, and doublet effects.
* Added native and optional Splatter technical backends, sparse 10x-style
  output, optional `SingleCellExperiment` export, donor-cell-type pseudobulk,
  cell-type eQTL benchmarking, and six-panel publication figures.

- Established the project as an installable R package.
- Added R-first genome, annotation, read, assembly, GWAS, phenotype, breeding,
  and SimuPOP interfaces.
- Added GWAS-linked RNA-seq simulation with cis/trans eQTLs, condition and
  batch effects, genotype-by-condition interactions, latent confounders,
  allele-specific expression, optional Rsubread FASTQ generation, and eQTL
  recovery benchmarks.
- Added RNA-seq/eQTL result summaries, six-panel publication figures, exported
  panel source data, and a reproducible end-to-end paper figure runner.
- Added lightweight package examples and on-demand reference download helpers.
