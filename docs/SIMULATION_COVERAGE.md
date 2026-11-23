# Can Simitall simulate Biodes' questions? (1,000 paper questions)

Questions: Biodes `training/paper_questions/splits/test.jsonl` (764) + the first
236 of `dev.jsonl`, one question per paper, stratified by species x technique.
Per-question results: [`data-raw/coverage/coverage_1000.tsv`](../data-raw/coverage/coverage_1000.tsv)
(`verdict`, `needs_status`, `blocker`, `reps_given`, and the rules-reader result).

Each question is run on its **reviewed spec** (`expected_spec`, as if Biodes read
it perfectly) -> Biodes design choice -> `studyforge::layout_from_spec()` ->
`write_design_package()` -> `simitall_from_biodes()`. 361 questions give no
replicate count (mostly "how many per group?"); they were run with 3 reps
(`reps_given = FALSE`, `REPS_IF_MISSING`), since the power sweep tries several.

## Current result (Biodes 568c9b5, Simitall main with methylation, 2026-11-22)

| verdict | n | before methylation | before templates | first run |
|---|---|---|---|---|
| simulated_all: every data type it needs | 661 | 594 | 437 | 277 |
| simulated_outcome_question: needs only a measured trait | 50 | 50 | 50 | 37 |
| simulated_part: some data types have no simulator | 60 | 95 | 89 | 59 |
| no_simulator_for_data_type | 117 | 149 | 313 | 189 |
| not_run: no level count for a factor | 54 | 54 | 58 | 33 |
| not_run: spec has no factors | 46 | 46 | 46 | 45 |
| not_run: other (genetics without genotype factor 3, Biodes error 2, ...) | 6 | 6 | 7 | 42 |

**711 / 1000 (71%) get everything they need; 60 more get part of it.**
Methylation (WGBS, RRBS, EM-seq, bisulfite, 450k/EPIC arrays) now simulates with
`simulate_methylation()`, which moved 67 questions; the 18 MeDIP questions still
wait for an enrichment simulator. Small RNA-seq and 16S/ITS microbiome simulate
from real templates (`simitall_templates()`).

## Second 1000 questions (2026-11-22)

Questions: the rest of `dev.jsonl` (92), `dev2` (155), `test_draw4` (206),
`test_draw5` (159), `test_heldout` (87), `test_newkinds` (100) and the first 201
of `train.jsonl` (`QUESTION_SET=second`). Per-question results:
[`data-raw/coverage/coverage_second_1000.tsv`](../data-raw/coverage/coverage_second_1000.tsv)
(from `summarize_coverage.R`). 466 had no replicate count and were run with 3.

| verdict | n |
|---|---|
| simulated_all: every data type it needs | 728 |
| simulated_outcome_question: needs only a measured trait | 13 |
| simulated_part: some data types not simulated | 26 |
| no_simulator_for_data_type | 51 |
| not_run: no level count for a factor (time 56, genotype 55, stage 9, ...) | 148 |
| not_run: spec has no factors | 31 |
| not_run: other (Biodes "non-character argument" 2, augmented design without checks 1) | 3 |

**741 / 1000 (74%) get everything they need; 26 more get part of it.** The 100
`test_newkinds` questions do worse (59 / 100): they hold most of the new gaps below.

New gap: **genetics with no lines** (24 questions). BSA / pool-seq bulks,
population resequencing (WGS of populations, evolve-and-resequence) and GWAS or
eQTL whose design factor is a phenotype group or a population, not the lines.
These used to stop Simitall; now the rest of the design is simulated and
genetics is listed in `sim_truth.json$not_simulated` as `needs_genotype_factor`.
They need their own simulators: bulks of individuals with allele frequencies by
phenotype (BSA), and populations with drift and selection (pool-seq, WGS).

Missing data types here (gain = only gap): m6a_seq +10, scatac_seq +9,
genetics with no lines +13 (wgs 7, gwas 4, eqtl 2) and +5 more for BSA/pool-seq, ribo_seq +6, hic +5,
medip +5, metabolomics +5, rip_clip +5, degradome, lipidomics, mnase_seq,
spatial +2 each.

## What we can simulate (by data type) (first 1000)

| data type | questions | simulated | blocked |
|---|---|---|---|
| rnaseq (+ 3' tag-seq) | 330 | 305 | 25 |
| chip_seq (+ CUT&Tag) | 74 | 66 | 8 |
| atac_seq | 60 | 55 | 5 |
| gwas | 55 | 41 | 14 |
| single_cell | 43 | 40 | 3 |
| qtl_mapping | 35 | 21 | 14 |
| long_read_rna | 25 | 24 | 1 |
| small_rnaseq (template / native NB) | 140 | 132 | 8 |
| amplicon (16S/ITS) + metagenomics (template) | 55 | 48 | 7 |
| methylation (template / native CG-CHG-CHH / array) | 71 | 67 | 4 |

"Blocked" = the simulator exists but the spec gives no factor or level count
(GWAS / QTL panels mostly have no genotype factor).

## What we cannot simulate (no simulator), in build order

Gain = questions that become fully covered when this data type is added.

| data type | questions | gain | nearest existing machinery |
|---|---|---|---|
| m6a_seq (MeRIP) | 41 | +40 | IP + input counts per peak, like ChIP |
| ribo_seq | 29 | +26 | RNA + footprint counts per gene, translational-efficiency truth |
| hic | 23 | +21 | contact matrices with planted loops/TAD changes (new) |
| degradome | 20 | +20 | cleavage-site counts (mostly asked with small RNA) |
| medip | 18 | +18 | enrichment counts per window, like ChIP IP/input |
| scatac_seq | 18 | +16 | ATAC peaks x cells (combine ATAC and scRNA code) |
| mnase_seq | 11 | +10 | nucleosome positions, like ATAC |
| spatial | 10 | +9 | scRNA counts on a grid |
| metabolomics, lipidomics | 12 | +10 | log-normal feature matrix (easy) |
| rip_clip | 6 | +5 | IP + input counts, like ChIP |

Gain counts questions where this is the only missing data type.

## Blockers left, by owner

**Biodes**

1. No level count (58), no factors (46: mostly panels and populations whose size
   is in the prompt). These are questions for the user or reader rules.
2. Two reviewed specs make `studyforge_ask` fail with "non-character argument".

**Simitall**

3. "How many per group?" questions (361) have no `design.reps`. They should go
   to `simitall_power_sweep()` over a range of reps instead of stopping.

## Fixed since the first run

- Biodes: factor-aware design choice and `layout_from_spec()` (every factor a
  column, factorial/split crossing, within-unit factors such as tissue).
- Biodes 568c9b5: joined levels go to `treatment_combination`; Simitall maps
  each model term through `analysis.terms` (a factor named `treatment` can be
  the second factor).
- Simitall: layouts come from Biodes; factors the model does not hold get
  effects and a unit's rows share its effect; amplicon is microbiome, not
  genotyping; CRISPR/overexpression are treatments; `suggestions.tsv` no longer
  invents factors (it gave 97 plant questions a made-up genotype x water design).

## Today's end-to-end path (Biodes rules reader, no Qwen here)

303 / 1000 reach simulation with the rules reader alone (`rules_reader_stage`).
Most stop at the layout because the reader leaves out replicate counts or factor
levels. Rerun with the Qwen reader (`engine = "llm"`) before judging the reader.

## Rerun

```sh
export BIODES_DIR=~/Documents/Github/BioDesignGPT SIMITALL_DIR=~/Documents/Github/simitall
Rscript data-raw/coverage/run_gold_specs.R 1 1000 gold_run.tsv      # reviewed specs
Rscript data-raw/coverage/run_rules_reader.R 1 1000 rules_run.tsv     # rules reader
Rscript data-raw/coverage/route_gold_specs.R                          # routes only
QUESTION_SET=second Rscript data-raw/coverage/run_gold_specs.R 1 1000 second.tsv   # the next 1000
Rscript data-raw/coverage/summarize_coverage.R coverage.tsv second.tsv   # verdict per question
IDS=ids.txt Rscript data-raw/coverage/run_gold_specs.R 1 89 some.tsv  # only the question ids in ids.txt
```
