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

## Current result (Biodes 568c9b5, Simitall main, 2026-11-22)

| verdict | n | first run |
|---|---|---|
| simulated_all: every data type it needs | 437 | 277 |
| simulated_outcome_question: needs only a measured trait | 50 | 37 |
| simulated_part: some data types have no simulator | 89 | 59 |
| no_simulator_for_data_type | 313 | 189 |
| not_run: no level count for a factor | 58 | 33 |
| not_run: spec has no factors | 46 | 45 |
| not_run: other (genetics without genotype factor 2, Biodes error 2, ...) | 7 | 42 |
| layout_drops_factor | 0 | 260 |

**487 / 1000 (49%) get everything they need; 89 more get part of it** (31% at
the first run). Among the 639 questions that give a replicate count: 322 all,
68 part. Every design that can be laid out now runs; what is left is missing
simulators (313) and questions for the user (111).

## What we can simulate (by data type)

| data type | questions | simulated | blocked |
|---|---|---|---|
| rnaseq (+ 3' tag-seq) | 330 | 305 | 25 |
| chip_seq (+ CUT&Tag) | 74 | 66 | 8 |
| atac_seq | 60 | 55 | 5 |
| gwas | 55 | 41 | 14 |
| single_cell | 43 | 40 | 3 |
| qtl_mapping | 35 | 21 | 14 |
| long_read_rna | 25 | 24 | 1 |

"Blocked" = the simulator exists but the spec gives no factor or level count
(GWAS / QTL panels mostly have no genotype factor).

## What we cannot simulate (no simulator), in build order

Gain = questions that become fully covered when this data type is added.

| data type | questions | gain | nearest existing machinery |
|---|---|---|---|
| small_rnaseq (miRNA, siRNA, lncRNA, circRNA) | 140 | +120 | count matrix like `simulate_rnaseq_experiment()`, small-RNA feature set |
| methylation + medip | 89 | +86 | rebuild `simulate_methylation()` / `call_dmrs()` (design in the handoff doc) |
| amplicon (16S/ITS) + metagenomics | 55 | +51 | taxa count matrix, Dirichlet-multinomial (route `assay_microbiome`, planned) |
| m6a_seq (MeRIP) | 41 | +40 | IP + input counts per peak, like ChIP |
| ribo_seq | 29 | +27 | RNA + footprint counts per gene, translational-efficiency truth |
| hic | 23 | +23 | contact matrices with planted loops/TAD changes (new) |
| scatac_seq | 18 | +18 | ATAC peaks x cells (combine ATAC and scRNA code) |
| degradome | 20 | +3 | cleavage-site counts (mostly asked with small RNA) |
| mnase_seq | 11 | +10 | nucleosome positions, like ATAC |
| spatial | 10 | +10 | scRNA counts on a grid |
| metabolomics, lipidomics | 12 | +8 | log-normal feature matrix (easy) |
| rip_clip | 6 | +6 | IP + input counts, like ChIP |

The first five take the ceiling (all data types needed have a simulator) from
57% to 90%.

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
```
