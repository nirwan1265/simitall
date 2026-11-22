# Can Simitall simulate Biodes' questions? (1,000 paper questions, 2026-11-22)

Questions: Biodes `training/paper_questions/splits/test.jsonl` (764) + the first
236 of `dev.jsonl`, one question per paper, stratified by species x technique.
Per-question results: [`data-raw/coverage/coverage_1000.tsv`](../data-raw/coverage/coverage_1000.tsv).
Nothing was fixed in this pass; this is the list of what works and what does not.

Each question was run twice:

1. **Reviewed spec** (`expected_spec`, as if Biodes read the question
   perfectly) -> Biodes design choice -> layout -> `simitall_from_biodes()`.
   This measures Simitall plus the design hand-off.
2. **Biodes rules reader** on the prompt (no Qwen here) -> same pipeline,
   without `use_suggestions`. This measures today's end-to-end path.

## Result (reviewed specs)

| verdict | n | meaning |
|---|---|---|
| simulated_all | 277 | every data type the question needs was simulated |
| simulated_outcome_question | 37 | question needs only a measured trait; simulated |
| simulated_part | 59 | some of its data types simulated, others have no simulator |
| no_simulator_for_data_type | 189 | ran, but its data type has no simulator (only a placeholder outcome was made) |
| layout_drops_factor | 260 | ran, but the layout kept only one of the spec's factors, so the data test a different experiment |
| not_run_design | 136 | no layout could be built from the spec |
| not_run_simitall_error | 42 | Simitall stopped (genetics asked for, no genotype factor) |

**314 / 1000 (31%) get everything they need today; 59 more get part of it.**

## What we can simulate (by data type, reviewed specs)

| data type | questions | simulated | blocked by design/layout |
|---|---|---|---|
| rnaseq (+ 3' tag-seq) | 330 | 168 | 162 |
| chip_seq (+ CUT&Tag) | 74 | 48 | 26 |
| atac_seq | 60 | 35 | 25 |
| gwas | 55 | 38 | 17 |
| single_cell | 43 | 28 | 15 |
| qtl_mapping | 35 | 17 | 18 |
| long_read_rna | 25 | 17 | 8 |
| gbs, wgs, bsa_seq, pool_seq, targeted_capture | 39 | 18 | 21 |
| nascent_rna, trap_seq, eqtl | 15 | 6 | 9 |

These simulators work; what stops them is the design hand-off (below).

## What we cannot simulate (no simulator)

Gain = questions that become fully covered when this data type is added (with
the design issues fixed), in priority order:

| data type | questions | gain | nearest existing machinery |
|---|---|---|---|
| small_rnaseq (miRNA, siRNA, lncRNA, circRNA) | 140 | +120 | count matrix like `simulate_rnaseq_experiment()`, small-RNA feature set |
| methylation + medip | 89 | +86 | rebuild the lost `simulate_methylation()` / `call_dmrs()` (design in the handoff doc) |
| amplicon (16S/ITS) + metagenomics | 55 | +51 | taxa count matrix, Dirichlet-multinomial, treatment shifts on a few taxa |
| m6a_seq (MeRIP) | 41 | +40 | IP + input counts per peak, like ChIP |
| ribo_seq | 29 | +27 | RNA + footprint counts per gene, translational-efficiency truth |
| hic | 23 | +23 | contact matrices with planted loops/TAD changes (new) |
| scatac_seq | 18 | +18 | ATAC peaks x cells (combine ATAC and scRNA code) |
| degradome | 20 | +3 | cleavage-site counts (mostly asked with small RNA) |
| mnase_seq | 11 | +10 | nucleosome positions, like ATAC |
| spatial | 10 | +10 | scRNA counts on a grid |
| metabolomics, lipidomics | 12 | +8 | log-normal feature matrix with treatment effects (easy) |
| rip_clip | 6 | +6 | IP + input counts, like ChIP |

With all of them and the design issues fixed, every one of the 1,000 is
covered. The first five (small RNA, methylation, 16S, m6A, Ribo-seq) take
coverage of needs from 57% to 90%.

## Design / hand-off blockers (fix before adding simulators)

1. **Multi-factor specs get one-factor designs: 260 questions.** Spec has
   genotype x time, treatment x tissue, genotype x temperature, ...; Biodes picks
   `rcbd`/`alpha` (model `y ~ treatment + block`) and the layout keeps only the
   first factor. Fix in Biodes design choice (factorial / split-plot when the
   spec has >1 factor), or have the layout cross all spec factors inside the
   blocks. Biggest single win: ~124 RNA-seq questions alone.
2. **Factor with no level count: 33.** ("several time points"). A real
   question for the user; keep asking.
3. **Only factor is `time`: 37.** The layout needs a treatment factor and does
   not take time as one. Time-course-only studies should use time as the
   treatment (or repeated-measures design).
4. **No factors at all: 45.** GWAS panels, mapping populations, genebank
   collections ("800 accessions", "4,300 rats", "192 F1 progeny"). The number
   of lines is in the prompt but not in a factor. Biodes should record it (e.g.
   `population.size` or a genotype factor) so Simitall gets its lines.
5. **Genetics with no genotype term: 42.** 36 are the amplicon misroute below;
   the other 6 are #4.
6. Unreplicated blocked design (3), population-plan designs with no layout
   rule (2), Biodes error on 2 specs ("non-character argument").

## Bugs found in Simitall (to fix)

- **`amplicon` is routed to genetics.** In these questions amplicon means
  16S/ITS microbiome (46 of 48) or edit amplicons; it must leave the genetics
  row (`simulators.tsv`) and become a microbiome route.
- **`inst/suite/suggestions.tsv` invents factors.** With `use_suggestions =
  TRUE`, any plant prompt with open factors gets `genotype:6, water:2`; 97 of
  1,000 questions got a made-up genotype x water experiment. It is also a
  word-matching reader. Remove the factor suggestions (keep only numeric
  defaults like residual SD, listed as assumptions), or move them to Biodes.
- `crispr` / `overexpression` are listed as unsupported assays; they are
  treatments (the edited lines), not data to simulate.
- Factor reference level follows natural sort, not Biodes' listed order.

## Today's end-to-end path (Biodes rules reader)

526 / 1000 reached simulation. 430 stopped before a layout: in 191 the rules
reader found no factor that the reviewed spec has, in 134 it found the factor
but not its level count. Rerun this with the Qwen reader (`engine = "llm"`)
before judging the reader; the rules engine is the fallback.

## Rerun

```sh
BIODES_DIR=~/Documents/Github/BioDesignGPT SIMITALL_DIR=~/Documents/Github/simitall \
  Rscript data-raw/coverage/run_gold_specs.R 1 1000 gold_run.tsv      # reviewed specs
Rscript data-raw/coverage/run_rules_reader.R 1 1000 rules_run.tsv     # rules reader
Rscript data-raw/coverage/route_gold_specs.R                          # routes only
```
