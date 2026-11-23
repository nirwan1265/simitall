# Can Simitall simulate Biodes' questions? (2,000 paper questions)

**First 1000:** Biodes `training/paper_questions/splits/test.jsonl` (764) + the
first 236 of `dev.jsonl`, one question per paper, stratified by species x
technique. **Second 1000** (`QUESTION_SET=second`): the rest of `dev.jsonl` (92),
`dev2` (155), `test_draw4` (206), `test_draw5` (159), `test_heldout` (87),
`test_newkinds` (100) and the first 201 of `train.jsonl`.

Per-question results: [`coverage_1000.tsv`](../data-raw/coverage/coverage_1000.tsv),
[`coverage_second_1000.tsv`](../data-raw/coverage/coverage_second_1000.tsv)
(`verdict`, `needs_status`, `blocker`, `defaults_filled`, `solve_for`,
`genetic_unit`, what was simulated).

Each question is run on its **reviewed spec** (`expected_spec`, as if Biodes read
it perfectly), then the steps Biodes' reader runs after extraction (panel
genotype factor, `spec$genetics`, `mark_asked()`), then
`fill_layout_defaults()` for level counts the user left open (as when the user
does not answer; each is an assumption) -> Biodes design choice ->
`studyforge::layout_from_spec()` -> `write_design_package()` ->
`simitall_from_biodes()`. Questions with no replicate count (`power.solve_for =
"reps"`, 996 of 2000) run with 3 reps here to keep the run fast;
`simitall_from_question()` chooses their count with a power sweep.

## Current result (Biodes 3675658 schema 0.5, Simitall main, 2026-11-22)

| verdict | first 1000 | second 1000 | both |
|---|---|---|---|
| simulated_all: every data type it needs | 768 | 895 | 1663 |
| simulated_outcome_question: needs only a measured trait | 52 | 17 | 69 |
| simulated_part: some data types not simulated | 65 | 21 | 86 |
| no_simulator_for_data_type | 85 | 47 | 132 |
| not_run: spec has no factors | 28 | 20 | 48 |
| not_run: other | 2 | 0 | 2 |

**1,732 / 2,000 (87%) get everything they need** (820 and 912); 86 more get
part of it. 233 of those use at least one flagged layout default (genotype count
125, time 75, stage 12, tissue 9, ...). Without the defaults: 747 and 752.

The two "other": pq_0105 (yeast amplicon: no microbiome/amplicon template for
yeast and no fallback) and pq_0180 (Biodes: the panel step adds a `genotype`
factor with no level count and `studyforge_ask()` stops in `.divisors()` with
"argument must be coercible to non-negative integer").

### History (fully simulated)

| run | first 1000 | second 1000 |
|---|---|---|
| first run | 314 | - |
| layouts from Biodes | 487 | - |
| + small RNA / microbiome templates | 644 | - |
| + methylation | 711 | 741 |
| + m6A | 751 | 751 |
| + Biodes 0.5 (panel genotype, defaults, genetics unit) | 820 | 912 |

## What we can simulate (by data type, both sets)

| data type | questions | simulated |
|---|---|---|
| RNA-seq (+ 3' tag, long-read, nascent, TRAP) | 1028 | 1025 |
| small RNA-seq (template / native NB) | 192 | 191 |
| GWAS / GBS / WGS / capture panels | 148 | 117 |
| ChIP-seq (+ CUT&Tag) | 110 | 106 |
| methylation (template / native CG-CHG-CHH / array) | 108 | 107 |
| amplicon (16S/ITS) + metagenomics (template) | 98 | 95 |
| ATAC-seq | 90 | 88 |
| single-cell RNA-seq | 84 | 78 |
| QTL mapping / BSA / pool-seq | 81 | 65 |
| m6A MeRIP-seq | 51 | 51 |
| eQTL | 10 | 9 |

Panels and mapping fall short where the genetic unit is a bulk or a population
(`needs_genotype_factor`, 34 questions): BSA / pool-seq bulks, population
resequencing. Biodes 0.5 writes `spec$genetics` (`unit`, `individuals_per_bulk`,
`bulk_fraction`, `generations`) for the simulators that will cover them.

## What we cannot simulate yet

Gain = questions fully covered once this data type is added (it is their only gap).

| data type | first | second | gain (both) | nearest existing machinery |
|---|---|---|---|---|
| ribo_seq | 28 | 6 | +33 | RNA + footprint counts per gene, translational-efficiency truth |
| scatac_seq | 18 | 11 | +29 | ATAC peaks x cells (combine ATAC and scRNA code) |
| hic | 23 | 5 | +28 | contact matrices with planted loops/TAD changes (new) |
| medip | 18 | 5 | +23 | enrichment per window, like `simulate_merip()` IP/input |
| degradome | 20 | 2 | +22 | cleavage-site counts (mostly asked with small RNA) |
| genetics with no lines (wgs, bsa_seq, pool_seq, gwas, ...) | 5 | 29 | +23 | bulks by phenotype tail; populations with drift/selection |
| metabolomics, lipidomics | 12 | 7 | +19 | log-normal feature matrix (easy) |
| mnase_seq | 11 | 2 | +12 | nucleosome positions, like ATAC |
| spatial | 10 | 2 | +12 | scRNA counts on a grid |
| rip_clip | 6 | 5 | +11 | IP + input counts, like `simulate_merip()` |

## Blockers left, by owner

**Biodes**

1. 48 specs have no factor at all (Biodes now asks "what do you compare or
   map?").
2. pq_0180: panel genotype factor with no level count crashes `.sf_recipe()`.
3. `practice.R:148` warns when `organism` has more than one value.

**Simitall**

4. The simulators in the table above; a fallback for amplicon data in
   organisms with no template.

## Fixed since the first run

- Biodes: factor-aware design choice and `layout_from_spec()`; joined levels go
  to `treatment_combination`; schema 0.5: flagged layout defaults, `power.solve_for`,
  panel genotype factor, per-environment layouts, `spec$genetics`, NA-safe organism.
- Simitall: layouts come from Biodes; factors the model does not hold get
  effects; amplicon is microbiome; CRISPR/overexpression are treatments; genetics
  with no genotype factor simulates the rest; "how many per group?" runs a power
  sweep; templates for small RNA, microbiome and methylation; m6A.

## Rerun

```sh
export BIODES_DIR=~/Documents/Github/BioDesignGPT SIMITALL_DIR=~/Documents/Github/simitall
Rscript data-raw/coverage/run_gold_specs.R 1 1000 first.tsv                      # first 1000
QUESTION_SET=second Rscript data-raw/coverage/run_gold_specs.R 1 1000 second.tsv  # second 1000
Rscript data-raw/coverage/summarize_coverage.R coverage_1000.tsv first.tsv       # verdict per question
IDS=ids.txt Rscript data-raw/coverage/run_gold_specs.R 1 89 some.tsv              # only these question ids
FILL_DEFAULTS=0 ...                                                               # stop at missing level counts instead
Rscript data-raw/coverage/run_rules_reader.R 1 1000 rules_run.tsv                 # Biodes rules reader end to end
```
