# Can Simitall simulate Biodes' questions? (all 4,695 paper questions)

**First 1000:** Biodes `training/paper_questions/splits/test.jsonl` (764) + the
first 236 of `dev.jsonl`, one question per paper, stratified by species x
technique. **Second 1000** (`QUESTION_SET=second`): the rest of `dev.jsonl` (92),
`dev2` (155), `test_draw4` (206), `test_draw5` (159), `test_heldout` (87),
`test_newkinds` (100) and the first 201 of `train.jsonl`. **Third set**
(`QUESTION_SET=third`): the rest of `train.jsonl` (2,695), so every question
Biodes has.

Per-question results: [`coverage_1000.tsv`](../data-raw/coverage/coverage_1000.tsv),
[`coverage_second_1000.tsv`](../data-raw/coverage/coverage_second_1000.tsv),
[`coverage_third_2695.tsv`](../data-raw/coverage/coverage_third_2695.tsv)
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

## Current result (Biodes 5523171 schema 0.5.1, Simitall main, 2026-11-23)

| verdict | first 1000 | second 1000 | third 2695 | all 4695 |
|---|---|---|---|---|
| simulated_all: every data type it needs | 776 | 896 | 2519 | 4191 |
| simulated_outcome_question: needs only a measured trait | 46 | 16 | 34 | 96 |
| simulated_part: some data types not simulated | 71 | 31 | 66 | 168 |
| no_simulator_for_data_type | 84 | 41 | 57 | 182 |
| not_run: spec has no factors | 22 | 16 | 19 | 57 |
| not_run: other | 1 | 0 | 0 | 1 |

**4,287 / 4,695 (91%) get everything they need** (822, 912, 2,553); 168 more
get part of it. 793 of those use at least one flagged layout default. With
Biodes 0.5.1, finer assay kinds are routed (`assay_kinds`: DNase/FAIRE as
ATAC, DAP-seq/SELEX as ChIP); amplicon genotyping, TCR/BCR repertoires, CRISPR
screens, proteomics and STARR-seq are now counted as not simulated instead of
being simulated as the wrong data, and survival outcomes (12) simulate.

The "other": pq_0105 (yeast amplicon: no microbiome/amplicon template for
yeast and no fallback; Biodes has no yeast amplicon study to learn from either).

### History (fully simulated)

| run | first 1000 | second 1000 |
|---|---|---|
| first run | 314 | - |
| layouts from Biodes | 487 | - |
| + small RNA / microbiome templates | 644 | - |
| + methylation | 711 | 741 |
| + m6A | 751 | 751 |
| + Biodes 0.5 (panel genotype, defaults, genetics unit) | 820 | 912 |
| + Biodes b78d893 (pq_0180 fixed) | 821 | 912 |
| third set (train), same code | - | - | 2555 / 2695 |
| + Biodes 0.5.1 (assay kinds, generations, survival) | 822 | 912 | 2553 |
| Biodes d39936b (cell_type / fraction within units) | 822 | 912 | 2553 |
| Biodes 7d98031 (occasion factors, repeated_factors) | 822 | 912 | 2553 |
| Biodes ab01def (repeated designs, time values) | 822 | 912 | 2553 |

## Biodes ab01def, schema 0.5.3 (2026-11-25)

Same 4,287 / 4,695. `repeated_factors` set on 40 reviewed specs (read from
the text): `[time]` 12, `[age]` 4, `[stage]` 2, `[age,time]` 1, `[]` 21.
Of the 19 repeated, 8 get `randomized_complete_block_repeated`; 11 have no
randomized factor (only the occasion, e.g. "the same 3 cows at 15 and 30
days") and stay `rcbd` with each animal a block (a paired design). Time
factors: 787, of which 199 carry `values` (107 evenly spaced, 92 not, e.g.
6/12/24/48/72 h). Simitall spaces time effects by those values: pq_0178 gives
0, 0.91, 2.73, 6.36, 10 for 6/12/24/48/72 h.

## Biodes 7d98031, schema 0.5.2 (2026-11-24)

Occasion factors in the reviewed specs: time 787, stage 230, age 126. None of
the reviewed specs carry `repeated_factors`; read from the text it is set on
17: `[time]` 10, `[age]` 4, `[stage]` 2, `[]` 1 (pq_0476); absent (asked) on
the rest, which simulate as randomized to units. All 17 simulate fully.
pq_0527 is a false `[stage]`: "a mix of 5 embryos from the same sow" at
E16/E18/E21/E28 is destructive (the sow is shared, not the embryos).

## Biodes d39936b (2026-11-24)

Same 4,287 / 4,695 on the reviewed specs: no verdict and no factor count
changed (the new count reading is in Biodes' reader, which this check skips).
`cell_type` and `fraction` are now within-unit factors: 43 reviewed specs have
one, 19 of them changed from a factorial layout to an RCBD of the treatment
with each unit repeated across the cell types or fractions; all 43 run
(41 fully, 2 with a missing data type). No layout errors.

Biodes' reader alone on all 4,695 prompts, 1cbf2f4 vs d39936b
(`read_all_prompts.R`, `score_reader_counts.R`;
`biodes_reader_1cbf2f4_vs_d39936b.tsv` lists the 258 prompts whose factors
changed), scored against the reviewed specs' 5,982 factor counts:

| reader | count right | count wrong | count left open | factor not found | specs with every count right |
|---|---|---|---|---|---|
| 1cbf2f4 | 2,906 | 472 | 454 | 2,150 | 1,666 / 4,302 |
| d39936b | 3,037 | 426 | 446 | 2,073 | 1,758 / 4,302 |

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

## What we cannot simulate yet (all 4,695)

Gain = questions fully covered once this data type is added (it is their only gap).

| data type | questions | gain | nearest existing machinery |
|---|---|---|---|
| metabolomics, lipidomics | 55 | +55 | log-normal feature matrix (easy) |
| genetics with no lines (wgs, bsa_seq, pool_seq, gwas, ...) | 70 | +45 | bulks by phenotype tail; populations with drift/selection (`spec$genetics`) |
| ribo_seq | 43 | +42 | RNA + footprint counts per gene, translational-efficiency truth |
| scatac_seq | 41 | +41 | ATAC peaks x cells (combine ATAC and scRNA code) |
| hic | 38 | +38 | contact matrices with planted loops/TAD changes (new) |
| medip | 31 | +31 | enrichment per window, like `simulate_merip()` IP/input |
| degradome | 22 | +22 | cleavage-site counts (mostly asked with small RNA) |
| spatial | 19 | +19 | scRNA counts on a grid |
| mnase_seq | 17 | +16 | nucleosome positions, like ATAC |
| rip_clip | 14 | +14 | IP + input counts, like `simulate_merip()` |

## Simulated, but not the right data (Simitall)

These run, but the data type is not what was asked:

- Amplicon that is not a microbial community (10): amplicon genotyping
  (CRISPR alleles, mitochondrial haplotypes), DSB-repair amplicons, immune
  repertoires. They route to the 16S template today.
- CRISPR / Tn-seq screens (6): guide or insertion counts are the data; `crispr`
  routes to "treatment" (a way to make lines) and nothing is simulated for it.
- Long-read RNA (13 isoform questions): gene counts only, no isoform truth.
- Proteomics (10), TCR/BCR repertoires (8), STARR-seq/MPRA (2): no simulator
  and no assay value from Biodes.
- Outcome types: binary, count, ordinal and survival simulators exist, but only
  2 of 4,695 specs set `primary_outcome.type`, so every outcome is continuous.

## Blockers left, by owner

**Biodes** (reading; details in the list sent to Biodes, 2026-11-22)

1. 70 specs have no factor at all.
2. Assays not read: 18 prompts name an assay the spec leaves out; DNase-seq,
   FAIRE-seq, DAP-seq, SELEX (7) have no assay value; amplicon kinds, CRISPR
   screens, proteomics, repertoires, reporter assays have none either.
3. `primary_outcome.type` is set in 2 of 4,695 specs.
4. Breeding notation (BC3F2, F2:3, BC1F3:5; 61 prompts) is not read into the
   number of backcross and selfing generations (pq_2667, BC1F3:5, got the NIL
   default of 3 backcrosses).
5. `practice.R:148` warns when `organism` has more than one value.

**Simitall**

6. The simulators in the two tables above; a fallback for amplicon data in
   organisms with no template (pq_0105).

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
QUESTION_SET=third Rscript data-raw/coverage/run_gold_specs.R 1 2695 third.tsv    # the rest (train)
Rscript data-raw/coverage/summarize_coverage.R coverage_1000.tsv first.tsv       # verdict per question
IDS=ids.txt Rscript data-raw/coverage/run_gold_specs.R 1 89 some.tsv              # only these question ids
FILL_DEFAULTS=0 ...                                                               # stop at missing level counts instead
Rscript data-raw/coverage/run_rules_reader.R 1 1000 rules_run.tsv                 # Biodes rules reader end to end
```
