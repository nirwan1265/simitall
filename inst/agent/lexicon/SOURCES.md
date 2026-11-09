# Organism and technique lexicon

Deterministic tables for reading a request: *which organism*, *which technique
or analysis*, *which settings and numbers*, and *what the user wants done with
which data*. Read by `R/agent_lexicon.R`; `.simitall_lexicon_read_prompt()`
combines them into a StudySpec-shaped view. Rebuild the generated tables with
the scripts in `data-raw/lexicon/`.

| Table | Files | Job |
|---|---|---|
| 1. organisms | `organisms.tsv`, `organism_aliases.tsv` | organism + biology facts |
| 2. techniques | `techniques.tsv`, `technique_aliases.tsv` | assay / platform / analysis -> catalog route or `not_supported` |
| 3. methods and design terms | `settings.tsv` (generated), `quantities.tsv` | words -> catalog parameter values ("brother-sister mating" -> `mating_scheme = sib_mating`, "EMMAX" -> `linear_mixed_model`); unit- or noun-anchored numbers ("h2 0.4", "30x", "1 Mb window", "4,000 cases") |
| 4. intent | `intent.tsv`, `intent_routes.tsv` | goal (simulate / analyze / both / plan), data source (provided / synthetic / demo / ask), "only", negation ("no eQTL", "besides"), and guard cues (real data requested, clinical use, mid-sentence corrections); ambiguous assay x goal -> analysis |

| File | Rows | Content |
|---|---|---|
| `organisms.tsv` | 88 | taxid, NCBI scientific name, group and family (from the NCBI lineage), ploidy, selfing, mating system, approximate generation time, common reference assembly, SimItAll knowledge file |
| `organism_aliases.tsv` | 912 | alias -> organism; type (scientific name/synonym, common name, misspelling, abbreviation, plural, strain, panel, breed, genome build, informal), source, strength, case sensitivity |
| `techniques.tsv` | 69 | assays, sequencing platforms and analyses -> catalog analysis (`routes_to`), candidate analyses when the technique alone is ambiguous (bulk RNA-seq, scRNA-seq, WGS), parameter set (e.g. `assay=atacseq`, `platform=NovaSeq_X`), support (catalog status or `not_supported`), ontology IDs |
| `technique_aliases.tsv` | 767 | alias -> technique; source (ontology ID, SRA vocabulary, analysis catalog, curated), strength, case sensitivity |
| `settings.tsv` | 367 | alias -> analysis, field, value; 91 rows are the catalog's own enum values, the rest curated synonyms; every value validated against the catalog |
| `quantities.tsv` | 68 | regex (one capture) -> field, analyses it applies to, transform (int, num, frac, bp, range, count, chr); number words and thousands separators normalised first |
| `intent.tsv` | 42 | regex cue -> signal (simulate, benchmark, analyze, plan, provided, demo, only, negation, real_data, clinical, correction) |
| `intent_routes.tsv` | 9 | ambiguous assay + goal -> analysis (simulate scRNA-seq -> `scrnaseq_simulation`) |

**Strength.** Weak aliases ("flies", "patients", "N2", "PCA", "annotation") count
only when no strong alias matched. **Case-sensitive** aliases ("DO", "CC", "B6",
"PBS", "ASE", "IBD") are matched exactly so "do", "pbs" or inflammatory-bowel
"IBD" in lower case do not trigger them.

## Sources and versions

| Source | Version | Used for | Licence |
|---|---|---|---|
| NCBI Taxonomy via OBO Foundry NCBITaxon (`ncbitaxon.obo.gz`, github.com/obophenotype/ncbitaxon) | release v2026-07-12 | taxids, scientific names, GenBank/common names, older Latin synonyms, misspellings, lineage | NCBI data: public domain; NCBITaxon OBO: CC0 |
| EFO, Experimental Factor Ontology (github.com/EBISPOT/efo) | v3.95.0 | assay and instrument names and exact synonyms | Apache 2.0 |
| OBI, Ontology for Biomedical Investigations (github.com/obi-ontology/obi) | 2026-07-27 | assay and instrument names and exact synonyms | CC BY 4.0 |
| EDAM (github.com/edamontology/edamontology, `releases/EDAM.tsv`) | current release | analysis labels (weak aliases) | CC BY-SA 4.0 |
| ENA/SRA schema (github.com/enasequence/schema) | current | library strategies, instrument models | Apache 2.0 |
| SimItAll analysis catalog | this repository | routes, parameters, support status, catalog aliases | package licence |
| Curated (`data-raw/lexicon/*_curated.tsv`, `organism_extra_aliases.tsv`) | 2026-11 | organism list and biology facts; strains, panels, breeds, cultivars, genome builds; technique routing and extra aliases | package licence |

Curated biology facts (ploidy, selfing, generation time, reference) are
summaries for planning, not authoritative; generation times are approximate.
NCBI renamed pea to *Lathyrus oleraceus*; *Pisum sativum* stays as a synonym.

## Rebuild

```sh
python3 data-raw/lexicon/build_organisms.py  /path/ncbitaxon.obo.gz
python3 data-raw/lexicon/build_techniques.py /path/efo.obo /path/obi.obo /path/EDAM.tsv
python3 data-raw/lexicon/build_settings.py
```
Both scripts check every taxid against NCBI, every ontology ID against
EFO/OBI/EDAM, every route and enum value against the analysis catalog, and
report aliases that would point to two targets.

## Measured on the prompt -> StudySpec examples (229, batches 1-2)

`Rscript inst/agent/prompt_eval/score_reader.R [test]` scores three readers with
the same scorer as `score_examples.R`: the current agent, the four tables alone,
and the tables with the agent filling fields the tables left empty.

| all 229 / frozen test 152 | current agent | tables only | tables + agent |
|---|---|---|---|
| reading recall (expected fields read) | 27.0 / 25.5 | 78.7 / 79.6 | 82.1 / 82.7 |
| goal and data source | 6.7 / 6.6 | 80.0 / 80.5 | 80.2 / 80.5 |
| primary analysis | 45.2 / 43.6 | 75.0 / 77.1 | 76.4 / 78.6 |
| secondary analyses exact | 69.9 / 67.8 | 83.0 / 81.6 | 75.1 / 74.3 |
| forbidden analyses chosen | 7 / 3 | 0 / 0 | 4 / 2 |
| stops when it should | 80.8 / 87.5 | 96.2 / 93.8 | 96.2 / 93.8 |
| examples with invented values | 9 / 5 | 1 / 1 | 10 / 6 |
| defaulted instead of asking | 5 / 4 | 0 / 0 | 5 / 4 |

Tables alone are the safer reader: adding the agent's values gains ~3 points of
reading but brings back its invented and defaulted values. The remaining misses
are mostly numbers that need role assignment ("6 healthy and 6 sick people" ->
12 donors, "control and drought" -> conditions, "10 founder pairs") — the job of
the trained prompt model.

Caveat: the examples, the curated aliases and the rules were written by the same
author, so these numbers are optimistic. Phrases copied from example prompts were
removed; a fair check needs prompts written after the tables were frozen.
