# Organism and technique lexicon

Deterministic tables for recognising *which organism* and *which technique or
analysis* a request is about. Read by `R/agent_lexicon.R`
(`.simitall_detect_organism()`, `.simitall_detect_techniques()`,
`.simitall_lexicon_checks()`). Rebuild with the scripts in `data-raw/lexicon/`.

| File | Rows | Content |
|---|---|---|
| `organisms.tsv` | 88 | taxid, NCBI scientific name, group and family (from the NCBI lineage), ploidy, selfing, mating system, approximate generation time, common reference assembly, SimItAll knowledge file |
| `organism_aliases.tsv` | 912 | alias -> organism; type (scientific name/synonym, common name, misspelling, abbreviation, plural, strain, panel, breed, genome build, informal), source, strength, case sensitivity |
| `techniques.tsv` | 69 | assays, sequencing platforms and analyses -> catalog analysis (`routes_to`), candidate analyses when the technique alone is ambiguous (bulk RNA-seq, scRNA-seq, WGS), parameter set (e.g. `assay=atacseq`, `platform=NovaSeq_X`), support (catalog status or `not_supported`), ontology IDs |
| `technique_aliases.tsv` | 762 | alias -> technique; source (ontology ID, SRA vocabulary, analysis catalog, curated), strength, case sensitivity |

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
```
Both scripts check every taxid against NCBI, every ontology ID against
EFO/OBI/EDAM, every route and enum value against the analysis catalog, and
report aliases that would point to two targets.

## Measured on the prompt -> StudySpec examples (229, batches 1-2)

`Rscript inst/agent/prompt_eval/score_lexicon.R`

| | current agent | lexicon only |
|---|---|---|
| named organism found | 40.7% | 98.6% |
| no organism invented | 97.8% | 96.6% |
| correct primary analysis | 45.2% | 68.8% |
| wrong analysis chosen | 61 | 30 |
| unsupported technique recognised | — | 92.9% |

Caveat: the examples and the curated aliases were written by the same author,
so these numbers are optimistic. Phrases copied from example prompts were
removed before scoring; a fair check needs prompts written after the lexicon
was frozen.
