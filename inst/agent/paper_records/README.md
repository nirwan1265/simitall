# Paper records

`method_records.tsv` holds method decisions extracted from papers, one row per
decision (for example "how many PCs to use in GWAS" or "how to normalize
single-cell counts"). Each row says what the paper recommends, where it says
it, and how simitall would use it.

This folder sits outside `inst/agent/knowledge/` on purpose: the retrieval
index never reads it, and the agent only uses rows a person has reviewed.

## How the agent uses rows

- `review_status = reviewed` — used. Relevant rows are added to spec plans
  (Limitations panel) for GWAS, RNA-seq, and genomic-selection workflows.
- `unreviewed` — ignored until someone checks it.
- `excluded` — never used (for example, a statement that is wrong in the source).

## How to review

1. Open the file in Excel or a spreadsheet.
2. For each row, check `guidance` against the paper at `page`.
   Page numbers are pages of the PDF file, which can differ from the journal's
   printed page numbers. Rows with `page = abstract` came from the abstract of
   an online open-access version.
3. Set `review_status` to `reviewed` or `excluded`, and put your initials in
   `reviewed_by`. Fix the wording if needed.
4. Save as tab-separated text, keeping the header row.

Read the `flag` column first: it marks a known error in a source, two papers
that disagree, unverified numbers, and single-study examples.

## Columns

| Column | Meaning |
|---|---|
| record_id | Stable id: domain plus number |
| domain | gwas, rnaseq, scrnaseq, scatac, atacseq, qtl, genomic_selection, sequencing, pangenome, microbiome, genomic_lm |
| decision | The choice the row is about |
| guidance | What the paper says, paraphrased (no copied text) |
| numbers | Key values, if any |
| applies_to | Organism, data type, or scope |
| simitall_use | default, guided_question, validation_check, design_rule, simulation_parameter, benchmark_idea, glossary |
| source, doi, page | Where it comes from |
| confidence | high (direct statement or table), medium (review summary), low (unclear or single study) |
| flag | Known problems or caveats |
| review_status, reviewed_by | Your review |

## Coverage and gaps

- 129 rows come from the 34 PDFs in the local paper folder (mostly reviews), and
  13 from open-access primary papers read online (Platt 2010, Kang 2008,
  Zhou & Stephens 2012, Liu 2016, Yin 2020, Love 2014, Hafemeister & Satija
  2019, Meuwissen 2001, Yu 2008, Gage 2020, Li 2011).
- Extraction covered passages that state a method choice or a number, not every
  sentence; a paper can contain more than its rows show.
- Not yet covered: most of the ~100 citation-only cards in
  `knowledge/papers/` (their bodies are templates without method content),
  and field-trial design (randomization, checks, border effects), which none of
  the current papers cover.
