# PmiREN3.1 FTP Download Guide

This directory provides downloadable miRNA annotations, sequences, expression matrices, predicted and degradome-supported miRNA-target relationships, target-gene functional annotations, and other multi-omics datasets from PmiREN3.1.

PmiREN website: https://pmiren.com

Last updated: 2026-09-26

## Quick navigation

| User need | Recommended location |
| --- | --- |
| Data for one species | `miRNA_info/<Species>/` |
| Basic miRNA annotations | `miRNA_info/<Species>/<Species>_basicinfo.tsv` |
| Mature, star, or precursor sequences | `miRNA_info/<Species>/*.fa` |
| Genomic coordinates | `miRNA_info/<Species>/*.bed` |
| miRNA expression matrices | `miRNA_info/<Species>/*mature.lib.exp.tsv` and `*mature.tissue.exp.tsv` |
| Annotation pipeline outputs and raw counts | `miRNA_annotation/<Species>/` |
| miRNA expression overview PDF | `miRNA_annotation/<Species>/*-miRNA-expression.pdf` |
| One library's mature miRNA RPM values (server-side supplementary files) | `../sup_function/sRNA_exp/<Species>/<Species>_<library>.tsv` |
| miRNA-target relationships for one species | `GSM_info/Target_predict/Target_gene_info/<Species>_targetGene.txt` |
| All miRNA-target relationships | `GSM_info/Target_predict/Target_gene_info/All_species_targetGene.txt` |
| Tool-specific target evidence | `GSM_info/Target_predict/psRNAtarget/`, `psRobot/`, and `Degradome/` |
| InterProScan annotations of target genes | `GSM_info/Target_annotation/All_target_interproscan.tsv` |
| KEGG annotations of target genes | `GSM_info/Target_annotation/All_target_KEGG.tsv` |
| DNA and RNA methylation records | `GSM_info/Methylation/` |
| Iso-Seq-supported miRNA loci | `GSM_info/iso_seq/Iso_seq.tsv` |
| Conserved family-level regulatory networks | `GSM_info/Core_GRN/` |

## Directory structure

```text
FTP_download/
├── README.md
├── miRNA_info/
│   └── <Species>/
│       ├── <Species>_basicinfo.tsv
│       ├── <Species>_mature.fa
│       ├── <Species>_mature.bed
│       ├── <Species>_star.fa
│       ├── <Species>_star.bed
│       ├── <Species>_stemloop.fa
│       ├── <Species>_stemloop.bed
│       ├── <Species>_mature.lib.exp.tsv
│       ├── <Species>_mature.tissue.exp.tsv
│       └── <Species>_targetGene.txt
├── miRNA_annotation/
│   └── <Species>/
│       ├── <Species>-all-mod_fp_prediction
│       ├── <Species>-basic-info
│       ├── <Species>-mature.count
│       ├── <Species>-mature.exp
│       ├── <Species>_annotation.log
│       └── <Species>-miRNA-expression.pdf
└── GSM_info/
    ├── Target_predict/
    │   ├── Target_gene_info/
    │   │   ├── All_species_targetGene.txt
    │   │   └── <Species>_targetGene.txt
    │   ├── psRNAtarget/
    │   │   └── <Species>_psRNAtarget.tsv
    │   ├── psRobot/
    │   │   └── <Species>_psRobot.tsv
    │   └── Degradome/
    │       └── <Species>_<library>.cleaveland.log
    ├── Target_annotation/
    │   ├── All_target_interproscan.tsv
    │   └── All_target_KEGG.tsv
    ├── Methylation/
    │   ├── All_DNA_methylation.tsv
    │   └── All_RNA_methylation.tsv
    ├── iso_seq/
    │   └── Iso_seq.tsv
    └── Core_GRN/
        ├── Brassicaceae_Core_GRN.tsv
        ├── Fabaceae_Core_GRN.tsv
        └── Poaceae_Core_GRN.tsv
```

File availability can vary among species because it depends on the available reference genome, sequencing libraries, and supporting datasets.

As of the update date, `miRNA_info/` contains 321 species/taxon basic-annotation tables and 227 library expression matrices. `miRNA_annotation/` contains 239 annotation batches across 237 species directories, with 239 expression PDFs. These are different data products; their row counts and taxon coverage need not agree.

## Species-level miRNA data

Each subdirectory under `miRNA_info/` is named after a species, subspecies, variety, or cultivar. Spaces in taxon names are represented by underscores.

### Basic annotations

`<Species>_basicinfo.tsv` is a tab-separated table containing the core miRNA annotation. Its fields include:

- `PmiREN_ID`
- `Species`
- `miRNA_Family`
- precursor chromosome, start, end, strand, sequence, and secondary structure
- mature miRNA name, coordinates, and sequence
- star miRNA name, coordinates, and sequence

### Sequence files

The FASTA files contain sequences for different miRNA features:

- `<Species>_mature.fa`: mature miRNA sequences
- `<Species>_star.fa`: miRNA star sequences
- `<Species>_stemloop.fa`: precursor stem-loop sequences

### Coordinate files

The corresponding `.bed` files provide genomic coordinates for mature, star, and precursor features:

- `<Species>_mature.bed`
- `<Species>_star.bed`
- `<Species>_stemloop.bed`

These files follow a BED-compatible, tab-delimited layout.

### Expression matrices

- `<Species>_mature.lib.exp.tsv`: mature miRNA expression organized by sequencing library or sample
- `<Species>_mature.tissue.exp.tsv`: mature miRNA expression organized by tissue or sample description

The first column, `MIR`, contains the mature miRNA identifier. The remaining columns represent libraries, samples, or tissue labels. Library expression values are RPM (reads per million). Tissue labels can repeat: a shared label does not imply that those columns have been combined into a single tissue average. Some species may not have expression data for every data category.

### Per-library RPM tables (supplementary server directory)

The library matrices have also been split into individual two-column tables at the sibling location `../sup_function/sRNA_exp/`, outside `FTP_download/`:

```text
PmiREN3.1_database/sup_function/sRNA_exp/
├── split_manifest.json
└── <Species>/
    └── <Species>_<library>.tsv
```

Each file has a header `miRNA\tRPM`. The first column is the mature miRNA name and the second is its RPM in that library. Names, values, row order, and zero values are retained from the corresponding `<Species>_mature.lib.exp.tsv`; no additional normalization or filtering is applied. There are 5,132 library files across 227 species directories. `split_manifest.json` records source files, library names, row counts, and checksums. This server-side directory is not included automatically when downloading only `FTP_download/`.

### Species-level target relationships

`<Species>_targetGene.txt` contains miRNA-target relationships and evidence fields. The current columns are:

`ID`, `miRNA_gene_ID`, `targetgeneID`, `psRNAtarget`, `psRobot`, `Degradome`, `DegradomePic`, `CreateTime`, and `basicID`.

The `psRNAtarget`, `psRobot`, and `Degradome` columns record the available evidence or score from the corresponding analysis. `NA` indicates that no value is recorded for that field. Users should consult the PmiREN3.1 methods and website documentation when interpreting prediction scores, degradome categories, and confidence levels.

## Annotation pipeline outputs and expression plots

`miRNA_annotation/<Species>/` provides annotation pipeline outputs and associated expression matrices. File prefixes use the containing species directory name. Where two batches are retained, an additional identifier distinguishes them, such as `Citrus_x_clementina__mirdeep3_86_new_525` or `Pyrus_x_bretschneideri__Pbr` and `Pyrus_x_bretschneideri__Pyrxbre`.

| File suffix | Content |
| --- | --- |
| `-basic-info` | Tab-separated annotation table without a header; 22 fields described below |
| `-mature.count` | Raw count matrix with a header: `MIR`, `MIR_ID`, then library columns |
| `-mature.exp` | Pipeline expression matrix with a header: `MIR`, `MIR_ID`, then library columns |
| `-all-mod_fp_prediction` | Original pipeline prediction output; retains its pipeline-specific format |
| `_annotation.log` | Pipeline processing messages and warnings; not an expression table |
| `-miRNA-expression.pdf` | Expression overview derived from the corresponding raw count matrix |

The 22 fields in `-basic-info`, in order, are:

1. miRNA gene name
2. Internal `MIR_ID`
3. miRNA family
4. Species
5. Chromosome
6. Strand
7. Flanked precursor start
8. Flanked precursor end
9. Flanked precursor sequence
10. Flanked precursor secondary structure
11. Stem-loop start
12. Stem-loop end
13. Stem-loop sequence
14. Stem-loop secondary structure
15. Mature miRNA name
16. Mature miRNA start
17. Mature miRNA end
18. Mature miRNA sequence
19. Star miRNA name
20. Star miRNA start
21. Star miRNA end
22. Star miRNA sequence

`MIR` in these pipeline expression matrices contains the gene/precursor name (for example, `Ath-MIR156a`), whereas the `miRNA_info/` library matrices use mature names (for example, `Ath-miR156a`). Join pipeline annotations and expression rows using species, batch, and internal `MIR_ID`; internal IDs are not globally unique. Use `miRNA_info/*/*_mature.lib.exp.tsv` or the supplementary per-library files when RPM values are required.

Names in the current pipeline tables were updated where species, coordinates, and sequences matched the supplied renamed annotations. Records confirmed absent from the supplied August 2026 reference table were removed from `-basic-info`, `-mature.count`, and `-mature.exp`. Unresolved records were retained. The prediction files retain their original rows; processing logs retain their historical messages, with unambiguous name substitutions where available. These files should not be assumed to contain identical record sets. The separate `miRNA_info/` exports were not regenerated by this pipeline-table update.

Each PDF shows all retained miRNA rows and libraries as a heatmap of `log2(raw count + 1)`, alongside up to 30 rows ranked by mean raw count. Bars show `log2(mean raw count + 1)`. No between-library normalization is applied in these plots. Rows are ranked by mean count, duplicate names remain separate, and large matrices display only selected library labels. An input matrix with a header but no miRNA rows produces a PDF stating `No miRNA records to plot`.

## Integrated and multi-species data

### Target prediction

`GSM_info/Target_predict/` provides a centralized location for miRNA-target relationship files.

- `Target_gene_info/<Species>_targetGene.txt`: relationships for an individual species
- `Target_gene_info/All_species_targetGene.txt`: combined relationships across all available species
- `psRNAtarget/<Species>_psRNAtarget.tsv`: tool-specific target prediction results
- `psRobot/<Species>_psRobot.tsv`: tool-specific target prediction results
- `Degradome/<Species>_<library>.cleaveland.log`: CleaveLand analysis logs for individual libraries

The tool-specific TSV files may not have a header, and their schemas differ. CleaveLand logs are textual analysis outputs rather than a uniform TSV table. Use the headered `Target_gene_info/` tables for a common evidence schema across methods.

The combined file is useful for comparative or bulk analyses, while the species-specific files are more convenient for targeted analyses.

### Target-gene functional annotations

`GSM_info/Target_annotation/` contains functional annotations linked to miRNA-target relationships.

#### InterProScan

`All_target_interproscan.tsv` contains compiled InterProScan results for target genes across the database. Key fields include:

`Organism`, `miRNA`, `Target_gene`, `Analysis`, `Signature_acce`, `Signature_desc`, `InterPro_acce`, and `InterPro_desc`.

This is a large file of approximately 1.1 GB. Users may wish to filter it by organism, miRNA, or target-gene identifier after downloading.

#### KEGG

`All_target_KEGG.tsv` contains KEGG orthology annotations for target genes. Key fields include:

`Organism`, `miRNA`, `target_gene`, `KO`, `Definition`, `Score`, `Second`, and `Best`.

This file is approximately 106 MB.

### DNA and RNA methylation

`GSM_info/Methylation/` contains:

- `All_DNA_methylation.tsv`
- `All_RNA_methylation.tsv`

These tables include modification type, species, miRNA or miRNA family, genomic coordinates, score, strand, source library, and sample description.

### Iso-Seq support

`GSM_info/iso_seq/Iso_seq.tsv` contains Iso-Seq-supported miRNA loci. Fields include organism, miRNA, chromosome, start, end, length, strand, source library, sequence, and the linked basic record identifier.

### Core gene regulatory networks

`GSM_info/Core_GRN/` contains conserved regulatory-network edge tables for three plant families:

- `Brassicaceae_Core_GRN.tsv`
- `Fabaceae_Core_GRN.tsv`
- `Poaceae_Core_GRN.tsv`

Each table contains `source`, `target`, `edge_type`, `weight`, and `species`. The `species` field may contain multiple semicolon-separated species.

## Recommended workflows

### Retrieve all data for one species

1. Open `miRNA_info/<Species>/`.
2. Use `<Species>_basicinfo.tsv` for the primary miRNA annotations.
3. Use the FASTA and BED files for sequence and coordinate analyses.
4. Use the expression matrices when expression data are available.
5. Use `<Species>_targetGene.txt` for the corresponding miRNA-target relationships.

### Perform a multi-species miRNA-target analysis

1. Download `GSM_info/Target_predict/Target_gene_info/All_species_targetGene.txt`.
2. Filter by species, miRNA, target gene, or evidence field.
3. Join target-gene identifiers to `All_target_interproscan.tsv` or `All_target_KEGG.tsv` for functional interpretation.

### Inspect evidence for a target relationship

1. Locate the relationship in a species-specific target file or `All_species_targetGene.txt`.
2. Review the `psRNAtarget`, `psRobot`, and `Degradome` fields.
3. Use the target-gene identifier to retrieve InterProScan or KEGG annotations.
4. Consult the PmiREN3.1 methods for the thresholds and confidence definitions used in the database.

## File-format and usage notes

- Tables with `.tsv` extensions are tab-separated.
- Files named `*_targetGene.txt` are also tab-separated.
- FASTA files use standard `>` headers followed by sequence lines.
- BED files are tab-separated coordinate tables.
- `NA` denotes a missing or unavailable value.
- Species names in paths use underscores, while species names inside tables may use spaces.
- Subspecies, varieties, cultivars, and assembly-specific entries may have separate directories.
- Very large files should be processed with streaming or command-line tools rather than opened directly in spreadsheet software.
- Data availability differs among species; the absence of a file or value does not necessarily indicate a failed analysis.

## Citation

When using these data, please cite the relevant PmiREN publication and report the PmiREN database version used for the analysis.
