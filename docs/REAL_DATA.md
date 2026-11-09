# Real-Data Import Examples

SimItAll simulations retain programmed truth. Public datasets exercise a
separate, complementary contract: the application must import authentic inputs,
preserve their provenance, apply the correct analysis, and report input/QC
checks without inventing a biological truth set.

The curated download records are in
`inst/agent/capabilities/real_data_catalog.json`. They describe sources,
expected file layouts, annotation requirements, and checks; they do **not**
download data or make an external dataset part of the R package.

## Storage policy

Keep downloaded public data outside the Git repository, for example:

```text
~/Documents/simitall_real_data/
  rnaseq_airway/
  rnaseq_maize_drought/
  chipseq_encode_k562/
  gwas_maize_nam/
  gwas_arabidopsis/
  scrnaseq_pbmc3k/
```

Copy `data-raw/real_data_manifests/manifest.template.json` into each dataset
folder as `manifest.json` and complete it after downloading. The manifest is
the provenance record: accession, date, source URL, exact reference assembly,
annotation release, identifiers, and any sample-ID mapping.

Do not commit downloaded counts, VCFs, BAMs, FASTQs, whole annotations, or GO
files. Commit only stable catalogue metadata, import code, small test fixtures,
and non-sensitive example manifests.

## Gene-information contract

A JSON manifest is deliberately not a replacement for a GFF3/GTF. It points to
the real annotation resource while preserving machine-readable compatibility
rules.

For a real analysis, retain:

- A compatible GFF3, GTF, or gene-coordinate TSV with `gene_id`, chromosome,
  start, end, and strand.
- The exact genome assembly/reference build used by the counts, peaks, or VCF.
- The gene-ID namespace and an explicit mapping, when the data and annotation
  use different identifiers.
- A dated organism-specific Gene Ontology Annotation File (GAF) plus a
  reviewed ID mapping if GO work is requested.

Before any gene annotation is reported, SimItAll must verify assembly,
chromosome-label, coordinate, and identifier compatibility. It must not assume
that an annotation from a different maize, Arabidopsis, or human reference is
interchangeable.

## Initial public examples

| Identifier | Analysis | First download | Annotation rule |
| --- | --- | --- | --- |
| `airway_dexamethasone` | Human bulk RNA-seq DE | Bioconductor `airway` counts and sample metadata | Confirm the original gene model/build before interval or GO interpretation. |
| `maize_drought_rnaseq_gse132113` | Maize drought RNA-seq | GEO supplementary data for `GSE132113` | Inspect whether values are counts or normalized expression, then obtain a matching maize annotation. |
| `encode_k562_h3k27ac_encsr000akp` | Human processed ChIP-seq peaks | ENCODE `ENCSR000AKP` narrowPeak, bigWig, metadata | Read the exact reference assembly from ENCODE file metadata before peak-to-gene annotation. |
| `maize_nam_panzea_gwas` | Maize NAM GWAS | Panzea phenotype plus MaizeGDB/SNPVersity VCF subset | Use an annotation matching the VCF reference build. |
| `arabidopsis_1001genomes_arapheno_gwas` | Arabidopsis GWAS | 1001 Genomes VCF subset plus AraPheno phenotype | Confirm TAIR/reference compatibility and accession-ID mappings. |
| `pbmc3k_10x` | Human single-cell RNA-seq | 10x `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz` | Features support clustering; add a compatible GTF for interval/GO interpretation. |

The catalogue records exact provider URLs and workflow-specific import checks.
Use processed data first. FASTQ/BAM imports are a separate raw-read route and
must not be represented as a completed analysis merely because a file was
downloaded.

## Reproducible import checks

Each real-data import should at least confirm:

1. All declared local files exist and their checksums/provenance are recorded.
2. Sample identifiers agree across all matrices, VCFs, metadata, and covariates.
3. The assay contains the expected data type: integer counts, sparse UMI matrix,
   genomic intervals, or valid genotype calls.
4. Reference assembly, contig names, and gene-ID namespace are compatible with
   the selected annotation.
5. The requested analysis produces its standard QC/results, while claims of
   simulated truth recovery are omitted for real data.
