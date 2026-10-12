# Compact Gene-Coordinate Catalogues

These compressed TSV files retain one row per GFF3 `gene` feature:

```text
seqname  start_bp  end_bp  strand  gene_id  source_gene_id  gene_name  gene_biotype
```

They support fast interval overlap, feature labelling, and gene-ID joins while
being much smaller than the original GFF3 files. They deliberately omit
transcripts, exons, CDS features, UTRs, parent-child structure, and regulatory
features. Use the original full annotation for transcript-aware work such as
isoform quantification, coding consequence prediction, or exon overlap.

The files are sourced from local annotations supplied for package development;
their assembly and release are recorded in `manifest.tsv`. They are not a
replacement for a freshly downloaded, fully versioned annotation in a real
study.

## GO Joins

`gene_id` is a join-friendly stable ID; for Ensembl GFF3 it removes the
technical `gene:` prefix. `source_gene_id` retains the exact identifier from
the source GFF3. `gene_id` is the field normally used to join a Gene Ontology
(GO) annotation table. This package does not infer GO assignments. Obtain a
dated GAF/GPAD file and an explicit ID mapping for the same organism and
assembly, then join the mapping to `gene_id`. See `../annotations/go_sources.tsv`
for source guidance.
