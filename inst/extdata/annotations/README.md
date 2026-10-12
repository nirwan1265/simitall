# Packaged Gene-Annotation Slices

This directory provides small, inspectable GFF3 examples for exercising
annotation-aware SimItAll workflows. They are **not** whole-genome annotation
sets and are not appropriate for biological interpretation, gene discovery, or
production enrichment analysis.

Each plant and fly file contains 25 `gene` records only. The records retain
their original sequence names, coordinates, strands, and source attributes but
intentionally omit transcript, exon, CDS, UTR, and regulatory records. The
human IRF6 demonstration annotation remains in `../human_irf6/` because it is
a region-specific toy study input.

Use `annotation_manifest.tsv` to identify the assembly, source release, and
the canonical source. A real study must supply a compatible full GFF3/GTF and
reference genome build.

## GO Annotations

Gene coordinates and Gene Ontology (GO) annotations are different resources:

- GFF3/GTF describes genomic features and their coordinates.
- A GO Annotation File (GAF) maps a gene product identifier to one or more GO
  terms, evidence codes, and supporting references.

Do not bundle a mutable whole-organism GO file in the package. Download a
dated GAF and the matching GO ontology for the organism and gene-ID namespace
in a real analysis. The `go_sources.tsv` catalogue records recommended
authoritative sources, but contains no inferred GO assignments.
