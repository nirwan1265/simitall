# Local GO Annotation Downloads

This directory is a source-controlled copy of the dated GAF downloads supplied
in `Downloads/GO_annotation` on 2026-10-07. These are full Gene Ontology (GO)
annotation resources, kept in the GitHub repository as raw study resources.
They are intentionally excluded from the built R package so that installing
`simitall` does not download approximately 88 MB of mutable external data.

Use a file directly with the package API:

```r
annotations <- read_go_gaf("resources/go_annotations/ARATH-mod.gaf.gz")
genes <- read_gene_coordinates("Arabidopsis_thaliana")
gene_go <- join_gene_go(genes, annotations)
```

`manifest.tsv` retains checksums and GAF header dates. Preserve those fields in
any analysis output. A GAF maps identifiers to GO terms, evidence codes, and
references; it does not establish a gene's experimental function by itself.

The source files remain subject to the Gene Ontology Consortium's terms and
release provenance. The bundled package demonstration subset lives separately
at `inst/extdata/go/`.
