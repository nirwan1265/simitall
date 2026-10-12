# Bundled GO Demonstration Data

`arabidopsis_tair10_go_demo.gaf.gz` is a small, verbatim subset of the local
GO Consortium `ARATH-mod.gaf.gz` download. It includes annotations for the 25
Arabidopsis genes in `../annotations/arabidopsis_tair10_gene_slice.gff3` and
exists only to demonstrate GAF parsing and coordinate-to-GO joins.

It is not a complete or current annotation set. For actual analysis, pass a
locally downloaded, dated organism-specific GAF to `read_go_gaf()` and record
the GAF header date, ontology release, and ID-mapping choices.

The supplied local GAF downloads can be used directly without copying their
full contents into this package. `read_go_gaf()` accepts `.gaf.gz` files.
