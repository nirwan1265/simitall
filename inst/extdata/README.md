# Bundled example data

This directory intentionally contains only lightweight manifests and synthetic
examples. Use `fetch_ref_files()` to download full reference sequences and
`fetch_marker_panel()` to download curated marker sequences.

Downloaded genomes and generated analysis outputs should be stored outside the
installed package directory.

`panels/demo_panel.fa` is a legacy single-chromosome founder panel.
`panels/demo_multichrom_panel.fa` contains four founders across three
chromosomes using `founder|chromosome` FASTA headers.

## Runnable Bundled Demos

The files below are deliberately small and entirely synthetic. They exist so a
new installation can execute an end-to-end example without downloading a
reference panel. They are not real organism data and must not be presented as
maize NAM founders, 1000 Genomes participants, or clinical evidence.

- `panels/demo_maize_nam_chr10.fa`: eight aligned synthetic founders labelled
  `hap1` through `hap8` on a toy `chr10` interval.
- `maps/demo_maize_nam_chr10_map.tsv`: compatible toy recombination map with a
  lower-recombination central interval.
- `human_irf6/demo_human_irf6_chr1.vcf.gz`: tiny synthetic GRCh38-like chr1
  marker panel near the IRF6-region demonstration coordinates.
- `human_irf6/demo_human_irf6_chr1.gff3`: matching synthetic annotation.
- `human_irf6/demo_human_irf6_chr1_recombination_map.tsv`: matching toy map.

## Gene-Annotation Examples

`annotations/` contains tiny coordinate-only GFF3 excerpts for maize,
Arabidopsis, wheat, rice, sorghum, and Drosophila, plus a manifest pointing to
the region-specific human IRF6 example. They make annotation parsing and
workflow demonstrations runnable without turning the package into a
multi-gigabyte genome-data download. See `annotations/README.md` for source
assemblies, boundaries, and the separate Gene Ontology (GO) download policy.

`gene_coordinates/` contains compact, genome-wide gene-only coordinate tables
for Arabidopsis, barley, rice, wheat, maize, and Drosophila. Each row retains
the chromosome/contig, start, end, strand, gene ID, gene name, and biotype but
not transcript/exon/CDS detail. These are appropriate for interval overlap and
gene-ID joins; see `gene_coordinates/README.md` for their assembly manifest and
GO join guidance.

`go/arabidopsis_tair10_go_demo.gaf.gz` is a small Gene Ontology Annotation File
(GAF) subset matching the Arabidopsis GFF3 excerpt. The functions
`read_go_gaf()`, `join_gene_go()`, and `summarize_go_terms()` also accept a
user's full local GAF download directly.

In the Shiny agent, open **Optional input details**, choose **Use package demo
only**, then ask for a maize NAM or human IRF6 pedigree demonstration. The
returned recipe is runnable with the **RUN** button.
