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

In the Shiny agent, open **Optional input details**, choose **Use package demo
only**, then ask for a maize NAM or human IRF6 pedigree demonstration. The
returned recipe is runnable with the **RUN** button.
