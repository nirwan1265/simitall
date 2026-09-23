# Bundled example data

This directory intentionally contains only lightweight manifests and synthetic
examples. Use `fetch_ref_files()` to download full reference sequences and
`fetch_marker_panel()` to download curated marker sequences.

Downloaded genomes and generated analysis outputs should be stored outside the
installed package directory.

`panels/demo_panel.fa` is a legacy single-chromosome founder panel.
`panels/demo_multichrom_panel.fa` contains four founders across three
chromosomes using `founder|chromosome` FASTA headers.
