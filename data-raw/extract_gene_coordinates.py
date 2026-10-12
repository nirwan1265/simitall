#!/usr/bin/env python3
"""Create compact, gene-only coordinate tables from GFF3/GTF annotations.

Usage:
  python3 data-raw/extract_gene_coordinates.py input.gff3 output.tsv.gz

The output is deliberately limited to stable identifiers plus genomic interval
metadata. It is useful for interval overlap and GO-ID mapping, but is not a
replacement for a full GFF3/GTF when transcript or exon structure is needed.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import re


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def attr_value(attributes: str, keys: tuple[str, ...]) -> str:
    values = {}
    for field in attributes.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
        else:
            continue
        values[key.strip()] = value.strip().strip('"')
    for key in keys:
        if key in values and values[key]:
            return values[key]
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_gff", type=Path)
    parser.add_argument("output_tsv_gz", type=Path)
    args = parser.parse_args()

    args.output_tsv_gz.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open_text(args.input_gff) as source, gzip.open(args.output_tsv_gz, "wt") as target:
        target.write(
            "seqname\tstart_bp\tend_bp\tstrand\tgene_id\tsource_gene_id\tgene_name\tgene_biotype\n"
        )
        for line in source:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attributes = parts[8]
            source_gene_id = attr_value(attributes, ("ID", "gene_id", "gene"))
            # Ensembl GFF3 commonly prefixes stable IDs with "gene:". Keep the
            # original identifier too, but provide the join-friendly form first.
            gene_id = re.sub(r"^gene:", "", source_gene_id)
            gene_name = attr_value(attributes, ("Name", "gene_name", "gene_symbol"))
            gene_biotype = attr_value(attributes, ("biotype", "gene_biotype", "gene_type"))
            target.write("\t".join((
                parts[0], parts[3], parts[4], parts[6], gene_id, source_gene_id,
                gene_name, gene_biotype
            )) + "\n")
            count += 1
    print(f"Wrote {count} gene intervals to {args.output_tsv_gz}")


if __name__ == "__main__":
    main()
