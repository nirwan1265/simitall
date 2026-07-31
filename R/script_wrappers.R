# Internal helpers ---------------------------------------------------------

.simitall_cli_args <- function(..., .args = list()) {
  args <- c(list(...), .args)
  out <- character()

  for (name in names(args)) {
    value <- args[[name]]
    if (is.null(value) || length(value) == 0L) {
      next
    }

    flag <- paste0("--", name)
    if (is.logical(value)) {
      if (isTRUE(value)) {
        out <- c(out, flag)
      }
    } else {
      pairs <- as.vector(rbind(rep(flag, length(value)), as.character(value)))
      out <- c(out, pairs)
    }
  }

  out
}

.simitall_run <- function(fun, args) {
  fun(args = args)
  invisible(TRUE)
}

.simitall_extdata <- function(...) {
  installed <- system.file("extdata", ..., package = "simitall")
  if (nzchar(installed)) {
    return(installed)
  }

  file.path("inst", "extdata", ...)
}

.simitall_append_lines <- function(text, path) {
  cat(paste0(text, "\n"), file = path, append = TRUE, sep = "")
}

#' Simulate a genome with repeats
#'
#' Generate a random genome or modify a reference FASTA with tandem repeats,
#' motif repeats, or both. The simulator can preserve local GC content, control
#' repeat spacing, generate multiple diverged genome copies, and write
#' coordinate truth in TSV, BED, GFF3, and JSON formats.
#'
#' @param out_fa Output FASTA path.
#' @param in_fa Optional input FASTA. Omit when `random_genome = TRUE`.
#' @param random_genome Generate a random genome instead of modifying `in_fa`.
#' @param ... Additional simulator options, such as `n_events`, `seg_len`,
#'   `copies`, `spacing_distribution`, `ploidy`, and `snp_rate`.
#'
#' @return Invisibly returns `TRUE` after writing the requested files.
#' @examples
#' out <- tempfile(fileext = ".fa")
#' simulate_genome(
#'   out_fa = out,
#'   random_genome = TRUE,
#'   random_length = 10000,
#'   n_events = 2,
#'   seg_len = 100,
#'   copies = 3,
#'   seed = 1
#' )
#' @export
simulate_genome <- function(out_fa, in_fa = NULL, random_genome = FALSE, ...) {
  extras <- list(...)
  prefix <- tools::file_path_sans_ext(out_fa)

  if (!is.null(extras$fixed_spacing) && is.null(extras$spacing_fixed)) {
    extras$spacing_fixed <- extras$fixed_spacing
    extras$fixed_spacing <- NULL
  }
  if (!is.null(extras$snp_rate) && is.null(extras$ploidy_snp_rate)) {
    extras$ploidy_snp_rate <- extras$snp_rate
    extras$snp_rate <- NULL
  }
  if (!is.null(extras$indel_rate) && is.null(extras$ploidy_indel_rate)) {
    extras$ploidy_indel_rate <- extras$indel_rate
    extras$indel_rate <- NULL
  }

  defaults <- list(
    coords_tsv = paste0(prefix, ".repeats.tsv"),
    coords_bed = paste0(prefix, ".repeats.bed"),
    coords_gff3 = paste0(prefix, ".repeats.gff3"),
    coords_tsv_per_copy = paste0(prefix, ".repeats.tsv"),
    coords_bed_per_copy = paste0(prefix, ".repeats.bed"),
    coords_gff3_per_copy = paste0(prefix, ".repeats.gff3"),
    summary_json = paste0(prefix, ".summary.json"),
    name_map_tsv = paste0(prefix, ".name_map.tsv"),
    label_blocks = TRUE
  )
  extras <- utils::modifyList(defaults, extras)

  args <- .simitall_cli_args(
    out_fa = out_fa,
    in_fa = in_fa,
    random_genome = random_genome,
    .args = extras
  )
  .simitall_run(main_01_make_tandem_repeats, args)
}

#' Generate genome annotations
#'
#' Generate gene models, CDS features, operons, TSSs, promoters, rRNA/tRNA
#' clusters, terminators, riboswitches, CRISPR arrays, GC-skew landmarks,
#' plasmids, and optional regulatory elements.
#'
#' @param in_fa Input genome FASTA.
#' @param out_prefix Prefix used to derive output paths.
#' @param ... Additional annotation options. Use `out_gff3` to override the
#'   default GFF3 path or `in_gff3` to import existing gene models.
#'
#' @return Invisibly returns `TRUE` after writing annotation outputs.
#' @examples
#' \dontrun{
#' simulate_annotations(
#'   in_fa = "genome.fa",
#'   out_prefix = "results/demo",
#'   n_genes = 100,
#'   rrna_clusters = 2,
#'   crispr_count = 1
#' )
#' }
#' @export
simulate_annotations <- function(in_fa, out_prefix, ...) {
  extras <- list(...)
  defaults <- list(
    out_gff3 = paste0(out_prefix, ".gff3"),
    genes_fa = paste0(out_prefix, ".genes.fa"),
    cds_fa = paste0(out_prefix, ".cds.fa"),
    promoters_fa = paste0(out_prefix, ".promoters.fa"),
    tss_tsv = paste0(out_prefix, ".tss.tsv"),
    operons_tsv = paste0(out_prefix, ".operons.tsv"),
    rrna_fa = paste0(out_prefix, ".rrna.fa"),
    trna_fa = paste0(out_prefix, ".trna.fa"),
    riboswitch_fa = paste0(out_prefix, ".riboswitch.fa"),
    crispr_fa = paste0(out_prefix, ".crispr.fa"),
    regulatory_fa_out = paste0(out_prefix, ".regulatory.fa"),
    plasmid_fa_out = paste0(out_prefix, ".plasmids.fa"),
    plasmid_gff3_out = paste0(out_prefix, ".plasmids.gff3")
  )
  extras <- utils::modifyList(defaults, extras)

  args <- .simitall_cli_args(
    genome_fa = in_fa,
    out_prefix = out_prefix,
    .args = extras
  )
  .simitall_run(main_01b_generate_annotations, args)
}

#' Fetch curated antibiotic marker sequences
#'
#' Download the marker accessions listed in the bundled marker manifest and
#' write them as a multi-record FASTA.
#'
#' @param out_fa Output FASTA path.
#' @param marker_tsv Marker manifest. By default, uses the manifest bundled
#'   with `simitall`.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' \dontrun{
#' fetch_marker_panel("marker_panel.fa")
#' }
#' @export
fetch_marker_panel <- function(
    out_fa,
    marker_tsv = .simitall_extdata("markers", "markers_curated.tsv")) {
  .simitall_run(main_01c_fetch_marker_panel, c(marker_tsv, out_fa))
}

#' Fetch example reference genomes
#'
#' Download one reference sequence each for E. coli K-12 MG1655, human
#' GRCh38 chromosome 1, maize B73 chromosome 1, and Arabidopsis chromosome 1.
#' Large reference genomes are downloaded on demand and are not bundled with
#' the R package.
#'
#' @param out_dir Directory in which reference FASTA files are written.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' \dontrun{
#' fetch_ref_files("references")
#' }
#' @export
fetch_ref_files <- function(out_dir = "references") {
  .simitall_run(main_01d_fetch_ref_files, out_dir)
}

#' Summarize QUAST reports
#'
#' Combine all `report.tsv` files beneath a QUAST result directory into one
#' tidy CSV containing coverage levels and key assembly metrics.
#'
#' @param quast_root Directory containing QUAST result folders.
#' @param out_csv Output CSV path.
#'
#' @return The summary data frame, invisibly.
#' @examples
#' \dontrun{
#' summarize_quast("results/quast", "results/quast_summary.csv")
#' }
#' @export
summarize_quast <- function(quast_root, out_csv) {
  reports <- list.files(
    quast_root,
    pattern = "report\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(reports)) {
    stop("No report.tsv files found under: ", quast_root)
  }

  rows <- lapply(reports, function(path) {
    report <- utils::read.delim(
      path,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    values <- stats::setNames(report[[2L]], report[[1L]])
    get_metric <- function(name) {
      if (name %in% names(values)) as.numeric(values[[name]]) else NA_real_
    }

    folder <- basename(dirname(path))
    match <- regexec("^ill([0-9]+)_pb([0-9]+)$", folder)
    coverage <- regmatches(folder, match)[[1L]]

    data.frame(
      ill_cov = if (length(coverage)) as.integer(coverage[2L]) else NA_integer_,
      pb_cov = if (length(coverage)) as.integer(coverage[3L]) else NA_integer_,
      total_length = get_metric("Total length"),
      n_contigs = get_metric("# contigs"),
      n50 = get_metric("N50"),
      largest_contig = get_metric("Largest contig"),
      gc = get_metric("GC (%)"),
      genome_fraction = get_metric("Genome fraction (%)"),
      mismatches_per_100kbp = get_metric("# mismatches per 100 kbp"),
      indels_per_100kbp = get_metric("# indels per 100 kbp"),
      misassemblies = get_metric("# misassemblies"),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  result <- result[order(result$ill_cov, result$pb_cov), , drop = FALSE]
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, out_csv, row.names = FALSE)
  invisible(result)
}

#' Simulate a GWAS cohort
#'
#' Simulate variants, LD blocks, recombination, population structure, sample
#' genotypes, and optional phenotypes from a reference genome.
#'
#' @param genome_fa Reference genome FASTA.
#' @param out_prefix Prefix for VCF, genotype, phenotype, and QC outputs.
#' @param ... Additional options such as `n_samples`, `snp_rate`,
#'   `ld_block_size`, `n_subpops`, and `recomb_map`.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' \dontrun{
#' simulate_gwas_cohort(
#'   genome_fa = "genome.fa",
#'   out_prefix = "results/gwas",
#'   n_samples = 200,
#'   n_subpops = 2
#' )
#' }
#' @export
simulate_gwas_cohort <- function(genome_fa, out_prefix, ...) {
  args <- .simitall_cli_args(
    genome_fa = genome_fa,
    out_prefix = out_prefix,
    .args = list(...)
  )
  .simitall_run(main_10_simulate_gwas_cohort, args)
}

#' Simulate phenotypes from genotypes
#'
#' Simulate quantitative, binary, multi-trait, dominance, epistatic, and
#' gene-by-environment phenotype architectures using `simplePHENOTYPES`.
#'
#' @param geno_file Input VCF, HapMap, or numeric genotype file.
#' @param out_prefix Prefix for phenotype and simulation-summary outputs.
#' @param ... Additional phenotype options such as `h2`, `n_add_qtn`,
#'   `n_traits`, `architecture`, and `binary`.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' \dontrun{
#' simulate_phenotypes(
#'   geno_file = "cohort.vcf",
#'   out_prefix = "results/traits",
#'   h2 = 0.5,
#'   n_add_qtn = 20
#' )
#' }
#' @export
simulate_phenotypes <- function(geno_file, out_prefix, ...) {
  args <- .simitall_cli_args(
    geno_file = geno_file,
    out_prefix = out_prefix,
    .args = list(...)
  )
  .simitall_run(main_10b_simulate_phenotypes, args)
}

#' Generate a random haplotype panel
#'
#' Create an aligned multi-FASTA panel by introducing SNPs and indels into a
#' shared synthetic founder sequence.
#'
#' @param out_fa Output multi-FASTA path.
#' @param n_haplotypes Number of founder haplotypes.
#' @param length Sequence length in base pairs.
#' @param ... Additional options such as `gc`, `snp_rate`, `indel_rate`, and
#'   `seed`.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' panel <- tempfile(fileext = ".fa")
#' generate_random_haplotype_panel(panel, n_haplotypes = 4, length = 2000)
#' @export
generate_random_haplotype_panel <- function(
    out_fa,
    n_haplotypes = 8,
    length = 100000,
    ...) {
  args <- .simitall_cli_args(
    out_fa = out_fa,
    n_haps = n_haplotypes,
    length = length,
    .args = list(...)
  )
  .simitall_run(main_11a_generate_random_haplotype_panel, args)
}

#' Simulate breeding populations
#'
#' Simulate F1, F2, backcross, selfing, RIL, NIL, doubled-haploid, NAM, and
#' MAGIC populations from an aligned haplotype FASTA panel. Optional models
#' include crossover interference, fixed loci, background selection,
#' segregation distortion, structural variants, missingness, and genotyping
#' error.
#'
#' @param haplotype_fa Aligned founder haplotypes in multi-FASTA format.
#' @param out_prefix Prefix for FASTA, metadata, VCF, graph, and QC outputs.
#' @param ... Additional breeding options such as `scheme`, `n_offspring`,
#'   `generations`, `parents`, `founders`, and `fix_locus`.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' \dontrun{
#' panel <- system.file("extdata", "panels", "demo_panel.fa",
#'                      package = "simitall")
#' simulate_breeding(
#'   haplotype_fa = panel,
#'   out_prefix = "results/ril",
#'   scheme = "RIL",
#'   n_offspring = 100,
#'   generations = 6
#' )
#' }
#' @export
simulate_breeding <- function(haplotype_fa, out_prefix, ...) {
  args <- .simitall_cli_args(
    haplotype_fa = haplotype_fa,
    out_prefix = out_prefix,
    .args = list(...)
  )
  .simitall_run(main_11_simulate_breeding, args)
}

#' Run a SimuPOP simulation from R
#'
#' Run selected SimuPOP mating schemes through `reticulate`, including
#' built-in F1, F2, NAM, and MAGIC presets and direct mating-scheme
#' configurations. Genotypes and sample metadata are exported automatically.
#'
#' @param out_prefix Prefix for genotype and metadata outputs.
#' @param config A configuration list or path to a JSON configuration file.
#'   If omitted, named values supplied through `...` form the configuration.
#' @param ... Configuration values added to or overriding `config`.
#'
#' @return Invisibly returns `TRUE`.
#' @examples
#' \dontrun{
#' simupop_api(
#'   out_prefix = "results/simupop_f2",
#'   config = list(
#'     population = list(size = 100, ploidy = 2, loci = 200),
#'     preset = "F2",
#'     generations = 2
#'   )
#' )
#' }
#' @export
simupop_api <- function(out_prefix, config = NULL, ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for simupop_api().")
  }

  overrides <- list(...)
  remove_config <- FALSE

  if (is.character(config) && length(config) == 1L) {
    if (!file.exists(config)) {
      stop("SimuPOP config file not found: ", config)
    }
    if (length(overrides)) {
      config <- utils::modifyList(jsonlite::fromJSON(config), overrides)
    } else {
      config_path <- config
    }
  } else {
    if (is.null(config)) {
      config <- list()
    }
    if (!is.list(config)) {
      stop("'config' must be a list or a JSON file path.")
    }
    config <- utils::modifyList(config, overrides)
  }

  if (!exists("config_path", inherits = FALSE)) {
    config_path <- tempfile("simitall-simupop-", fileext = ".json")
    jsonlite::write_json(config, config_path, auto_unbox = TRUE, pretty = TRUE)
    remove_config <- TRUE
  }
  if (remove_config) {
    on.exit(unlink(config_path), add = TRUE)
  }

  args <- c("--config", config_path, "--out_prefix", out_prefix)
  .simitall_run(main_12_simupop_api, args)
}
