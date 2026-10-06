#' List validated sequencing-platform profiles
#'
#' The profiles connect plain-language platform names to simulator settings.
#' ART does not provide current error models for every current Illumina system,
#' so its system profile is an explicitly labelled compatible proxy rather than
#' a claim that it reproduces an instrument's current chemistry.
#'
#' @return A data frame describing supported Illumina and PacBio profiles.
#' @export
simitall_sequencing_profiles <- function() {
  data.frame(
    platform = c(
      "MiSeq", "MiSeq i100", "NextSeq 1000/2000", "NovaSeq X",
      "Sequel IIe HiFi", "Revio HiFi", "Vega HiFi"
    ),
    technology = c(rep("Illumina short read", 4L), rep("PacBio HiFi", 3L)),
    simulator = c(rep("ART", 4L), rep("PBSIM/PBSIM3", 3L)),
    simulator_profile = c("MSv3", "MSv3", "NS50", "HS25", rep("HIFI", 3L)),
    default_read_length_bp = c(250L, 150L, 150L, 150L, NA_integer_, NA_integer_, NA_integer_),
    maximum_read_length_bp = c(300L, 300L, 300L, 150L, NA_integer_, NA_integer_, NA_integer_),
    simulation_note = c(
      "ART MiSeq v3 error-profile proxy.",
      "ART has no MiSeq i100 model; uses the MiSeq v3 proxy.",
      "ART NextSeq 500 profile proxy for NextSeq 1000/2000.",
      "ART has no NovaSeq X model; uses the HiSeq 2500 proxy.",
      "HiFi read type; supply an empirical PBSIM model when instrument-specific errors matter.",
      "HiFi read type; supply an empirical PBSIM model when instrument-specific errors matter.",
      "HiFi read type; supply an empirical PBSIM model when instrument-specific errors matter."
    ),
    stringsAsFactors = FALSE
  )
}

.simitall_normalize_sequencing_platform <- function(platform) {
  key <- gsub("[^a-z0-9]+", "", tolower(platform))
  aliases <- c(
    generic = "generic_illumina", genericillumina = "generic_illumina", illumina = "generic_illumina",
    miseq = "miseq", miseqv3 = "miseq",
    miseqi100 = "miseq_i100", i100 = "miseq_i100", i100plus = "miseq_i100",
    nextseq = "nextseq_1000_2000", nextseq1000 = "nextseq_1000_2000", nextseq2000 = "nextseq_1000_2000", nextseq10002000 = "nextseq_1000_2000",
    novaseq = "novaseq_x", novaseqx = "novaseq_x", novaseqxplus = "novaseq_x",
    sequeliie = "sequel_iie_hifi", sequeliiehifi = "sequel_iie_hifi", revio = "revio_hifi", reviohifi = "revio_hifi", vega = "vega_hifi", vegahifi = "vega_hifi",
    hifi = "generic_hifi", generichifi = "generic_hifi", clr = "clr"
  )
  if (!key %in% names(aliases)) {
    stop("Unsupported sequencing platform: ", platform,
         ". Use simitall_sequencing_profiles() to view supported profiles.")
  }
  unname(aliases[[key]])
}

.simitall_illumina_profile <- function(platform) {
  platform <- .simitall_normalize_sequencing_platform(platform)
  profiles <- list(
    generic_illumina = list(art_system = "HS25", readlen = 150L, max_readlen = 150L, insert = 350L, insert_sd = 50L),
    miseq = list(art_system = "MSv3", readlen = 250L, max_readlen = 300L, insert = 550L, insert_sd = 60L),
    miseq_i100 = list(art_system = "MSv3", readlen = 150L, max_readlen = 300L, insert = 350L, insert_sd = 50L),
    nextseq_1000_2000 = list(art_system = "NS50", readlen = 150L, max_readlen = 300L, insert = 350L, insert_sd = 50L),
    novaseq_x = list(art_system = "HS25", readlen = 150L, max_readlen = 150L, insert = 350L, insert_sd = 50L)
  )
  if (!platform %in% names(profiles)) stop("platform is not an Illumina profile")
  c(list(key = platform), profiles[[platform]])
}

.simitall_pacbio_profile <- function(platform) {
  platform <- .simitall_normalize_sequencing_platform(platform)
  profiles <- list(
    clr = list(type = "CLR", label = "CLR"),
    generic_hifi = list(type = "HIFI", label = "HiFi"),
    sequel_iie_hifi = list(type = "HIFI", label = "Sequel_IIe_HiFi"),
    revio_hifi = list(type = "HIFI", label = "Revio_HiFi"),
    vega_hifi = list(type = "HIFI", label = "Vega_HiFi")
  )
  if (!platform %in% names(profiles)) stop("platform is not a PacBio profile")
  c(list(key = platform), profiles[[platform]])
}

#' Simulate Illumina reads with ART
#'
#' Wrapper around the `art_illumina` command for paired-end read simulation.
#' This function runs ART and compresses output FASTQ files with pigz if available.
#'
#' @param ref_fa Character. Reference FASTA path.
#' @param outprefix Character. Output prefix (ART will write `<outprefix>1.fq` and `<outprefix>2.fq`).
#' @param cov Numeric. Fold coverage.
#' @param platform Illumina profile: `"generic_illumina"`, `"MiSeq"`,
#'   `"MiSeq i100"`, `"NextSeq 1000/2000"`, or `"NovaSeq X"`.
#' @param readlen Optional read length. Defaults depend on `platform` and are
#'   checked against the profile's supported maximum.
#' @param ins Optional mean insert size; defaults depend on `platform`.
#' @param sd Optional insert-size standard deviation; defaults depend on
#'   `platform`.
#'
#' @return A list with `r1` and `r2` file paths to gzipped FASTQ files.
#' @examples
#' \dontrun{
#' sim_illumina_art("01_simref/ecoli_repMed.fa", "02_reads/ecoli/illumina/ecoli.ill_cov30_", 30)
#' }
#' @export
sim_illumina_art <- function(
    ref_fa,
    outprefix,
    cov,
    readlen = NULL,
    ins = NULL,
    sd = NULL,
    platform = "generic_illumina") {
  if (!file.exists(ref_fa)) stop("Reference FASTA not found: ", ref_fa)
  profile <- .simitall_illumina_profile(platform)
  readlen <- if (is.null(readlen)) profile$readlen else as.integer(readlen)
  ins <- if (is.null(ins)) profile$insert else as.integer(ins)
  sd <- if (is.null(sd)) profile$insert_sd else as.integer(sd)
  if (!is.finite(readlen) || readlen < 1L || readlen > profile$max_readlen) {
    stop("readlen must be between 1 and ", profile$max_readlen,
         " bp for platform ", platform)
  }
  if (!is.finite(ins) || !is.finite(sd) || ins < 1L || sd < 0L) {
    stop("ins must be positive and sd must be non-negative")
  }
  if (!nzchar(Sys.which("art_illumina"))) {
    stop("ART was not found in PATH. Install the Bioconda package 'art'.")
  }
  dir.create(dirname(outprefix), showWarnings = FALSE, recursive = TRUE)
  cmd <- sprintf(
    "art_illumina -ss %s -i %s -p -l %s -f %s -m %s -s %s -o %s",
    profile$art_system, shQuote(ref_fa), readlen, cov, ins, sd, shQuote(outprefix)
  )
  status <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (status != 0) stop("ART simulation failed")

  r1 <- paste0(outprefix, "1.fq")
  r2 <- paste0(outprefix, "2.fq")
  if (file.exists(r1)) system(sprintf("pigz -f %s", shQuote(r1)))
  if (file.exists(r2)) system(sprintf("pigz -f %s", shQuote(r2)))
  list(r1 = paste0(r1, ".gz"), r2 = paste0(r2, ".gz"))
}

#' Simulate PacBio reads with PBSIM/PBSIM3
#'
#' Wrapper around `pbsim` or `pbsim3` for CLR/HiFi simulation.
#'
#' @param ref_fa Character. Reference FASTA path.
#' @param outdir Character. Output directory.
#' @param cov Numeric. Fold coverage.
#' @param type Deprecated compatibility argument: `"CLR"` or `"HIFI"`.
#' @param platform PacBio profile: `"CLR"`, `"HiFi"`, `"Sequel IIe HiFi"`,
#'   `"Revio HiFi"`, or `"Vega HiFi"`.
#' @param seed Integer. Random seed (default 1).
#'
#' @return Character. Output directory path.
#' @examples
#' \dontrun{
#' sim_pacbio("01_simref/ecoli_repMed.fa", "02_reads/ecoli/pacbio_HIFI/cov20", 20, type = "HIFI")
#' }
#' @export
sim_pacbio <- function(ref_fa, outdir, cov, type = NULL, seed = 1, platform = NULL) {
  if (!file.exists(ref_fa)) stop("Reference FASTA not found: ", ref_fa)
  ref_fa <- normalizePath(ref_fa, mustWork = TRUE)
  outdir <- normalizePath(
    outdir,
    mustWork = FALSE
  )
  if (is.null(platform)) platform <- if (is.null(type)) "HiFi" else type
  profile <- .simitall_pacbio_profile(platform)
  type <- profile$type
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  has_pbsim3 <- nzchar(Sys.which("pbsim3"))
  has_pbsim <- nzchar(Sys.which("pbsim"))

  if (has_pbsim3) {
    prefix <- paste0("pb", profile$label)
    cmd <- sprintf("pbsim3 --depth %s --seed %s --prefix %s_cov%s %s", cov, seed, prefix, cov, shQuote(ref_fa))
  } else if (has_pbsim) {
    dtype <- ifelse(type == "HIFI", "CCS", "CLR")
    prefix <- paste0("pb", profile$label)
    cmd <- sprintf("pbsim --depth %s --seed %s --data-type %s --prefix %s_cov%s %s", cov, seed, dtype, prefix, cov, shQuote(ref_fa))
  } else {
    stop("pbsim3 or pbsim not found in PATH")
  }

  old <- getwd()
  setwd(outdir)
  on.exit(setwd(old), add = TRUE)
  status <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (status != 0) stop("PacBio simulation failed")

  fq <- list.files(".", pattern = "\\.(fastq|fq)$", full.names = TRUE)
  if (length(fq) > 0) {
    for (f in fq) system(sprintf("pigz -f %s", shQuote(f)))
  }
  invisible(outdir)
}

#' Simulate Oxford Nanopore reads with Badread
#'
#' Simulate long reads using predefined identity profiles for older R9.4,
#' R9.4.1, modern R10.4.1, or perfect reads.
#'
#' @param ref_fa Reference FASTA path.
#' @param outdir Output directory.
#' @param cov Fold coverage.
#' @param model Identity preset: `"nanopore2018"`, `"nanopore2020"`,
#'   `"nanopore2023"`, or `"perfect"`.
#' @param mean_len Mean read length.
#' @param len_sd Read-length standard deviation.
#' @param identity Optional Badread identity string in
#'   `"mean,max,standard_deviation"` format.
#' @param seed Random seed.
#'
#' @return The gzipped FASTQ path, invisibly.
#' @examples
#' \dontrun{
#' sim_nanopore("genome.fa", "results/nanopore", cov = 20)
#' }
#' @export
sim_nanopore <- function(
    ref_fa,
    outdir,
    cov,
    model = "nanopore2023",
    mean_len = 15000,
    len_sd = 13000,
    identity = NULL,
    seed = 1) {
  if (!file.exists(ref_fa)) {
    stop("Reference FASTA not found: ", ref_fa)
  }
  if (!nzchar(Sys.which("badread"))) {
    stop("Badread was not found in PATH. Install the package 'badread'.")
  }

  presets <- c(
    nanopore2018 = "87,95,5",
    nanopore2020 = "92,98,3",
    nanopore2023 = "98,100,1",
    perfect = "100,100,0"
  )
  model <- match.arg(tolower(model), names(presets))
  if (is.null(identity)) {
    identity <- unname(presets[[model]])
  }

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out_fq <- file.path(outdir, sprintf("nanopore_%s_cov%s.fastq", model, cov))
  command <- sprintf(
    "badread simulate --reference %s --quantity %sx --length %s,%s --identity %s --seed %s > %s",
    shQuote(normalizePath(ref_fa, mustWork = TRUE)),
    cov,
    mean_len,
    len_sd,
    identity,
    seed,
    shQuote(out_fq)
  )
  status <- system(command)
  if (status != 0L) {
    stop("Nanopore simulation failed.")
  }

  if (nzchar(Sys.which("pigz"))) {
    status <- system2("pigz", c("-f", out_fq))
  } else {
    status <- system2("gzip", c("-f", out_fq))
  }
  if (status != 0L) {
    stop("Could not compress Nanopore FASTQ output.")
  }

  invisible(paste0(out_fq, ".gz"))
}

#' Run Illumina + PacBio grid simulation
#'
#' Simulates read sets across fixed Illumina and PacBio coverage grids.
#'
#' @param simref_fa Character. Reference FASTA.
#' @param tag Character. Tag used for output subfolders.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_grid_both("genome.fa", "ecoli_demo")
#' }
#' @param ill_covs Illumina coverage values.
#' @param pb_covs PacBio coverage values.
#' @param output_dir Root output directory.
#' @param seed PacBio random seed.
#' @export
run_grid_both <- function(
    simref_fa,
    tag,
    ill_covs = c(10, 20, 30, 40, 60, 80),
    pb_covs = c(5, 10, 15, 20, 30, 40),
    output_dir = "reads",
    seed = 1) {

  dir.create(file.path(output_dir, tag, "illumina"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, tag, "pacbio_CLR"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, tag, "pacbio_HIFI"), showWarnings = FALSE, recursive = TRUE)

  for (ic in ill_covs) {
    outp <- file.path(output_dir, tag, "illumina", paste0(tag, ".ill_cov", ic, "_"))
    if (!file.exists(paste0(outp, "1.fq.gz"))) {
      sim_illumina_art(simref_fa, outp, ic, 150, 350, 50)
    }
  }

  for (pc in pb_covs) {
    outd <- file.path(output_dir, tag, "pacbio_CLR", paste0("cov", pc))
    if (!length(list.files(outd, pattern = paste0("pbCLR_cov", pc, ".*fastq.gz"), full.names = TRUE))) {
      sim_pacbio(simref_fa, outd, pc, "CLR", seed)
    }
  }

  for (pc in pb_covs) {
    outd <- file.path(output_dir, tag, "pacbio_HIFI", paste0("cov", pc))
    if (!length(list.files(outd, pattern = paste0("pbHIFI_cov", pc, ".*fastq.gz"), full.names = TRUE))) {
      sim_pacbio(simref_fa, outd, pc, "HIFI", seed)
    }
  }

  TRUE
}

#' Run Unicycler hybrid assemblies across a grid
#'
#' @param tag Character. Tag used for input/output folders.
#' @param threads Integer. Number of threads.
#' @param reads_dir Root directory containing simulated reads.
#' @param output_dir Root assembly output directory.
#' @param ill_covs Illumina coverage values.
#' @param pb_covs PacBio coverage values.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_unicycler_grid("ecoli_repMed", threads = 8)
#' }
#' @export
run_unicycler_grid <- function(
    tag,
    threads = 8,
    reads_dir = "reads",
    output_dir = "assemblies",
    ill_covs = c(10, 20, 30, 40, 60, 80),
    pb_covs = c(5, 10, 15, 20, 30, 40)) {
  if (!nzchar(Sys.which("unicycler"))) {
    stop("Unicycler was not found in PATH.")
  }

  for (ic in ill_covs) {
    r1 <- file.path(reads_dir, tag, "illumina", paste0(tag, ".ill_cov", ic, "_1.fq.gz"))
    r2 <- file.path(reads_dir, tag, "illumina", paste0(tag, ".ill_cov", ic, "_2.fq.gz"))

    for (pc in pb_covs) {
      pb_dir <- file.path(reads_dir, tag, "pacbio_CLR", paste0("cov", pc))
      pb_fastq <- list.files(pb_dir, pattern = paste0("pbCLR_cov", pc, ".*fastq.gz"), full.names = TRUE)
      if (length(pb_fastq) > 0) {
        out_dir <- file.path(output_dir, tag, "CLR", paste0("ill", ic, "_pb", pc))
        if (!file.exists(file.path(out_dir, "assembly.fasta"))) {
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          cmd <- sprintf("unicycler -1 %s -2 %s -l %s -o %s -t %s --mode normal",
                         shQuote(r1), shQuote(r2), shQuote(pb_fastq[1]), shQuote(out_dir), threads)
          if (system(cmd) != 0L) stop("Unicycler failed for ", basename(out_dir))
        }
      }

      pb_dir <- file.path(reads_dir, tag, "pacbio_HIFI", paste0("cov", pc))
      pb_fastq <- list.files(pb_dir, pattern = paste0("pbHIFI_cov", pc, ".*fastq.gz"), full.names = TRUE)
      if (length(pb_fastq) > 0) {
        out_dir <- file.path(output_dir, tag, "HIFI", paste0("ill", ic, "_pb", pc))
        if (!file.exists(file.path(out_dir, "assembly.fasta"))) {
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          cmd <- sprintf("unicycler -1 %s -2 %s -l %s -o %s -t %s --mode normal",
                         shQuote(r1), shQuote(r2), shQuote(pb_fastq[1]), shQuote(out_dir), threads)
          if (system(cmd) != 0L) stop("Unicycler failed for ", basename(out_dir))
        }
      }
    }
  }
  TRUE
}

#' Run QUAST evaluation across assembly grid
#'
#' @param tag Character. Tag used for input/output folders.
#' @param ref_fa Truth reference FASTA.
#' @param assemblies_dir Root assembly directory.
#' @param output_dir Root QUAST output directory.
#' @param threads Number of QUAST threads.
#' @param ill_covs Illumina coverage values.
#' @param pb_covs PacBio coverage values.
#'
#' @return TRUE on success.
#' @examples
#' \dontrun{
#' run_quast_grid("ecoli_demo", ref_fa = "genome.fa")
#' }
#' @export
run_quast_grid <- function(
    tag,
    ref_fa,
    assemblies_dir = "assemblies",
    output_dir = "evaluation",
    threads = 4,
    ill_covs = c(10, 20, 30, 40, 60, 80),
    pb_covs = c(5, 10, 15, 20, 30, 40)) {
  quast <- Sys.which("quast.py")
  if (!nzchar(quast)) {
    quast <- Sys.which("quast")
  }
  if (!nzchar(quast)) {
    stop("QUAST was not found in PATH.")
  }
  if (!file.exists(ref_fa)) {
    stop("Reference FASTA not found: ", ref_fa)
  }

  for (mode in c("CLR", "HIFI")) {
    outroot <- file.path(output_dir, tag, mode)
    dir.create(outroot, showWarnings = FALSE, recursive = TRUE)

    for (ic in ill_covs) {
      for (pc in pb_covs) {
        asm <- file.path(assemblies_dir, tag, mode, paste0("ill", ic, "_pb", pc), "assembly.fasta")
        if (!file.exists(asm)) next
        outdir <- file.path(outroot, paste0("ill", ic, "_pb", pc))
        if (file.exists(file.path(outdir, "report.tsv"))) next

        cmd <- sprintf("%s -r %s -o %s --threads %s %s",
                       shQuote(quast), shQuote(ref_fa), shQuote(outdir),
                       threads, shQuote(asm))
        if (system(cmd) != 0L) stop("QUAST failed for ", basename(outdir))
      }
    }
  }
  TRUE
}
