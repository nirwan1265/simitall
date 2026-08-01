.simitall_dependency_catalog <- function(profile) {
  core <- data.frame(
    component = c("jsonlite", "Matrix", "reticulate"),
    kind = "R package",
    candidates = c("jsonlite", "Matrix", "reticulate"),
    required_for = "minimal",
    install_hint = "Installed by the minimal profile",
    stringsAsFactors = FALSE
  )
  population <- data.frame(
    component = c("simplePHENOTYPES", "SimuPOP"),
    kind = c("R package", "Python module"),
    candidates = c("simplePHENOTYPES", "simuPOP"),
    required_for = "population",
    install_hint = c(
      "remotes::install_github('samuelbfernandes/simplePHENOTYPES')",
      "conda install -c conda-forge simupop"
    ),
    stringsAsFactors = FALSE
  )
  omics <- data.frame(
    component = c(
      "Rsubread", "ChIPsim", "SingleCellExperiment",
      "SummarizedExperiment", "splatter", "ggplot2", "patchwork"
    ),
    kind = c(rep("R package", 7L)),
    candidates = c(
      "Rsubread", "ChIPsim", "SingleCellExperiment",
      "SummarizedExperiment", "splatter", "ggplot2", "patchwork"
    ),
    required_for = "omics",
    install_hint = c(
      rep("Install with BiocManager::install()", 5L),
      rep("Install from CRAN", 2L)
    ),
    stringsAsFactors = FALSE
  )
  sequencing <- data.frame(
    component = c(
      "ART", "PBSIM/PBSIM3", "Badread", "Unicycler", "QUAST",
      "samtools", "seqkit", "pigz/gzip"
    ),
    kind = "command",
    candidates = c(
      "art_illumina", "pbsim|pbsim3", "badread", "unicycler",
      "quast.py|quast", "samtools", "seqkit", "pigz|gzip"
    ),
    required_for = "sequencing",
    install_hint = c(
      rep("Install through Conda/Bioconda", 4L),
      "Install QUAST through pip or Bioconda",
      rep("Install through Conda/Bioconda", 3L)
    ),
    stringsAsFactors = FALSE
  )
  selected <- switch(
    profile,
    minimal = list(core),
    population = list(core, population),
    omics = list(core, omics),
    sequencing = list(core, sequencing),
    full = list(core, population, omics, sequencing)
  )
  unique(do.call(rbind, selected))
}

.simitall_check_python_module <- function(module) {
  python <- Sys.which(c("python", "python3"))
  python <- unname(python[nzchar(python)][1L])
  if (!length(python) || is.na(python)) {
    return(list(available = FALSE, location = NA_character_))
  }
  status <- suppressWarnings(system2(
    python,
    c("-c", shQuote(paste0("import ", module))),
    stdout = FALSE,
    stderr = FALSE
  ))
  list(available = identical(status, 0L), location = python)
}

#' Check simitall software dependencies
#'
#' Inspect the R packages, Python modules, and external commands required by a
#' selected installation profile. This function does not install or modify
#' anything and is safe to use in scripts, support requests, and bug reports.
#'
#' @param profile Dependency profile: `"minimal"`, `"population"`,
#'   `"omics"`, `"sequencing"`, or `"full"`.
#' @param verbose Print a compact availability table and installation summary.
#'
#' @return Invisibly returns a data frame with component, dependency type,
#'   profile, availability, detected location/version, and installation hint.
#' @examples
#' check_simitall_dependencies("minimal")
#' \dontrun{
#' check_simitall_dependencies("full")
#' }
#' @export
check_simitall_dependencies <- function(
    profile = c("minimal", "population", "omics", "sequencing", "full"),
    verbose = TRUE) {
  profile <- match.arg(profile)
  result <- .simitall_dependency_catalog(profile)
  result$available <- FALSE
  result$detected <- NA_character_
  for (i in seq_len(nrow(result))) {
    kind <- result$kind[i]
    candidates <- strsplit(result$candidates[i], "|", fixed = TRUE)[[1L]]
    if (kind == "R package") {
      package <- candidates[1L]
      available <- requireNamespace(package, quietly = TRUE)
      result$available[i] <- available
      if (available) {
        result$detected[i] <- as.character(utils::packageVersion(package))
      }
    } else if (kind == "Python module") {
      check <- .simitall_check_python_module(candidates[1L])
      result$available[i] <- check$available
      result$detected[i] <- check$location
    } else {
      locations <- Sys.which(candidates)
      locations <- unname(locations[nzchar(locations)])
      result$available[i] <- length(locations) > 0L
      if (length(locations)) result$detected[i] <- locations[1L]
    }
  }
  result$candidates <- NULL
  result <- result[, c(
    "component", "kind", "required_for", "available", "detected",
    "install_hint"
  )]
  if (isTRUE(verbose)) {
    display <- data.frame(
      status = ifelse(result$available, "[OK]", "[MISSING]"),
      component = result$component,
      type = result$kind,
      detected = ifelse(is.na(result$detected), "", result$detected),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    print(display, row.names = FALSE)
    missing <- sum(!result$available)
    message(
      "Profile '", profile, "': ", nrow(result) - missing, "/",
      nrow(result), " dependencies available."
    )
    if (missing) {
      message(
        "Run install_simitall_dependencies('", profile,
        "') or use the platform bootstrap installer."
      )
    }
  }
  invisible(result)
}

#' Install a simitall dependency profile
#'
#' Run the platform bootstrap installer shipped with the package. The installer
#' creates or updates a Conda environment and installs only the dependencies
#' selected by the profile. On native Windows, the `sequencing` and `full`
#' profiles delegate to WSL2 because ART, PBSIM, and Unicycler are distributed
#' through Linux/macOS Bioconda channels.
#'
#' @param profile Installation profile: `"minimal"`, `"population"`,
#'   `"omics"`, `"sequencing"`, or `"full"`.
#' @param envname Conda environment name.
#' @param ask Ask for confirmation before running in an interactive session.
#' @param dry_run Return the command without executing it.
#'
#' @return Invisibly returns the installer exit status. With `dry_run = TRUE`,
#'   returns a list containing the command and arguments.
#' @examples
#' \dontrun{
#' install_simitall_dependencies("population")
#' install_simitall_dependencies("full", envname = "simitall")
#' }
#' @export
install_simitall_dependencies <- function(
    profile = c("minimal", "population", "omics", "sequencing", "full"),
    envname = "simitall",
    ask = interactive(),
    dry_run = FALSE) {
  profile <- match.arg(profile)
  if (length(envname) != 1L || is.na(envname) ||
      !grepl("^[A-Za-z0-9_.-]+$", envname)) {
    stop("envname may contain only letters, numbers, dots, underscores, and hyphens")
  }
  windows <- identical(.Platform$OS.type, "windows")
  installer_name <- if (windows) {
    "install_simitall_windows.ps1"
  } else {
    "install_simitall_unix.sh"
  }
  installer <- system.file("installers", installer_name, package = "simitall")
  if (!nzchar(installer) || !file.exists(installer)) {
    source_candidate <- file.path("inst", "installers", installer_name)
    if (file.exists(source_candidate)) installer <- normalizePath(source_candidate)
  }
  if (!nzchar(installer) || !file.exists(installer)) {
    stop("The platform installer could not be located")
  }
  if (windows) {
    command <- Sys.which("powershell.exe")
    if (!nzchar(command)) stop("powershell.exe was not found")
    arguments <- c(
      "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
      shQuote(installer), "-Profile", profile, "-EnvName", envname
    )
  } else {
    command <- Sys.which("bash")
    if (!nzchar(command)) stop("bash was not found")
    arguments <- c(shQuote(installer), profile, envname)
  }
  specification <- list(
    command = unname(command), arguments = arguments,
    profile = profile, envname = envname
  )
  if (isTRUE(dry_run)) return(invisible(specification))
  if (isTRUE(ask)) {
    answer <- readline(paste0(
      "Install the '", profile, "' profile into Conda environment '",
      envname, "'? [y/N] "
    ))
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      message("Installation cancelled; no changes were made.")
      return(invisible(1L))
    }
  }
  status <- system2(command, arguments)
  if (!identical(status, 0L)) {
    stop("Dependency installer exited with status ", status)
  }
  invisible(status)
}
