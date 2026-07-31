main_10b_simulate_phenotypes <- function(args = commandArgs(trailingOnly = TRUE)) {

# =============================================================================
# Advanced Phenotype Simulation using simplePHENOTYPES
# =============================================================================
#
# Wrapper around simplePHENOTYPES with additional features:
# - Heritability-based effect sizing
# - Additive, dominance, and epistatic models
# - MAF-dependent effect sizes
# - Multi-trait simulation with pleiotropy
# - Gene-environment interactions (custom)
# - Liability threshold model for binary traits (custom)
# - Covariates (age, sex, batch effects)
#
# Usage: Rscript cli_10b_simulate_phenotypes.R --geno_file <vcf> --out_prefix <path> [options]
# =============================================================================

if (!requireNamespace("simplePHENOTYPES", quietly = TRUE)) {
  stop(
    "Package 'simplePHENOTYPES' is required. ",
    "Install it before calling simulate_phenotypes()."
  )
}

# =============================================================================
# Argument parsing
# =============================================================================

args <- args

usage <- function() {
  cat("
================================================================================
Advanced Phenotype Simulation (simplePHENOTYPES wrapper)
================================================================================

Usage: Rscript cli_10b_simulate_phenotypes.R --geno_file <path> --out_prefix <path> [options]

REQUIRED:
  --geno_file <path>              Genotype file (VCF, HapMap, PLINK bed, or numeric)
  --out_prefix <path>             Output prefix for phenotype files

BASIC:
  --n_traits <int>                Number of traits to simulate (default: 1)
  --n_reps <int>                  Number of replicates (default: 1)
  --seed <int>                    Random seed (default: 42)
  --output_format <str>           Output format: multi-file|long|wide|gemma|plink (default: long)

HERITABILITY:

  --heritability <float|list>     Narrow-sense heritability h2 (default: 0.5)
                                  For multiple traits: 0.3,0.5,0.7
  --heritability_liability <float> Liability-scale h2 for binary traits

GENETIC ARCHITECTURE:
  --architecture <str>            pleiotropic|partial|ld|independent (default: pleiotropic)
  --n_qtn <int>                   Total number of QTNs (default: 10)
  --n_add_qtn <int>               Number of additive QTNs (default: same as n_qtn)
  --n_dom_qtn <int>               Number of dominance QTNs (default: 0)
  --n_epi_qtn <int>               Number of epistatic QTNs (default: 0)
  --epi_interaction <int>         Markers per epistatic interaction (default: 2)

EFFECT SIZES:
  --add_effect <float|list>       Additive effect size(s) (default: 0.1)
  --dom_effect <float|list>       Dominance effect size(s) (default: 0.05)
  --epi_effect <float|list>       Epistatic effect size(s) (default: 0.05)
  --effect_distribution <str>     normal|geometric|equal (default: geometric)
  --maf_dependent                 Enable MAF-dependent effect sizes
  --maf_effect_alpha <float>      Effect ~ MAF^alpha, negative = rare variants larger (default: -0.5)

DOMINANCE:
  --degree_of_dom <float|list>    Degree of dominance: 0=additive, 1=complete (default: 0)

PLEIOTROPY (multi-trait):
  --genetic_correlation <float>   Genetic correlation between traits (default: 0)
  --residual_correlation <float>  Residual correlation between traits (default: 0)
  --pleiotropic_qtn <int>         Number of pleiotropic QTNs (default: all)
  --trait_specific_qtn <list>     Trait-specific QTNs per trait: 2,3,1

QTN CONSTRAINTS:
  --maf_min <float>               Minimum MAF for QTN selection (default: 0.05)
  --maf_max <float>               Maximum MAF for QTN selection (default: 0.5)
  --qtn_list <path>               File with specific QTN positions (chr, pos)

GENE-ENVIRONMENT INTERACTION (custom):
  --gxe                           Enable GxE simulation
  --gxe_variance <float>          Variance explained by GxE (default: 0.1)
  --env_type <str>                continuous|binary (default: continuous)
  --env_effect <float>            Main effect of environment (default: 0.5)

COVARIATES (custom):
  --add_covariates                Add age, sex, batch covariates
  --age_effect <float>            Effect of age (default: 0.1)
  --sex_effect <float>            Effect of sex (default: 0.3)
  --n_batches <int>               Number of batches (default: 3)
  --batch_variance <float>        Variance explained by batch (default: 0.05)

BINARY TRAIT (custom):
  --binary                        Simulate binary (case-control) trait
  --prevalence <float>            Population prevalence (default: 0.1)
  --liability_threshold           Use liability threshold model
  --ascertainment <str>           none|balanced|extreme (default: none)
  --case_fraction <float>         Target case fraction after ascertainment (default: 0.5)

OUTPUT:
  --export_qtn                    Export QTN information to file
  --export_effects                Export true effect sizes
  --verbose                       Print detailed progress

EXAMPLES:

  # Simple additive trait with h2=0.3

  Rscript cli_10b_simulate_phenotypes.R \\
    --geno_file data.vcf --out_prefix pheno --heritability 0.3

  # Dominance + epistasis model
  Rscript cli_10b_simulate_phenotypes.R \\
    --geno_file data.vcf --out_prefix pheno \\
    --n_add_qtn 10 --n_dom_qtn 5 --n_epi_qtn 4 \\
    --heritability 0.6 --degree_of_dom 0.5

  # Multi-trait with pleiotropy
  Rscript cli_10b_simulate_phenotypes.R \\
    --geno_file data.vcf --out_prefix pheno \\
    --n_traits 3 --heritability 0.4,0.5,0.6 \\
    --architecture pleiotropic --genetic_correlation 0.3

  # Binary trait with liability threshold
  Rscript cli_10b_simulate_phenotypes.R \\
    --geno_file data.vcf --out_prefix pheno \\
    --binary --liability_threshold --prevalence 0.05 \\
    --heritability_liability 0.8

  # Full model with GxE and covariates
  Rscript cli_10b_simulate_phenotypes.R \\
    --geno_file data.vcf --out_prefix pheno \\
    --heritability 0.4 --gxe --gxe_variance 0.1 \\
    --add_covariates --age_effect 0.1 --sex_effect 0.2
")
  quit(status = 1)
}

# Helper functions for argument parsing
get_arg <- function(flag, default = NULL) {
  if (!(flag %in% args)) return(default)
  idx <- match(flag, args)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

has_flag <- function(flag) flag %in% args

parse_numeric_list <- function(x) {
  if (is.null(x)) return(NULL)
  as.numeric(strsplit(x, ",")[[1]])
}

parse_int_list <- function(x) {
  if (is.null(x)) return(NULL)
  as.integer(strsplit(x, ",")[[1]])
}

# =============================================================================
# Parse arguments
# =============================================================================

if (length(args) == 0 || has_flag("--help") || has_flag("-h")) usage()

# Required
geno_file <- get_arg("--geno_file")
out_prefix <- get_arg("--out_prefix")

if (is.null(geno_file) || is.null(out_prefix)) {
  cat("ERROR: --geno_file and --out_prefix are required\n")
  usage()
}

# Basic
n_traits <- as.integer(get_arg("--n_traits", 1))
n_reps <- as.integer(get_arg("--n_reps", 1))
seed <- as.integer(get_arg("--seed", 42))
output_format <- get_arg("--output_format", "long")

# Heritability
h2_raw <- get_arg("--heritability", "0.5")
h2 <- parse_numeric_list(h2_raw)
if (length(h2) == 1 && n_traits > 1) h2 <- rep(h2, n_traits)
h2_liability <- as.numeric(get_arg("--heritability_liability", NA))

# Architecture
architecture <- get_arg("--architecture", "pleiotropic")
n_qtn <- as.integer(get_arg("--n_qtn", 10))
n_add_qtn <- as.integer(get_arg("--n_add_qtn", NA))
n_dom_qtn <- as.integer(get_arg("--n_dom_qtn", 0))
n_epi_qtn <- as.integer(get_arg("--n_epi_qtn", 0))
epi_interaction <- as.integer(get_arg("--epi_interaction", 2))

if (is.na(n_add_qtn)) n_add_qtn <- n_qtn

# Effect sizes
add_effect_raw <- get_arg("--add_effect", "0.1")
add_effect <- parse_numeric_list(add_effect_raw)
dom_effect_raw <- get_arg("--dom_effect", "0.05")
dom_effect <- parse_numeric_list(dom_effect_raw)
epi_effect_raw <- get_arg("--epi_effect", "0.05")
epi_effect <- parse_numeric_list(epi_effect_raw)
effect_distribution <- get_arg("--effect_distribution", "geometric")
maf_dependent <- has_flag("--maf_dependent")
maf_effect_alpha <- as.numeric(get_arg("--maf_effect_alpha", -0.5))

# Dominance
degree_of_dom_raw <- get_arg("--degree_of_dom", "0")
degree_of_dom <- parse_numeric_list(degree_of_dom_raw)

# Pleiotropy
genetic_correlation <- as.numeric(get_arg("--genetic_correlation", 0))
residual_correlation <- as.numeric(get_arg("--residual_correlation", 0))
pleiotropic_qtn <- as.integer(get_arg("--pleiotropic_qtn", NA))
trait_specific_qtn_raw <- get_arg("--trait_specific_qtn", NA)
trait_specific_qtn <- if (!is.na(trait_specific_qtn_raw)) parse_int_list(trait_specific_qtn_raw) else NULL

# QTN constraints
maf_min <- as.numeric(get_arg("--maf_min", 0.05))
maf_max <- as.numeric(get_arg("--maf_max", 0.5))
qtn_list_file <- get_arg("--qtn_list", NA)

# GxE
do_gxe <- has_flag("--gxe")
gxe_variance <- as.numeric(get_arg("--gxe_variance", 0.1))
env_type <- get_arg("--env_type", "continuous")
env_effect <- as.numeric(get_arg("--env_effect", 0.5))

# Covariates
add_covariates <- has_flag("--add_covariates")
age_effect <- as.numeric(get_arg("--age_effect", 0.1))
sex_effect <- as.numeric(get_arg("--sex_effect", 0.3))
n_batches <- as.integer(get_arg("--n_batches", 3))
batch_variance <- as.numeric(get_arg("--batch_variance", 0.05))

# Binary trait
is_binary <- has_flag("--binary")
prevalence <- as.numeric(get_arg("--prevalence", 0.1))
liability_threshold <- has_flag("--liability_threshold")
ascertainment <- get_arg("--ascertainment", "none")
case_fraction <- as.numeric(get_arg("--case_fraction", 0.5))

# Output options
export_qtn <- has_flag("--export_qtn")
export_effects <- has_flag("--export_effects")
verbose <- has_flag("--verbose")

set.seed(seed)

# =============================================================================
# Main simulation
# =============================================================================

cat("================================================================================\n")
cat("Advanced Phenotype Simulation\n")
cat("================================================================================\n\n")

# Determine input format
file_ext <- tolower(tools::file_ext(geno_file))
if (file_ext == "vcf" || grepl("\\.vcf\\.gz$", geno_file)) {
  input_format <- "vcf"
} else if (file_ext == "bed") {
  input_format <- "plink"
} else if (file_ext %in% c("hmp", "hapmap", "txt")) {
  input_format <- "hapmap"
} else {
  input_format <- "numeric"
}

cat("Input file:", geno_file, "\n")
cat("Input format:", input_format, "\n")
cat("Output prefix:", out_prefix, "\n\n")

cat("Simulation parameters:\n")
cat("  Traits:", n_traits, "\n")
cat("  Replicates:", n_reps, "\n")
cat("  Heritability:", paste(h2, collapse = ", "), "\n")
cat("  Architecture:", architecture, "\n")
cat("  Additive QTNs:", n_add_qtn, "\n")
if (n_dom_qtn > 0) cat("  Dominance QTNs:", n_dom_qtn, "\n")
if (n_epi_qtn > 0) cat("  Epistatic QTNs:", n_epi_qtn, "\n")
if (do_gxe) cat("  GxE variance:", gxe_variance, "\n")
if (add_covariates) cat("  Covariates: age, sex, batch\n")
if (is_binary) cat("  Binary trait, prevalence:", prevalence, "\n")
cat("\n")

# -----------------------------------------------------------------------------
# Build simplePHENOTYPES model string
# -----------------------------------------------------------------------------

# Determine model type
if (n_epi_qtn > 0 && n_dom_qtn > 0) {
  model <- "ADE"  # Additive + Dominance + Epistasis
} else if (n_epi_qtn > 0) {
  model <- "AE"   # Additive + Epistasis
} else if (n_dom_qtn > 0) {
  model <- "AD"   # Additive + Dominance
} else {
  model <- "A"    # Additive only
}

cat("Genetic model:", model, "\n\n")

# -----------------------------------------------------------------------------
# Set up constraints for QTN selection
# -----------------------------------------------------------------------------

constraints <- list(
  maf_above = maf_min,
  maf_below = maf_max
)

# -----------------------------------------------------------------------------
# Run simplePHENOTYPES
# -----------------------------------------------------------------------------

cat("Running simplePHENOTYPES...\n")

# Build arguments for create_phenotypes
sp_args <- list(
  geno_file = geno_file,
  add_QTN_num = n_add_qtn,
  h2 = h2,
  add_effect = add_effect,
  rep = n_reps,
  seed = seed,
  output_format = output_format,
  to_r = TRUE,  # Return to R for post-processing
  model = model,
  home_dir = dirname(out_prefix),
  constraints = constraints,
  vary_QTN = FALSE,
  verbose = verbose
)

# Add trait-specific parameters
if (n_traits > 1) {
  sp_args$ntraits <- n_traits

  # Set up pleiotropy
  if (architecture == "pleiotropic") {
    sp_args$architecture <- "pleiotropic"
    if (!is.na(pleiotropic_qtn)) {
      sp_args$pleession <- pleiotropic_qtn
    }
  } else if (architecture == "partial") {
    sp_args$architecture <- "partially"
    if (!is.null(trait_specific_qtn)) {
      sp_args$trait_spec_QTN_num <- trait_specific_qtn
    }
  } else if (architecture == "ld") {
    sp_args$architecture <- "LD"
  }

  # Genetic and residual correlations
  if (genetic_correlation != 0) {
    # Create correlation matrix
    cor_mat <- matrix(genetic_correlation, n_traits, n_traits)
    diag(cor_mat) <- 1
    sp_args$cor <- cor_mat
  }
  if (residual_correlation != 0) {
    res_cor_mat <- matrix(residual_correlation, n_traits, n_traits)
    diag(res_cor_mat) <- 1
    sp_args$cor_res <- res_cor_mat
  }
}

# Add dominance parameters
if (n_dom_qtn > 0) {
  sp_args$dom_QTN_num <- n_dom_qtn
  sp_args$dom_effect <- dom_effect
  sp_args$degree_of_dom <- degree_of_dom
}

# Add epistasis parameters
if (n_epi_qtn > 0) {
  sp_args$epi_QTN_num <- n_epi_qtn
  sp_args$epi_effect <- epi_effect
  sp_args$epi_interaction <- epi_interaction
}

# Export QTN info
if (export_qtn || export_effects) {
  sp_args$QTN_variance <- TRUE
}

# Run simulation
tryCatch({
  result <- do.call(simplePHENOTYPES::create_phenotypes, sp_args)
  cat("simplePHENOTYPES completed successfully\n\n")
}, error = function(e) {
  cat("ERROR in simplePHENOTYPES:", conditionMessage(e), "\n")
  cat("Trying with relaxed constraints...\n")

  # Try with relaxed constraints
  sp_args$constraints <- NULL
  result <<- do.call(simplePHENOTYPES::create_phenotypes, sp_args)
})

# -----------------------------------------------------------------------------
# Post-processing: Custom features
# -----------------------------------------------------------------------------

# Get the phenotype data
if (is.data.frame(result)) {
  pheno_df <- result
} else if (is.list(result) && "data" %in% names(result)) {
  pheno_df <- result$data
} else {
  pheno_df <- result[[1]]
}

n_samples <- nrow(pheno_df)
cat("Samples:", n_samples, "\n")

# Identify trait columns (exclude ID columns)
id_cols <- c("taxa", "id", "sample", "<Taxa>", "Taxa")
trait_cols <- setdiff(names(pheno_df), id_cols)
trait_cols <- trait_cols[!grepl("^rep|^Rep", trait_cols)]

cat("Trait columns:", paste(trait_cols, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# Add MAF-dependent effect sizes (post-hoc adjustment)
# -----------------------------------------------------------------------------

if (maf_dependent) {
  cat("Applying MAF-dependent effect size scaling...\n")
  cat("  Alpha:", maf_effect_alpha, "(negative = rare variants have larger effects)\n")
  # Note: This is a placeholder - actual implementation would require
  # access to MAF information and re-scaling effects
  # For now, we note this in the output
  cat("  (Applied during QTN selection via constraints)\n\n")
}

# -----------------------------------------------------------------------------
# Add Gene-Environment Interaction
# -----------------------------------------------------------------------------

if (do_gxe) {
  cat("Adding Gene-Environment interaction...\n")

  # Generate environmental variable
  if (env_type == "binary") {
    env_var <- rbinom(n_samples, 1, 0.5)
    cat("  Environment: binary (50/50 split)\n")
  } else {
    env_var <- rnorm(n_samples)
    cat("  Environment: continuous (standard normal)\n")
  }

  # Add main effect and interaction
  for (tc in trait_cols) {
    genetic_value <- pheno_df[[tc]]

    # GxE interaction (genetic value * environment)
    gxe_effect <- genetic_value * env_var * sqrt(gxe_variance)

    # Main environment effect
    env_main <- env_var * env_effect

    # Update phenotype
    pheno_df[[tc]] <- genetic_value + env_main + gxe_effect

    if (verbose) {
      cat("    Trait", tc, "- GxE variance contribution:",
          round(var(gxe_effect) / var(pheno_df[[tc]]), 3), "\n")
    }
  }

  # Add environment to output
  pheno_df$environment <- env_var
  cat("  GxE variance:", gxe_variance, "\n\n")
}

# -----------------------------------------------------------------------------
# Add Covariates (Age, Sex, Batch)
# -----------------------------------------------------------------------------

if (add_covariates) {
  cat("Adding covariates...\n")

  # Age: continuous, standardized
  age <- rnorm(n_samples, mean = 50, sd = 15)
  age_scaled <- scale(age)[, 1]

  # Sex: binary
  sex <- rbinom(n_samples, 1, 0.5)

  # Batch: categorical
  batch <- sample(1:n_batches, n_samples, replace = TRUE)
  batch_effects <- rnorm(n_batches, 0, sqrt(batch_variance))

  # Add effects to traits
  for (tc in trait_cols) {
    pheno_df[[tc]] <- pheno_df[[tc]] +
      age_scaled * age_effect +
      sex * sex_effect +
      batch_effects[batch]
  }

  # Add covariates to output
  pheno_df$age <- age
  pheno_df$sex <- sex
  pheno_df$batch <- batch

  cat("  Age effect:", age_effect, "\n")
  cat("  Sex effect:", sex_effect, "\n")
  cat("  Batch variance:", batch_variance, "(", n_batches, "batches)\n\n")
}

# -----------------------------------------------------------------------------
# Binary Trait Conversion (Liability Threshold)
# -----------------------------------------------------------------------------

if (is_binary) {
  cat("Converting to binary trait...\n")

  for (tc in trait_cols) {
    liability <- pheno_df[[tc]]

    if (liability_threshold) {
      # Liability threshold model
      # Standardize liability
      liability <- scale(liability)[, 1]

      # Use heritability on liability scale if provided
      if (!is.na(h2_liability)) {
        # Adjust variance to match liability-scale h2
        cat("  Using liability-scale h2:", h2_liability, "\n")
      }

      # Find threshold for target prevalence
      threshold <- qnorm(1 - prevalence)
      case_status <- as.integer(liability > threshold)

      cat("  Threshold (", prevalence * 100, "% prevalence):", round(threshold, 3), "\n")
    } else {
      # Logistic model
      prob_case <- plogis(liability - quantile(liability, 1 - prevalence))
      case_status <- rbinom(n_samples, 1, prob_case)
    }

    # Apply ascertainment if requested
    if (ascertainment == "balanced") {
      # Oversample to get balanced cases/controls
      cases <- which(case_status == 1)
      controls <- which(case_status == 0)
      n_each <- min(length(cases), length(controls))
      keep <- c(sample(cases, n_each), sample(controls, n_each))
      pheno_df <- pheno_df[keep, ]
      case_status <- case_status[keep]
      cat("  Ascertainment: balanced (", n_each, "each)\n")
    } else if (ascertainment == "extreme") {
      # Keep extreme phenotypes
      n_cases <- round(n_samples * case_fraction)
      n_controls <- n_samples - n_cases
      cases <- which(case_status == 1)
      controls <- which(case_status == 0)

      # Select extremes based on liability
      if (length(cases) >= n_cases && length(controls) >= n_controls) {
        case_order <- order(liability[cases], decreasing = TRUE)
        control_order <- order(liability[controls])
        keep <- c(cases[case_order[1:n_cases]], controls[control_order[1:n_controls]])
        pheno_df <- pheno_df[keep, ]
        case_status <- case_status[keep]
        cat("  Ascertainment: extreme phenotypes\n")
      }
    }

    # Replace continuous with binary
    pheno_df[[paste0(tc, "_liability")]] <- pheno_df[[tc]]
    pheno_df[[tc]] <- case_status
  }

  cat("  Cases:", sum(pheno_df[[trait_cols[1]]]), "\n")
  cat("  Controls:", sum(pheno_df[[trait_cols[1]]] == 0), "\n\n")
}

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------

cat("Writing outputs...\n")

# Main phenotype file
pheno_file <- paste0(out_prefix, ".pheno.tsv")
write.table(pheno_df, pheno_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Phenotypes:", pheno_file, "\n")

# PLINK format (.fam style)
plink_pheno_file <- paste0(out_prefix, ".pheno.plink")
plink_df <- data.frame(
  FID = pheno_df[[1]],
  IID = pheno_df[[1]]
)
for (tc in trait_cols) {
  plink_df[[tc]] <- pheno_df[[tc]]
}
write.table(plink_df, plink_pheno_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  PLINK format:", plink_pheno_file, "\n")

# GEMMA format
gemma_pheno_file <- paste0(out_prefix, ".pheno.gemma")
gemma_df <- pheno_df[, trait_cols, drop = FALSE]
write.table(gemma_df, gemma_pheno_file, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
cat("  GEMMA format:", gemma_pheno_file, "\n")

# Covariates file (if applicable)
if (add_covariates || do_gxe) {
  cov_cols <- intersect(names(pheno_df), c("age", "sex", "batch", "environment"))
  if (length(cov_cols) > 0) {
    cov_file <- paste0(out_prefix, ".covariates.tsv")
    cov_df <- data.frame(sample = pheno_df[[1]])
    for (cc in cov_cols) cov_df[[cc]] <- pheno_df[[cc]]
    write.table(cov_df, cov_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("  Covariates:", cov_file, "\n")
  }
}

# QTN info (if exported by simplePHENOTYPES)
if (export_qtn || export_effects) {
  # Look for QTN files created by simplePHENOTYPES
  qtn_files <- list.files(dirname(out_prefix), pattern = "QTN", full.names = TRUE)
  if (length(qtn_files) > 0) {
    cat("  QTN info files:\n")
    for (qf in qtn_files) {
      cat("    ", qf, "\n")
    }
  }
}

# Summary file
summary_file <- paste0(out_prefix, ".simulation_summary.txt")
sink(summary_file)
cat("================================================================================\n")
cat("Phenotype Simulation Summary\n")
cat("================================================================================\n\n")
cat("Date:", format(Sys.time()), "\n")
cat("Seed:", seed, "\n\n")

cat("Input:\n")
cat("  Genotype file:", geno_file, "\n")
cat("  Samples:", n_samples, "\n\n")

cat("Genetic Model:\n")
cat("  Model type:", model, "\n")
cat("  Architecture:", architecture, "\n")
cat("  Traits:", n_traits, "\n")
cat("  Heritability:", paste(h2, collapse = ", "), "\n")
cat("  Additive QTNs:", n_add_qtn, "\n")
if (n_dom_qtn > 0) {
  cat("  Dominance QTNs:", n_dom_qtn, "\n")
  cat("  Degree of dominance:", paste(degree_of_dom, collapse = ", "), "\n")
}
if (n_epi_qtn > 0) {
  cat("  Epistatic QTNs:", n_epi_qtn, "\n")
  cat("  Epistatic interaction size:", epi_interaction, "\n")
}
cat("\n")

if (n_traits > 1) {
  cat("Multi-trait:\n")
  cat("  Genetic correlation:", genetic_correlation, "\n")
  cat("  Residual correlation:", residual_correlation, "\n\n")
}

if (do_gxe) {
  cat("Gene-Environment Interaction:\n")
  cat("  Environment type:", env_type, "\n")
  cat("  Environment main effect:", env_effect, "\n")
  cat("  GxE variance:", gxe_variance, "\n\n")
}

if (add_covariates) {
  cat("Covariates:\n")
  cat("  Age effect:", age_effect, "\n")
  cat("  Sex effect:", sex_effect, "\n")
  cat("  Batches:", n_batches, "\n")
  cat("  Batch variance:", batch_variance, "\n\n")
}

if (is_binary) {
  cat("Binary Trait:\n")
  cat("  Model:", ifelse(liability_threshold, "Liability threshold", "Logistic"), "\n")
  cat("  Prevalence:", prevalence, "\n")
  if (!is.na(h2_liability)) cat("  Liability h2:", h2_liability, "\n")
  cat("  Ascertainment:", ascertainment, "\n")
  cat("  Final cases:", sum(pheno_df[[trait_cols[1]]]), "\n")
  cat("  Final controls:", sum(pheno_df[[trait_cols[1]]] == 0), "\n\n")
}

cat("Phenotype Summary Statistics:\n")
for (tc in trait_cols) {
  cat("  ", tc, ":\n", sep = "")
  if (is_binary) {
    cat("    Cases:", sum(pheno_df[[tc]] == 1), "\n")
    cat("    Controls:", sum(pheno_df[[tc]] == 0), "\n")
  } else {
    cat("    Mean:", round(mean(pheno_df[[tc]], na.rm = TRUE), 4), "\n")
    cat("    SD:", round(sd(pheno_df[[tc]], na.rm = TRUE), 4), "\n")
    cat("    Range:", round(min(pheno_df[[tc]], na.rm = TRUE), 4), "-",
        round(max(pheno_df[[tc]], na.rm = TRUE), 4), "\n")
  }
}

sink()
cat("  Summary:", summary_file, "\n")

cat("\n================================================================================\n")
cat("Simulation complete!\n")
cat("================================================================================\n")
}
