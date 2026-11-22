# Intent extraction and workflow compilation are deliberately separate from
# knowledge retrieval. Retrieval explains decisions; this registry controls
# which tested simulation modules may be assembled into runnable code.

.simitall_agent_extract_intent <- function(question) {
  text <- tolower(question)
  terms <- .simitall_agent_terms(question)
  relationship_groups <- character()
  if (grepl("unrelated|random mating", text)) relationship_groups <- c(relationship_groups, "unrelated")
  if (grepl("distantly[ -]?related|distant relation|second[ -]?cousin", text)) {
    relationship_groups <- c(relationship_groups, "distantly_related")
  }
  if (grepl("first[ -]?cousin", text)) relationship_groups <- c(relationship_groups, "first_cousin_descendant")
  analyses <- c(
    if (grepl("gwas|association", text)) "gwas",
    if (grepl("qtl", text)) "qtl",
    if (grepl("rna.?seq", text)) "rnaseq",
    if (grepl("chip.?seq", text)) "chipseq",
    if (grepl("genomic.?selection|gblup", text)) "genomic_selection",
    if (grepl("gxe|genotype.by.environment|irrigated", text)) "gxe"
  )
  figures <- c(
    if (grepl("manhattan", text)) "manhattan",
    if (grepl("qq", text)) "qq",
    if (grepl("pca", text)) "pca",
    if (grepl("volcano", text)) "volcano",
    if (grepl("probability", text)) "probability_comparison",
    if (grepl("homozygos", text)) "homozygosity",
    if (grepl("distribution", text)) "trait_distribution"
  )
  is_human_pedigree <- identical(.simitall_agent_species(question), "human") &&
    grepl("family|pedigree|unrelated|related|cousin|consanguin", text)
  list(
    organism = .simitall_agent_species(question),
    region = if ("irf6" %in% terms) "IRF6" else NA_character_,
    relationship_groups = unique(relationship_groups),
    analyses = unique(analyses),
    figures = unique(figures),
    trait_model = if (grepl("probability|liability", text)) "synthetic_liability" else NA_character_,
    workflow = if (is_human_pedigree) "human_pedigree" else NA_character_
  )
}

.simitall_agent_workflow_registry <- function() {
  data.frame(
    workflow = c("human_pedigree", "nam_quantitative_gwas", "nam_binary_gwas", "nam_gxe_multiomics"),
    required_modules = c(
      "marker_frequencies;pedigree_templates;synthetic_liability;pedigree_figures",
      "breeding;quantitative_trait;gwas;gwas_figures",
      "breeding;binary_trait;binary_gwas;gwas_figures",
      "breeding;gxe;rnaseq;eqtl;gwas;genomic_selection"
    ),
    stringsAsFactors = FALSE
  )
}

.simitall_agent_compile_workflow <- function(question, preflight) {
  intent <- .simitall_agent_extract_intent(question)
  if (identical(intent$workflow, "human_pedigree") && identical(intent$region, "IRF6")) {
    return(.simitall_agent_human_pedigree_recipe(question, preflight))
  }
  ""
}
