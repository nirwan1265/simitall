#' List optional R packages for SimItAll result figures
#'
#' Return a curated catalogue of optional R and Bioconductor packages that can
#' extend SimItAll's built-in result figures. The catalogue is deliberately
#' advisory: a plot should be selected because its input data and scientific
#' question support it, not because a package happens to be installed.
#'
#' SimItAll's core plotting functions only require `ggplot2`. Packages in this
#' catalogue are optional analysis and visualization backends. Some figures
#' also require inputs that simulation alone does not create, such as aligned
#' BAM files, bigWig tracks, peak calls, or a reference annotation.
#'
#' @param workflow One or more workflow names: `"gwas"`, `"rnaseq"`,
#'   `"eqtl"`, `"chipseq"`, `"scrnaseq"`, `"breeding"`,
#'   `"genomic_selection"`, or `"cross_cutting"`. The default returns every
#'   catalogue entry.
#' @param available_only If `TRUE`, keep only packages installed in the current
#'   R library.
#'
#' @return A data frame with package, repository, workflow, role, expected
#'   inputs, and example outputs. `available` reports whether the package is
#'   installed locally.
#' @examples
#' simitall_figure_tool_catalog("gwas")
#' simitall_figure_tool_catalog(c("rnaseq", "chipseq"), available_only = TRUE)
#' @export
simitall_figure_tool_catalog <- function(
    workflow = c(
      "gwas", "rnaseq", "eqtl", "chipseq", "scrnaseq", "breeding",
      "genomic_selection", "cross_cutting"
    ),
    available_only = FALSE) {
  catalogue <- data.frame(
    workflow = c(
      rep("gwas", 10L), rep("rnaseq", 11L), rep("eqtl", 4L),
      rep("chipseq", 8L), rep("scrnaseq", 11L), rep("breeding", 8L),
      rep("genomic_selection", 5L), rep("cross_cutting", 7L)
    ),
    package = c(
      "qqman", "CMplot", "SNPRelate", "GENESIS", "GWASTools",
      "locuszoomr", "susieR", "coloc", "bigsnpr", "PRROC",
      "DESeq2", "edgeR", "limma", "EnhancedVolcano", "PCAtools",
      "pcaExplorer", "ComplexHeatmap", "clusterProfiler", "fgsea",
      "WGCNA", "DEXSeq",
      "Qvalue", "susieR", "coloc", "ComplexUpset",
      "ChIPQC", "DiffBind", "csaw", "ChIPseeker", "ChIPpeakAnno",
      "motifmatchr", "universalmotif", "ComplexHeatmap",
      "Seurat", "scater", "scDblFinder", "scran", "scuttle",
      "batchelor", "clustree", "slingshot", "tradeSeq", "CellChat",
      "infercnv",
      "qtl", "qtl2", "LDheatmap", "vcfR", "adegenet", "pegas",
      "hierfstat", "sommer",
      "rrBLUP", "BGLR", "sommer", "yardstick", "pROC",
      "ComplexUpset", "circlize", "ggalluvial", "ggraph", "igraph",
      "patchwork", "ggplot2"
    ),
    repository = c(
      "CRAN", "CRAN", "Bioconductor", "Bioconductor", "Bioconductor",
      "CRAN", "CRAN", "CRAN", "CRAN", "CRAN",
      "Bioconductor", "Bioconductor", "Bioconductor", "Bioconductor",
      "Bioconductor", "Bioconductor", "Bioconductor", "Bioconductor",
      "Bioconductor", "CRAN", "Bioconductor",
      "Bioconductor", "CRAN", "CRAN", "CRAN",
      "Bioconductor", "Bioconductor", "Bioconductor", "Bioconductor",
      "Bioconductor", "Bioconductor", "CRAN", "Bioconductor",
      "CRAN", "Bioconductor", "Bioconductor", "Bioconductor", "Bioconductor",
      "Bioconductor", "CRAN", "Bioconductor", "Bioconductor", "GitHub",
      "Bioconductor",
      "CRAN", "GitHub", "Bioconductor", "CRAN", "CRAN", "CRAN",
      "CRAN", "CRAN",
      "CRAN", "CRAN", "CRAN", "CRAN", "CRAN",
      "CRAN", "CRAN", "CRAN", "CRAN", "CRAN", "CRAN", "CRAN"
    ),
    role = c(
      "Manhattan and QQ plots", "Circular and multi-track association plots",
      "Genotype PCA, kinship, and QC", "Relatedness-aware association and PCs",
      "GWAS QC and QQ diagnostics", "Regional association and LD plots",
      "Fine-mapping posterior summaries", "GWAS-eQTL colocalization",
      "Large genotype matrices and PRS", "ROC and precision-recall curves",
      "Negative-binomial differential expression", "Negative-binomial differential expression",
      "Voom/linear-model differential expression", "Volcano plots with labels",
      "RNA-seq PCA diagnostics", "Interactive RNA-seq exploratory analysis",
      "Heatmaps and annotation tracks", "GO/pathway enrichment plots",
      "Fast gene-set enrichment", "Co-expression modules and trait heatmaps",
      "Exon usage and differential splicing",
      "FDR and q-value diagnostics", "Fine-mapped eQTL credible sets",
      "GWAS-eQTL colocalization", "eQTL-sharing UpSet plots",
      "ChIP-seq quality metrics", "Differential binding and replicate plots",
      "Window-based differential binding", "Peak annotation, profiles, and TSS plots",
      "Peak annotation and overlap", "Motif matching", "Motif import and plotting",
      "Peak or signal heatmaps",
      "Single-cell QC, embedding, and markers", "Single-cell QC diagnostics",
      "Doublet detection", "Normalization and feature selection", "Single-cell utilities",
      "Batch correction", "Clustering-resolution diagnostics", "Trajectory inference",
      "Pseudotime differential expression", "Cell-cell communication", "CNV heatmaps",
      "QTL scans and effects", "Modern multi-parent QTL mapping", "LD heatmaps",
      "VCF import and population genetics", "Population structure and PCA",
      "Diversity and population-genetic statistics", "F statistics and structure",
      "Mixed models, GxE, and BLUPs",
      "GBLUP/RR-BLUP prediction", "Bayesian genomic prediction",
      "Mixed models, GxE, and BLUPs", "Prediction-performance metrics",
      "ROC, calibration, and AUC", "Set overlap diagrams", "Circos plots",
      "Alluvial/Sankey plots", "Network graphs", "Network data structures",
      "Multi-panel figure composition", "Core plotting"
    ),
    expected_inputs = c(
      "GWAS marker table", "GWAS marker table", "Genotypes and sample metadata",
      "Genotypes, pedigree/kinship, phenotype", "Genotypes and association results",
      "Peak locus, LD, and annotation", "Summary statistics or LD matrix",
      "GWAS and eQTL summary statistics", "Genotypes and phenotype", "Truth labels and scores",
      "Counts and design matrix", "Counts and design matrix", "Counts and design matrix",
      "Differential-expression table", "Normalized counts", "Counts and sample metadata",
      "Matrix plus row/column annotations", "Gene list and identifier mapping",
      "Ranked gene statistics", "Normalized expression and traits", "Exon-level counts",
      "eQTL p-values", "eQTL summary statistics", "GWAS and eQTL summaries",
      "Sets of eGenes or SNP-gene pairs",
      "BAM and control/sample sheet", "BAM/peak sets and design", "BAM and design",
      "Peak BED and genome annotation", "Peak BED and annotation", "Motif database and sequences",
      "Motif matrices", "Signal matrix or peak matrix",
      "Single-cell count matrix and metadata", "Single-cell count matrix and metadata",
      "Single-cell count matrix", "Single-cell count matrix", "Single-cell count matrix",
      "Integrated expression matrix", "Cluster assignments", "Reduced dimensions/cluster labels",
      "Pseudotime and counts", "Cell metadata and expression", "Expression and genomic positions",
      "Cross genotypes and phenotype", "Multi-parent genotypes and phenotype", "Genotypes and marker map",
      "VCF and metadata", "Genotypes and labels", "Genotypes and population labels",
      "Genotypes and population labels", "Phenotype, pedigree, environment",
      "Genotypes and phenotype", "Genotypes and phenotype", "Phenotype and environment",
      "Observed/predicted values", "Binary truth and scores", "Named sets", "Genomic links",
      "Long-format flow table", "Node/edge tables", "Node/edge tables", "ggplot objects", "Tidy data"
    ),
    example_outputs = c(
      "Manhattan, QQ", "Miami/circular association", "PCA, kinship heatmap, MAF QC",
      "PC-AiR/kinship-aware diagnostics", "HWE/QQ diagnostics", "LocusZoom-style regional plot",
      "Credible-set plot", "Colocalization posterior", "PRS and large-scale QC", "ROC/PR",
      "MA/dispersion and DE table", "MA/dispersion and DE table", "MA plot and batch diagnostics",
      "Labelled volcano", "PCA/scree plot", "Interactive PCA/QC", "Annotated heatmap",
      "GO dotplot", "GSEA curve", "Module-trait heatmap", "PSI/exon plots",
      "p-value histogram/FDR", "Credible set", "Colocalization", "UpSet plot",
      "FRiP/fragment QC", "Differential-binding MA/volcano", "Differential windows",
      "TSS profile/annotation", "Peak annotation", "Motif enrichment", "Motif logo", "Signal heatmap",
      "UMAP/marker plots", "QC violins/PCA", "Doublet score", "Normalization diagnostics",
      "Library/QC plots", "Before/after integration", "Resolution sweep", "Pseudotime trajectory",
      "Pseudotime gene trends", "Communication heatmap", "inferCNV heatmap",
      "LOD/effect plots", "Founder-effect and QTL plots", "LD heatmap", "MAF/heterozygosity",
      "PCA/admixture", "Tajima's D/diversity", "Fst", "GxE reaction norms/BLUPs",
      "Predicted-vs-observed", "Marker-effect/shrinkage plots", "Variance/GxE plots",
      "Accuracy and calibration", "ROC/calibration", "UpSet", "Circos", "Sankey", "Network", "Network", "Multi-panel report", "General QC figures"
    ),
    stringsAsFactors = FALSE
  )
  workflow <- unique(match.arg(workflow, unique(catalogue$workflow), several.ok = TRUE))
  catalogue <- catalogue[catalogue$workflow %in% workflow, , drop = FALSE]
  catalogue$available <- vapply(
    catalogue$package, requireNamespace, logical(1), quietly = TRUE
  )
  if (isTRUE(available_only)) catalogue <- catalogue[catalogue$available, , drop = FALSE]
  rownames(catalogue) <- NULL
  catalogue
}
