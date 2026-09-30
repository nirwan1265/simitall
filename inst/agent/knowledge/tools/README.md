---
id: tool_and_tutorial_map
document_type: tool_collection
status: active
domains: [software, tutorials, vignettes, workflows]
source_priority: package_docs_and_user_supplied_tutorials
last_reviewed: 2026-09-29
---

# Tool and Tutorial Map

Each tool card records what a dependency does, where SimItAll uses it, and its
boundary. This map connects those cards with supplied tutorials/vignettes and
reproducible SimItAll workflows. It is navigation material, not a substitute
for the original manual or paper.

| Tool or package | SimItAll role | Tool card | Tutorial or documentation basis | Related workflow |
| --- | --- | --- | --- | --- |
| TASSEL | Plant-oriented association, LD, kinship, PCA, and GBS workflows | `TASSEL.Rmd` | User-supplied TASSEL guide plus current official manual | Plant/crop GWAS and linkage analysis |
| GEMMA | External relatedness-aware LMM/mvLMM/BSLMM association | `GEMMA.Rmd` | Official manual and primary paper | Animal/human GWAS benchmark comparisons |
| GAPIT/FarmCPU | Comparative GWAS models including FarmCPU | `GAPIT_FarmCPU.Rmd` | User-supplied GAPIT and FarmCPU manuals | GWAS benchmark comparisons |
| maSigPro | Bulk RNA-seq time-course differential expression | `maSigPro.Rmd` | User-supplied maSigPro tutorial | Future time-course RNA-seq workflow |
| powsimR | Bulk/single-cell RNA-seq power analysis | `powsimR.Rmd` | User-supplied powsimR vignette | RNA-seq and single-cell study design |
| RNA-seq resources | Analysis-stage decision support | `rnaseq_analysis_resources.Rmd` | User-supplied review and navigation page | Bulk RNA-seq workflows |
| SimuPOP | Forward-time populations and mating schemes via `reticulate` | `SimuPOP.Rmd` | User-supplied Peng tutorial/paper plus package documentation | Breeding populations, NAM, MAGIC, RIL, pedigree designs |
| simplePHENOTYPES | Causal phenotype and trait architecture simulation | `phenotype_packages.Rmd` | User-supplied reference/documentation | GWAS and breeding phenotypes |
| simstudy | Family, maternal, and study-design effects | `phenotype_packages.Rmd` | User-supplied package-index and simulation tutorial | Breeding phenotypes |
| pedtricks | Pedigree utilities | `phenotype_packages.Rmd` | User-supplied package help | Pedigree-aware breeding |
| rrBLUP | Mixed-model GWAS and genomic prediction | `rrBLUP.Rmd` | User-supplied rrBLUP documentation | GWAS and genomic selection |
| BGLR | Bayesian genomic prediction | `BGLR.Rmd` | User-supplied BGLR tutorial/PDF | Genomic selection |
| Splatter/splatPop | Population-aware single-cell simulation | `Seurat_Splatter.Rmd` | User-supplied splatPop tutorial | Single-cell eQTL |
| Seurat | Single-cell QC, embedding, clustering, and interpretation | `Seurat_Splatter.Rmd` | User-supplied Seurat tutorials | Single-cell eQTL |
| ChIPsim | ChIP-seq signal/read simulation | `ChIPsim.Rmd` | User-supplied ChIPsim introduction | ChIP-seq from annotations |
| ART | Illumina read simulation | `ART.Rmd` | Official documentation and cited paper | Hybrid assembly benchmark |
| PBSIM3 | PacBio long-read simulation | `PBSIM3.Rmd` | Official documentation and cited paper | Hybrid assembly benchmark |
| Badread | Nanopore read simulation | `Badread.Rmd` | Official documentation and cited paper | Hybrid assembly benchmark |
| Unicycler | Bacterial hybrid assembly | `Unicycler.Rmd` | Official documentation and cited paper | Hybrid assembly benchmark |
| QUAST | Assembly evaluation | `QUAST.Rmd` | Official documentation and cited paper | Hybrid assembly benchmark |
| RFMix | Local ancestry inference boundary | `RFMix.Rmd` | Official documentation and cited paper | Ancestry-aware GWAS |
| Tractor | Local-ancestry-aware association boundary | `Tractor.Rmd` | Official documentation and cited paper | Ancestry-aware GWAS |
| GENESIS | Relatedness/population-aware association tools | `GENESIS.Rmd` | Official documentation | Ancestry-aware GWAS |
| GenomeAdmixR | Admixture simulation/visualization support | `GenomeAdmixR.Rmd` | Installed package documentation | Ancestry-aware GWAS |
| Figure and diagnostic package map | Optional R backends for GWAS, RNA-seq, ChIP-seq, single-cell, breeding, and genomic-selection figures | `figure_visualization_packages.Rmd` | Official package documentation | Result reporting and validation |

## How the Agent Uses This Map

1. Retrieve the tool card and the matching workflow before suggesting code.
2. State whether SimItAll runs the tool directly, wraps it through R, or only
   documents an external-analysis boundary.
3. Prefer the supplied tutorial/vignette for tool-specific parameters.
4. Never claim a tool is bundled or executable unless the local installation
   and workflow explicitly support it.

## Adding a New Tutorial

For each new tutorial, add the original file outside the package source if it
is copyrighted, then add or update one tool card with: tool version, source
link/citation, core commands/functions, inputs, outputs, failure modes,
SimItAll wrapper function, and a runnable small demo.
