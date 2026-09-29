---
id: paper_map_by_analysis
document_type: paper_collection
status: active
domains: [GWAS, RNAseq, ChIPseq, single_cell, breeding, sequencing, population_genetics, ancestry]
source_priority: navigation_only
last_reviewed: 2026-09-29
---

# Paper Map by Analysis Area

This index maps the individual paper summaries in this directory to their main
analysis role. A paper can appear in more than one area. The map is for
navigation and retrieval; use the individual `.Rmd` summary and original source
for methods, numbers, and claims.

## GWAS, QTL, Population Structure, and Ancestry

- **GWAS/QTL in plants and maize:** `altaf_2024_qtl_gwas_plant_improvement.Rmd`,
  `sahito_2024_gwas_maize.Rmd`, `dinanty_2026_qtl_analysis_sorghum.Rmd`,
  `sharma_2023_meta_qtl_wheat.Rmd`, `singh_2026_qtl_mapping_plants_genomic_era.Rmd`.
- **Association models and structure:** `yu_2006_a_unified_mixed_model_method_for_association.Rmd`,
  `kang_2008_efficient_control_of_population_structure_in.Rmd`,
  `zhou_2012_genome_wide_efficient_mixed_model_analysis_fo.Rmd`,
  `lipka_2012_gapit_genome_association_and_prediction_integ.Rmd`,
  `bradbury_2007_tassel_software_for_association_mapping_of_complex.Rmd`,
  `zhou_2012_genome_wide_efficient_mixed_model_analysis_fo.Rmd`,
  `price_2006_principal_components_analysis_corrects_for_st.Rmd`,
  `pritchard_2000_association_mapping_in_structured_populations.Rmd`.
- **Human/admixture-aware association:** `atkinson_2021_tractor_uses_local_ancestry_to_enable_the_inclusion.Rmd`,
  `maples_2013_rfmix_a_discriminative_modeling_approach_for_rapid_an.Rmd`,
  `alexander_2009_fast_model_based_estimation_of_ancestry_in_un.Rmd`,
  `sun_2025_opportunities_and_challenges_of_local_ancestry_in_geneti.Rmd`.
- **Population resources:** `romay_2013_comprehensive_genotyping_of_the_usa_national_maize_inb.Rmd`,
  `flintgarcia_2005_maize_association_population_a_highresolution_pl.Rmd`,
  `yu_2008_maize_nam.Rmd`, `1001_genomes_2016_arabidopsis.Rmd`,
  `3k_rice_genomes_2014.Rmd`, `paper_2015_a_global_reference_for_human_genetic_variatio.Rmd`.

## Phenotypes, Breeding, and Genomic Selection

- **Trait architecture and phenotype simulation:** `fernandes_2020_simplephenotypes.Rmd`,
  `panigrahi_2024_cattle_reproduction_qtls.Rmd`, `nguyen_2024_aquaculture_infectious_disease_genomics.Rmd`.
- **Genomic prediction and selection:** `meuwissen_2001_genomic_selection.Rmd`,
  `vanraden_2008_efficient_methods_to_compute_genomic_predicti.Rmd`,
  `endelman_2011_rrblup.Rmd`, `perez_2014_genome_wide_regression_and_prediction_with_th.Rmd`,
  `gianola_2006_genomic_assisted_prediction_of_genetic_value.Rmd`,
  `heffner_2009_genomic_selection_for_crop_improvement.Rmd`,
  `alemu_2024_genomic_selection_plant_breeding.Rmd`,
  `sinha_2023_integrated_genomic_selection_cereals.Rmd`,
  `yanez_2023_gwas_genomic_selection_aquaculture.Rmd`,
  `liu_2016_farmcpu.Rmd`.
- **Breeding/population simulation:** `peng_2005_simupop.Rmd`,
  `gaynor_2021_alphasimr.Rmd`, `haller_2019_slim3.Rmd`,
  `baumdicker_2022_msprime.Rmd`, `kelleher_2016_efficient_coalescent_simulation_and_genealogi.Rmd`,
  `kelleher_2019_inferring_whole_genome_histories_in_large_pop.Rmd`.

## Bulk RNA-seq, Differential Expression, and eQTLs

- **RNA-seq simulation and quantification:** `frazee_2015_polyester.Rmd`,
  `patro_2017_salmon_provides_fast_and_bias_aware_quantific.Rmd`,
  `dobin_2013_star_ultrafast_universal_rna_seq_aligner.Rmd`,
  `kim_2015_hisat_a_fast_spliced_aligner_with_low_memory.Rmd`,
  `liao_2014_featurecounts_an_efficient_general_purpose_pr.Rmd`.
- **RNA-seq analysis and differential expression:** `conesa_2016_a_survey_of_best_practices_for_rna_seq_data_a.Rmd`,
  `love_2014_moderated_estimation_of_fold_change_and_dispe.Rmd`,
  `robinson_2010_tt_edger_tt_a_bioconductor_package_for_differ.Rmd`,
  `law_2014_voom_precision_weights_unlock_linear_model_an.Rmd`,
  `trapnell_2013_differential_analysis_of_gene_regulation_at_t.Rmd`,
  `geraci_2020_rnaseq_editorial.Rmd`,
  `dawadi_2025_rnaseq_practical_guide.Rmd`.

## Single-Cell RNA-seq and Multi-Omics

- **Simulation:** `zappia_2017_splatter.Rmd`, `azodi_2021_splatpop_simulating_population_scale_single_cell_rna_s.Rmd`,
  `song_2023_scdesign3.Rmd`, `vieth_2017_powsimr.Rmd`.
- **Analysis, integration, and QC:** `nesari_2026_scrna_analysis_review.Rmd`,
  `su_2022_scrna_biomedical_guidelines.Rmd`, `butler_2018_integrating_single_cell_transcriptomic_data_a.Rmd`,
  `stuart_2019_comprehensive_integration_of_single_cell_data.Rmd`,
  `hafemeister_2019_normalization_and_variance_stabilization_of_s.Rmd`,
  `korsunsky_2019_fast_sensitive_and_accurate_integration_of_si.Rmd`,
  `jiang_2022_scrna_zero_inflation.Rmd`.
- **Joint chromatin/transcriptome:** `cao_2018_joint_profiling_of_chromatin_accessibility_an.Rmd`.

## ChIP-seq, ATAC-seq, and Regulatory Genomics

- **ChIP-seq and peak calling:** `zhang_2008_model_based_analysis_of_chip_seq_macs.Rmd`.
- **ATAC-seq and accessibility:** `buenrostro_2013_transposition_of_native_chromatin_for_fast_an.Rmd`,
  `schep_2017_chromvar_inferring_transcription_factor_associated_acc.Rmd`.
- **Methylation:** `krueger_2011_bismark_a_flexible_aligner_and_methylation_ca.Rmd`.

## Sequencing, Assembly, and Variant Processing

- **Read simulation:** `huang_2012_art.Rmd`, `ono_2013_pbsim.Rmd`,
  `wick_2019_badread.Rmd`, `scarano_2024_third_generation_sequencing.Rmd`.
- **Assembly and evaluation:** `wick_2017_unicycler.Rmd`,
  `bankevich_2012_spades_a_new_genome_assembly_algorithm_and_it.Rmd`,
  `kolmogorov_2019_assembly_of_long_error_prone_reads_using_repe.Rmd`,
  `koren_2017_canu_scalable_and_accurate_long_read_assembly.Rmd`,
  `vaser_2017_fast_and_accurate_de_novo_genome_assembly_fro.Rmd`,
  `gurevich_2013_quast_quality_assessment_tool_for_genome_asse.Rmd`.
- **Variant/data formats and QC:** `li_2009_the_sequence_alignment_map_format_and_samtool.Rmd`,
  `danecek_2011_the_variant_call_format_and_vcftools.Rmd`,
  `danecek_2021_twelve_years_of_samtools_and_bcftools.Rmd`,
  `mckenna_2010_the_genome_analysis_toolkit_a_mapreduce_frame.Rmd`,
  `cingolani_2012_a_program_for_annotating_and_predicting_the_e.Rmd`,
  `chen_2018_fastp_an_ultra_fast_all_in_one_fastq_preproce.Rmd`.

## Species and Pangenomes

- **Crop genomes and pangenomes:** `huang_2022_integrated_genomics_crop_domestication.Rmd`,
  `shi_2023_plant_pan_genomics.Rmd`, `walkowiak_2020_multiple_wheat_genomes_reveal_global_variation_in.Rmd`,
  `allen_2017_characterization_of_a_wheat_breeders_array_suitable_fo.Rmd`,
  `xu_2022_smart_breeding_igep.Rmd`,
  `petit_2016_tomato_gpat6.Rmd`, `khan_2026_cuticle_drought_tolerance.Rmd`,
  `jung_2006_rice_wda1.Rmd`.
- **Bacterial analysis context:** `bradbury_2007_tassel_software_for_association_mapping_of_complex.Rmd` is not bacterial-specific; bacterial-GWAS questions require a dedicated bacterial workflow/tool boundary rather than blindly applying plant or human GWAS defaults.
