#!/usr/bin/env python3
"""Batch 1 of the SimItAll prompt -> StudySpec evaluation set (draft, needs review).

Covers breeding/cross design, pedigree inbreeding, phenotype outcome models,
GWAS (quantitative, binary, multi-trait, rare-variant, relatedness calibration),
linkage QTL, genomic prediction, multi-step combinations, traps and
not-yet-supported techniques.

Run from the repo root:
    python3 data-raw/agent_eval/batch1_genetics.py
It (re)writes the batch-1 rows of inst/agent/prompt_eval/study_spec_examples_draft.jsonl.

Field conventions (see inst/agent/prompt_eval/README.md):
  expected_spec  only the fields the prompt determines; the scorer checks just these.
  must_ask       dotted fields the agent must NOT fill silently (ask, or disclose as an assumption).
  must_not       analyses that must not appear anywhere in the plan.
"""
import json, os

ROWS = []


def ex(id, workflow, situation, style, prompt, spec, must_ask=(), must_not=(), source="", notes=""):
    s = {"spec_version": "1.0"}
    s.update(spec)
    s.setdefault("secondary_analyses", [])
    ROWS.append({
        "id": id, "batch": 1, "review_status": "draft", "workflow": workflow,
        "situation": situation, "style": style, "source": source, "prompt": prompt,
        "expected_spec": s, "must_ask": list(must_ask), "must_not": list(must_not), "notes": notes,
    })


SIM, AN, BOTH = "simulate", "analyze", "simulate_and_analyze"
SYN, PROV, ASK, DEMO = "synthetic", "provided", "ask", "package_demo"
T4 = "docs/TUTORIAL.md sec 4"
T3 = "docs/TUTORIAL.md sec 3"

# ------------------------------------------------------------------ cross design (tutorial 4.3-4.6)
ex("cross_ril_8self", "cross_design", "simulate", "expert",
   "Simulate 200 RILs from two founders with eight generations of selfing.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"scheme": "RIL", "founder_count": 2, "selfing_generations": 8, "population_size": 200}},
   must_ask=["organism"], source=T4 + ".4")
ex("cross_f2_maize", "cross_design", "simulate", "beginner",
   "can you make me an F2 population of 250 maize plants from a cross of two inbreds?",
   {"goal": SIM, "organism": "maize", "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"scheme": "F2", "founder_count": 2, "population_size": 250}}, source=T4 + ".3")
ex("cross_nam_5founders", "cross_design", "simulate", "expert",
   "Build a synthetic NAM population: one common parent crossed to four other founders, 400 offspring in total.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"scheme": "NAM", "founder_count": 5, "population_size": 400}},
   must_ask=["organism"], source=T4 + ".6", notes="1 common + 4 = 5 founders, as in the tutorial hap1..hap5 example.")
ex("cross_magic8", "cross_design", "simulate", "terse",
   "MAGIC, 8 founders, 300 lines, 6 rounds of selfing, rice.",
   {"goal": SIM, "organism": "rice", "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"scheme": "MAGIC", "founder_count": 8, "selfing_generations": 6, "population_size": 300}}, source=T4 + ".6")
ex("cross_dh", "cross_design", "simulate", "chatty",
   "We're thinking about doubled haploids instead of RILs. Could you simulate 200 DH lines from a biparental wheat cross so I can see what they look like?",
   {"goal": SIM, "organism": "wheat", "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"scheme": "DH", "founder_count": 2, "population_size": 200}},
   notes="Catalog scheme enum lacks DH; tutorial supports scheme='DH'. Reviewer: extend enum or map to biparental.", source=T4 + ".4")
ex("cross_bc2s7_chr10", "cross_design", "simulate", "expert",
   "Simulate a synthetic biparental cross for chromosome 10 with 2 backcrosses and 7 selfing generations, then show final allele frequencies.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "cross_design", "region": "chr10",
    "design": {"scheme": "biparental", "founder_count": 2, "backcross_generations": 2, "selfing_generations": 7}},
   source="docs/TUTORIAL.md 'Agent worked example: biparental BC2S7'", notes="Tutorial's own agent worked example; expected donor frequency 0.125.")
ex("cross_nil", "cross_design", "simulate", "expert",
   "Make 100 near-isogenic lines: four backcrosses to the recurrent parent followed by three selfing generations.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"scheme": "NIL", "founder_count": 2, "backcross_generations": 4, "selfing_generations": 3, "population_size": 100}},
   must_ask=["organism"], source=T4 + ".5")
ex("cross_f23", "cross_design", "simulate", "nonnative",
   "I need population of F2:3 families, 200 plant, from soybean two parent. please simulate.",
   {"goal": SIM, "organism": "soybean", "data_source": SYN, "primary_analysis": "cross_design",
    "design": {"founder_count": 2, "selfing_generations": 2, "population_size": 200}}, source=T4 + ".3",
   notes="Tutorial: F2:3-style = sequence F1,SELF:2.")
ex("cross_missing_scheme", "cross_design", "ask", "beginner",
   "Simulate a breeding population of sorghum for me.",
   {"goal": SIM, "organism": "sorghum", "data_source": SYN, "primary_analysis": "cross_design"},
   must_ask=["design.scheme", "design.population_size", "design.founder_count"],
   notes="No scheme or size given: must ask, not default to NAM.")
ex("cross_my_founders", "cross_design", "provided", "expert",
   "Using my phased founder FASTA (4 founders, 3 chromosomes), simulate an F2 with 100 offspring per family.",
   {"goal": SIM, "data_source": PROV, "primary_analysis": "cross_design",
    "design": {"scheme": "F2", "founder_count": 4}},
   must_ask=["organism"], source=T4 + ".2")

# ------------------------------------------------------------------ pedigree inbreeding (planned)
ex("ped_selfing_10gen", "pedigree_inbreeding", "simulate", "expert",
   "Simulate 10 generations of selfing in a population of 500 Arabidopsis plants and track the inbreeding coefficient.",
   {"goal": SIM, "organism": "arabidopsis", "data_source": SYN, "primary_analysis": "pedigree_inbreeding",
    "parameters": {"mating_scheme": "selfing", "generations": 10, "population_size": 500}})
ex("ped_sib_mice", "pedigree_inbreeding", "simulate", "chatty",
   "Our mouse colony is kept by brother-sister mating. Can you simulate 20 generations of that with 50 breeding pairs and show how heterozygosity drops?",
   {"goal": SIM, "organism": "mouse", "data_source": SYN, "primary_analysis": "pedigree_inbreeding",
    "parameters": {"mating_scheme": "sib_mating", "generations": 20}},
   must_ask=["parameters.population_size"], notes="50 pairs = 100 individuals; reviewer to decide whether population_size should be 100.")
ex("ped_first_cousin_human", "pedigree_inbreeding", "simulate", "expert",
   "Simulate a synthetic human chr10 pedigree with first-cousin marriages for 10 generations and compute pedigree F and marker homozygosity.",
   {"goal": SIM, "organism": "human", "data_source": SYN, "region": "chr10", "primary_analysis": "pedigree_inbreeding",
    "parameters": {"mating_scheme": "first_cousin", "generations": 10}},
   must_ask=["parameters.population_size"], source="AGENTS.md planning-model example")
ex("ped_no_scheme", "pedigree_inbreeding", "ask", "beginner",
   "Simulate inbreeding in a small cattle herd over 15 generations.",
   {"goal": SIM, "organism": "cattle", "data_source": SYN, "primary_analysis": "pedigree_inbreeding",
    "parameters": {"generations": 15}},
   must_ask=["parameters.mating_scheme", "parameters.population_size"], source="AGENTS.md: 'inbreeding' must not imply one mating scheme")
ex("ped_chr14_roh", "pedigree_inbreeding", "provided", "expert",
   "For my phased human chr14 VCF, recombination map and GRCh38 GTF, simulate ten generations of inbreeding and rank the most affected genes.",
   {"goal": BOTH, "organism": "human", "data_source": PROV, "region": "chr14", "primary_analysis": "pedigree_inbreeding",
    "secondary_analyses": ["gene_coordinate_annotation"], "parameters": {"generations": 10}, "execution_status": "unsupported", "do_not_run": True},
   must_ask=["parameters.mating_scheme"], must_not=["quantitative_gwas", "binary_gwas"],
   source="AGENTS.md required regression case", notes="Must not emit generic chr10 family-risk recipe.")
ex("ped_custom", "pedigree_inbreeding", "provided", "terse",
   "Here's my pedigree file for 300 dogs. Estimate inbreeding coefficients.",
   {"goal": AN, "organism": "dog", "data_source": PROV, "primary_analysis": "pedigree_inbreeding",
    "parameters": {"mating_scheme": "custom_pedigree"}})

# ------------------------------------------------------------------ phenotype outcome models (tutorial 4.9)
ex("pheno_h2_gauss", "phenotype_outcome_model", "simulate", "expert",
   "Simulate a continuous trait with heritability 0.4 and 20 additive QTN on my maize genotype file.",
   {"goal": SIM, "organism": "maize", "data_source": PROV, "primary_analysis": "phenotype_outcome_model",
    "parameters": {"outcome_family": "gaussian", "heritability": 0.4, "n_causal": 20}}, source=T3)
ex("pheno_binary_liability", "phenotype_outcome_model", "simulate", "chatty",
   "I'd like a disease-like yes/no trait on top of simulated genotypes, about 30% affected, liability heritability 0.5.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "phenotype_outcome_model",
    "parameters": {"outcome_family": "liability", "heritability": 0.5, "prevalence": 0.3}}, must_ask=["organism"], source="docs/TUTORIAL.md sec 4.9")
ex("pheno_ordinal", "phenotype_outcome_model", "simulate", "terse",
   "Ordinal disease score 1-5 for 300 sorghum lines, h2 0.3.",
   {"goal": SIM, "organism": "sorghum", "data_source": SYN, "primary_analysis": "phenotype_outcome_model",
    "design": {"population_size": 300}, "parameters": {"outcome_family": "ordinal", "heritability": 0.3, "n_categories": 5}}, source="docs/TUTORIAL.md sec 4.9")
ex("pheno_counts", "phenotype_outcome_model", "simulate", "expert",
   "Simulate tiller counts with overdispersion (negative binomial) for 250 rice lines.",
   {"goal": SIM, "organism": "rice", "data_source": SYN, "primary_analysis": "phenotype_outcome_model",
    "design": {"population_size": 250}, "parameters": {"outcome_family": "negative_binomial"}},
   must_ask=["parameters.heritability"], source="docs/TUTORIAL.md sec 4.9")
ex("pheno_maternal", "phenotype_outcome_model", "simulate", "expert",
   "Simulate seed weight with a maternal effect and additive genetics in F2 families; h2 = 0.35, maternal variance 20%.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "phenotype_outcome_model", "secondary_analyses": ["cross_design"],
    "design": {"scheme": "F2"}, "parameters": {"outcome_family": "gaussian", "heritability": 0.35, "maternal_variance": 0.2}},
   must_ask=["organism"], source="docs/TUTORIAL.md sec 4.9")
ex("pheno_no_h2", "phenotype_outcome_model", "ask", "beginner",
   "make up a phenotype for my vcf",
   {"goal": SIM, "data_source": PROV, "primary_analysis": "phenotype_outcome_model"},
   must_ask=["parameters.outcome_family", "parameters.heritability"])

# ------------------------------------------------------------------ quantitative GWAS (tutorial 3, 3.1)
ex("gwas_sim_run_bench", "quantitative_gwas", "simulate", "expert",
   "Simulate a GWAS cohort of 200 samples from 2 subpopulations with 10 causal SNPs, run a mixed-model GWAS with 2 PCs, and benchmark recovery of the causal loci.",
   {"goal": BOTH, "data_source": SYN, "primary_analysis": "quantitative_gwas",
    "design": {"population_size": 200, "n_subpopulations": 2},
    "parameters": {"association_model": "linear_mixed_model", "n_PCs": 2, "n_causal": 10, "trait_type": "continuous"}},
   must_ask=["organism", "parameters.heritability"], must_not=["bulk_counts_simulation", "genomic_prediction"], source=T3 + ".1")
ex("gwas_my_vcf_farmcpu", "quantitative_gwas", "provided", "expert",
   "Run FarmCPU on my sorghum VCF and plant height phenotypes, with Bonferroni correction.",
   {"goal": AN, "organism": "sorghum", "data_source": PROV, "primary_analysis": "quantitative_gwas",
    "parameters": {"association_model": "FarmCPU", "multiple_testing": "Bonferroni", "trait_type": "continuous"}},
   must_not=["genomic_prediction", "bulk_counts_simulation"])
ex("gwas_find_gene", "quantitative_gwas", "provided", "beginner",
   "I have genotypes and flowering time for 350 maize inbreds. Which genes control flowering?",
   {"goal": AN, "organism": "maize", "data_source": PROV, "primary_analysis": "quantitative_gwas",
    "design": {"population_size": 350}, "parameters": {"trait_type": "continuous"}},
   must_ask=["parameters.association_model"], must_not=["bulk_counts_simulation", "differential_expression"],
   notes="Gene hunting = GWAS; must not add RNA-seq.")
ex("gwas_nam_quant", "quantitative_gwas", "simulate", "terse",
   "NAM GWAS on simulated maize, quantitative trait, h2 0.6, 15 QTL.",
   {"goal": BOTH, "organism": "maize", "data_source": SYN, "primary_analysis": "quantitative_gwas", "secondary_analyses": ["cross_design"],
    "design": {"scheme": "NAM"}, "parameters": {"heritability": 0.6, "n_causal": 15, "trait_type": "continuous"}},
   source="inst/agent/knowledge/workflows/maize_nam_quantitative_gwas.Rmd")
ex("gwas_lmm_kinship_pcs", "quantitative_gwas", "provided", "expert",
   "GWAS with an LMM, kinship plus 3 PCs, FDR 5%, on my rice panel of 413 accessions.",
   {"goal": AN, "organism": "rice", "data_source": PROV, "primary_analysis": "quantitative_gwas", "design": {"population_size": 413},
    "parameters": {"association_model": "linear_mixed_model", "structure_adjustment": ["PCs", "kinship"], "n_PCs": 3, "multiple_testing": "FDR"}})
ex("gwas_demo", "quantitative_gwas", "demo", "beginner",
   "Show me a GWAS example using the package demo data.",
   {"goal": AN, "data_source": DEMO, "primary_analysis": "quantitative_gwas"})
ex("gwas_nonnative", "quantitative_gwas", "provided", "nonnative",
   "I am have the vcf of wheat 500 line and grain yield, want association study please with population structure correction.",
   {"goal": AN, "organism": "wheat", "data_source": PROV, "primary_analysis": "quantitative_gwas", "design": {"population_size": 500},
    "parameters": {"trait_type": "continuous"}}, must_ask=["parameters.structure_adjustment"])
ex("gwas_no_trait", "quantitative_gwas", "ask", "terse",
   "run gwas",
   {"goal": AN, "data_source": ASK, "primary_analysis": "quantitative_gwas"},
   must_ask=["organism", "parameters.trait_type"], notes="Trait type unknown: quantitative vs binary must be asked.")

# ------------------------------------------------------------------ binary GWAS
ex("gwas_binary_nam", "binary_gwas", "simulate", "terse",
   "binomial phenotype GWAS in NAM",
   {"goal": BOTH, "data_source": SYN, "primary_analysis": "binary_gwas", "secondary_analyses": ["cross_design"], "design": {"scheme": "NAM"}},
   must_ask=["organism"], must_not=["quantitative_gwas"],
   source="AGENTS.md product decision: binary-trait GWAS, not a generic F1/F2/F3 demo",
   notes="AGENTS.md: means binary-trait GWAS; no unrelated extra analyses.")
ex("gwas_case_control", "binary_gwas", "provided", "expert",
   "Case-control GWAS of resistance (1) vs susceptibility (0) in my 280 barley lines using Firth logistic regression.",
   {"goal": AN, "organism": "barley", "data_source": PROV, "primary_analysis": "binary_gwas", "design": {"population_size": 280},
    "parameters": {"association_model": "Firth_logistic"}})
ex("gwas_binary_sim", "binary_gwas", "simulate", "chatty",
   "Could you simulate 1,000 cases and 1,000 controls with 5 causal variants and then run a logistic mixed model GWAS?",
   {"goal": BOTH, "data_source": SYN, "primary_analysis": "binary_gwas",
    "design": {"cases": 1000, "controls": 1000}, "parameters": {"association_model": "logistic_mixed_model", "n_causal": 5}},
   must_ask=["organism"])
ex("gwas_binary_maize_nam_wf", "binary_gwas", "simulate", "expert",
   "Simulate maize NAM families with a binary disease trait (prevalence 0.2) and run a binary GWAS.",
   {"goal": BOTH, "organism": "maize", "data_source": SYN, "primary_analysis": "binary_gwas", "secondary_analyses": ["cross_design"],
    "design": {"scheme": "NAM"}, "parameters": {"prevalence": 0.2}}, source="inst/agent/knowledge/workflows/maize_nam_binary_gwas.Rmd")

# ------------------------------------------------------------------ other GWAS (planned / knowledge)
ex("gwas_multi_env", "multitrait_longitudinal_gwas", "provided", "expert",
   "GxE GWAS of yield measured in 4 environments on 300 soybean lines.",
   {"goal": AN, "organism": "soybean", "data_source": PROV, "primary_analysis": "multitrait_longitudinal_gwas",
    "design": {"population_size": 300, "environments": 4}, "parameters": {"design_type": "GxE"}})
ex("gwas_longitudinal", "multitrait_longitudinal_gwas", "provided", "chatty",
   "We measured plant height every week for 8 weeks on the same 200 lines. Can we do a GWAS on the growth curves?",
   {"goal": AN, "data_source": PROV, "primary_analysis": "multitrait_longitudinal_gwas",
    "design": {"population_size": 200, "time_points": 8}, "parameters": {"design_type": "longitudinal"}}, must_ask=["organism"])
ex("gwas_rare_gene", "rare_variant_and_gene_based_gwas", "provided", "expert",
   "Gene-based burden test of rare variants (MAF < 1%) from my exome VCF against a binary phenotype.",
   {"goal": AN, "data_source": PROV, "primary_analysis": "rare_variant_and_gene_based_gwas",
    "parameters": {"variant_set_definition": "gene", "MAF_threshold": 0.01}}, must_ask=["organism"])
ex("gwas_relatedness_calib", "population_relatedness_gwas_calibration", "simulate", "expert",
   "Show how cryptic relatedness inflates GWAS p-values and how kinship correction fixes it, using simulated human data.",
   {"goal": BOTH, "organism": "human", "data_source": SYN, "primary_analysis": "population_relatedness_gwas_calibration",
    "parameters": {"adjustment": "kinship"}}, source="inst/agent/knowledge/workflows/human_relatedness_gwas_calibration.Rmd")

# ------------------------------------------------------------------ linkage QTL (planned)
ex("qtl_ril_reviewed", "linkage_qtl", "simulate", "expert",
   "Simulate 200 maize RILs from two founders, selfed for six generations, then map a drought trait.",
   {"goal": BOTH, "organism": "maize", "data_source": SYN, "primary_analysis": "linkage_qtl",
    "secondary_analyses": ["cross_design", "phenotype_outcome_model"],
    "design": {"scheme": "RIL", "founder_count": 2, "selfing_generations": 6, "population_size": 200}, "parameters": {"outcome_family": "gaussian"}},
   source="inst/agent/capabilities/study_spec_examples.jsonl (reviewed)", notes="Copied from the reviewed set as an anchor.")
ex("qtl_f2_cim", "linkage_qtl", "provided", "expert",
   "Composite interval mapping on my F2 tomato population (genotypes, map and fruit weight), 1000 permutations.",
   {"goal": AN, "organism": "tomato", "data_source": PROV, "primary_analysis": "linkage_qtl", "design": {"scheme": "F2"},
    "parameters": {"mapping_method": "composite_interval_mapping", "permutations": 1000}})
ex("qtl_no_map", "linkage_qtl", "provided", "beginner",
   "I have 150 RILs genotyped and phenotyped for height. Find the QTL.",
   {"goal": AN, "data_source": PROV, "primary_analysis": "linkage_qtl", "design": {"scheme": "RIL", "population_size": 150}},
   must_ask=["organism", "parameters.mapping_method"], notes="Genetic map is a required input: must request it.")
ex("qtl_dh_barley", "linkage_qtl", "simulate", "terse",
   "Sim 180 barley DH lines + QTL scan for heading date.",
   {"goal": BOTH, "organism": "barley", "data_source": SYN, "primary_analysis": "linkage_qtl", "secondary_analyses": ["cross_design"],
    "design": {"scheme": "DH", "population_size": 180}})
ex("qtl_vs_gwas_ambiguous", "linkage_qtl", "ask", "beginner",
   "I want to map genes for drought tolerance in sorghum.",
   {"goal": AN, "organism": "sorghum", "data_source": ASK}, must_ask=["primary_analysis"],
   notes="Linkage (biparental population) vs association (diversity panel) depends on the population: ask.")

# ------------------------------------------------------------------ genomic prediction (tutorial 3.2)
ex("gs_reviewed", "genomic_prediction", "provided", "expert",
   "Use my genotype and yield data to predict GEBVs and choose the top 10 parents.",
   {"goal": AN, "data_source": PROV, "primary_analysis": "genomic_prediction", "parameters": {"model": "GBLUP", "selection_count": 10}},
   source="inst/agent/capabilities/study_spec_examples.jsonl (reviewed)", notes="Reviewed example sets GBLUP although unstated; reviewer: confirm or move to must_ask.")
ex("gs_cv_bayesb", "genomic_prediction", "provided", "expert",
   "Compare GBLUP and BayesB prediction accuracy with 5-fold cross-validation on my wheat data.",
   {"goal": AN, "organism": "wheat", "data_source": PROV, "primary_analysis": "genomic_prediction",
    "parameters": {"validation_design": "random_kfold", "folds": 5}}, notes="Two models requested; model field should hold both.")
ex("gs_cycle_sim", "genomic_prediction", "simulate", "expert",
   "Simulate one round of genomic selection: fit GBLUP, pick 10 parents, make 20 crosses with minimum-kinship mating.",
   {"goal": BOTH, "data_source": SYN, "primary_analysis": "genomic_prediction", "secondary_analyses": ["cross_design"],
    "parameters": {"model": "GBLUP", "selection_count": 10, "n_crosses": 20, "mating": "minimum_kinship"}},
   must_ask=["organism"], source=T3 + ".2", notes="tutorial run_genomic_selection(n_parents=10, n_crosses=20)")
ex("gs_forward", "genomic_prediction", "provided", "chatty",
   "We have three breeding cycles of maize data. Can you train on cycles 1-2 and see how well it predicts cycle 3?",
   {"goal": AN, "organism": "maize", "data_source": PROV, "primary_analysis": "genomic_prediction",
    "parameters": {"validation_design": "forward_generation"}}, must_ask=["parameters.model"])
ex("gs_env_held_out", "genomic_prediction", "provided", "terse",
   "rrBLUP, leave-one-environment-out, sorghum biomass.",
   {"goal": AN, "organism": "sorghum", "data_source": PROV, "primary_analysis": "genomic_prediction",
    "parameters": {"model": "rrBLUP", "validation_design": "environment_held_out"}})
ex("gs_which_parents", "genomic_prediction", "provided", "beginner",
   "which lines should I cross next season? I have genotypes and yield for 400 lines",
   {"goal": AN, "data_source": PROV, "primary_analysis": "genomic_prediction", "design": {"population_size": 400}},
   must_ask=["organism", "parameters.model"], must_not=["quantitative_gwas"], source="agent_spec.R genomic_selection phrasing")

# ------------------------------------------------------------------ combinations
ex("combo_cross_gwas_gs", "combination", "simulate", "expert",
   "Simulate a 500-line maize NAM, run GWAS for yield, then fit GBLUP and select the top 20 lines.",
   {"goal": BOTH, "organism": "maize", "data_source": SYN, "primary_analysis": "quantitative_gwas",
    "secondary_analyses": ["cross_design", "genomic_prediction"], "design": {"scheme": "NAM", "population_size": 500},
    "parameters": {"model": "GBLUP", "selection_count": 20}}, must_ask=["parameters.heritability"])
ex("combo_gwas_only_not_rnaseq", "combination", "provided", "expert",
   "Run a GWAS on my VCF and phenotype file. Don't do anything else.",
   {"goal": AN, "data_source": PROV, "primary_analysis": "quantitative_gwas"},
   must_ask=["parameters.trait_type"], must_not=["bulk_counts_simulation", "genomic_prediction", "chipseq_atacseq", "single_cell_eqtl"],
   source="AGENTS.md: GWAS must not silently run other workflows")
ex("combo_pheno_then_qtl", "combination", "simulate", "chatty",
   "First make an F2 of 300 tomato plants, then simulate fruit weight with h2 0.5, and finally do interval mapping.",
   {"goal": BOTH, "organism": "tomato", "data_source": SYN, "primary_analysis": "linkage_qtl",
    "secondary_analyses": ["cross_design", "phenotype_outcome_model"], "design": {"scheme": "F2", "population_size": 300},
    "parameters": {"heritability": 0.5, "mapping_method": "interval_mapping"}})
ex("combo_gwas_then_annot", "combination", "provided", "expert",
   "GWAS on my Arabidopsis VCF for flowering time, then annotate significant SNPs with nearby genes and GO terms.",
   {"goal": AN, "organism": "arabidopsis", "data_source": PROV, "primary_analysis": "quantitative_gwas",
    "secondary_analyses": ["gene_coordinate_annotation"], "parameters": {"trait_type": "continuous"}})

# ------------------------------------------------------------------ traps
ex("trap_real_population", "trap", "ask", "expert",
   "Simulate the real maize NAM population with the actual 25 founders.",
   {"goal": SIM, "organism": "maize", "primary_analysis": "cross_design", "design": {"scheme": "NAM", "founder_count": 25}, "do_not_run": True},
   must_ask=["data_source"], source="docs/TUTORIAL.md: named maize NAM founders are not bundled",
   notes="Must ask for real founder files; must not present synthetic founders as the real NAM.")
ex("trap_h2_impossible", "trap", "ask", "expert",
   "Simulate a trait with heritability 1.4 on 200 lines.",
   {"goal": SIM, "primary_analysis": "phenotype_outcome_model", "design": {"population_size": 200}, "do_not_run": True},
   must_ask=["parameters.heritability"], notes="h2 must be in [0,1]: flag and ask.")
ex("trap_conflicting_n", "trap", "ask", "chatty",
   "Simulate 200 RILs, actually make it 300, from two founders.",
   {"goal": SIM, "data_source": SYN, "primary_analysis": "cross_design", "design": {"scheme": "RIL", "founder_count": 2}},
   must_ask=["design.population_size"], notes="Conflicting sizes: confirm 300 rather than silently picking one.")
ex("trap_clinical", "trap", "ask", "chatty",
   "Use my family's genotypes to tell me my child's risk of cleft lip.",
   {"goal": AN, "organism": "human", "data_source": PROV, "do_not_run": True},
   notes="Clinical prediction is out of scope: synthetic research demonstrations only (AGENTS.md).")
ex("trap_selfing_mice", "trap", "ask", "expert",
   "Simulate 10 generations of selfing in mice.",
   {"goal": SIM, "organism": "mouse", "primary_analysis": "pedigree_inbreeding", "do_not_run": True},
   must_ask=["parameters.mating_scheme"], notes="Mice cannot self: flag biological infeasibility (AGENTS species rules).")
ex("trap_gwas_10_samples", "trap", "ask", "beginner",
   "GWAS with my 10 plants",
   {"goal": AN, "data_source": PROV, "primary_analysis": "quantitative_gwas", "design": {"population_size": 10}},
   notes="Should warn that 10 samples has essentially no power; still a plan, not silent run.")

# ------------------------------------------------------------------ known but not yet supported techniques
ex("unsup_riboseq", "unsupported", "unsupported", "expert",
   "Simulate Ribo-seq footprints for heat vs control with matched RNA-seq.",
   {"goal": SIM, "data_source": SYN, "execution_status": "unsupported", "do_not_run": True},
   must_not=["bulk_counts_simulation"], notes="Ribo-seq is not RNA-seq: must not silently route to bulk counts simulation.")
ex("unsup_hic", "unsupported", "unsupported", "terse",
   "simulate Hi-C contact maps for chr1",
   {"goal": SIM, "data_source": SYN, "execution_status": "unsupported", "do_not_run": True},
   must_not=["chipseq_atacseq"])
ex("unsup_methylation", "unsupported", "unsupported", "chatty",
   "Could you simulate bisulfite sequencing data so I can test DMR calling?",
   {"goal": SIM, "data_source": SYN, "execution_status": "unsupported", "do_not_run": True},
   must_not=["read_simulation", "chipseq_atacseq"])
ex("unsup_proteomics", "unsupported", "unsupported", "expert",
   "Simulate label-free proteomics intensities for 3 conditions with 4 replicates.",
   {"goal": SIM, "data_source": SYN, "execution_status": "unsupported", "do_not_run": True},
   must_not=["bulk_counts_simulation"])

# ------------------------------------------------------------------ splits and output
# Deterministic split: every third row goes to dev, the rest to the frozen test set.
for i, r in enumerate(ROWS):
    r["split"] = "dev" if i % 3 == 0 else "test"

out = "inst/agent/prompt_eval/study_spec_examples_draft.jsonl"
os.makedirs(os.path.dirname(out), exist_ok=True)
keep = []
if os.path.exists(out):
    keep = [l for l in open(out) if json.loads(l).get("batch") != 1]
ids = [r["id"] for r in ROWS]
assert len(ids) == len(set(ids)), "duplicate ids"
with open(out, "w") as f:
    f.writelines(keep)
    for r in ROWS:
        f.write(json.dumps(r, ensure_ascii=False) + "\n")
print(f"batch 1: {len(ROWS)} rows ({sum(r['split'] == 'dev' for r in ROWS)} dev, {sum(r['split'] == 'test' for r in ROWS)} test) -> {out}")
