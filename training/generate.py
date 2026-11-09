#!/usr/bin/env python3
"""Synthetic training data for the SimItAll prompt model (labels by construction).

Each sentence is assembled from labelled parts, so every token tag and every
sentence label (goal, data source, primary analysis) is known exactly. The
vocabulary (organisms, strains/panels/breeds, technique names, method words)
comes from the lexicon tables in inst/agent/lexicon, so the model learns the
same names the lexicon recognises, in many more contexts.

Sentences that look too much like any prompt in the evaluation set
(inst/agent/prompt_eval/) are dropped, so the frozen test split stays unseen.

    python training/generate.py --n 60000 --seed 1 --out training/data/train.jsonl
    python training/generate.py --n 3000  --seed 2 --out training/data/dev.jsonl
"""
import argparse, collections, csv, json, os, random, re, sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "inst", "python"))
from prompt_model import ANALYSES, LABELS, tokenize  # noqa: E402

LEX = os.path.join(ROOT, "inst", "agent", "lexicon")
R = random.Random(0)


def tsv(name):
    with open(os.path.join(LEX, name)) as f:
        return list(csv.DictReader(f, delimiter="\t"))


# ---------------------------------------------------------------- vocabulary
ORG = {o["organism"]: o for o in tsv("organisms.tsv")}
ORG_NAMES = collections.defaultdict(list)    # organism -> names usable as "in {name}" / "{name} {unit}"
ORG_PLURAL = collections.defaultdict(list)
ORG_PANEL = collections.defaultdict(list)    # strains, panels, breeds, cultivars
for a in tsv("organism_aliases.tsv"):
    if a["strength"] != "strong":
        continue
    t, o, al = a["alias_type"], a["organism"], a["alias"]
    if t in ("common_name", "scientific_name", "abbreviation") or (t == "canonical_id" and "_" not in al):
        ORG_NAMES[o].append(al)
    elif t in ("plural", "informal") and al.lower().endswith("s"):
        ORG_PLURAL[o].append(al)
    elif t == "informal":
        ORG_NAMES[o].append(al)
    elif t in ("strain", "panel", "breed"):
        ORG_PANEL[o].append(al)
GROUP = collections.defaultdict(list)
for o, r in ORG.items():
    if ORG_NAMES[o] and o != "bacteria":
        GROUP[r["group"]].append(o)
PLANTS, ANIMALS = GROUP["plant"], GROUP["animal"]
MICROBES = GROUP["bacterium"] + GROUP["fungus"] + GROUP["protist"]
SELFERS = [o for o in PLANTS if ORG[o]["can_self"] == "yes"] + ["c_elegans"]
LAB_ANIMALS = [o for o in ("mouse", "rat", "drosophila", "zebrafish", "medaka", "c_elegans") if o in ORG]
LIVESTOCK = [o for o in ("cattle", "pig", "sheep", "goat", "chicken", "dog", "horse", "cat") if o in ORG]
CROPS = [o for o in ("maize", "rice", "sorghum", "wheat", "barley", "soybean", "tomato", "arabidopsis", "potato",
                     "cotton", "rapeseed", "common_bean", "pearl_millet", "foxtail_millet", "oat", "chickpea") if o in ORG]

TECH = {t["technique"]: t for t in tsv("techniques.tsv")}
TECH_ALIAS = collections.defaultdict(list)
for a in tsv("technique_aliases.tsv"):
    if a["strength"] == "strong" and a["source"] in ("curated", "analysis_catalog") and len(a["alias"]) <= 40:
        TECH_ALIAS[a["technique"]].append(a["alias"])
UNSUPPORTED = [t for t, r in TECH.items() if r["support"] == "not_supported" and TECH_ALIAS[t]]
SETTING = collections.defaultdict(list)    # (analysis, field) -> [(alias, value)]
for s in tsv("settings.tsv"):
    if s["strength"] == "strong":
        SETTING[(s["analysis"], s["field"])].append(s["alias"])


def org_units(o):
    g = ORG[o]["group"]
    if o == "human":
        return ["people", "individuals", "participants", "samples", "patients"]
    if g == "plant":
        return ["lines", "plants", "accessions", "genotypes", "inbreds", "individuals", "entries"]
    if g == "animal":
        return ORG_PLURAL[o][:3] + ["animals", "individuals"] if ORG_PLURAL[o] else ["animals", "individuals"]
    return ["strains", "isolates", "samples"]


# ---------------------------------------------------------------- sentence builder
class S:
    def __init__(self):
        self.p = []

    def __call__(self, text, tag=None):
        if text:
            self.p.append((text, tag))
        return self

    def build(self):
        text = " ".join(t for t, _ in self.p)
        text = re.sub(r"\s+([,.?;:!)])", r"\1", text)
        text = re.sub(r"\(\s+", "(", text)
        toks, tags = [], []
        for t, tag in self.p:
            tt = tokenize(t)
            toks += tt
            tags += (["B-" + tag] + ["I-" + tag] * (len(tt) - 1)) if tag else ["O"] * len(tt)
        return text, toks, tags


WORDNUM = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
           10: "ten", 11: "eleven", 12: "twelve", 15: "fifteen", 20: "twenty"}


def fmt(n, words=True):
    if words and n in WORDNUM and R.random() < 0.3:
        return WORDNUM[n]
    if n >= 1000 and R.random() < 0.5:
        return f"{n:,}"
    return str(n)


def ch(*xs):
    return R.choice(xs)


def maybe(p=0.5):
    return R.random() < p


def org_name(o):
    return R.choice(ORG_NAMES[o])


def mention(technique):
    return R.choice(TECH_ALIAS[technique]) if TECH_ALIAS[technique] else technique.replace("_", " ")


def setting(analysis, field):
    xs = SETTING.get((analysis, field))
    return R.choice(xs) if xs else None


# intent openers ------------------------------------------------------------
def sim_open(s):
    s(ch("Simulate", "Please simulate", "Can you simulate", "Could you simulate", "I want to simulate", "Help me simulate",
         "I need to simulate", "Let's simulate", "sim", "Simulate for me", "Generate", "Please generate", "Make me",
         "Can you make", "Create", "I'd like to simulate", "We want to simulate", "Could you generate"))


def prov_open(s, noun):
    op = ch("I have", "We have", "Here is my", "Here's our", "I've got", "We collected", "Using my", "With our", "From my",
            "I already have", "Attached is my", "Our lab has", "We generated", "I uploaded")
    if op.split()[-1] in ("my", "our"):
        noun = re.sub(r"^(?:a|an)\s+", "", noun)
    s(op)(noun)


def ask_end(s):
    s(ch(".", "?", ". Thanks!", ", please.", "", ". What do you suggest?", ". Can you help?", "."))


def distractor(s):
    if maybe(0.25):
        s(ch(", seed 42", f", on chromosome {R.randint(1, 12)}", f" (version {R.randint(2, 6)})", f", deadline in {R.randint(2, 9)} weeks",
             f", ~{R.randint(2, 40)} GB of data", ", output as CSV", f", random seed {R.randint(1, 999)}"))


def org_phrase(s, o, p=0.8):
    if o and maybe(p):
        s(ch("in", "for", "using", "of", "with")).__call__(org_name(o), "ORG")


def pop_phrase(s, o, lo=50, hi=2000):
    s(fmt(R.randint(lo, hi)), "POP")
    if maybe(0.4) and ORG_NAMES[o]:
        s(org_name(o), "ORG")
    s(ch(*org_units(o)))


def conds(s, k=2):
    pairs = [("control", "drought"), ("control", "heat"), ("wild type", "mutant"), ("treated", "untreated"),
             ("control", "salt"), ("mock", "infected"), ("healthy", "disease"), ("control", "treatment"),
             ("low nitrogen", "high nitrogen"), ("vehicle", "drug"), ("wild-type", "knockout"), ("virgin", "mated"),
             ("control", "cold"), ("well-watered", "water-limited"), ("day", "night")]
    a, b = R.choice(pairs)
    if k == 3:
        s(a, "COND")(",")(b, "COND")(ch("and", ","))(ch("recovery", "combined stress", "rewatered"), "COND")
    else:
        s(a, "COND")(ch("vs", "versus", "and", "or", "compared with"))(b, "COND")


CELLTYPES = ["T cells", "B cells", "monocytes", "NK cells", "neurons", "astrocytes", "microglia", "hepatocytes",
             "epidermis", "cortex", "endodermis", "xylem", "phloem", "guard cells", "mesophyll", "enterocytes",
             "stem cells", "fibroblasts", "macrophages", "root hair cells"]
TRAITS_Q = ["plant height", "flowering time", "grain yield", "kernel weight", "days to anthesis", "seed protein",
            "leaf area", "root length", "body weight", "fat mass", "blood glucose", "BMI", "LDL cholesterol",
            "milk yield", "egg weight", "starvation resistance", "wing length", "tiller number", "oil content",
            "lipid content", "growth rate", "chlorophyll content"]
TRAITS_PLANT = ["plant height", "flowering time", "grain yield", "kernel weight", "days to anthesis", "seed protein",
                "leaf area", "root length", "tiller number", "oil content", "lipid content", "chlorophyll content", "fruit weight"]
TRAITS_ANIMAL = ["body weight", "fat mass", "blood glucose", "BMI", "LDL cholesterol", "milk yield", "egg weight",
                 "starvation resistance", "wing length", "growth rate", "litter size", "lifespan"]


def trait(o):
    return R.choice(TRAITS_PLANT if ORG.get(o, {}).get("group") == "plant" else TRAITS_ANIMAL)


TRAITS_B = ["disease status", "resistance vs susceptibility", "hip dysplasia", "mastitis", "type 2 diabetes",
            "lodging", "infection status", "survival", "presence of the phenotype", "early flowering yes/no"]


# ---------------------------------------------------------------- families
# each returns (S, goal, source, primary)
def f_cross():
    o = R.choice(CROPS); s = S()
    scheme = ch("RIL", "F2", "NAM", "MAGIC", "backcross", "biparental", "F1")
    sim_open(s); s(ch("a", "an", ""));
    if maybe(0.5):
        pop_phrase(s, o, 80, 1000)
    s(ch(f"{scheme} population", f"{scheme} lines", f"{scheme}s", f"{scheme} families", f"{scheme} cross"), "PRIMARY")
    if maybe(0.7):
        s(ch("from", "with", "derived from"))(fmt(R.choice([2, 2, 2, 4, 5, 8, 16, 26])), "FOUNDERS")(ch("founders", "parents", "parental lines", "inbred parents"))
    if maybe(0.6):
        s(ch(",", "and", "with"))(fmt(R.randint(2, 10)), "SELFGEN")(ch("generations of selfing", "selfing generations", "rounds of selfing"))
    if maybe(0.25):
        s(ch(",", "after"))(fmt(R.randint(1, 5)), "BCGEN")(ch("backcrosses", "backcross generations", "rounds of backcrossing"))
    org_phrase(s, o, 0.6); distractor(s); ask_end(s)
    return s, "simulate", "synthetic", "cross_design"


def f_cross_animal():
    o = R.choice(LAB_ANIMALS); s = S()
    sim_open(s)
    s(ch("a", "an"))(ch("F2 intercross", "recombinant inbred panel", "advanced intercross", "backcross"), "PRIMARY")
    s(ch("of", "with"))(fmt(R.randint(100, 800)), "POP")(ch(*org_units(o)))
    if maybe(0.6):
        s("from")(fmt(R.choice([2, 2, 8])), "FOUNDERS")(ch("founder strains", "founders", "parental strains"))
    org_phrase(s, o, 0.7); ask_end(s)
    return s, "simulate", "synthetic", "cross_design"


def f_pedigree():
    s = S()
    kind = R.random()
    if kind < 0.3:
        o = R.choice(SELFERS); scheme = ch("selfing", "self-pollination", "self-fertilization")
    elif kind < 0.7:
        o = R.choice(LAB_ANIMALS + LIVESTOCK); scheme = ch("sib mating", "full-sib mating", "brother-sister mating", "sibling mating", "full sib mating")
    else:
        o = "human"; scheme = ch("first-cousin marriages", "first cousin mating", "second-cousin marriages", "consanguineous marriages")
    sim_open(s)
    g = R.randint(3, 30)
    if maybe(0.5):
        s(fmt(g), "GENS")(ch("generations of", "gens of"))(scheme)
    else:
        s(scheme)(ch("for", "over"))(fmt(g), "GENS")(ch("generations", "gens"))
    if maybe(0.6):
        s(ch("in a population of", "with", "keeping"))(fmt(R.randint(20, 3000)), "POP")(ch(*org_units(o)) + ch("", " per generation", ""))
    org_phrase(s, o, 0.8)
    if maybe(0.7):
        s(ch("and track", "then estimate", "and show", "and plot"))(ch("inbreeding coefficients", "pedigree F", "runs of homozygosity", "heterozygosity loss", "the inbreeding coefficient"), "PRIMARY")
    ask_end(s)
    return s, "simulate", "synthetic", "pedigree_inbreeding"


def f_pedigree_provided():
    o = R.choice(LIVESTOCK + ["human", "mouse"]); s = S()
    prov_open(s, ch("a pedigree of", "a pedigree file with", "pedigree records for"))
    s(fmt(R.randint(100, 5000)), "POP")(ch(*org_units(o)))
    if maybe(0.5):
        s(ch("spanning", "over"))(fmt(R.randint(4, 20)), "GENS")("generations")
    org_phrase(s, o, 0.4)
    s(ch(". Estimate", ". Compute", ". What are the", ". Can you calculate"))(ch("inbreeding coefficients", "pedigree inbreeding", "kinship and inbreeding"), "PRIMARY")
    ask_end(s)
    return s, "analyze", "provided", "pedigree_inbreeding"


def f_pheno():
    o = R.choice(CROPS + LIVESTOCK); s = S()
    sim_open(s)
    fam = ch("continuous", "binary", "ordinal", "count", "negative binomial", "liability-scale")
    s(ch("a", ""))(fam)(ch("trait", "phenotype", "phenotypes"), "PRIMARY")
    if maybe(0.8):
        s(ch("with heritability", "with h2 =", "with h2", "h2"))(ch("0.2", "0.3", "0.35", "0.4", "0.5", "0.6", "0.75", "0.8"), "H2")
    if maybe(0.6):
        s(ch("and", ","))(fmt(R.randint(1, 100)), "NCAUSAL")(ch("causal variants", "QTL", "QTNs", "causal loci", "additive QTL"))
    if maybe(0.5):
        s("for")(fmt(R.randint(100, 1000)), "POP")(ch(*org_units(o)))
    if maybe(0.4):
        s(ch("on my genotype file", "using my VCF", "from our genotypes"))
        src = "provided"
    else:
        src = "synthetic"
    org_phrase(s, o, 0.5); ask_end(s)
    return s, "simulate", src, "phenotype_outcome_model"


def f_gwas():
    o = R.choice(CROPS + LIVESTOCK + LAB_ANIMALS + ["human"]); s = S(); binary = maybe(0.35)
    m = ch("GWAS", "genome-wide association study", "association mapping", "association study", "mixed-model GWAS",
           "genome-wide association scan")
    if binary and maybe(0.5) and o != "human":
        tr = R.choice(TRAITS_B)
    else:
        tr = R.choice(TRAITS_B) if binary else trait(o)
    mode = R.random()
    if mode < 0.55:      # analyze provided data
        prov_open(s, ch("genotypes and phenotypes for", "a VCF and", "SNP data on"))
        if binary and maybe(0.6):
            s(fmt(R.randint(100, 8000)), "CASES")(ch("cases", "affected animals" if o != "human" else "cases"))(ch("and", "plus"))(fmt(R.randint(100, 9000)), "CONTROLS")("controls")
        else:
            if o in ORG_PANEL and maybe(0.3):
                s(fmt(R.randint(80, 3000)), "POP")(R.choice(ORG_PANEL[o]), "ORG")(ch(*org_units(o)))
            else:
                pop_phrase(s, o, 80, 20000)
        s(ch("with", "measured for", "scored for"))(tr)
        s(ch(". Run a", ". Can you do a", ". Please run a", ". I want a", ". Do"))(m, "PRIMARY")
        if maybe(0.5):
            s(ch("with", "using", "including"))(ch("kinship", "PCs", "a kinship matrix", "principal components", "FarmCPU", "an LMM", "Firth correction" if binary else "Bonferroni"))
        goal, src = "analyze", "provided"
    else:                # simulate then analyze
        sim_open(s)
        s(ch("a cohort of", "a population of", ""))(fmt(R.randint(200, 5000)), "POP")(ch(*org_units(o)))
        if binary:
            s(ch("with a binary trait", "with a disease trait", "with a case-control phenotype"))
        if maybe(0.7):
            s(ch("with heritability", "h2"))(ch("0.3", "0.4", "0.5", "0.6", "0.25"), "H2")
        if maybe(0.6):
            s(ch("and", ","))(fmt(R.randint(1, 60)), "NCAUSAL")(ch("causal SNPs", "QTL", "causal variants"))
        s(ch(", then run a", "and run a", ", then do", "and test with a"))(m, "PRIMARY")
        goal, src = "simulate_and_analyze", "synthetic"
    org_phrase(s, o, 0.6); distractor(s); ask_end(s)
    return s, goal, src, "binary_gwas" if binary else "quantitative_gwas"


def f_gwas_panel():
    o = R.choice([x for x in ORG_PANEL if x in CROPS + LAB_ANIMALS + LIVESTOCK] or ["maize"]); s = S()
    s(ch("We phenotyped", "We scored", "I measured"))(fmt(R.randint(60, 2000)), "POP")(R.choice(ORG_PANEL[o]), "ORG")(ch(*org_units(o)))
    s(ch("for", "on"))(trait(o))(ch(". Can we", ". How do we", ". Please"))(ch("map it with a", "run a", "do a"))(ch("GWAS", "genome-wide association study", "association scan"), "PRIMARY")
    ask_end(s)
    return s, "analyze", "provided", "quantitative_gwas"


def f_qtl():
    o = R.choice(CROPS + LAB_ANIMALS); s = S()
    if maybe(0.5):
        prov_open(s, ""); s(fmt(R.randint(80, 600)), "POP")(ch("RILs", "F2 plants", "DH lines", "backcross progeny", "F2 mice"))(ch("genotyped and phenotyped for", "with genotypes and"))(trait(o))
        s(ch(". Find the", ". Map the", ". Run"))(ch("QTL", "QTL mapping", "interval mapping", "composite interval mapping"), "PRIMARY"); goal, src = "analyze", "provided"
    else:
        sim_open(s); s(fmt(R.randint(80, 600)), "POP")(ch("RILs", "F2 plants", "DH lines"))
        if maybe(0.5):
            s("from")(fmt(2), "FOUNDERS")(ch("parents", "founders"))
        s(ch("and then", ", then", "and"))(ch("map QTL", "run QTL mapping", "do interval mapping"), "PRIMARY")
        s(ch("", "for"))(trait(o)); goal, src = "simulate_and_analyze", "synthetic"
    org_phrase(s, o, 0.6); ask_end(s)
    return s, goal, src, "linkage_qtl"


def f_gs():
    o = R.choice(CROPS + LIVESTOCK); s = S()
    prov_open(s, ""); pop_phrase(s, o, 200, 5000); s(ch("with genotypes and", "genotyped and phenotyped for"))(trait(o))
    s(ch(". Fit", ". Run", ". Train", ". Use"))(ch("GBLUP", "rrBLUP", "genomic prediction", "genomic selection", "BayesB"), "PRIMARY")
    if maybe(0.5):
        s(ch("with 5-fold cross-validation", "and predict the next cycle", "with leave-one-family-out validation", ""))
    org_phrase(s, o, 0.5); ask_end(s)
    return s, "analyze", "provided", "genomic_prediction"


def f_multitrait():
    o = R.choice(CROPS + LIVESTOCK + LAB_ANIMALS); s = S()
    prov_open(s, ""); pop_phrase(s, o, 100, 3000)
    s(ch("measured for", "phenotyped for"))(trait(o))("and")(trait(o))
    s(ch(". Run a", ". Can we do a", ". I want a"))(ch("multi-trait GWAS", "multivariate GWAS", "longitudinal GWAS", "repeated measures GWAS"), "PRIMARY")
    ask_end(s)
    return s, "analyze", "provided", "multitrait_longitudinal_gwas"


def f_rare():
    s = S(); prov_open(s, ""); s(fmt(R.randint(500, 50000)), "POP")(ch("exomes", "whole genomes", "sequenced people"))
    s(ch(". Run", ". Do", ". Please run"))(ch("SKAT", "a burden test", "gene-based tests", "SKAT-O"), "PRIMARY")(ch("for", "on"))(trait("human"))
    ask_end(s)
    return s, "analyze", "provided", "rare_variant_and_gene_based_gwas"


def f_relcal():
    o = R.choice(["human", "mouse", "dog", "cattle"]); s = S(); sim_open(s)
    s(ch("related", "closely related", "family-structured"))(ch(*org_units(o)))
    if maybe(0.6):
        s("(")(fmt(R.randint(200, 3000)), "POP")(ch("in total)", "individuals)"))
    s(ch("and compare", "and check"))(ch("GWAS inflation", "genomic inflation", "lambda GC"), "PRIMARY")(ch("with and without a kinship matrix", "under PCs vs kinship correction", ""))
    ask_end(s)
    return s, "simulate_and_analyze", "synthetic", "population_relatedness_gwas_calibration"


def f_bulk_sim():
    o = R.choice(CROPS + LAB_ANIMALS + ["human"]); s = S()
    sim_open(s); s(ch("bulk RNA-seq counts", "RNA-seq counts", "RNA-seq data", "expression counts", "bulk RNA-seq", "a count matrix"), "PRIMARY")
    s(ch("for", "with"))
    if maybe(0.8):
        conds(s, ch(2, 2, 3))
    else:
        s(fmt(R.randint(2, 6)), "NCOND")(ch("conditions", "treatments", "groups"))
    if maybe(0.8):
        s(ch(",", "with", "and"))(fmt(R.randint(2, 12)), "REPS")(ch("biological replicates", "replicates", "reps", "biological reps"))(ch("each", "per condition", "per group", ""))
    if maybe(0.25):
        s("and")(fmt(R.randint(2, 4)), "TECHREPS")("technical replicates")
    if maybe(0.3):
        s(ch("across", "in"))(fmt(R.randint(2, 4)), "BATCHES")("batches")
    if maybe(0.3):
        s(ch("using a", "with a"))(ch("negative binomial", "zero-inflated", "ZINB"))("model")
    org_phrase(s, o, 0.6); distractor(s); ask_end(s)
    return s, "simulate", "synthetic", "bulk_counts_simulation"


def f_de():
    o = R.choice(CROPS + LAB_ANIMALS + ["human"]); s = S()
    prov_open(s, ch("a count matrix", "RNA-seq counts", "a counts table", "featureCounts output"))
    org_phrase(s, o, 0.5)
    if maybe(0.5):
        s("with")(fmt(R.randint(2, 8)), "REPS")(ch("replicates per group", "biological replicates each"))
    s(ch(". Find", ". Run", ". I need", ". Which genes are"))(ch("DEGs", "differential expression", "DESeq2", "edgeR", "limma-voom", "differentially expressed genes", "differentially expressed"), "PRIMARY")
    s(ch("between", "for", "comparing"))
    conds(s)
    ask_end(s)
    return s, "analyze", "provided", "differential_expression"


def f_timeseries():
    o = R.choice(CROPS + LAB_ANIMALS); s = S(); sim = maybe(0.4)
    if sim:
        sim_open(s); s(ch("a time course", "a time-series RNA-seq experiment", "time-course expression"), "PRIMARY")
    else:
        prov_open(s, ch("RNA-seq from", "bulk RNA-seq across", "expression data at"))
    if maybe(0.6):
        s(fmt(R.randint(3, 10)), "TIMEPTS")(ch("time points", "timepoints", "stages", "sampling times"))
    else:
        s(ch("at", "sampled at"))("0, 2, 6, 12 and 24 hours")
    if maybe(0.5):
        s(",")(fmt(R.randint(2, 6)), "REPS")(ch("replicates each", "reps per time point"))
    if not sim:
        s(ch(". Which genes change", ". Model expression", ". Find genes that respond"))(ch("over time", "across the time course", "with a spline"), "PRIMARY")
    org_phrase(s, o, 0.6); ask_end(s)
    return s, ("simulate" if sim else "analyze"), ("synthetic" if sim else "provided"), "time_series_rnaseq"


def f_splicing():
    o = R.choice(CROPS + LAB_ANIMALS + ["human"]); s = S()
    prov_open(s, ch("RNA-seq BAMs and a GTF", "aligned reads", "BAM files")); org_phrase(s, o, 0.5)
    s(ch(". Test", ". Look for", ". Run"))(ch("differential exon usage", "alternative splicing", "isoform switching", "LeafCutter", "DEXSeq", "differential splicing"), "PRIMARY")
    s("between"); conds(s); ask_end(s)
    return s, "analyze", "provided", "alternative_splicing"


def f_ase():
    o = R.choice(["maize", "mouse", "human", "rice", "drosophila"]); s = S()
    prov_open(s, ch("RNA-seq and phased SNPs for", "F1 hybrid RNA-seq from", "RNA-seq of"))
    s(fmt(R.randint(3, 60)), "POP")(ch("F1 hybrids", "individuals", "samples"))
    s(ch(". Measure", ". Estimate", ". Test"))(ch("allele-specific expression", "allelic imbalance", "ASE"), "PRIMARY"); org_phrase(s, o, 0.6); ask_end(s)
    return s, "analyze", "provided", "allele_specific_expression"


def f_power():
    o = R.choice(LAB_ANIMALS + CROPS + ["human"]); s = S()
    s(ch("How many", "how many"))(ch(*org_units(o), "replicates", "samples"))(ch("do I need", "should we use", "per group do we need"))
    s(ch("for", "in an"))(ch("RNA-seq", "an RNA-seq experiment", "bulk RNA-seq"), "PRIMARY")(ch("to detect a 2-fold change", "to find 1.5-fold changes", "", "at 80% power"))
    org_phrase(s, o, 0.5); ask_end(s)
    return s, "plan", "ask", "rnaseq_power"


def f_eqtl():
    o = R.choice(CROPS + ["human", "mouse", "rat", "cattle"]); s = S()
    prov_open(s, ch("genotypes and expression for", "a VCF and RNA-seq from", "matched genotypes and RNA-seq for"))
    s(fmt(R.randint(40, 1000)), "DONORS")(ch("donors", "individuals", "samples", "lines", "accessions"))
    s(ch(". Map", ". Run", ". Find"))(ch("cis-eQTLs", "eQTLs", "cis and trans eQTLs", "trans-eQTLs", "expression QTLs"), "PRIMARY")
    org_phrase(s, o, 0.6); ask_end(s)
    return s, "analyze", "provided", "bulk_eqtl"


def f_sqtl():
    o = R.choice(CROPS + ["human"]); s = S()
    prov_open(s, ch("junction counts and genotypes for", "RNA-seq and SNPs for"))
    s(fmt(R.randint(50, 800)), "DONORS")(ch("accessions", "donors", "individuals"))
    s(ch(". Map", ". Run", ". Find"))(ch("splicing QTLs", "sQTLs", "splice QTL"), "PRIMARY"); org_phrase(s, o, 0.6); ask_end(s)
    return s, "analyze", "provided", "splicing_qtl"


def celltypes(s, k):
    cts = R.sample(CELLTYPES, k)
    for i, c in enumerate(cts):
        if i:
            s("and" if i == k - 1 else ",")
        s(c, "CELLTYPE")


def f_sceqtl():
    s = S(); sim = maybe(0.5)
    if sim:
        sim_open(s); s(ch("single-cell data", "scRNA-seq"))(ch("for", "from"))(fmt(R.randint(20, 400)), "DONORS")("donors")
        s(ch("with", ","))(fmt(R.randint(50, 2000)), "CELLS")(ch("cells per donor", "cells each"))
        s(ch("and map", ", then run", "and test"))(ch("cell-type eQTLs", "single-cell eQTLs", "pseudobulk eQTLs", "cell-type-specific eQTLs"), "PRIMARY")
        goal, src, prim = "simulate_and_analyze", "synthetic", "scrnaseq_simulation"
    else:
        prov_open(s, ch("pseudobulk matrices and genotypes for", "scRNA-seq and genotypes from"))
        s(fmt(R.randint(20, 400)), "DONORS")("donors")
        s(ch(". Map", ". Run"))(ch("cell-type eQTLs", "single-cell eQTLs", "pseudobulk eQTLs"), "PRIMARY")
        if maybe(0.5):
            s(ch("in", "for"))
            celltypes(s, R.randint(2, 3))
        goal, src, prim = "analyze", "provided", "single_cell_eqtl"
    ask_end(s)
    return s, goal, src, prim


def f_scsim():
    o = R.choice(["human", "mouse", "arabidopsis", "zebrafish", "maize"]); s = S()
    sim_open(s); s(ch("scRNA-seq", "single-cell RNA-seq", "single-cell data", "a single-cell dataset"), "PRIMARY")
    s(ch("with", "for"))
    if maybe(0.5):
        s(fmt(R.randint(2, 12)), "NCELLTYPES")("cell types")
    else:
        celltypes(s, R.randint(2, 4))
    if maybe(0.7):
        s(",")(fmt(R.randint(2, 40)), "DONORS")(ch("donors", "samples", "individuals"))
    if maybe(0.5):
        s(ch(",", "and"))(fmt(R.randint(50, 2000)), "CELLS")(ch("cells per donor", "cells each"))
    if maybe(0.4):
        s("and"); conds(s)
    org_phrase(s, o, 0.6); ask_end(s)
    return s, "simulate", "synthetic", "scrnaseq_simulation"


def f_scan():
    o = R.choice(["human", "mouse", "arabidopsis", "drosophila", "zebrafish"]); s = S()
    prov_open(s, ch("10x data from", "a Seurat object with", "scRNA-seq from"))
    s(fmt(R.randint(2, 40)), "DONORS")(ch("donors", "samples", "animals"))
    v, obj = R.choice([(". Cluster", "the cells"), (". Annotate", "the cell types"), (". Run", "pseudobulk DE"),
                       (". Find", "marker genes"), (". Do", "clustering and annotation"), (". Integrate and cluster", "the samples")])
    s(v)(obj, "PRIMARY")
    org_phrase(s, o, 0.6); ask_end(s)
    return s, "analyze", "provided", "scrnaseq_analysis"


def f_chip():
    o = R.choice(CROPS + LAB_ANIMALS + ["human"]); s = S(); atac = maybe(0.5); sim = maybe(0.5)
    assay = ch("ATAC-seq", "ATAC", "chromatin accessibility") if atac else ch("ChIP-seq", "TF ChIP-seq", "H3K27ac ChIP-seq", "histone ChIP")
    if sim:
        sim_open(s); s(assay, "PRIMARY")
    else:
        prov_open(s, ""); s(assay, "PRIMARY")(ch("data", "BAMs", "reads"))
    s(ch("for", "with"))
    conds(s)
    if maybe(0.7):
        s(",")(fmt(R.randint(1, 5)), "REPS")(ch("replicates each", "biological replicates", "reps"))
    if not sim:
        s(ch(". Call peaks and test", ". Find"))(ch("differential binding", "differential accessibility", "differential peaks"))
    org_phrase(s, o, 0.6); ask_end(s)
    return s, ("simulate" if sim else "analyze"), ("synthetic" if sim else "provided"), "chipseq_atacseq"


def f_reads():
    s = S(); o = R.choice(MICROBES + CROPS)
    sim_open(s); s(ch("NovaSeq X", "MiSeq", "PacBio HiFi", "nanopore", "Revio HiFi", "NextSeq 2000", "Illumina"))
    s(ch("reads", "sequencing reads", "FASTQ files"), "PRIMARY")(ch("at", ""))(f"{R.choice([5, 10, 15, 20, 30, 40, 50, 60, 100])}x", "COV")(ch("coverage", ""))
    s(ch("from my reference", "from a synthetic genome", "from this FASTA", "from the genome"))
    org_phrase(s, o, 0.4); ask_end(s)
    return s, "simulate", "synthetic", "read_simulation"


def f_assembly():
    s = S(); o = R.choice(MICROBES + CROPS)
    prov_open(s, ch("HiFi reads", "Illumina reads", "nanopore reads", "Illumina and HiFi reads"))
    org_phrase(s, o, 0.6)
    s(ch(". Assemble", ". Do a", ". Run a"))(ch("the genome", "de novo assembly", "hybrid assembly", "genome assembly"), "PRIMARY")(ch("and check it with QUAST", "", "and report N50"))
    ask_end(s)
    return s, "analyze", "provided", "genome_assembly"


def f_variants():
    o = R.choice(CROPS + LIVESTOCK + ["human"]); s = S()
    prov_open(s, ch("FASTQs for", "BAMs for", "WGS data for"))
    s(fmt(R.randint(3, 1000)), "POP")(ch(*org_units(o)))
    s(ch(". Call", ". Joint-call", ". Find"))(ch("SNPs and indels", "variants", "structural variants", "SNPs"), "PRIMARY")
    org_phrase(s, o, 0.6); ask_end(s)
    return s, "analyze", "provided", "variant_calling"


def f_annot():
    o = R.choice(CROPS + ["human", "mouse"]); s = S()
    prov_open(s, ch("a list of", "my"))
    s(fmt(R.randint(20, 900)))(ch("candidate genes", "genes", "SNPs"))
    org_phrase(s, o, 0.7)
    s(ch(". Which", ". Map", ". Annotate"))(ch("GO terms are enriched", "GO terms", "the nearest genes", "variant effects", "gene coordinates"), "PRIMARY")
    ask_end(s)
    return s, "analyze", "provided", "gene_coordinate_annotation"


def f_multiomics():
    o = R.choice(CROPS); s = S(); sim = maybe(0.5)
    if sim:
        sim_open(s); s(ch("genotypes, phenotypes and RNA-seq", "matched genotype and expression data"))(ch("in", "across"))(fmt(R.randint(2, 6)), "ENVS")("environments")
        s(ch("and test which genes explain", "and model"))(ch("the GxE", "genotype by environment interaction", "multi-omics GxE"), "PRIMARY")
        goal, src = "simulate_and_analyze", "synthetic"
    else:
        prov_open(s, ch("genotypes, RNA-seq and phenotypes for", "multi-omics data on"))
        pop_phrase(s, o, 50, 600)
        s(ch(". Does expression", ". Can we", ". I want to"))(ch("mediate", "integrate", "predict"), "PRIMARY")(ch("the SNP effect on yield", "the layers", "the phenotype"))
        goal, src = "analyze", "provided"
    org_phrase(s, o, 0.6); ask_end(s)
    return s, goal, src, "gxe_multiomics"


def f_popgen():
    o = R.choice(CROPS + LIVESTOCK + ["human"]); s = S()
    a = R.choice(["population_structure", "linkage_disequilibrium", "selection_scan"])
    prov_open(s, ch("a VCF of", "SNP chip data for", "genotypes for"))
    s(fmt(R.randint(50, 3000)), "POP")(ch(*org_units(o)))
    m = {"population_structure": ch("PCA", "ADMIXTURE", "FST", "kinship", "population structure"),
         "linkage_disequilibrium": ch("LD decay", "linkage disequilibrium", "haplotype blocks"),
         "selection_scan": ch("a selection scan", "iHS", "XP-EHH", "selective sweeps")}[a]
    s(ch(". Run", ". Compute", ". Look for", ". Estimate"))(m, "PRIMARY"); org_phrase(s, o, 0.6); ask_end(s)
    return s, "analyze", "provided", a


def f_microbiome():
    o = R.choice(["mouse", "human", "chicken", "pig", "cattle"]); s = S()
    prov_open(s, ch("16S data from", "shotgun metagenomes from", "16S amplicon reads from"))
    s(fmt(R.randint(10, 400)), "POP")(ch("samples", "animals", "birds", "people"))
    s(ch(". Compare", ". Run", ". Which taxa"))(ch("alpha and beta diversity", "differential abundance", "taxonomic profiles", "differ"), "PRIMARY")
    org_phrase(s, o, 0.5); ask_end(s)
    return s, "analyze", "provided", "microbiome_analysis"


# combinations, inputs, exclusions ------------------------------------------
def f_combo():
    s = S(); kind = R.random()
    if kind < 0.25:   # simulate -> analyze; simulation step is primary
        sim_open(s); s(ch("RNA-seq counts", "bulk RNA-seq"), "PRIMARY")(ch("for", "from"))(fmt(R.randint(50, 500)), "DONORS")("donors")
        s(ch("and then map", ", then run", "and test for"))(ch("eQTLs", "cis-eQTLs", "differential expression"), "SECONDARY")
        return s, "simulate_and_analyze", "synthetic", "bulk_counts_simulation"
    if kind < 0.45:   # analyze then annotate
        prov_open(s, ch("RNA-seq counts", "a count matrix")); s(ch(". Run", ". Do"))(ch("DESeq2", "differential expression", "edgeR"), "PRIMARY")
        s(ch("and then", ", then"))(ch("GO enrichment", "GO term enrichment", "annotate the genes"), "SECONDARY")
        return s, "analyze", "provided", "differential_expression"
    if kind < 0.65:   # input, not a step
        prov_open(s, "")
        s(ch("GWAS hits", "GWAS results", "eQTL results", "DE genes", "QTL peaks"), "INPUT")
        o = R.choice(CROPS + ["human"]); org_phrase(s, o, 0.6)
        s(ch(". Which genes are within", ". Find genes near them, within", ". List genes within"))(f"{R.choice([5, 10, 20, 50, 100])} kb")
        s(ch(".", "?"))
        if maybe(0.5):
            s(ch("Then", "Also"))(ch("annotate them with GO terms", "check GO enrichment"), "PRIMARY")
        else:
            s(ch("Annotate", "Map"))(ch("the nearest genes", "gene coordinates"), "PRIMARY")
        return s, "analyze", "provided", "gene_coordinate_annotation"
    if kind < 0.8:    # exclusion
        sim_open(s); s(ch("RNA-seq counts", "expression counts"), "PRIMARY")(ch("for", "with"))
        conds(s)
        s(ch(", but no", "- don't run", ", without", ", skip"))(ch("eQTL mapping", "GWAS", "DE analysis", "eQTLs"), "NEG")
        s(ch(", I'll do the analysis myself", "", ", just the counts"))
        return s, "simulate", "synthetic", "bulk_counts_simulation"
    # cross + GWAS: GWAS on the simulated population is the main analysis
    o = R.choice(CROPS); sim_open(s); s(fmt(R.randint(200, 2000)), "POP")
    s(ch("NAM", "MAGIC", "RIL"), "SECONDARY")(ch("lines", "families"))
    s(ch("and run a", ", then", "and do a"))(ch("GWAS", "association mapping"), "PRIMARY")(ch("for", "on"))(trait(o))
    org_phrase(s, o, 0.7); ask_end(s)
    return s, "simulate_and_analyze", "synthetic", "quantitative_gwas"


def f_unsupported():
    t = R.choice(UNSUPPORTED); s = S(); o = R.choice(CROPS + LAB_ANIMALS + ["human"])
    if maybe(0.5):
        sim_open(s); s(mention(t))(ch("data", "reads", "profiles", ""))
        goal, src = "simulate", "synthetic"
    else:
        prov_open(s, ""); s(mention(t))(ch("data", "results", "samples"))(ch(". Analyze it", ". What changes between", ". Can you process it"))
        goal, src = "analyze", "provided"
    if maybe(0.4):
        s("for"); conds(s)
    org_phrase(s, o, 0.5); ask_end(s)
    return s, goal, src, "unsupported"


def f_vague():
    s = S(); kind = R.random()
    if kind < 0.4:
        s(ch("RNA-seq", "rna seq", "single cell", "scRNA-seq", "ATAC-seq?", "help with my sequencing data", "bulk RNA-seq please"))
        ask_end(s); return s, "none", "ask", "none"
    if kind < 0.7:
        s(ch("run gwas", "do a GWAS", "gwas please", "association study"), "PRIMARY"); ask_end(s)
        return s, "analyze", "ask", "quantitative_gwas"
    s(ch("make single cell data", "simulate some data", "simulate a population", "generate a dataset"))
    ask_end(s)
    return s, "simulate", "synthetic", "none"


def f_demo():
    s = S(); a, m = R.choice([("quantitative_gwas", "GWAS"), ("cross_design", "NAM population"), ("chipseq_atacseq", "ChIP-seq"),
                              ("scrnaseq_simulation", "single-cell simulation"), ("bulk_counts_simulation", "RNA-seq counts")])
    s(ch("Show me a", "Run the", "Give me a", "Can I see a"))(m, "PRIMARY")(ch("example using the package demo data", "demo", "with the built-in example data"))
    ask_end(s)
    return s, ("analyze" if a in ("quantitative_gwas", "chipseq_atacseq") else "simulate"), "package_demo", a


FAMILIES = {
    "cross": (f_cross, 6), "cross_animal": (f_cross_animal, 2), "pedigree": (f_pedigree, 6), "pedigree_prov": (f_pedigree_provided, 2),
    "pheno": (f_pheno, 3), "gwas": (f_gwas, 8), "gwas_panel": (f_gwas_panel, 2), "qtl": (f_qtl, 3), "gs": (f_gs, 3),
    "multitrait": (f_multitrait, 2), "rare": (f_rare, 1), "relcal": (f_relcal, 1), "bulk_sim": (f_bulk_sim, 5), "de": (f_de, 4),
    "timeseries": (f_timeseries, 3), "splicing": (f_splicing, 2), "ase": (f_ase, 1), "power": (f_power, 2), "eqtl": (f_eqtl, 3),
    "sqtl": (f_sqtl, 1), "sceqtl": (f_sceqtl, 3), "scsim": (f_scsim, 4), "scan": (f_scan, 3), "chip": (f_chip, 4),
    "reads": (f_reads, 3), "assembly": (f_assembly, 2), "variants": (f_variants, 2), "annot": (f_annot, 2),
    "multiomics": (f_multiomics, 2), "popgen": (f_popgen, 3), "microbiome": (f_microbiome, 2), "combo": (f_combo, 7),
    "unsupported": (f_unsupported, 4), "vague": (f_vague, 2), "demo": (f_demo, 1),
}


def noise(text):
    if maybe(0.12):
        text = text.lower()
    if maybe(0.05):
        text = text.rstrip(".?!")
    return text


def eval_prompts():
    path = os.path.join(ROOT, "inst", "agent", "prompt_eval", "study_spec_examples_draft.jsonl")
    out = []
    if os.path.exists(path):
        for l in open(path):
            if l.strip():
                out.append(set(w.lower() for w in tokenize(json.loads(l)["prompt"]) if w.isalnum()))
    return out


def main():
    global R
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20000); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--weights", help="JSON {family: multiplier} from mine_errors.py")
    ap.add_argument("--max-overlap", type=float, default=0.6, help="drop sentences this similar (Jaccard) to any eval prompt")
    a = ap.parse_args()
    R = random.Random(a.seed)
    w = json.load(open(a.weights)) if a.weights else {}
    names = list(FAMILIES)
    weights = [FAMILIES[n][1] * w.get(n, 1.0) for n in names]
    guard = eval_prompts()
    rows, dropped, bad, seen = [], 0, 0, set()
    while len(rows) < a.n:
        fam = R.choices(names, weights)[0]
        s, goal, src, prim = FAMILIES[fam][0]()
        text, toks, tags = s.build()
        if tokenize(text) != toks:
            bad += 1
            continue
        text = noise(text)
        new = tokenize(text)
        if [x.lower() for x in new] != [x.lower() for x in toks[:len(new)]] or any(g != "O" for g in tags[len(new):]):
            bad += 1
            continue
        toks, tags = new, tags[:len(new)]
        if text in seen:
            continue
        ws = set(x.lower() for x in toks if x.isalnum())
        if any(len(ws & g) / max(len(ws | g), 1) >= a.max_overlap for g in guard):
            dropped += 1
            continue
        seen.add(text)
        assert prim in ANALYSES and all(t in LABELS for t in tags), (fam, prim)
        rows.append({"text": text, "tokens": toks, "tags": tags, "goal": goal, "data_source": src, "primary": prim, "family": fam})
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    c = collections.Counter(r["family"] for r in rows)
    print(f"wrote {len(rows)} sentences to {a.out}; {dropped} dropped as too close to eval prompts; {bad} tokenisation mismatches")
    print("families:", dict(c.most_common()))


if __name__ == "__main__":
    main()
