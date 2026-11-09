#!/usr/bin/env python3
"""Build inst/agent/lexicon/organisms.tsv and organism_aliases.tsv.

Sources
  * NCBI Taxonomy, via the OBO Foundry NCBITaxon release (ncbitaxon.obo.gz from
    https://github.com/obophenotype/ncbitaxon/releases): scientific names, GenBank
    common names, layperson names, equivalent names, older Latin synonyms,
    misspellings, lineage (group and family).
  * data-raw/lexicon/organisms_curated.tsv: the organism list and biology facts
    (ploidy, selfing, mating system, generation time, common reference). Curated.
  * data-raw/lexicon/organism_extra_aliases.tsv: strains, panels, breeds, cultivars,
    genome builds and informal words NCBI does not carry. Curated.

    python3 data-raw/lexicon/build_organisms.py /path/to/ncbitaxon.obo.gz
"""
import csv, gzip, re, sys
from collections import defaultdict

OBO = sys.argv[1]
HERE = "data-raw/lexicon"
OUT = "inst/agent/lexicon"

KEEP_TYPES = {  # NCBI synonym type -> alias_type
    "genbank_common_name": "common_name", "OMO:0003003": "common_name", "equivalent_name": "scientific_synonym",
    "synonym": "scientific_synonym", "genbank_synonym": "scientific_synonym", "OMO:0003006": "misspelling",
    "anamorph": "scientific_synonym", "teleomorph": "scientific_synonym", "genbank_anamorph": "scientific_synonym",
    "OMO:0003012": "acronym", "genbank_acronym": "acronym",
}
# Last words that take a regular English plural (others are mass nouns or irregular: mouse, sheep, rice...).
COUNTABLE = {"fly", "mosquito", "bee", "rat", "cat", "cow", "goat", "horse", "pig", "dog", "boar", "chicken", "bantam",
             "frog", "platanna", "monkey", "macaque", "worm", "silkworm", "moth", "stickleback", "bean", "pea", "potato",
             "tomato", "pepper", "sunflower", "cucumber", "banana", "apple", "grape", "strawberry", "poplar",
             "cottonwood", "tree", "medic", "nematode", "hen", "bull", "peanut", "chickpea", "soybean", "medaka"}
# NCBI names reviewed and removed: ambiguous, offensive in other senses, or food words that are not organism cues.
DROP = {"rape", "chocolate", "tapioca", "cocoa", "red rice", "bee", "noble cane"}
WEAKEN = {"wild boar", "eubacteria", "brewer's yeast"}
GROUPS = [(33208, "animal"), (33090, "plant"), (4751, "fungus"), (2, "bacterium"), (2157, "archaeon"),
          (10239, "virus"), (2759, "protist")]

curated = list(csv.DictReader(open(f"{HERE}/organisms_curated.tsv"), delimiter="\t"))
extra = list(csv.DictReader(open(f"{HERE}/organism_extra_aliases.tsv"), delimiter="\t"))
targets = {int(r["taxid"]) for r in curated}

# pass 1: parent map --------------------------------------------------------------
parent, cur = {}, None
with gzip.open(OBO, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("id: NCBITaxon:"):
            cur = int(line[14:].strip()) if line[14:].strip().isdigit() else None
        elif cur is not None and line.startswith("is_a: NCBITaxon:"):
            p = line[16:].split(" ")[0].strip()
            if p.isdigit():
                parent[cur] = int(p)
need = set()
for t in targets:
    x = t
    while x in parent and x not in need:
        need.add(x)
        x = parent[x]
    need.add(x)

# pass 2: names, ranks, synonyms for targets and their ancestors ------------------
info, cur = {}, None
with gzip.open(OBO, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("id: NCBITaxon:"):
            tid = line[14:].strip()
            cur = info.setdefault(int(tid), {"name": "", "rank": "", "syn": []}) if tid.isdigit() and int(tid) in need else None
        elif cur is None:
            continue
        elif line.startswith("name: "):
            cur["name"] = line[6:].strip()
        elif line.startswith("property_value: has_rank NCBITaxon:"):
            cur["rank"] = line.split("NCBITaxon:")[1].strip()
        elif line.startswith("synonym: "):
            m = re.match(r'synonym: "((?:[^"\\]|\\.)*)" (\w+) ([\w:]+)', line)
            if m:
                cur["syn"].append((m.group(1), m.group(3)))


def lineage(t):
    out = [t]
    while out[-1] in parent:
        out.append(parent[out[-1]])
    return out


def plural(w):
    if re.search(r"(s|sh|ch|x|z)$", w):
        return w + "es"
    if re.search(r"[^aeiou]y$", w):
        return w[:-1] + "ies"
    if w.endswith("o"):
        return w + "es"
    return w + "s"


problems, organisms, aliases = [], [], []
for r in curated:
    t = int(r["taxid"])
    if t not in info or not info[t]["name"]:
        problems.append(f"taxid {t} ({r['organism']}) not found in NCBITaxon")
        continue
    name = info[t]["name"]
    if name.lower() != r["label"].lower():
        problems.append(f"{r['organism']}: curated label '{r['label']}' != NCBI name '{name}'")
    lin = lineage(t)
    group = next((g for gid, g in GROUPS if gid in lin), "")
    family = next((info[x]["name"] for x in lin if x in info and info[x]["rank"] == "family"), "")
    organisms.append({**r, "scientific_name": name, "group": group, "family": family})
    org = r["organism"]
    aliases.append((name, org, "scientific_name", "NCBITaxon", "strong", "no"))
    aliases.append((org.replace("_", " "), org, "canonical_id", "curated", "strong", "no"))
    common = []
    for s, typ in info[t]["syn"]:
        if typ in KEEP_TYPES:
            at = KEEP_TYPES[typ]
            aliases.append((s, org, at, "NCBITaxon", "weak" if at == "acronym" else "strong", "yes" if at == "acronym" else "no"))
            if at == "common_name":
                common.append(s)
    # generated: "D. melanogaster", "D.melanogaster" for binomials
    parts = name.split()
    if len(parts) >= 2 and parts[0][0].isupper() and parts[1].islower():
        aliases.append((f"{parts[0][0]}. {' '.join(parts[1:])}", org, "abbreviation", "generated", "strong", "no"))
        aliases.append((f"{parts[0][0]}.{' '.join(parts[1:])}", org, "abbreviation", "generated", "strong", "no"))
    for c in common:  # plural of countable common names only: "rats", "fruit flies"
        words = c.split()
        if words[-1].lower() in COUNTABLE:
            aliases.append((" ".join(words[:-1] + [plural(words[-1])]), org, "plural", "generated", "strong", "no"))
for e in extra:
    aliases.append((e["alias"], e["organism"], e["alias_type"], "curated", e["strength"], e["case_sensitive"]))

known = {o["organism"] for o in organisms}
for a in aliases:
    if a[1] not in known:
        problems.append(f"alias '{a[0]}' points to unknown organism '{a[1]}'")

# de-duplicate; resolve collisions (curated wins, otherwise drop the alias) ---------
def norm(a, cs):
    a = re.sub(r"\s+", " ", a.strip())
    return a if cs == "yes" else a.lower()

by = defaultdict(list)
aliases = [a for a in aliases if a[0].lower() not in DROP]
aliases = [a[:4] + ("weak",) + a[5:] if a[0].lower() in WEAKEN else a for a in aliases]
for a in aliases:
    if len(a[0].strip()) < 3 and a[5] != "yes":
        continue
    by[norm(a[0], a[5])].append(a)
final, dropped = [], []
rank_src = {"curated": 0, "NCBITaxon": 1, "generated": 2}
for k, lst in by.items():
    orgs = {a[1] for a in lst}
    cur_orgs = {a[1] for a in lst if a[3] == "curated"}
    if len(orgs) > 1:
        if len(cur_orgs) == 1:
            lst = [a for a in lst if a[1] in cur_orgs]
        else:
            dropped.append(f"{k}: {sorted(orgs)}")
            continue
    best = sorted(lst, key=lambda a: (rank_src[a[3]], a[4] != "strong"))[0]
    final.append(best)
final.sort(key=lambda a: (a[1], a[2], a[0].lower()))

with open(f"{OUT}/organisms.tsv", "w", newline="") as f:
    cols = ["organism", "taxid", "scientific_name", "group", "family", "ploidy", "ploidy_note", "can_self",
            "mating_system", "generation_time", "common_reference", "simitall_knowledge"]
    w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    w.writerows(organisms)
with open(f"{OUT}/organism_aliases.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["alias", "organism", "alias_type", "source", "strength", "case_sensitive"])
    w.writerows(final)
print(f"{len(organisms)} organisms, {len(final)} aliases ({sum(a[3] == 'NCBITaxon' for a in final)} from NCBI)")
print("dropped (ambiguous):", dropped)
print("PROBLEMS:" if problems else "no problems", *problems, sep="\n  ")
