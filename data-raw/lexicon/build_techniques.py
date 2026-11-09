#!/usr/bin/env python3
"""Build inst/agent/lexicon/techniques.tsv and technique_aliases.tsv.

Sources
  * EFO (Experimental Factor Ontology, EBI; efo.obo from github.com/EBISPOT/efo releases)
  * OBI (Ontology for Biomedical Investigations; views/obi.obo from github.com/obi-ontology/obi)
  * EDAM (bioinformatics operations/topics; releases/EDAM.tsv from github.com/edamontology/edamontology)
  * ENA/SRA controlled vocabularies (library strategy, instrument model) from github.com/enasequence/schema
  * inst/agent/capabilities/analysis_catalog.json (subanalysis aliases, parameters, execution status)
  * data-raw/lexicon/techniques_curated.tsv: technique list, routing and extra aliases (curated).
    In extra_aliases, a leading "~" marks a weak alias (used only when nothing stronger matches)
    and "!" marks a case-sensitive one.

    python3 data-raw/lexicon/build_techniques.py EFO_OBO OBI_OBO EDAM_TSV
"""
import csv, json, re, sys
from collections import defaultdict

sys.path.insert(0, "data-raw/lexicon")
from obo import parse

efo_path, obi_path, edam_path = sys.argv[1:4]
OUT = "inst/agent/lexicon"
catalog = json.load(open("inst/agent/capabilities/analysis_catalog.json"))
subs = {s["id"]: s for f in catalog["families"] for s in f["subanalyses"]}

onto = {}
for path in (efo_path, obi_path):
    for k, v in parse(path).items():
        if v["name"] and not v["obsolete"]:
            onto[k.replace("efo:EFO_", "EFO:")] = v
edam = {}
with open(edam_path) as f:
    r = csv.reader(f, delimiter="\t")
    h = next(r)
    iC, iL, iS = h.index("Class ID"), h.index("Preferred Label"), h.index("Synonyms")
    for row in r:
        if len(row) > iS:
            edam["EDAM:" + row[iC].rsplit("/", 1)[-1]] = (row[iL], [s for s in row[iS].split("|") if s])

problems, techs, aliases = [], [], []


def add(alias, tech, source, strength="strong", cs="no"):
    alias = re.sub(r"\s+", " ", alias.strip())
    if alias.startswith("!"):
        alias, cs = alias[1:], "yes"
    if alias.startswith("~"):
        alias, strength = alias[1:], "weak"
    if len(alias) >= 2:
        aliases.append((alias, tech, source, strength, cs))


for r in csv.DictReader(open("data-raw/lexicon/techniques_curated.tsv"), delimiter="\t"):
    t = r["technique"]
    for a in [r["routes_to"]] + [c for c in r["candidates"].split(";") if c]:
        if a and a not in subs:
            problems.append(f"{t}: unknown catalog analysis {a}")
    if r["set_parameter"]:
        k, v = r["set_parameter"].split("=")
        params = {p["id"]: p for p in subs[r["routes_to"]].get("parameters", [])}
        if k not in params:
            problems.append(f"{t}: {r['routes_to']} has no parameter {k}")
        elif params[k].get("allowed_values") and v not in params[k]["allowed_values"]:
            problems.append(f"{t}: {k}={v} not in {params[k]['allowed_values']}")
    if r["routes_to"]:
        support = subs[r["routes_to"]].get("execution_status", "")
    elif r["candidates"]:
        support = "ask_which_analysis"
    else:
        support = "not_supported"
    techs.append({**r, "support": support})
    add(r["label"], t, "curated")
    for oid in [o for o in r["ontology_ids"].split(";") if o]:
        if oid.startswith("EDAM:"):
            if oid not in edam:
                problems.append(f"{t}: {oid} not in EDAM")
                continue
            add(edam[oid][0], t, oid, "weak")  # EDAM labels are broad ("Genetic mapping"); weak
            continue
        if oid not in onto:
            problems.append(f"{t}: {oid} not found in EFO/OBI")
            continue
        add(re.sub(r" assay$", "", onto[oid]["name"]), t, oid)
        for s, scope, _ in onto[oid]["syn"]:
            if scope == "EXACT" and len(s) <= 60:
                add(s, t, oid)
    for s in [x for x in r["sra_terms"].split(";") if x]:
        add(s, t, "ENA/SRA vocabulary")
    for s in [x for x in r["extra_aliases"].split("|") if x]:
        add(s, t, "curated")
    if r["routes_to"] and r["kind"] == "analysis":
        for s in subs[r["routes_to"]].get("aliases", []) or []:
            add(s, t, "analysis_catalog")


def norm(a, cs):
    a = re.sub(r"[-_/]", " ", a)
    a = re.sub(r"\s+", " ", a).strip()
    return a if cs == "yes" else a.lower()


by = defaultdict(list)
for a in aliases:
    by[norm(a[0], a[4])].append(a)
prio = lambda a: (0 if a[2] == "curated" else 1 if a[2] == "analysis_catalog" else 2, a[3] != "strong")
final, dropped = [], []
for k, lst in by.items():
    ts = {a[1] for a in lst}
    if len(ts) > 1:
        cur = {a[1] for a in lst if a[2] == "curated"}
        if len(cur) == 1:
            lst = [a for a in lst if a[1] in cur]
        else:
            dropped.append(f"{k}: {sorted(ts)}")
            continue
    final.append(sorted(lst, key=prio)[0])
final.sort(key=lambda a: (a[1], a[0].lower()))

cols = ["technique", "kind", "label", "routes_to", "candidates", "set_parameter", "support", "ontology_ids", "notes"]
with open(f"{OUT}/techniques.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
    w.writeheader()
    w.writerows(techs)
with open(f"{OUT}/technique_aliases.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["alias", "technique", "source", "strength", "case_sensitive"])
    w.writerows(final)
print(f"{len(techs)} techniques, {len(final)} aliases "
      f"({sum(not a[2] in ('curated', 'analysis_catalog') for a in final)} from EFO/OBI/EDAM/SRA)")
print("dropped (ambiguous):", dropped)
print("PROBLEMS:" if problems else "no problems", *problems, sep="\n  ")
