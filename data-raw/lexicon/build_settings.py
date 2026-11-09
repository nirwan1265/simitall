#!/usr/bin/env python3
"""Build inst/agent/lexicon/settings.tsv (table 3: methods and design terms -> settings).

Each row says: when this alias appears and the request routes to one of these
analyses, set this field to this value. Sources:
  * inst/agent/capabilities/analysis_catalog.json: every enum/set value of every
    typed parameter becomes an alias of itself ("Firth_logistic" -> "Firth logistic").
  * data-raw/lexicon/settings_curated.tsv: synonyms people actually write
    ("brother-sister mating" -> sib_mating, "EMMAX" -> linear_mixed_model).
Every analysis, field and value is checked against the catalog.

    python3 data-raw/lexicon/build_settings.py
"""
import csv, json, re

catalog = json.load(open("inst/agent/capabilities/analysis_catalog.json"))
subs = {s["id"]: s for f in catalog["families"] for s in f["subanalyses"]}
# StudySpec field for each catalog parameter (design keys follow prompt_eval/README.md)
DESIGN_KEYS = {("cross_design", "scheme"): "design.scheme", ("cross_design", "founder_count"): "design.founder_count",
               ("cross_design", "selfing_generations"): "design.selfing_generations"}
STOP = {"none", "gene", "cell", "linear", "cis", "trans", "continuous", "categorical", "cohort", "pathway", "binomial",
        "trio", "single_sample", "coordinates", "taxonomy", "diversity", "assembly", "survival", "ordinal", "liability",
        "gaussian", "bernoulli", "selfing", "backcross", "introgression", "bottleneck", "hybrid", "short_read",
        "long_read", "prediction", "mediation", "association", "multitrait", "longitudinal", "spline", "biparental"}


def field(analysis, pid):
    return DESIGN_KEYS.get((analysis, pid), f"parameters.{pid}")


rows, problems = [], []
for sid, s in subs.items():
    for p in s.get("parameters", []) or []:
        if p.get("type") == "array":  # lists of inputs (e.g. modalities), not a choice
            continue
        for v in p.get("allowed_values") or []:
            if v in STOP or len(v) < 3:
                continue
            strong = bool(re.search(r"[A-Z0-9_]", v))
            for alias in {v, v.replace("_", " ")}:
                rows.append({"alias": alias, "analysis": sid, "field": field(sid, p["id"]), "value": v,
                             "strength": "strong" if strong else "weak", "case_sensitive": "no", "source": "analysis_catalog"})

for r in csv.DictReader(open("data-raw/lexicon/settings_curated.tsv"), delimiter="\t"):
    for a in r["analysis"].split(";"):
        params = {field(a, p["id"]): p for p in (subs.get(a, {}).get("parameters") or [])}
        if a not in subs:
            problems.append(f"{r['alias']}: unknown analysis {a}")
        elif r["field"] not in params:
            problems.append(f"{r['alias']}: {a} has no {r['field']}")
        else:
            p = params[r["field"]]
            ok = (p.get("type") == "boolean" and r["value"] in ("true", "false")) or \
                 (p.get("allowed_values") and r["value"] in p["allowed_values"])
            if not ok:
                problems.append(f"{r['alias']}: {r['value']} not allowed for {a} {r['field']}")
        rows.append({**r, "analysis": a, "source": "curated"})

# curated rows override catalog rows for the same (alias, analysis, field)
seen, final = set(), []
for r in sorted(rows, key=lambda r: r["source"] != "curated"):
    k = (r["alias"].lower(), r["analysis"], r["field"])
    if k not in seen:
        seen.add(k)
        final.append(r)
final.sort(key=lambda r: (r["analysis"], r["field"], r["alias"].lower()))
with open("inst/agent/lexicon/settings.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, ["alias", "analysis", "field", "value", "strength", "case_sensitive", "source"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(final)
print(f"{len(final)} settings rows ({sum(r['source'] == 'analysis_catalog' for r in final)} from catalog enums) "
      f"covering {len({r['analysis'] for r in final})} analyses")
print("PROBLEMS:" if problems else "no problems", *problems, sep="\n  ")
