#!/usr/bin/env python3
"""Write a human review sheet for one batch and refresh coverage_map.json counts.

    python3 data-raw/agent_eval/review_sheet.py 2
If inst/agent/prompt_eval/results/all_by_example.tsv exists (from score_examples.R),
the sheet also shows what the current local agent picked, so reviewers can see bugs.
"""
import csv, json, os, sys
from collections import OrderedDict

batch = int(sys.argv[1]) if len(sys.argv) > 1 else 2
base = "inst/agent/prompt_eval"
rows = [json.loads(l) for l in open(f"{base}/study_spec_examples_draft.jsonl") if l.strip()]

# coverage counts --------------------------------------------------------------
cov_path = f"{base}/coverage_map.json"
cov = json.load(open(cov_path))
counts = {}
for r in rows:
    counts[r["workflow"]] = counts.get(r["workflow"], 0) + 1
known = {w["workflow"] for w in cov["workflows"]}
for w in cov["workflows"]:
    w["have"] = counts.get(w["workflow"], 0)
for w in sorted(set(counts) - known):
    cov["workflows"].append({"workflow": w, "family": "", "catalog_status": "", "target": 0, "have": counts[w], "tutorial_source": ""})
cov["totals"]["have"] = len(rows)
json.dump(cov, open(cov_path, "w"), indent=2)
open(cov_path, "a").write("\n")

# current agent output -----------------------------------------------------------
got = {}
tsv = f"{base}/results/all_by_example.tsv"
if os.path.exists(tsv):
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        flag = []
        if r["primary_ok"] == "FALSE":
            flag.append("wrong analysis")
        if r["must_not_violations"] not in ("0", ""):
            flag.append("ran a forbidden analysis")
        if r["stop_ok"] == "FALSE":
            flag.append("did not stop")
        if r["invented_fields"]:
            flag.append("invented " + r["invented_fields"])
        got[r["id"]] = (r["got_primary"] if r["got_primary"] != "NA" else "nothing") + (" (" + "; ".join(flag) + ")" if flag else "")


def md(x):
    return str(x).replace("|", "\\|").replace("\n", " ")


keys = ["goal", "organism", "data_source", "region", "primary_analysis", "secondary_analyses", "design", "parameters",
        "execution_status", "do_not_run"]
mine = [r for r in rows if r["batch"] == batch]
groups = OrderedDict()
for r in mine:
    groups.setdefault(r["workflow"], []).append(r)

out = [f"# Batch {batch} review sheet (draft)", "",
       f"{len(mine)} examples. For each row: is the prompt realistic, and is the expected spec what SimItAll *should* produce?",
       "Tell Claude which ids to fix, or set `review_status` to `reviewed` in `study_spec_examples_draft.jsonl`.",
       "The **SimItAll today** column is what the current local agent picked (from `score_examples.R`); it shows bugs, not the answer key.", ""]
out.append("## Contents\n")
for w, rs in groups.items():
    out.append(f"- [{w}](#{w.replace('_', '-')}) ({len(rs)})")
out.append("")
for w, rs in groups.items():
    out += [f"## {w}", "", "| id | split | prompt | expected (key fields) | must ask | must not | SimItAll today | note |",
            "|---|---|---|---|---|---|---|---|"]
    for r in rs:
        s = {k: r["expected_spec"][k] for k in keys if k in r["expected_spec"] and r["expected_spec"][k] not in ([], {}, None)}
        out.append("| `{}` | {} | {} | `{}` | {} | {} | {} | {} |".format(
            r["id"], r["split"], md(r["prompt"]), md(json.dumps(s, ensure_ascii=False)),
            md(", ".join(r["must_ask"])), md(", ".join(r["must_not"])), md(got.get(r["id"], "")), md(r.get("notes", ""))))
    out.append("")
open(f"{base}/review_batch{batch}.md", "w").write("\n".join(out))
print(f"wrote {base}/review_batch{batch}.md ({len(mine)} rows); coverage have = {len(rows)}")
