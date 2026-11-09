#!/usr/bin/env python3
"""Paraphrase generated prompts with a free local LLM (Ollama) while keeping every label exact.

Labelled spans are wrapped in [[double brackets]] and the LLM is asked to keep them.
A paraphrase is kept only if every labelled span can be found again, unchanged
(brackets may be dropped and case may change), and no number appears that was
not in the source sentence. Sentence labels (goal, data source, primary analysis)
are copied, because the meaning must not change.

    ollama pull llama3.1                     # or qwen2.5:7b (follows rules better), llama3.2:3b (faster)
    OLLAMA_NUM_PARALLEL=4 ollama serve       # in another terminal, lets 4 requests run at once
    python training/paraphrase_ollama.py --inp training/data/train.jsonl --out training/data/para.jsonl \
        --n 10000 --variants 2 --workers 4 --max-minutes 360

Resumable: every attempt (kept or rejected) is logged in <out>.tried, so a restart
continues where it stopped. --selftest checks the relabelling logic without Ollama.
"""
import argparse, concurrent.futures as cf, json, os, random, re, sys, time, urllib.request

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "inst", "python"))
from prompt_model import spans, tokenize  # noqa: E402

STYLES = [
    "as a busy PhD student typing quickly: casual, short, maybe lowercase",
    "as a careful postdoc writing to a statistician: complete sentences",
    "as a non-native English speaker: simple words, slightly unusual word order",
    "as a terse lab note: fragments, commas, no greeting",
    "as a PI asking for advice in a friendly tone",
    "as someone explaining the project in chat, with a short bit of context first",
    "as a question starting with 'Can you' or 'Is it possible to'",
]
PROMPT = """Rewrite the scientist's request {style}.
Rules:
- Copy every [[bracketed]] phrase exactly, brackets included, each exactly once. Do not translate, pluralise or split them.
- Do not add numbers, organisms, methods, analyses or facts. Do not drop any request.
- You may reorder clauses, add filler words, or turn it into a question.
- Reply with the rewritten request only.

Example
Request: Simulate [[200]] [[RIL lines]] from [[two]] founders with [[eight]] generations of selfing.
Rewritten: could you make [[RIL lines]] for me? [[200]] of them, [[two]] founders, [[eight]] generations of selfing

Request: {src}
Rewritten:"""


def bracketed(d):
    out, i = [], 0
    sp = spans(d["tokens"], d["tags"])
    for lab, text, start in sp:
        out += d["tokens"][i:start] + [f"[[{text}]]"]
        i = start + len(text.split(" "))
    return " ".join(out + d["tokens"][i:]), sp


def _norm(x):
    return re.sub(r"\s+", " ", " ".join(tokenize(x))).strip().lower()


def relabel(para, d):
    """Return (tokens, tags) for a paraphrase, or None if any label cannot be recovered exactly."""
    sp = spans(d["tokens"], d["tags"])
    plain = re.sub(r"\[\[|\]\]", "", para).strip()
    toks = tokenize(plain)
    low = [t.lower() for t in toks]
    tags = ["O"] * len(toks)
    used = [False] * len(toks)
    # find each span (in source order first, anywhere otherwise), case-insensitive, token-exact
    for lab, text, _ in sp:
        want = [t.lower() for t in tokenize(text)]
        n, hit = len(want), None
        for i in range(len(low) - n + 1):
            if low[i:i + n] == want and not any(used[i:i + n]):
                hit = i
                break
        if hit is None:
            return None
        tags[hit] = "B-" + lab
        for j in range(hit + 1, hit + n):
            tags[j] = "I-" + lab
        for j in range(hit, hit + n):
            used[j] = True
    # numbers outside labelled spans must have been in the source sentence too
    src_nums = sorted(t for t, g in zip(d["tokens"], d["tags"]) if g == "O" and re.match(r"\d", t))
    new_nums = sorted(t for t, u in zip(toks, used) if not u and re.match(r"\d", t))
    pool = list(src_nums)
    for x in new_nums:
        if x in pool:
            pool.remove(x)
        else:
            return None
    if _norm(plain) == _norm(d["text"]) or len(toks) < 3:
        return None
    return toks, tags


def ask(model, endpoint, prompt, temperature, retries=3):
    body = json.dumps({"model": model, "prompt": prompt, "stream": False, "keep_alive": "2h",
                       "options": {"temperature": temperature, "num_predict": 120, "num_ctx": 1024}}).encode()
    for k in range(retries):
        try:
            req = urllib.request.Request(endpoint, data=body, headers={"Content-Type": "application/json"})
            return json.loads(urllib.request.urlopen(req, timeout=180).read())["response"].strip()
        except Exception:
            if k == retries - 1:
                raise
            time.sleep(2 * (k + 1))


def work(job, a):
    idx, v, d = job
    src, _ = bracketed(d)
    style = STYLES[(idx * 3 + v) % len(STYLES)]
    try:
        para = ask(a.model, a.endpoint, PROMPT.format(style=style, src=src), a.temperature)
    except Exception:
        return idx, v, None, "error"
    para = (para or "").strip().strip('"').split("\n")[0].strip()
    r = relabel(para, d)
    if not r:
        return idx, v, None, "rejected"
    toks, tags = r
    text = re.sub(r"\[\[|\]\]", "", para).strip()
    return idx, v, {**d, "text": text, "tokens": toks, "tags": tags, "family": d["family"] + "+para",
                    "para_of": idx, "para_variant": v, "para_style": style}, "kept"


def selftest():
    d = {"text": "Simulate 200 RIL lines from two founders , seed 42 .",
         "tokens": ["Simulate", "200", "RIL", "lines", "from", "two", "founders", ",", "seed", "42", "."],
         "tags": ["O", "B-POP", "B-PRIMARY", "I-PRIMARY", "O", "B-FOUNDERS", "O", "O", "O", "O", "O"], "family": "x"}
    ok = [
        "could you make [[RIL lines]] for me? [[200]] of them from [[two]] founders",
        "Can you simulate 200 ril lines with two founders (seed 42)?",          # brackets dropped, case changed
    ]
    bad = [
        "Simulate 300 RIL lines from two founders.",                             # number changed
        "Simulate 200 RILs from two founders.",                                  # span altered
        "Simulate 200 RIL lines from two founders, 5 reps.",                      # new number
    ]
    for p in ok:
        assert relabel(p, d), p
    for p in bad:
        assert relabel(p, d) is None, p
    t, g = relabel(ok[1], d)
    assert [x for x, y in zip(t, g) if y != "O"] == ["200", "ril", "lines", "two"], (t, g)
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp"); ap.add_argument("--out")
    ap.add_argument("--n", type=int, default=10000, help="source sentences to paraphrase")
    ap.add_argument("--variants", type=int, default=2, help="paraphrases per source sentence")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--temperature", type=float, default=0.8)
    ap.add_argument("--max-minutes", type=float, default=0, help="stop after this long (0 = no limit)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--model", default=os.environ.get("SIMITALL_OLLAMA_MODEL", "llama3.1"))
    ap.add_argument("--endpoint", default="http://localhost:11434/api/generate")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    src = [json.loads(l) for l in open(a.inp) if l.strip()]
    random.Random(a.seed).shuffle(src)
    src = src[:a.n]
    tried_path = a.out + ".tried"
    tried = set()
    if os.path.exists(tried_path):
        tried = {tuple(json.loads(l)) for l in open(tried_path)}
    jobs = [(i, v, d) for i, d in enumerate(src) for v in range(a.variants) if (i, v) not in tried]
    try:
        ask(a.model, a.endpoint, "Say OK.", 0.0, retries=1)
    except Exception as e:
        sys.exit(f"Ollama not reachable at {a.endpoint} ({e}). Start it with `ollama serve` and `ollama pull {a.model}`.")
    print(f"{len(jobs)} paraphrase jobs left ({len(tried)} already tried) with {a.model}", flush=True)

    counts = {"kept": 0, "rejected": 0, "error": 0}
    t0 = time.time()
    deadline = t0 + a.max_minutes * 60 if a.max_minutes else None
    with open(a.out, "a") as f, open(tried_path, "a") as ft, cf.ThreadPoolExecutor(a.workers) as ex:
        pending = set()
        it = iter(jobs)
        for _ in range(a.workers * 2):
            j = next(it, None)
            if j:
                pending.add(ex.submit(work, j, a))
        while pending:
            done, pending = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                idx, v, row, status = fut.result()
                counts[status] += 1
                if status != "error":
                    ft.write(json.dumps([idx, v]) + "\n"); ft.flush()
                if row:
                    f.write(json.dumps(row) + "\n"); f.flush()
                n = sum(counts.values())
                if n % 50 == 0:
                    rate = n / max(time.time() - t0, 1) * 3600
                    print(f"{n}/{len(jobs)}: kept {counts['kept']} ({counts['kept'] / n:.0%}), rejected {counts['rejected']}, "
                          f"errors {counts['error']}, {rate:.0f}/hour", flush=True)
                if not (deadline and time.time() > deadline):
                    j = next(it, None)
                    if j:
                        pending.add(ex.submit(work, j, a))
    print(f"done: {counts} -> {a.out}" + (" (time limit reached; rerun to continue)" if deadline and time.time() > deadline else ""))


if __name__ == "__main__":
    main()
