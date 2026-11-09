"""SimItAll prompt model: a small BERT that tags the roles of words in a request.

The lexicon (inst/agent/lexicon) already recognises organisms, techniques,
settings and intent. This model learns what the lexicon cannot: which number
belongs to which field ("6 healthy and 6 sick people" -> donors), which
condition names were given, and which mentioned analysis is the requested
(primary) step, a later (secondary) step, an input ("my GWAS hits") or excluded
("no eQTL"). It also predicts goal, data source and primary analysis for the
whole request as a fallback.

    python inst/python/prompt_model.py --model inst/models/simitall-tagger --text "..."
    python inst/python/prompt_model.py --model inst/models/simitall-tagger --examples X.jsonl --out preds.json
"""
import argparse, json, os, re

TOKEN_RE = re.compile(r"\d+(?:[.,]\d+)*|[A-Za-z]+(?:[-'/][A-Za-z0-9]+)*|[^\sA-Za-z\d]")

# token tags: numbers by field, named spans, and analysis roles
NUM_TAGS = ["POP", "DONORS", "FOUNDERS", "GENS", "SELFGEN", "BCGEN", "REPS", "TECHREPS", "TIMEPTS", "CASES",
            "CONTROLS", "ENVS", "BATCHES", "NCAUSAL", "CELLS", "NCELLTYPES", "NCOND", "COV", "H2"]
SPAN_TAGS = ["COND", "CELLTYPE", "ORG"]
ROLE_TAGS = ["PRIMARY", "SECONDARY", "INPUT", "NEG"]
TAGS = NUM_TAGS + SPAN_TAGS + ROLE_TAGS
LABELS = ["O"] + [p + t for t in TAGS for p in ("B-", "I-")]

GOALS = ["simulate", "analyze", "simulate_and_analyze", "plan", "none"]
SOURCES = ["provided", "synthetic", "package_demo", "ask"]
ANALYSES = ["pedigree_inbreeding", "population_relatedness_gwas_calibration", "cross_design", "phenotype_outcome_model",
            "quantitative_gwas", "binary_gwas", "multitrait_longitudinal_gwas", "rare_variant_and_gene_based_gwas",
            "linkage_qtl", "bulk_counts_simulation", "differential_expression", "time_series_rnaseq",
            "alternative_splicing", "allele_specific_expression", "rnaseq_power", "bulk_eqtl", "splicing_qtl",
            "single_cell_eqtl", "scrnaseq_simulation", "scrnaseq_analysis", "chipseq_atacseq", "genomic_prediction",
            "population_structure", "linkage_disequilibrium", "selection_scan", "read_simulation", "genome_assembly",
            "variant_calling", "gene_coordinate_annotation", "gxe_multiomics", "microbiome_analysis",
            "none", "unsupported"]


def tokenize(text):
    return TOKEN_RE.findall(text)


def spans(tokens, tags):
    """[(tag, text, start_index)] from BIO tags."""
    out, cur = [], None
    for i, (tok, tag) in enumerate(zip(tokens, tags)):
        if tag.startswith("B-") or (tag.startswith("I-") and (cur is None or cur[0] != tag[2:])):
            if cur:
                out.append(cur)
            cur = [tag[2:], [tok], i]
        elif tag.startswith("I-") and cur:
            cur[1].append(tok)
        else:
            if cur:
                out.append(cur)
            cur = None
    if cur:
        out.append(cur)
    return [(t, " ".join(w), s) for t, w, s in out]


def _torch():
    import torch
    import torch.nn as nn
    return torch, nn


def make_model(encoder, hidden):
    torch, nn = _torch()

    class Tagger(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = encoder
            self.drop = nn.Dropout(0.1)
            self.tok_head = nn.Linear(hidden, len(LABELS))
            self.goal_head = nn.Linear(hidden, len(GOALS))
            self.source_head = nn.Linear(hidden, len(SOURCES))
            self.primary_head = nn.Linear(hidden, len(ANALYSES))

        def forward(self, input_ids, attention_mask):
            h = self.drop(self.encoder(input_ids=input_ids, attention_mask=attention_mask).last_hidden_state)
            c = h[:, 0]
            return self.tok_head(h), self.goal_head(c), self.source_head(c), self.primary_head(c)

    return Tagger()


def encode(tokenizer, toks, labels=None, max_len=160):
    enc = tokenizer(toks, is_split_into_words=True, truncation=True, max_length=max_len)
    wid = enc.word_ids()
    lab, prev = [], None
    for w in wid:
        if w is None or w == prev:
            lab.append(-100)
        else:
            lab.append(LABELS.index(labels[w]) if labels else 0)
        prev = w
    return enc["input_ids"], enc["attention_mask"], lab, wid


def pick_device(name="auto"):
    torch, _ = _torch()
    if name != "auto":
        return name
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


class PromptTagger:
    def __init__(self, model_dir, device="cpu"):
        torch, _ = _torch()
        from transformers import AutoConfig, AutoModel, AutoTokenizer, PreTrainedTokenizerFast
        meta = json.load(open(os.path.join(model_dir, "meta.json")))
        self.tok = PreTrainedTokenizerFast.from_pretrained(model_dir) if meta.get("custom_tokenizer") \
            else AutoTokenizer.from_pretrained(model_dir)
        cfg = AutoConfig.from_pretrained(model_dir)
        enc = AutoModel.from_config(cfg, add_pooling_layer=False) if cfg.model_type == "bert" else AutoModel.from_config(cfg)
        self.model = make_model(enc, cfg.hidden_size)
        self.model.load_state_dict(torch.load(os.path.join(model_dir, "tagger.pt"), map_location=device))
        self.model.eval()
        self.lower = meta.get("lowercase", True)

    def predict(self, text):
        torch, _ = _torch()
        toks = tokenize(text)
        inp = [t.lower() for t in toks] if self.lower else toks
        ids, mask, _, wid = encode(self.tok, inp)
        with torch.no_grad():
            tl, gl, sl, pl = self.model(torch.tensor([ids]), torch.tensor([mask]))
        pred = tl[0].argmax(-1).tolist()
        tags, prev = ["O"] * len(toks), None
        for i, w in enumerate(wid):
            if w is not None and w != prev:
                tags[w] = LABELS[pred[i]]
            prev = w

        def head(logits, names, k=1):
            p = logits[0].softmax(-1)
            top = p.topk(k)
            return [[names[int(i)], round(float(v), 3)] for v, i in zip(top.values, top.indices)]

        return {"spans": [list(s) for s in spans(toks, tags)],
                "goal": head(gl, GOALS)[0], "data_source": head(sl, SOURCES)[0],
                "primary": head(pl, ANALYSES)[0], "primary_top3": head(pl, ANALYSES, 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--text"); ap.add_argument("--examples"); ap.add_argument("--out")
    a = ap.parse_args()
    tg = PromptTagger(a.model)
    if a.text:
        print(json.dumps(tg.predict(a.text)))
        return
    prompts = [json.loads(l)["prompt"] for l in open(a.examples) if l.strip()]
    preds = {p: tg.predict(p) for p in prompts}
    json.dump(preds, open(a.out, "w"), indent=1)
    print(f"wrote {len(preds)} predictions to {a.out}")


if __name__ == "__main__":
    main()
