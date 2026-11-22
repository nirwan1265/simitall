#!/usr/bin/env python3
"""Train the SimItAll prompt model (token roles + goal + data source + primary analysis).

Tiny BERT from scratch (fast, ~8 MB, ships with the package):
    python training/train.py --train training/data/train.jsonl --dev training/data/dev.jsonl --out inst/models/simitall-tagger
Fine-tune a pretrained biomedical BERT (better on unusual wording, ~440 MB, downloads once):
    python training/train.py --base microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract --epochs 4 --out models/simitall-biomedbert ...
Device is chosen automatically: CUDA, then Apple MPS, then CPU (--device to force).
The checkpoint with the best dev score is kept.
"""
import argparse, json, os, random, sys, time

import torch
from torch.nn.utils.rnn import pad_sequence

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "inst", "python"))
from prompt_model import ANALYSES, GOALS, LABELS, SOURCES, encode, make_model, pick_device  # noqa: E402


def load(path):
    return [json.loads(l) for l in open(path) if l.strip()]


def build_tiny(train, vocab_size, hidden, layers, heads):
    from tokenizers import Tokenizer, models, normalizers, pre_tokenizers, trainers
    from tokenizers.processors import TemplateProcessing
    from transformers import BertConfig, BertModel, PreTrainedTokenizerFast
    tk = Tokenizer(models.WordPiece(unk_token="[UNK]"))
    tk.normalizer = normalizers.BertNormalizer(lowercase=True)
    tk.pre_tokenizer = pre_tokenizers.BertPreTokenizer()
    trainer = trainers.WordPieceTrainer(vocab_size=vocab_size, special_tokens=["[PAD]", "[UNK]", "[CLS]", "[SEP]", "[MASK]"])
    # vocabulary from training text plus every lexicon alias, so unseen names stay in-vocab
    lex = os.path.join(ROOT, "inst", "agent", "lexicon")
    extra = []
    for f in ("organism_aliases.tsv", "technique_aliases.tsv", "settings.tsv"):
        extra += [l.split("\t")[0] for l in open(os.path.join(lex, f))][1:]
    tk.train_from_iterator([d["text"] for d in train] + extra * 3, trainer)
    tk.post_processor = TemplateProcessing(single="[CLS] $A [SEP]", special_tokens=[
        ("[CLS]", tk.token_to_id("[CLS]")), ("[SEP]", tk.token_to_id("[SEP]"))])
    hf = PreTrainedTokenizerFast(tokenizer_object=tk, unk_token="[UNK]", pad_token="[PAD]", cls_token="[CLS]",
                                 sep_token="[SEP]", mask_token="[MASK]")
    cfg = BertConfig(vocab_size=hf.vocab_size, hidden_size=hidden, num_hidden_layers=layers, num_attention_heads=heads,
                     intermediate_size=hidden * 4, max_position_embeddings=256)
    return hf, BertModel(cfg, add_pooling_layer=False), cfg


def batches(data, tok, bs, shuffle, lower):
    idx = list(range(len(data)))
    if shuffle:
        random.shuffle(idx)
    for i in range(0, len(idx), bs):
        items = [data[j] for j in idx[i:i + bs]]
        ids, mask, lab, g, s, p = [], [], [], [], [], []
        for d in items:
            toks = [t.lower() for t in d["tokens"]] if lower else d["tokens"]
            a, m, l, _ = encode(tok, toks, d["tags"])
            ids.append(torch.tensor(a)); mask.append(torch.tensor(m)); lab.append(torch.tensor(l))
            g.append(GOALS.index(d["goal"])); s.append(SOURCES.index(d["data_source"])); p.append(ANALYSES.index(d["primary"]))
        yield (pad_sequence(ids, batch_first=True, padding_value=tok.pad_token_id), pad_sequence(mask, batch_first=True),
               pad_sequence(lab, batch_first=True, padding_value=-100), torch.tensor(g), torch.tensor(s), torch.tensor(p))


@torch.no_grad()
def evaluate(model, data, tok, lower, dev):
    model.eval()
    tp = fp = fn = 0
    ok = {"goal": 0, "source": 0, "primary": 0}
    for ids, mask, lab, g, s, p in batches(data, tok, 64, False, lower):
        tl, gl, sl, pl = model(ids.to(dev), mask.to(dev))
        pred = tl.argmax(-1).cpu()
        keep = lab != -100
        pr, gt = pred[keep], lab[keep]
        tp += ((pr == gt) & (gt != 0)).sum().item()
        fp += ((pr != gt) & (pr != 0)).sum().item()
        fn += ((pr != gt) & (gt != 0)).sum().item()
        ok["goal"] += (gl.argmax(-1).cpu() == g).sum().item()
        ok["source"] += (sl.argmax(-1).cpu() == s).sum().item()
        ok["primary"] += (pl.argmax(-1).cpu() == p).sum().item()
    model.train()
    prec, rec = tp / max(tp + fp, 1), tp / max(tp + fn, 1)
    out = {"token_f1": round(2 * prec * rec / max(prec + rec, 1e-9), 4)}
    out.update({k + "_acc": round(v / len(data), 4) for k, v in ok.items()})
    out["score"] = round((out["token_f1"] * 2 + out["goal_acc"] + out["source_acc"] + out["primary_acc"]) / 5, 4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True); ap.add_argument("--dev", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--base", default=None, help="HF model name/path; omit to train a tiny BERT from scratch")
    ap.add_argument("--epochs", type=int, default=8); ap.add_argument("--bs", type=int, default=32)
    ap.add_argument("--lr", type=float, default=None); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--hidden", type=int, default=256); ap.add_argument("--layers", type=int, default=4)
    ap.add_argument("--heads", type=int, default=4); ap.add_argument("--vocab", type=int, default=6000)
    a = ap.parse_args()
    random.seed(a.seed); torch.manual_seed(a.seed)
    train, devset = load(a.train), load(a.dev)
    dev = pick_device(a.device)
    print(f"device {dev}; {len(train)} train / {len(devset)} dev sentences", flush=True)

    if a.base:
        from transformers import AutoConfig, AutoModel, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(a.base)
        cfg = AutoConfig.from_pretrained(a.base)
        enc = AutoModel.from_pretrained(a.base, add_pooling_layer=False) if cfg.model_type == "bert" else AutoModel.from_pretrained(a.base)
        custom, lower, lr = False, getattr(tok, "do_lower_case", True), a.lr or 3e-5
    else:
        tok, enc, cfg = build_tiny(train, a.vocab, a.hidden, a.layers, a.heads)
        custom, lower, lr = True, True, a.lr or 5e-4

    model = make_model(enc, cfg.hidden_size).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    steps = a.epochs * ((len(train) + a.bs - 1) // a.bs)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=lr, total_steps=steps, pct_start=0.1)
    ce = torch.nn.CrossEntropyLoss(ignore_index=-100)
    best, t0 = None, time.time()
    os.makedirs(a.out, exist_ok=True)
    for ep in range(a.epochs):
        tot = 0.0
        for ids, mask, lab, g, s, p in batches(train, tok, a.bs, True, lower):
            ids, mask, lab, g, s, p = (x.to(dev) for x in (ids, mask, lab, g, s, p))
            tl, gl, sl, pl = model(ids, mask)
            loss = ce(tl.reshape(-1, len(LABELS)), lab.reshape(-1)) + 0.5 * (ce(gl, g) + ce(sl, s) + ce(pl, p))
            opt.zero_grad(); loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step(); sched.step()
            tot += loss.item()
        m = evaluate(model, devset, tok, lower, dev)
        print(f"epoch {ep + 1}/{a.epochs}: loss {tot:.1f}  dev {m}  ({time.time() - t0:.0f}s)", flush=True)
        if best is None or m["score"] > best["score"]:
            best = dict(m, epoch=ep + 1)
            tok.save_pretrained(a.out); cfg.save_pretrained(a.out)
            torch.save({k: v.cpu() for k, v in model.state_dict().items()}, os.path.join(a.out, "tagger.pt"))
            json.dump({"labels": LABELS, "goals": GOALS, "sources": SOURCES, "analyses": ANALYSES, "custom_tokenizer": custom,
                       "lowercase": lower, "base": a.base, "dev": best, "n_train": len(train)},
                      open(os.path.join(a.out, "meta.json"), "w"), indent=2)
    print(f"best epoch {best['epoch']}: {best}; saved {a.out}")


if __name__ == "__main__":
    main()
