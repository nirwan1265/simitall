# Training the SimItAll prompt model

The lexicon (`inst/agent/lexicon/`, four tables) reads organisms, techniques,
settings, numbers with units and intent. This small BERT model learns what the
tables cannot: **which number belongs to which field** ("6 healthy and 6 sick
people" -> 12 donors, "from two parents" -> founders), **condition and cell-type
names**, and **the role of each analysis mentioned** (primary step, later step,
input such as "my GWAS hits", or excluded such as "no eQTL"). It also predicts
goal, data source and primary analysis for the whole request.

Rules first: `R/agent_model.R` (`.simitall_model_merge`) only fills fields the
lexicon left empty, reorders analyses the lexicon already found, and adds
nothing when the lexicon decided to stop (unsupported technique, impossible
value, clinical use, real data missing, mid-sentence correction).

## Files

| File | Purpose |
|---|---|
| `generate.py` | Synthetic sentences with labels by construction, 35 families covering all 31 catalog analyses, combinations, inputs, exclusions, unsupported techniques, vague requests and demo data. Vocabulary comes from the lexicon tables. Sentences too similar to any evaluation prompt are dropped. |
| `paraphrase_ollama.py` | Free local-LLM rewrites in 7 writing styles. Labels survive if every labelled phrase is found again (brackets may be dropped, case may change) and no new number appears. Resumable (rejects are remembered), time-limited, `--selftest` |
| `train.py` | Multi-task model: token roles + goal + data source + primary analysis. Tiny BERT from scratch, or `--base` a pretrained model. CUDA, Apple MPS or CPU. Keeps the best dev checkpoint |
| `run_overnight.sh` | One command for the whole night: generate -> paraphrase -> train -> predict the 229 evaluation prompts -> score |
| `../inst/python/prompt_model.py` | Labels, model, inference (`--text` or `--examples`) |
| `../inst/agent/prompt_eval/score_model.R` | Lexicon vs lexicon + model with the same scorer as `score_examples.R` |

## Overnight run (Mac)

```bash
cd simitall && git checkout agent-lexicon && git pull
source ~/sfenv/bin/activate            # the venv with torch + transformers from StudyForge
python -m pip install --upgrade torch transformers tokenizers
ollama pull llama3.1                   # once
BASE=microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract bash training/run_overnight.sh
```

The script keeps the Mac awake (`caffeinate`), starts `ollama serve` if needed,
paraphrases for at most `PARA_MINUTES` (default 360), trains the tiny model
(always) and the BiomedBERT model (when `BASE` is set), then scores both.
Rerunning the same command resumes. Logs and scores: `training/logs/`; models:
`training/models/` (both git-ignored). Other knobs: `MODEL`, `N_PARA`,
`N_TRAIN`, `WORKERS`, `EPOCHS_TINY`, `EPOCHS_BASE`, `DEVICE`, `SKIP_PARA=1`.

In the morning look at `training/logs/score_simitall-tagger-tiny_test.txt` (and
`..._base_test.txt`). If the hybrid column beats the lexicon column on the test
split without more invented values, ship the tiny model:

```bash
rm -rf inst/models/simitall-tagger && cp -r training/models/simitall-tagger-tiny inst/models/simitall-tagger
```

## Results so far (tiny model, 30k template sentences, 4 epochs on 2 CPUs, no paraphrases)

| | lexicon | lexicon + model |
|---|---|---|
| reading recall, all 229 / test 152 | 78.7 / 79.6 | 84.0 / 83.9 |
| primary analysis | 75.0 / 77.1 | 81.2 / 81.4 |
| goal and data source | 80.0 / 80.5 | 81.4 / 81.5 |
| invented values (examples) | 1 / 1 | 2 / 2 |
| improved / made worse (all) | | 34 / 0 |

Caveat: the merge guards were adjusted after reading failures on all 229 examples
(test split included), and the same author wrote the examples, lexicon and
generator, so these are optimistic. The honest check is new prompts written by
someone else after the model is frozen.
