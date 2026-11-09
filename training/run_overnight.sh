#!/usr/bin/env bash
# Overnight training run for the SimItAll prompt model.
#   generate -> paraphrase with Ollama (time-limited) -> train tiny model (+ optional BiomedBERT)
#   -> predict the 229 evaluation prompts -> score lexicon vs lexicon+model.
#
# From the SimItAll repo root:
#   bash training/run_overnight.sh
# Options (environment variables):
#   PY=python3            Python with torch, transformers, tokenizers
#   MODEL=llama3.1        Ollama model for paraphrases (qwen2.5:7b follows rules better; llama3.2:3b is faster)
#   PARA_MINUTES=360      stop paraphrasing after this many minutes, then train
#   N_PARA=12000          source sentences sent to Ollama (x VARIANTS paraphrases each)
#   N_TRAIN=60000         template sentences
#   BASE=microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract   also fine-tune a pretrained model (optional)
#   SKIP_PARA=1           skip Ollama entirely (templates only)
#   DEVICE=auto           auto | cpu | mps | cuda
# Everything is resumable: rerun the same command and finished steps are reused.
set -euo pipefail

# keep a Mac awake for the whole run
if [[ "$(uname)" == "Darwin" && -z "${CAFFEINATED:-}" ]]; then
  CAFFEINATED=1 exec caffeinate -dimsu bash "$0" "$@"
fi

PY=${PY:-python3}
MODEL=${MODEL:-llama3.1}
PARA_MINUTES=${PARA_MINUTES:-360}
N_PARA=${N_PARA:-12000}
VARIANTS=${VARIANTS:-2}
WORKERS=${WORKERS:-4}
N_TRAIN=${N_TRAIN:-60000}
N_DEV=${N_DEV:-3000}
EPOCHS_TINY=${EPOCHS_TINY:-10}
EPOCHS_BASE=${EPOCHS_BASE:-4}
BASE=${BASE:-}
DEVICE=${DEVICE:-auto}
D=training/data
L=training/logs
M=training/models
EVAL=inst/agent/prompt_eval/study_spec_examples_draft.jsonl
mkdir -p "$D" "$L" "$M"
log() { echo "[$(date '+%a %H:%M:%S')] $*" | tee -a "$L/overnight.log"; }

[[ -f DESCRIPTION && -d inst/agent/lexicon ]] || { echo "Run from the SimItAll repo root."; exit 1; }
$PY -c "import torch, transformers, tokenizers" 2>/dev/null || {
  echo "Install first:  $PY -m pip install --upgrade torch transformers tokenizers"; exit 1; }
log "start: python $($PY -c 'import torch; print(torch.__version__)'), device=$DEVICE, Ollama model=$MODEL"

# 1. template sentences ------------------------------------------------------
[[ -s $D/train.jsonl ]] || $PY training/generate.py --n "$N_TRAIN" --seed 1 --out $D/train.jsonl | tee -a "$L/overnight.log"
[[ -s $D/dev.jsonl ]]   || $PY training/generate.py --n "$N_DEV"   --seed 2 --out $D/dev.jsonl   | tee -a "$L/overnight.log"

# 2. paraphrases with Ollama (time-limited, resumable) ------------------------
if [[ -z "${SKIP_PARA:-}" ]]; then
  if command -v ollama >/dev/null 2>&1; then
    if ! ollama list >/dev/null 2>&1; then
      log "starting ollama serve (OLLAMA_NUM_PARALLEL=$WORKERS)"
      OLLAMA_NUM_PARALLEL=$WORKERS nohup ollama serve > "$L/ollama.log" 2>&1 &
      sleep 8
    fi
    ollama pull "$MODEL" >> "$L/ollama.log" 2>&1 || true
    log "paraphrasing up to $N_PARA x $VARIANTS sentences for at most $PARA_MINUTES minutes"
    $PY training/paraphrase_ollama.py --inp $D/train.jsonl --out $D/para.jsonl --n "$N_PARA" --variants "$VARIANTS" \
      --workers "$WORKERS" --model "$MODEL" --max-minutes "$PARA_MINUTES" 2>&1 | tee -a "$L/paraphrase.log" | tail -n 3
  else
    log "ollama not found: skipping paraphrases (install from https://ollama.com to add them)"
  fi
fi

# 3. mix: templates + paraphrases; hold out paraphrases of 1 in 20 source sentences as a harder dev set
$PY - <<'EOF'
import json, random
rows = [l for l in open("training/data/train.jsonl")]
para_train, para_dev = [], []
try:
    for l in open("training/data/para.jsonl"):
        (para_dev if json.loads(l)["para_of"] % 20 == 0 else para_train).append(l)
except FileNotFoundError:
    pass
mix = rows + para_train * 2          # paraphrases are scarcer and more varied: count them twice
random.Random(0).shuffle(mix)
open("training/data/train_mix.jsonl", "w").writelines(mix)
open("training/data/dev_para.jsonl", "w").writelines(para_dev or open("training/data/dev.jsonl").readlines()[:500])
print(f"mixed training rows: {len(mix)} ({len(para_train)} paraphrases x2); paraphrase dev rows: {len(para_dev)}")
EOF
cat $D/dev.jsonl $D/dev_para.jsonl > $D/dev_all.jsonl

# 4. train ---------------------------------------------------------------------
log "training tiny model ($EPOCHS_TINY epochs)"
$PY training/train.py --train $D/train_mix.jsonl --dev $D/dev_all.jsonl --out $M/simitall-tagger-tiny \
  --epochs "$EPOCHS_TINY" --device "$DEVICE" 2>&1 | { grep -v -i warn || true; } | tee "$L/train_tiny.log"
MODELS=(simitall-tagger-tiny)
if [[ -n "$BASE" ]]; then
  log "fine-tuning $BASE ($EPOCHS_BASE epochs)"
  $PY training/train.py --train $D/train_mix.jsonl --dev $D/dev_all.jsonl --out $M/simitall-tagger-base --base "$BASE" \
    --epochs "$EPOCHS_BASE" --bs 32 --device "$DEVICE" 2>&1 | { grep -v -i warn || true; } | tee "$L/train_base.log"
  MODELS+=(simitall-tagger-base)
fi

# 5. predict the evaluation prompts and score -----------------------------------
for m in "${MODELS[@]}"; do
  $PY inst/python/prompt_model.py --model $M/$m --examples $EVAL --out $L/preds_$m.json
  for s in all test; do
    arg=""; [[ $s == test ]] && arg=test
    Rscript inst/agent/prompt_eval/score_model.R $L/preds_$m.json $arg 2>&1 | tee "$L/score_${m}_$s.txt"
  done
done
log "done. Scores: $L/score_*_test.txt  (frozen test split).  Models: $M/"
log "If hybrid beats lexicon on the test split with no more invented values, copy the tiny model:"
log "  rm -rf inst/models/simitall-tagger && cp -r $M/simitall-tagger-tiny inst/models/simitall-tagger"
