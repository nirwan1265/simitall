# old_model: the prompt reader and Shiny agent (before Biodes)

Saved copy of everything removed on 2026-11-22 when Biodes took over reading
the user's question. Nothing here is built, installed, loaded or tested: the
folder is in `.Rbuildignore`, `devtools::load_all()` only reads `R/`, and
testthat only runs `tests/testthat/`.

Paths mirror the package layout:

| here | what it was |
| --- | --- |
| `R/agent_spec.R` | keyword spec layer: extract -> validate -> resolve -> compile plan + code; guided questions; option comparison |
| `R/agent_api.R` | `simitall_ask()`, `launch_simitall_agent()`, knowledge search, code validator |
| `R/agent_workflow_compiler.R` | legacy router |
| `R/agent_lexicon.R`, `R/agent_reader.R`, `inst/agent/lexicon/` | lexicon reader |
| `R/agent_model.R`, `inst/models/`, `inst/python/prompt_model.py`, `training/` | BERT prompt tagger and its training |
| `R/agent_capabilities.R`, `inst/agent/capabilities/` | capability catalog and StudySpec schema |
| `R/agent_evaluation.R`, `inst/agent/prompt_eval/`, `data-raw/agent_eval/` | prompt evaluation sets |
| `inst/agent/aliases.tsv`, `data-raw/lexicon/` | aliases and lexicon builders |
| `inst/shiny/simitall-agent/` | the Shiny agent app (Seed2Grow demo) |
| `tests/testthat/test-agent-*.R`, `man/` | their tests and help pages |

## Running the old agent

The code needs the package around it, so install the last version that had it
(commit `4bfddcf`):

```r
remotes::install_github("nirwan1265/simitall@4bfddcf")
simitall::launch_simitall_agent()
```

To bring a piece back into the package, copy it from here to the same path
outside `old_model/` and run `devtools::document()`.
