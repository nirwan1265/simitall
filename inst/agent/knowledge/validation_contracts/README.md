# SimItAll Validation Contracts

A validation contract is the four-part, workflow-specific definition of what a
simulation or analysis must retain and check:

1. **Simulation truth**: the programmed inputs and hidden truth used to create
   synthetic data. This is not an expected biological result.
2. **Required results**: files the workflow must write for reuse, auditing, and
   downstream analysis.
3. **Required figures**: reproducible figures generated from saved source data.
4. **Pass criteria**: structural or statistical checks that show the configured
   mechanism ran as intended. A pass does not establish real biological or
   clinical validity.

These cards are training and retrieval material for the local agent. They also
serve as the specification for fast demo tests in `tests/testthat/`. Tests use
bundled toy data and small dimensions; larger benchmark runs remain separate.

## Figure Rules

- Generate a figure only when the workflow has the data needed to support it.
- Save panel source data next to the figure and record the random seed.
- Label synthetic, bundled-demo, and user-supplied inputs clearly.
- Do not turn a simulation pass criterion into a claim about a real species,
  human trait, or named population.
