# Showcase: simulated data for three designed experiments

Each folder takes one design written by [Biodes](https://github.com/nirwan1265/BioDesignGPT/tree/main/examples/showcase)
from a plain-language question, simulates data for that exact layout, and records everything it planted in
`sim_truth.json`: the true effects an analysis should find. [Bioana](https://github.com/nirwan1265/Bioana/tree/main/examples/showcase)
then analyses the same folder and scores itself against that truth. The folder names are the same in every repository.

| Example | Design (from Biodes) | What Simitall planted |
|---|---|---|
| [01 GWAS, sorghum panel](01_gwas_sorghum) | 300 accessions, alpha lattice, 2 replicates; `y ~ treatment + rep + (1 \| rep:block)` | 300 inbred lines x 1,905 markers with population structure, 10 additive QTNs, heritability 0.6 |
| [02 RNA-seq, genotype x drought](02_rnaseq_drought) | wild type / mutant x well-watered / drought, 4 plants each; `y ~ treatment * treatment2 + block` | 5,000 genes: 500 respond to drought in both genotypes, 500 respond differently in the mutant (the goal) |
| [03 Seed priming, time to germination](03_germination_survival) | 2 lines x primed / not, 12 trays, 10 seeds per cell; checked every 2 days to day 21 | hazard ratio 0.6 that starts after a delay, seed rot as a competing event, tray and cell frailty, censoring at day 21 |

Each folder's README gives the question, what was planted in words, one summary figure and the links on to Bioana.

Rebuild all three (seed 2026), from the repository root:

```sh
Rscript examples/showcase/make_showcase.R                       # reads ../BioDesignGPT/examples/showcase, else GitHub
BIODES_SHOWCASE=/path/to/BioDesignGPT/examples/showcase Rscript examples/showcase/make_showcase.R
```
