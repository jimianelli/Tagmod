# Atka mackerel tag model

This directory contains the revived 2018 AD Model Builder (ADMB) tag model,
its input data, archived results, and R analysis helpers.

## Requirements

- ADMB (tested with ADMB 12.3)
- POSIX shell (`sh`)
- R for result summaries and plots

## Quick start

From this directory:

```sh
./scripts/build.sh
./scripts/run-model.sh both
./scripts/smoke-test.sh
```

`run-model.sh` accepts the stem of any `.dat` file in `data/`. It performs a
mode fit by default and writes a new, timestamped result directory under
`runs/generated/`, leaving the historical results untouched.

To run MCMC:

```sh
./scripts/run-model.sh both --mcmc
```

The historical settings are one million iterations with every 200th sample
saved. They can be overridden for a diagnostic run:

```sh
MCMC_ITERATIONS=10000 MCMC_SAVE=20 ./scripts/run-model.sh both --mcmc
```

Summarize archived MCMC results:

```sh
Rscript R/plotTM.r runs/results runs/generated/archived-mcmc-summary.pdf
```

Compare newly generated mode fits with the archived parameter files:

```sh
Rscript R/validate-archive.R .
```

Rerun every archived mode-fit dataset and then validate it:

```sh
./scripts/reproduce-archive.sh
```

The validator exits nonzero when any fit differs. The currently known
`Seguam_M` provenance mismatch is documented in the workspace-level
`VALIDATION_REPORT.md`.

## Layout

- `src/tm.tpl`: current one-area ADMB model
- `src/tm2.tpl`, `src/tm3.tpl`: historical model variants
- `data/`: model inputs
- `runs/results/`: archived 2018 parameter files used as validation evidence
- `runs/generated/`: newly generated results (created on demand)
- `R/read.admb.R`: readers for ADMB output formats
- `R/plotTM.r`: reproducible MCMC summary plots
- `R/plotTM-legacy.r`: original exploratory analysis, retained for provenance

## Reproducibility notes

The original executable was an Intel-only macOS binary and is intentionally not
tracked here. `build.sh` rebuilds it for the current machine. The repository
root retains the original Windows batch files as historical records.

The model input convention uses a small `tm.dat` control file whose first line
is the path to the actual dataset. The revived runner creates this control file
inside an isolated working directory for every run.

## Known scientific checks

The revival corrects a reporting-only issue where `Biomass` was left at zero
during an ordinary fit. It does not alter the likelihood, priors, or population
dynamics. Before producing new scientific advice, compare regenerated estimates
against archived `.par`, `.std`, and `.rep` files and document which model
variant and input dataset are authoritative.
