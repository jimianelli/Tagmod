# Revival assessment and provenance

Date: 2026-07-29

## Authoritative runnable baseline

`2018 tag model/Main` is the current runnable baseline. Its one-area model
(`src/tm.tpl`) builds with ADMB 12.3 on Apple Silicon and reproduces the archived
`both` fit.

The regenerated and archived fits have:

- objective function: `810.919633885051`
- `lnNinit`: `6.90766133005`
- maximum-gradient difference only at floating-point roundoff

The archived report's zero biomass was caused by a derived output being assigned
only in ADMB standard-deviation and MCMC phases. The revived template assigns
derived outputs in all phases. The regenerated biomass is `652.939`.

## Historical Git repository

The original local repository at `NPRB Final Report Analysis/Tagmod` contains
macOS `Icon` metadata files inside `.git/refs` and `.git/objects`. These produce
invalid Git ref names. Its working tree also contains deleted tracked source
files and numerous untracked ADMB outputs.

A clean clone is stored at:

`NPRB Final Report Analysis/Tagmod-recovered`

No remote branches were modified. The remote branches preserved in the clone
are:

- `Dev`
- `master`
- `sex`

The `sex` branch contains the most recent historical commits and the June 2018
sex-structured model work. It should be evaluated as a separate model lineage,
not silently merged into the revived one-area model.

## Directory roles

- `2018 tag model/Main`: executable model and archived 2018 results
- `2018 tag model/misc files`: raw/intermediate data, GIS material, and notes
- `NPRB Final Report Analysis`: final-report inputs and analysis provenance
- `NPRB Final Report Analysis/Tagmod`: damaged local Git working copy; preserve
- `NPRB Final Report Analysis/Tagmod-recovered`: clean remote history

## Remaining work before scientific reuse

1. Recover or reconstruct the missing generating input/template for the
   `runs/results/Seguam_M` result.
2. Identify which one-area datasets and model variants support each reported
   conclusion.
3. Run short MCMC diagnostics before committing resources to the historical
   one-million-iteration runs.
4. Review recruitment-factor and natural-mortality sensitivity assumptions.
5. Reconcile the sex-structured Petrel input/revision discrepancy and treat
   `tm3.tpl` as a distinct candidate requiring
   its own validation.
6. Build a data dictionary connecting spreadsheet/GIS source fields to each
   `.dat` input.

See `VALIDATION_REPORT.md` for the completed mode-fit comparison.
