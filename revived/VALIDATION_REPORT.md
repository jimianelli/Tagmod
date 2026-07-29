# Atka tag-model validation

Date: 2026-07-29

## One-area 2018 model

All 13 archived result names have a case-insensitive match to a surviving input
in `2018 tag model/Main/data`. Every input completed an ordinary ADMB fit.

Twelve fits reproduce their archived objective and parameter estimates to
numerical precision:

- `both`
- `both_no_RF`
- `nearshore`
- `nearshore_no_RF`
- `seamounts`
- `seamounts_no_RF`
- `Seguam`
- `Seguam_no_RF`
- `Seguam_no_RF_M`
- `Seguam1`
- `Seguam1b`
- `Seguam34`

The detailed machine-readable comparison is in:

`2018 tag model/Main/runs/generated/mode-validation.csv`

The formatted comparison is in:

`2018 tag model/Main/runs/generated/mode-validation.md`

### Unresolved `Seguam_M` result

`Seguam_M` does not reproduce the parameter file in `runs/results`:

- archived objective: `658.809295259956`
- regenerated objective: `692.186888823867`
- absolute difference: `33.377593563911`

The regenerated objective matches the separate
`runs/res_seguam_sensitivity/Seguam_M.par` objective exactly. This shows the
surviving `data/Seguam_M.dat` belongs to the sensitivity result rather than the
file archived under `runs/results/Seguam_M.par`.

The latter result therefore lacks its exact generating input or template. It
should not be presented as reproducible until that version is found in backups
or reconstructed from source spreadsheets.

`Seguam1In_34Out2.dat` is not a substitute: it fails the current template's
`icheck` data-layout guard.

## Sex-structured lineage

The latest remote `sex` branch was recovered and its `tm3.tpl`:

- compiles successfully with ADMB 12.3 on Apple Silicon;
- runs successfully with its tracked `Pet_All.dat`;
- converges with 27 parameters and a maximum gradient of approximately
  `2.26e-05`.

The recovered combination gives an objective of approximately `1656.466`.
The surviving historical `misc files/src/arc/Pet_All.par` reports `1655.83`.
Consequently, the sex-structured code is operational but its archived Petrel
fit is not yet exactly reproducible. The discrepancy is consistent with input
or revision drift among the loose working copies.

No archived parameter outputs are committed on the recovered remote branch, so
the loose `misc files/src/arc` outputs are the only comparison evidence found.

## Interpretation

The one-area model is strongly revived: 12 of 13 archived fits reproduce. The
exception is isolated to the provenance of one sensitivity dataset. The
sex-structured lineage is executable but remains a research candidate requiring
data/version reconciliation before reuse.
