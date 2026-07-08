# replicateAGCApaper

Author: Alberto Quaini

`replicateAGCApaper` is the replication compendium for the paper on anchored
geodesic component analysis (AGCA) for multivariate extremes.

This repository is intentionally separate from `AGCA4extremes`.

- `AGCA4extremes` is the reusable, CRAN-oriented methods package.
- `replicateAGCApaper` stores the paper-specific raw data, transformed data,
  scripts, generated results, and figures.

The compendium is large because it includes the raw Open Source Asset Pricing
daily portfolio archive and the generated empirical outputs. It is not intended
for CRAN submission.

## Development Assistance

OpenAI Codex was used as a programming assistant during development, mainly for
code scaffolding, refactoring, documentation, and tests. All methodological
choices, validation, final code, and responsibility for the package remain with
the author.

## Structure

```text
data-raw/
  empirics/
    ff/       # original Fama-French zip archives
    osap/     # original Open Source Asset Pricing zip archive

data/
  empirics/
    ff/       # transformed Fama-French RDS/CSV files
    osap/     # transformed OSAP RDS/CSV files

inst/
  empirics/
    scripts/  # paper-specific empirical scripts
    results/  # generated empirical CSV/PDF outputs
  simulations/
    scripts/  # paper-specific simulation scripts
    results/  # generated simulation CSV/PDF outputs

R/
  GeodesicExtreme.R  # legacy helper used by the paper scripts
  paths.R
  run.R

scripts/
  run_data.R
  run_empirics.R
  run_simulations.R
  run_all.R
```

## Dependencies

Install the companion methods package first:

```r
devtools::install("../AGCA4extremes")
```

The historical paper scripts also use standard plotting and data-manipulation
packages available in the project R library. If a script reports a missing
package, install it before rerunning the full workflow.

## Reproducing the Data

Small raw downloaded archives are stored in `data-raw/`. The large Open Source
Asset Pricing archive is downloaded by the data script when missing because it
exceeds GitHub's file-size limit. To regenerate the transformed empirical
datasets:

```sh
Rscript scripts/run_data.R
```

The data scripts write transformed files to `data/empirics/ff` and
`data/empirics/osap`. Use `--refresh` when running the underlying scripts if
you want to redownload the raw archives.

## Reproducing the Simulations

```sh
Rscript scripts/run_simulations.R
```

The default simulation scripts include bootstrap and finite-sample loops and
can take time. The underlying scripts accept command-line options such as
`--skip-bootstrap`, `--bootstrap-reps=`, `--reps=`, `--n=`, and `--k=`.

## Reproducing the Empirical Figures and Tables

```sh
Rscript scripts/run_empirics.R
```

Outputs are written to `inst/empirics/results`.

## Full Workflow

```sh
Rscript scripts/run_all.R
```

The full workflow regenerates transformed data, simulations, and empirical
outputs. Existing results are included so that the paper can be inspected
without rerunning the complete computation.
