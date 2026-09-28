# Subcortical brain iron development in psychosis-spectrum youth

Analysis code and de-identified data for Larsen et al. (in preparation), examining the
developmental trajectory of subcortical magnetic susceptibility — an MRI marker of brain iron
content (quantitative susceptibility mapping, QSM) — in youth with and without subthreshold
psychosis-spectrum symptoms, and its association with cognition.

This repository contains the code and released data used for all analyses, tables, and figures in
the manuscript. Running `QSM_Analysis_Submission_public.Rmd` reproduces them from the included
de-identified dataset.

## File organization

```
qsm-iron-psychosis/
  |- QSM_Analysis_Submission_public.Rmd   # Main analysis — reproduces all figures and tables
  |- data_public/                         # De-identified analysis dataset
  |     |- qsm_iron_psychosis_analysis_data.csv
  |     |- DATA_DICTIONARY.md             # Column definitions and de-identification notes
  |- R/                                   # Helper functions sourced by the Rmd
  |- LICENSE                              # MIT
  |- README.md
```

`figs/` and `tables/data/` are created automatically when the Rmd is run, and collect the figure
(SVG) and table (CSV) outputs.

## System requirements

- **Operating system:** macOS 12+, Ubuntu 20.04+, or Windows 10/11
- **R version:** R 4.3.1 (2023-06-16) or later
- **Platform tested:** aarch64-apple-darwin20 (macOS, Apple Silicon)
- **Bayesian backend:** the mediation models use [brms](https://paulbuerkner.com/brms/) with the
  **cmdstanr** backend, which requires a working **CmdStan** installation (see below).
- **R package dependencies:**
  `tidyverse`, `data.table`, `mgcv`, `gratia`, `broom`, `broom.mixed`, `lmerTest`, `kableExtra`,
  `patchwork`, `effectsize`, `ggsignif`, `ggeffects`, `brms`, `tidybayes`, `ggtext`, `scico`,
  `svglite`, `gtsummary` (plus `MASS`, `digest`, `knitr`, `rmarkdown`).

## Installation guide

1. Install R from [CRAN](https://cran.r-project.org/) (and, optionally, [RStudio](https://posit.co/download/rstudio-desktop/)).
2. Clone this repository:
   ```
   git clone https://github.com/Larsen-Lab/qsm-iron-psychosis.git
   cd qsm-iron-psychosis
   ```
3. Install the required R packages:
   ```r
   install.packages(c("tidyverse", "data.table", "mgcv", "gratia", "broom", "broom.mixed",
                      "lmerTest", "kableExtra", "patchwork", "effectsize", "ggsignif",
                      "ggeffects", "brms", "tidybayes", "ggtext", "scico", "svglite",
                      "gtsummary", "MASS", "digest", "knitr", "rmarkdown"))
   ```
4. Install the CmdStan toolchain used by the Bayesian mediation models:
   ```r
   install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
   cmdstanr::install_cmdstan()
   ```

Typical install time: ~15 minutes on a normal desktop computer (longer if CmdStan compiles from
source).

## Data

The **de-identified analysis dataset** required to reproduce all reported results is included in
`data_public/qsm_iron_psychosis_analysis_data.csv` (263 imaging sessions from 183 participants);
see `data_public/DATA_DICTIONARY.md` for column definitions. Participant and session identifiers
were replaced with fresh, non-linkable IDs, and variables not used in the analyses (e.g., detailed
demographics, raw neurocognitive subtests) were removed.

## Running the code

Open `QSM_Analysis_Submission_public.Rmd` in RStudio and knit it (or run
`rmarkdown::render("QSM_Analysis_Submission_public.Rmd")`). It loads the included
dataset and reproduces all manuscript figures (`figs/`) and tables (`tables/data/`).

Expected run time: roughly 30–60 minutes on a laptop, dominated by compiling and sampling the
Bayesian mediation models; the non-Bayesian analyses complete in a few minutes.

## Reproducibility note

Deterministic results (mixed-model coefficients, F/t statistics, p-values, FDR) reproduce exactly.
Simulation-based quantities — the posterior divergence-onset ages and the Bayesian mediation
estimates — reproduce to within Monte Carlo error and may differ from the manuscript in the last
reported digit depending on platform, random-number stream, and Stan/CmdStan version.

## License

Released under the MIT License (see `LICENSE`).
