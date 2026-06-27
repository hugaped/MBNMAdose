# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`MBNMAdose` is an R package for Bayesian dose-response Model-Based Network Meta-Analysis (MBNMA). It fits models in JAGS (an external dependency, ≥4.3.0) via `R2jags`/`rjags`. Models synthesise relative effects across agents modelled through dose-response functions, operating on arm-level aggregate data (binomial, normal, or Poisson).

## Development Commands

All standard package operations run via `devtools` from within an R session:

```r
devtools::document()      # Regenerate documentation (roxygen2)
devtools::check()         # Full R CMD check
devtools::test()          # Run all tests
devtools::build()         # Build package tarball
devtools::install()       # Install locally
```

Run a single test file:
```r
testthat::test_file("tests/testthat/test_run.functions.R")
```

Most JAGS-dependent tests are skipped on CI (`skip_on_ci()`, `skip_on_cran()`). Only `test_write.functions.R`, `test_prepare.functions.R`, and `test_dose.function.R` run without JAGS.

## Architecture

### User-facing workflow
1. `mbnma.network(data.ab)` — validates and indexes arm-level data → `mbnma.network` S3 object
2. `mbnma.run(network, fun=...)` — writes JAGS code, runs MCMC → `mbnma` S3 object
3. `predict(mbnma, ...)` / `rank(mbnma, ...)` / `plot(mbnma, ...)` — downstream analysis

Supporting: `nma.run()` (treatment-level NMA), `nma.nodesplit()` (consistency checking).

### S3 Classes
| Class | Constructed by | Key contents |
|---|---|---|
| `mbnma.network` | `mbnma.network()` | `data.ab`, `agents`, `treatments`, `classes` |
| `mbnma` | `mbnma.run()` | Inherits `rjags`; adds `model.arg` (fun, jagscode, jagsdata, priors, etc.), `network`, `type="dose"` |
| `dosefun` | `demax()`, `dloglin()`, etc. | `name`, `params`, `apool`, `jags` (JAGS expression string), `fun` (R formula), `bname` |
| `nma` | `nma.run()` | `jagsresult`, `trt.labs`, `UME` |
| `mbnma.predict` | `predict.mbnma()` | Prediction samples per agent/dose |
| `mbnma.rank` | `rank.mbnma()` | Ranking probabilities and cumulative rank |

### Key source files
- `R/prepare.functions.R` — `mbnma.network()`, `mbnma.validate.data()`, `add_index()`, data indexing helpers
- `R/run.functions.R` — `mbnma.run()` (orchestrator), `mbnma.jags()` (internal JAGS runner), `nma.run()`, `check.likelink()`, `check.fun()`
- `R/dose.functions.R` — All `dosefun` constructors: `demax()`, `dloglin()`, `dexp()`, `dpoly()`, `dfpoly()`, `dspline()`, `dnonparam()`, `duser()`, `dmulti()`, `ditp()`
- `R/write.jags.R` — `mbnma.write()` generates JAGS model code as a character vector from a `dosefun` + model options
- `R/write.functions.R` — Lower-level JAGS code fragments, `get.prior()`, `replace.prior()`, `write.nma()`
- `R/predict.functions.R` — `predict.mbnma()`, `get.model.vals()`, placebo synthesis
- `R/rank.functions.R` — `rank.mbnma()`, `rank.nma()`, `calcprob()`, `sumrank()`
- `R/inconsistency.functions.R` — `nma.nodesplit()`, `inconsistency.loops()`
- `R/plot.functions.R` — `plot.mbnma()`, `plot.mbnma.network()`, `plot.mbnma.predict()`, `plot.mbnma.rank()`
- `R/mbnma-class.R`, `R/mbnma.network-class.R`, etc. — S3 methods (print, summary, plot) for each class

### JAGS model generation pipeline
`mbnma.run()` calls `mbnma.write()` → returns a character vector of JAGS model lines. The model is written to a `tempfile()` and passed to `R2jags::jags()`. If `model.file` is supplied, the auto-generated model is bypassed entirely.

The `dosefun` object's `jags` field contains the JAGS expression string for the dose-response (e.g. `"s.beta.1[agent[i,k]] * (dose[i,k] / (s.beta.2[agent[i,k]] + dose[i,k]))"`). `mbnma.write()` injects this into a standard NMA loop template.

### Dose-response parameter pooling
Each parameter in a `dosefun` has an `apool` value:
- `"rel"` — agent-specific relative effects (indexed `[agent[i,k]]`)
- `"common"` — single pooled value across all agents
- `"random"` — exchangeable across agents (adds `sd.*` parameter)
- numeric — fixed constant, not estimated

At least one parameter must use `"rel"`.

### Data format
Input data (`data.ab`) is long-format with one row per study arm. Required columns depend on likelihood: `{r, n}` (binomial), `{y, se}` (normal), `{r, E}` (Poisson). Column `dose=0` arms are automatically relabelled as Placebo. `studyID`, `agent`, `dose` are always required.

### Datasets included
`triptans` (binomial), `gout` (normal, disconnected agent), `osteopain` (normal, class structure), `ssri` (normal), `psoriasis75`, `psoriasis90`, `alog_pcfb` — used extensively in tests and vignettes.
