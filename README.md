# lctmr

Latent Class Trajectory Modeling (LCTM) for pediatric growth data in R.

## Overview

**lctmr** identifies hidden subgroups (latent classes) of individuals with distinct growth patterns over time. It is built on the [`lcmm`](https://cran.r-project.org/package=lcmm) package's `hlme` fitter and implements the framework from Lennon et al. (2018) as an *investigator-driven* workflow: rather than letting automation pick everything, you inspect residual diagnostics first and decide the model's functional form, then hand the K-search and adequacy checks to the package.

For example, in pediatric growth data, LCTM might discover:

- Children with consistently high weight
- Children with rapid early growth that plateaus
- Children with slow but steady growth
- Children with declining growth trajectories

The package exposes a small, layered API: four workflow functions for the end-to-end process, an adequacy evaluator, two plotting functions, and four low-level metric/utility functions for advanced use.

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("evanxliu1/lctmr")
```

Requires R and the `lcmm` and `ggplot2` packages (plus `stats`/`grDevices`, which ship with R). `haven`, `splines`, `knitr`, `rmarkdown`, and `testthat` are suggested for reading data, splines, vignettes, and tests.

## Quick Start

```r
library(lctmr)
data(sample_growth)

# Step 1: Clean data (optional but recommended)
cleaned <- lctm_clean(
  sample_growth,
  outcome  = "weight_raw",
  time_var = "anthroage",
  id_var   = "childid"
)

# Step 2: Fit an initial model and review diagnostics
initial <- lctm_initial(cleaned, k = 3, degree = 2)
plot(initial)            # spaghetti, LOESS, residuals, and a guide panel

# Step 3: Refine — choose random effects from the residual patterns,
#         then let the package search K and Model A/B
result <- lctm_refine(initial, random = ~ 1 + anthroage, k_range = 2:5)

# Step 4: Interpret
summary(result)          # best model, BIC table, search history
plot(result)             # mean trajectories for the selected model
```

## The 4-Step Workflow

```
lctm_clean()        Prepare data, drop NAs, flag implausible values
      |
lctm_initial()      Fit initial model (random intercept only),
      |             emit residual diagnostics
      |
   [Investigator reviews residual patterns to choose random effects]
      |
lctm_refine()       Automated BIC search across K + Model A/B + adequacy
      |
plot() / summary()  Visualize and interpret the selected model
```

1. **Prepare data** (`lctm_clean`) — drop missing values, remove biologically implausible measurements using WHO/CDC z-score criteria, and enforce a minimum number of observations per subject.
2. **Fit initial model** (`lctm_initial`) — fit with a *random intercept only* so residual patterns reflect mean-trajectory misspecification rather than slope variability. Examine the plots to decide which random-effects terms to add.
3. **Refine** (`lctm_refine`) — sweep K, and at each K compare Model A (common variance) and Model B (proportional variance), stopping at the first model that meets all adequacy thresholds.
4. **Interpret** — mean-trajectory and spaghetti plots, class assignments, BIC table, and the full search history.

### Reading residual diagnostics to choose random effects

The diagnostic panel from `lctm_initial()` is the heart of the investigator-driven workflow. The residual pattern suggests which random-effects structure to pass to `lctm_refine(random = ...)`:

| Residual pattern | Suggested random effects |
|---|---|
| Flat horizontal band | `~ 1` (random intercept only) |
| Diagonal / linear trend | `~ 1 + time` |
| Curved (U or ∩) | `~ 1 + time + I(time^2)` |
| S-shaped / multiple bends | cubic terms or natural splines (`knots`) |

A common convention is to match the random-effects polynomial to the mean-trajectory `degree` term-for-term, so each subject gets their own coefficients on top of the class mean.

### How the trajectory is specified: `degree`, `knots`, `spline_degree`, and `random`

These four arguments are easy to confuse because they sound related, but they act on **different layers of the model**. For a subject *i* in latent class *k*:

```
outcome_ij  =  [ class-k MEAN trajectory ]  +  [ subject i's RANDOM deviation ]  +  error
```

| Argument | Layer | What it controls |
|---|---|---|
| `degree` | class mean (fixed) | Polynomial degree of the mean trajectory **when not using splines**: 1 = linear, 2 = quadratic, 3 = cubic. When `knots` are supplied it is **inert** — the spline defines the mean — and is automatically set to 1. |
| `knots` | class mean (fixed) | Switches the mean trajectory from a global polynomial to a **piecewise spline** with bends at the given time values (e.g. `c(6, 12)`). Supplying it forces `degree = 1`. |
| `spline_degree` | class mean (fixed) | Order of the spline *pieces* when `knots` are set: `3` = natural cubic spline (`splines::ns()`, C² smooth), `2` = quadratic B-spline (`splines::bs(degree = 2)`, C¹, fewer parameters per knot). Ignored without `knots`. |
| `random` | per-subject (random) | How each individual **deviates from their class mean over time**. Contains the time variable and its powers only, e.g. `~ 1 + age + I(age^2)`. Independent of the mean's form — see below. |

Two consequences that trip people up:

- **`degree` and `random` are different things.** `degree = 1, random = ~ 1 + age + I(age^2)` is *not* contradictory: `degree` describes the class *mean* (and is inert once you use knots), while the `random` formula describes per-subject *deviations*. A quadratic random structure on top of a spline mean is perfectly valid.
- **The random part doesn't have to match the fixed part.** Even with a spline mean trajectory, you typically keep `random` as a low-order polynomial (`~ 1 + age + I(age^2)`). Putting the full spline basis in `random` means estimating a covariance over every basis coefficient *per subject*, which is heavy and usually won't converge. A quadratic random effect is the stable, conventional way to let individuals curve around the class mean.

> ⚠️ **WARNING — these two arguments fail quietly if misused:**
>
> - **`spline_degree` does nothing without `knots`.** Setting `spline_degree` (to either `2` or `3`) but forgetting `knots` does **not** give you a spline — you get a plain polynomial of whatever `degree` is. The package emits a warning whenever you pass `spline_degree` without `knots`, but the model still fits, so read your warnings. **Always pass `knots` when you want a spline.**
> - **`degree` is overridden to 1 whenever `knots` are present.** If you supply `knots`, any `degree` you set (or that is inherited from the initial model) is silently replaced by 1, with a message — the spline defines the shape. If you actually wanted a polynomial, *remove* `knots`.

Baseline covariates (e.g. socioeconomic status, prematurity) are a *fifth*, separate thing — they go in neither `degree` nor `random` but in the `covariates` argument of `lctm_refine()`, which adds them to the fixed formula only (shifting the overall mean, the same across classes; not the class-specific shape).

### Splines

`lctm_initial()` and `lctm_refine()` accept a `knots` argument for splines, useful when residual diagnostics show non-polynomial patterns (e.g., distinct early vs. late growth phases). Supplying `knots` replaces the polynomial trajectory: `degree` is automatically set to 1 (with a message) because the spline carries the trajectory shape.

```r
init   <- lctm_initial(cleaned, k = 3, knots = c(6, 12))
result <- lctm_refine(init, knots = c(6, 12), k_range = 2:5)
```

The `spline_degree` argument controls the polynomial order of the spline pieces:

- `spline_degree = 3` (default) — **natural cubic spline** via `splines::ns()`. Piecewise cubic with continuous curvature (C²) and linear tails beyond the boundary knots. The conventional choice.
- `spline_degree = 2` — **quadratic B-spline** via `splines::bs(degree = 2)`. Piecewise quadratic with continuous slope but discontinuous curvature (C¹). It uses one fewer parameter per knot than cubic, so it is more parsimonious and tends to estimate more stably in a multi-class mixture — useful when a cubic spline is over-flexible or fails to converge.

```r
# Quadratic (piecewise) spline instead of the default cubic
init   <- lctm_initial(cleaned, k = 3, knots = c(6, 12), spline_degree = 2)
result <- lctm_refine(init, knots = c(6, 12), spline_degree = 2, k_range = 2:5)
```

`spline_degree` is inherited by `lctm_refine()` from the initial model unless overridden, and is ignored when `knots` is NULL.

## Functions

| Layer | Function | Purpose |
|-------|----------|---------|
| Workflow | `lctm_clean()` | Validate input, drop NAs, apply WHO/CDC outlier cutoffs, enforce min observations per subject. Optional — raw data can go straight into `lctm_initial()`. |
| Workflow | `lctm_initial()` | Fit a single-K model with random intercept only; produce spaghetti/LOESS/residual plots and a guide panel. |
| Workflow | `lctm_refine()` | Automated refinement: BIC sweep across K, Model A/B comparison at each K, adequacy gating. |
| Adequacy | `lctm_adequacy()` | Evaluate APPA, OCC, relative entropy, and the minimum class proportion for any fitted model; used internally by `lctm_refine()` but exposed for standalone use. |
| Adequacy | `lctm_filter_small_classes()` | Investigator-driven companion to the `min_prop` floor: remove subjects assigned to below-floor classes and return a new dataset (plus the removed IDs) to refit on. |
| Visualization | `lctm_plot_trajectories()` | Mean and/or spaghetti trajectory plots (legend shows n and %). |
| Visualization | `lctm_plot_residuals()` | Residual diagnostic plots, optionally faceted by class. |
| Low-level | `calc_appa()` | Average Posterior Probability of Assignment, per class. |
| Low-level | `calc_occ()` | Odds of Correct Classification, per class. |
| Low-level | `calc_relative_entropy()` | Relative (normalized) entropy of the posterior probabilities. |
| Utility | `assign_classes()` | Hard class assignment from a posterior probability matrix. |

### Key arguments and defaults

**`lctm_clean(data, outcome, time_var, id_var, sex_var = NULL, type = NULL, standards = "WHO", zscore_col = NULL, zscore_col_cdc = NULL, birth_weight_col = NULL, check_decrease = FALSE, verbose = TRUE)`**

- `type` — outcome type for outlier detection: `"weight"`, `"height"`, `"hc"`, or `NULL` (no z-score filtering). When set, `zscore_col` is required.
- `standards` — `"WHO"` (default), `"CDC"`, or `"both"` (then `zscore_col_cdc` is required).
- `birth_weight_col` — birth weight in grams; rows with ≤ 500 g or ≥ 5000 g are dropped.
- `check_decrease` — for `"height"`/`"hc"`, flags implausible decreases > 3 cm between visits.
- Returns an `lctm_cleaned` object (cleaned `data` plus the column names and counts of rows/subjects removed).

**`lctm_initial(data, outcome = NULL, time_var = NULL, id_var = NULL, k = 2, degree = 2, knots = NULL, spline_degree = NULL, sex_var = NULL, save_pdf = NULL, verbose = TRUE)`**

- Accepts an `lctm_cleaned` object *or* a raw data frame (in which case `outcome`/`time_var`/`id_var` are required).
- `degree` — 1 = linear, 2 = quadratic, 3 = cubic; mutually exclusive with `knots`.
- `knots` / `spline_degree` — knot positions and the spline piece degree (3 = natural cubic via `ns()`, 2 = quadratic via `bs()`); see [Splines](#splines).
- `save_pdf` — path to write the diagnostic panel.
- Returns an `lctm_initial` object holding the fitted random-intercept model, a base (single-class) model for start values, and a named list of `plots`.

**`lctm_refine(initial, random = NULL, k_range = 2:7, degree = NULL, knots = NULL, spline_degree = NULL, covariates = NULL, models = c("A", "B"), adequacy_thresholds = list(appa = 0.70, occ = 5.0, entropy = 0.5, min_prop = 0.05), start_simple = FALSE, save_pdf = NULL, verbose = TRUE)`**

- `random` — random-effects formula chosen from the diagnostics; defaults to the full polynomial matching `degree`.
- `degree`/`knots`/`spline_degree` — inherited from `initial` if `NULL`.
- `covariates` — names added to the *fixed* part of the formula only.
- `models` — which parameterizations to try and in what order (see below).
- `start_simple` — use a single-class model for starting values, which speeds up large datasets.
- Returns an `lctm_result` with `best_model`, `best_k`, `best_model_type`, the `adequacy` object, the full `bic_table`, `all_models`, `class_assignments`, and the `search_history`.

### Model A vs. Model B

At each K, `lctm_refine()` can fit two residual-variance parameterizations:

- **Model A** — common variance across classes (`hlme(nwg = FALSE)`).
- **Model B** — proportional (class-specific) variance (`hlme(nwg = TRUE)`).

The default `models = c("A", "B")` tries the simpler shared-variance model first and escalates to proportional variance only if needed. Reverse with `models = c("B", "A")`, or restrict to one.

## Objects and S3 methods

Objects flow from one stage to the next:

```
data.frame ── lctm_clean() ──▶ lctm_cleaned
                                    │
                                    ▼
                              lctm_initial() ──▶ lctm_initial   (+ diagnostic plots)
                                    │
                          [choose random effects]
                                    │
                                    ▼
                              lctm_refine()  ──▶ lctm_result    (best lctm_model + adequacy)
                                                      │
                              lctm_adequacy() ◀───────┘         lctm_adequacy (APPA / OCC / entropy)
```

The package defines `print` methods for `lctm_cleaned`, `lctm_initial`, `lctm_model`, `lctm_result`, and `lctm_adequacy`; `plot` methods for `lctm_initial`, `lctm_model`, and `lctm_result`; and `summary`/`coef` methods for the model and result objects. So `print(result)`, `summary(result)`, and `plot(result)` all work directly on the output of `lctm_refine()`.

- `plot(initial, which = NULL)` shows all diagnostic panels, or one of `"spaghetti"`, `"loess"`, `"residuals"`, `"guide"`.
- `plot(result, type = "mean")` delegates to `lctm_plot_trajectories()`; pass `type = "spaghetti"` or `"both"`.

## Model Adequacy Metrics

| Metric | Threshold | Description |
|--------|-----------|-------------|
| APPA | ≥ 0.70 | Average Posterior Probability of Assignment — mean confidence of assignment within each class. |
| OCC | ≥ 5.0 | Odds of Correct Classification — assignment certainty relative to class prevalence (`Inf` = perfect). |
| Entropy | ≥ 0.5 | Relative (normalized) entropy — overall class separation, 0 to 1. |
| Min. class proportion | ≥ 0.05 | Smallest class as a share of subjects. **Reported only** — see below. |

A model **passes** (`overall_pass = TRUE`) when APPA, OCC, and entropy all meet their thresholds across non-empty classes *and* the model is not degenerate (no empty classes). Thresholds are configurable via `lctm_refine(adequacy_thresholds = ...)` or `lctm_adequacy(thresholds = ...)`. The returned `lctm_adequacy` object reports per-class values, per-metric pass flags, an `overall_pass`, an `is_degenerate` flag, the `class_proportions`, and a `min_prop_pass` flag.

### Minimum class proportion and the filter/refit loop

APPA, OCC, and entropy all reward *confident* assignment, which a model can achieve by isolating a handful of unusual subjects into their own tiny class. The `min_prop` floor (default 0.05) flags any class smaller than that share of subjects.

It is reported via `class_proportions` and `min_prop_pass` but is **deliberately not folded into `overall_pass`**: removing subjects is an investigator decision, not an automatic one. To act on the flag, use `lctm_filter_small_classes()`, which never modifies the original data and never refits — it just hands back a filtered dataset and the full list of removed IDs (so exclusions can be reported in, e.g., a participant-flow diagram):

```r
result <- lctm_refine(initial, random = ~ 1 + anthroage, k_range = 2:5)
print(result$adequacy)                      # see which classes fall below the floor

filt   <- lctm_filter_small_classes(result, min_prop = 0.05)  # remove them, get new data
init2  <- lctm_initial(filt$data, k = result$best_k, degree = 2)
result2 <- lctm_refine(init2, random = ~ 1 + anthroage,
                       k_range = result$best_k)                # refit (same K first is usual)
```

Repeat as many rounds as you judge appropriate, or stop.

## Bundled data

`data(sample_growth)` loads an anonymized longitudinal growth dataset (935 observations of 290 children) with columns:

| Column | Description |
|---|---|
| `childid` | Child identifier |
| `timepoint` | Measurement timepoint label |
| `anthroage` | Age in months at measurement |
| `weight_raw` | Raw weight (kg) |
| `waz` | Weight-for-age z-score (WHO standard) |
| `month` | Month of measurement |

## Vignettes

Two vignettes ship with the package:

- `vignette("introduction", package = "lctmr")` — what LCTM is, when to use it, and the difference between the investigator-driven and automated workflows.
- `vignette("workflow", package = "lctmr")` — an end-to-end walkthrough on `sample_growth`, including the residual-pattern interpretation guide.

Build them locally with `devtools::build_vignettes()` or read them on the GitHub repo.

## Dependencies

- [lcmm](https://cran.r-project.org/package=lcmm) — core latent class mixed model fitting (`hlme`)
- [ggplot2](https://cran.r-project.org/package=ggplot2) — visualization

## References

- **Lennon et al. (2018)** — the LCTM framework this package adapts.
- **Gee (2014)** — APPA / OCC adequacy criteria.
- **Proust-Lima et al.** — the `lcmm` package and its `hlme` fitter.

## License

MIT
