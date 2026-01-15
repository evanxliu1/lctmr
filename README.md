# lctmr

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/evanliu/lctmr)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Latent Class Trajectory Modeling (LCTM) for pediatric growth data in R.

## Overview

**lctmr** identifies hidden subgroups (latent classes) of individuals with distinct growth patterns over time. For example, in pediatric growth data, LCTM might discover:

- Children with consistently high weight
- Children with rapid early growth that plateaus
- Children with slow but steady growth
- Children with declining growth trajectories

The package provides three levels of abstraction:

| Level | Function | Use Case |
|-------|----------|----------|
| High | `lctm_auto()` | Fully automated analysis |
| Mid | `lctm_setup()` / `lctm_fit()` / `lctm_adequacy()` | Step-by-step with control |
| Low | `calc_appa()` / `calc_occ()` / `calc_relative_entropy()` | Advanced customization |

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("evanliu/lctmr")
```

## Quick Start

```r
library(lctmr)

# Load sample data
data(sample_growth)

# Option A: Fully automated (recommended for most users)
result <- lctm_auto(
  sample_growth,
  outcome = "weight_raw",
  time_var = "anthroage",
  id_var = "childid"
)
summary(result)
plot(result)

# Option B: Step-by-step (for iteration and customization)
setup <- lctm_setup(sample_growth, "weight_raw", "anthroage", "childid")
print(setup$bic_table)  # See ranked K values

model <- lctm_fit(setup, k = 4, model = "F")
lctm_adequacy(model)    # Check if model passes adequacy criteria

# If adequacy fails, try different K
model <- lctm_fit(setup, k = 3, model = "F")
lctm_adequacy(model)

# Visualize results
lctm_plot_trajectories(model)
lctm_plot_trajectories(model, type = "spaghetti")
```

## The 5-Step Workflow

1. **Determine random effects structure** - Examine residual plots to choose model complexity
2. **Compare BIC across K values** - Fit models with K=1 to K=7 classes, rank by BIC
3. **Fit candidate models** - Test Models E (common variance) and F (proportional variance)
4. **Check adequacy** - Verify APPA > 0.70, OCC > 5.0, Entropy > 0.5
5. **Visualize** - Create mean trajectory and spaghetti plots

See `vignette("workflow")` for detailed guidance.

## Model Adequacy Metrics

| Metric | Threshold | Description |
|--------|-----------|-------------|
| APPA | > 0.70 | Average Posterior Probability of Assignment |
| OCC | > 5.0 | Odds of Correct Classification |
| Entropy | > 0.5 | Relative Entropy (class separation) |

## Dependencies

- [lcmm](https://cran.r-project.org/package=lcmm) - Core latent class mixed model fitting
- [ggplot2](https://cran.r-project.org/package=ggplot2) - Visualization


## License

MIT
