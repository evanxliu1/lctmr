# lctmr

Latent Class Trajectory Modeling (LCTM) for pediatric growth data in R.

## Overview

**lctmr** identifies hidden subgroups (latent classes) of individuals with distinct growth patterns over time. Built on the `lcmm` package, it implements the framework from Lennon et al. (2018) with an investigator-driven workflow that places researcher judgment at the center of model specification.

For example, in pediatric growth data, LCTM might discover:

- Children with consistently high weight
- Children with rapid early growth that plateaus
- Children with slow but steady growth
- Children with declining growth trajectories

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("evanxliu1/lctmr")
```

## Quick Start

```r
library(lctmr)
data(sample_growth)

# Step 1: Clean data
cleaned <- lctm_clean(sample_growth, "weight_raw", "anthroage", "childid")

# Step 2: Initial model + diagnostic plots
initial <- lctm_initial(cleaned, k = 3, degree = 2)
plot(initial)  # Review residual patterns

# Step 3: Refine (choose random effects based on residual patterns)
result <- lctm_refine(initial, random = ~ 1 + anthroage, k_range = 2:5)

# Step 4: Interpret
summary(result)
plot(result)
```

## The 4-Step Workflow

```
lctm_clean()        Prepare data, flag implausible values
      |
lctm_initial()      Fit initial model, examine residual diagnostics
      |
   [Investigator reviews residual patterns to decide random effects]
      |
lctm_refine()       Automated BIC search across K + Model A/B + adequacy
      |
lctm_plot_*()       Visualize and interpret results
```

1. **Prepare data** (`lctm_clean`) - Remove biologically implausible values using WHO/CDC criteria
2. **Fit initial model** (`lctm_initial`) - Fit without random effects; examine residual patterns to guide model specification
3. **Refine** (`lctm_refine`) - Automates BIC comparison across K values, tests Model A (common variance) and Model B (proportional variance), checks adequacy criteria
4. **Interpret results** - Mean trajectory plots, spaghetti plots, class assignments

### Automated Alternative

For quick exploratory analyses, `lctm_auto()` handles all decisions automatically:

```r
result <- lctm_auto(
  sample_growth,
  outcome = "weight_raw",
  time_var = "anthroage",
  id_var = "childid"
)
```

## Functions

| Layer | Function | Purpose |
|-------|----------|---------|
| Workflow | `lctm_clean()` | Data preparation and outlier detection |
| Workflow | `lctm_initial()` | Initial model with diagnostic plots and residual guide |
| Workflow | `lctm_refine()` | Automated refinement: BIC search + Model A/B + adequacy |
| Automated | `lctm_auto()` | Fully automated analysis |
| Legacy | `lctm_setup()` / `lctm_fit()` | Step-by-step model fitting |
| Adequacy | `lctm_adequacy()` | Evaluate APPA, OCC, and relative entropy |
| Low-level | `calc_appa()` / `calc_occ()` / `calc_relative_entropy()` | Individual adequacy metrics |
| Utility | `assign_classes()` / `compare_models()` | Helper functions |
| Visualization | `lctm_plot_trajectories()` | Mean and spaghetti trajectory plots |
| Visualization | `lctm_plot_residuals()` | Residual diagnostic plots |

## Model Adequacy Metrics

| Metric | Threshold | Description |
|--------|-----------|-------------|
| APPA | > 0.70 | Average Posterior Probability of Assignment |
| OCC | > 5.0 | Odds of Correct Classification |
| Entropy | > 0.5 | Relative Entropy (class separation) |

A model passes adequacy when all three metrics meet their thresholds across all classes.

## Dependencies

- [lcmm](https://cran.r-project.org/package=lcmm) - Core latent class mixed model fitting
- [ggplot2](https://cran.r-project.org/package=ggplot2) - Visualization

## License

MIT
