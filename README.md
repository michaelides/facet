# facet

<!-- badges: start -->
[![R-CMD-check](https://github.com/michaelides/facet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/michaelides/facet/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/facet)](https://CRAN.R-project.org/package=facet)
<!-- badges: end -->

Modern Generalizability Theory with Variance Component Models

## Overview

**facet** provides a comprehensive framework for conducting Generalizability Theory (G-theory) analyses using variance component models. It provides a bridge between classical test theory and standard linear mixed-effects models, offering:

- **Multiple backends**: **lme4** (frequentist), **brms** (Bayesian), and **mom** (Method of Moments)
- **Univariate and multivariate analyses** using `mvbind()`
- **Automatic backend selection** based on your design
- **Built-in datasets**: Includes classic G-theory datasets like `brennan` and `rajaratnam`
- **Rich Visualization**: Built-in support for variance plots, D-study sweeps, and random effects caterpillar plots

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("michaelides/facet")
```

## Usage

### G-Study

A G-study estimates variance components for each facet in your measurement design. The **facet** package uses standard **lme4** formula syntax:

```r
library(facet)

# Load included Brennan dataset (10 persons, 3 tasks, 12 raters)
data(brennan)

# Conduct a G-study with Persons, Tasks, and Raters
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + 
            (1 | Person:Task) + (1 | Person:Rater), 
            data = brennan)

# View variance components
summary(g)
```

### D-Study

A D-study computes reliability coefficients for specific measurement designs:

```r
# Conduct a D-study with 5 tasks and 2 raters
d <- dstudy(g, n = list(Task = 5, Rater = 2))

# View Generalizability (G) and Dependability (D) coefficients
summary(d)
```

### Exploring Sample Sizes (Sweep)

Rapidly evaluate reliability across multiple sample size combinations:

```r
# Explore different numbers of Tasks and Raters
d_sweep <- dstudy(g, n = list(Task = c(2, 5, 10), Rater = c(1, 2, 5)))

# Plot the results
plot(d_sweep, type = "sweep")
```

### Visualization

Visualize random effects and uncertainty using caterpillar plots:

```r
# Plot random effects for Persons
plot_ranef(g, which = "Person")
```

## Features

### Analysis Capabilities
- **Automatic backend selection** — intelligently chooses between lme4, brms, and method of moments based on your formula structure
- **Univariate and multivariate G-studies** — analyze multiple outcomes simultaneously using `mvbind()`
- **Support for Complex Designs** — full support for crossed, nested, and mixed designs
- **Sample size exploration** — rapidly evaluate reliability across multiple designs using simple list-based input
- **Advanced D-study options** — support for custom universe scores, error specifications, and score aggregation

### Variance Component Analysis
- **Confidence/credible intervals** — all backends provide interval estimates: brms offers Bayesian credible intervals, lme4/mom provide bootstrap or analytical intervals
- **Flexible display scales** — view variance components as variances or standard deviations using the `scale` argument
- **Proportion of variance** — each component expressed as percentage of total variance for easy interpretation

### Bayesian Features (brms backend)
- **Full Posterior Control** — full access to MCMC samples and convergence diagnostics (Rhat, ESS)
- **Custom priors** — specify informative priors on variance components
- **Multivariate Models** — estimated residual correlations and cross-component variances

### Output and Integration
- **Tidyverse Compatible** — fully compatible with `broom` methods like `tidy()` and `glance()`
- **Model inspection** — easily access underlying model objects (`g$model`)
- **Rich Visualization** — built-in `ggplot2`-based plotting for G-studies, D-studies, and random effects

## Backends

### lme4 (Default)
- Fast frequentist estimation for univariate designs
- Best for routine analysis and large datasets

### brms (Bayesian)
- Bayesian estimation using Stan
- Required for multivariate models and custom priors
- Provides full posterior distributions

### Method of Moments (MOM)
- Fast closed-form solution for balanced designs
- Analytical confidence intervals without bootstrapping
- Useful for large balanced datasets where lme4 might be slow

## Output Structure

### G-Study Object

```r
# Variance components
tidy(g)

# Summary statistics
glance(g)

# Plot
plot(g)
```

### D-Study Object

```r
# Coefficients
tidy(d)

# Summary
summary(d)

# Plot
plot(d)
```

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.

Shavelson, R. J., & Webb, N. M. (1991). *Generalizability theory: A primer*. Sage Publications.

## License

MIT
