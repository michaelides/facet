# Caterpillar Plot for Random Effects

Creates a caterpillar (dot plot with error bars) to display random
effects from a gstudy object. Works with both lme4 and brms estimators,
providing a unified interface regardless of estimator used.

## Usage

``` r
plot_ranef(
  x,
  which = NULL,
  ci_level = 0.95,
  sort = TRUE,
  colors = NULL,
  ncol = NULL,
  ...
)
```

## Arguments

- x:

  A gstudy object.

- which:

  Which random effect(s) to plot. If NULL (default), plots all random
  effects in a multi-panel plot. Can be a character vector of specific
  facet names, or a single facet name.

- ci_level:

  Confidence/credible interval level. Default is 0.95.

- sort:

  Logical; if TRUE (default), sort effects by estimate size.

- colors:

  Optional vector of two colors for points and error bars.

- ncol:

  Number of columns for multi-panel plot when which is NULL. Default is
  NULL (automatic layout).

- ...:

  Additional arguments (passed to ggplot2 theme functions).

## Value

A ggplot object (invisibly).

## Examples

``` r
# Fit a G-study with lme4 (default) using the brennan dataset
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan
)

# Plot all random effects
plot_ranef(g)


# Plot only one specific random effect
plot_ranef(g, which = "Person")


# \donttest{
# Fit G-study with brms for Bayesian credible intervals
g_brms <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan,
  estimator = "brms",
  iter = 2000, cores = 4, refresh = 1000
)
#> Compiling Stan program...
#> Start sampling
#> Warning: There were 90 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# Plot with brms (uses posterior credible intervals)
plot_ranef(g_brms)

# }
```
