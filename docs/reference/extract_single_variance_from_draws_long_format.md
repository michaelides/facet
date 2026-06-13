# Extract Single Variance Component from Long-Format brms Model

Helper function to extract a single variance component from posterior
draws for long-format multivariate models.

## Usage

``` r
extract_single_variance_from_draws_long_format(
  draws,
  grp,
  resp,
  type,
  random_summary,
  spec_pars,
  is_mv,
  is_log_link = FALSE
)
```

## Arguments

- draws:

  Posterior draws matrix.

- grp:

  Group name (facet name or "residual\_\_").

- resp:

  Response dimension name.

- type:

  Component type: "main", "interaction", or "residual".

- random_summary:

  Random effects summary from brms.

- spec_pars:

  Special parameters summary from brms.

- is_mv:

  Whether this is a multivariate model.

- is_log_link:

  Whether sigma uses log link.

## Value

A data.frame with variance component estimates.
