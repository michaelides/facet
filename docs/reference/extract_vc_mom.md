# Extract Variance Components from Method of Moments Fit

Extracts variance components from a MoM fit with confidence intervals
using asymptotic approximations.

## Usage

``` r
extract_vc_mom(model, conf_level = 0.95, formula = NULL)
```

## Arguments

- model:

  A momfit object from
  [`fit_mom()`](https://github.com/yourorg/facet/reference/fit_mom.md).

- conf_level:

  Confidence level for intervals (default 0.95).

- formula:

  The original formula (used to extract response name). Optional.

## Value

A tibble with columns:

- component:

  Name of the variance component

- facet:

  Associated facet name

- type:

  Type: "main", "interaction", or "residual"

- var:

  Point estimate of variance

- pct:

  Percentage of total variance

- lower:

  Lower confidence interval bound

- upper:

  Upper confidence interval bound

- se:

  Standard error of the variance estimate

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  score = rnorm(100),
  person = factor(rep(1:20, each = 5)),
  item = factor(rep(1:5, times = 20))
)
model <- fit_mom(score ~ (1 | person) + (1 | item), data)
extract_vc_mom(model)
} # }
```
