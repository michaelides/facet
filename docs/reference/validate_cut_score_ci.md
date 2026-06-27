# Validate Cut Score and Credible Interval Parameters

Extracts grand mean if cut_score provided, validates ci parameter
against estimator, validates probs, and warns if phi-cut CI requested
without cut_score.

## Usage

``` r
validate_cut_score_ci(
  gstudy_obj,
  cut_score = NULL,
  ci = NULL,
  probs = c(0.025, 0.975)
)
```

## Arguments

- gstudy_obj:

  A gstudy or mgstudy object

- cut_score:

  Optional cutoff score

- ci:

  Credible interval specification

- probs:

  Probability levels for credible intervals

## Value

A list with mu_y, ci (possibly modified), probs
