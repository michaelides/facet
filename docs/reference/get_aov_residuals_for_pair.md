# Get AOv Residuals for Pairwise Correlation

Extracts residuals from an aov model, handling the case where the model
was fit on a different subset of data than the pairwise data.

## Usage

``` r
get_aov_residuals_for_pair(
  aov_model,
  pairwise_data,
  resp,
  orig_data,
  resp_aov_formula = NULL
)
```

## Arguments

- aov_model:

  Fitted aov model

- pairwise_data:

  Data frame with pairwise complete cases

- resp:

  Response variable name

- orig_data:

  Original data used to fit the model

- resp_aov_formula:

  The original aov formula string (with Error terms) used to fit the
  model. When provided, the model is refit on pairwise data using this
  formula to preserve the Error() structure.

## Value

Vector of residuals aligned with pairwise_data
