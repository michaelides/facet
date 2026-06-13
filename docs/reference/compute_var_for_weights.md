# Compute VAR for Given Weights

Computes the Value Added Ratio for all subscales given a weight vector.
Used as the objective function for subscale optimization.

## Usage

``` r
compute_var_for_weights(weights, dstudy_obj, type = "rel")
```

## Arguments

- weights:

  Numeric vector of weights (will be normalized)

- dstudy_obj:

  A dstudy object

- type:

  "rel" for var_rel or "abs" for var_abs

## Value

Named vector of VAR values for each subscale
