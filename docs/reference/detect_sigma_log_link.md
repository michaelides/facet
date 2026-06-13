# Detect if Sigma Uses Log Link in brms Model

Checks whether the sigma parameters in a brms model use a log link
(indicated by 'b_sigma\_' prefix in parameter names).

## Usage

``` r
detect_sigma_log_link(model)
```

## Arguments

- model:

  A brmsfit object.

## Value

Logical; TRUE if sigma uses log link, FALSE otherwise.
