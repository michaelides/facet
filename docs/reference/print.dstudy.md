# Print Method for dstudy Objects

Print Method for dstudy Objects

## Usage

``` r
# S3 method for class 'dstudy'
print(x, digits = 3, scale = c("variance", "sd"), sem = FALSE, ...)
```

## Arguments

- x:

  A dstudy object.

- digits:

  Number of digits to display.

- scale:

  Scale for displaying results: "variance" (default) or "sd".

- sem:

  Logical; if TRUE, include standard errors of measurement (sem_rel and
  sem_abs) in the output. Default is FALSE.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns x.
