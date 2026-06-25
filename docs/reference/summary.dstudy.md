# Summary Method for dstudy Objects

Summary Method for dstudy Objects

## Usage

``` r
# S3 method for class 'dstudy'
summary(object, scale = c("variance", "sd"), digits = 4, sem = FALSE, ...)
```

## Arguments

- object:

  A dstudy object.

- scale:

  Scale for displaying results: "variance" (default) or "sd".

- digits:

  Number of digits to display (default 4).

- sem:

  Logical; if TRUE, include standard errors of measurement (sem_rel and
  sem_abs) in the output. Default is FALSE.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns object.
