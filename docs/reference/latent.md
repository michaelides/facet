# Extract Latent Scores from G-Study Objects

Creates latent scores (universe score estimates) for the object of
measurement based on random effects extracted from a gstudy or mgstudy
object. The latent scores are formed on the basis of the universe
defined in the gstudy/mgstudy, or a user-specified universe.

## Usage

``` r
latent(x, universe = NULL, ...)

# Default S3 method
latent(x, universe = NULL, ...)

# S3 method for class 'gstudy'
latent(x, universe = NULL, ...)

# S3 method for class 'mgstudy'
latent(x, universe = NULL, ...)
```

## Arguments

- x:

  A gstudy or mgstudy object.

- universe:

  Specification for components that contribute to the universe score.
  Can be:

  - NULL (default): universe includes only the object of measurement

  - A character vector: c("Person", "Person:Rater")

  - A formula: ~ Person + Person:Rater

  The object of measurement is always included. Interaction terms in the
  universe create separate variables for each level combination.

- ...:

  Additional arguments (currently unused).

## Value

A data frame with:

- A column for the object of measurement levels

- A `latent` column (univariate) or `latent_<dim>` columns
  (multivariate)

- Additional columns for interaction combinations if specified in
  universe

## Details

### Universe Specification

By default, the latent score is based only on the object of
measurement's random effect (e.g., Person). You can expand the universe
to include interactions with the object (e.g., Person:Rater).

When interaction terms are included in the universe, separate latent
score variables are created for each level of the non-object facet. For
example, if universe includes "Person:Rater" with raters A, B, C, the
output will have columns `latent_RaterA`, `latent_RaterB`,
`latent_RaterC` in addition to the base `latent` column.

### Multivariate Models

For mgstudy objects, latent scores are computed separately for each
dimension and returned in wide format with separate columns per
dimension.

## See also

[`gstudy()`](https://github.com/yourorg/facet/reference/gstudy.md),
[`ranef()`](https://github.com/yourorg/facet/reference/ranef.md)

## Examples

``` r
# Basic usage - object of measurement only
g <- gstudy(Score ~ (1 | Person) + (1 | Rater), data = brennan)
latent_scores <- latent(g)
head(latent_scores)
#>   Person        latent
#> 1      1  3.563932e-15
#> 2      2  7.292988e-01
#> 3      3 -6.077490e-01
#> 4      4 -6.077490e-01
#> 5      5  6.685239e-01
#> 6      6 -9.116235e-01

# With interaction in universe
if (FALSE) { # \dontrun{
latent_scores <- latent(g, universe = c("Person", "Person:Rater"))
} # }
```
