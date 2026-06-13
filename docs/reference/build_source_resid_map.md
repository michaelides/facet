# Build Residual MS/DF Mapping for Each Source

Builds a mapping from each ANOVA source to the appropriate residual MS
and df, using stratum information to select the correct residual when
multiple strata contain "Residuals" entries (common in nested and
crossed designs with Error() terms).

## Usage

``` r
build_source_resid_map(anova_results)
```
