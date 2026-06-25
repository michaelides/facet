# Caterpillar Plot for Random Effects

Creates a caterpillar (dot plot with error bars) to display random
effects from a gstudy object. Works with both lme4 and brms backends,
providing a unified interface regardless of backend used.

## Usage

``` r
plot_ranef(
  x,
  which = NULL,
  ci_level = 0.95,
  sort = TRUE,
  colors = NULL,
  ncol = NULL,
  ...
)
```

## Arguments

- x:

  A gstudy object.

- which:

  Which random effect(s) to plot. If NULL (default), plots all random
  effects in a multi-panel plot. Can be a character vector of specific
  facet names, or a single facet name.

- ci_level:

  Confidence/credible interval level. Default is 0.95.

- sort:

  Logical; if TRUE (default), sort effects by estimate size.

- colors:

  Optional vector of two colors for points and error bars.

- ncol:

  Number of columns for multi-panel plot when which is NULL. Default is
  NULL (automatic layout).

- ...:

  Additional arguments (passed to ggplot2 theme functions).

## Value

A ggplot object (invisibly).

## Examples

``` r
# Fit a G-study with lme4 (default) using the brennan dataset
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan
)

# Plot all random effects
plot_ranef(g)


# Plot only one specific random effect
plot_ranef(g, which = "Person")


# \donttest{
# Fit G-study with brms for Bayesian credible intervals
g_brms <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
  (1 | Person:Task),
  data = brennan,
  backend = "brms",
  iter = 2000, cores = 4, refresh = 1000
)
#> Compiling Stan program...
#> Trying to compile a simple C file
#> Running /Library/Frameworks/R.framework/Resources/bin/R CMD SHLIB foo.c
#> using C compiler: ‘Apple clang version 21.0.0 (clang-2100.1.1.101)’
#> using SDK: ‘MacOSX26.5.sdk’
#> clang -arch arm64 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I"/Users/afp18axu/Library/R/arm64/4.6/library/Rcpp/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/unsupported"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/BH/include" -I"/Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/src/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/RcppParallel/include/"  -I"/Users/afp18axu/Library/R/arm64/4.6/library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c foo.c -o foo.o
#> In file included from <built-in>:1:
#> In file included from /Users/afp18axu/Library/R/arm64/4.6/library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22:
#> In file included from /Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/Eigen/Dense:1:
#> In file included from /Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/Eigen/Core:19:
#> /Users/afp18axu/Library/R/arm64/4.6/library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: 'cmath' file not found
#>   679 | #include <cmath>
#>       |          ^~~~~~~
#> 1 error generated.
#> make: *** [foo.o] Error 1
#> Start sampling
#> Warning: There were 27 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems

# Plot with brms (uses posterior credible intervals)
plot_ranef(g_brms)

# }
```
