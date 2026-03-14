## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(facet)
library(dplyr)
library(tidyr)
library(broom)

# Load Brennan dataset
data(brennan)

# Conduct G-study with default (lme4) backend
# We specify a person-by-task design.
g <- gstudy(Score ~ (1 | Person) + (1 | Task), 
            data = brennan)

# View variance components
print(g)

## -----------------------------------------------------------------------------
# Summary with variance estimates (default)
summary(g)

# Summary with standard deviation estimates
summary(g, scale = "sd")

## -----------------------------------------------------------------------------
tidy(g)

## -----------------------------------------------------------------------------
glance(g)

## -----------------------------------------------------------------------------
# G-study with profile likelihood confidence intervals
g_profile <- gstudy(
  Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  ci_method = "profile"
)

summary(g_profile)

## -----------------------------------------------------------------------------
# G-study with bootstrap confidence intervals
# Using fewer iterations for speed in this example
g_boot <- gstudy(
  Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  ci_method = "boot",
  nsim = 50,  # Number of bootstrap samples
  boot.type = "perc"  # Percentile bootstrap
)

summary(g_boot)

## -----------------------------------------------------------------------------
# Method of moments G-study
g_mom <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + 
          (1 | Person:Task),
  data = brennan,
  backend = "mom"
)

# Summary shows confidence intervals
summary(g_mom)

## -----------------------------------------------------------------------------
# Example with nested data (Method of Moments)
# Using rajaratnam dataset: items nested within subtests
data(rajaratnam)

g_nested <- gstudy(
  Score ~ (1 | Person) + (1 | Subtest) + (1 | Item:Subtest) + (1 | Person:Subtest),
  data = rajaratnam,
  backend = "mom"
)

print(g_nested)

## ----eval=FALSE---------------------------------------------------------------
# # Bayesian G-study with brms
# # Note: Requires brms and Stan to be installed
# g_bayes <- gstudy(
#   Score ~ (1 | Person) + (1 | Task),
#   data = brennan,
#   backend = "brms"
# )
# 
# print(g_bayes)
# summary(g_bayes)

## -----------------------------------------------------------------------------
# Conduct D-study with 3 tasks and 4 raters
# (using our previously fitted lme4 g-study 'g')
d <- dstudy(g, n = list(Task = 3, Rater = 4))

print(d)
summary(d)

## -----------------------------------------------------------------------------
# Explore different numbers of tasks and raters
d_sweep <- dstudy(g, n = list(Task = c(3, 5, 10), Rater = c(2, 4, 8)))

# View coefficients for each design
# Using tidy() to get a nice data frame of coefficients
library(broom)
tidy(d_sweep)

## -----------------------------------------------------------------------------
# D-study with a cut-score of 7 (on a 0-10 scale)
# First, let's check the grand mean from our G-study
mu_y <- extract_grand_mean(g)
cat("Grand mean:", mu_y, "\n")

# Now conduct D-study with cut_score = 7
d_cut <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 7)

print(d_cut)

## -----------------------------------------------------------------------------
# Compare coefficients with and without cut-score
d_standard <- dstudy(g, n = list(Task = 3, Rater = 4))

# Standard phi (no cut-score)
cat("Standard phi:", d_standard$coefficients$phi, "\n")

# Phi-cut with cut-score of 7
cat("Phi-cut (c=7):", d_cut$coefficients$phi_cut, "\n")

# Phi-cut with different cut-scores
d_cut_5 <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 5)
d_cut_9 <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 9)

cat("Phi-cut (c=5):", d_cut_5$coefficients$phi_cut, "\n")
cat("Phi-cut (c=9):", d_cut_9$coefficients$phi_cut, "\n")

## -----------------------------------------------------------------------------
# Sweep with cut-score
d_sweep_cut <- dstudy(
  g, 
  n = list(Task = c(3, 5, 10), Rater = c(2, 4, 8)),
  cut_score = 7
)

# View results
tidy(d_sweep_cut)

## -----------------------------------------------------------------------------
# Access the grand mean used
d_cut$mu_y
d_cut$cut_score

## -----------------------------------------------------------------------------
# D-study with custom universe specification
# The object of measurement is Person
d_custom <- dstudy(g, n = list(Task = 3, Rater = 4), 
                   universe = c("Person", "Person:Task"))

## -----------------------------------------------------------------------------
# D-study with custom error specification
# Using only interaction terms as error sources
d_error <- dstudy(
  g,
  n = list(Task = 3, Rater = 4),
  error = c("Person:Task", "Person:Rater")
)

summary(d_error)

## -----------------------------------------------------------------------------
# D-study with aggregation over raters
# This models the case where we average scores across Tasks or Raters
d_agg <- dstudy(
  g,
  n = list(Task = 3, Rater = 4),
  aggregation = "Rater",
  residual_is = "Person:Task:Rater"
)

summary(d_agg)

## ----eval=FALSE---------------------------------------------------------------
# # D-study with posterior estimation (requires brms backend)
# # Note: This uses full posterior to compute coefficients with uncertainty
# d_post <- dstudy(
#   g_bayes,
#   n = list(rater = 3),
#   estimation = "posterior"
# )
# 
# print(d_post)
# 
# # Access posterior distributions
# names(d_post$posterior)

## ----eval=FALSE---------------------------------------------------------------
# # D-study with posterior from lme4 gstudy
# d_post2 <- dstudy(
#   g,
#   n = list(rater = 3),
#   estimation = "posterior"
# )

## ----eval=FALSE---------------------------------------------------------------
# # D-study with 95% credible intervals
# d_ci <- dstudy(
#   g_bayes,
#   n = list(Task = 3, Rater = 4),
#   ci = c("g", "phi")
# )
# 
# print(d_ci)

## ----eval=FALSE---------------------------------------------------------------
# # 90% credible interval
# d_ci_90 <- dstudy(
#   g_bayes,
#   n = list(Task = 3, Rater = 4),
#   ci = "g",
#   probs = c(0.05, 0.95)
# )

## ----eval=FALSE---------------------------------------------------------------
# d_cut_ci <- dstudy(
#   g_bayes,
#   n = list(Task = 3, Rater = 4),
#   cut_score = 7,
#   ci = c("g", "phi", "phi-cut")
# )

## ----eval=FALSE---------------------------------------------------------------
# # This will produce a warning
# g_lme4 <- gstudy(Score ~ (1|Person) + (1|Task), data = brennan)
# d <- dstudy(g_lme4, n = list(Task = 3), ci = "g")
# # Warning: Credible intervals for 'mom' and 'lme4' backends are not yet implemented.

## -----------------------------------------------------------------------------
# Bar plot of variance components (proportions)
plot(g, type = "proportion")

## -----------------------------------------------------------------------------
# Bar plot on variance scale
plot(g, type = "variance")

## -----------------------------------------------------------------------------
# Bar plot of G and D coefficients
plot(d, type = "coefficients")

## -----------------------------------------------------------------------------
# Sweep plot (for sample size exploration)
plot(d_sweep, type = "sweep")

## -----------------------------------------------------------------------------
# Plot all random effects
plot_ranef(g)

# Plot a specific facet
plot_ranef(g, which = "Person")

## -----------------------------------------------------------------------------
# Get variance components as data frame
vc <- tidy(g)
vc

## -----------------------------------------------------------------------------
# Filter specific components
vc %>%
  filter(component != "Residual") %>%
  mutate(pct_total = var / sum(var) * 100)

## -----------------------------------------------------------------------------
# Access the fitted model
model <- g$model
print(summary(model))

## ----eval=FALSE---------------------------------------------------------------
# # Extract variance-covariance matrices (requires loading lme4)
# library(lme4)
# vc_obj <- VarCorr(g)
# print(vc_obj)
# 
# # Extract random effects (BLUPs/conditional modes)
# re <- ranef(g)
# print(re)

## ----eval=FALSE---------------------------------------------------------------
# # Extract random effects from Bayesian model
# re_bayes <- ranef(g_bayes)
# print(re_bayes)
# 
# # Parameter pairs plot
# # pairs(g_bayes)

## ----eval=FALSE---------------------------------------------------------------
# # Extract posterior draws
# draws <- brms::as_draws_matrix(g_bayes$model)

## -----------------------------------------------------------------------------
# Extract latent scores from a G-study
# Using the default universe (Person only)
latent_scores <- latent(g)

head(latent_scores)

## -----------------------------------------------------------------------------
# Include Person:Task interaction in the universe
latent_with_interaction <- latent(g, universe = c("Person", "Person:Task"))

head(latent_with_interaction)

## -----------------------------------------------------------------------------
# Using formula syntax
latent_formula <- latent(g, universe = ~ Person + Person:Task)

head(latent_formula)

## ----eval=FALSE---------------------------------------------------------------
# # Create multivariate data
# brennan_wide <- brennan %>%
#   pivot_wider(names_from = Task, values_from = Score, names_prefix = "Task")
# 
# # Multivariate G-study
# g_mv <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = brennan_wide
# )
# 
# # Extract latent scores for each dimension
# latent_mv <- latent(g_mv)
# head(latent_mv)

## ----eval=FALSE---------------------------------------------------------------
# # Create multivariate data by spreading the tasks
# brennan_wide <- brennan %>%
#   pivot_wider(names_from = Task, values_from = Score, names_prefix = "Task")
# 
# # Multivariate G-study (automatically uses brms)
# g_mv <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = brennan_wide
# )
# 
# summary(g_mv)
# 
# # D-study for multivariate models
# d_mv <- dstudy(g_mv, n = list(Rater = 4))
# print(d_mv)

## ----eval=FALSE---------------------------------------------------------------
# data(rajaratnam)
# head(rajaratnam)

## ----eval=FALSE---------------------------------------------------------------
# library(brms)
# 
# # Long-format multivariate G-study
# # - Subtest is the dimension variable (3 levels: A, B, C)
# # - (0+Subtest|r|Person) - correlated random effects for Person
# # - (0+Subtest||Item) - uncorrelated random effects for Item
# # - sigma ~ 0 + Subtest - separate residual variance per subtest
# 
# g_long <- gstudy(
#   bf(Score ~ 0 + Subtest + (0+Subtest|r|Person) + (0+Subtest||Item),
#      sigma ~ 0 + Subtest),
#   data = rajaratnam,
#   backend = "brms"
# )
# 
# # Verify it was detected as long-format
# g_long$long_format_multivariate  # TRUE
# g_long$dimension_var             # "Subtest"

## ----eval=FALSE---------------------------------------------------------------
# g_long$sample_size_tibble

## ----eval=FALSE---------------------------------------------------------------
# print(g_long)

## ----eval=FALSE---------------------------------------------------------------
# g_long$correlations$random_effect_cor

## ----eval=FALSE---------------------------------------------------------------
# # D-study with user-specified sample sizes (applies to all dimensions)
# d_long <- dstudy(g_long, n = list(Person = 10))
# 
# # Coefficients are computed per dimension
# d_long$coefficients

## ----eval=FALSE---------------------------------------------------------------
# # The residual variances shown are on the original scale
# # (already exponentiated from log scale)
# g_long$variance_components %>%
#   filter(component == "Residual") %>%
#   select(dim, var, lower, upper)

## -----------------------------------------------------------------------------
# G-study with three facets (see Basic G-Study section)
# D-study: what if we used 5 tasks rated by 2 raters?
d_3f <- dstudy(
  g,
  n = list(Task = 5, Rater = 2)
)

summary(d_3f)

## -----------------------------------------------------------------------------
library(facet)
library(dplyr)

# 1. Load data
data(brennan)

# 2. Conduct G-study with profile confidence intervals
g <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
          (1 | Person:Task),
  data = brennan,
  ci_method = "profile"
)

# 3. View results
summary(g)

# 4. Conduct D-study with specific design (e.g., 5 tasks, 3 raters)
d <- dstudy(g, n = list(Task = 5, Rater = 3))

# 5. View coefficients
summary(d)

# 6. Explore alternative designs (Sample Size Sweep)
d_alt <- dstudy(g, n = list(Task = c(2, 5, 10), Rater = c(1, 2, 5)))

# 7. Plot sweep results
plot(d_alt, type = "sweep")

# 8. Extract and analyze variance components
tidy(g) %>%
  filter(component != "Residual") %>%
  arrange(desc(var))

## -----------------------------------------------------------------------------
is.gstudy(g)
is.dstudy(d)

## -----------------------------------------------------------------------------
# Display on SD scale
print(g, scale = "sd")
summary(g, scale = "sd")

