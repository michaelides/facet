## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center"
)

## -----------------------------------------------------------------------------
library(facet)
library(broom)
library(dplyr)

# Brennan (2001) dataset: Person crossed with Task and Rater
data(brennan)
head(brennan)

# Rajaratnam, Cronbach, & Gleser (1965): Person crossed with Subtest, 
# Items nested within Subtest
data(rajaratnam)
head(rajaratnam)

## -----------------------------------------------------------------------------
# Simple G-study with default REML backend
g <- gstudy(Score ~ (1 | Person) + (1 | Task),
            data = brennan)

# Print basic results
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
# Components ordered by variance
tidy(g) %>%
  arrange(desc(var)) %>%
  mutate(pct_cumulative = cumsum(pct))

## -----------------------------------------------------------------------------
# G-study with profile likelihood confidence intervals
# Note: Rater is nested within Task, so we use Rater:Task notation
g_profile <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater:Task) +
    (1 | Person:Task),
  data = brennan,
  ci_method = "profile"
)

summary(g_profile)

## -----------------------------------------------------------------------------
# G-study with bootstrap confidence intervals (50 iterations for speed)
g_boot <- gstudy(
  Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  ci_method = "boot",
  nsim = 50,
  boot.type = "perc"
)

summary(g_boot)

## -----------------------------------------------------------------------------
# Method of moments G-study with nested design (Rater nested in Task)
g_mom <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater:Task) +
    (1 | Person:Task),
  data = brennan,
  backend = "mom"
)

summary(g_mom)

## -----------------------------------------------------------------------------
# Items nested within subtests (rajaratnam dataset)
g_nested <- gstudy(
  Score ~ (1 | Person) + (1 | Subtest) + (1 | Item:Subtest) +
    (1 | Person:Subtest),
  data = rajaratnam,
  backend = "mom"
)

summary(g_nested)

## ----eval=FALSE---------------------------------------------------------------
# # Bayesian G-study with default priors
# g_bayes <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater),
#   data = brennan,
#   backend = "brms"
# )
# 
# print(g_bayes)
# summary(g_bayes)

## ----eval=FALSE---------------------------------------------------------------
# # View default priors for a model
# library(brms)
# default_prior(Score ~ (1 | Person) + (1 | Task),
#               data = brennan, backend = "brms")
# 
# # Define custom priors
# my_prior <- c(
#   set_prior("exponential(1)", class = "sd", group = "Person"),
#   set_prior("exponential(0.5)", class = "sd", group = "Task"),
#   set_prior("student_t(3, 0, 2)", class = "sd", group = "Rater"),
#   set_prior("exponential(1)", class = "sigma")
# )
# 
# g_custom <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater),
#   data = brennan,
#   backend = "brms",
#   prior = my_prior
# )

## ----eval=FALSE---------------------------------------------------------------
# # Extract and check diagnostics
# library(posterior)
# 
# # View Rhat values (should be near 1.0)
# g_bayes$variance_components %>%
#   select(component, Rhat, Bulk_ESS, Tail_ESS)
# 
# # Extract posterior draws
# draws <- as_draws_matrix(g_bayes$model)
# 
# # Plot posterior distributions
# plot(g_bayes, type = "variance")

## -----------------------------------------------------------------------------
# Conduct D-study with 3 tasks and 4 raters
d <- dstudy(g, n = list(Task = 3, Rater = 4))

print(d)

## -----------------------------------------------------------------------------
# Access coefficients directly
d$coefficients

# Universe score variance (true variance)
cat("Universe score variance:", d$coefficients$uni, "\n")

# Relative error variance (affects rank-ordering)
cat("Relative error variance:", d$coefficients$sigma2_delta, "\n")

# Absolute error variance (affects absolute scores)
cat("Absolute error variance:", d$coefficients$sigma2_delta_abs, "\n")

# G coefficient
cat("G coefficient:", d$coefficients$g, "\n")

# Phi coefficient
cat("Phi coefficient:", d$coefficients$phi, "\n")

## -----------------------------------------------------------------------------
# Explore different numbers of tasks and raters
d_sweep <- dstudy(g, n = list(Task = c(3, 5, 10), Rater = c(2, 4, 8)))

# View all combinations
tidy(d_sweep)

## -----------------------------------------------------------------------------
# Plot sweep results
plot(d_sweep, type = "sweep")

## -----------------------------------------------------------------------------
# Extract grand mean from G-study
mu_y <- mean(brennan$Score)
cat("Grand mean:", mu_y, "\n")

# D-study with cut-score of 7
d_cut <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 7)

print(d_cut)

## -----------------------------------------------------------------------------
# Compare different cut-scores
d_cut_5 <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 5)
d_cut_7 <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 7)
d_cut_9 <- dstudy(g, n = list(Task = 3, Rater = 4), cut_score = 9)

cat("Phi-cut (c=5):", d_cut_5$coefficients$phi_cut, "\n")
cat("Phi-cut (c=7):", d_cut_7$coefficients$phi_cut, "\n")
cat("Phi-cut (c=9):", d_cut_9$coefficients$phi_cut, "\n")

## -----------------------------------------------------------------------------
# Include Person:Task interaction in universe
d_universe <- dstudy(
  g, 
  n = list(Task = 3, Rater = 4),
  universe = c("Person", "Person:Task")
)

summary(d_universe)

## -----------------------------------------------------------------------------
# Only use specific interactions as error sources
d_error <- dstudy(
  g,
  n = list(Task = 3, Rater = 4),
  error = c("Person:Task", "Person:Rater")
)

summary(d_error)

## -----------------------------------------------------------------------------
# Average scores across Raters
d_agg <- dstudy(
  g,
  n = list(Task = 3, Rater = 4),
  aggregation = "Rater",
  residual_is = "Person:Task:Rater"
)

summary(d_agg)

## ----eval=FALSE---------------------------------------------------------------
# # Bayesian G-study
# g_brms <- gstudy(
#   Score ~ (1 | Person) + (1 | Task),
#   data = brennan,
#   backend = "brms"
# )
# 
# # D-study with 95% credible intervals
# d_ci <- dstudy(
#   g_brms,
#   n = list(Task = 3),
#   ci = c("g", "phi")
# )
# 
# print(d_ci)

## ----eval=FALSE---------------------------------------------------------------
# # 90% credible interval
# d_ci_90 <- dstudy(
#   g_brms,
#   n = list(Task = 3),
#   ci = "g",
#   probs = c(0.05, 0.95)
# )

## ----eval=FALSE---------------------------------------------------------------
# # Prepare multivariate data (wide format)
# brennan_wide <- brennan %>%
#   tidyr::pivot_wider(
#     names_from = Task,
#     values_from = Score,
#     names_prefix = "Task"
#   )
# 
# head(brennan_wide)
# 
# # Multivariate G-study
# g_mv <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = brennan_wide,
#   backend = "brms"
# )
# 
# summary(g_mv)

## ----eval=FALSE---------------------------------------------------------------
# # View variance components by dimension
# tidy(g_mv)
# 
# # View correlations between dimensions
# g_mv$correlations$residual_cor
# 
# # Random effect correlations
# g_mv$correlations$random_effect_cor

## ----eval=FALSE---------------------------------------------------------------
# # D-study for multivariate model
# d_mv <- dstudy(g_mv, n = list(Rater = 4))
# 
# # Coefficients for each dimension
# tidy(d_mv)

## ----eval=FALSE---------------------------------------------------------------
# # Multivariate G-study
# g_mv <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = brennan_wide,
#   backend = "brms"
# )
# 
# # D-study with default equal weights
# d_mv <- dstudy(g_mv, n = list(Rater = 4))
# 
# # View coefficients (includes Composite row)
# print(d_mv)

## ----eval=FALSE---------------------------------------------------------------
# # Custom weights: Task1=50%, Task2=30%, Task3=20%
# d_weighted <- dstudy(
#   g_mv,
#   n = list(Rater = 4),
#   weights = c(0.5, 0.3, 0.2)
# )
# 
# # Compare composite coefficients
# cat(
#   "Equal weights G:",
#   d_mv$coefficients$g[d_mv$coefficients$dim == "Composite"],
#   "\n"
# )
# cat(
#   "Custom weights G:",
#   d_weighted$coefficients$g[d_weighted$coefficients$dim == "Composite"],
#   "\n"
# )

## ----eval=FALSE---------------------------------------------------------------
# library(facet)
# library(brms)
# 
# # Multivariate G-study with items nested within subtests
# # Note: This uses long-format specification
# g_mv <- gstudy(
#   bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person) + (0 + Subtest || Item),
#      sigma ~ 0 + Subtest),
#   data = rajaratnam,
#   backend = "brms",
#   chains = 2,
#   iter = 1000
# )
# 
# # D-study computing VAR
# d_mv <- dstudy(
#   g_mv,
#   n = list(Person = 5)
# )
# 
# # Access PRMSE and VAR results using prmse()
# prmse(d_mv)

## ----eval=FALSE---------------------------------------------------------------
# # Get a tidy data frame with all PRMSE and VAR results
# prmse(d_mv)
# 
# # Access posterior means directly (use [[ to avoid partial matching with variance_components)
# d_mv[["var"]]$var_rel     # VAR (relative) means
# d_mv[["var"]]$var_abs     # VAR (absolute) means
# d_mv[["var"]]$prmse_c_rel # PRMSE(C→S_i) (relative) means
# d_mv[["var"]]$prmse_c_abs # PRMSE(C→S_i) (absolute) means
# 
# # Access full posterior draws (for custom analysis)
# dim(d_mv[["var"]]$var_rel_draws)     # [n_draws, n_dims]
# dim(d_mv[["var"]]$prmse_c_rel_draws) # [n_draws, n_dims]

## ----eval=FALSE---------------------------------------------------------------
# # Default 95% credible intervals
# prmse(d_mv)
# 
# # Custom credible intervals (e.g., 90%)
# prmse(d_mv, probs = c(0.05, 0.95))

## ----eval=FALSE---------------------------------------------------------------
# # Check if VAR results are available (use [[ to avoid partial matching)
# if (!is.null(d_mv[["var"]])) {
#   prmse(d_mv)
# } else {
#   message("VAR not available for this model")
# }

## ----eval=FALSE---------------------------------------------------------------
# # Long-format specification
# library(brms)
# 
# g_long <- gstudy(
#   bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person) + (0 + Subtest || Item),
#      sigma ~ 0 + Subtest),
#   data = rajaratnam,
#   backend = "brms"
# )
# 
# # Check detection
# g_long$long_format_multivariate  # TRUE
# g_long$dimension_var             # "Subtest"
# g_long$dimensions                # c("A", "B", "C")
# 
# summary(g_long)

## -----------------------------------------------------------------------------
# Basic latent scores (Person only)
latent_scores <- latent(g)

head(latent_scores)

## -----------------------------------------------------------------------------
# Include Person:Task interaction
latent_interaction <- latent(g, universe = c("Person", "Person:Task"))

head(latent_interaction)

## ----eval=FALSE---------------------------------------------------------------
# # Latent scores from Bayesian G-study
# latent_bayes <- latent(g_bayes)
# 
# head(latent_bayes)

## -----------------------------------------------------------------------------
# Variance scale
plot(g, type = "variance")

# Proportion of variance
plot(g, type = "proportion")

## -----------------------------------------------------------------------------
# Forest plot (requires CIs)
g_ci <- gstudy(
  Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  ci_method = "profile"
)

plot(g_ci, type = "forest")

## -----------------------------------------------------------------------------
# All random effects
plot_ranef(g)

# Specific facet
plot_ranef(g, which = "Person")

## -----------------------------------------------------------------------------
# Create sweep
d_sweep <- dstudy(g, n = list(Task = c(2, 5, 10), Rater = c(2, 4, 8)))

# Plot sweep
plot(d_sweep, type = "sweep")

# Single coefficient
plot(d_sweep, type = "sweep", coefficient = "g")

## -----------------------------------------------------------------------------
# Access the underlying model object
model <- g$model

# For lme4 models
class(model)

# Extract random effects
re <- ranef(g)
re

# Extract variance-covariance
VarCorr(g)

## ----eval=FALSE---------------------------------------------------------------
# # Extract posterior draws
# draws <- brms::as_draws_matrix(g_bayes$model)
# 
# # Examine parameter names
# colnames(draws)
# 
# # Compute custom summaries
# posterior::summarise_draws(draws)

## -----------------------------------------------------------------------------
# Check for warnings
g$model@optinfo$conv$lme4

# Singular fits indicate near-zero variance components
lme4::isSingular(g$model)

## -----------------------------------------------------------------------------
# Check if any components were truncated
g_mom$variance_components %>%
  filter(var == 0)

## ----eval=FALSE---------------------------------------------------------------
# # Check Rhat values
# g_bayes$variance_components %>%
#   select(component, Rhat, Bulk_ESS, Tail_ESS) %>%
#   mutate(
#     converged = Rhat < 1.01,
#     sufficient_ess = Bulk_ESS > 400 & Tail_ESS > 400
#   )
# 
# # View estimation issues
# g_bayes$estimation_issues

## -----------------------------------------------------------------------------
library(facet)
library(dplyr)
library(tidyr)
library(broom)

# 1. Load and explore data
data(brennan)
summary(brennan)

# Check design structure (Raters nested within Task)
brennan %>%
  group_by(Task) %>%
  summarise(
    raters = paste(sort(unique(Rater)), collapse = ", "),
    n_raters = n_distinct(Rater)
  )

# 2. Conduct G-study with profile CIs (using nested design)
g_full <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater:Task) +
    (1 | Person:Task),
  data = brennan,
  ci_method = "profile"
)

# 3. Examine variance components
summary(g_full)

# Identify largest sources of variance
tidy(g_full) %>%
  arrange(desc(var)) %>%
  mutate(
    cum_pct = cumsum(pct),
    component = factor(component, levels = component)
  )

# 4. Conduct D-study to optimize design
# Current design: 3 tasks, 4 raters per task
d_current <- dstudy(g_full, n = list(Task = 3, `Rater:Task` = 4))
cat("Current design G:", d_current$coefficients$g, "\n")
cat("Current design Phi:", d_current$coefficients$phi, "\n")

# 5. Explore alternative designs
d_explore <- dstudy(
  g_full,
  n = list(
    Task = c(1, 2, 3, 5, 10),
    `Rater:Task` = c(1, 2, 4, 8, 12)
  )
)

# Find designs achieving G >= 0.80
tidy(d_explore) %>%
  filter(g >= 0.80) %>%
  arrange(g)

# 6. Visualize
plot(d_explore, type = "sweep")

# 7. Make recommendation
# Optimal design: minimum resources for G >= 0.80
optimal <- tidy(d_explore) %>%
  filter(g >= 0.80) %>%
  arrange(Task * `Rater:Task`) %>%
  slice(1)

print(optimal)

## -----------------------------------------------------------------------------
# G-study with nested design (Raters nested within Task)
g_3f <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater:Task) +
    (1 | Person:Task),
  data = brennan
)

summary(g_3f)

# D-study: what if we used 5 tasks with 2 raters each?
# Note: Use backticks for nested facet names in dstudy
d_3f <- dstudy(g_3f, n = list(Task = 5, `Rater:Task` = 2))

summary(d_3f)

# Compare with original design
d_orig <- dstudy(g_3f, n = list(Task = 3, `Rater:Task` = 4))

cat("Original (3 tasks, 4 raters/task): G =", d_orig$coefficients$g, "\n")
cat("Alternative (5 tasks, 2 raters/task): G =", d_3f$coefficients$g, "\n")

## ----eval=FALSE---------------------------------------------------------------
# # Optimize weights for maximum composite reliability
# result <- prmse(d_mv, optimize = "composite")
# 
# # View optimal weights
# result$weights
# 
# # View composite reliability achieved
# result$composite_reliability

## ----eval=FALSE---------------------------------------------------------------
# # Optimize weights for maximum subscale VAR
# result <- prmse(d_mv, optimize = "subscale")
# 
# # View minimax results (maximizes minimum VAR)
# result$minimax$weights
# result$minimax$min_var
# 
# # View per-subscale results
# result$per_subscale
# 
# # View comparison table
# result$comparison

## ----eval=FALSE---------------------------------------------------------------
# # Grid search with 0.2 resolution
# result <- prmse(d_mv, optimize = "tuning", grid_resolution = 0.2)
# 
# # View best solution
# result$best$weights
# 
# # View full grid results
# head(result$grid_results)
# 
# # Number of viable solutions (all VAR > 1)
# result$n_viable

## ----eval=FALSE---------------------------------------------------------------
# # Compare all three methods
# comp <- prmse(d_mv, optimize = "composite")
# sub <- prmse(d_mv, optimize = "subscale")
# tun <- prmse(d_mv, optimize = "tuning", grid_resolution = 0.1)
# 
# # Compare weights
# cat("Composite optimization weights:", round(comp$weights, 3), "\n")
# cat("Subscale minimax weights:", round(sub$minimax$weights, 3), "\n")
# cat("Tuning best weights:", round(tun$best$weights, 3), "\n")
# 
# # Compare metrics
# cat("Composite reliability:", round(comp$composite_reliability, 3), "\n")
# cat("Minimum VAR (minimax):", round(sub$minimax$min_var, 3), "\n")
# cat("Minimum VAR (tuning):", round(tun$best$min_var, 3), "\n")

## ----eval=FALSE---------------------------------------------------------------
# # Include composite reliability row
# prmse(d_mv, include_composite = TRUE)

