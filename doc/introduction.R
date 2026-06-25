## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.align = "center",
  cache = TRUE,
  cache.lazy = FALSE,
  message = FALSE,
  warning = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
# # Fully crossed design with REPLICATION per cell (e.g. each rater rates
# # every person on every task on multiple occasions). The three-way
# # Person:Task:Rater interaction IS estimable here because each cell has
# # more than one observation.
# g_crossed <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task) + (1 | Person:Rater) + (1 | Task:Rater) +
#     (1 | Person:Task:Rater),
#   data = crossed_data,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )

## ----eval=FALSE---------------------------------------------------------------
# g_nested <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) + (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )

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
# Canonical G-study for the brennan data: all main effects plus
# the only estimable interaction (Person:Task). The Person x Rater:Task
# interaction is confounded with the residual and absorbed into the
# error term.
g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
              (1 | Person:Task),
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
# Note: (1 | Rater) and (1 | Rater:Task) are equivalent in the brennan data
# because the 12 Rater labels are unique within each Task. We use the
# un-nested form so the formula reads as "all main effects + the only
# estimable interaction".
g_profile <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
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
# Fit the same model with both methods
g_profile <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
    (1 | Person:Task),
  data = brennan,
  ci_method = "profile"
)

# Compare the CIs
profile_ci <- tidy(g_profile) %>%
  select(component, var, lower, upper) %>%
  mutate(method = "profile")

# Bootstrap CIs were computed above in g_boot
boot_ci <- tidy(g_boot) %>%
  select(component, var, lower, upper) %>%
  mutate(method = "boot (nsim=50)")

bind_rows(profile_ci, boot_ci) %>%
  select(method, component, var, lower, upper) %>%
  mutate(across(c(var, lower, upper), ~ sprintf("%.3f", .x)))

## -----------------------------------------------------------------------------
# Method of moments G-study with nested design (Rater nested in Task;
# (1 | Rater) and (1 | Rater:Task) are equivalent here)
g_mom <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
    (1 | Person:Task),
  data = brennan,
  backend = "mom"
)

summary(g_mom)

## -----------------------------------------------------------------------------
# Items nested within subtests (rajaratnam dataset)
# (1 | ItemId) and (1 | Item:Subtest) are exactly equivalent here.
g_nested <- gstudy(
  Score ~ (1 | Person) + (1 | Subtest) + (1 | ItemId) +
    (1 | Person:Subtest),
  data = rajaratnam,
  backend = "mom"
)

summary(g_nested)

## ----eval=FALSE---------------------------------------------------------------
# library(facet)
# data(brennan)
# 
# # Examine structure
# head(brennan)
# #   Person Task Rater Score
# # 1      1    1     1     5
# # 2      2    1     1     9
# # 3      3    1     1     3
# # 4      4    1     1     7
# # 5      5    1     1     9
# # 6      6    1     1     3
# 
# str(brennan)

## ----eval=FALSE---------------------------------------------------------------
# # Bayesian G-study with default priors
# g_bayes <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )
# 
# print(g_bayes)
# summary(g_bayes)

## ----eval=FALSE---------------------------------------------------------------
# g_full <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   chains = 4,        # Number of chains
#   iter = 2000,       # Total iterations per chain
#   warmup = 1000,     # Warmup iterations
#   thin = 1,          # Thinning interval
#   cores = 4,         # Parallel processing
#   refresh = 1000,    # Progress message frequency
#   seed = 123         # Reproducibility
# )

## ----eval=FALSE---------------------------------------------------------------
# # View default priors for a model
# library(brms)
# default_prior(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#                 (1 | Person:Task),
#               data = brennan)
# 
# # Define custom priors
# my_prior <- c(
#   # Person variance: half-Cauchy with scale 1
#   set_prior("cauchy(0, 1)", class = "sd", group = "Person"),
# 
#   # Task variance: exponential with rate 1
#   set_prior("exponential(1)", class = "sd", group = "Task"),
# 
#   # Rater variance: exponential with rate 0.5 (more spread)
#   set_prior("exponential(0.5)", class = "sd", group = "Rater"),
# 
#   # Person:Task interaction: exponential with rate 1
#   set_prior("exponential(1)", class = "sd", group = "Person:Task"),
# 
#   # Residual: half-Cauchy with scale 1
#   set_prior("cauchy(0, 1)", class = "sigma")
# )
# 
# g_custom <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   prior = my_prior,
#   iter = 2000, cores = 4, refresh = 1000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Weakly informative prior
# prior_weak <- set_prior("cauchy(0, 5)", class = "sd")
# g_weak <- gstudy(
#   Score ~ (1 | Person) + (1 | Task),
#   data = brennan,
#   backend = "brms",
#   prior = prior_weak,
#   iter = 2000, cores = 4, refresh = 1000
# )
# 
# # More informative prior
# prior_info <- set_prior("exponential(1)", class = "sd")
# g_info <- gstudy(
#   Score ~ (1 | Person) + (1 | Task),
#   data = brennan,
#   backend = "brms",
#   prior = prior_info,
#   iter = 2000, cores = 4, refresh = 1000
# )
# 
# # Compare posterior means
# g_weak$variance_components$estimate
# g_info$variance_components$estimate

## ----eval=TRUE----------------------------------------------------------------
# Fit two models with contrasting priors
prior_weak_run <- set_prior("cauchy(0, 5)", class = "sd")
g_weak_run <- gstudy(
  Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  backend = "brms",
  prior = prior_weak_run,
  chains = 2, iter = 2000, warmup = 500, seed = 123, cores = 4, refresh = 1000
)

prior_info_run <- set_prior("exponential(1)", class = "sd")
g_info_run <- gstudy(
  Score ~ (1 | Person) + (1 | Task),
  data = brennan,
  backend = "brms",
  prior = prior_info_run,
  chains = 2, iter = 2000, warmup = 500, seed = 123, cores = 4, refresh = 1000
)

# Side-by-side comparison
bind_rows(
  tidy(g_weak_run)  %>% mutate(prior = "cauchy(0, 5)"),
  tidy(g_info_run)  %>% mutate(prior = "exponential(1)")
) %>%
  select(prior, component, var, lower, upper) %>%
  mutate(across(c(var, lower, upper), ~ sprintf("%.3f", .x)))

## ----eval=FALSE---------------------------------------------------------------
# # Extract and check diagnostics
# library(posterior)
# 
# # View Rhat values (should be near 1.0)
# g_bayes$variance_components %>%
#   select(component, Rhat, Bulk_ESS, Tail_ESS)
# 
# # Check for convergence problems
# any(g_bayes$variance_components$Rhat > 1.01)
# any(g_bayes$variance_components$Bulk_ESS < 400)
# 
# # Extract posterior draws
# draws <- as_draws_matrix(g_bayes$model)
# 
# # Plot posterior distributions
# plot(g_bayes, type = "variance")

## ----eval=FALSE---------------------------------------------------------------
# library(posterior)
# 
# # Extract draws matrix
# draws <- as_draws_matrix(g_bayes$model)
# 
# # View parameter names
# colnames(draws)
# 
# # Summary statistics
# summarise_draws(draws)

## ----eval=FALSE---------------------------------------------------------------
# # Variance components with posterior summaries
# g_bayes$variance_components

## ----eval=FALSE---------------------------------------------------------------
# # Different credible interval levels
# g_bayes$variance_components %>%
#   dplyr::mutate(
#     q50 = posterior::quantile2(draws[, "sd_Person__Intercept"], probs = 0.5),
#     q90_lower = posterior::quantile2(draws[, "sd_Person__Intercept"], probs = 0.05),
#     q90_upper = posterior::quantile2(draws[, "sd_Person__Intercept"], probs = 0.95)
#   )

## ----eval=FALSE---------------------------------------------------------------
# # Variance component distributions
# plot(g_bayes, type = "variance")
# 
# # Using bayesplot
# library(bayesplot)
# 
# # Posterior density plots
# mcmc_dens(draws, pars = c("sd_Person__Intercept", "sd_Task__Intercept"))
# 
# # Trace plots (should look like "fuzzy caterpillars")
# mcmc_trace(draws, pars = c("sd_Person__Intercept", "sd_Task__Intercept"))
# 
# # Pairs plot (check for correlations)
# mcmc_pairs(draws, pars = c("sd_Person__Intercept", "sd_Task__Intercept", "sigma"))

## ----eval=FALSE---------------------------------------------------------------
# # Posterior predictive check
# pp_check(g_bayes$model)
# 
# # Group-specific checks
# pp_check(g_bayes$model, type = "dens_overlay", group = "Person")

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
# Default universe: only Person is in the universe
d_default <- dstudy(g, n = list(Task = 3, Rater = 4))

# Custom universe: Person + Person:Task
d_custom <- dstudy(
  g,
  n = list(Task = 3, Rater = 4),
  universe = c("Person", "Person:Task")
)

# Side-by-side comparison
tribble(
  ~universe,                            ~G,    ~Phi,
  "default (Person only)",              round(d_default$coefficients$g, 3),
                                         round(d_default$coefficients$phi, 3),
  "custom (Person + Person:Task)",      round(d_custom$coefficients$g, 3),
                                         round(d_custom$coefficients$phi, 3)
)

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
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
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
# # Synthetic fully crossed design: 6 persons x 3 tasks x 4 raters
# set.seed(42)
# n_p <- 6; n_t <- 3; n_r <- 4
# crossed_data <- expand.grid(
#   Person = paste0("P", seq_len(n_p)),
#   Task   = seq_len(n_t),
#   Rater  = paste0("R", seq_len(n_r)),
#   stringsAsFactors = FALSE
# )
# mu <- 5
# task_eff  <- c(0.5, 0, -0.5)
# rater_eff <- c(0.2, 0.1, -0.1, -0.2)
# p_eff <- setNames(rnorm(n_p, 0, 1), paste0("P", seq_len(n_p)))
# crossed_data$Score <- with(crossed_data,
#   mu + p_eff[Person] +
#     task_eff[Task] +
#     rater_eff[as.integer(sub("R", "", Rater))] +
#     rnorm(nrow(crossed_data), 0, 0.8)   # residual + interactions
# )
# 
# # Pivot to wide format: one row per Person x Rater, one column per Task
# crossed_wide <- crossed_data %>%
#   tidyr::pivot_wider(
#     names_from  = Task,
#     values_from = Score,
#     names_prefix = "Task"
#   )
# 
# head(crossed_wide)
# #   Person Rater Task1 Task2 Task3
# # 1     P1    R1   8.3   5.5   4.1
# # 2     P2    R1   5.1   4.4   5.2
# # 3     P3    R1   7.7   5.5   4.8
# # 4     P4    R1   6.3   6.3   3.9
# # 5     P5    R1   7.2   5.4   5.0
# # 6     P6    R1   7.4   3.0   5.6
# 
# # Multivariate G-study
# g_mv <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = crossed_wide,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
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
# # Multivariate G-study on the fully crossed synthetic data
# g_mv <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = crossed_wide,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
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
# # Fully crossed synthetic data (defined in the Multivariate G-Study section)
# library(tidyr)
# crossed_wide <- crossed_data %>%
#   pivot_wider(
#     names_from = Task,
#     values_from = Score,
#     names_prefix = "Task"
#   )
# 
# head(crossed_wide)
# #   Person Rater Task1 Task2 Task3
# # 1     P1    R1   8.3   5.5   4.1
# # 2     P2    R1   5.1   4.4   5.2
# # 3     P3    R1   7.7   5.5   4.8
# # 4     P4    R1   6.3   6.3   3.9
# # 5     P5    R1   7.2   5.4   5.0
# # 6     P6    R1   7.4   3.0   5.6
# 
# # Multivariate G-study
# g_mv_wide <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = crossed_wide,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Use rajaratnam dataset (nested: items within subtests)
# data(rajaratnam)
# 
# head(rajaratnam)
# #   Person Subtest Item Score
# # 1      1       1    1     4
# # 2      1       1    2     5
# # ...
# 
# # Multivariate G-study with long format
# g_mv_long <- gstudy(
#   bf(Score ~ 0 + Subtest +
#          (0 + Subtest | r | Person) +
#          (0 + Subtest || ItemId)),
#   sigma ~ 0 + Subtest,
#   data = rajaratnam,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Variance components by dimension
# g_mv_wide$variance_components
# 
# # Residual correlations
# g_mv_wide$correlations$residual_cor
# 
# # Random effect correlations (e.g., Person effects)
# g_mv_wide$correlations$random_effect_cor$Person

## ----eval=FALSE---------------------------------------------------------------
# d_mv <- dstudy(g_mv_wide, n = list(Rater = 4), ci = c("g", "phi"))
# 
# print(d_mv)

## ----eval=FALSE---------------------------------------------------------------
# # Custom weights: Task1=50%, Task2=30%, Task3=20%
# d_weighted <- dstudy(
#   g_mv_wide,
#   n = list(Rater = 4),
#   weights = c(0.5, 0.3, 0.2),
#   ci = c("g", "phi")
# )
# 
# print(d_weighted)

## ----eval=FALSE---------------------------------------------------------------
# # Wide format: one column per dimension
# crossed_wide <- crossed_data %>%
#   tidyr::pivot_wider(
#     names_from = Task,
#     values_from = Score,
#     names_prefix = "Task"
#   )
# 
# # Long format on a multivariate variable: each row is a (person, task) score
# # and the multivariate response is built from the columns of crossed_data.
# # (The brms bf() syntax is more complex; see Multivariate G-Study section
# # for the full multivariate model.)
# 
# # Wide-format G-study (default for multivariate)
# g_wide_run <- gstudy(
#   mvbind(Task1, Task2, Task3) ~ (1 | Person) + (1 | Rater),
#   data = crossed_wide,
#   backend = "brms",
#   chains = 2, iter = 2000, warmup = 500, seed = 123, cores = 4, refresh = 1000
# )
# 
# # Print the wide-format variance components
# tidy(g_wide_run)

## ----eval=FALSE---------------------------------------------------------------
# library(facet)
# set.seed(123)
# n_persons <- 1500
# n_items <- 10
# 
# # True data-generating parameters:
# #   Var(tau_A) = Var(tau_B) = Var(tau_C) = 1.0
# #   Cov(tau_A, tau_B) = Cov(tau_A, tau_C) = Cov(tau_B, tau_C) = 0.5
# #   Item variance = 0.1 per subscale
# #   Person:Item variance = 0.9 per subscale
# Sigma_person <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), 3, 3,
#   dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
# 
# # Generate Person effects and build the data
# L <- chol(Sigma_person)
# tau <- matrix(rnorm(n_persons * 3), n_persons, 3) %*% L
# data <- do.call(rbind, lapply(seq_len(n_persons), function(p) {
#   data.frame(
#     person = factor(p),
#     item   = factor(seq_len(n_items)),
#     A = tau[p, 1] + rnorm(n_items, 0, sqrt(0.1)) + rnorm(n_items, 0, sqrt(0.9)),
#     B = tau[p, 2] + rnorm(n_items, 0, sqrt(0.1)) + rnorm(n_items, 0, sqrt(0.9)),
#     C = tau[p, 3] + rnorm(n_items, 0, sqrt(0.1)) + rnorm(n_items, 0, sqrt(0.9))
#   )
# }))
# 
# # Fit and run prmse
# g <- gstudy(mvbind(A, B, C) ~ (1 | person) + (1 | item), data = data, backend = "mom")
# d <- dstudy(g, n = list(item = 10), weights = c(A = 1, B = 1, C = 1))
# result <- suppressWarnings(prmse(d))
# 
# # Closed-form reference (Haberman 2008)
# n_i <- 10
# Sigma_tau <- Sigma_person
# Sigma_obs_rel <- Sigma_tau
# diag(Sigma_obs_rel) <- diag(Sigma_obs_rel) + 0.9 / n_i
# w <- rep(1/3, 3)
# var_tau_C <- as.numeric(t(w) %*% Sigma_tau %*% w)
# var_obs_C_rel <- as.numeric(t(w) %*% Sigma_obs_rel %*% w)
# Rel_C <- var_tau_C / var_obs_C_rel
# cov_tau_C <- as.numeric(Sigma_tau %*% w)
# 
# # PRMSE_S reference = G coefficient
# ref_s <- diag(Sigma_tau) / (diag(Sigma_tau) + 0.9 / n_i)
# 
# # PRMSE_C reference = Haberman closed form
# ref_c <- (cov_tau_C^2) / (diag(Sigma_tau) * Rel_C * var_obs_C_rel)
# 
# # Compare
# data.frame(
#   dim = c("A", "B", "C"),
#   ref_prmse_s = ref_s,
#   pkg_prmse_s = result$prmse_s_rel,
#   ref_prmse_c = ref_c,
#   pkg_prmse_c = result$prmse_c_rel
# )

## ----eval=FALSE---------------------------------------------------------------
# library(facet)
# library(brms)
# 
# # Multivariate G-study with items nested within subtests
# # Note: This uses long-format specification with the ItemId column
# g_mv <- gstudy(
#   bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person) + (0 + Subtest || ItemId),
#      sigma ~ 0 + Subtest),
#   data = rajaratnam,
#   backend = "brms",
#   chains = 2,
#   iter = 2000,
#   cores = 4,
#   refresh = 1000
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
# d_mv[["var"]]$prmse_c_rel # PRMSE(C->S_i) (relative) means
# d_mv[["var"]]$prmse_c_abs # PRMSE(C->S_i) (absolute) means
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
#   bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person) + (0 + Subtest || ItemId),
#      sigma ~ 0 + Subtest),
#   data = rajaratnam,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
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

# 2. Conduct G-study with profile CIs (canonical formula)
g_full <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
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
d_current <- dstudy(g_full, n = list(Task = 3, Rater = 4))
cat("Current design G:", d_current$coefficients$g, "\n")
cat("Current design Phi:", d_current$coefficients$phi, "\n")

# 5. Explore alternative designs
d_explore <- dstudy(
  g_full,
  n = list(
    Task = c(1, 2, 3, 5, 10),
    Rater = c(1, 2, 4, 8, 12)
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
  arrange(Task * Rater) %>%
  slice(1)

print(optimal)

## -----------------------------------------------------------------------------
# G-study with the canonical brennan formula
g_3f <- gstudy(
  Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
    (1 | Person:Task),
  data = brennan
)

summary(g_3f)

# D-study: what if we used 5 tasks with 2 raters each?
d_3f <- dstudy(g_3f, n = list(Task = 5, Rater = 2))

summary(d_3f)

# Compare with original design
d_orig <- dstudy(g_3f, n = list(Task = 3, Rater = 4))

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
# # Fit all three methods on the same multivariate D-study
# comp <- prmse(d_mv, optimize = "composite")
# sub  <- prmse(d_mv, optimize = "subscale")
# tun  <- prmse(d_mv, optimize = "tuning", grid_resolution = 0.1)
# 
# # Build a small comparison tibble
# tribble(
#   ~method,        ~weights,                            ~objective_value,    ~objective_name,
#   "composite",    round(comp$weights, 3),             round(comp$composite_reliability, 3),  "composite reliability",
#   "subscale (minimax)", round(sub$minimax$weights, 3), round(sub$minimax$min_var, 3),       "min VAR across subscales",
#   "tuning",       round(tun$best$weights, 3),         round(tun$best$min_var, 3),           "min VAR (grid search)"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Include composite reliability row
# prmse(d_mv, include_composite = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# # Compare models with LOO-CV
# library(loo)
# 
# # Fit alternative models
# g_simple <- gstudy(Score ~ (1 | Person) + (1 | Task),
#                    data = brennan, backend = "brms",
#                    iter = 2000, cores = 4, refresh = 1000)
# g_full <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#                    (1 | Person:Task),
#                  data = brennan, backend = "brms",
#                  iter = 2000, cores = 4, refresh = 1000)
# 
# # LOO comparison
# loo_simple <- loo(g_simple$model)
# loo_full <- loo(g_full$model)
# 
# loo_compare(loo_simple, loo_full)

## ----eval=FALSE---------------------------------------------------------------
# # Bayes factor approximation
# bayes_factor(g_full$model, g_simple$model)

## ----eval=FALSE---------------------------------------------------------------
# # Missing data is automatically handled by brms
# # using Bayesian imputation during model fitting
# 
# g_missing <- gstudy(
#   Score ~ (1 | Person) + (1 | Task),
#   data = data_with_missing,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )

## ----eval=FALSE---------------------------------------------------------------
# # Student-t likelihood for robustness
# g_robust <- gstudy(
#   bf(Score ~ (1 | Person) + (1 | Task),
#      family = student()),
#   data = brennan,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )

## ----eval=FALSE---------------------------------------------------------------
# tidy(g) %>%
#   mutate(
#     Variance = sprintf("%.3f", var),
#     `95% CI` = if("lower" %in% names(.)) {
#       sprintf("[%.3f, %.3f]", lower, upper)
#     } else {
#       "---"
#     },
#     `% of Total` = sprintf("%.1f%%", pct)
#   ) %>%
#   select(Component = component, Variance, `95% CI`, `% of Total`) %>%
#   knitr::kable()

## ----eval=FALSE---------------------------------------------------------------
# # Function to create APA-style table for variance components
# make_apa_variance_table <- function(gstudy_obj) {
#   vc <- gstudy_obj$variance_components
# 
#   # Create formatted table
#   data.frame(
#     Source = vc$component,
#     Estimate = sprintf("%.3f", vc$estimate),
#     SD = sprintf("%.3f", vc$sd),
#     `95% CI` = sprintf("[%.3f, %.3f]", vc$lower, vc$upper),
#     Rhat = sprintf("%.3f", vc$Rhat),
#     ESS = round(vc$Bulk_ESS)
#   )
# }

## ----eval=FALSE---------------------------------------------------------------
# library(facet)
# data(brennan)
# 
# # a) Bayesian G-study
# g_bayes <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   iter = 2000, cores = 4, refresh = 1000
# )
# 
# # b) Posterior summaries
# print(g_bayes)
# 
# # c) Compare to REML
# g_reml <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "lme4"
# )
# 
# # Compare estimates
# g_bayes$variance_components$estimate
# g_reml$variance_components$var

## ----eval=FALSE---------------------------------------------------------------
# library(brms)
# 
# # a) Custom prior encoding prior information
# # "Moderate" variance with SD around 1.0 suggests half-Cauchy(0, 1)
# # or exponential(1) for the SD
# custom_prior <- c(
#   set_prior("cauchy(0, 1)", class = "sd", group = "Person"),
#   set_prior("cauchy(0, 2)", class = "sd", group = "Task"),
#   set_prior("cauchy(0, 2)", class = "sd", group = "Rater"),
#   set_prior("cauchy(0, 2)", class = "sd", group = "Person:Task"),
#   set_prior("cauchy(0, 2)", class = "sigma")
# )
# 
# # b) Fit with custom priors
# g_custom <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   prior = custom_prior,
#   iter = 2000, cores = 4, refresh = 1000
# )
# 
# # Compare
# g_bayes$variance_components$estimate
# g_custom$variance_components$estimate
# 
# # c) Sensitivity analysis
# scales <- c(0.5, 1, 2, 5)
# results <- lapply(scales, function(s) {
#   prior_s <- set_prior(paste0("cauchy(0, ", s, ")"), class = "sd")
#   g_s <- gstudy(
#     Score ~ (1 | Person) + (1 | Task),
#     data = brennan,
#     backend = "brms",
#     prior = prior_s
#   )
#   g_s$variance_components$estimate[1]  # Person variance
# })

## ----eval=FALSE---------------------------------------------------------------
# # a) G-study
# g_ex3 <- gstudy(
#   Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#     (1 | Person:Task),
#   data = brennan,
#   backend = "brms",
#   chains = 4,
#   iter = 2000,
#   cores = 4,
#   refresh = 1000
# )
# 
# # b) D-study sweep
# d_ex3 <- dstudy(
#   g_ex3,
#   n = list(Task = 2:5, Rater = 1:3),
#   ci = c("g", "phi")
# )
# 
# # c) Designs achieving G > 0.80
# library(dplyr)
# d_ex3$coefficients %>%
#   filter(g_lower > 0.80) %>%
#   arrange(g_lower)
# 
# # d) Visualization
# plot(d_ex3, type = "sweep")

## ----eval=FALSE---------------------------------------------------------------
# data(rajaratnam)
# 
# # a) Multivariate G-study
# g_mv_ex4 <- gstudy(
#   bf(Score ~ 0 + Subtest +
#          (0 + Subtest | r | Person) +
#          (0 + Subtest || ItemId)),
#   sigma ~ 0 + Subtest,
#   data = rajaratnam,
#   backend = "brms",
#   chains = 2,
#   iter = 2000,
#   cores = 4,
#   refresh = 1000
# )
# 
# # b) VAR with credible intervals
# d_mv_ex4 <- dstudy(
#   g_mv_ex4,
#   n = list(Person = 5),
#   ci = c("var_rel", "var_abs")
# )
# 
# prmse(d_mv_ex4)
# 
# # c) Interpretation
# # VAR > 1 indicates subscale adds value
# # VAR < 1 indicates composite is better
# # VAR $\approx$ 1 means either approach acceptable

## ----eval=FALSE---------------------------------------------------------------
# # a) Three prior specifications
# prior_weak <- set_prior("cauchy(0, 5)", class = "sd")
# prior_mod <- set_prior("cauchy(0, 1)", class = "sd")
# prior_strong <- set_prior("exponential(2)", class = "sd")
# 
# g_weak <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#                  (1 | Person:Task),
#                  data = brennan, backend = "brms", prior = prior_weak,
#                  iter = 2000, cores = 4, refresh = 1000)
# g_mod <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#                 (1 | Person:Task),
#                 data = brennan, backend = "brms", prior = prior_mod,
#                 iter = 2000, cores = 4, refresh = 1000)
# g_strong <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater) +
#                    (1 | Person:Task),
#                    data = brennan, backend = "brms", prior = prior_strong,
#                    iter = 2000, cores = 4, refresh = 1000)
# 
# # b) Compare posteriors
# library(ggplot2)
# draws_weak <- as_draws_df(g_weak$model)
# draws_mod <- as_draws_df(g_mod$model)
# draws_strong <- as_draws_df(g_strong$model)
# 
# # c) Prior sensitivity matters most when:
# # - Sample size is small
# # - Variance components are near zero
# # - Priors strongly conflict with data

