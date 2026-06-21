#' Tidy PRMSE and VAR Results
#'
#' Extract PRMSE(C→S_i) and VAR results as a data frame from a dstudy object.
#' For univariate models, returns the G and Phi coefficients as PRMSE metrics.
#'
#' @section Experimental:
#' This function is experimental. Testing has been limited, and results
#' should be verified before being cited in published work.
#'
#' @param dstudy_obj A dstudy object from \code{\link{dstudy()}}
#' @param include_composite Logical. If TRUE, adds a Composite row with summary
#' metrics. Default FALSE. The Composite row shows composite reliability (G/Phi)
#' in prmse_s_rel/abs columns, with prmse_c_rel/abs and var_rel/abs set to 1.0
#' (trivially, the composite perfectly predicts itself).
#' @param include_profile Logical. If TRUE (default), includes prmse_p_rel and
#' prmse_p_abs columns showing how well the entire profile predicts each subscale.
#' If FALSE, these columns are omitted for a more compact output.
#' @param ci Character vector specifying which metrics to compute CIs for.
#' Options: "prmse", "var", or both. Default NULL (no CIs).
#' **Note:** Credible intervals are only available for brms backend with posterior
#' estimation. For mom backend, this parameter is ignored with a warning.
#' @param probs Numeric vector of length 2 specifying the quantile probabilities
#' for credible intervals (brms backend only). Default is \code{c(0.025, 0.975)} for 95% CI.
#' @param weights Numeric vector of weights for computing composite coefficients.
#' Length must match the number of dimensions. Default NULL uses weights from dstudy.
#' When alternative weights are provided, all calculations are performed using
#' posterior draws to properly propagate uncertainty.
#' @param optimize Character string specifying optimization method. One of:
#' \describe{
#' \item{"composite"}{Maximize composite reliability via generalized eigenvalue}
#' \item{"subscale"}{Maximize VAR for subscales (per-subscale and minimax)}
#' \item{"tuning"}{Grid search guided by VAR diagnostics}
#' }
#' Default NULL (no optimization).
#' @param optimize_target Character string specifying which VAR to optimize.
#' Either "rel" (G-based) or "abs" (Phi-based). Default "rel".
#' @param grid_resolution Numeric step size for grid search in "tuning" method.
#' Must be between 0.01 and 0.5. Default 0.1.
#' @param subscale Character string specifying target subscale for "subscale"
#' method. If NULL (default), returns all per-subscale results plus minimax.
#'
#' @return For univariate models, a tibble with:
#' \describe{
#' \item{dim}{The response variable name}
#' \item{prmse_rel}{The G coefficient (relative reliability)}
#' \item{prmse_abs}{The Phi coefficient (absolute reliability)}
#' \item{prmse_rel_LL, prmse_rel_UL, prmse_abs_LL, prmse_abs_UL}{Confidence intervals (if ci specified)}
#' }
#'
#' For multivariate models, a tibble with one row per dimension and:
#' \describe{
#' \item{dim}{Dimension/subscale name}
#' \item{prmse_s_rel, prmse_s_abs}{
#' Univariate subscale reliabilities (G and Phi coefficients).
#' These are baseline measures that ignore information from other dimensions.
#' }
#' \item{prmse_c_rel, prmse_c_abs}{
#' PRMSE when predicting subscale true scores from the weighted composite.
#' Used for subscale viability analysis: if PRMSE_S > PRMSE_C,
#' the subscale provides unique information beyond the composite.
#' }
#' \item{prmse_p_rel, prmse_p_abs}{
#' PRMSE when projecting from the entire profile (all dimensions).
#' Can exceed PRMSE_S when dimensions are correlated by borrowing
#' predictive information from related dimensions.
#' Only included when \code{include_profile = TRUE} (default).
#' }
#' \item{var_rel, var_abs}{Value Added Ratio (PRMSE_S / PRMSE_C).}
#' \item{..._LL, ..._UL}{Confidence interval limits for requested metrics.}
#' }
#'
#' When \code{include_composite = TRUE}, a Composite row is added with:
#' \itemize{
#' \item prmse_s_rel/abs: Composite reliability (G/Phi coefficients)
#' \item prmse_c_rel/abs: 1.0 (composite trivially predicts itself)
#' \item prmse_p_rel/abs: Multivariate reliability (if include_profile = TRUE)
#' \item var_rel/abs: 1.0
#' }
#'
#' When \code{optimize} is specified, returns a list with optimization results
#' containing the optimal weights, eigenvalue/value, and recomputed metrics.
#'
#' @details
#' ## Univariate Models
#'
#' For univariate models, PRMSE_rel and PRMSE_abs correspond directly to the
#' G and Phi coefficients, respectively:
#' \itemize{
#' \item \code{prmse_rel} = G coefficient (relative decision reliability)
#' \item \code{prmse_abs} = Phi coefficient (absolute decision reliability)
#' }
#'
#' These are the same metrics but labeled differently to maintain consistency
#' with the multivariate output structure.
#'
#' ## Multivariate Models: Three Types of PRMSE
#'
#' For multivariate models, three types of PRMSE metrics are computed,
#' each answering a different question about prediction accuracy:
#'
#' ### PRMSE_S (Subscale Reliability)
#'
#' The univariate G coefficient for each dimension. This measures how well
#' the subscale's own observed scores predict its true scores, ignoring
#' information from other dimensions.
#'
#' \deqn{\text{PRMSE}_S = G = \frac{\sigma^2_{\text{universe}}}{\sigma^2_{\text{universe}} + \sigma^2_{\text{relative error}}}}
#'
#' This is the baseline reliability for each dimension when considered in isolation.
#'
#' ### PRMSE_C (Composite Prediction)
#'
#' How well the weighted composite score predicts each subscale's true scores.
#' This is the Haberman (2008) formula:
#'
#' \deqn{\mathrm{PRMSE}(C \to S_i) = \frac{[\mathrm{Cov}(\tau_i, C)]^2}
#'       {\mathrm{Var}(\tau_i) \cdot \mathrm{Rel}(C) \cdot \mathrm{Var}(C)}}
#'
#' where
#' \itemize{
#'   \item \eqn{\mathrm{Cov}(\tau_i, C) = \sum_j w_j \, \mathrm{Cov}(\tau_i, \tau_j)}
#'     is the weighted covariance between subscale \eqn{i}'s true score and the
#'     composite true score, computed from the universe-score (disattenuated)
#'     covariance matrix, not from observed-score sample covariances.
#'   \item \eqn{\mathrm{Var}(\tau_i) = \Sigma_{\tau,ii}} is the universe-score
#'     variance of subscale \eqn{i}.
#'   \item \eqn{\mathrm{Rel}(C) = \mathrm{Var}(\tau_C) / \mathrm{Var}(C)} is
#'     the composite reliability (G coefficient for relative, \eqn{\Phi} for
#'     absolute).
#'   \item \eqn{\mathrm{Var}(C) = \mathbf{w}^\top \Sigma_{\mathrm{obs}} \mathbf{w}}
#'     is the D-study observed variance of the composite (scaled for the
#'     decision study sample sizes).
#' }
#'
#' **Note on Vispoel et al. (2023) Equation 38:** Vispoel et al. (2023, p. 9)
#' label this quantity as "Equation 38", but the printed equation contains a
#' typesetting error (the denominator uses \eqn{\hat\sigma^2(S_i) \cdot
#' \mathrm{Rel}(C)} where it should be \eqn{\mathrm{Var}(\tau_i) \cdot
#' \mathrm{Rel}(C)}). This implementation uses the corrected Haberman (2008)
#' form, which matches the numerical results reported in Vispoel et al.'s
#' Table 8.
#'
#' **Viability Analysis:** Subscale viability is supported when
#' \eqn{\text{PRMSE}_S > \text{PRMSE}_C}, meaning the subscale predicts its own
#' true scores better than the composite does. The Value Added Ratio (VAR)
#' formalizes this: \eqn{\text{VAR} = \text{PRMSE}_S / \text{PRMSE}_C}.
#'
#' ### PRMSE_P (Profile Projection)
#'
#' How well the entire profile (all dimensions simultaneously) predicts each
#' subscale's true scores. This uses the multivariate structure to potentially
#' improve prediction accuracy by borrowing information from correlated dimensions:
#'
#' \deqn{\text{PRMSE}_P[d] = \frac{[\Sigma_\tau \Sigma_{\text{obs}}^{-1} \Sigma_\tau]_{dd}}{\Sigma_{\tau,dd}}}
#'
#' where \eqn{\Sigma_\tau} is the universe-score covariance matrix and
#' \eqn{\Sigma_{\text{obs}}} is the D-study observed covariance matrix.
#'
#' **Diagonal slice vs. trace:** \code{prmse_p_rel} returns the \emph{diagonal
#' slice} of the matrix \eqn{\Sigma_\tau \Sigma_{\text{obs}}^{-1} \Sigma_\tau}
#' (one value per dimension). The companion attribute
#' \code{attr(result, "prmse_mv_rel")} returns the \emph{trace} of the same
#' matrix, normalized by the trace of \eqn{\Sigma_\tau}, which is a single
#' global multivariate profile reliability:
#'
#' \deqn{\text{PRMSE}_{MV} = \frac{\mathrm{tr}(\Sigma_\tau \Sigma_{\text{obs}}^{-1} \Sigma_\tau)}{\mathrm{tr}(\Sigma_\tau)}}
#'
#' This trace equals the matrix analog of a squared canonical correlation and
#' always satisfies \eqn{\text{PRMSE}_{MV} \geq \text{PRMSE}_P[d]} for any
#' dimension \eqn{d} (by Cauchy-Schwarz on the Cholesky factor of
#' \eqn{\Sigma_{\text{obs}}^{-1/2}\Sigma_\tau}).
#'
#' **Relationship to PRMSE_S:** When dimensions are uncorrelated (all off-diagonal
#' covariances = 0), the matrices become diagonal and PRMSE_P simplifies to PRMSE_S.
#' In this case, \code{prmse_p_rel} will approximately equal \code{prmse_s_rel}.
#'
#' These metrics will only differ meaningfully when:
#' \itemize{
#' \item Dimensions are correlated with each other
#' \item The correlations are substantial enough to provide predictive value
#' }
#'
#' A large difference between PRMSE_P and PRMSE_S indicates that the multivariate
#' structure provides additional predictive information beyond what each dimension
#' offers alone. This can occur when one dimension has high reliability and is
#' correlated with a less reliable dimension, allowing the profile to "borrow
#' strength" for improved prediction.
#'
#' ### Summary Comparison
#'
#' | Metric | Question Answered |
#' |--------|-------------------|
#' | PRMSE_S | "How reliable is this subscale alone?" |
#' | PRMSE_C | "How well does the composite predict this subscale?" |
#' | PRMSE_P | "How well does the full profile predict this subscale?" |
#'
#' ## Composite Row
#'
#' When \code{include_composite = TRUE}, a Composite row is added to the output.
#' This row summarizes the composite score's reliability:
#' \itemize{
#' \item \code{prmse_s_rel/abs}: The composite's G and Phi coefficients
#' \item \code{prmse_c_rel/abs}: Always 1.0 (the composite perfectly predicts itself)
#' \item \code{prmse_p_rel/abs}: Multivariate reliability (when include_profile = TRUE)
#' \item \code{var_rel/abs}: Always 1.0
#' }
#'
#' ## Weight Optimization
#'
#' The three optimization methods address a fundamental tension in multivariate
#' generalizability theory:
#'
#' - **"composite"**: Finds weights that maximize composite reliability (PRMSE).
#' Uses analytical solution via generalized eigenvalue decomposition.
#' Best for selection/norm-referencing decisions.
#' - **"subscale"**: Finds weights that maximize VAR for subscales. Returns both
#' per-subscale optimizations and minimax optimization.
#' Best for profile interpretation where all subscales matter.
#' - **"tuning"**: Exhaustive grid search over weight combinations.
#' Useful for diagnostic purposes and understanding weight-VAR relationships.
#'
#' ## Confidence Intervals
#'
#' Credible intervals are only available when using the brms backend with posterior
#' estimation. These are computed directly from the posterior draws of the variance
#' components and PRMSE metrics.
#'
#' For the mom backend (method of moments), credible intervals are not available
#' because the method produces only point estimates without a posterior distribution.
#' If credible intervals are needed, consider using the brms backend:
#'
#' \preformatted{
#' g_mv <- gstudy(formula, data, backend = "brms")
#' d_mv <- dstudy(g_mv, n = ...)
#' prmse(d_mv, ci = "prmse") # CIs from posterior draws
#' }
#'
#' @seealso
#' \code{\link{dstudy}()} for computing D-study results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Univariate model
#' g_uni <- gstudy(Score ~ (1 | Person) + (1 | Item), data = mydata)
#' d_uni <- dstudy(g_uni, n = list(Item = 5))
#' prmse(d_uni) # Returns G and Phi as prmse_rel and prmse_abs
#'
#' # Multivariate model with mom backend (no CIs available)
#' g_mv <- gstudy(Score ~ 0 + Subtest + (0 + Subtest | Person),
#'   data = data, dimension_var = "Subtest", backend = "mom")
#' d_mv <- dstudy(g_mv, n = list(Person = 5))
#' prmse(d_mv)
#'
#' # Include Composite row
#' prmse(d_mv, include_composite = TRUE)
#'
#' # Exclude profile projection metrics for compact output
#' prmse(d_mv, include_profile = FALSE)
#'
#' # Understanding PRMSE_S vs PRMSE_P:
#' # When dimensions are uncorrelated, prmse_p_rel ≈ prmse_s_rel
#' # When correlated, prmse_p_rel may exceed prmse_s_rel
#' result <- prmse(d_mv)
#' diff <- result$prmse_s_rel - result$prmse_p_rel
#' if (all(abs(diff) < 0.01)) {
#'   message("Dimensions appear uncorrelated: PRMSE_P ≈ PRMSE_S")
#' } else {
#'   message("Dimensions are correlated: PRMSE_P differs from PRMSE_S")
#' }
#'
#' # Weight optimization
#' prmse(d_mv, optimize = "composite")
#' prmse(d_mv, optimize = "subscale")
#'
#' # Multivariate model with brms backend (credible intervals available)
#' library(brms)
#' g_mv <- gstudy(
#'   bf(Score ~ 0 + Subtest + (0 + Subtest | r | Person)),
#'   data = data, backend = "brms"
#' )
#' d_mv <- dstudy(g_mv, n = list(Person = 5))
#' prmse(d_mv, ci = "prmse") # Returns point estimates with credible intervals
#' }
prmse <- function(dstudy_obj,
  include_composite = FALSE,
  include_profile = TRUE,
  ci = NULL,
  probs = c(0.025, 0.975),
  weights = NULL,
  optimize = NULL,
  optimize_target = "rel",
  grid_resolution = 0.1,
  subscale = NULL) {

  if (!inherits(dstudy_obj, "dstudy")) {
    stop("'dstudy_obj' must be a dstudy object", call. = FALSE)
  }

  warning(
    "This is an experimental function. Testing has been limited.\n",
    "Results should be verified, ideally before you cite them.",
    call. = FALSE
  )

  if (!is.null(ci)) {
    ci <- match.arg(ci, c("prmse", "var"), several.ok = TRUE)
  }

  if (length(probs) != 2 || probs[1] >= probs[2]) {
    stop("'probs' must be a numeric vector of length 2 with probs[1] < probs[2]", call. = FALSE)
  }

  if (!is.null(optimize)) {
    optimize <- match.arg(optimize, c("composite", "subscale", "tuning"))
    optimize_target <- match.arg(optimize_target, c("rel", "abs"))

    if (grid_resolution < 0.01) {
      warning("'grid_resolution' (", grid_resolution, ") is below the minimum 0.01; using 0.01.",
              call. = FALSE)
      grid_resolution <- 0.01
    } else if (grid_resolution > 0.5) {
      warning("'grid_resolution' (", grid_resolution, ") is above the maximum 0.5; using 0.5.",
              call. = FALSE)
      grid_resolution <- 0.5
    }

    if (!isTRUE(dstudy_obj$is_multivariate)) {
      stop("Weight optimization is only available for multivariate models.", call. = FALSE)
    }

    result <- switch(optimize,
      "composite" = optimize_composite_weights_internal(dstudy_obj, optimize_target),
      "subscale" = optimize_subscale_weights_internal(dstudy_obj, subscale, optimize_target),
      "tuning" = optimize_tuning_weights_internal(dstudy_obj, grid_resolution, optimize_target)
    )

    dims <- dstudy_obj$gstudy$dimensions
    if (!is.null(result$metrics)) {
      result$metrics <- recalculate_var_with_weights(dstudy_obj, result$weights, dims, ci, probs, n = NULL)
    }
    if (optimize == "subscale" && is.null(subscale)) {
      if (!is.null(result$minimax$metrics)) {
        result$minimax$metrics <- recalculate_var_with_weights(dstudy_obj, result$minimax$weights, dims, ci, probs, n = NULL)
      }
    }
    if (optimize == "tuning") {
      if (!is.null(result$best$metrics)) {
        result$best$metrics <- recalculate_var_with_weights(dstudy_obj, result$best$weights, dims, ci, probs, n = NULL)
      }
    }

    return(result)
  }

  if (!isTRUE(dstudy_obj$is_multivariate)) {
    message(
      "For univariate models, PRMSE_rel and PRMSE_abs correspond to the G and Phi coefficients:\n",
      " - prmse_rel = G coefficient (relative decision reliability)\n",
      " - prmse_abs = Phi coefficient (absolute decision reliability)\n",
      "These are returned in a tibble for consistency with multivariate output."
    )

    if (!is.null(dstudy_obj$coefficients)) {
      coef_df <- dstudy_obj$coefficients
    } else if (!is.null(dstudy_obj$summary)) {
      coef_df <- dstudy_obj$summary
    } else {
      stop("No coefficients or summary found in dstudy object", call. = FALSE)
    }

    res <- tibble::tibble(
      dim = coef_df$dim,
      prmse_rel = coef_df$g,
      prmse_abs = coef_df$phi
    )

    if (!is.null(ci) && "prmse" %in% ci) {
      if (!is.null(coef_df$g_LL)) {
        res$prmse_rel_LL <- coef_df$g_LL
        res$prmse_rel_UL <- coef_df$g_UL
        res$prmse_abs_LL <- coef_df$phi_LL
        res$prmse_abs_UL <- coef_df$phi_UL
      } else {
        warning(
          "Credible intervals are only available for brms backend with posterior estimation. ",
          "For mom backend, consider using a Bayesian approach (backend = 'brms') if CIs are needed.",
          call. = FALSE
        )
      }
    }

    return(res)
  }

  if (is.null(dstudy_obj$var)) {
    if (!isTRUE(dstudy_obj$is_sweep)) {
      stop("No VAR results available. VAR is only computed for multivariate models.", call. = FALSE)
    }
    coef_df <- dstudy_obj$coefficients
    dims <- unique(coef_df$dim)
    dims <- dims[dims != "Composite"]
    w <- dstudy_obj$weights
    if (is.null(w) || length(w) != length(dims)) {
      w <- rep(1, length(dims))
    }
    names(w) <- dims
    return(recalculate_var_with_weights(dstudy_obj, w, dims, ci, probs, n = NULL))
  }

  coef_df <- dstudy_obj$coefficients
  dims <- unique(coef_df$dim)
  dims <- dims[dims != "Composite"]
  n_dims <- length(dims)

  if (length(dims) == 1) {
    message(
      "Only one dimension found. Treating as univariate:\n",
      " - PRMSE_P is not defined for single dimension\n",
      " - Use prmse_s and prmse_c metrics"
    )
  }

  if (!is.null(weights)) {
    if (length(weights) != n_dims) {
      stop("weights must have length ", n_dims, " (number of dimensions), got ", length(weights), call. = FALSE)
    }
    names(weights) <- dims
    return(recalculate_var_with_weights(dstudy_obj, weights, dims, ci, probs, n = NULL))
  }

  var <- dstudy_obj$var
  has_draws <- !is.null(var$prmse_s_rel_draws) || !is.null(var$var_rel_draws)

  if (has_draws) {
    res <- tibble::tibble(dim = dims)

    add_metric <- function(df, draws, name) {
      if (is.null(draws)) {
        df[[name]] <- NA_real_
        df[[paste0(name, "_LL")]] <- NA_real_
        df[[paste0(name, "_UL")]] <- NA_real_
        return(df)
      }
      if (is.matrix(draws)) {
        df[[name]] <- colMeans(draws, na.rm = TRUE)
        df[[paste0(name, "_LL")]] <- apply(draws, 2, quantile, probs = probs[1], na.rm = TRUE)
        df[[paste0(name, "_UL")]] <- apply(draws, 2, quantile, probs = probs[2], na.rm = TRUE)
      } else {
        df[[name]] <- rep(mean(draws, na.rm = TRUE), nrow(df))
        df[[paste0(name, "_LL")]] <- rep(quantile(draws, probs = probs[1], na.rm = TRUE), nrow(df))
        df[[paste0(name, "_UL")]] <- rep(quantile(draws, probs = probs[2], na.rm = TRUE), nrow(df))
      }
      df
    }

    res <- add_metric(res, var$prmse_s_rel_draws, "prmse_s_rel")
    res <- add_metric(res, var$prmse_s_abs_draws, "prmse_s_abs")
    res <- add_metric(res, var$prmse_c_rel_draws, "prmse_c_rel")
    res <- add_metric(res, var$prmse_c_abs_draws, "prmse_c_abs")
    if (include_profile) {
      res <- add_metric(res, var$prmse_p_rel_draws, "prmse_p_rel")
      res <- add_metric(res, var$prmse_p_abs_draws, "prmse_p_abs")
    }
    res <- add_metric(res, var$var_rel_draws, "var_rel")
    res <- add_metric(res, var$var_abs_draws, "var_abs")

    if (is.null(ci)) {
      ci_cols <- grep("_LL$|_UL$", names(res), value = TRUE)
      res <- res[, !(names(res) %in% ci_cols), drop = FALSE]
    } else if (!("prmse" %in% ci)) {
      prmse_ci_cols <- grep("^prmse_s_|^prmse_c_|^prmse_p_", grep("_LL$|_UL$", names(res), value = TRUE), value = TRUE)
      res <- res[, !(names(res) %in% prmse_ci_cols), drop = FALSE]
    } else if (!("var" %in% ci)) {
      var_ci_cols <- grep("^var_rel_|^var_abs_", grep("_LL$|_UL$", names(res), value = TRUE), value = TRUE)
      res <- res[, !(names(res) %in% var_ci_cols), drop = FALSE]
    }

    prmse_mv_rel_val <- if (!is.null(var$prmse_mv_rel_draws)) {
      c(mean = mean(var$prmse_mv_rel_draws, na.rm = TRUE),
        LL = quantile(var$prmse_mv_rel_draws, probs = probs[1], na.rm = TRUE),
        UL = quantile(var$prmse_mv_rel_draws, probs = probs[2], na.rm = TRUE))
    } else {
      c(mean = NA_real_, LL = NA_real_, UL = NA_real_)
    }

    prmse_mv_abs_val <- if (!is.null(var$prmse_mv_abs_draws)) {
      c(mean = mean(var$prmse_mv_abs_draws, na.rm = TRUE),
        LL = quantile(var$prmse_mv_abs_draws, probs = probs[1], na.rm = TRUE),
        UL = quantile(var$prmse_mv_abs_draws, probs = probs[2], na.rm = TRUE))
    } else {
      c(mean = NA_real_, LL = NA_real_, UL = NA_real_)
    }

    attr(res, "prmse_mv_rel") <- prmse_mv_rel_val
    attr(res, "prmse_mv_abs") <- prmse_mv_abs_val

    if (include_composite) {
      comp_row <- coef_df[coef_df$dim == "Composite", ]
      comp_g <- if (nrow(comp_row) > 0) comp_row$g[1] else NA_real_
      comp_phi <- if (nrow(comp_row) > 0) comp_row$phi[1] else NA_real_

      comp_cols <- list(dim = "Composite")
      comp_cols$prmse_s_rel <- unname(comp_g)
      comp_cols$prmse_s_abs <- unname(comp_phi)
      comp_cols$prmse_c_rel <- unname(1.0)
      comp_cols$prmse_c_abs <- unname(1.0)
      if (include_profile) {
        comp_cols$prmse_p_rel <- unname(prmse_mv_rel_val["mean"])
        comp_cols$prmse_p_abs <- unname(prmse_mv_abs_val["mean"])
      }
      comp_cols$var_rel <- unname(1.0)
      comp_cols$var_abs <- unname(1.0)

      ci_cols <- grep("_LL$|_UL$", names(res), value = TRUE)
      needs_draws_ci <- any(grepl("^prmse_s_", ci_cols)) ||
        (include_profile && any(grepl("^prmse_p_", ci_cols)))

      if (needs_draws_ci) {
        w <- dstudy_obj$weights
        if (is.null(w) || length(w) != n_dims) {
          w <- rep(1, n_dims)
        }
        names(w) <- dims
        draws_ci <- recalculate_var_with_weights(
          dstudy_obj, w, dims, ci, probs, n = NULL
        )
        comp_draws <- draws_ci[draws_ci$dim == "Composite", , drop = FALSE]
        if (nrow(comp_draws) == 1) {
          for (col in ci_cols) {
            if (grepl("^var_", col) || grepl("_c_", col)) {
              comp_cols[[col]] <- unname(1.0)
            } else if (col %in% names(comp_draws)) {
              val <- comp_draws[[col]][1]
              if (length(val) == 1 && !is.na(val)) {
                comp_cols[[col]] <- unname(val)
              } else {
                comp_cols[[col]] <- unname(comp_cols[[sub("_(LL|UL)$", "", col)]])
              }
            } else {
              comp_cols[[col]] <- unname(comp_cols[[sub("_(LL|UL)$", "", col)]])
            }
          }
        }
      } else {
        for (col in ci_cols) {
          if (grepl("^var_", col) || grepl("_c_", col)) {
            comp_cols[[col]] <- unname(1.0)
          } else if (grepl("_LL$", col)) {
            base_name <- sub("_LL$", "", col)
            comp_cols[[col]] <- unname(comp_cols[[base_name]])
          } else {
            base_name <- sub("_UL$", "", col)
            comp_cols[[col]] <- unname(comp_cols[[base_name]])
          }
        }
      }

      for (cn in names(res)) {
        names(res[[cn]]) <- NULL
      }
      for (cn in names(comp_cols)) {
        names(comp_cols[[cn]]) <- NULL
      }
      res <- rbind(res, comp_cols)
    }

    return(res)
  } else {
    prmse_s_rel <- setNames(coef_df$g[coef_df$dim != "Composite"],
      coef_df$dim[coef_df$dim != "Composite"])
    prmse_s_abs <- setNames(coef_df$phi[coef_df$dim != "Composite"],
      coef_df$dim[coef_df$dim != "Composite"])

    prmse_c_rel <- sapply(dims, function(d) {
      val <- if (!is.null(var[[d]]) && !is.null(var[[d]]$prmse_c_rel)) var[[d]]$prmse_c_rel else NA_real_
      unname(val)
    })
    prmse_c_abs <- sapply(dims, function(d) {
      val <- if (!is.null(var[[d]]) && !is.null(var[[d]]$prmse_c_abs)) var[[d]]$prmse_c_abs else NA_real_
      unname(val)
    })
    prmse_p_rel <- sapply(dims, function(d) {
      val <- if (!is.null(var[[d]]) && !is.null(var[[d]]$prmse_p_rel)) var[[d]]$prmse_p_rel else NA_real_
      unname(val)
    })
    prmse_p_abs <- sapply(dims, function(d) {
      val <- if (!is.null(var[[d]]) && !is.null(var[[d]]$prmse_p_abs)) var[[d]]$prmse_p_abs else NA_real_
      unname(val)
    })
    var_rel <- sapply(dims, function(d) {
      val <- if (!is.null(var[[d]]) && !is.null(var[[d]]$var_rel)) var[[d]]$var_rel else NA_real_
      unname(val)
    })
    var_abs <- sapply(dims, function(d) {
      val <- if (!is.null(var[[d]]) && !is.null(var[[d]]$var_abs)) var[[d]]$var_abs else NA_real_
      unname(val)
    })

    base_cols <- list(
      dim = dims,
      prmse_s_rel = unname(prmse_s_rel[dims]),
      prmse_s_abs = unname(prmse_s_abs[dims]),
      prmse_c_rel = unname(prmse_c_rel[dims]),
      prmse_c_abs = unname(prmse_c_abs[dims])
    )

    if (include_profile) {
      base_cols$prmse_p_rel <- unname(prmse_p_rel[dims])
      base_cols$prmse_p_abs <- unname(prmse_p_abs[dims])
    }

    base_cols$var_rel <- unname(var_rel[dims])
    base_cols$var_abs <- unname(var_abs[dims])

    res <- tibble::as_tibble(base_cols)

    if (!is.null(ci)) {
      warning(
        "Credible intervals are only available for brms backend with posterior estimation. ",
        "For mom backend, consider using a Bayesian approach (backend = 'brms') if CIs are needed.",
        call. = FALSE
      )
    }

    attr(res, "prmse_mv_rel") <- c(mean = NA_real_, LL = NA_real_, UL = NA_real_)
    attr(res, "prmse_mv_abs") <- c(mean = NA_real_, LL = NA_real_, UL = NA_real_)

    if (include_composite) {
      comp_row <- coef_df[coef_df$dim == "Composite", ]
      comp_g <- if (nrow(comp_row) > 0) comp_row$g[1] else NA_real_
      comp_phi <- if (nrow(comp_row) > 0) comp_row$phi[1] else NA_real_

      comp_cols <- list(dim = "Composite")
      comp_cols$prmse_s_rel <- comp_g
      comp_cols$prmse_s_abs <- comp_phi
      comp_cols$prmse_c_rel <- 1.0
      comp_cols$prmse_c_abs <- 1.0
      if (include_profile) {
        comp_cols$prmse_p_rel <- NA_real_
        comp_cols$prmse_p_abs <- NA_real_
      }
      comp_cols$var_rel <- 1.0
      comp_cols$var_abs <- 1.0

      res <- rbind(res, comp_cols)
    }

    return(res)
  }
}

#' Compute VAR Using Haberman Formula (Per Draw)
#'
#' Computes the Value Added Ratio (VAR) for each subscale using the
#' Haberman (2008) formula, which evaluates whether a subscale's own scores
#' better estimate its universe score than the composite does.
#'
#' The formula is:
#' \deqn{VAR(S_i) = \frac{PRMSE(S_i)}{PRMSE(C \rightarrow S_i)}}
#'
#' Where PRMSE(C -> S_i) is the Haberman (2008) correct form:
#' \deqn{PRMSE(C \rightarrow S_i) = \frac{[\mathrm{Cov}(\tau_i, C)]^2}{\mathrm{Var}(\tau_i) \cdot \mathrm{Rel}(C) \cdot \mathrm{Var}(C)}}
#' with \eqn{\mathrm{Cov}(\tau_i, C) = \sum_j w_j \, \mathrm{Cov}(\tau_i, \tau_j)}
#' computed from the universe-score covariance matrix.
#' This is the corrected form of what Vispoel et al. (2023) labeled as
#' "Equation 38" (the printed equation contains a typesetting error in the
#' denominator; see \code{\link{prmse}} for details).
#'
#' @param uni_cov_draws 3D array (n_draws x n_dims x n_dims) of universe score covariances.
#' @param total_rel_cov_draws 3D array (n_draws x n_dims x n_dims) of total observed covariances for relative error.
#' @param total_abs_cov_draws 3D array (n_draws x n_dims x n_dims) of total observed covariances for absolute error.
#' @param dimensions Character vector of dimension names.
#' @param weights Named numeric vector of weights for the composite scores.
#'
#' @return A list containing:
#' \describe{
#'   \item{prmse_s_rel, prmse_s_abs}{Univariate subscale reliabilities.}
#'   \item{prmse_c_rel, prmse_c_abs}{PRMSE when projecting from weighted composite.}
#'   \item{prmse_p_rel, prmse_p_abs}{PRMSE when projecting from entire profile (diagonal slice).}
#'   \item{var_rel, var_abs}{Value Added Ratios (PRMSE_S / PRMSE_C).}
#'   \item{prmse_mv_rel, prmse_mv_abs}{Global multivariate PRMSE (trace formula).}
#' }
#'
#' @references
#' Haberman, S. J. (2008). When can subscores have value? \emph{Journal of
#' Educational and Behavioral Statistics}, 33(2), 204-229.
#'
#' Vispoel, W. P., Lee, H., Hong, H., & Chen, T. (2023). Applying multivariate
#' generalizability theory to psychological assessments. \emph{Psychological
#' Methods}, 28(4), 847-870.
#'
#' @keywords internal
compute_var_haberman_draws <- function(
  uni_cov_draws,   # [n_draws, n_dims, n_dims]
  total_rel_cov_draws, # [n_draws, n_dims, n_dims]
  total_abs_cov_draws, # [n_draws, n_dims, n_dims]
  dimensions,
  weights) {

  n_draws <- dim(uni_cov_draws)[1]
  k <- length(dimensions)

  # Metrics for subscales
  prmse_s_rel_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(prmse_s_rel_draws) <- dimensions
  prmse_s_abs_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(prmse_s_abs_draws) <- dimensions

  prmse_c_rel_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(prmse_c_rel_draws) <- dimensions
  prmse_c_abs_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(prmse_c_abs_draws) <- dimensions

  prmse_p_rel_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(prmse_p_rel_draws) <- dimensions
  prmse_p_abs_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(prmse_p_abs_draws) <- dimensions

  var_rel_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(var_rel_draws) <- dimensions
  var_abs_draws <- matrix(NA_real_, nrow = n_draws, ncol = k)
  colnames(var_abs_draws) <- dimensions

  # Profile-level overall metrics
  prmse_mv_rel <- numeric(n_draws)
  prmse_mv_abs <- numeric(n_draws)

  w <- weights[dimensions]

  for (i in seq_len(n_draws)) {
    Sigma_tau <- uni_cov_draws[i, , ]
    Sigma_rel <- total_rel_cov_draws[i, , ]
    Sigma_abs <- total_abs_cov_draws[i, , ]

    # 1. Global Profile Reliability (Trace Formula)
    # PRMSE_MV = tr(Sigma_tau * Sigma_obs^-1 * Sigma_tau) / tr(Sigma_tau)

    # Pre-calculate common matrix products
    S_rel_inv <- tryCatch(solve(Sigma_rel), error = function(e) matrix(0, k, k))
    S_abs_inv <- tryCatch(solve(Sigma_abs), error = function(e) matrix(0, k, k))

    T_S_rel_T <- Sigma_tau %*% S_rel_inv %*% Sigma_tau
    T_S_abs_T <- Sigma_tau %*% S_abs_inv %*% Sigma_tau

    tr_tau <- sum(diag(Sigma_tau))
    prmse_mv_rel[i] <- if (tr_tau > 0) sum(diag(T_S_rel_T)) / tr_tau else NA_real_
    prmse_mv_abs[i] <- if (tr_tau > 0) sum(diag(T_S_abs_T)) / tr_tau else NA_real_

    # 2. Composite prediction (for current weights)
    var_C_rel <- as.numeric(t(w) %*% Sigma_rel %*% w)
    var_C_abs <- as.numeric(t(w) %*% Sigma_abs %*% w)
    var_tau_C <- as.numeric(t(w) %*% Sigma_tau %*% w)
    Rel_C_rel <- if (var_C_rel > 0) var_tau_C / var_C_rel else NA_real_
    Rel_C_abs <- if (var_C_abs > 0) var_tau_C / var_C_abs else NA_real_
    cov_tau_C <- as.numeric(Sigma_tau %*% w)

    for (d_idx in seq_len(k)) {
      var_tau_d <- Sigma_tau[d_idx, d_idx]

      # 3. Univariate metrics (PRMSE_S = Rel(S))
      rel_S_rel <- if (Sigma_rel[d_idx, d_idx] > 0) var_tau_d / Sigma_rel[d_idx, d_idx] else NA_real_
      rel_S_abs <- if (Sigma_abs[d_idx, d_idx] > 0) var_tau_d / Sigma_abs[d_idx, d_idx] else NA_real_

      prmse_s_rel_draws[i, d_idx] <- pmin(1, pmax(0, rel_S_rel))
      prmse_s_abs_draws[i, d_idx] <- pmin(1, pmax(0, rel_S_abs))

      # 4. Composite-based metrics (PRMSE_C)
      # PRMSE_C = (Cov(tau_d, C))^2 / (Var(tau_d) * Rel_C * Var(C))
      den_c_rel <- var_tau_d * Rel_C_rel * var_C_rel
      den_c_abs <- var_tau_d * Rel_C_abs * var_C_abs

      prmse_c_rel <- if (!is.na(den_c_rel) && den_c_rel > 0) (cov_tau_C[d_idx]^2) / den_c_rel else NA_real_
      prmse_c_abs <- if (!is.na(den_c_abs) && den_c_abs > 0) (cov_tau_C[d_idx]^2) / den_c_abs else NA_real_

      prmse_c_rel_draws[i, d_idx] <- pmin(1, pmax(0, prmse_c_rel))
      prmse_c_abs_draws[i, d_idx] <- pmin(1, pmax(0, prmse_c_abs))

      # 5. Profile-based metrics (PRMSE_P = Diagonal Slice)
      # PRMSE_P_i = [T Sigma_obs^-1 T]_ii / [T]_ii
      prmse_p_rel <- if (var_tau_d > 0) T_S_rel_T[d_idx, d_idx] / var_tau_d else NA_real_
      prmse_p_abs <- if (var_tau_d > 0) T_S_abs_T[d_idx, d_idx] / var_tau_d else NA_real_

      prmse_p_rel_draws[i, d_idx] <- pmin(1, pmax(0, prmse_p_rel))
      prmse_p_abs_draws[i, d_idx] <- pmin(1, pmax(0, prmse_p_abs))

      # 6. Value Added Ratio (VAR = PRMSE_S / PRMSE_C)
      # Use clipped values for consistency with stored PRMSE values
      var_rel_draws[i, d_idx] <- if (!is.na(prmse_c_rel) && prmse_c_rel > 0) {
        pmin(1, pmax(0, rel_S_rel)) / pmin(1, pmax(0, prmse_c_rel))
      } else {
        NA_real_
      }

      var_abs_draws[i, d_idx] <- if (!is.na(prmse_c_abs) && prmse_c_abs > 0) {
        pmin(1, pmax(0, rel_S_abs)) / pmin(1, pmax(0, prmse_c_abs))
      } else {
        NA_real_
      }
    }
  }

  list(
    prmse_s_rel = prmse_s_rel_draws,
    prmse_s_abs = prmse_s_abs_draws,
    prmse_c_rel = prmse_c_rel_draws,
    prmse_c_abs = prmse_c_abs_draws,
    prmse_p_rel = prmse_p_rel_draws,
    prmse_p_abs = prmse_p_abs_draws,
    var_rel = var_rel_draws,
    var_abs = var_abs_draws,
    prmse_mv_rel = prmse_mv_rel,
    prmse_mv_abs = prmse_mv_abs
  )
}

