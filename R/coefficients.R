#' Generalizability and Dependability Coefficients
#'
#' Functions for computing G (generalizability) and D (dependability) coefficients
#' from variance components in generalizability theory.
#'
#' @name coefficients
#' @keywords internal
#' @importFrom magrittr %>%
NULL

#' Calculate D-Study Variance Components
#'
#' Scales variance components according to D-study sample sizes.
#' For each facet, the variance is divided by the number of levels in the D-study.
#'
#' @param vc G-study variance components tibble.
#' @param n Named list of sample sizes for each facet.
#' @param object Character vector naming the object(s) of measurement.
#' @param aggregation Character vector of facets to aggregate over. When specified,
#' components containing aggregation facets are divided only by the aggregation
#' facet's n (not by the full scale factor).
#' @param n_provided Logical indicating if n was explicitly provided by the user.
#' When FALSE (default), no rescaling is applied.
#' @param residual_is Character string specifying which facets make up the residual.
#' @param facet_n Named vector of facet sample sizes from G-study. Used to fill
#' in missing aggregation facets when they are not specified in `n`.
#' @return A tibble of D-study variance components with columns:
#' \itemize{
#' \item component: Name of variance component
#' \item var_unscaled: Unscaled variance estimate (from G-study)
#' \item pct_unscaled: Percentage of total unscaled variance
#' \item var_scaled: Scaled variance (divided by D-study sample sizes)
#' \item pct_scaled: Percentage of total scaled variance
#' }
#' When aggregation is specified, returns both scaled and unscaled estimates.
#'
#' @keywords internal
calculate_dstudy_variance <- function(vc, n, object, aggregation = NULL, n_provided = FALSE, residual_is = NULL, facet_n = NULL) {
  d_vc <- vc

  object_spec <- parse_specification(object)

  agg_facets <- if (!is.null(aggregation)) parse_specification(aggregation) else NULL
  residual_facets <- if (!is.null(residual_is) && residual_is != "") parse_component_facets(residual_is) else NULL

  has_aggregation <- !is.null(aggregation) && length(aggregation) > 0

  # Store original n for scale_factor computation
  # scale_factor should only include facets that were originally in n (user-specified)
  n_original <- n

  # Fill in missing aggregation facets from facet_n (G-study sample sizes)
  if (has_aggregation && !is.null(facet_n)) {
    for (f in agg_facets) {
      if (!(f %in% names(n)) && f %in% names(facet_n)) {
        n[[f]] <- facet_n[[f]]
      }
    }
  }

  n_facets <- if (n_provided && length(n) > 0) names(n) else character(0)
  n_original_facets <- if (n_provided && length(n_original) > 0) names(n_original) else character(0)

  if (has_aggregation) {
    effective_agg_facets <- unique(c(agg_facets, n_facets))
    additional_agg_facets <- setdiff(agg_facets, n_original_facets)
  } else {
    effective_agg_facets <- n_facets
    additional_agg_facets <- character(0)
  }

  d_vc <- d_vc %>%
    mutate(
      facets_list = purrr::map(component, function(c) {
        if (c == "Residual" && !is.null(residual_facets)) {
          return(residual_facets)
        }
        parse_component_facets(c)
      }),
      is_object = component %in% object_spec,
      is_residual = (component == "Residual"),
      is_interaction = purrr::map_lgl(facets_list, ~ length(.x) > 1),
      scale_factor = purrr::map_dbl(facets_list, compute_scale_factor_from_facets, n = n_original),
      n_facet = purrr::map_dbl(facets_list, function(flist) {
        if (length(flist) == 1 && flist %in% names(n)) {
          n[[flist]]
        } else {
          1
        }
      }),
      has_additional_agg = purrr::map_lgl(component, function(c) {
        if (length(additional_agg_facets) == 0) return(FALSE)
        if (c == "Residual" && !is.null(residual_facets)) {
          return(any(additional_agg_facets %in% residual_facets))
        }
        facets <- parse_component_facets(c)
        any(facets %in% additional_agg_facets)
      }),
      additional_agg_factor = purrr::map_dbl(component, function(c) {
        if (length(additional_agg_facets) == 0) return(1)
        if (c == "Residual" && !is.null(residual_facets)) {
          agg_in_comp <- residual_facets[residual_facets %in% additional_agg_facets]
        } else {
          facets <- parse_component_facets(c)
          agg_in_comp <- facets[facets %in% additional_agg_facets]
        }
        if (length(agg_in_comp) > 0) {
          agg_factor <- 1
          for (f in agg_in_comp) {
            if (f %in% names(n)) agg_factor <- agg_factor * n[[f]]
          }
          return(agg_factor)
        }
        return(1)
      }),
      has_agg = purrr::map_lgl(component, function(c) {
        if (is.null(effective_agg_facets) || length(effective_agg_facets) == 0) return(FALSE)
        if (c == "Residual" && !is.null(residual_facets)) {
          return(any(effective_agg_facets %in% residual_facets))
        }
        facets <- parse_component_facets(c)
        any(facets %in% effective_agg_facets)
      }),
      agg_n = purrr::map_dbl(component, function(c) {
        if (is.null(effective_agg_facets) || length(effective_agg_facets) == 0) return(1)
        if (c == "Residual" && !is.null(residual_facets)) {
          agg_in_comp <- residual_facets[residual_facets %in% effective_agg_facets]
        } else {
          facets <- parse_component_facets(c)
          agg_in_comp <- facets[facets %in% effective_agg_facets]
        }
        if (length(agg_in_comp) > 0) {
          agg_factor <- 1
          for (f in agg_in_comp) {
            if (f %in% names(n)) agg_factor <- agg_factor * n[[f]]
          }
          return(agg_factor)
        }
        return(1)
      })
    )

  do_rescaling <- !is.null(aggregation) || (n_provided && length(n) > 0)

  if (!do_rescaling) {
    d_vc <- d_vc %>%
      mutate(
        var_unscaled = .data$var,
        var_scaled = .data$var,
        pct_unscaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100,
        pct_scaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100
      ) %>%
      select(-dplyr::any_of(c("var", "pct", "facets_list", "is_object", "is_residual", "is_interaction",
        "scale_factor", "n_facet", "has_additional_agg", "additional_agg_factor",
        "has_agg", "agg_n", "error", "se", "lower", "upper", "sd", "Rhat", "Bulk_ESS", "Tail_ESS")))
  } else {
    total_n <- if (length(n) > 0) prod(unlist(n)) else 1

    if (has_aggregation) {
      d_vc <- d_vc %>%
        mutate(
          var_unscaled = .data$var,
          var = case_when(
            is_object ~ .data$var,
            has_additional_agg ~ .data$var / (scale_factor * additional_agg_factor),
            is_residual & has_agg ~ .data$var / (scale_factor * additional_agg_factor),
            is_residual ~ .data$var / total_n,
            has_agg & !is_interaction ~ .data$var / (scale_factor * additional_agg_factor),
            has_agg & is_interaction ~ .data$var / (scale_factor * additional_agg_factor),
            is_interaction ~ .data$var / scale_factor,
            n_facet > 1 ~ .data$var / n_facet,
            TRUE ~ .data$var
          ),
          var_scaled = .data$var,
          pct_scaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100,
          pct_unscaled = (var_unscaled / sum(var_unscaled, na.rm = TRUE)) * 100
        ) %>%
        select(-dplyr::any_of(c("pct", "facets_list", "is_object", "is_residual", "is_interaction",
          "scale_factor", "n_facet", "has_additional_agg", "additional_agg_factor",
          "has_agg", "agg_n", "error", "se", "lower", "upper", "sd", "Rhat", "Bulk_ESS", "Tail_ESS")))
    } else {
      d_vc <- d_vc %>%
        mutate(
          var_unscaled = .data$var,
          var = case_when(
            is_object ~ .data$var,
            has_agg ~ .data$var / scale_factor,
            is_residual ~ .data$var / total_n,
            is_interaction ~ .data$var / scale_factor,
            n_facet > 1 ~ .data$var / n_facet,
            TRUE ~ .data$var
          ),
          var_scaled = .data$var,
          pct_unscaled = (var_unscaled / sum(var_unscaled, na.rm = TRUE)) * 100,
          pct_scaled = (.data$var / sum(.data$var, na.rm = TRUE)) * 100
        ) %>%
        select(-dplyr::any_of(c("pct", "facets_list", "is_object", "is_residual", "is_interaction",
          "scale_factor", "n_facet", "has_additional_agg", "additional_agg_factor",
          "has_agg", "agg_n", "error", "se", "lower", "upper", "sd", "Rhat", "Bulk_ESS", "Tail_ESS")))
    }
  }

  # Ensure dim column comes after component for multivariate models
  if ("dim" %in% names(d_vc)) {
    d_vc <- d_vc %>%
      dplyr::relocate(dim, .after = component)
  }

  d_vc
}

#' Calculate Divided Variance Components
#'
#' Divides each variance component by the sample sizes of non-object facets
#' that appear in that component. This provides an alternative scaling where:
#' - Object component is NOT scaled (same as standard D-study)
#' - For non-object components, divides only by n for non-object facets in component
#' - For Residual, divides by product of all non-object facet sample sizes
#'
#' This differs from standard D-study variance scaling in how interactions are handled.
#' Standard D-study scales interactions by all facets in the component, while this
#' method scales only by the non-object facets.
#'
#' @param vc Variance components tibble (unscaled, from G-study)
#' @param n Named list of sample sizes for each facet
#' @param object Specification for object of measurement
#' @param residual_is Character string specifying which facets make up the residual
#' @return Tibble with additional 'var_divided' column
#'
#' @keywords internal
calculate_divided_variance <- function(vc, n, object, residual_is = NULL) {
  if (is.null(n) || length(n) == 0) {
    d_vc <- vc %>%
      mutate(var_divided = var)
    return(d_vc)
  }

  object_spec <- parse_specification(object)

  all_facets <- names(n)
  non_object_facets <- setdiff(all_facets, object_spec)

  residual_facets <- if (!is.null(residual_is) && residual_is != "") {
    parse_component_facets(residual_is)
  } else {
    NULL
  }

  d_vc <- vc %>%
    mutate(
      comp_facets = purrr::map(component, function(c) {
        if (c == "Residual") {
          if (!is.null(residual_facets)) {
            residual_facets
          } else {
            non_object_facets
          }
        } else {
          parse_component_facets(c)
        }
      }),
      is_object = component %in% object_spec,
      divisor = purrr::map_dbl(comp_facets, function(facets) {
        if (length(facets) == 0) return(1)

        facets_to_use <- facets[facets %in% non_object_facets]
        if (length(facets_to_use) == 0) return(1)

        div <- 1
        for (f in facets_to_use) {
          if (f %in% names(n)) {
            div <- div * n[[f]]
          }
        }
        div
      }),
      var_divided = ifelse(is_object, var, var / divisor)
    ) %>%
    select(-comp_facets, -is_object, -divisor)

  d_vc
}

#' Calculate Coefficients Using Divided Variance Components
#'
#' Computes coefficients using the "divided" variance estimates where each
#' component is divided only by the sample sizes of non-object facets.
#'
#' @param vc_divided Variance components tibble with 'var_divided' column
#' @param object Specification for object of measurement
#' @param error Specification for error components (optional)
#' @param aggregation Aggregation facets (optional)
#' @param residual_is Residual composition specification
#' @param universe Universe specification (optional)
#' @return Data frame with coefficient estimates
#'
#' @keywords internal
calculate_divided_coefficients <- function(vc_divided, object, error = NULL,
  aggregation = NULL, residual_is = NULL, universe = NULL,
  cut_score = NULL, mu_y = NULL) {
  vc_for_calc <- vc_divided %>%
    mutate(var = var_divided) %>%
    select(-var_divided)

  calculate_coefficients(vc_for_calc, n = NULL, object, error, aggregation, residual_is, universe, cut_score, mu_y)
}

#' Calculate G and D Coefficients
#'
#' Computes generalizability (G) and dependability (D/Phi) coefficients
#' from D-study variance components.
#'
#' @param vc D-study variance components tibble.
#' @param n Named list of sample sizes for each facet (required for aggregation).
#' @param object Specification for object of measurement. Can be:
#'   - NULL (default): uses first component in vc
#'   - A character string: "p"
#'   - A character vector: c("p", "p:d")
#'   - A formula: obj ~ p + p:d
#' @param error Specification for error components. Can be:
#'   - NULL (default): all components not in object specification
#'   - A character string: "p:i"
#'   - A character vector: c("p:i", "p:r")
#'   - A formula: ~ p:i + p:r
#'   Note: If the same facet is specified in both `object` and `error`, an error is raised.
#' @param aggregation Character vector of facets to aggregate over.
#'   Components containing these facets will be divided by the sample size.
#'   Default is NULL (no aggregation).
#'   Note: If the same facet is specified in both `object` and `aggregation`, an error is raised.
#' @param residual_is Character string specifying which facets make up the residual.
#'   For example, "p:i:r" means the residual represents the p:i:r interaction.
#'   Default is NULL (residual is not rescaled).
#' @return A data frame with columns:
#' \item{uni}{Universe score variance. Represents the true individual differences
#' in the trait being measured. It reflects the systematic variance attributable
#' to the object of measurement (e.g., persons) after accounting for measurement
#' error. In classical test theory terms, this corresponds to the "true score"
#' variance. A larger universe score variance indicates greater heterogeneity
#' among objects of measurement, which generally improves the precision of
#' relative comparisons.}
#' \item{sigma2_delta}{Relative error variance. Error variance affecting relative
#' decisions (rankings). Includes only variance components that interact with the
#' object of measurement (i.e., components that contain the object facet). For
#' example, in a person x rater design with person as object, this includes
#' person:rater and residual components.}
#' \item{sigma2_delta_abs}{Absolute error variance. Error variance affecting absolute
#' decisions. Includes ALL variance components except the object of measurement.
#' For example, in a person x rater design with person as object, this includes
#' rater, person:rater, and residual components.}
#' \item{g}{Generalizability coefficient. Quantifies the reliability of
#' measurements for relative decisions—situations where the focus is on
#' rank-ordering or comparing objects of measurement (e.g., identifying high
#' vs. low performers). Calculated as: g = uni / (uni + sigma2_delta).
#' The G coefficient ranges from 0 to 1, with higher values indicating better
#' reliability. Values >= 0.70 are typically considered acceptable for group-level
#' decisions, while >= 0.80 is preferred for individual-level decisions.}
#' \item{phi}{Dependability coefficient (Phi coefficient). Quantifies the
#' reliability of measurements for absolute decisions—situations where the
#' focus is on the absolute level of a measurement (e.g., determining if a
#' score exceeds a cutoff). Calculated as: phi = uni / (uni + sigma2_delta_abs).
#' The Phi coefficient is always <= G coefficient because absolute error includes
#' all error sources, not just those affecting relative rankings.}
#' \item{sem_rel}{Standard error of measurement for relative decisions.
#' Calculated as sqrt(sigma2_delta). Indicates the typical error when
#' interpreting relative standing or rank-order among objects of measurement.}
#' \item{sem_abs}{Standard error of measurement for absolute decisions.
#' Calculated as sqrt(sigma2_delta_abs). Indicates the typical error when
#' interpreting absolute score levels against criterion-referenced standards.}
#'
#' @keywords internal
calculate_coefficients <- function(vc, n = NULL, object = NULL, error = NULL,
  aggregation = NULL, residual_is = NULL, universe = NULL,
  cut_score = NULL, mu_y = NULL) {
  if (is.null(object) && is.null(error)) {
    default_object <- vc$component[1]
    warning(
      "No object specified. Using first component '", default_object,
      "' as the object of measurement.",
      call. = FALSE
    )
    object_spec <- default_object
  } else if (is.null(object)) {
    object_spec <- character(0)
  } else {
    object_spec <- parse_specification(object)
  }

  # Process universe specification
  # Default: universe = object only
  if (is.null(universe) || length(universe) == 0) {
    universe_spec <- object_spec
  } else {
    universe_spec <- parse_specification(universe)
  }

  if (!is.null(error)) {
    error_spec <- parse_specification(error)
  } else {
    error_spec <- NULL
  }

  if (!is.null(aggregation)) {
    agg_facets <- parse_specification(aggregation)
  } else {
    agg_facets <- NULL
  }

  if (!is.null(error_spec) && length(error_spec) > 0) {
    overlap <- intersect(universe_spec, error_spec)
    if (length(overlap) > 0) {
      stop(
        "Component(s) '", paste(overlap, collapse = "', '"), "' specified as both ",
        "in the universe and in error. ",
        "The same component cannot be in both the universe and error.",
        call. = FALSE
      )
    }
  }

  has_dim <- "dim" %in% names(vc)
  dims <- if (has_dim) unique(vc$dim) else "response"

  # Check if we have unscaled variance (from calculate_dstudy_variance with aggregation)
  has_var_unscaled <- "var_unscaled" %in% names(vc)

  results <- lapply(dims, function(d) {
    if (has_dim) {
      vc_dim <- vc[vc$dim == d, , drop = FALSE]
    } else {
      vc_dim <- vc
    }

    vc_modified <- vc_dim

    # Calculate universe score variance
    # Use unscaled variance if available for proper scaling
    uni <- 0
    for (comp in universe_spec) {
      # Check if component exists (might be named differently)
      actual_comp <- comp
      if (!(comp %in% vc_modified$component)) {
        # Check if this is a residual specification
        if (!is.null(residual_is) && comp == residual_is && "Residual" %in% vc_modified$component) {
          actual_comp <- "Residual"
        }
      }

      if (actual_comp %in% vc_modified$component) {
        # Use unscaled variance if available, otherwise use var
        if (has_var_unscaled) {
          comp_var <- vc_modified$var_unscaled[vc_modified$component == actual_comp]
        } else {
          comp_var <- vc_modified$var[vc_modified$component == actual_comp]
        }

        if (comp %in% object_spec) {
          # Object component: NOT scaled
          uni <- uni + comp_var
        } else {
          # Non-object universe component: scaled by sample sizes
          # Need to determine scale factor based on component facets
          scale_factor <- get_universe_scale_factor_for_component(comp, n, agg_facets, residual_is, actual_comp)
          uni <- uni + comp_var / scale_factor
        }
      }
    }

    # For relative error, exclude universe components (not just object)
    sigma2_delta <- compute_relative_error(vc_modified, universe_spec, error_spec, agg_facets, residual_is)

    # For absolute error, exclude universe components (not just object)
    sigma2_delta_abs <- compute_absolute_error(vc_modified, universe_spec, error_spec, agg_facets, residual_is)

    g <- uni / (uni + sigma2_delta)

    phi <- uni / (uni + sigma2_delta_abs)

    phi_cut <- NA_real_
    if (!is.null(cut_score) && !is.null(mu_y)) {
      mu_y_val <- if (is.list(mu_y)) mu_y[[d]] else mu_y
      if (!is.null(mu_y_val) && !is.na(mu_y_val)) {
        adjustment <- (mu_y_val - cut_score)^2
        phi_cut <- (uni + adjustment) / (uni + sigma2_delta_abs + adjustment)
      }
    }

    sem_rel <- sqrt(sigma2_delta)
    sem_abs <- sqrt(sigma2_delta_abs)

    if (has_dim) {
      result_df <- data.frame(
        dim = d,
        uni = uni,
        sigma2_delta = sigma2_delta,
        sigma2_delta_abs = sigma2_delta_abs,
        g = g,
        phi = phi,
        sem_rel = sem_rel,
        sem_abs = sem_abs,
        stringsAsFactors = FALSE
      )
    } else {
      result_df <- data.frame(
        uni = uni,
        sigma2_delta = sigma2_delta,
        sigma2_delta_abs = sigma2_delta_abs,
        g = g,
        phi = phi,
        sem_rel = sem_rel,
        sem_abs = sem_abs
      )
    }

    if (!is.null(cut_score)) {
      result_df$phi_cut <- phi_cut
    }

    result_df
  })

  result_df <- do.call(rbind, results)

  result_df
}

#' Build Variance-Covariance Matrix for a Component
#'
#' Constructs a variance-covariance matrix for a specific variance component
#' in a multivariate design. The diagonal elements are the variances for each
#' dimension, and off-diagonal elements are covariances (when available).
#'
#' @param vc A variance components tibble with columns: component, dim, var
#' @param component Character string naming the component (e.g., "Person", "Residual")
#' @param dimensions Character vector of dimension names
#' @param correlations List containing covariance information from the G-study:
#'   - `residual_cov`: Tibble of residual covariances
#'   - `random_effect_cov`: Named list of random effect covariances per facet
#'
#' @return A symmetric variance-covariance matrix with dimensions as row/column names
#'
#' @keywords internal
build_component_covariance_matrix <- function(vc, component, dimensions, correlations = NULL, scale_factor = 1) {
  n_dim <- length(dimensions)
  cov_matrix <- matrix(0, n_dim, n_dim)
  rownames(cov_matrix) <- dimensions
  colnames(cov_matrix) <- dimensions

  for (i in seq_along(dimensions)) {
    dim_var <- vc$var[vc$component == component & vc$dim == dimensions[i]]
    if (length(dim_var) > 0) {
      cov_matrix[i, i] <- dim_var[1] / scale_factor
    }
  }

  if (!is.null(correlations) && component == "Residual") {
    if (!is.null(correlations$residual_cov) && nrow(correlations$residual_cov) > 0) {
      for (i in seq_len(nrow(correlations$residual_cov))) {
        row <- correlations$residual_cov[i, ]
        idx1 <- match(row$dim1, dimensions)
        idx2 <- match(row$dim2, dimensions)
        if (!is.na(idx1) && !is.na(idx2)) {
          # Scale the covariance by the scale_factor
          scaled_cov <- row$estimate / scale_factor
          cov_matrix[idx1, idx2] <- scaled_cov
          cov_matrix[idx2, idx1] <- scaled_cov
        }
      }
    }
  } else if (!is.null(correlations) && !is.null(correlations$random_effect_cov)) {
    facet_name <- component
    if (facet_name %in% names(correlations$random_effect_cov)) {
      facet_cov <- correlations$random_effect_cov[[facet_name]]
      if (!is.null(facet_cov) && nrow(facet_cov) > 0) {
        for (i in seq_len(nrow(facet_cov))) {
          row <- facet_cov[i, ]
          idx1 <- match(row$dim1, dimensions)
          idx2 <- match(row$dim2, dimensions)
          if (!is.na(idx1) && !is.na(idx2)) {
            # Scale the covariance by the scale_factor
            scaled_cov <- row$estimate / scale_factor
            cov_matrix[idx1, idx2] <- scaled_cov
            cov_matrix[idx2, idx1] <- scaled_cov
          }
        }
      }
    }
  }

  cov_matrix
}

#' Compute Weighted Variance from Covariance Matrix
#'
#' Computes the variance of a weighted composite using the quadratic form
#' w' * Sigma * w, where w is the weight vector and Sigma is the
#' variance-covariance matrix.
#'
#' @param cov_matrix A symmetric variance-covariance matrix
#' @param weights Numeric vector of weights (one per dimension)
#'
#' @return The weighted variance as a numeric scalar
#'
#' @keywords internal
compute_weighted_variance <- function(cov_matrix, weights) {
  as.numeric(t(weights) %*% cov_matrix %*% weights)
}

#' Calculate Composite G and Phi Coefficients
#'
#' Computes composite reliability coefficients for multivariate designs by
#' combining variance-covariance information across dimensions using weights.
#'
#' @param vc A variance components tibble (scaled for D-study)
#' @param n Named list of sample sizes for each facet
#' @param weights Named numeric vector of weights (names = dimensions)
#' @param object The object of measurement specification
#' @param error Optional error component specification
#' @param aggregation Optional aggregation specification
#' @param residual_is Optional residual specification
#' @param universe Optional universe specification
#' @param correlations List of covariance information from the G-study
#' @param cut_score Optional cut-score for phi-cut calculation
#' @param mu_y Optional grand mean(s) for phi-cut calculation
#'
#' @return A data.frame with one row containing:
#'   - dim: "Composite"
#'   - uni: Universe score variance (composite)
#'   - sigma2_delta: Relative error variance (composite)
#'   - sigma2_delta_abs: Absolute error variance (composite)
#'   - g: G coefficient (composite)
#'   - phi: Phi coefficient (composite)
#'   - sem_rel: Standard error of measurement (relative)
#'   - sem_abs: Standard error of measurement (absolute)
#'   - phi_cut: Phi-cut coefficient (if cut_score provided)
#'
#' @keywords internal
calculate_composite_coefficients <- function(vc, n, weights, object, error = NULL, aggregation = NULL, residual_is = NULL, universe = NULL, correlations = NULL, cut_score = NULL, mu_y = NULL, gstudy_data = NULL, dimension_var = NULL) {
  dimensions <- names(weights)

  if (is.null(object)) {
    object_spec <- vc$component[1]
  } else {
    object_spec <- parse_specification(object)
  }

  if (is.null(universe) || length(universe) == 0) {
    universe_spec <- object_spec
  } else {
    universe_spec <- parse_specification(universe)
  }

  if (!is.null(error)) {
    error_spec <- parse_specification(error)
  } else {
    error_spec <- NULL
  }

  components <- unique(vc$component)

  # Universe score variance must use G-study (unscaled) variance components.
  # When vc contains d_vc (D-study scaled), we temporarily swap var/unscaled.
  vc_uni <- vc
  if ("var_unscaled" %in% names(vc)) {
    vc_uni$var <- vc$var_unscaled
  }

  uni_composite <- 0
  for (comp in universe_spec) {
    if (comp %in% components) {
      cov_mat <- build_component_covariance_matrix(vc_uni, comp, dimensions, correlations, scale_factor = 1)
      uni_composite <- uni_composite + compute_weighted_variance(cov_mat, weights)
    }
  }

  sigma2_delta_composite <- compute_composite_error_variance(
    vc, dimensions, weights, correlations, universe_spec, error_spec, aggregation, residual_is, "relative", n = n, object_spec = object_spec
  )

  sigma2_delta_abs_composite <- compute_composite_error_variance(
    vc, dimensions, weights, correlations, universe_spec, error_spec, aggregation, residual_is, "absolute", n = n, object_spec = object_spec
  )

  g_composite <- uni_composite / (uni_composite + sigma2_delta_composite)
  phi_composite <- uni_composite / (uni_composite + sigma2_delta_abs_composite)
  sem_rel_composite <- sqrt(sigma2_delta_composite)
  sem_abs_composite <- sqrt(sigma2_delta_abs_composite)

  var_rel <- NA_real_
  var_abs <- NA_real_
  per_subscale_var <- NULL

  if ("dim" %in% names(vc) && length(dimensions) > 1) {
    # Universe score variance uses G-study (unscaled) variance
    uni_var_col <- if ("var_unscaled" %in% names(vc)) "var_unscaled" else "var"
    # Error variance uses D-study (scaled) variance
    err_var_col <- if ("var_scaled" %in% names(vc)) "var_scaled" else "var"

    dim_g <- vapply(dimensions, function(d) {
      vc_d <- vc[vc$dim == d, , drop = FALSE]
      uni_d <- sum(vc_d[[uni_var_col]][vc_d$component %in% universe_spec], na.rm = TRUE)
      rel_comps <- vc_d$component[sapply(vc_d$component, function(c) {
        if (c %in% universe_spec) return(FALSE)
        if (c == "Residual") return(TRUE)
        facets <- parse_component_facets(c)
        any(universe_spec %in% facets)
      })]
      err_rel_d <- sum(vc_d[[err_var_col]][vc_d$component %in% rel_comps], na.rm = TRUE)
      g_d <- if ((uni_d + err_rel_d) > 0) uni_d / (uni_d + err_rel_d) else NA_real_
      g_d
    }, numeric(1))
    names(dim_g) <- dimensions

    dim_phi <- vapply(dimensions, function(d) {
      vc_d <- vc[vc$dim == d, , drop = FALSE]
      uni_d <- sum(vc_d[[uni_var_col]][vc_d$component %in% universe_spec], na.rm = TRUE)
      all_comps <- vc_d$component[!(vc_d$component %in% universe_spec)]
      err_abs_d <- sum(vc_d[[err_var_col]][vc_d$component %in% all_comps], na.rm = TRUE)
      phi_d <- if ((uni_d + err_abs_d) > 0) uni_d / (uni_d + err_abs_d) else NA_real_
      phi_d
    }, numeric(1))
    names(dim_phi) <- dimensions

    # Build model-based covariance matrices from variance components
    # This ensures consistency with the draws-based path and Equation 38's
    # requirement for disattenuated (universe-score) covariances
    k <- length(dimensions)

    # Sigma_tau: universe-score covariance matrix from model VCs
    # Uses G-study (unscaled) variance - vc_uni created earlier has unscaled values
    Sigma_tau <- matrix(0, k, k, dimnames = list(dimensions, dimensions))
    for (comp in universe_spec) {
      if (comp %in% components) {
        cov_mat <- build_component_covariance_matrix(
          vc_uni, comp, dimensions, correlations, scale_factor = 1
        )
        Sigma_tau <- Sigma_tau + cov_mat
      }
    }

    # Sigma_obs_rel: model-based observed covariance matrix (relative decisions)
    # = universe + relative error components (scaled by D-study sample sizes)
    Sigma_obs_rel <- Sigma_tau
    # Sigma_obs_abs: model-based observed covariance matrix (absolute decisions)
    # = universe + all error components (scaled by D-study sample sizes)
    Sigma_obs_abs <- Sigma_tau

    for (comp in components) {
      if (comp %in% universe_spec) next
      # vc$var already contains D-study scaled variance, so use scale_factor = 1
      cov_mat <- build_component_covariance_matrix(
        vc, comp, dimensions, correlations, scale_factor = 1
      )
      # All non-universe components contribute to absolute error
      Sigma_obs_abs <- Sigma_obs_abs + cov_mat
      # Only relative error components contribute to relative error
      is_rel <- comp == "Residual" || {
        facets <- parse_component_facets(comp)
        any(universe_spec %in% facets)
      }
      if (is_rel) {
        Sigma_obs_rel <- Sigma_obs_rel + cov_mat
      }
    }

    # Compute PRMSE_P for all dimensions (profile projection)
    # PRMSE_P[d] = [Sigma_tau * solve(Sigma_obs) * Sigma_tau][d,d] / Sigma_tau[d,d]
    prmse_p_rel_all <- rep(NA_real_, k)
    prmse_p_abs_all <- rep(NA_real_, k)
    names(prmse_p_rel_all) <- dimensions
    names(prmse_p_abs_all) <- dimensions

    S_rel_inv <- tryCatch(solve(Sigma_obs_rel), error = function(e) NULL)
    S_abs_inv <- tryCatch(solve(Sigma_obs_abs), error = function(e) NULL)

    if (!is.null(S_rel_inv)) {
      T_S_rel_T <- Sigma_tau %*% S_rel_inv %*% Sigma_tau
      for (i in seq_along(dimensions)) {
        d <- dimensions[i]
        var_tau_d <- Sigma_tau[i, i]
        if (is.finite(T_S_rel_T[i, i]) && var_tau_d > 0) {
          prmse_p_rel_all[d] <- T_S_rel_T[i, i] / var_tau_d
        }
      }
    }

    if (!is.null(S_abs_inv)) {
      T_S_abs_T <- Sigma_tau %*% S_abs_inv %*% Sigma_tau
      for (i in seq_along(dimensions)) {
        d <- dimensions[i]
        var_tau_d <- Sigma_tau[i, i]
        if (is.finite(T_S_abs_T[i, i]) && var_tau_d > 0) {
          prmse_p_abs_all[d] <- T_S_abs_T[i, i] / var_tau_d
        }
      }
    }

    # Compute PRMSE_C for all dimensions using model-based Sigma_tau %*% w
    # PRMSE_C = [Cov(tau_d, C)]^2 / [Var(tau_d) * Rel(C) * Var(C)]
    # This is equivalent to Equation 38 when expanded, using universe-score
    # (disattenuated) covariances as the paper specifies.
    w <- weights[dimensions]
    cov_tau_C <- as.numeric(Sigma_tau %*% w)  # Cov(tau_d, C) for each d
    var_C_rel <- as.numeric(t(w) %*% Sigma_obs_rel %*% w)
    var_C_abs <- as.numeric(t(w) %*% Sigma_obs_abs %*% w)
    var_tau_C <- as.numeric(t(w) %*% Sigma_tau %*% w)
    Rel_C_rel <- if (var_C_rel > 0) var_tau_C / var_C_rel else NA_real_
    Rel_C_abs <- if (var_C_abs > 0) var_tau_C / var_C_abs else NA_real_

    per_subscale_var <- lapply(dimensions, function(d) {
      d_idx <- which(dimensions == d)
      var_tau_d <- Sigma_tau[d_idx, d_idx]

      # PRMSE_C numerator: (Cov(tau_d, C))^2
      # where Cov(tau_d, C) = [Sigma_tau %*% w][d]
      #   = sigma^2_tau_d * w_d + sum_{j != d} sigma_tau_{d,j} * w_j
      # This matches Equation 38's numerator using disattenuated covariances
      num_c <- cov_tau_C[d_idx]^2

      # PRMSE_C denominator: Var(tau_d) * Rel(C) * Var(C)
      den_rel <- var_tau_d * Rel_C_rel * var_C_rel
      den_abs <- var_tau_d * Rel_C_abs * var_C_abs

      prmse_c_rel <- if (!is.na(den_rel) && is.finite(den_rel) && den_rel > 0) {
        num_c / den_rel
      } else NA_real_

      prmse_c_abs <- if (!is.na(den_abs) && is.finite(den_abs) && den_abs > 0) {
        num_c / den_abs
      } else NA_real_

      var_rel_d <- if (!is.na(prmse_c_rel) && prmse_c_rel > 0) dim_g[d] / prmse_c_rel else NA_real_
      var_abs_d <- if (!is.na(prmse_c_abs) && prmse_c_abs > 0) dim_phi[d] / prmse_c_abs else NA_real_

      list(
        prmse_c_rel = prmse_c_rel,
        prmse_c_abs = prmse_c_abs,
        prmse_p_rel = prmse_p_rel_all[d],
        prmse_p_abs = prmse_p_abs_all[d],
        var_rel = var_rel_d,
        var_abs = var_abs_d
      )
    })
    names(per_subscale_var) <- dimensions
  }

  result <- data.frame(
    dim = "Composite",
    uni = uni_composite,
    sigma2_delta = sigma2_delta_composite,
    sigma2_delta_abs = sigma2_delta_abs_composite,
    g = g_composite,
    phi = phi_composite,
    sem_rel = sem_rel_composite,
    sem_abs = sem_abs_composite,
    stringsAsFactors = FALSE
  )

  if (!is.null(cut_score) && !is.null(mu_y)) {
    mu_y_weighted <- sum(weights * sapply(dimensions, function(d) {
      if (is.list(mu_y)) mu_y[[d]] else mu_y
    }))
    adjustment <- (mu_y_weighted - cut_score)^2
    result$phi_cut <- (uni_composite + adjustment) / (uni_composite + sigma2_delta_abs_composite + adjustment)
  }

  list(
    summary = result,
    var_results = per_subscale_var
  )
}

#' Calculate Composite Variance Draws
#'
#' Computes posterior draws of composite variance for each component type
#' using weighted variance-covariance matrices.
#'
#' @param vc_draws Named list (by dimension) of named lists (by component) of variance draws
#' @param cov_draws Named list from extract_covariance_draws()
#' @param weights Named numeric vector of weights (names = dimensions)
#' @param scale_factor Factor to scale variances and covariances (for D-study scaling)
#'
#' @return Named list (by component) of posterior draws for composite variance
#'
#' @keywords internal
calculate_composite_variance_draws <- function(vc_draws, cov_draws, weights, scale_factor = 1) {
  # Input validation
  if (!is.list(vc_draws) || length(vc_draws) == 0) {
    stop("vc_draws must be a non-empty list", call. = FALSE)
  }
  if (!is.numeric(weights) || length(weights) == 0) {
    stop("weights must be a numeric vector", call. = FALSE)
  }

  dimensions <- names(vc_draws)
  n_dim <- length(dimensions)

  if (!all(names(weights) %in% dimensions)) {
    stop("weights names must match dimensions in vc_draws", call. = FALSE)
  }

  # Check vc_draws structure
  first_dim <- vc_draws[[1]]
  if (!is.list(first_dim) || length(first_dim) == 0) {
    stop("vc_draws must be a nested list structure", call. = FALSE)
  }

  n_draws <- length(first_dim[[1]])

  # Get all component names (from first dimension)
  components <- names(vc_draws[[1]])

  composite_draws <- list()

  for (comp in components) {
    comp_draws <- numeric(n_draws)

    # Get variance draws for each dimension
    var_draws_list <- list()
    for (d in dimensions) {
      var_draws_list[[d]] <- vc_draws[[d]][[comp]] / scale_factor
    }

    # For each posterior draw, compute w' * Sigma * w
    for (draw_idx in 1:n_draws) {
      # Build covariance matrix for this draw
      sigma_mat <- matrix(0, n_dim, n_dim)
      rownames(sigma_mat) <- dimensions
      colnames(sigma_mat) <- dimensions

      for (i in 1:n_dim) {
        sigma_mat[i, i] <- var_draws_list[[dimensions[i]]][draw_idx]
      }

      # Add off-diagonal covariances
      cov_list <- NULL
      if (comp == "Residual") {
        cov_list <- cov_draws$residual
      } else if (comp %in% names(cov_draws$random_effect)) {
        cov_list <- cov_draws$random_effect[[comp]]
      }

      if (!is.null(cov_list)) {
        for (pair_name in names(cov_list)) {
          dims_pair <- strsplit(pair_name, "_")[[1]]
          if (length(dims_pair) == 2 && all(dims_pair %in% dimensions)) {
            idx1 <- match(dims_pair[1], dimensions)
            idx2 <- match(dims_pair[2], dimensions)
            cov_val <- cov_list[[pair_name]][draw_idx] / scale_factor
            sigma_mat[idx1, idx2] <- cov_val
            sigma_mat[idx2, idx1] <- cov_val
          }
        }
      }

      # Compute weighted variance: t(weights) %*% sigma_mat %*% weights
      comp_draws[draw_idx] <- as.numeric(t(weights) %*% sigma_mat %*% weights)
    }

    composite_draws[[comp]] <- comp_draws
  }

  composite_draws
}

#' Calculate D-Study Variance Components with Composite Row
#'
#' Adds composite variance component rows for multivariate models using
#' posterior draws for full uncertainty propagation.
#'
#' @param vc_draws Variance draws (from extract_variance_draws)
#' @param cov_draws Covariance draws (from extract_covariance_draws)
#' @param weights Named vector of weights
#' @param n Sample sizes (named list)
#' @param object Object specification
#' @param n_provided Whether n was provided
#'
#' @return List with:
#'   \describe{
#'     \item{vc_table}{Tibble with variance components including composite rows}
#'     \item{posterior}{List of composite draws}
#'   }
#'
#' @keywords internal
calculate_dstudy_variance_composite <- function(vc_draws, cov_draws, weights, n, object, n_provided) {
  dimensions <- names(weights)
  components <- names(vc_draws[[1]])

  object_spec <- parse_specification(object)

  composite_draws <- list()
  composite_summaries <- list()

  for (comp in components) {
    scale_factor <- compute_component_scale_factor(
      comp, n, object_spec, n_provided
    )

    comp_draws <- calculate_composite_variance_draws(
      vc_draws, cov_draws, weights, scale_factor
    )

    composite_draws[[comp]] <- comp_draws[[comp]]

    var_scaled_mean <- mean(comp_draws[[comp]], na.rm = TRUE)
    var_se <- sd(comp_draws[[comp]], na.rm = TRUE)

    # Unscaled = scaled * scale_factor (reverse the D-study scaling)
    var_unscaled_mean <- var_scaled_mean * scale_factor

    composite_summaries[[comp]] <- list(
      var_scaled = var_scaled_mean,
      var_unscaled = var_unscaled_mean,
      se = var_se
    )
  }

  total_var_scaled <- sum(sapply(composite_summaries, function(x) x$var_scaled))
  total_var_unscaled <- sum(sapply(composite_summaries, function(x) x$var_unscaled))

  composite_rows <- lapply(components, function(comp) {
    data.frame(
      component = comp,
      dim = "Composite",
      var_unscaled = composite_summaries[[comp]]$var_unscaled,
      pct_unscaled = (composite_summaries[[comp]]$var_unscaled / total_var_unscaled) * 100,
      var_scaled = composite_summaries[[comp]]$var_scaled,
      pct_scaled = (composite_summaries[[comp]]$var_scaled / total_var_scaled) * 100,
      stringsAsFactors = FALSE
    )
  })

  list(
    vc_table = do.call(rbind, composite_rows),
    posterior = composite_draws
  )
}

#' Compute Composite Error Variance
#'
#' Calculates the composite relative or absolute error variance by summing
#' weighted variance-covariance matrices for all error components.
#'
#' @param vc A variance components tibble
#' @param dimensions Character vector of dimension names
#' @param weights Named numeric vector of weights
#' @param correlations List of covariance information from the G-study
#' @param universe_spec Character vector of universe components
#' @param error_spec Optional character vector of error components
#' @param aggregation Optional aggregation specification
#' @param residual_is Optional residual specification
#' @param error_type Either "relative" or "absolute"
#' @param scale_factor Factor by which to scale covariances
#'
#' @return The composite error variance as a numeric scalar
#'
#' @keywords internal
compute_composite_error_variance <- function(vc, dimensions, weights, correlations, universe_spec, error_spec, aggregation, residual_is, error_type, n = NULL, object_spec = NULL) {
  components <- unique(vc$component)

  if (!is.null(error_spec)) {
    error_components <- error_spec
  } else {
    if (error_type == "relative") {
      error_components <- vc$component[
        sapply(vc$component, function(comp) {
          if (comp == "Residual") return(TRUE)
          if (comp %in% universe_spec) return(FALSE)
          facets <- parse_component_facets(comp)
          any(universe_spec %in% facets)
        })
      ]
    } else {
      error_components <- setdiff(components, universe_spec)
    }
  }

  total_error <- 0
  for (comp in error_components) {
    if (comp %in% components) {
      # vc$var already contains D-study scaled variance, so use scale_factor = 1
      cov_mat <- build_component_covariance_matrix(vc, comp, dimensions, correlations, scale_factor = 1)
      total_error <- total_error + compute_weighted_variance(cov_mat, weights)
    }
  }

  total_error
}

#' Calculate Composite Coefficients from Posterior Draws
#'
#' Computes composite reliability coefficients from posterior draws of
#' variance components and covariances. This ensures proper uncertainty
#' propagation for multivariate D-studies with brms backend.
#'
#' @param vc_draws Named list (by dimension) of named lists (by component) of variance draws
#' @param cov_draws Named list from extract_covariance_draws()
#' @param weights Named numeric vector of weights (names = dimensions)
#' @param scale_factors Named list of scale factors per component
#' @param object_spec Parsed object specification
#' @param universe_spec Parsed universe specification
#' @param error_spec Parsed error specification (or NULL)
#' @param agg_facets Parsed aggregation specification (or NULL)
#' @param residual_is Residual composition specification
#' @param cut_score Optional cut score for phi-cut
#' @param mu_y Named list of grand means (or single value)
#' @param ci Optional CI specification (character vector)
#' @param probs Probability levels for credible intervals
#'
#' @return List with:
#'   \describe{
#'     \item{summary}{Data frame with coefficient point estimates}
#'     \item{distributions}{Named list of posterior draws for each coefficient}
#'   }
#'
#' @keywords internal
calculate_composite_coefficients_from_draws <- function(vc_draws, cov_draws,
  weights, scale_factors,
  object_spec, universe_spec,
  error_spec, agg_facets,
  residual_is, cut_score = NULL,
  mu_y = NULL, ci = NULL,
  probs = c(0.025, 0.975),
  gstudy_data = NULL,
  dimension_var = NULL) {
  dimensions <- names(weights)
  n_draws <- length(vc_draws[[1]][[1]])
  components <- unique(unlist(lapply(vc_draws, names)))

  error_info <- identify_error_components_for_draws(components, universe_spec, error_spec, agg_facets)

  composite_variance_draws <- list()
  for (comp in components) {
    scale_factor <- scale_factors[[comp]]
    composite_variance_draws[[comp]] <- calculate_composite_variance_draws(
      vc_draws, cov_draws, weights, scale_factor
    )[[comp]]
  }

  # Universe score variance must use G-study (unscaled) variance.
  # composite_variance_draws[[comp]] = t(w) * (vc_draws/sf) * w, already scaled.
  # To get unscaled: multiply by scale_factor to undo the division.
  uni_draws <- rep(0, n_draws)
  for (comp in universe_spec) {
    if (comp %in% components) {
      scale_factor <- scale_factors[[comp]]
      uni_draws <- uni_draws + composite_variance_draws[[comp]] * scale_factor
    }
  }

  sigma2_delta_draws <- rep(0, n_draws)
  for (i in which(error_info$is_relative_error)) {
    comp <- components[i]
    sigma2_delta_draws <- sigma2_delta_draws + composite_variance_draws[[comp]]
  }

  sigma2_delta_abs_draws <- rep(0, n_draws)
  for (i in which(error_info$is_absolute_error)) {
    comp <- components[i]
    sigma2_delta_abs_draws <- sigma2_delta_abs_draws + composite_variance_draws[[comp]]
  }

  g_draws <- uni_draws / (uni_draws + sigma2_delta_draws)
  phi_draws <- uni_draws / (uni_draws + sigma2_delta_abs_draws)
  sem_rel_draws <- sqrt(sigma2_delta_draws)
  sem_abs_draws <- sqrt(sigma2_delta_abs_draws)

  g_draws[is.nan(g_draws) | is.infinite(g_draws)] <- NA
  phi_draws[is.nan(phi_draws) | is.infinite(phi_draws)] <- NA

  # Compute per-dimension G/Phi draws for VAR calculation using Haberman formula
  # Each dim_g_draws[[d]] is a list with g, uni, err_abs
  var_rel_draws <- NULL
  var_abs_draws <- NULL
  var_rel_draws_matrix <- NULL
  var_abs_draws_matrix <- NULL

  if (length(dimensions) > 1 && !is.null(gstudy_data)) {
    dim_g_draws <- lapply(dimensions, function(d) {
      uni_d <- rep(0, n_draws)
      err_rel_d <- rep(0, n_draws)
      err_abs_d <- rep(0, n_draws)
      for (comp in components) {
        sf <- scale_factors[[comp]]
        if (comp %in% universe_spec) {
          # Universe components: use G-study (unscaled) variance
          uni_d <- uni_d + vc_draws[[d]][[comp]]
        } else {
          # Error components: use D-study (scaled) variance
          v <- vc_draws[[d]][[comp]] / sf
          err_abs_d <- err_abs_d + v
          is_rel <- comp == "Residual" || grepl("^Residual", comp) || {
            facets <- parse_component_facets(comp)
            any(universe_spec %in% facets)
          }
          if (is_rel) err_rel_d <- err_rel_d + v
        }
      }
      denom_rel <- uni_d + err_rel_d
      g_d <- ifelse(denom_rel > 0, uni_d / denom_rel, NA_real_)
      denom_abs <- uni_d + err_abs_d
      phi_d <- ifelse(denom_abs > 0, uni_d / denom_abs, NA_real_)
      list(g = g_d, phi = phi_d, uni = uni_d, err_abs = err_abs_d)
    })
    names(dim_g_draws) <- dimensions

    # Construct covariance arrays for Haberman formula
    # uni_cov_draws: universe-score covariance (G-study, UNSCALED)
    # total_rel_cov_draws: observed covariance (D-study, SCALED)
    # total_abs_cov_draws: observed covariance (D-study, SCALED)
    uni_cov_draws <- array(0, dim = c(n_draws, length(dimensions), length(dimensions)))
    total_rel_cov_draws <- array(0, dim = c(n_draws, length(dimensions), length(dimensions)))
    total_abs_cov_draws <- array(0, dim = c(n_draws, length(dimensions), length(dimensions)))

    for (i in seq_along(dimensions)) {
      d1 <- dimensions[i]
      for (j in seq_along(dimensions)) {
        d2 <- dimensions[j]

        for (comp in components) {
          sf <- scale_factors[[comp]]

          # Determine covariance list for off-diagonal elements
          cov_list <- NULL
          if (i != j) {
            pair_name <- paste0(pmin(d1, d2), "_", pmax(d1, d2))
            if (comp == "Residual" || grepl("^Residual", comp)) {
              cov_list <- cov_draws$residual
            } else if (!is.null(cov_draws$random_effect) && comp %in% names(cov_draws$random_effect)) {
              cov_list <- cov_draws$random_effect[[comp]]
            }
          }

          # Universe components: UNSCALED (G-study variance)
          if (comp %in% universe_spec) {
            val_uni <- if (i == j) {
              if (!is.null(vc_draws[[d1]][[comp]])) vc_draws[[d1]][[comp]] else numeric(n_draws)
            } else {
              if (!is.null(cov_list) && pair_name %in% names(cov_list)) {
                cov_list[[pair_name]]
              } else {
                numeric(n_draws)
              }
            }
            uni_cov_draws[, i, j] <- uni_cov_draws[, i, j] + val_uni
          }

          # Observed covariance: SCALED (D-study variance)
          val_obs <- if (i == j) {
            if (!is.null(vc_draws[[d1]][[comp]])) vc_draws[[d1]][[comp]] / sf else numeric(n_draws)
          } else {
            if (!is.null(cov_list) && pair_name %in% names(cov_list)) {
              cov_list[[pair_name]] / sf
            } else {
              numeric(n_draws)
            }
          }

          if (error_info$is_relative_error[which(components == comp)] || comp %in% universe_spec) {
            total_rel_cov_draws[, i, j] <- total_rel_cov_draws[, i, j] + val_obs
          }
          if (error_info$is_absolute_error[which(components == comp)] || comp %in% universe_spec) {
            total_abs_cov_draws[, i, j] <- total_abs_cov_draws[, i, j] + val_obs
          }
        }
      }
    }

    var_result <- compute_var_haberman_draws(
      uni_cov_draws = uni_cov_draws,
      total_rel_cov_draws = total_rel_cov_draws,
      total_abs_cov_draws = total_abs_cov_draws,
      dimensions = dimensions,
      weights = weights
    )

    var_rel_draws_matrix <- var_result$var_rel
    var_abs_draws_matrix <- var_result$var_abs

    var_rel_draws_matrix <- var_result$var_rel
    var_abs_draws_matrix <- var_result$var_abs

    # Summary values (means) for Composite row - these will be NA since VAR is per-subscale
    var_rel_draws <- NA_real_
    var_abs_draws <- NA_real_
  } else if (length(dimensions) > 1) {
    # No data available - cannot compute VAR
    var_rel_draws <- NA_real_
    var_abs_draws <- NA_real_
  }

  phi_cut_draws <- NULL
  if (!is.null(cut_score) && !is.null(mu_y)) {
    mu_y_weighted <- sum(weights * sapply(dimensions, function(d) {
      if (is.list(mu_y)) mu_y[[d]] else mu_y
    }))
    adjustment <- (mu_y_weighted - cut_score)^2
    phi_cut_draws <- (uni_draws + adjustment) / (uni_draws + sigma2_delta_abs_draws + adjustment)
    phi_cut_draws[is.nan(phi_cut_draws) | is.infinite(phi_cut_draws)] <- NA
  }

  summary_df <- data.frame(
    dim = "Composite",
    uni = mean(uni_draws, na.rm = TRUE),
    sigma2_delta = mean(sigma2_delta_draws, na.rm = TRUE),
    sigma2_delta_abs = mean(sigma2_delta_abs_draws, na.rm = TRUE),
    g = mean(g_draws, na.rm = TRUE),
    phi = mean(phi_draws, na.rm = TRUE),
    sem_rel = mean(sem_rel_draws, na.rm = TRUE),
    sem_abs = mean(sem_abs_draws, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (!is.null(cut_score)) {
    summary_df$phi_cut <- mean(phi_cut_draws, na.rm = TRUE)
  }

  if (!is.null(ci)) {
    if ("g" %in% ci) {
      summary_df$g_LL <- quantile(g_draws, probs[1], na.rm = TRUE)
      summary_df$g_UL <- quantile(g_draws, probs[2], na.rm = TRUE)
    }
    if ("phi" %in% ci) {
      summary_df$phi_LL <- quantile(phi_draws, probs[1], na.rm = TRUE)
      summary_df$phi_UL <- quantile(phi_draws, probs[2], na.rm = TRUE)
    }
    if ("phi-cut" %in% ci && !is.null(cut_score)) {
      summary_df$phi_cut_LL <- quantile(phi_cut_draws, probs[1], na.rm = TRUE)
      summary_df$phi_cut_UL <- quantile(phi_cut_draws, probs[2], na.rm = TRUE)
    }
  }

  distributions <- list(
    uni = uni_draws,
    sigma2_delta = sigma2_delta_draws,
    sigma2_delta_abs = sigma2_delta_abs_draws,
    g = g_draws,
    phi = phi_draws,
    sem_rel = sem_rel_draws,
    sem_abs = sem_abs_draws
  )

  if (!is.null(cut_score)) {
    distributions$phi_cut <- phi_cut_draws
  }

  var_results <- NULL
  if (!is.null(var_rel_draws_matrix)) {
    var_results <- list(
      prmse_s_rel_draws = var_result$prmse_s_rel,
      prmse_s_abs_draws = var_result$prmse_s_abs,
      prmse_c_rel_draws = var_result$prmse_c_rel,
      prmse_c_abs_draws = var_result$prmse_c_abs,
      prmse_p_rel_draws = var_result$prmse_p_rel,
      prmse_p_abs_draws = var_result$prmse_p_abs,
      var_rel_draws = var_rel_draws_matrix,
      var_abs_draws = var_abs_draws_matrix,
      prmse_mv_rel_draws = var_result$prmse_mv_rel,
      prmse_mv_abs_draws = var_result$prmse_mv_abs,
      prmse_s_rel = colMeans(var_result$prmse_s_rel, na.rm = TRUE),
      prmse_s_abs = colMeans(var_result$prmse_s_abs, na.rm = TRUE),
      prmse_c_rel = colMeans(var_result$prmse_c_rel, na.rm = TRUE),
      prmse_c_abs = colMeans(var_result$prmse_c_abs, na.rm = TRUE),
      prmse_p_rel = colMeans(var_result$prmse_p_rel, na.rm = TRUE),
      prmse_p_abs = colMeans(var_result$prmse_p_abs, na.rm = TRUE),
      var_rel = colMeans(var_rel_draws_matrix, na.rm = TRUE),
      var_abs = colMeans(var_abs_draws_matrix, na.rm = TRUE),
      prmse_mv_rel = mean(var_result$prmse_mv_rel, na.rm = TRUE),
      prmse_mv_abs = mean(var_result$prmse_mv_abs, na.rm = TRUE),
      prmse_s_rel_LL = apply(var_result$prmse_s_rel, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_s_rel_UL = apply(var_result$prmse_s_rel, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_s_abs_LL = apply(var_result$prmse_s_abs, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_s_abs_UL = apply(var_result$prmse_s_abs, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_c_rel_LL = apply(var_result$prmse_c_rel, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_c_rel_UL = apply(var_result$prmse_c_rel, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_c_abs_LL = apply(var_result$prmse_c_abs, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_c_abs_UL = apply(var_result$prmse_c_abs, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_p_rel_LL = apply(var_result$prmse_p_rel, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_p_rel_UL = apply(var_result$prmse_p_rel, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_p_abs_LL = apply(var_result$prmse_p_abs, 2, quantile, probs = probs[1], na.rm = TRUE),
      prmse_p_abs_UL = apply(var_result$prmse_p_abs, 2, quantile, probs = probs[2], na.rm = TRUE),
      var_rel_LL = apply(var_rel_draws_matrix, 2, quantile, probs = probs[1], na.rm = TRUE),
      var_rel_UL = apply(var_rel_draws_matrix, 2, quantile, probs = probs[2], na.rm = TRUE),
      var_abs_LL = apply(var_abs_draws_matrix, 2, quantile, probs = probs[1], na.rm = TRUE),
      var_abs_UL = apply(var_abs_draws_matrix, 2, quantile, probs = probs[2], na.rm = TRUE),
      prmse_mv_rel_LL = quantile(var_result$prmse_mv_rel, probs = probs[1], na.rm = TRUE),
      prmse_mv_rel_UL = quantile(var_result$prmse_mv_rel, probs = probs[2], na.rm = TRUE),
      prmse_mv_abs_LL = quantile(var_result$prmse_mv_abs, probs = probs[1], na.rm = TRUE),
    prmse_mv_abs_UL = quantile(var_result$prmse_mv_abs, probs = probs[2], na.rm = TRUE)
    )
  }

  list(
    summary = summary_df,
    distributions = distributions,
    var_results = var_results
  )
}

#' Compute Relative Error Variance
#'
#' Calculates the relative error variance for a D-study design.
#' Relative error includes variance components that affect relative rankings
#' (interactions with the universe components).
#'
#' @param vc Variance components tibble.
#' @param universe_spec Character vector of universe component names.
#' @param error_spec Character vector of error component names, or NULL for default.
#' @return The relative error variance.
#'
#' @keywords internal
compute_relative_error <- function(vc, universe_spec, error_spec = NULL, agg_facets = NULL, residual_is = NULL) {
  if (!is.null(error_spec)) {
    vc_filtered <- vc %>%
      filter(component %in% error_spec)

    sum(vc_filtered$var, na.rm = TRUE)
  } else {
    residual_facets <- if (!is.null(residual_is)) parse_component_facets(residual_is) else NULL

    vc_filtered <- vc %>%
      mutate(
        is_universe = component %in% universe_spec,
        is_residual = (component == "Residual"),
        is_interaction = purrr::map_lgl(component, is_interaction),
        has_universe = purrr::map_lgl(component, function(c) {
          facets <- parse_component_facets(c)
          any(universe_spec %in% facets)
        }),
        is_agg_main = if (!is.null(agg_facets)) {
          purrr::map_lgl(component, function(c) {
            if (c == "Residual") return(FALSE)
            facets <- parse_component_facets(c)
            length(facets) == 1 && facets %in% agg_facets
          })
        } else {
          FALSE
        }
      ) %>%
      filter(
        !is_universe,
        !is_agg_main,
        is_residual | has_universe
      )

    sum(vc_filtered$var, na.rm = TRUE)
  }
}

#' Compute Absolute Error Variance
#'
#' Calculates the absolute error variance for a D-study design.
#' Absolute error includes all variance components that affect absolute scores
#' (all components except the universe components).
#'
#' @param vc Variance components tibble.
#' @param universe_spec Character vector of universe component names.
#' @param error_spec Character vector of error component names, or NULL for default.
#' @return The absolute error variance.
#'
#' @keywords internal
compute_absolute_error <- function(vc, universe_spec, error_spec = NULL, agg_facets = NULL, residual_is = NULL) {
  if (!is.null(error_spec)) {
    vc_filtered <- vc %>%
      filter(component %in% error_spec)

    sum(vc_filtered$var, na.rm = TRUE)
  } else {
    vc_filtered <- vc %>%
      mutate(
        is_agg_main = if (!is.null(agg_facets)) {
          purrr::map_lgl(component, function(c) {
            if (c == "Residual") return(FALSE)
            facets <- parse_component_facets(c)
            length(facets) == 1 && facets %in% agg_facets
          })
        } else {
          FALSE
        }
      ) %>%
      filter(
        !(component %in% universe_spec),
        !is_agg_main
      ) %>%
      select(-is_agg_main)

    sum(vc_filtered$var, na.rm = TRUE)
  }
}

#' Compute Generalizability Coefficient
#'
#' Calculates the generalizability coefficient (G coefficient), which represents
#' the proportion of observed variance attributable to universe score variance
#' for relative decisions.
#'
#' @param variance_components A tibble of variance components.
#' @param n A named list specifying the number of levels for each facet.
#' @param object Specification for object of measurement (character, vector, or formula).
#' @param error Specification for error components, or NULL for default.
#' @param aggregation Character vector of facets to aggregate over.
#' @param residual_is Character string specifying which facets make up the residual.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @return The G coefficient as a numeric value.
#'
#' @keywords internal
compute_g_coefficient <- function(variance_components, n, object, error = NULL,
                                  aggregation = NULL, residual_is = NULL, n_provided = FALSE) {
  d_vc <- calculate_dstudy_variance(variance_components, n, object, aggregation, n_provided, residual_is)

  coefs <- calculate_coefficients(d_vc, n, object, error, aggregation, residual_is)

  coefs$g
}

#' Compute Dependability Coefficient
#'
#' Calculates the dependability coefficient (Phi coefficient), which represents
#' the proportion of observed variance attributable to universe score variance
#' for absolute decisions.
#'
#' @param variance_components A tibble of variance components.
#' @param n A named list specifying the number of levels for each facet.
#' @param object Specification for object of measurement (character, vector, or formula).
#' @param error Specification for error components, or NULL for default.
#' @param aggregation Character vector of facets to aggregate over.
#' @param residual_is Character string specifying which facets make up the residual.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @return The Phi coefficient as a numeric value.
#'
#' @keywords internal
compute_phi_coefficient <- function(variance_components, n, object, error = NULL,
                                    aggregation = NULL, residual_is = NULL, n_provided = FALSE) {
  d_vc <- calculate_dstudy_variance(variance_components, n, object, aggregation, n_provided, residual_is)

  coefs <- calculate_coefficients(d_vc, n, object, error, aggregation, residual_is)

  coefs$phi
}

#' Compute Standard Error of Measurement
#'
#' Calculates the standard error of measurement (SEM) for a D-study design.
#'
#' @param error_variance The error variance (relative or absolute).
#' @return The standard error of measurement.
#'
#' @keywords internal
compute_sem <- function(error_variance) {
  sqrt(error_variance)
}

#' Compute Value Added Ratio (VAR)
#'
#' Calculates the Value Added Ratio (VAR) for a multivariate composite design.
#' VAR expresses how much the composite reliability improves over that of a
#' Compute Observed Score Variances and Covariances for VAR Calculation
#'
#' Computes observed score variances and covariances from the original data,
#' which are required for the Haberman (2008) formula in VAR calculation.
#' These are data quantities (not model estimates) and are constant across
#' posterior draws.
#'
#' @param data A data frame containing the subscale scores.
#' @param dimensions Character vector of dimension/subscale names.
#' @param weights Named numeric vector of weights for each dimension.
#'
#' @return A list with:
#' \describe{
#' \item{obs_var}{Named vector of observed score variances per dimension}
#' \item{obs_cov}{Matrix of observed score covariances between dimensions}
#' \item{obs_var_C}{Observed variance of the weighted composite}
#' }
#'
#' @keywords internal
compute_observed_covariances <- function(data, dimensions, weights, dimension_var = NULL) {
  # Handle long-format data: reshape to wide format
  if (!is.null(dimension_var) && dimension_var %in% names(data)) {
    # Find the score column (first numeric column that isn't dimension_var)
    score_col <- NULL
    for (col in names(data)) {
      if (col != dimension_var && is.numeric(data[[col]])) {
        score_col <- col
        break
      }
    }

    if (!is.null(score_col)) {
      # Reshape from long to wide format
      data <- tidyr::pivot_wider(
        data,
        id_cols = setdiff(names(data), c(dimension_var, score_col)),
        names_from = all_of(dimension_var),
        values_from = all_of(score_col)
      )
    }
  }

  n_dims <- length(dimensions)
  obs_var <- vapply(dimensions, function(d) {
    var(data[[d]], na.rm = TRUE)
  }, numeric(1))
  names(obs_var) <- dimensions

  obs_cov <- matrix(0, nrow = n_dims, ncol = n_dims)
  rownames(obs_cov) <- colnames(obs_cov) <- dimensions
  for (i in seq_len(n_dims)) {
    for (j in seq_len(n_dims)) {
      if (i == j) {
        obs_cov[i, j] <- obs_var[i]
      } else {
        obs_cov[i, j] <- cov(data[[dimensions[i]]], data[[dimensions[j]]], use = "complete.obs")
      }
    }
  }

  w <- weights[dimensions]
  obs_var_C <- as.numeric(t(w) %*% obs_cov %*% w)

  list(
    obs_var = obs_var,
    obs_cov = obs_cov,
    obs_var_C = obs_var_C
  )
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
#' Where PRMSE(C -> S_i) is:
#' \deqn{PRMSE(C \rightarrow S_i) = \frac{(\sigma^2_{S_i} \cdot Rel(S_i) + \sum_{j \neq i} \sigma_{S_i, S_j})^2}{\sigma^2_{S_i} \cdot Rel(C) \cdot \sigma^2_C}}
#'
#' @param uni_cov_draws 3D array [n_draws, n_dims, n_dims] of universe score covariances.
#' @param total_rel_cov_draws 3D array [n_draws, n_dims, n_dims] of total observed covariances for relative error.
#' @param total_abs_cov_draws 3D array [n_draws, n_dims, n_dims] of total observed covariances for absolute error.
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

#' Compute PRMSE Confidence Intervals via Delta Method
#'
#' Computes confidence intervals for PRMSE metrics using the delta method,
#' propagating uncertainty from variance component estimates.
#'
#' @param dstudy_obj A dstudy object
#' @param gstudy_obj The associated gstudy object
#' @param metrics Character vector: "prmse", "var", or both
#' @param probs Numeric vector of length 2 for CI bounds
#' @param n Named list of D-study sample sizes
#' @param weights Named numeric vector of dimension weights
#' @param object Object of measurement specification
#' @param universe Universe components specification
#' @param error Error components specification
#' @param aggregation Aggregation specification
#' @param residual_is Residual specification
#'
#' @return A list with CI bounds for each requested metric
#'
#' @details
#' The delta method approximates the variance of a function using:
#' Var(f(X)) = sum((df/dx_i)^2 * Var(x_i))
#'
#' For PRMSE metrics, this involves computing partial derivatives with
#' respect to each variance component and propagating their SEs.
#'
#' @keywords internal
compute_prmse_delta_ci <- function(dstudy_obj, gstudy_obj, metrics, probs,
                                   n, weights, object, universe, error,
                                   aggregation, residual_is) {
  # Extract variance components with CIs
  vc <- gstudy_obj$variance_components

  # Check if SEs are available
  if (!"se" %in% names(vc) && !"lower" %in% names(vc)) {
    warning("Variance component SEs not available. Cannot compute delta method CIs.", call. = FALSE)
    # Return dimension-specific NA values
    dims <- unique(dstudy_obj$coefficients$dim)
    dims <- dims[dims != "Composite"]
    ci_result <- list()
    for (d in dims) {
      ci_result[[paste0("prmse_s_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_abs_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_abs_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_abs_UL_", d)]] <- NA_real_
      ci_result[[paste0("var_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("var_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("var_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("var_abs_UL_", d)]] <- NA_real_
    }
    return(ci_result)
  }

  # Get point estimates from dstudy coefficients
  coefs <- dstudy_obj$coefficients
  dims <- unique(coefs$dim)
  dims <- dims[dims != "Composite"]

  z_crit <- qnorm(1 - (1 - (probs[2] - probs[1])) / 2)

  # Initialize result list
  ci_result <- list()

  # Get variance component SEs per dimension
  vc_se <- vc
  if ("se" %in% names(vc_se)) {
    # Use SE column
  } else if ("lower" %in% names(vc_se) && "upper" %in% names(vc_se)) {
    # Approximate SE from CI width
    vc_se$se <- (vc_se$upper - vc_se$lower) / (2 * z_crit)
  }

  # For each dimension, compute CI for G and Phi coefficients
  for (dim in dims) {
    # Get variance components for this dimension
    vc_dim <- vc_se[vc_se$dim == dim, ]

    if (nrow(vc_dim) == 0) next

    # Compute universe and error variance
    universe_spec <- if (!is.null(universe)) parse_specification(universe) else parse_specification(object)

    uni_var <- sum(vc_dim$var[vc_dim$component %in% universe_spec], na.rm = TRUE)
    err_rel_var <- sum(vc_dim$var[vc_dim$component == "Residual" |
      grepl("^Residual", vc_dim$component) |
      sapply(vc_dim$component, function(c) {
        if (c %in% universe_spec) return(FALSE)
        facets <- parse_component_facets(c)
        any(universe_spec %in% facets)
      })], na.rm = TRUE)

    # SEs for universe and error variance (using sum of variances)
    uni_se <- sqrt(sum(vc_dim$se[vc_dim$component %in% universe_spec]^2, na.rm = TRUE))
    err_rel_se <- sqrt(sum(vc_dim$se[vc_dim$component == "Residual" |
      grepl("^Residual", vc_dim$component) |
      sapply(vc_dim$component, function(c) {
        if (c %in% universe_spec) return(FALSE)
        facets <- parse_component_facets(c)
        any(universe_spec %in% facets)
      })]^2, na.rm = TRUE))

    # G coefficient and its SE
    g <- if ((uni_var + err_rel_var) > 0) uni_var / (uni_var + err_rel_var) else NA_real_

    # Delta method for ratio: Var(a/b) = (a/b)^2 * (Var(a)/a^2 + Var(b)/b^2 - 2Cov(a,b)/(a*b))
    # For simplicity, assume independent estimates (Cov = 0)
    if (!is.na(g) && g > 0 && !is.na(uni_se) && !is.na(err_rel_se)) {
      var_num <- uni_se^2
      var_den <- uni_se^2 + err_rel_se^2
      se_g <- abs(g) * sqrt(var_num / uni_var^2 + var_den / (uni_var + err_rel_var)^2)

      ci_result[[paste0("prmse_s_rel_LL_", dim)]] <- max(0, g - z_crit * se_g)
      ci_result[[paste0("prmse_s_rel_UL_", dim)]] <- min(1, g + z_crit * se_g)
    } else {
      ci_result[[paste0("prmse_s_rel_LL_", dim)]] <- NA_real_
      ci_result[[paste0("prmse_s_rel_UL_", dim)]] <- NA_real_
    }

    # Phi coefficient (absolute error)
    all_err_var <- sum(vc_dim$var[!vc_dim$component %in% universe_spec], na.rm = TRUE)
    all_err_se <- sqrt(sum(vc_dim$se[!vc_dim$component %in% universe_spec]^2, na.rm = TRUE))

    phi <- if ((uni_var + all_err_var) > 0) uni_var / (uni_var + all_err_var) else NA_real_

    if (!is.na(phi) && phi > 0 && !is.na(uni_se) && !is.na(all_err_se)) {
      var_num <- uni_se^2
      var_den <- uni_se^2 + all_err_se^2
      se_phi <- abs(phi) * sqrt(var_num / uni_var^2 + var_den / (uni_var + all_err_var)^2)

      ci_result[[paste0("prmse_s_abs_LL_", dim)]] <- max(0, phi - z_crit * se_phi)
      ci_result[[paste0("prmse_s_abs_UL_", dim)]] <- min(1, phi + z_crit * se_phi)
    } else {
      ci_result[[paste0("prmse_s_abs_LL_", dim)]] <- NA_real_
      ci_result[[paste0("prmse_s_abs_UL_", dim)]] <- NA_real_
    }
  }

  # For PRMSE_C and VAR, the formulas are more complex
  # We'll use a simplified approximation based on G coefficient SE
  # Full implementation would require numerical derivatives of Equation 38

  # Get existing VAR values
  if (!is.null(dstudy_obj$var)) {
    var_data <- dstudy_obj$var
    for (dim in dims) {
      if (!is.null(var_data[[dim]])) {
        # Approximate CI using delta method through VAR
        # VAR = G / PRMSE_C, so SE(VAR) is complex
        # For now, use a conservative approximation

        # Get prmse_c values
        prmse_c_rel <- var_data[[dim]]$prmse_c_rel
        prmse_c_abs <- var_data[[dim]]$prmse_c_abs

        # Approximate SE using coefficient of variation
        g_val <- coefs$g[coefs$dim == dim][1]
        g_se <- if (!is.na(ci_result[[paste0("prmse_s_rel_LL_", dim)]])) {
          (ci_result[[paste0("prmse_s_rel_UL_", dim)]] - ci_result[[paste0("prmse_s_rel_LL_", dim)]]) / (2 * z_crit)
        } else NA_real_

        if (!is.na(prmse_c_rel) && prmse_c_rel > 0 && !is.na(g_se)) {
          # Approximate SE for prmse_c (conservative)
          se_prmse_c <- g_se * prmse_c_rel / g_val
          ci_result[[paste0("prmse_c_rel_LL_", dim)]] <- max(0, prmse_c_rel - z_crit * se_prmse_c)
          ci_result[[paste0("prmse_c_rel_UL_", dim)]] <- min(1, prmse_c_rel + z_crit * se_prmse_c)
        } else {
          ci_result[[paste0("prmse_c_rel_LL_", dim)]] <- NA_real_
          ci_result[[paste0("prmse_c_rel_UL_", dim)]] <- NA_real_
        }

        if (!is.na(prmse_c_abs) && prmse_c_abs > 0 && !is.na(g_se)) {
          se_prmse_c <- g_se * prmse_c_abs / g_val
          ci_result[[paste0("prmse_c_abs_LL_", dim)]] <- max(0, prmse_c_abs - z_crit * se_prmse_c)
          ci_result[[paste0("prmse_c_abs_UL_", dim)]] <- min(1, prmse_c_abs + z_crit * se_prmse_c)
        } else {
          ci_result[[paste0("prmse_c_abs_LL_", dim)]] <- NA_real_
          ci_result[[paste0("prmse_c_abs_UL_", dim)]] <- NA_real_
        }

        # VAR CI (propagate through ratio)
        var_rel <- var_data[[dim]]$var_rel
        var_abs <- var_data[[dim]]$var_abs

        if (!is.na(var_rel) && !is.na(g_se) && !is.na(prmse_c_rel) && prmse_c_rel > 0) {
          # Delta method for ratio: SE(VAR) = VAR * sqrt((SE_g/g)^2 + (SE_c/c)^2)
          se_var <- var_rel * sqrt((g_se / g_val)^2 + (g_se / prmse_c_rel)^2)
          ci_result[[paste0("var_rel_LL_", dim)]] <- max(0, var_rel - z_crit * se_var)
          ci_result[[paste0("var_rel_UL_", dim)]] <- var_rel + z_crit * se_var
        } else {
          ci_result[[paste0("var_rel_LL_", dim)]] <- NA_real_
          ci_result[[paste0("var_rel_UL_", dim)]] <- NA_real_
        }

        if (!is.na(var_abs) && !is.na(g_se)) {
          se_var <- var_abs * sqrt((g_se / g_val)^2)
          ci_result[[paste0("var_abs_LL_", dim)]] <- max(0, var_abs - z_crit * se_var)
          ci_result[[paste0("var_abs_UL_", dim)]] <- var_abs + z_crit * se_var
        } else {
          ci_result[[paste0("var_abs_LL_", dim)]] <- NA_real_
          ci_result[[paste0("var_abs_UL_", dim)]] <- NA_real_
        }

        # PRMSE_P (use similar approximation)
        prmse_p_rel <- var_data[[dim]]$prmse_p_rel
        prmse_p_abs <- var_data[[dim]]$prmse_p_abs

        if (!is.na(prmse_p_rel) && !is.na(g_se)) {
          se_p <- g_se * prmse_p_rel / g_val
          ci_result[[paste0("prmse_p_rel_LL_", dim)]] <- max(0, prmse_p_rel - z_crit * se_p)
          ci_result[[paste0("prmse_p_rel_UL_", dim)]] <- min(1, prmse_p_rel + z_crit * se_p)
        } else {
          ci_result[[paste0("prmse_p_rel_LL_", dim)]] <- NA_real_
          ci_result[[paste0("prmse_p_rel_UL_", dim)]] <- NA_real_
        }

        if (!is.na(prmse_p_abs) && !is.na(g_se)) {
          se_p <- g_se * prmse_p_abs / g_val
          ci_result[[paste0("prmse_p_abs_LL_", dim)]] <- max(0, prmse_p_abs - z_crit * se_p)
          ci_result[[paste0("prmse_p_abs_UL_", dim)]] <- min(1, prmse_p_abs + z_crit * se_p)
        } else {
          ci_result[[paste0("prmse_p_abs_LL_", dim)]] <- NA_real_
          ci_result[[paste0("prmse_p_abs_UL_", dim)]] <- NA_real_
        }
      }
    }
  }

  ci_result
}

#' Compute PRMSE Confidence Intervals via Bootstrap
#'
#' Computes confidence intervals for PRMSE metrics using parametric bootstrap,
#' resampling variance components from their estimated distributions.
#'
#' @param dstudy_obj A dstudy object
#' @param gstudy_obj The associated gstudy object
#' @param metrics Character vector: "prmse", "var", or both
#' @param probs Numeric vector of length 2 for CI bounds
#' @param n_bootstrap Number of bootstrap samples
#' @param n Named list of D-study sample sizes
#' @param weights Named numeric vector of dimension weights
#' @param object Object of measurement specification
#' @param universe Universe components specification
#' @param error Error components specification
#' @param aggregation Aggregation specification
#' @param residual_is Residual specification
#'
#' @return A list with CI bounds for each requested metric
#'
#' @details
#' For each bootstrap iteration:
#' 1. Resample variance components from N(VC, SE^2)
#' 2. Recompute PRMSE metrics
#' 3. Use percentiles for CI bounds
#'
#' @keywords internal
compute_prmse_bootstrap_ci <- function(dstudy_obj, gstudy_obj, metrics, probs,
                                       n_bootstrap, n, weights, object, universe,
                                       error, aggregation, residual_is) {
  # Extract variance components with SEs
  vc <- gstudy_obj$variance_components

  # Check if SEs are available
  if (!"se" %in% names(vc) && !"lower" %in% names(vc)) {
    warning("Variance component SEs not available. Cannot compute bootstrap CIs.", call. = FALSE)
    # Return dimension-specific NA values
    dims <- unique(dstudy_obj$coefficients$dim)
    dims <- dims[dims != "Composite"]
    ci_result <- list()
    for (d in dims) {
      ci_result[[paste0("prmse_s_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_abs_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_abs_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_p_abs_UL_", d)]] <- NA_real_
      ci_result[[paste0("var_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("var_rel_UL_", d)]] <- NA_real_
      ci_result[[paste0("var_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("var_abs_UL_", d)]] <- NA_real_
    }
    return(ci_result)
  }

  # Get SE estimates
  z_crit <- qnorm(1 - (1 - (probs[2] - probs[1])) / 2)
  if (!"se" %in% names(vc) && "lower" %in% names(vc)) {
    vc$se <- (vc$upper - vc$lower) / (2 * z_crit)
  }

  # Get dimensions
  coefs <- dstudy_obj$coefficients
  dims <- unique(coefs$dim)
  dims <- dims[dims != "Composite"]

  # Initialize bootstrap storage
  n_dims <- length(dims)
  boot_prmse_s_rel <- matrix(NA_real_, n_bootstrap, n_dims)
  boot_prmse_s_abs <- matrix(NA_real_, n_bootstrap, n_dims)
  boot_prmse_c_rel <- matrix(NA_real_, n_bootstrap, n_dims)
  boot_prmse_c_abs <- matrix(NA_real_, n_bootstrap, n_dims)
  boot_var_rel <- matrix(NA_real_, n_bootstrap, n_dims)
  boot_var_abs <- matrix(NA_real_, n_bootstrap, n_dims)
  colnames(boot_prmse_s_rel) <- dims
  colnames(boot_prmse_s_abs) <- dims
  colnames(boot_prmse_c_rel) <- dims
  colnames(boot_prmse_c_abs) <- dims
  colnames(boot_var_rel) <- dims
  colnames(boot_var_abs) <- dims

  # Get specification
  universe_spec <- if (!is.null(universe)) parse_specification(universe) else parse_specification(object)

  # Bootstrap loop
  for (b in seq_len(n_bootstrap)) {
    # Resample variance components
    vc_boot <- vc
    for (i in seq_len(nrow(vc_boot))) {
      if (!is.na(vc_boot$se[i]) && vc_boot$se[i] > 0) {
        vc_boot$var[i] <- max(0, rnorm(1, mean = vc_boot$var[i], sd = vc_boot$se[i]))
      }
    }

    # Recompute coefficients for each dimension
    for (j in seq_along(dims)) {
      d <- dims[j]
      vc_d <- vc_boot[vc_boot$dim == d, ]

      if (nrow(vc_d) == 0) next

      # Universe variance
      uni_var <- sum(vc_d$var[vc_d$component %in% universe_spec], na.rm = TRUE)

      # Relative error variance
      rel_comps <- vc_d$component[sapply(vc_d$component, function(c) {
        if (c %in% universe_spec) return(FALSE)
        if (c == "Residual") return(TRUE)
        facets <- parse_component_facets(c)
        any(universe_spec %in% facets)
      })]
      err_rel_var <- sum(vc_d$var[vc_d$component %in% rel_comps], na.rm = TRUE)

      # Absolute error variance
      err_abs_var <- sum(vc_d$var[!vc_d$component %in% universe_spec], na.rm = TRUE)

      # G and Phi
      g_b <- if ((uni_var + err_rel_var) > 0) uni_var / (uni_var + err_rel_var) else NA_real_
      phi_b <- if ((uni_var + err_abs_var) > 0) uni_var / (uni_var + err_abs_var) else NA_real_

      boot_prmse_s_rel[b, j] <- g_b
      boot_prmse_s_abs[b, j] <- phi_b

      # Get existing PRMSE_C and VAR (approximate - would need full recomputation)
      if (!is.null(dstudy_obj$var) && !is.null(dstudy_obj$var[[d]])) {
        # Scale by bootstrap G relative to original G
        g_orig <- coefs$g[coefs$dim == d][1]
        if (!is.na(g_b) && !is.na(g_orig) && g_orig > 0) {
          ratio <- g_b / g_orig
          boot_prmse_c_rel[b, j] <- dstudy_obj$var[[d]]$prmse_c_rel * ratio
          boot_prmse_c_abs[b, j] <- dstudy_obj$var[[d]]$prmse_c_abs * ratio
          # VAR = G / PRMSE_C
          if (!is.na(dstudy_obj$var[[d]]$prmse_c_rel) && dstudy_obj$var[[d]]$prmse_c_rel > 0) {
            boot_var_rel[b, j] <- g_b / (dstudy_obj$var[[d]]$prmse_c_rel * ratio)
          }
          if (!is.na(dstudy_obj$var[[d]]$prmse_c_abs) && dstudy_obj$var[[d]]$prmse_c_abs > 0) {
            boot_var_abs[b, j] <- g_b / (dstudy_obj$var[[d]]$prmse_c_abs * ratio)
          }
        }
      }
    }
  }

  # Compute CI from percentiles
  ci_result <- list()

  for (j in seq_along(dims)) {
    d <- dims[j]

    # prmse_s_rel CI
    if (sum(!is.na(boot_prmse_s_rel[, j])) > 10) {
      ci_result[[paste0("prmse_s_rel_LL_", d)]] <- quantile(boot_prmse_s_rel[, j], probs = probs[1], na.rm = TRUE)
      ci_result[[paste0("prmse_s_rel_UL_", d)]] <- quantile(boot_prmse_s_rel[, j], probs = probs[2], na.rm = TRUE)
    } else {
      ci_result[[paste0("prmse_s_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_rel_UL_", d)]] <- NA_real_
    }

    # prmse_s_abs CI
    if (sum(!is.na(boot_prmse_s_abs[, j])) > 10) {
      ci_result[[paste0("prmse_s_abs_LL_", d)]] <- quantile(boot_prmse_s_abs[, j], probs = probs[1], na.rm = TRUE)
      ci_result[[paste0("prmse_s_abs_UL_", d)]] <- quantile(boot_prmse_s_abs[, j], probs = probs[2], na.rm = TRUE)
    } else {
      ci_result[[paste0("prmse_s_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_s_abs_UL_", d)]] <- NA_real_
    }

    # prmse_c_rel CI
    if (sum(!is.na(boot_prmse_c_rel[, j])) > 10) {
      ci_result[[paste0("prmse_c_rel_LL_", d)]] <- quantile(boot_prmse_c_rel[, j], probs = probs[1], na.rm = TRUE)
      ci_result[[paste0("prmse_c_rel_UL_", d)]] <- quantile(boot_prmse_c_rel[, j], probs = probs[2], na.rm = TRUE)
    } else {
      ci_result[[paste0("prmse_c_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_rel_UL_", d)]] <- NA_real_
    }

    # prmse_c_abs CI
    if (sum(!is.na(boot_prmse_c_abs[, j])) > 10) {
      ci_result[[paste0("prmse_c_abs_LL_", d)]] <- quantile(boot_prmse_c_abs[, j], probs = probs[1], na.rm = TRUE)
      ci_result[[paste0("prmse_c_abs_UL_", d)]] <- quantile(boot_prmse_c_abs[, j], probs = probs[2], na.rm = TRUE)
    } else {
      ci_result[[paste0("prmse_c_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("prmse_c_abs_UL_", d)]] <- NA_real_
    }

    # var_rel CI
    if (sum(!is.na(boot_var_rel[, j])) > 10) {
      ci_result[[paste0("var_rel_LL_", d)]] <- quantile(boot_var_rel[, j], probs = probs[1], na.rm = TRUE)
      ci_result[[paste0("var_rel_UL_", d)]] <- quantile(boot_var_rel[, j], probs = probs[2], na.rm = TRUE)
    } else {
      ci_result[[paste0("var_rel_LL_", d)]] <- NA_real_
      ci_result[[paste0("var_rel_UL_", d)]] <- NA_real_
    }

    # var_abs CI
    if (sum(!is.na(boot_var_abs[, j])) > 10) {
      ci_result[[paste0("var_abs_LL_", d)]] <- quantile(boot_var_abs[, j], probs = probs[1], na.rm = TRUE)
      ci_result[[paste0("var_abs_UL_", d)]] <- quantile(boot_var_abs[, j], probs = probs[2], na.rm = TRUE)
    } else {
      ci_result[[paste0("var_abs_LL_", d)]] <- NA_real_
      ci_result[[paste0("var_abs_UL_", d)]] <- NA_real_
    }

    # prmse_p CIs - not computed in bootstrap (requires full matrix recomputation)
    # Set to NA for now
    ci_result[[paste0("prmse_p_rel_LL_", d)]] <- NA_real_
    ci_result[[paste0("prmse_p_rel_UL_", d)]] <- NA_real_
    ci_result[[paste0("prmse_p_abs_LL_", d)]] <- NA_real_
    ci_result[[paste0("prmse_p_abs_UL_", d)]] <- NA_real_
  }

  ci_result
}

#' Identify Error Components for D-Study (Internal)
#'
#' Determines which variance components belong to relative and absolute error.
#' This consolidates logic used by both point-estimate and draw-based computations.
#'
#' @param components Character vector of component names
#' @param universe_spec Character vector of universe components
#' @param error_spec Character vector of explicit error components (or NULL)
#' @param agg_facets Character vector of aggregation facets (or NULL)
#'
#' @return List with `is_relative_error` and `is_absolute_error` logical vectors
#'
#' @keywords internal
identify_error_components_for_draws <- function(components, universe_spec, error_spec = NULL,
                                                agg_facets = NULL) {
  n_comp <- length(components)
  is_universe <- components %in% universe_spec

  is_agg_main <- logical(n_comp)
  if (!is.null(agg_facets)) {
    for (i in seq_len(n_comp)) {
      comp <- components[i]
      if (comp == "Residual") {
        is_agg_main[i] <- FALSE
      } else {
        comp_facets <- parse_component_facets(comp)
        is_agg_main[i] <- length(comp_facets) == 1 && all(comp_facets %in% agg_facets)
      }
    }
  }

  if (!is.null(error_spec)) {
    is_error_component <- components %in% error_spec
    return(list(
      is_relative_error = is_error_component,
      is_absolute_error = is_error_component,
      error_spec_used = TRUE
    ))
  }

  is_relative_error <- logical(n_comp)
  is_absolute_error <- logical(n_comp)

  for (i in seq_len(n_comp)) {
    comp <- components[i]

    if (is_universe[i] || is_agg_main[i]) {
      next
    }

    is_absolute_error[i] <- TRUE

    if (comp == "Residual" || grepl("^Residual", comp)) {
      is_relative_error[i] <- TRUE
    } else {
      comp_facets <- parse_component_facets(comp)
      if (any(universe_spec %in% comp_facets)) {
        is_relative_error[i] <- TRUE
      }
    }
  }

  list(
    is_relative_error = is_relative_error,
    is_absolute_error = is_absolute_error,
    error_spec_used = FALSE
  )
}

#' Calculate Coefficients Using Posterior Draws
#'
#' Computes generalizability and dependability coefficients using full posterior
#' distributions from a brms model. This provides proper uncertainty quantification
#' for D-study coefficients.
#'
#' @param gstudy_obj A gstudy object fitted with backend = "brms".
#' @param n Named list of sample sizes for each facet.
#' @param object Specification for object of measurement.
#' @param error Specification for error components.
#'   Note: If the same facet is specified in both `object` and `error`, an error is raised.
#' @param aggregation Character vector of facets to aggregate over.
#'   Note: If the same facet is specified in both `object` and `aggregation`, an error is raised.
#' @param residual_is Character string specifying residual composition.
#' @param is_sweep Logical indicating if this is a sweep over multiple sample sizes.
#' @param n_grid If is_sweep, a data frame with all combinations of sample sizes.
#' @param cut_score Optional cut score for phi-cut calculation.
#' @param mu_y Optional grand mean for phi-cut calculation.
#'
#' @return A list with:
#'   \item{coefficients}{Data frame with coefficient means (and sweep combinations if applicable)}
#'   \item{posterior}{List of posterior distribution vectors (or list of lists for sweep)}
#'
#' @keywords internal
calculate_coefficients_posterior <- function(gstudy_obj, n, object = NULL, universe = NULL, error = NULL,
                                             aggregation = NULL, residual_is = NULL,
                                             is_sweep = FALSE, n_grid = NULL, n_provided = FALSE, use_scaled = TRUE,
                                             cut_score = NULL, mu_y = NULL, ci = NULL, probs = c(0.025, 0.975),
                                             weights = NULL, n_per_dim = NULL) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required for posterior estimation.", call. = FALSE)
  }

  # Extract posterior draws
  draws <- brms::as_draws_matrix(gstudy_obj$model)

  # Get variance component names and their draws
  vc_draws <- extract_variance_draws(gstudy_obj, draws)

  # Extract covariance draws for multivariate composite calculation
  cov_draws <- NULL
  if (gstudy_obj$is_multivariate && length(gstudy_obj$dimensions) > 1) {
    cov_draws <- extract_covariance_draws(
      gstudy_obj$model,
      gstudy_obj$dimensions
    )
  }

  # Get object of measurement (always from G-study)
  if (is.null(object)) {
    object <- gstudy_obj$object
  }
  object_spec <- parse_specification(object)

  # Process universe specification
  # Default: universe = object only
  if (is.null(universe) || length(universe) == 0) {
    universe_spec <- object_spec
  } else {
    universe_spec <- parse_specification(universe)
  }

  # Parse error specification
  if (!is.null(error)) {
    error_spec <- parse_specification(error)
  } else {
    error_spec <- NULL
  }

  # Parse aggregation specification
  if (!is.null(aggregation)) {
    agg_facets <- parse_specification(aggregation)
  } else {
    agg_facets <- NULL
  }

  # Validate universe vs error
  if (!is.null(error_spec) && length(error_spec) > 0) {
    overlap <- intersect(universe_spec, error_spec)
    if (length(overlap) > 0) {
      stop(
        "Component(s) '", paste(overlap, collapse = "', '"), "' specified as both ",
        "in the universe and in error. ",
        "The same component cannot be in both the universe and error.",
        call. = FALSE
      )
    }
  }

  is_mv <- gstudy_obj$is_multivariate
  mv_dims <- if (is_mv) gstudy_obj$dimensions else character(0)

  if (is_sweep) {
    # Calculate for each sweep combination
    results <- vector("list", nrow(n_grid))
    posterior_list <- vector("list", nrow(n_grid))

    for (i in seq_len(nrow(n_grid))) {
      n_current <- as.list(n_grid[i, , drop = FALSE])
      n_current <- lapply(n_current, as.numeric)

      if (is_mv) {
        # Multivariate: iterate over dimensions
        dim_results <- vector("list", length(mv_dims))
        names(dim_results) <- mv_dims
        dim_posteriors <- vector("list", length(mv_dims))
        names(dim_posteriors) <- mv_dims

        for (d in mv_dims) {
          combo_result <- calculate_single_posterior(
            vc_draws = vc_draws[[d]],
            n = n_current,
            object_spec = object_spec,
            universe_spec = universe_spec,
            error_spec = error_spec,
            agg_facets = agg_facets,
            residual_is = residual_is,
            gstudy_obj = gstudy_obj,
            n_provided = n_provided,
            use_scaled = use_scaled,
            cut_score = cut_score,
            mu_y = if (is.list(mu_y)) mu_y[[d]] else mu_y,
            ci = ci,
            probs = probs
          )
          dim_results[[d]] <- combo_result$summary
          dim_posteriors[[d]] <- combo_result$distributions
        }

        # Combine results with dimension column
        results[[i]] <- cbind(
          data.frame(n_current, stringsAsFactors = FALSE, check.names = FALSE),
          do.call(rbind, lapply(names(dim_results), function(d) {
            cbind(dim = d, dim_results[[d]])
          }))
        )
        posterior_list[[i]] <- dim_posteriors
      } else {
        # Univariate: single dimension
        combo_result <- calculate_single_posterior(
          vc_draws = vc_draws,
          n = n_current,
          object_spec = object_spec,
          universe_spec = universe_spec,
          error_spec = error_spec,
          agg_facets = agg_facets,
          residual_is = residual_is,
          gstudy_obj = gstudy_obj,
          n_provided = n_provided,
          use_scaled = use_scaled,
          cut_score = cut_score,
          mu_y = mu_y,
          ci = ci,
          probs = probs
        )

        results[[i]] <- cbind(
          data.frame(n_current, stringsAsFactors = FALSE, check.names = FALSE),
          combo_result$summary
        )
        posterior_list[[i]] <- combo_result$distributions
      }
    }

    coefficients <- do.call(rbind, results)
    coefficients <- tibble::as_tibble(coefficients)

    return(list(
      coefficients = coefficients,
      posterior = posterior_list
    ))
  } else {
    # Single sample size
    if (is_mv) {
      dim_results <- vector("list", length(mv_dims))
      names(dim_results) <- mv_dims
      dim_posteriors <- vector("list", length(mv_dims))
      names(dim_posteriors) <- mv_dims

      for (d in mv_dims) {
        result <- calculate_single_posterior(
          vc_draws = vc_draws[[d]],
          n = n,
          object_spec = object_spec,
          universe_spec = universe_spec,
          error_spec = error_spec,
          agg_facets = agg_facets,
          residual_is = residual_is,
          gstudy_obj = gstudy_obj,
          n_provided = n_provided,
          use_scaled = use_scaled,
          cut_score = cut_score,
          mu_y = if (is.list(mu_y)) mu_y[[d]] else mu_y,
          ci = ci,
          probs = probs
        )
        dim_results[[d]] <- cbind(dim = d, result$summary)
        dim_posteriors[[d]] <- result$distributions
      }

      coefficients <- tibble::as_tibble(do.call(rbind, lapply(mv_dims, function(d) dim_results[[d]])))

      composite_result <- NULL
      composite_posterior <- NULL

      if (!is.null(cov_draws) && !is.null(weights)) {
        vc_result <- calculate_dstudy_variance_composite(
          vc_draws = vc_draws,
          cov_draws = cov_draws,
          weights = weights,
          n = n,
          object = object,
          n_provided = n_provided
        )

        components <- names(vc_draws[[1]])
        scale_factors <- list()
        for (comp in components) {
          scale_factors[[comp]] <- compute_component_scale_factor(comp, n, object_spec, n_provided)
        }

        coef_result <- calculate_composite_coefficients_from_draws(
          vc_draws = vc_draws,
          cov_draws = cov_draws,
          weights = weights,
          scale_factors = scale_factors,
          object_spec = object_spec,
          universe_spec = universe_spec,
          error_spec = error_spec,
          agg_facets = agg_facets,
          residual_is = residual_is,
          cut_score = cut_score,
          mu_y = mu_y,
          ci = ci,
          probs = probs,
          gstudy_data = gstudy_obj$data,
          dimension_var = gstudy_obj$dimension_var
        )

        coefficients <- dplyr::bind_rows(coefficients, coef_result$summary)

        var_results <- coef_result$var_results

        composite_posterior <- coef_result$distributions

        return(list(
          coefficients = coefficients,
          posterior = dim_posteriors,
          composite_vc = vc_result$vc_table,
          composite_posterior = composite_posterior,
          var_results = var_results
        ))
      }

      return(list(
        coefficients = coefficients,
        posterior = dim_posteriors,
        composite_vc = NULL,
        composite_posterior = NULL,
        var_results = NULL
      ))
    } else {
      result <- calculate_single_posterior(
        vc_draws = vc_draws,
        n = n,
        object_spec = object_spec,
        universe_spec = universe_spec,
        error_spec = error_spec,
        agg_facets = agg_facets,
        residual_is = residual_is,
        gstudy_obj = gstudy_obj,
        n_provided = n_provided,
        use_scaled = use_scaled,
        cut_score = cut_score,
        mu_y = mu_y,
        ci = ci,
        probs = probs
      )

      return(list(
        coefficients = tibble::as_tibble(result$summary),
        posterior = result$distributions
      ))
    }
  }
}

#' Extract Variance Component Draws from brms Model
#'
#' Extracts posterior draws of variance components (SD^2) from a brms model.
#' For multivariate models, returns a nested list structure by dimension.
#'
#' @param gstudy_obj A gstudy object with brms backend.
#' @param draws A posterior draws_matrix object (e.g., from brms::as_draws_matrix()).
#'
#' @return A named list. For univariate: list of variance draws per component.
#'   For multivariate: nested list with dimensions as names, each containing
#'   variance draws for that dimension.
#'
#' @keywords internal
extract_variance_draws <- function(gstudy_obj, draws) {
  model <- gstudy_obj$model
  vc <- gstudy_obj$variance_components
  is_mv <- gstudy_obj$is_multivariate
  is_long_format <- isTRUE(gstudy_obj$long_format_multivariate)
  dimension_var <- gstudy_obj$dimension_var

  if (is_mv) {
    dims <- unique(vc$dim)
    vc_draws <- vector("list", length(dims))
    names(vc_draws) <- dims

    for (d in dims) {
      vc_dim <- vc[vc$dim == d, ]
      vc_draws[[d]] <- list()

      for (comp in vc_dim$component) {
        sd_draws <- extract_sd_draws_multivariate(
          draws = draws,
          comp = comp,
          dim = d,
          vc_dim = vc_dim,
          is_long_format = is_long_format,
          dimension_var = dimension_var
        )
        vc_draws[[d]][[comp]] <- sd_draws^2
      }
    }
  } else {
    vc_draws <- list()

    for (comp in vc$component) {
      sd_draws <- extract_sd_draws_univariate(
        draws = draws,
        comp = comp,
        vc = vc
      )
      vc_draws[[comp]] <- sd_draws^2
    }
  }

  vc_draws
}

extract_sd_draws_multivariate <- function(draws, comp, dim, vc_dim, is_long_format = FALSE, dimension_var = NULL) {
  param_candidates <- character()

  if (is_long_format && !is.null(dimension_var)) {
    dim_in_params <- paste0(dimension_var, dim)
  } else {
    dim_in_params <- dim
  }

  if (comp == "Residual") {
    param_candidates <- c(
      paste0("sigma_", dim_in_params),
      paste0("sigma_", gsub(" ", "_", dim_in_params)),
      paste0("b_sigma_", dim_in_params),
      paste0("sigma_", dim),
      paste0("b_sigma_", dim)
    )
  } else {
    clean_comp <- gsub(":", "__", comp)
    param_candidates <- c(
      paste0("sd_", clean_comp, "__", dim_in_params, "_Intercept"),
      paste0("sd_", clean_comp, "__", dim_in_params),
      paste0("sd_", comp, "__", dim_in_params, "_Intercept"),
      paste0("sd_", comp, "__", dim_in_params),
      paste0("sd_", clean_comp, "__", dim, "_Intercept"),
      paste0("sd_", clean_comp, "__", dim),
      paste0("sd_", comp, "__", dim, "_Intercept"),
      paste0("sd_", comp, "__", dim)
    )
  }

  for (param_name in param_candidates) {
    if (param_name %in% colnames(draws)) {
      sd_draws <- posterior::extract_variable(draws, param_name)
      if (comp == "Residual" && grepl("^b_sigma_", param_name)) {
        sd_draws <- exp(sd_draws)
      }
      return(sd_draws)
    }
  }

  base_pattern <- if (comp == "Residual") {
    paste0("sigma_", dim_in_params)
  } else {
    paste0("sd_", gsub(":", "__", comp), "__", dim_in_params)
  }
  matching_params <- grep(paste0("^", base_pattern), colnames(draws), value = TRUE)
  if (length(matching_params) > 0) {
    sd_draws <- posterior::extract_variable(draws, matching_params[1])
    return(sd_draws)
  }

  base_pattern_fallback <- if (comp == "Residual") {
    paste0("sigma_", dim)
  } else {
    paste0("sd_", gsub(":", "__", comp), "__", dim)
  }
  matching_params_fallback <- grep(paste0("^", base_pattern_fallback), colnames(draws), value = TRUE)
  if (length(matching_params_fallback) > 0) {
    sd_draws <- posterior::extract_variable(draws, matching_params_fallback[1])
    return(sd_draws)
  }

  var_estimate <- vc_dim$var[vc_dim$component == comp]
  n_draws <- if (is.matrix(draws)) nrow(draws) else length(draws)
  rep(sqrt(var_estimate), n_draws)
}

extract_sd_draws_univariate <- function(draws, comp, vc) {
  if (comp == "Residual") {
    param_name <- "sigma"
    param_name_alt <- "b_sigma"
    matching_pattern <- "^sigma"
  } else {
    clean_comp <- gsub(":", "__", comp)
    param_name <- paste0("sd_", clean_comp, "__Intercept")
    param_name_alt <- paste0("sd_", clean_comp)
    matching_pattern <- paste0("^sd_", clean_comp)
  }

  if (param_name %in% colnames(draws)) {
    sd_draws <- posterior::extract_variable(draws, param_name)
    if (comp == "Residual" && startsWith(param_name, "b_sigma")) {
      sd_draws <- exp(sd_draws)
    }
    return(sd_draws)
  }

  if (param_name_alt %in% colnames(draws)) {
    sd_draws <- posterior::extract_variable(draws, param_name_alt)
    if (comp == "Residual" && startsWith(param_name_alt, "b_sigma")) {
      sd_draws <- exp(sd_draws)
    }
    return(sd_draws)
  }

  matching_params <- grep(matching_pattern, colnames(draws), value = TRUE)
  if (length(matching_params) > 0) {
    sd_draws <- posterior::extract_variable(draws, matching_params[1])
    if (comp == "Residual" && startsWith(matching_params[1], "b_sigma")) {
      sd_draws <- exp(sd_draws)
    }
    return(sd_draws)
  }

  var_estimate <- vc$var[vc$component == comp]
  n_draws <- if (is.matrix(draws)) nrow(draws) else length(draws)
  rep(sqrt(var_estimate), n_draws)
}

#' Calculate Posterior Coefficients for a Single Sample Size
#'
#' @param vc_draws Named list of variance component draws.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets (or NULL).
#' @param residual_is Character string for residual composition.
#' @param gstudy_obj Original gstudy object.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @param use_scaled Logical indicating whether to use scaled variance draws.
#'   When FALSE, uses original variance draws and does not scale universe components.
#'   When TRUE, uses scaled variance draws and scales non-object universe components.
#' @param cut_score Optional cut score for phi-cut calculation.
#' @param mu_y Optional grand mean for phi-cut calculation.
#'
#' @return List with 'summary' (data frame of means) and 'distributions' (list of vectors).
#'
#' @keywords internal
calculate_single_posterior <- function(vc_draws, n, object_spec, universe_spec,
                                       error_spec, agg_facets, residual_is, gstudy_obj,
                                       n_provided = FALSE, use_scaled = TRUE,
                                       cut_score = NULL, mu_y = NULL,
                                       ci = NULL, probs = c(0.025, 0.975)) {
  n_draws <- length(vc_draws[[1]])

  # Default universe to object if not specified
  if (is.null(universe_spec) || length(universe_spec) == 0) {
    universe_spec <- object_spec
  }

  # Scale variance components for D-study (only if use_scaled)
  if (use_scaled) {
    scaled_draws <- scale_variance_draws(vc_draws, n, object_spec, agg_facets,
      residual_is, gstudy_obj, n_provided)
  } else {
    scaled_draws <- vc_draws
  }

  # Calculate universe score variance (uni)
  # Object components are NEVER scaled, non-object universe components are scaled only if use_scaled
  uni <- rep(0, n_draws)
  for (comp in universe_spec) {
    if (comp %in% names(scaled_draws)) {
      if (comp %in% object_spec) {
        # Object component: NEVER scaled
        uni <- uni + scaled_draws[[comp]]
      } else {
        # Non-object universe component
        if (use_scaled) {
          # Scale by sample sizes
          scale_factor <- get_universe_scale_factor_for_component(comp, n, agg_facets, residual_is, comp)
          uni <- uni + scaled_draws[[comp]] / scale_factor
        } else {
          # No scaling for unscaled estimate
          uni <- uni + scaled_draws[[comp]]
        }
      }
    }
  }

  # Calculate relative error variance (sigma2_delta)
  sigma2_delta <- compute_relative_error_draws(scaled_draws, universe_spec,
    error_spec, agg_facets, residual_is)

  # Calculate absolute error variance (sigma2_delta_abs)
  sigma2_delta_abs <- compute_absolute_error_draws(scaled_draws, universe_spec,
    error_spec, agg_facets, residual_is)

  # Calculate coefficients
  g <- uni / (uni + sigma2_delta)
  phi <- uni / (uni + sigma2_delta_abs)

  phi_cut <- NA_real_
  if (!is.null(cut_score) && !is.null(mu_y)) {
    adjustment <- (mu_y - cut_score)^2
    phi_cut <- (uni + adjustment) / (uni + sigma2_delta_abs + adjustment)
  }

  sem_rel <- sqrt(sigma2_delta)
  sem_abs <- sqrt(sigma2_delta_abs)

  # Handle infinite/NaN cases
  g[is.nan(g) | is.infinite(g)] <- NA
  phi[is.nan(phi) | is.infinite(phi)] <- NA
  if (!all(is.na(phi_cut))) {
    phi_cut[is.nan(phi_cut) | is.infinite(phi_cut)] <- NA
  }

  ci_cols <- list()
  if (!is.null(ci)) {
    if ("g" %in% ci) {
      ci_cols$g_LL <- quantile(g, probs[1], na.rm = TRUE)
      ci_cols$g_UL <- quantile(g, probs[2], na.rm = TRUE)
    }
    if ("phi" %in% ci) {
      ci_cols$phi_LL <- quantile(phi, probs[1], na.rm = TRUE)
      ci_cols$phi_UL <- quantile(phi, probs[2], na.rm = TRUE)
    }
    if ("phi-cut" %in% ci && !is.null(cut_score)) {
      ci_cols$phi_cut_LL <- quantile(phi_cut, probs[1], na.rm = TRUE)
      ci_cols$phi_cut_UL <- quantile(phi_cut, probs[2], na.rm = TRUE)
    }
  }

  summary_df <- data.frame(
    uni = mean(uni, na.rm = TRUE),
    sigma2_delta = mean(sigma2_delta, na.rm = TRUE),
    sigma2_delta_abs = mean(sigma2_delta_abs, na.rm = TRUE),
    g = mean(g, na.rm = TRUE),
    phi = mean(phi, na.rm = TRUE),
    sem_rel = mean(sem_rel, na.rm = TRUE),
    sem_abs = mean(sem_abs, na.rm = TRUE)
  )

  if (!is.null(cut_score)) {
    summary_df$phi_cut <- mean(phi_cut, na.rm = TRUE)
  }

  if (length(ci_cols) > 0) {
    summary_df <- cbind(summary_df, as.data.frame(ci_cols))
  }

  distributions <- list(
    uni = uni,
    sigma2_delta = sigma2_delta,
    sigma2_delta_abs = sigma2_delta_abs,
    g = g,
    phi = phi,
    sem_rel = sem_rel,
    sem_abs = sem_abs
  )

  if (!is.null(cut_score)) {
    distributions$phi_cut <- phi_cut
  }

  list(
    summary = summary_df,
    distributions = distributions
  )
}

#' Scale Variance Draws for D-Study
#'
#' Applies D-study scaling to posterior draws of variance components.
#'
#' @param vc_draws Named list of variance component draws.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param agg_facets Character vector of aggregation facets.
#' @param residual_is Character string for residual composition.
#' @param gstudy_obj Original gstudy object.
#' @param n_provided Logical indicating if n was explicitly provided.
#'
#' @return Named list of scaled variance draws.
#'
#' @keywords internal
scale_variance_draws <- function(vc_draws, n, object_spec, agg_facets,
  residual_is, gstudy_obj, n_provided = FALSE) {
  scaled_draws <- vc_draws

  do_rescaling <- !is.null(agg_facets) || (n_provided && length(n) > 0)

  if (!do_rescaling) {
    return(scaled_draws)
  }

  total_n <- if (length(n) > 0) prod(unlist(n)) else 1

  for (comp in names(scaled_draws)) {
    if (comp %in% object_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      if (!is.null(agg_facets)) {
        agg_n <- 1
        for (f in agg_facets) {
          if (f %in% names(n)) agg_n <- agg_n * n[[f]]
        }
        scaled_draws[[comp]] <- scaled_draws[[comp]] / agg_n
      } else {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / total_n
      }
      next
    }

    comp_facets <- parse_component_facets(comp)

    has_agg <- !is.null(agg_facets) && any(agg_facets %in% comp_facets)

    if (has_agg) {
      agg_in_comp <- comp_facets[comp_facets %in% agg_facets]
      agg_n <- 1
      for (f in agg_in_comp) {
        if (f %in% names(n)) agg_n <- agg_n * n[[f]]
      }
      if (agg_n > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / agg_n
      }
    } else {
      scale_factor <- 1
      for (facet in comp_facets) {
        if (facet %in% names(n)) {
          scale_factor <- scale_factor * n[[facet]]
        }
      }
      if (scale_factor > 1) {
        scaled_draws[[comp]] <- scaled_draws[[comp]] / scale_factor
      }
    }
  }

  if (!is.null(agg_facets) && length(agg_facets) > 0 && !is.null(residual_is)) {
    for (comp in names(scaled_draws)) {
      if (comp == "Residual" || grepl("^Residual", comp)) {
        residual_facets <- parse_component_facets(residual_is)
        intersect_facets <- agg_facets[agg_facets %in% residual_facets]
        if (length(intersect_facets) > 0) {
          agg_factor <- 1
          for (facet in intersect_facets) {
            if (facet %in% names(n)) {
              agg_factor <- agg_factor * n[[facet]]
            }
          }
          if (agg_factor > 1) {
            scaled_draws[[comp]] <- scaled_draws[[comp]] * agg_factor
          }
        }
      }
    }
  }

  scaled_draws
}

#' Compute Relative Error Variance from Draws
#'
#' @param scaled_draws Named list of scaled variance draws.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets.
#'
#' @return Numeric vector of relative error variances.
#'
#' @keywords internal
compute_relative_error_draws <- function(scaled_draws, universe_spec, error_spec,
  agg_facets, residual_is = NULL) {
  n_draws <- length(scaled_draws[[1]])

  if (!is.null(error_spec)) {
    result <- rep(0, n_draws)
    for (comp in error_spec) {
      if (comp %in% names(scaled_draws)) {
        result <- result + scaled_draws[[comp]]
      }
    }
    return(result)
  }

  result <- rep(0, n_draws)

  for (comp in names(scaled_draws)) {
    if (comp %in% universe_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      result <- result + scaled_draws[[comp]]
      next
    }

    comp_facets <- parse_component_facets(comp)
    is_agg_main <- length(comp_facets) == 1 && comp_facets %in% agg_facets
    if (is_agg_main) {
      next
    }

    if (any(universe_spec %in% comp_facets)) {
      result <- result + scaled_draws[[comp]]
    }
  }

  result
}

#' Compute Absolute Error Variance from Draws
#'
#' @param scaled_draws Named list of scaled variance draws.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of error components (or NULL).
#' @param agg_facets Character vector of aggregation facets.
#'
#' @return Numeric vector of absolute error variances.
#'
#' @keywords internal
compute_absolute_error_draws <- function(scaled_draws, universe_spec, error_spec,
  agg_facets, residual_is = NULL) {
  n_draws <- length(scaled_draws[[1]])

  if (!is.null(error_spec)) {
    result <- rep(0, n_draws)
    for (comp in error_spec) {
      if (comp %in% names(scaled_draws)) {
        result <- result + scaled_draws[[comp]]
      }
    }
    return(result)
  }

  result <- rep(0, n_draws)

  for (comp in names(scaled_draws)) {
    if (comp %in% universe_spec) {
      next
    }

    if (comp == "Residual" || grepl("^Residual", comp)) {
      result <- result + scaled_draws[[comp]]
      next
    }

    comp_facets <- parse_component_facets(comp)
    is_agg_main <- length(comp_facets) == 1 && comp_facets %in% agg_facets
    if (is_agg_main) {
      next
    }

    result <- result + scaled_draws[[comp]]
  }

  result
}

#' @export
extract_grand_mean <- function(gstudy_obj) {
  UseMethod("extract_grand_mean")
}

#' @export
extract_grand_mean.gstudy <- function(gstudy_obj) {
  model <- gstudy_obj$model
  backend <- gstudy_obj$backend

  if (backend == "lme4") {
    fe <- lme4::fixef(model)
    if ("(Intercept)" %in% names(fe)) {
      return("(Intercept)" = fe["(Intercept)"])
    } else if ("Intercept" %in% names(fe)) {
      return(Intercept = fe["Intercept"])
    }
  } else if (backend == "brms") {
    fe <- brms::fixef(model)
    if ("Intercept" %in% rownames(fe)) {
      return(Intercept = fe["Intercept", "Estimate"])
    }
  } else if (backend == "mom") {
    response <- gstudy_obj$response
    if (!is.null(response) && response %in% names(gstudy_obj$data)) {
      return(Intercept = mean(gstudy_obj$data[[response]], na.rm = TRUE))
    }
  }

  stop("Could not extract grand mean from model", call. = FALSE)
}

#' @export
extract_grand_mean.mgstudy <- function(gstudy_obj) {
  model <- gstudy_obj$model
  dimensions <- gstudy_obj$dimensions
  backend <- gstudy_obj$backend

  mu_y <- list()

  if (backend == "brms") {
    fe <- brms::fixef(model)
    for (d in dimensions) {
      param_name <- paste0(d, "_Intercept")
      if (param_name %in% rownames(fe)) {
        mu_y[[d]] <- fe[param_name, "Estimate"]
      } else {
        mu_y[[d]] <- mean(gstudy_obj$data[[d]], na.rm = TRUE)
      }
    }
  } else if (backend == "mom") {
    for (d in dimensions) {
      mu_y[[d]] <- mean(gstudy_obj$data[[d]], na.rm = TRUE)
    }
  } else {
    stop("Multivariate models require brms or mom backend", call. = FALSE)
  }

  return(mu_y)
}
