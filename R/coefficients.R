#' Generalizability and Dependability Coefficients
#'
#' Functions for computing G (generalizability) and D (dependability) coefficients
#' from variance components in generalizability theory.
#'
#' @name coefficients
#' @family decision studies
#' @keywords internal
#' @importFrom magrittr %>%
NULL

#' Compute Scale Factor for D-Study Variance Component
#'
#' Calculates the scaling factor for a variance component based on
#' the D-study sample sizes. For main effects, divides by n_facet.
#' For interactions, divides by the product of n for each facet.
#'
#' @param facets Character vector of facet names in the component.
#' @param n Named list of sample sizes for each facet.
#' @return The scale factor (numeric).
#'
#' @keywords internal
compute_scale_factor_from_facets <- function(facets, n) {
  if (length(facets) == 0 || identical(facets, "Residual")) {
    return(1)
  }

  scale_factor <- 1
  for (facet in facets) {
    if (facet %in% names(n)) {
      scale_factor <- scale_factor * n[[facet]]
    }
  }
  scale_factor
}

#' Compute Scale Factor for a Component
#'
#' Computes the scale factor for a variance component.
#' Used in composite variance calculations.
#'
#' @param comp The variance component name.
#' @param n Named list of sample sizes.
#' @param object_spec Character vector of object components.
#' @param n_provided Logical indicating if n was explicitly provided.
#' @return The scale factor for this component.
#'
#' @keywords internal
compute_component_scale_factor <- function(comp, n, object_spec, n_provided) {
  if (comp %in% object_spec) {
    return(1)
  }

  total_n <- if (length(n) > 0) prod(unlist(n)) else 1

  if (n_provided && length(n) > 0) {
    if (comp == "Residual") {
      return(total_n)
    }

    comp_facets <- parse_component_facets(comp)
    scale_factor <- 1
    for (f in comp_facets) {
      if (f %in% names(n)) {
        scale_factor <- scale_factor * n[[f]]
      }
    }
    return(scale_factor)
  }

  return(1)
}

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
    # Universe score is always estimated from unscaled (G-study) variance components
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

        uni <- uni + comp_var
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
  components <- vc$component
  classification <- classify_error_components(components, universe_spec, error_spec, agg_facets)
  sum(vc$var[classification$is_relative_error], na.rm = TRUE)
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
  components <- vc$component
  classification <- classify_error_components(components, universe_spec, error_spec, agg_facets)
  sum(vc$var[classification$is_absolute_error], na.rm = TRUE)
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
  classify_error_components(components, universe_spec, error_spec, agg_facets)
}

#' Classify Error Components for D-Study
#'
#' Unified function that determines which variance components belong to
#' relative and absolute error, consolidating the logic previously
#' duplicated across `compute_relative_error()`, `compute_absolute_error()`,
#' and `identify_error_components_for_draws()`.
#'
#' @param components Character vector of component names.
#' @param universe_spec Character vector of universe components.
#' @param error_spec Character vector of explicit error components (or NULL for default).
#' @param agg_facets Character vector of aggregation facets (or NULL).
#'
#' @return A list with:
#' \item{is_relative_error}{Logical vector: TRUE for components in relative error}
#' \item{is_absolute_error}{Logical vector: TRUE for components in absolute error}
#' \item{error_spec_used}{Logical: TRUE if an explicit error_spec was provided}
#'
#' @keywords internal
classify_error_components <- function(components, universe_spec, error_spec = NULL,
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
