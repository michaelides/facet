#' Extract Latent Scores from G-Study Objects
#'
#' Creates latent scores (universe score estimates) for the object of measurement
#' based on random effects extracted from a gstudy or mgstudy object. The latent
#' scores are formed on the basis of the universe defined in the gstudy/mgstudy,
#' or a user-specified universe.
#'
#' @param x A gstudy or mgstudy object.
#' @param universe Specification for components that contribute to the universe score.
#' Can be:
#' \itemize{
#'   \item NULL (default): universe includes only the object of measurement
#'   \item A character vector: c("Person", "Person:Rater")
#'   \item A formula: ~ Person + Person:Rater
#' }
#' The object of measurement is always included. Interaction terms in the universe
#' create separate variables for each level combination.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with:
#' \itemize{
#'   \item A column for the object of measurement levels
#'   \item A `latent` column (univariate) or `latent_<dim>` columns (multivariate)
#'   \item Additional columns for interaction combinations if specified in universe
#' }
#'
#' @details
#' ## Universe Specification
#'
#' By default, the latent score is based only on the object of measurement's
#' random effect (e.g., Person). You can expand the universe to include
#' interactions with the object (e.g., Person:Rater).
#'
#' When interaction terms are included in the universe, separate latent score
#' variables are created for each level of the non-object facet. For example,
#' if universe includes "Person:Rater" with raters A, B, C, the output will
#' have columns `latent_RaterA`, `latent_RaterB`, `latent_RaterC` in addition
#' to the base `latent` column.
#'
#' ## Multivariate Models
#'
#' For mgstudy objects, latent scores are computed separately for each dimension
#' and returned in wide format with separate columns per dimension.
#'
#' @seealso [gstudy()], [ranef()]
#'
#' @export
#'
#' @examples
#' # Basic usage - object of measurement only
#' g <- gstudy(Score ~ (1 | Person) + (1 | Rater), data = brennan)
#' latent_scores <- latent(g)
#' head(latent_scores)
#'
#' # With interaction in universe
#' \dontrun{
#' latent_scores <- latent(g, universe = c("Person", "Person:Rater"))
#' }
latent <- function(x, universe = NULL, ...) {
  UseMethod("latent")
}

#' @rdname latent
#' @export
latent.gstudy <- function(x, universe = NULL, ...) {
  if (!inherits(x, "gstudy") && !inherits(x, "mgstudy")) {
    stop("'x' must be a gstudy or mgstudy object", call. = FALSE)
  }

  object <- x$object
  if (is.null(object)) {
    stop("Cannot determine object of measurement from gstudy object", call. = FALSE)
  }

  universe_spec <- parse_universe_spec(universe, object)

  if (x$backend == "brms") {
    re <- ranef(x, summary = FALSE)
  } else {
    re <- ranef(x)
  }

  is_multivariate <- inherits(x, "mgstudy")

  if (is_multivariate) {
    build_latent_multivariate(x, re, object, universe_spec)
  } else {
    build_latent_univariate(x, re, object, universe_spec)
  }
}

#' @rdname latent
#' @export
latent.mgstudy <- latent.gstudy

parse_universe_spec <- function(universe, object) {
  if (is.null(universe)) {
    return(object)
  }

  spec <- parse_specification(universe)

  if (!object %in% spec) {
    warning(
      "The universe specification did not include the object of measurement '",
      object, "'. Adding '", object, "' to the universe.",
      call. = FALSE
    )
    spec <- c(object, spec)
  }

  spec
}

build_latent_univariate <- function(x, re, object, universe_spec) {
  re_object <- extract_object_random_effects(re, object, is_multivariate = FALSE, backend = x$backend)

  if (is.null(re_object) || length(re_object) == 0) {
    warning("No random effects found for object '", object, "'", call. = FALSE)
    return(data.frame())
  }

  if (is.list(re_object) && !is.data.frame(re_object)) {
    object_levels <- names(re_object)
    latent_scores <- unlist(re_object)
    names(latent_scores) <- NULL
  } else if (is.data.frame(re_object)) {
    object_levels <- rownames(re_object)
    latent_scores <- re_object[, 1]
  } else {
    object_levels <- names(re_object)
    latent_scores <- as.numeric(re_object)
  }

  result <- data.frame(
    object = object_levels,
    latent = latent_scores,
    stringsAsFactors = FALSE
  )
  names(result)[1] <- object

  interaction_specs <- universe_spec[grepl(":", universe_spec)]

  for (spec in interaction_specs) {
    result <- add_interaction_latent(result, re, spec, object, x$data, is_multivariate = FALSE, dimension = NULL, backend = x$backend)
  }

  result
}

build_latent_multivariate <- function(x, re, object, universe_spec) {
  dimensions <- x$dimensions

  if (is.null(dimensions) || length(dimensions) == 0) {
    warning("No dimensions found in mgstudy object", call. = FALSE)
    return(data.frame())
  }

  results_by_dim <- list()
  object_levels <- NULL

  for (dim in dimensions) {
    re_dim <- extract_object_random_effects(re, object, is_multivariate = TRUE, dimension = dim, backend = x$backend)

    if (is.null(re_dim) || length(re_dim) == 0) {
      warning("No random effects found for object '", object, "' in dimension '", dim, "'", call. = FALSE)
      next
    }

    if (is.list(re_dim) && !is.data.frame(re_dim)) {
      if (is.null(object_levels)) {
        object_levels <- names(re_dim)
      }
      latent_scores <- unlist(re_dim)
      names(latent_scores) <- NULL
    } else if (is.data.frame(re_dim)) {
      if (is.null(object_levels)) {
        object_levels <- rownames(re_dim)
      }
      latent_scores <- re_dim[, 1]
    } else {
      if (is.null(object_levels)) {
        object_levels <- names(re_dim)
      }
      latent_scores <- as.numeric(re_dim)
    }

    results_by_dim[[dim]] <- latent_scores
  }

  if (length(results_by_dim) == 0) {
    return(data.frame())
  }

  result <- data.frame(object = object_levels, stringsAsFactors = FALSE)
  names(result)[1] <- object

  for (dim in names(results_by_dim)) {
    col_name <- paste0("latent_", dim)
    result[[col_name]] <- results_by_dim[[dim]]
  }

  interaction_specs <- universe_spec[grepl(":", universe_spec)]

  for (spec in interaction_specs) {
    for (dim in dimensions) {
      result <- add_interaction_latent(result, re, spec, object, x$data, is_multivariate = TRUE, dimension = dim, backend = x$backend)
    }
  }

  result
}

extract_object_random_effects <- function(re, object, is_multivariate, dimension = NULL, backend = "lme4") {
  if (backend == "brms") {
    if (!object %in% names(re)) {
      return(NULL)
    }
    re_array <- re[[object]]

    if (is.null(re_array) || length(dim(re_array)) < 3) {
      return(NULL)
    }

    if (is_multivariate && !is.null(dimension)) {
      param_name <- paste0(dimension, "_Intercept")
      if (!param_name %in% dimnames(re_array)[[3]]) {
        return(NULL)
      }
      param_idx <- which(dimnames(re_array)[[3]] == param_name)
      estimates <- apply(re_array[, , param_idx], 2, mean)
      names(estimates) <- dimnames(re_array)[[2]]
      return(estimates)
    } else {
      estimates <- apply(re_array[, , 1], 2, mean)
      names(estimates) <- dimnames(re_array)[[2]]
      return(estimates)
    }
  }

  if (is_multivariate && !is.null(dimension)) {
    if (!dimension %in% names(re)) {
      return(NULL)
    }
    re_dim <- re[[dimension]]
    if (!object %in% names(re_dim)) {
      return(NULL)
    }
    return(re_dim[[object]])
  }

  if (!object %in% names(re)) {
    return(NULL)
  }

  re[[object]]
}

add_interaction_latent <- function(result, re, spec, object, data, is_multivariate, dimension, backend = "lme4") {
  parts <- strsplit(spec, ":")[[1]]
  parts <- trimws(parts)

  if (!object %in% parts) {
    warning(
      "Interaction '", spec, "' does not include object '", object, "'. Skipping.",
      call. = FALSE
    )
    return(result)
  }

  other_facet <- parts[parts != object]
  if (length(other_facet) != 1) {
    warning(
      "Interaction '", spec, "' has unexpected structure. Skipping.",
      call. = FALSE
    )
    return(result)
  }
  other_facet <- other_facet[1]

  if (is.null(data) || !other_facet %in% names(data)) {
    warning(
      "Cannot find facet '", other_facet, "' in data. Skipping interaction '",
      spec, "'.",
      call. = FALSE
    )
    return(result)
  }

  other_levels <- unique(data[[other_facet]])
  other_levels <- sort(other_levels)

  re_interaction <- NULL

  if (backend == "brms") {
    if (spec %in% names(re)) {
      re_interaction <- re[[spec]]
    } else {
      alt_spec <- paste(rev(parts), collapse = ":")
      if (alt_spec %in% names(re)) {
        re_interaction <- re[[alt_spec]]
      }
    }
  } else if (is_multivariate && !is.null(dimension)) {
    if (dimension %in% names(re)) {
      re_dim <- re[[dimension]]
      if (spec %in% names(re_dim)) {
        re_interaction <- re_dim[[spec]]
      } else {
        alt_spec <- paste(rev(parts), collapse = ":")
        if (alt_spec %in% names(re_dim)) {
          re_interaction <- re_dim[[alt_spec]]
        }
      }
    }
  } else {
    if (spec %in% names(re)) {
      re_interaction <- re[[spec]]
    } else {
      alt_spec <- paste(rev(parts), collapse = ":")
      if (alt_spec %in% names(re)) {
        re_interaction <- re[[alt_spec]]
      }
    }
  }

  if (is.null(re_interaction)) {
    warning(
      "No random effects found for interaction '", spec, "'. Skipping.",
      call. = FALSE
    )
    return(result)
  }

  object_levels <- result[[object]]

  if (backend == "brms") {
    if (is_multivariate && !is.null(dimension)) {
      param_name <- paste0(dimension, "_Intercept")
      if (!param_name %in% dimnames(re_interaction)[[3]]) {
        warning(
          "No parameter '", param_name, "' found for interaction '", spec, "'. Skipping.",
          call. = FALSE
        )
        return(result)
      }
      param_idx <- which(dimnames(re_interaction)[[3]] == param_name)
      estimates <- apply(re_interaction[, , param_idx], 2, mean)
      names(estimates) <- dimnames(re_interaction)[[2]]
      re_interaction <- estimates
    } else {
      interaction_names <- dimnames(re_interaction)[[2]]
      estimates <- apply(re_interaction[, , 1], 2, mean)
      names(estimates) <- interaction_names
      re_interaction <- estimates
    }
  } else if (is.data.frame(re_interaction)) {
    interaction_names <- rownames(re_interaction)
  } else {
    interaction_names <- names(re_interaction)
  }

  for (level in other_levels) {
    level_str <- as.character(level)
    col_name <- if (is_multivariate && !is.null(dimension)) {
      paste0("latent_", level_str, "_", dimension)
    } else {
      paste0("latent_", level_str)
    }

    latent_values <- numeric(nrow(result))

    for (i in seq_along(object_levels)) {
      obj_level <- as.character(object_levels[i])

      patterns_to_try <- c(
        paste0(obj_level, ".", level_str),
        paste0(level_str, ".", obj_level),
        paste0(obj_level, ":", level_str),
        paste0(level_str, ":", obj_level)
      )

      matched <- FALSE
      for (pattern in patterns_to_try) {
        if (pattern %in% interaction_names) {
          if (is.data.frame(re_interaction)) {
            latent_values[i] <- re_interaction[pattern, 1]
          } else {
            latent_values[i] <- re_interaction[pattern]
          }
          matched <- TRUE
          break
        }
      }

      if (!matched) {
        latent_values[i] <- NA_real_
      }
    }

    result[[col_name]] <- latent_values
  }

  result
}
