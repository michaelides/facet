#' Simulate Data for Generalizability Theory Analysis
#'
#' Generate simulated data for univariate and multivariate generalizability theory
#' (G-theory) analyses. Supports both wide-format and long-format multivariate
#' designs, with optional correlation structures for residual and random effects.
#'
#' @param facets A named list specifying the number of levels for each facet.
#'   For example, `list(person = 50, item = 10)` creates a crossed design with
#'   50 persons and 10 items. Names should be valid R variable names.
#' @param formula An optional lme4-style formula from which to derive the facet
#'   structure. If provided, overrides `facets`. The formula should use
#'   random effects syntax, e.g., `y ~ (1|person) + (1|item)`.
#' @param vc A named list of variance components (as variances, not SDs).
#'   Names should match facet names or interactions (use `:` for interactions).
#'   For example, `list(person = 1, item = 0.5, "person:item" = 0.3)`.
#' @param sd_residual Standard deviation of the residual error term. Default is 1.
#' @param multivariate Logical. If TRUE, generates multivariate data. Default is FALSE.
#' @param n_dims Number of dimensions (outcome variables) for multivariate data.
#'   Default is 2. Ignored if `multivariate = FALSE`.
#' @param dim_names Character vector of dimension names. If NULL, dimensions
#'   are named `y1`, `y2`, etc. for wide format.
#' @param format Character string specifying the output format: "wide" or "long".
#'   Wide format has multiple outcome columns; long format has a dimension
#'   indicator variable. Default is "wide".
#' @param dim_var Name of the dimension indicator variable for long format.
#'   Default is "dimension".
#' @param response_var Name of the response variable. Default is "y" for
#'   univariate and "score" for multivariate long format.
#' @param residual_cor A correlation matrix for residual errors across dimensions.
#'   Must be a symmetric positive-definite matrix with dimensions matching `n_dims`.
#'   NULL (default) implies independent residuals (identity matrix).
#' @param re_cor A named list of correlation matrices for random effects.
#'   Names should match facet names. Each matrix should be symmetric positive-definite
#'   with dimensions matching `n_dims`. NULL (default) implies independent random effects.
#' @param nested A named list specifying nesting relationships. Names are nested
#'   facets, values are the facets they are nested within. For example,
#'   `list(item = "person")` means items are nested within persons.
#' @param seed Optional random seed for reproducibility.
#' @param prefix Prefix for factor level names. Default is empty string,
#'   producing levels like "1", "2", etc.
#'
#' @return A data frame with:
#'   \itemize{
#'     \item Columns for each facet (as factors)
#'     \item One or more outcome columns (univariate: `y`, wide multivariate: `y1`, `y2`, ...,
#'           long multivariate: `score` with `dimension` column)
#'     \item For long format: a dimension indicator column
#'   }
#'
#' @details
#' ## Variance Component Specification
#'
#' Variance components can be specified using either the `vc` parameter or
#' individual `sd_*` parameters. The `vc` parameter takes a named list where
#' names correspond to facet names or interactions (use `:` to separate).
#'
#' ## Univariate Simulation
#'
#' For univariate designs, the function generates data according to:
#' \deqn{y_{ijk...} = \mu + \epsilon_p + \epsilon_i + \epsilon_{pi} + \epsilon_{residual}}
#' where each \eqn{\epsilon} term is drawn from a normal distribution with
#' the specified variance component.
#'
#' ## Multivariate Simulation
#'
#' For multivariate designs, correlations between dimensions can be specified:
#' - **Residual correlations**: Correlations between residual errors across dimensions
#' - **Random effect correlations**: Correlations between random effects for a specific facet
#'
#' Correlations are implemented using Cholesky decomposition to ensure
#' proper multivariate normal distributions.
#'
#' ## Format Options
#'
#' - **Wide format**: Multiple outcome columns (`y1`, `y2`, ...) - suitable for
#'   `mvbind()` syntax in `gstudy()`
#' - **Long format**: Single outcome column with dimension indicator - suitable for
#'   brms long-format models like `score ~ 0 + dimension + (0+dimension|facet)`
#'
#' @seealso [gstudy()] for conducting G-studies on the simulated data
#'
#' @examples
#' # Univariate p x i design
#' data_uni <- simulate_gtheory_data(
#'   facets = list(person = 20, item = 10),
#'   vc = list(person = 1, item = 0.5, "person:item" = 0.3),
#'   sd_residual = 1
#' )
#'
#' # Wide-format multivariate (2 dimensions)
#' data_wide <- simulate_gtheory_data(
#'   facets = list(person = 50, item = 10),
#'   vc = list(person = 1, item = 0.5),
#'   multivariate = TRUE,
#'   n_dims = 2,
#'   format = "wide"
#' )
#'
#' # Long-format multivariate with correlations
#' cor_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' data_long <- simulate_gtheory_data(
#'   facets = list(person = 50, item = 10),
#'   vc = list(person = 1, item = 0.5),
#'   multivariate = TRUE,
#'   n_dims = 2,
#'   format = "long",
#'   residual_cor = cor_matrix
#' )
#'
#' @export
simulate_gtheory_data <- function(
  facets = NULL,
  formula = NULL,
  vc = NULL,
  sd_residual = 1,
  multivariate = FALSE,
  n_dims = 2,
  dim_names = NULL,
  format = c("wide", "long"),
  dim_var = "dimension",
  response_var = NULL,
  residual_cor = NULL,
  re_cor = NULL,
  nested = NULL,
  seed = NULL,
  prefix = ""
) {

  format <- match.arg(format)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(facets) && is.null(formula)) {
    stop("Either 'facets' or 'formula' must be provided", call. = FALSE)
  }

  if (!is.null(formula)) {
    facets <- extract_facets_from_formula(formula)
  }

  facet_names <- names(facets)
  facet_levels <- unlist(facets)

  if (is.null(vc)) {
    vc <- setNames(as.list(rep(0, length(facet_names))), facet_names)
  }

  if (multivariate) {
    if (is.null(dim_names)) {
      if (format == "wide") {
        dim_names <- paste0("y", seq_len(n_dims))
      } else {
        dim_names <- paste0("dim", seq_len(n_dims))
      }
    }

    if (is.null(response_var)) {
      response_var <- if (format == "long") "score" else "y"
    }

    if (is.null(residual_cor)) {
      residual_cor <- diag(n_dims)
    } else {
      residual_cor <- validate_correlation_matrix(residual_cor, n_dims, "residual_cor")
    }

    if (!is.null(re_cor)) {
      for (facet_name in names(re_cor)) {
        re_cor[[facet_name]] <- validate_correlation_matrix(
          re_cor[[facet_name]], n_dims,
          paste0("re_cor$", facet_name)
        )
      }
    }

    if (format == "wide") {
      data <- generate_multivariate_wide(
        facets = facets,
        vc = vc,
        n_dims = n_dims,
        dim_names = dim_names,
        sd_residual = sd_residual,
        residual_cor = residual_cor,
        re_cor = re_cor,
        nested = nested,
        prefix = prefix
      )
    } else {
      data <- generate_multivariate_long(
        facets = facets,
        vc = vc,
        n_dims = n_dims,
        dim_names = dim_names,
        dim_var = dim_var,
        response_var = response_var,
        sd_residual = sd_residual,
        residual_cor = residual_cor,
        re_cor = re_cor,
        nested = nested,
        prefix = prefix
      )
    }
  } else {
    if (is.null(response_var)) {
      response_var <- "y"
    }

    data <- generate_univariate(
      facets = facets,
      vc = vc,
      sd_residual = sd_residual,
      response_var = response_var,
      nested = nested,
      prefix = prefix
    )
  }

  data
}

#' Extract Facet Names and Levels from Formula
#'
#' @param formula An lme4-style formula
#' @return Named list of facet levels (all set to 10 as default)
#' @keywords internal
extract_facets_from_formula <- function(formula) {
  formula_char <- deparse(formula)

  random_pattern <- "\\(\\s*1\\s*\\|\\s*([^)]+)\\)"
  random_terms <- regmatches(formula_char, gregexpr(random_pattern, formula_char))[[1]]

  if (length(random_terms) == 0) {
    stop("No random effect terms found in formula", call. = FALSE)
  }

  facets <- list()
  for (term in random_terms) {
    facet_spec <- sub(random_pattern, "\\1", term)
    facet_spec <- trimws(facet_spec)

    if (grepl("\\|", facet_spec)) {
      parts <- strsplit(facet_spec, "\\|")[[1]]
      facet_spec <- trimws(parts[length(parts)])
    }

    facets[[facet_spec]] <- 10
  }

  facets
}

#' Validate Correlation Matrix
#'
#' @param cor_matrix A matrix to validate
#' @param n_dims Expected number of dimensions
#' @param param_name Parameter name for error messages
#' @return The validated correlation matrix
#' @keywords internal
validate_correlation_matrix <- function(cor_matrix, n_dims, param_name) {
  if (!is.matrix(cor_matrix)) {
    cor_matrix <- as.matrix(cor_matrix)
  }

  if (nrow(cor_matrix) != n_dims || ncol(cor_matrix) != n_dims) {
    stop(
      param_name, " must be a ", n_dims, "x", n_dims, " matrix",
      call. = FALSE
    )
  }

  if (!isSymmetric(cor_matrix)) {
    stop(param_name, " must be a symmetric matrix", call. = FALSE)
  }

  eigen_values <- eigen(cor_matrix, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigen_values <= 0)) {
    stop(param_name, " must be positive definite", call. = FALSE)
  }

  if (!all(diag(cor_matrix) == 1)) {
    warning("Diagonal of ", param_name, " should be 1 for correlation matrices", call. = FALSE)
  }

  cor_matrix
}

#' Generate Design Structure
#'
#' @param facets Named list of facet levels
#' @param nested Named list of nesting relationships
#' @param prefix Prefix for level names
#' @return Data frame with facet columns
#' @keywords internal
generate_design_structure <- function(facets, nested = NULL, prefix = "") {
  facet_names <- names(facets)
  facet_levels <- unlist(facets)

  n_total <- prod(facet_levels)

  design <- data.frame(row_id = seq_len(n_total))

  ordered_facets <- facet_names
  if (!is.null(nested)) {
    nested_names <- names(nested)
    nesting_names <- unlist(nested)
    top_level <- setdiff(facet_names, c(nested_names, nesting_names))
    middle_level <- setdiff(nesting_names, nested_names)
    bottom_level <- nested_names

    ordered_facets <- c(top_level, middle_level, bottom_level)
    ordered_facets <- ordered_facets[ordered_facets %in% facet_names]
  }

  for (facet_name in ordered_facets) {
    n_levels <- facets[[facet_name]]

    if (!is.null(nested) && facet_name %in% names(nested)) {
      nesting_facet <- nested[[facet_name]]
      if (nesting_facet %in% names(design)) {
        n_nesting <- facets[[nesting_facet]]
        design[[facet_name]] <- paste0(
          prefix,
          rep(rep(seq_len(n_levels), each = n_total / (n_levels * n_nesting)), n_nesting)
        )
      }
    } else {
      cycle_length <- n_levels
      for (prev_facet in ordered_facets[1:which(ordered_facets == facet_name)]) {
        if (prev_facet != facet_name && !(prev_facet %in% names(nested))) {
          cycle_length <- cycle_length * facets[[prev_facet]]
        }
      }

      n_reps <- n_total / (n_levels * (n_total / cycle_length))
      design[[facet_name]] <- paste0(
        prefix,
        rep(seq_len(n_levels), each = n_total / cycle_length, length.out = n_total)
      )
    }
  }

  design$row_id <- NULL

  for (facet_name in facet_names) {
    if (!is.na(prefix) && prefix != "") {
      design[[facet_name]] <- paste0(prefix, design[[facet_name]])
    }
  }

  design
}

#' Generate Univariate Data
#'
#' @param facets Named list of facet levels
#' @param vc Named list of variance components
#' @param sd_residual Residual standard deviation
#' @param response_var Name of response variable
#' @param nested Named list of nesting relationships
#' @param prefix Prefix for level names
#' @return Data frame with simulated data
#' @keywords internal
generate_univariate <- function(facets, vc, sd_residual, response_var, nested = NULL, prefix = "") {
  facet_names <- names(facets)
  facet_levels <- unlist(facets)
  n_total <- prod(facet_levels)

  design <- generate_design_structure(facets, nested, prefix)

  y <- rnorm(n_total, 0, sd_residual)

  for (facet_name in facet_names) {
    vc_name <- facet_name
    vc_value <- vc[[vc_name]]

    if (is.null(vc_value)) {
      vc_value <- 0
    }

    if (vc_value > 0) {
      n_levels <- facets[[facet_name]]
      facet_effects <- rnorm(n_levels, 0, sqrt(vc_value))

      y <- y + facet_effects[as.integer(factor(design[[facet_name]]))]
    }
  }

  interaction_names <- names(vc)[grepl(":", names(vc))]
  for (int_name in interaction_names) {
    vc_value <- vc[[int_name]]
    if (vc_value > 0) {
      facet_parts <- strsplit(int_name, ":")[[1]]
      facet_parts <- trimws(facet_parts)

      interaction_term <- interaction(
        design[facet_parts],
        drop = TRUE
      )

      n_interaction_levels <- nlevels(interaction_term)
      interaction_effects <- rnorm(n_interaction_levels, 0, sqrt(vc_value))

      y <- y + interaction_effects[interaction_term]
    }
  }

  design[[response_var]] <- y

  for (facet_name in facet_names) {
    design[[facet_name]] <- factor(design[[facet_name]])
  }

  design
}

#' Generate Wide-Format Multivariate Data
#'
#' @param facets Named list of facet levels
#' @param vc Named list of variance components
#' @param n_dims Number of dimensions
#' @param dim_names Names of dimensions
#' @param sd_residual Residual standard deviation (or vector for heterogeneous)
#' @param residual_cor Residual correlation matrix
#' @param re_cor Named list of random effect correlation matrices
#' @param nested Named list of nesting relationships
#' @param prefix Prefix for level names
#' @return Data frame with simulated multivariate data
#' @keywords internal
generate_multivariate_wide <- function(facets, vc, n_dims, dim_names, sd_residual,
                                       residual_cor, re_cor, nested, prefix) {
  facet_names <- names(facets)
  n_total <- prod(unlist(facets))

  design <- generate_design_structure(facets, nested, prefix)

  if (length(sd_residual) == 1) {
    sd_residual <- rep(sd_residual, n_dims)
  }

  y_matrix <- matrix(0, nrow = n_total, ncol = n_dims)

  residual_effects <- generate_correlated_effects(n_total, residual_cor, sd_residual)
  y_matrix <- y_matrix + residual_effects

  for (facet_name in facet_names) {
    vc_name <- facet_name
    vc_value <- vc[[vc_name]]

    if (is.null(vc_value)) {
      vc_value <- 0
    }

    if (vc_value > 0) {
      n_levels <- facets[[facet_name]]

      cor_matrix <- if (!is.null(re_cor) && facet_name %in% names(re_cor)) {
        re_cor[[facet_name]]
      } else {
        diag(n_dims)
      }

      facet_effects <- generate_correlated_effects(n_levels, cor_matrix, rep(sqrt(vc_value), n_dims))

      y_matrix <- y_matrix + facet_effects[as.integer(factor(design[[facet_name]])), ]
    }
  }

  interaction_names <- names(vc)[grepl(":", names(vc))]
  for (int_name in interaction_names) {
    vc_value <- vc[[int_name]]
    if (vc_value > 0) {
      facet_parts <- strsplit(int_name, ":")[[1]]
      facet_parts <- trimws(facet_parts)

      interaction_term <- interaction(
        design[facet_parts],
        drop = TRUE
      )

      n_interaction_levels <- nlevels(interaction_term)
      interaction_effects <- generate_correlated_effects(
        n_interaction_levels,
        diag(n_dims),
        rep(sqrt(vc_value), n_dims)
      )

      y_matrix <- y_matrix + interaction_effects[interaction_term, ]
    }
  }

  colnames(y_matrix) <- dim_names

  for (facet_name in facet_names) {
    design[[facet_name]] <- factor(design[[facet_name]])
  }

  cbind(design, as.data.frame(y_matrix))
}

#' Generate Long-Format Multivariate Data
#'
#' @param facets Named list of facet levels
#' @param vc Named list of variance components
#' @param n_dims Number of dimensions
#' @param dim_names Names of dimensions
#' @param dim_var Name of dimension indicator variable
#' @param response_var Name of response variable
#' @param sd_residual Residual standard deviation (or vector for heterogeneous)
#' @param residual_cor Residual correlation matrix
#' @param re_cor Named list of random effect correlation matrices
#' @param nested Named list of nesting relationships
#' @param prefix Prefix for level names
#' @return Data frame with simulated multivariate data in long format
#' @keywords internal
generate_multivariate_long <- function(facets, vc, n_dims, dim_names, dim_var,
                                       response_var, sd_residual, residual_cor,
                                       re_cor, nested, prefix) {
  facet_names <- names(facets)
  n_base <- prod(unlist(facets))

  design_wide <- generate_design_structure(facets, nested, prefix)

  if (length(sd_residual) == 1) {
    sd_residual <- rep(sd_residual, n_dims)
  }

  y_matrix <- matrix(0, nrow = n_base, ncol = n_dims)

  residual_effects <- generate_correlated_effects(n_base, residual_cor, sd_residual)
  y_matrix <- y_matrix + residual_effects

  for (facet_name in facet_names) {
    vc_name <- facet_name
    vc_value <- vc[[vc_name]]

    if (is.null(vc_value)) {
      vc_value <- 0
    }

    if (vc_value > 0) {
      n_levels <- facets[[facet_name]]

      cor_matrix <- if (!is.null(re_cor) && facet_name %in% names(re_cor)) {
        re_cor[[facet_name]]
      } else {
        diag(n_dims)
      }

      facet_effects <- generate_correlated_effects(n_levels, cor_matrix, rep(sqrt(vc_value), n_dims))

      y_matrix <- y_matrix + facet_effects[as.integer(factor(design_wide[[facet_name]])), ]
    }
  }

  interaction_names <- names(vc)[grepl(":", names(vc))]
  for (int_name in interaction_names) {
    vc_value <- vc[[int_name]]
    if (vc_value > 0) {
      facet_parts <- strsplit(int_name, ":")[[1]]
      facet_parts <- trimws(facet_parts)

      interaction_term <- interaction(
        design_wide[facet_parts],
        drop = TRUE
      )

      n_interaction_levels <- nlevels(interaction_term)
      interaction_effects <- generate_correlated_effects(
        n_interaction_levels,
        diag(n_dims),
        rep(sqrt(vc_value), n_dims)
      )

      y_matrix <- y_matrix + interaction_effects[interaction_term, ]
    }
  }

  design_long <- data.frame(
    design_wide[rep(seq_len(n_base), n_dims), ],
    dim_var = rep(dim_names, each = n_base),
    response_var = as.vector(y_matrix),
    stringsAsFactors = FALSE
  )
  names(design_long)[names(design_long) == "dim_var"] <- dim_var
  names(design_long)[names(design_long) == "response_var"] <- response_var

  for (facet_name in facet_names) {
    design_long[[facet_name]] <- factor(design_long[[facet_name]])
  }
  design_long[[dim_var]] <- factor(design_long[[dim_var]], levels = dim_names)

  rownames(design_long) <- NULL

  design_long
}

#' Generate Correlated Effects Using Cholesky Decomposition
#'
#' @param n Number of observations
#' @param cor_matrix Correlation matrix
#' @param sds Standard deviations (vector of length ncol(cor_matrix))
#' @return Matrix of correlated effects (n x ncol(cor_matrix))
#' @keywords internal
generate_correlated_effects <- function(n, cor_matrix, sds) {
  n_dims <- ncol(cor_matrix)

  if (n_dims == 1) {
    return(matrix(rnorm(n, 0, sds[1]), ncol = 1))
  }

  Z <- matrix(rnorm(n * n_dims), nrow = n, ncol = n_dims)

  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)

  chol_decomp <- chol(cov_matrix)

  Y <- Z %*% chol_decomp

  Y
}

#' Simulate Data with Individual SD Parameters (Legacy Style)
#'
#' A convenience wrapper that uses individual `sd_*` parameters similar to
#' the `makedata()` function style. Useful for complex designs with many
#' variance components.
#'
#' @param n_p Number of persons
#' @param n_i Number of items
#' @param n_d Number of domains (optional)
#' @param n_l Number of levels (optional)
#' @param sd_res Residual standard deviation
#' @param sd_p Person standard deviation
#' @param sd_i Item standard deviation
#' @param sd_d Domain standard deviation
#' @param sd_l Level standard deviation
#' @param sd_pd Person-by-domain interaction SD
#' @param sd_li Level-by-item interaction SD
#' @param sd_ld Level-by-domain interaction SD
#' @param multivariate Logical. If TRUE, generates multivariate data
#' @param n_dims Number of dimensions for multivariate data
#' @param seed Random seed
#'
#' @return A data frame with simulated data
#'
#' @examples
#' # Simulate data similar to the original makedata() example
#' data <- simulate_gtheory_legacy(
#'   n_p = 200, n_i = 25, n_d = 5, n_l = 20,
#'   sd_res = 1, sd_p = 1, sd_i = .5, sd_d = .3, sd_l = .5,
#'   sd_pd = .3, sd_li = .3, sd_ld = .3
#' )
#'
#' @export
simulate_gtheory_legacy <- function(
  n_p = 200,
  n_i = 25,
  n_d = 5,
  n_l = 20,
  sd_res = 1,
  sd_p = 1,
  sd_i = 0.5,
  sd_d = 0.3,
  sd_l = 0.5,
  sd_pd = 0.3,
  sd_li = 0.3,
  sd_ld = 0.3,
  multivariate = FALSE,
  n_dims = 2,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_res <- n_p * n_i

  p <- paste0("p", rep(1:n_p, n_i))
  i <- paste0("i", rep(1:n_i, each = n_p))
  d <- paste0("d", rep(1:n_d, each = n_res / n_d))
  l <- rep(paste0("l", rep(1:n_l, each = n_p / n_l)), n_i)

  p_factor <- factor(p, levels = unique(p))
  i_factor <- factor(i, levels = unique(i))
  d_factor <- factor(d, levels = unique(d))
  l_factor <- factor(l, levels = unique(l))

  if (!multivariate) {
    y <- rnorm(n_res, 0, sd_res)

    y <- y + rnorm(nlevels(p_factor), 0, sd_p)[p_factor]
    y <- y + rnorm(nlevels(i_factor), 0, sd_i)[i_factor]
    y <- y + rnorm(nlevels(d_factor), 0, sd_d)[d_factor]
    y <- y + rnorm(nlevels(l_factor), 0, sd_l)[l_factor]

    pd_factor <- interaction(p, d, drop = TRUE)
    y <- y + rnorm(nlevels(pd_factor), 0, sd_pd)[pd_factor]

    li_factor <- interaction(l, i, drop = TRUE)
    y <- y + rnorm(nlevels(li_factor), 0, sd_li)[li_factor]

    ld_factor <- interaction(l, d, drop = TRUE)
    y <- y + rnorm(nlevels(ld_factor), 0, sd_ld)[ld_factor]

    data.frame(p = p, i = i, d = d, l = l, y = y)
  } else {
    y_matrix <- matrix(0, nrow = n_res, ncol = n_dims)

    residual_cor <- diag(n_dims)
    y_matrix <- y_matrix + generate_correlated_effects(n_res, residual_cor, rep(sd_res, n_dims))

    y_matrix <- y_matrix + generate_correlated_effects(nlevels(p_factor), diag(n_dims), rep(sd_p, n_dims))[p_factor, ]

    y_matrix <- y_matrix + generate_correlated_effects(nlevels(i_factor), diag(n_dims), rep(sd_i, n_dims))[i_factor, ]

    y_matrix <- y_matrix + generate_correlated_effects(nlevels(d_factor), diag(n_dims), rep(sd_d, n_dims))[d_factor, ]

    y_matrix <- y_matrix + generate_correlated_effects(nlevels(l_factor), diag(n_dims), rep(sd_l, n_dims))[l_factor, ]

    pd_factor <- interaction(p, d, drop = TRUE)
    y_matrix <- y_matrix + generate_correlated_effects(nlevels(pd_factor), diag(n_dims), rep(sd_pd, n_dims))[pd_factor, ]

    li_factor <- interaction(l, i, drop = TRUE)
    y_matrix <- y_matrix + generate_correlated_effects(nlevels(li_factor), diag(n_dims), rep(sd_li, n_dims))[li_factor, ]

    ld_factor <- interaction(l, d, drop = TRUE)
    y_matrix <- y_matrix + generate_correlated_effects(nlevels(ld_factor), diag(n_dims), rep(sd_ld, n_dims))[ld_factor, ]

    colnames(y_matrix) <- paste0("y", seq_len(n_dims))

    result <- data.frame(p = p, i = i, d = d, l = l)
    cbind(result, as.data.frame(y_matrix))
  }
}
