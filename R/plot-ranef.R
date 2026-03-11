#' Caterpillar Plot for Random Effects
#'
#' Creates a caterpillar (dot plot with error bars) to display random effects
#' from a gstudy object. Works with both lme4 and brms backends, providing a
#' unified interface regardless of backend used.
#'
#' @param x A gstudy object.
#' @param which Which random effect(s) to plot. If NULL (default), plots all
#'   random effects in a multi-panel plot. Can be a character vector of
#'   specific facet names, or a single facet name.
#' @param ci_level Confidence/credible interval level. Default is 0.95.
#' @param sort Logical; if TRUE (default), sort effects by estimate size.
#' @param colors Optional vector of two colors for points and error bars.
#' @param ncol Number of columns for multi-panel plot when which is NULL.
#'   Default is NULL (automatic layout).
#' @param ... Additional arguments (passed to ggplot2 theme functions).
#' @return A ggplot object (invisibly).
#' @export
#' @examples
#' # Fit a G-study with lme4 (default) using the brennan dataset
#' g <- gstudy(Score ~ (1 | Person) + (1 | Task) + (1 | Rater),
#'   data = brennan
#' )
#'
#' # Plot all random effects
#' plot_ranef(g)
#'
#' # Plot only one specific random effect
#' plot_ranef(g, which = "Person")
#'
#' \donttest{
#' # Fit G-study with brms for Bayesian credible intervals
#' g_brms <- gstudy(Score ~ (1 | Person) + (1 | Task),
#'   data = brennan,
#'   backend = "brms"
#' )
#'
#' # Plot with brms (uses posterior credible intervals)
#' plot_ranef(g_brms)
#' }
plot_ranef <- function(x, which = NULL, ci_level = 0.95,
                       sort = TRUE, colors = NULL, ncol = NULL, ...) {
  UseMethod("plot_ranef")
}

#' @export
plot_ranef.gstudy <- function(x, which = NULL, ci_level = 0.95,
                                sort = TRUE, colors = NULL, ncol = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (!inherits(x, "gstudy")) {
    stop("x must be a gstudy object", call. = FALSE)
  }

  if (x$backend == "mom") {
    message("Caterpillar plots are not available for the method of moments backend.")
    return(invisible(NULL))
  }

  model <- x$model

  if (x$backend == "lme4") {
    ranef_data <- prepare_ranef_lme4(model, ci_level = ci_level)
  } else if (x$backend == "brms") {
    ranef_data <- prepare_ranef_brms(model, ci_level = ci_level)
  } else {
    stop("Unknown backend: ", x$backend, call. = FALSE)
  }

  if (is.null(ranef_data) || nrow(ranef_data) == 0) {
    message("No random effects to plot.")
    return(invisible(NULL))
  }

  available_facets <- unique(ranef_data$facet)

  if (!is.null(which)) {
    which <- match.arg(which, available_facets, several.ok = TRUE)
    ranef_data <- ranef_data[ranef_data$facet %in% which, ]
  }

  if (sort) {
    ranef_data <- ranef_data[order(ranef_data$estimate), ]
  }

  ranef_data$level <- factor(ranef_data$level, levels = unique(ranef_data$level))

  if (is.null(colors)) {
    colors <- c("steelblue", "darkblue")
  }

  p <- ggplot2::ggplot(ranef_data, ggplot2::aes(x = level, y = estimate)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = lower, ymax = upper),
      color = colors[1],
      size = 0.5,
      fatten = 2
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Caterpillar Plot of Random Effects",
      x = "Level",
      y = "Estimate"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_text(size = 8),
      ...
    )

  if (!is.null(which) || length(available_facets) == 1) {
    if (!is.null(which) && length(which) == 1) {
      p <- p + ggplot2::labs(subtitle = paste("Facet:", which))
    }
  } else {
    p <- p + ggplot2::facet_wrap(~facet, ncol = ncol, scales = "free_y")
  }

  print(p)
  invisible(p)
}

prepare_ranef_lme4 <- function(model, ci_level = 0.95) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required.", call. = FALSE)
  }

  rr <- lme4::ranef(model, condVar = TRUE)

  if (length(rr) == 0) {
    return(NULL)
  }

  df_all <- lapply(names(rr), function(facet_name) {
    re_mat <- rr[[facet_name]]
    if (is.null(re_mat)) {
      return(NULL)
    }

    post_var <- attr(re_mat, "postVar")
    if (is.null(post_var)) {
      stop("Conditional variances not available. Ensure condVar = TRUE.",
           call. = FALSE)
    }

    n_levels <- nrow(re_mat)
    n_params <- ncol(re_mat)

    results <- lapply(seq_len(n_params), function(j) {
      estimates <- re_mat[, j]
      ses <- sqrt(squeeze(post_var[, j, j], .Machine$double.eps))

      z <- qnorm(1 - (1 - ci_level) / 2)

      param_name <- colnames(re_mat)[j]

      data.frame(
        level = rownames(re_mat),
        estimate = estimates,
        lower = estimates - z * ses,
        upper = estimates + z * ses,
        facet = facet_name,
        parameter = param_name,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, results)
  })

  do.call(rbind, df_all)
}

prepare_ranef_brms <- function(model, ci_level = 0.95) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.", call. = FALSE)
  }

  rr <- brms::ranef(model, summary = FALSE)

  if (length(rr) == 0) {
    return(NULL)
  }

  quantiles <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)

  df_all <- lapply(names(rr), function(facet_name) {
    re_array <- rr[[facet_name]]

    if (is.null(re_array) || length(dim(re_array)) < 3) {
      return(NULL)
    }

    dim_info <- dim(re_array)
    n_levels <- dim_info[2]
    n_params <- dim_info[3]

    results <- lapply(seq_len(n_params), function(j) {
      level_names <- dimnames(re_array)[[2]]
      param_name <- dimnames(re_array)[[3]][j]

      estimate <- apply(re_array[, , j], 2, mean)
      lower <- apply(re_array[, , j], 2, quantile, probs = quantiles[1])
      upper <- apply(re_array[, , j], 2, quantile, probs = quantiles[2])

      data.frame(
        level = level_names,
        estimate = estimate,
        lower = lower,
        upper = upper,
        facet = facet_name,
        parameter = param_name,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, results)
  })

  do.call(rbind, df_all)
}

squeeze <- function(x, eps) {
  pmax(pmin(x, 1 - eps), eps)
}
