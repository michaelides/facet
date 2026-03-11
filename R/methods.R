#' S3 Methods for mgt Classes
#'
#' This file defines S3 methods for the gstudy and dstudy classes,
#' including print, summary, plot, and tidy methods.
#'
#' @name mgt-methods
NULL

#' Print Method for gstudy Objects
#'
#' @param x A gstudy object.
#' @param digits Number of digits to display.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
#' @rdname print.gstudy
print.gstudy <- function(x, digits = 3, scale = c("variance", "sd"), ...) {
scale <- match.arg(scale)

cat("Generalizability Study (G-Study)\n")
cat("================================\n\n")

cat("Backend:", x$backend, "\n")
cat("Formula:", deparse(x$formula), "\n")
cat("Number of observations:", x$n_obs, "\n")
cat("Multivariate:", if (x$is_multivariate) "Yes" else "No", "\n\n")

cat("Object of measurement:", x$object, "\n")
cat("Facets:", paste(x$facets, collapse = ", "), "\n")

  cat("\n")
  
  ssi <- x$sample_size_info
  ss_tibble <- x$sample_size_tibble
  
  if (!is.null(ss_tibble) && nrow(ss_tibble) > 0) {
    cat("Sample Sizes:\n")
    print(ss_tibble, row.names = FALSE)
    
    # Add nested footnote if exists
    if (!is.null(ssi$nested) && length(ssi$nested) > 0) {
      cat("\nNested details:\n")
      footnote_lines <- format_nested_footnote(ssi$nested)
      for (line in footnote_lines) {
        cat("  ", line, "\n", sep = "")
      }
    }
    cat("\n")
  } else {
    if (!is.null(ssi)) {
      cat("Sample Sizes:\n")
      
      if (length(ssi$main) > 0) {
        cat(" Main Effects:\n")
        for (i in seq_along(ssi$main)) {
          cat(sprintf("  %s: %.0f\n", names(ssi$main)[i], ssi$main[i]))
        }
      }
      
      if (length(ssi$interactions) > 0) {
        cat(" Interactions:\n")
        for (i in seq_along(ssi$interactions)) {
          cat(sprintf("  %s: %.0f\n", names(ssi$interactions)[i], ssi$interactions[i]))
        }
      }
      
      if (!is.null(ssi$residual) && ssi$residual$facets != "" && !is.na(ssi$residual$n)) {
        cat(" Residual:\n")
        cat(sprintf("  %s: %.0f\n", ssi$residual$facets, ssi$residual$n))
      }
      
      if (!is.null(ssi$nested) && length(ssi$nested) > 0) {
        cat(" Nested Effects:\n")
        nested_lines <- format_nested_info(ssi$nested, indent = 2)
        for (line in nested_lines) {
          cat(line, "\n")
        }
      }
      cat("\n")
    } else if (!is.null(x$facet_n) && length(x$facet_n) > 0) {
      cat("Facet Sample Sizes:\n")
      for (i in seq_along(x$facet_n)) {
        cat(sprintf(" %s: %.0f\n", names(x$facet_n)[i], x$facet_n[i]))
      }
      cat("\n")
    }
  }

cat("Variance Components:\n")
vc_summary <- summarize_vc(x$variance_components, digits = digits, scale = scale)
print(vc_summary, row.names = FALSE, ...)

  invisible(x)
}

print_correlations <- function(correlations, digits = 3, format = c("long", "matrix")) {
  format <- match.arg(format)

  if (format == "matrix") {
    if (!is.null(correlations$random_effect_cor_matrix) && length(correlations$random_effect_cor_matrix) > 0) {
      for (facet_name in names(correlations$random_effect_cor_matrix)) {
        cat(sprintf("Correlated Random Effects (%s):\n", facet_name))
        cor_mat <- round(correlations$random_effect_cor_matrix[[facet_name]], digits)
        print(cor_mat, quote = FALSE)
        cat("\n")
      }
    }

    if (!is.null(correlations$residual_cor_matrix)) {
      cat("Residual Correlations:\n")
      cor_mat <- round(correlations$residual_cor_matrix, digits)
      print(cor_mat, quote = FALSE)
      cat("\n")
    }
  } else {
    if (!is.null(correlations$random_effect_cor) && length(correlations$random_effect_cor) > 0) {
      for (facet_name in names(correlations$random_effect_cor)) {
        if (!is.null(correlations$random_effect_cor[[facet_name]]) &&
            nrow(correlations$random_effect_cor[[facet_name]]) > 0) {
          cat(sprintf("Correlated Random Effects (%s):\n", facet_name))
          cor_summary <- summarize_cor(correlations$random_effect_cor[[facet_name]], digits)
          print(cor_summary, row.names = FALSE)
          cat("\n")
        }
      }
    }

    if (!is.null(correlations$residual_cor) && nrow(correlations$residual_cor) > 0) {
      cat("Correlated Residuals:\n")
      cor_summary <- summarize_cor(correlations$residual_cor, digits)
      print(cor_summary, row.names = FALSE)
      cat("\n")
    }
  }
}

#' Print Method for mgstudy Objects
#'
#' @param x An mgstudy object.
#' @param digits Number of digits to display.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @param cor_format Format for displaying correlations: "long" (default) or "matrix".
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
#' @rdname print.mgstudy
print.mgstudy <- function(x, digits = 3, scale = c("variance", "sd"), cor_format = c("long", "matrix"), ...) {
  scale <- match.arg(scale)
  cor_format <- match.arg(cor_format)

  cat("Multivariate Generalizability Study (MG-Study)\n")
  cat("=============================================\n\n")

  cat("Backend:", x$backend, "\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Number of observations:", x$n_obs, "\n")
  cat("Dimensions:", paste(x$dimensions, collapse = ", "), "\n\n")

  cat("Object of measurement:", x$object, "\n")
  cat("Facets:", paste(x$facets, collapse = ", "), "\n")

  ssi <- x$sample_size_info
  ss_tibble <- x$sample_size_tibble
  
  if (!is.null(ss_tibble) && nrow(ss_tibble) > 0) {
    cat("\nSample Sizes:\n")
    print(ss_tibble, row.names = FALSE)
    
    # Add nested footnote if exists
    if (!is.null(ssi$nested) && length(ssi$nested) > 0) {
      cat("\nNested details:\n")
      footnote_lines <- format_nested_footnote(ssi$nested)
      for (line in footnote_lines) {
        cat("  ", line, "\n", sep = "")
      }
    }
    cat("\n")
  } else {
    if (!is.null(ssi)) {
      cat("\nSample Sizes:\n")
      
      if (length(ssi$main) > 0) {
        cat(" Main Effects:\n")
        for (i in seq_along(ssi$main)) {
          cat(sprintf("  %s: %.0f\n", names(ssi$main)[i], ssi$main[i]))
        }
      }
      cat("\n")
    } else if (!is.null(x$facet_n) && length(x$facet_n) > 0) {
      cat("Facet Sample Sizes:\n")
      for (i in seq_along(x$facet_n)) {
        cat(sprintf("  %s: %.0f\n", names(x$facet_n)[i], x$facet_n[i]))
      }
      cat("\n")
    }
  }

  cat("Variance Components (by dimension):\n")
  for (dim in x$dimensions) {
    cat(sprintf("\nDimension: %s\n", dim))
    vc_dim <- x$variance_components[x$variance_components$dim == dim, ]
    vc_summary <- summarize_vc(vc_dim, digits = digits, scale = scale)
    print(vc_summary, row.names = FALSE, ...)
  }

  if (!is.null(x$correlations)) {
    cat("\n")
    print_correlations(x$correlations, digits = digits, format = cor_format)
  }

  invisible(x)
}

#' Print Method for dstudy Objects
#'
#' @param x A dstudy object.
#' @param digits Number of digits to display.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @param sem Logical; if TRUE, include standard errors of measurement
#'   (sem_rel and sem_abs) in the output. Default is FALSE.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns x.
#' @export
#' @rdname print.dstudy
print.dstudy <- function(x, digits = 3, scale = c("variance", "sd"), sem = FALSE, ...) {
  scale <- match.arg(scale)

cat("Decision Study (D-Study)\n")
cat("========================\n\n")

cat("Based on G-Study with", x$gstudy$backend, "backend\n")
cat("Object of measurement:", x$object, "\n")

# Display universe components
if (!is.null(x$universe) && length(x$universe) > 0) {
cat("Universe components:", paste(x$universe, collapse = ", "), "\n")
}

error_components <- identify_error_components(
x$variance_components, x$universe, x$error,
x$residual_composition
)
  cat("Error components for relative error (sigma2_delta):",
    paste(error_components$relative, collapse = ", "), "\n")
  cat("Error components for absolute error (sigma2_delta_abs):",
    paste(error_components$absolute, collapse = ", "), "\n\n")

if (x$is_sweep) {
    has_dim <- "dim" %in% names(x$coefficients)
    has_estimate <- "estimate" %in% names(x$coefficients)

    if (has_dim) {
      dims <- unique(x$coefficients$dim)
      cat("Sample Size Sweep (by dimension):\n\n")
      for (d in dims) {
        cat(sprintf("Dimension: %s\n", d))
        cat(paste(rep("-", 40), collapse = ""), "\n")
        coefs <- x$coefficients[x$coefficients$dim == d, ]

        if (!sem && all(c("sem_rel", "sem_abs") %in% names(coefs))) {
          coefs <- coefs[, !names(coefs) %in% c("sem_rel", "sem_abs"), drop = FALSE]
        }
        print(coefs, digits = digits, row.names = FALSE, ...)
        cat("\n")
      }
    } else {
      cat("Sample Size Sweep:\n")
      coefs <- x$coefficients

      if (!sem && all(c("sem_rel", "sem_abs") %in% names(coefs))) {
        coefs <- coefs[, !names(coefs) %in% c("sem_rel", "sem_abs")]
      }
      print(coefs, digits = digits, ...)
    }
 } else {
 cat("Sample Sizes:\n")
 for (name in names(x$n)) {
 cat(sprintf(" %s: %.0f\n", name, x$n[[name]]))
 }
 cat("\n")

 cat("Variance Components:\n")
 vc_summary <- summarize_vc(x$variance_components, digits = digits, scale = scale)
 print(vc_summary, row.names = FALSE, ...)
 cat("\n")

 cat("Coefficients:\n")
 coefs <- x$coefficients

 if (!sem && all(c("sem_rel", "sem_abs") %in% names(coefs))) {
 coefs <- coefs[, !names(coefs) %in% c("sem_rel", "sem_abs")]
 }
 print(coefs, digits = digits, row.names = FALSE, ...)

 invisible(x)
 }
 }

identify_error_components <- function(vc, universe_spec, error_spec, residual_composition = NULL) {
universe_parsed <- parse_specification_internal(universe_spec)

if (!is.null(error_spec)) {
error_parsed <- parse_specification_internal(error_spec)
relative_components <- error_parsed
absolute_components <- error_parsed
} else {
relative_components <- vc$component[
sapply(vc$component, function(comp) {
is_rel_error_component(comp, universe_parsed)
})
]

absolute_components <- vc$component[!(vc$component %in% universe_parsed)]
}

  if (!is.null(residual_composition) && residual_composition != "") {
    relative_components <- sapply(relative_components, function(comp) {
      if (comp == "Residual") {
        paste0(residual_composition, " (Residual)")
      } else {
        comp
      }
    })
    absolute_components <- sapply(absolute_components, function(comp) {
      if (comp == "Residual") {
        paste0(residual_composition, " (Residual)")
      } else {
        comp
      }
    })
  }

  list(
    relative = relative_components,
    absolute = absolute_components
  )
}

parse_specification_internal <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  if (inherits(x, "formula")) {
    if (length(x) == 2) {
      rhs <- x[[2]]
    } else {
      rhs <- x[[3]]
    }

    components <- character(0)

    if (is.call(rhs)) {
      parse_formula_terms_internal <- function(expr) {
        result <- character(0)
        if (is.call(expr)) {
          if (expr[[1]] == as.symbol("+")) {
            result <- c(result, parse_formula_terms_internal(expr[[2]]))
            result <- c(result, parse_formula_terms_internal(expr[[3]]))
          } else if (expr[[1]] == as.symbol(":")) {
            result <- c(result, deparse(expr))
          } else {
            result <- c(result, parse_formula_terms_internal(expr[[2]]))
          }
        } else if (is.symbol(expr)) {
          result <- c(result, as.character(expr))
        }
        result
      }
      components <- parse_formula_terms_internal(rhs)
    } else if (is.symbol(rhs)) {
      components <- as.character(rhs)
    }

    components <- components[components != "."]
    return(components)
  }

  if (is.character(x)) {
    return(x)
  }

  character(0)
}

is_rel_error_component <- function(component, universe_spec) {
if (component %in% universe_spec) {
return(FALSE)
}

if (component == "Residual") {
return(TRUE)
}

facets <- strsplit(component, ":")[[1]]

any(universe_spec %in% facets)
}

#' Summary Method for gstudy Objects
#'
#' @param object A gstudy object.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @param digits Number of digits to display (default 3).
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns object.
#' @export
#' @rdname summary.gstudy
summary.gstudy <- function(object, scale = c("variance", "sd"), digits = 3, ...) {
scale <- match.arg(scale)

cat("=== G-Study Summary ===\n\n")

cat("Design Information:\n")
cat(" Backend:", object$backend, "\n")
cat(" Formula:", deparse(object$formula), "\n")
cat(" Observations:", object$n_obs, "\n")
cat(" Multivariate:", if (object$is_multivariate) "Yes" else "No", "\n\n")

cat("Facet Information:\n")
cat(" Object of measurement:", object$object, "\n")
cat(" Facets:", paste(object$facets, collapse = ", "), "\n")

  ssi <- object$sample_size_info
  ss_tibble <- object$sample_size_tibble
  
  if (!is.null(ss_tibble) && nrow(ss_tibble) > 0) {
    cat("\nSample Sizes:\n")
    print(ss_tibble, row.names = FALSE)
    
    # Add nested footnote if exists
    if (!is.null(ssi$nested) && length(ssi$nested) > 0) {
      cat("\nNested details:\n")
      footnote_lines <- format_nested_footnote(ssi$nested)
      for (line in footnote_lines) {
        cat("  ", line, "\n", sep = "")
      }
    }
    cat("\n")
  } else {
    if (!is.null(ssi)) {
      cat("\nSample Sizes:\n")
      
      if (length(ssi$main) > 0) {
        cat(" Main Effects:\n")
        for (i in seq_along(ssi$main)) {
          cat(sprintf("  %s: %.0f\n", names(ssi$main)[i], ssi$main[i]))
        }
      }
      
      if (length(ssi$interactions) > 0) {
        cat(" Interactions:\n")
        for (i in seq_along(ssi$interactions)) {
          cat(sprintf("  %s: %.0f\n", names(ssi$interactions)[i], ssi$interactions[i]))
        }
      }
      
      if (!is.null(ssi$residual) && ssi$residual$facets != "" && !is.na(ssi$residual$n)) {
        cat(" Residual:\n")
        cat(sprintf("  %s: %.0f\n", ssi$residual$facets, ssi$residual$n))
      }
      
      if (!is.null(ssi$nested) && length(ssi$nested) > 0) {
        cat(" Nested Effects:\n")
        nested_lines <- format_nested_info(ssi$nested, indent = 2)
        for (line in nested_lines) {
          cat(line, "\n")
        }
      }
      cat("\n")
    } else if (!is.null(object$facet_n) && length(object$facet_n) > 0) {
      cat(" Facet Sample Sizes:\n")
      for (i in seq_along(object$facet_n)) {
        cat(sprintf("  %s: %.0f\n", names(object$facet_n)[i], object$facet_n[i]))
      }
      cat("\n")
    }
  }

cat("Variance Components:\n")
vc_summary <- summarize_vc(object$variance_components, digits = digits, scale = scale)
print(vc_summary, row.names = FALSE, ...)

cat("\nTotal variance:", sum(object$variance_components$var, na.rm = TRUE), "\n")

  invisible(object)
}

#' Summary Method for mgstudy Objects
#'
#' @param object An mgstudy object.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @param digits Number of digits to display (default 3).
#' @param cor_format Format for displaying correlations: "long" (default) or "matrix".
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns object.
#' @export
#' @rdname summary.mgstudy
summary.mgstudy <- function(object, scale = c("variance", "sd"), digits = 3, cor_format = c("long", "matrix"), ...) {
  scale <- match.arg(scale)
  cor_format <- match.arg(cor_format)

  cat("=== Multivariate G-Study Summary ===\n\n")

  cat("Design Information:\n")
  cat(" Backend:", object$backend, "\n")
  cat(" Formula:", deparse(object$formula), "\n")
  cat(" Observations:", object$n_obs, "\n")
  cat(" Dimensions:", paste(object$dimensions, collapse = ", "), "\n\n")

  cat("Facet Information:\n")
  cat(" Object of measurement:", object$object, "\n")
  cat(" Facets:", paste(object$facets, collapse = ", "), "\n")

  ssi <- object$sample_size_info
  ss_tibble <- object$sample_size_tibble
  
  if (!is.null(ss_tibble) && nrow(ss_tibble) > 0) {
    cat("\nSample Sizes:\n")
    print(ss_tibble, row.names = FALSE)
    
    # Add nested footnote if exists
    if (!is.null(ssi$nested) && length(ssi$nested) > 0) {
      cat("\nNested details:\n")
      footnote_lines <- format_nested_footnote(ssi$nested)
      for (line in footnote_lines) {
        cat("  ", line, "\n", sep = "")
      }
    }
    cat("\n")
  } else {
    if (!is.null(ssi)) {
      cat("\nSample Sizes:\n")
      
      if (length(ssi$main) > 0) {
        cat(" Main Effects:\n")
        for (i in seq_along(ssi$main)) {
          cat(sprintf("  %s: %.0f\n", names(ssi$main)[i], ssi$main[i]))
        }
      }
      cat("\n")
    } else if (!is.null(object$facet_n) && length(object$facet_n) > 0) {
      cat(" Facet Sample Sizes:\n")
      for (i in seq_along(object$facet_n)) {
        cat(sprintf("  %s: %.0f\n", names(object$facet_n)[i], object$facet_n[i]))
      }
      cat("\n")
    }
  }

  cat("Variance Components (by dimension):\n")
  for (dim in object$dimensions) {
    cat(sprintf("\nDimension: %s\n", dim))
    vc_dim <- object$variance_components[object$variance_components$dim == dim, ]
    vc_summary <- summarize_vc(vc_dim, digits = digits, scale = scale)
    print(vc_summary, row.names = FALSE, ...)
  }

  if (!is.null(object$correlations)) {
    cat("\n")
    print_correlations(object$correlations, digits = digits, format = cor_format)
  }

  invisible(object)
}

#' Summary Method for dstudy Objects
#'
#' @param object A dstudy object.
#' @param scale Scale for displaying results: "variance" (default) or "sd".
#' @param digits Number of digits to display.
#' @param sem Logical; if TRUE, include standard errors of measurement
#'   (sem_rel and sem_abs) in the output. Default is FALSE.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns object.
#' @export
#' @rdname summary.dstudy
summary.dstudy <- function(object, scale = c("variance", "sd"), digits = 3, sem = FALSE, ...) {
  scale <- match.arg(scale)

  cat("=== D-Study Summary ===\n\n")

  cat("Object of measurement:", object$object, "\n")
  cat("Based on G-Study:", deparse(object$gstudy$formula), "\n\n")

  if (object$is_sweep) {
    has_dim <- "dim" %in% names(object$coefficients)
    if (has_dim) {
      dims <- unique(object$coefficients$dim)
      cat("Sample Size Sweep Results (by dimension):\n\n")
      for (d in dims) {
        cat(sprintf("Dimension: %s\n", d))
        cat(paste(rep("-", 40), collapse = ""), "\n")
        coefs <- object$coefficients[object$coefficients$dim == d, ]
        coefs <- coefs[, names(coefs) != "dim", drop = FALSE]
        if (!sem && all(c("sem_rel", "sem_abs") %in% names(coefs))) {
          coefs <- coefs[, !names(coefs) %in% c("sem_rel", "sem_abs")]
        }
        print(coefs, row.names = FALSE, ...)
        
        if ("g" %in% names(coefs)) {
          best_g <- coefs[which.max(coefs$g), ]
          cat("\nHighest G coefficient for", d, ":\n")
          print(best_g, row.names = FALSE)
        }
        cat("\n")
      }
    } else {
      cat("Sample Size Sweep Results:\n")
      coefs <- object$coefficients
      if (!sem && all(c("sem_rel", "sem_abs") %in% names(coefs))) {
        coefs <- coefs[, !names(coefs) %in% c("sem_rel", "sem_abs")]
      }
      print(coefs, ...)

      # Find optimal configurations
      if ("g" %in% names(object$coefficients)) {
        best_g <- object$coefficients[which.max(object$coefficients$g), ]
        cat("\nHighest G coefficient:\n")
        print(best_g)
      }
    }
  } else {
    cat("Sample Sizes:\n")
 for (name in names(object$n)) {
 cat(sprintf(" %s: %.0f\n", name, object$n[[name]]))
 }
 cat("\n")

 cat("Variance Components:\n")
 vc_summary <- summarize_vc(object$variance_components, digits = digits, scale = scale)
 print(vc_summary, row.names = FALSE, ...)
 cat("\n")

 cat("Coefficients:\n")
 coefs <- object$coefficients
 if (!sem && all(c("sem_rel", "sem_abs") %in% names(coefs))) {
 coefs <- coefs[, !names(coefs) %in% c("sem_rel", "sem_abs")]
 }
 print(coefs, ...)
 }

 invisible(object)
}

#' Plot Method for gstudy Objects
#'
#' @param x A gstudy object.
#' @param type Type of plot: "variance", "proportion", or "forest".
#' @param ... Additional arguments passed to plotting functions.
#' @return A ggplot object (invisibly).
#' @export
#' @rdname plot.gstudy
plot.gstudy <- function(x, type = c("variance", "proportion", "forest"), ...) {
  type <- match.arg(type)

  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  vc <- x$variance_components

  p <- switch(type,
    "variance" = {
      ggplot2::ggplot(vc, ggplot2::aes(x = reorder(component, -var), y = var)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Variance Components",
          x = "Component",
          y = "Variance Estimate"
        ) +
        ggplot2::theme_minimal()
    },
    "proportion" = {
      ggplot2::ggplot(vc, ggplot2::aes(x = reorder(component, -pct), y = pct)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Proportion of Total Variance",
          x = "Component",
          y = "Percentage (%)"
        ) +
        ggplot2::theme_minimal()
    },
    "forest" = {
      # Forest plot with confidence intervals (if available)
      if (!all(c("lower", "upper") %in% names(vc))) {
        warning("Confidence intervals not available. Using 'variance' plot instead.",
          call. = FALSE
        )
        return(plot.gstudy(x, type = "variance", ...))
      }
      ggplot2::ggplot(vc, ggplot2::aes(x = reorder(component, -var), y = var)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Variance Components with Confidence Intervals",
          x = "Component",
          y = "Variance Estimate"
        ) +
        ggplot2::theme_minimal()
    }
  )

  print(p)
  invisible(p)
}

#' Plot Method for mgstudy Objects
#'
#' @param x An mgstudy object.
#' @param type Type of plot: "variance", "proportion", or "forest".
#' @param ... Additional arguments passed to plotting functions.
#' @return A ggplot object (invisibly).
#' @export
#' @rdname plot.mgstudy
plot.mgstudy <- function(x, type = c("variance", "proportion", "forest"), ...) {
  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  vc <- x$variance_components
  dimensions <- x$dimensions

  p <- switch(type,
    "variance" = {
      ggplot2::ggplot(vc, ggplot2::aes(x = reorder(component, -var), y = var)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Variance Components by Dimension",
          x = "Component",
          y = "Variance Estimate"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::facet_wrap(~dim, scales = "free_x")
    },
    "proportion" = {
      ggplot2::ggplot(vc, ggplot2::aes(x = reorder(component, -pct), y = pct)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Proportion of Total Variance by Dimension",
          x = "Component",
          y = "Percentage (%)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::facet_wrap(~dim)
    },
    "forest" = {
      if (!all(c("lower", "upper") %in% names(vc))) {
        warning("Confidence intervals not available. Using 'variance' plot instead.",
          call. = FALSE
        )
        return(plot.mgstudy(x, type = "variance", ...))
      }
      ggplot2::ggplot(vc, ggplot2::aes(x = reorder(component, -var), y = var)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.2) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Variance Components with Confidence Intervals by Dimension",
          x = "Component",
          y = "Variance Estimate"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::facet_wrap(~dim, scales = "free_x")
    }
  )

  print(p)
  invisible(p)
}

#' Plot Method for dstudy Objects
#'
#' @param x A dstudy object.
#' @param type Type of plot: "coefficients" or "sweep".
#' @param coefficient Which coefficient to plot: "g", "phi", or "both".
#' @param ... Additional arguments passed to plotting functions.
#' @return A ggplot object (invisibly).
#' @export
#' @rdname plot.dstudy
plot.dstudy <- function(x, type = c("coefficients", "sweep"), coefficient = c("both", "g", "phi"), ...) {
  type <- match.arg(type)
  coefficient <- match.arg(coefficient)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.", call. = FALSE)
  }

  if (type == "sweep") {
    if (!x$is_sweep) {
      warning(
        "Plot is not defined for individual g and phi coefficients. ",
        "n must be a vector (or list of vectors) to create a sweep plot.",
        call. = FALSE
      )
      return(invisible(NULL))
    }

    coefs <- x$coefficients
    facet_names <- names(x$n)

    if (length(facet_names) == 0) {
      stop("No facet information available for sweep plot.", call. = FALSE)
    }

    sweeping_facets <- facet_names[sapply(x$n, length) > 1]

    if (length(sweeping_facets) == 0) {
      warning(
        "No sweeping facets found. n contains only single values.",
        call. = FALSE
      )
      return(invisible(NULL))
    }

    x_var <- sweeping_facets[1]
    group_var <- if (length(sweeping_facets) > 1) sweeping_facets[2] else NULL
    
    is_multivariate <- x$is_multivariate && "dim" %in% names(coefs)

    if (coefficient == "both") {
      coef_long <- tidyr::pivot_longer(
        coefs,
        cols = c("g", "phi"),
        names_to = "coefficient_type",
        values_to = "value"
      )
      coef_long$coefficient_type <- factor(
        coef_long$coefficient_type,
        levels = c("g", "phi"),
        labels = c("G (Relative)", "Phi (Absolute)")
      )

      if (!is.null(group_var)) {
        p <- ggplot2::ggplot(coef_long, ggplot2::aes(x = .data[[x_var]], y = value, color = .data[[group_var]], group = .data[[group_var]])) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::geom_point(size = 2) +
          ggplot2::ylim(0, 1) +
          ggplot2::labs(
            title = "Coefficients Across Sample Sizes",
            x = x_var,
            y = "Coefficient Value",
            color = group_var
          ) +
          ggplot2::theme_minimal()
      } else {
        p <- ggplot2::ggplot(coef_long, ggplot2::aes(x = .data[[x_var]], y = value, color = coefficient_type)) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::geom_point(size = 2) +
          ggplot2::ylim(0, 1) +
          ggplot2::labs(
            title = "Coefficients Across Sample Sizes",
            x = x_var,
            y = "Coefficient Value",
            color = "Coefficient"
          ) +
          ggplot2::theme_minimal()
      }
      
      if (is_multivariate) {
        p <- p + ggplot2::facet_grid(coefficient_type ~ dim)
      } else {
        p <- p + ggplot2::facet_wrap(~coefficient_type)
      }
    } else {
      y_var <- coefficient
      y_label <- if (coefficient == "g") "G Coefficient" else "Phi Coefficient"

      if (!is.null(group_var)) {
        p <- ggplot2::ggplot(coefs, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[group_var]], group = .data[[group_var]])) +
          ggplot2::geom_line(linewidth = 1) +
          ggplot2::geom_point(size = 2) +
          ggplot2::ylim(0, 1) +
          ggplot2::labs(
            title = paste(y_label, "Across Sample Sizes"),
            x = x_var,
            y = y_label,
            color = group_var
          ) +
          ggplot2::theme_minimal()
      } else {
        p <- ggplot2::ggplot(coefs, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
          ggplot2::geom_line(linewidth = 1, color = "steelblue") +
          ggplot2::geom_point(size = 2, color = "steelblue") +
          ggplot2::ylim(0, 1) +
          ggplot2::labs(
            title = paste(y_label, "Across Sample Sizes"),
            x = x_var,
            y = y_label
          ) +
          ggplot2::theme_minimal()
      }
      
      if (is_multivariate) {
        p <- p + ggplot2::facet_wrap(~dim)
      }
    }

    print(p)
    invisible(p)
  } else {
    coefs <- x$coefficients
    if (!"g" %in% names(coefs)) {
      stop("Coefficients not available for plotting.", call. = FALSE)
    }

    coef_df <- data.frame(
      coefficient = c("G (Relative)", "Phi (Absolute)"),
      value = c(coefs$g[1], coefs$phi[1])
    )

    p <- ggplot2::ggplot(coef_df, ggplot2::aes(x = coefficient, y = value)) +
      ggplot2::geom_bar(stat = "identity", fill = c("steelblue", "coral")) +
      ggplot2::ylim(0, 1) +
      ggplot2::labs(
        title = "Generalizability and Dependability Coefficients",
        x = "Coefficient Type",
        y = "Value"
      ) +
      ggplot2::theme_minimal()

    print(p)
    invisible(p)
  }
}

#' Tidy Method for gstudy Objects
#'
#' Returns a tidy tibble of variance components.
#'
#' @param x A gstudy object.
#' @param ... Additional arguments (ignored).
#' @return A tibble.
#' @export
#' @rdname tidy.gstudy
tidy.gstudy <- function(x, ...) {
  x$variance_components
}

#' Tidy Method for dstudy Objects
#'
#' Returns a tidy tibble of coefficients and variance components.
#'
#' @param x A dstudy object.
#' @param ... Additional arguments (ignored).
#' @return A tibble.
#' @export
#' @rdname tidy.dstudy
tidy.dstudy <- function(x, ...) {
  x$coefficients
}

#' Tidy Method for mgstudy Objects
#'
#' Returns a tidy tibble of variance components for multivariate G-studies.
#'
#' @param x An mgstudy object.
#' @param ... Additional arguments (ignored).
#' @return A tibble.
#' @export
#' @rdname tidy.gstudy
tidy.mgstudy <- function(x, ...) {
  x$variance_components
}

#' Glance Method for gstudy Objects
#'
#' Returns a one-row summary of model statistics.
#'
#' @param x A gstudy object.
#' @param ... Additional arguments (ignored).
#' @return A one-row tibble.
#' @export
#' @rdname glance.gstudy
glance.gstudy <- function(x, ...) {
  tibble::tibble(
    n_obs = x$n_obs,
    n_facets = length(x$facets),
    backend = x$backend,
    is_multivariate = x$is_multivariate,
    object = x$object
  )
}

#' Glance Method for dstudy Objects
#'
#' Returns a one-row summary of D-study results.
#'
#' @param x A dstudy object.
#' @param ... Additional arguments (ignored).
#' @return A one-row tibble.
#' @export
#' @rdname glance.dstudy
glance.dstudy <- function(x, ...) {
  if (x$is_sweep) {
    tibble::tibble(
      is_sweep = TRUE,
      n_combinations = nrow(x$coefficients),
      object = x$object
    )
  } else {
    tibble::tibble(
      is_sweep = FALSE,
      g = x$coefficients$g[1],
      phi = x$coefficients$phi[1],
      object = x$object
    )
  }
}

#' Glance Method for mgstudy Objects
#'
#' Returns a one-row summary of model statistics for multivariate G-studies.
#'
#' @param x An mgstudy object.
#' @param ... Additional arguments (ignored).
#' @return A one-row tibble.
#' @export
#' @rdname glance.gstudy
glance.mgstudy <- function(x, ...) {
  tibble::tibble(
    n_obs = x$n_obs,
    n_facets = length(x$facets),
    n_dimensions = length(x$dimensions),
    dimensions = list(x$dimensions),
    backend = x$backend,
    is_multivariate = x$is_multivariate,
    object = x$object
  )
}

#' Variance-Covariance Matrix Extractor Generic
#'
#' Generic function for extracting variance-covariance matrices from fitted models.
#' Dispatches to the appropriate method based on the class of the object.
#'
#' @param x An object from which to extract variance components.
#' @param ... Additional arguments passed to methods.
#' @return Variance-covariance matrices or posterior summaries, depending on the method.
#' @export
VarCorr <- function(x, ...) {
  UseMethod("VarCorr")
}

#' Extract Variance-Covariance Matrices from gstudy Objects
#'
#' Extracts the variance-covariance matrices of random effects from a gstudy object.
#' This is a wrapper that calls the appropriate VarCorr method based on the backend
#' used to fit the model.
#'
#' @param x A gstudy object.
#' @param ... Additional arguments passed to the underlying VarCorr method
#' (e.g., \code{\link[lme4:VarCorr]{lme4::VarCorr}} or
#' \code{\link[brms:VarCorr]{brms::VarCorr}}).
#' @return For lme4 backend, returns a list of variance-covariance matrices
#' for each random effect term. For brms backend, returns a list with
#' posterior summaries.
#' @method VarCorr gstudy
#' @export
VarCorr.gstudy <- function(x, ...) {
  if (!inherits(x, "gstudy") && !inherits(x, "mgstudy")) {
    stop("x must be a gstudy object", call. = FALSE)
  }

  model <- x$model

  if (x$backend == "lme4") {
    lme4::VarCorr(model, ...)
  } else if (x$backend == "brms") {
    if (x$is_multivariate) {
      resp_names <- NULL
      if (!is.null(model$formula) && inherits(model$formula, "mvbrmsformula")) {
        resp_names <- model$formula$responses
      }
      if (is.null(resp_names) && !is.null(model$formula) && !is.null(model$formula$formula)) {
        bform <- model$formula$formula
        if (inherits(bform, "list") && length(bform) > 0) {
          resp_names <- names(bform)
        }
      }
      if (is.null(resp_names) || length(resp_names) <= 1) {
        brms::VarCorr(model, ...)
      } else {
        result <- list()
        for (resp in resp_names) {
          result[[resp]] <- brms::VarCorr(model, resp = resp, ...)
        }
        result
      }
    } else {
      brms::VarCorr(model, ...)
    }
  } else if (x$backend == "mom") {
    # VarCorr for method of moments backend
    model <- x$model
    
    # Extract variance components
    vc <- model$variance_components
    
    # Check if multivariate
    if (x$is_multivariate) {
      # Multivariate method of moments
      dimensions <- x$dimensions
      
      result <- list()
      
      for (dim in dimensions) {
        # Filter VC for this dimension
        vc_dim <- vc[vc$dim == dim, ]
        
        # Build data frame with Group, Std.Dev., Var.
        groups <- character()
        std_dev <- numeric()
        var_est <- numeric()
        
        # Add random effects
        random_effects <- vc_dim[vc_dim$component != "Residual", ]
        for (i in seq_len(nrow(random_effects))) {
          groups <- c(groups, random_effects$component[i])
          var_val <- max(0, random_effects$var[i])
          var_est <- c(var_est, var_val)
          std_dev <- c(std_dev, sqrt(var_val))
        }
        
        # Add residual
        resid_vc <- vc_dim[vc_dim$component == "Residual", ]
        if (nrow(resid_vc) > 0) {
          groups <- c(groups, "Residual")
          resid_var <- max(0, resid_vc$var[1])
          var_est <- c(var_est, resid_var)
          std_dev <- c(std_dev, sqrt(resid_var))
        }
        
        result[[dim]] <- data.frame(
          Group = groups,
          "Std.Dev." = std_dev,
          "Var." = var_est,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
      
      # Add correlation info if available
      if (!is.null(model$correlations)) {
        result$random_effect_cor <- model$correlations$random_effect_cor
        result$residual_cor <- model$correlations$residual_cor
      }
      
      return(result)
    } else {
      # Univariate method of moments
      # Build data frame with Group, Std.Dev., Var.
      groups <- character()
      std_dev <- numeric()
      var_est <- numeric()
      
      # Add random effects
      random_effects <- vc[vc$component != "Residual", ]
      for (i in seq_len(nrow(random_effects))) {
        groups <- c(groups, random_effects$component[i])
        var_val <- max(0, random_effects$var[i])
        var_est <- c(var_est, var_val)
        std_dev <- c(std_dev, sqrt(var_val))
      }
      
      # Add residual variance
      resid_vc <- vc[vc$component == "Residual", ]
      if (nrow(resid_vc) > 0) {
        groups <- c(groups, "Residual")
        resid_var <- max(0, resid_vc$var[1])
        var_est <- c(var_est, resid_var)
        std_dev <- c(std_dev, sqrt(resid_var))
      }
      
      result <- data.frame(
        Group = groups,
        "Std.Dev." = std_dev,
        "Var." = var_est,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      return(result)
    }
  } else {
    stop("Unknown backend: ", x$backend, call. = FALSE)
  }
}

#' @method VarCorr mgstudy
#' @export
VarCorr.mgstudy <- VarCorr.gstudy

#' Random Effects Extractor Generic
#'
#' Generic function for extracting random effects (BLUPs/conditional modes) from fitted models.
#' Dispatches to the appropriate method based on the class of the object.
#'
#' @param object An object from which to extract random effects.
#' @param ... Additional arguments passed to methods.
#' @return Random effects from the fitted model.
#' @export
ranef <- function(object, ...) {
UseMethod("ranef")
}

#' Extract Random Effects from gstudy Objects
#'
#' Extracts the random effects (BLUPs/conditional modes) from a gstudy object.
#' This is a wrapper that calls the appropriate ranef method based on the backend
#' used to fit the model.
#'
#' @param object A gstudy object.
#' @param ... Additional arguments passed to the underlying ranef method
#' (e.g., \code{\link[lme4:ranef.merMod]{lme4::ranef}} or
#' \code{\link[brms:ranef]{brms::ranef}}).
#' @return The random effects from the fitted model. For lme4 backend, returns
#' a list of matrices. For brms backend, returns a list of arrays with
#' posterior summaries.
#' @method ranef gstudy
#' @export
ranef.gstudy <- function(object, ...) {
if (!inherits(object, "gstudy") && !inherits(object, "mgstudy")) {
stop("object must be a gstudy object", call. = FALSE)
}

  model <- object$model

  if (object$backend == "lme4") {
    lme4::ranef(model, ...)
  } else if (object$backend == "brms") {
    brms::ranef(model, ...)
  } else if (object$backend == "mom") {
    # Empirical BLUPs for method of moments
    model <- object$model
    data <- object$data
    
    # Check if multivariate
    if (object$is_multivariate) {
      # Multivariate method of moments - empirical BLUPs
      responses <- model$responses
      random_facet_specs <- model$random_facet_specs
      
      result <- list()
      
      for (resp in responses) {
        # Compute grand mean for this response
        grand_mean <- mean(data[[resp]], na.rm = TRUE)
        
        resp_result <- list()
        
        for (facet_spec in random_facet_specs) {
          # Check if this is an interaction term (contains ":")
          if (grepl(":", facet_spec)) {
            # Interaction term: create combined factor
            # e.g., "rater:item" -> paste rater and item values
            facet_parts <- strsplit(facet_spec, ":")[[1]]
            facet_parts <- trimws(facet_parts)
            
            # Check if all parts exist in data
            if (!all(facet_parts %in% names(data))) next
            
            # Create combined factor
            combined_factor <- interaction(data[, facet_parts, drop = FALSE], drop = TRUE)
            
            # Compute means by combined factor level
            resp_data <- data[[resp]]
            facet_means <- aggregate(
              resp_data,
              by = list(facet = combined_factor),
              FUN = mean, na.rm = TRUE
            )
          } else {
            # Simple facet
            if (!facet_spec %in% names(data)) next
            
            # Compute means by facet level
            facet_means <- aggregate(
              data[[resp]],
              by = list(facet = data[[facet_spec]]),
              FUN = mean, na.rm = TRUE
            )
          }
          
          # Empirical BLUPs = level mean - grand mean
          eblup <- setNames(
            facet_means$x - grand_mean,
            facet_means$facet
          )
          
          resp_result[[facet_spec]] <- eblup
        }
        
        result[[resp]] <- resp_result
      }
      
      return(result)
    } else {
      # Univariate method of moments
      response <- model$response
      random_facet_specs <- model$random_facet_specs
      
      # Compute grand mean
      grand_mean <- mean(data[[response]], na.rm = TRUE)
      
      result <- list()
      
      for (facet_spec in random_facet_specs) {
        # Check if this is an interaction term (contains ":")
        if (grepl(":", facet_spec)) {
          # Interaction term: create combined factor
          facet_parts <- strsplit(facet_spec, ":")[[1]]
          facet_parts <- trimws(facet_parts)
          
          # Check if all parts exist in data
          if (!all(facet_parts %in% names(data))) next
          
          # Create combined factor
          combined_factor <- interaction(data[, facet_parts, drop = TRUE], drop = TRUE)
          
          # Compute means by combined factor level
          facet_means <- aggregate(
            data[[response]],
            by = list(facet = combined_factor),
            FUN = mean, na.rm = TRUE
          )
        } else {
          # Simple facet
          if (!facet_spec %in% names(data)) next
          
          # Compute means by facet level
          facet_means <- aggregate(
            data[[response]],
            by = list(facet = data[[facet_spec]]),
            FUN = mean, na.rm = TRUE
          )
        }
        
        # Empirical BLUPs = level mean - grand mean
        eblup <- setNames(
          facet_means$x - grand_mean,
          facet_means$facet
        )
        
        result[[facet_spec]] <- eblup
      }
      
      return(result)
    }
  } else {
    stop("Unknown backend: ", object$backend, call. = FALSE)
  }
}

#' @method ranef mgstudy
#' @export
ranef.mgstudy <- ranef.gstudy

#' Pairs Method for gstudy Objects
#'
#' Wrapper for pairs.brmsfit to visualize parameter pairs from a gstudy
#' object that was fit with the brms backend.
#'
#' @param x A gstudy object.
#' @param ... Additional arguments passed to pairs.brmsfit.
#' @return The result of pairs.brmsfit.
#' @export
#' @rdname pairs.gstudy
pairs.gstudy <- function(x, ...) {
if (x$backend != "brms") {
    stop("pairs.gstudy only works with gstudy objects fit with the brms backend", call. = FALSE)
  }

  if (!inherits(x$model, "brmsfit")) {
    stop("Model in gstudy object is not a brmsfit object", call. = FALSE)
  }

  pairs(x$model, ...)
}

#' Extract Draws from gstudy Objects
#'
#' Extract posterior draws from a gstudy object that was fit with the brms backend.
#' This is a wrapper around brms::as_draws_matrix.
#'
#' @param object A gstudy object.
#' @param ... Additional arguments passed to brms::as_draws_matrix.
#' @return A posterior draws_matrix object containing posterior draws for all model parameters.
#' @export
#' @rdname extract_draws.gstudy
extract_draws.gstudy <- function(object, ...) {
  if (!inherits(object, "gstudy")) {
    stop("object must be a gstudy object", call. = FALSE)
  }

  if (object$backend != "brms") {
    stop("extract_draws.gstudy only works with gstudy objects fit with the brms backend", call. = FALSE)
  }

  if (!inherits(object$model, "brmsfit")) {
    stop("Model in gstudy object is not a brmsfit object", call. = FALSE)
  }

  brms::as_draws_matrix(object$model, ...)
}
