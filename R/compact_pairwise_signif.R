#' Make and Filter Pairwise Statistical Comparisons
#'
#' Builds the observed factor cells, generates the biologically relevant
#' pairwise contrasts, performs a single family of pairwise tests, applies the
#' requested multiple-comparison correction, and filters on the adjusted
#' p-value threshold. This function performs no bracket positioning and does
#' not require ggplot2 or ggsignif.
#'
#' @inheritParams compact_pairwise_signif
#'
#' @return An object of class `pairwise_comparison_result`, containing:
#'   \describe{
#'     \item{comparisons}{A data frame of comparisons passing `alpha`.}
#'     \item{all_comparisons}{Every biologically retained comparison, including
#'       raw p-values, adjusted p-values, and test errors.}
#'     \item{groups}{Observed factor cells with plotting coordinates, sample
#'       sizes, and response maxima.}
#'     \item{normality}{Cell-level Shapiro-Wilk results used by `test = "auto"`.}
#'     \item{plot_data}{A copy of `data` with standardized plotting columns.}
#'     \item{observed_y}{The finite response values used for testing and layout.}
#'     \item{observed_data}{The response values together with their facet
#'       membership, when `facet` is supplied.}
#'     \item{method}{The test family actually used.}
#'     \item{settings}{The selected x layout and dodge width.}
#'   }
#'
#' @seealso [position_pairwise_significance()], [compact_pairwise_signif()]
#' @export
#'
#' @examples
#' staged_data <- transform(ToothGrowth, dose = factor(dose))
#' staged_tests <- make_pairwise_comparisons(
#'   staged_data,
#'   response = "len",
#'   factor1 = "dose",
#'   test = "t"
#' )
#' staged_bars <- position_pairwise_significance(
#'   staged_tests,
#'   annotation = "stars",
#'   panel_height_mm = 85,
#'   y_min = 0
#' )
#' staged_bars[, c("f1_1", "f1_2", "padj", "y_position")]
#' @md
make_pairwise_comparisons <- function(
    data,
    response,
    factor1,
    factor2 = NULL,
    facet = NULL,
    test = c("t", "wilcox", "auto"),
    test_args = list(),
    normality_alpha = 0.05,
    min_n = 2L,
    keep = NULL,
    p_adjust_method = "fdr",
    alpha = 0.05,
    x_layout = c("dodge", "interaction"),
    facet_scales = c("free_y", "fixed", "free", "free_x"),
    dodge_width = 0.75,
    interaction_sep = "_") {

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.null(facet) &&
      (!is.character(facet) || length(facet) != 1L || is.na(facet))) {
    stop("`facet` must be NULL or one column name.")
  }
  nm <- c(response, factor1, factor2, facet)
  nm <- nm[!is.na(nm)]
  if (!all(nm %in% names(data))) {
    stop("Unknown column(s): ", paste(setdiff(nm, names(data)), collapse = ", "))
  }
  if (!is.numeric(data[[response]])) stop("`response` must name a numeric column.")
  if (!is.null(factor2) && identical(factor1, factor2)) {
    stop("`factor1` and `factor2` must be different columns.")
  }
  if (!is.list(test_args)) stop("`test_args` must be a list.")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a number between 0 and 1.")
  }
  if (length(normality_alpha) != 1L || !is.finite(normality_alpha) ||
      normality_alpha < 0 || normality_alpha > 1) {
    stop("`normality_alpha` must be a number between 0 and 1.")
  }
  if (length(min_n) != 1L || !is.finite(min_n) || min_n < 1) {
    stop("`min_n` must be a positive number.")
  }

  x_layout <- match.arg(x_layout)
  facet_scales <- match.arg(facet_scales)
  if (is.null(factor2)) x_layout <- "interaction"

  if (!is.null(facet)) {
    facet_data <- data[!is.na(data[[facet]]), , drop = FALSE]
    if (nrow(facet_data) == 0L) stop("No non-missing facet values remain.")
    facet_values <- unique(as.character(facet_data[[facet]]))

    facet_results <- lapply(facet_values, function(facet_value) {
      panel_data <- facet_data[
        as.character(facet_data[[facet]]) == facet_value, , drop = FALSE
      ]
      make_pairwise_comparisons(
        data = panel_data,
        response = response,
        factor1 = factor1,
        factor2 = factor2,
        facet = NULL,
        test = test,
        test_args = test_args,
        normality_alpha = normality_alpha,
        min_n = min_n,
        keep = keep,
        p_adjust_method = p_adjust_method,
        alpha = alpha,
        x_layout = x_layout,
        facet_scales = facet_scales,
        dodge_width = dodge_width,
        interaction_sep = interaction_sep
      )
    })

    add_facet <- function(x, facet_value) {
      x[[facet]] <- rep(facet_value, nrow(x))
      x
    }
    bind_component <- function(component) {
      pieces <- Map(function(result, facet_value) {
        add_facet(result[[component]], facet_value)
      }, facet_results, facet_values)
      answer <- do.call(rbind, pieces)
      rownames(answer) <- NULL
      answer
    }

    comparisons <- bind_component("comparisons")
    all_comparisons <- bind_component("all_comparisons")
    groups <- bind_component("groups")
    normality <- bind_component("normality")
    plot_data <- bind_component("plot_data")

    observed_data <- do.call(rbind, Map(function(result, facet_value) {
      answer <- data.frame(.observed_y = result$observed_y)
      answer[[facet]] <- rep(facet_value, nrow(answer))
      answer
    }, facet_results, facet_values))
    rownames(observed_data) <- NULL

    ## Fixed x scales use one coordinate system across every panel. Free x
    ## scales retain the panel-local coordinates calculated by the recursive
    ## calls above.
    if (facet_scales %in% c("fixed", "free_y")) {
      complete_global <- stats::complete.cases(
        facet_data[, c(response, factor1, factor2), drop = FALSE]
      ) & is.finite(facet_data[[response]])
      global_data <- facet_data[complete_global, , drop = FALSE]

      level_values <- function(x) {
        if (is.factor(x)) levels(droplevels(x)) else unique(as.character(x))
      }
      f1_levels <- level_values(global_data[[factor1]])
      f2_levels <- if (is.null(factor2)) NULL else level_values(global_data[[factor2]])

      coordinate <- function(f1, f2 = NULL) {
        f1_index <- match(as.character(f1), f1_levels)
        if (is.null(factor2)) return(f1_index)
        f2_index <- match(as.character(f2), f2_levels)
        if (x_layout == "dodge") {
          n2 <- length(f2_levels)
          return(f1_index + (f2_index - (n2 + 1) / 2) * (dodge_width / n2))
        }
        grid <- unique(data.frame(
          f1 = as.character(global_data[[factor1]]),
          f2 = as.character(global_data[[factor2]]),
          stringsAsFactors = FALSE
        ))
        grid$f1_index <- match(grid$f1, f1_levels)
        grid$f2_index <- match(grid$f2, f2_levels)
        grid <- grid[order(grid$f1_index, grid$f2_index), , drop = FALSE]
        lookup <- stats::setNames(
          seq_len(nrow(grid)), paste(grid$f1, grid$f2, sep = "\r")
        )
        unname(lookup[paste(as.character(f1), as.character(f2), sep = "\r")])
      }

      remap_comparisons <- function(x) {
        x$x1 <- coordinate(x$f1_1, if (is.null(factor2)) NULL else x$f2_1)
        x$x2 <- coordinate(x$f1_2, if (is.null(factor2)) NULL else x$f2_2)
        x$xmin <- pmin(x$x1, x$x2)
        x$xmax <- pmax(x$x1, x$x2)
        x
      }
      groups$x <- coordinate(
        groups$factor1,
        if (is.null(factor2)) NULL else groups$factor2
      )
      comparisons <- remap_comparisons(comparisons)
      all_comparisons <- remap_comparisons(all_comparisons)

      plot_data$.signif_fill <- if (is.null(factor2)) {
        factor(rep("(one factor)", nrow(plot_data)))
      } else {
        factor(plot_data[[factor2]], levels = f2_levels)
      }
      if (is.null(factor2) || x_layout == "dodge") {
        plot_data$.signif_x <- factor(plot_data[[factor1]], levels = f1_levels)
      } else {
        grid <- unique(data.frame(
          f1 = as.character(global_data[[factor1]]),
          f2 = as.character(global_data[[factor2]]),
          stringsAsFactors = FALSE
        ))
        grid$f1_index <- match(grid$f1, f1_levels)
        grid$f2_index <- match(grid$f2, f2_levels)
        grid <- grid[order(grid$f1_index, grid$f2_index), , drop = FALSE]
        group_labels <- make.unique(paste(grid$f1, grid$f2, sep = interaction_sep))
        label_lookup <- stats::setNames(
          group_labels, paste(grid$f1, grid$f2, sep = "\r")
        )
        plot_labels <- unname(label_lookup[paste(
          as.character(plot_data[[factor1]]),
          as.character(plot_data[[factor2]]),
          sep = "\r"
        )])
        plot_data$.signif_x <- factor(plot_labels, levels = group_labels)
      }
    }

    methods <- stats::setNames(
      vapply(facet_results, function(x) x$method, character(1)),
      facet_values
    )
    return(structure(
      list(
        comparisons = comparisons,
        all_comparisons = all_comparisons,
        groups = groups,
        normality = normality,
        plot_data = plot_data,
        observed_y = observed_data$.observed_y,
        observed_data = observed_data,
        method = methods,
        settings = list(
          x_layout = x_layout,
          dodge_width = dodge_width,
          facet = facet,
          facet_scales = facet_scales
        )
      ),
      class = "pairwise_comparison_result"
    ))
  }

  as_plot_factor <- function(x) {
    if (is.factor(x)) droplevels(x) else factor(x, levels = unique(x[!is.na(x)]))
  }

  needed <- c(response, factor1, factor2)
  ok <- stats::complete.cases(data[, needed, drop = FALSE]) &
    is.finite(data[[response]])
  d <- data[ok, , drop = FALSE]
  if (nrow(d) == 0L) stop("No complete, finite observations remain.")

  d$.sig_f1 <- as_plot_factor(d[[factor1]])
  if (nlevels(d$.sig_f1) < 2L && is.null(factor2)) {
    stop("At least two observed groups are required.")
  }
  if (!is.null(factor2)) d$.sig_f2 <- as_plot_factor(d[[factor2]])

  f1_code <- as.integer(d$.sig_f1)
  f2_code <- if (is.null(factor2)) rep.int(1L, nrow(d)) else as.integer(d$.sig_f2)
  row_key <- paste(f1_code, f2_code, sep = ":")
  key_table <- unique(data.frame(f1_code = f1_code, f2_code = f2_code))
  key_table <- key_table[order(key_table$f1_code, key_table$f2_code), , drop = FALSE]
  group_key <- paste(key_table$f1_code, key_table$f2_code, sep = ":")
  d$.sig_group <- match(row_key, group_key)

  groups <- data.frame(
    group = seq_len(nrow(key_table)),
    factor1 = levels(d$.sig_f1)[key_table$f1_code],
    stringsAsFactors = FALSE
  )
  if (!is.null(factor2)) {
    groups$factor2 <- levels(d$.sig_f2)[key_table$f2_code]
  }

  if (is.null(factor2)) {
    groups$x <- key_table$f1_code
  } else if (x_layout == "dodge") {
    n2 <- nlevels(d$.sig_f2)
    groups$x <- key_table$f1_code +
      (key_table$f2_code - (n2 + 1) / 2) * (dodge_width / n2)
  } else {
    groups$x <- seq_len(nrow(groups))
  }

  values <- split(d[[response]], d$.sig_group)
  values <- values[as.character(seq_len(nrow(groups)))]
  groups$n <- vapply(values, length, integer(1))
  groups$y_max <- vapply(values, max, numeric(1), na.rm = TRUE)
  if (nrow(groups) < 2L) stop("At least two observed groups are required.")

  plot_data <- data
  plot_f1 <- factor(plot_data[[factor1]], levels = levels(d$.sig_f1))
  plot_f2 <- if (is.null(factor2)) {
    factor(rep("(one factor)", nrow(plot_data)))
  } else {
    factor(plot_data[[factor2]], levels = levels(d$.sig_f2))
  }
  plot_key <- paste(
    as.integer(plot_f1),
    if (is.null(factor2)) 1L else as.integer(plot_f2),
    sep = ":"
  )
  plot_gid <- match(plot_key, group_key)

  if (is.null(factor2) || x_layout == "dodge") {
    plot_data$.signif_x <- plot_f1
  } else {
    raw_labels <- paste(groups$factor1, groups$factor2, sep = interaction_sep)
    group_labels <- make.unique(raw_labels)
    plot_data$.signif_x <- factor(group_labels[plot_gid], levels = group_labels)
  }
  plot_data$.signif_y <- plot_data[[response]]
  plot_data$.signif_fill <- plot_f2
  plot_data$.signif_group <- plot_gid

  pair_index <- utils::combn(seq_len(nrow(groups)), 2L)
  comp <- data.frame(
    group1 = pair_index[1L, ],
    group2 = pair_index[2L, ],
    stringsAsFactors = FALSE
  )
  comp$f1_1 <- groups$factor1[comp$group1]
  comp$f1_2 <- groups$factor1[comp$group2]
  if (!is.null(factor2)) {
    comp$f2_1 <- groups$factor2[comp$group1]
    comp$f2_2 <- groups$factor2[comp$group2]
  }
  comp$x1 <- groups$x[comp$group1]
  comp$x2 <- groups$x[comp$group2]
  comp$xmin <- pmin(comp$x1, comp$x2)
  comp$xmax <- pmax(comp$x1, comp$x2)

  if (is.null(keep)) {
    keep_row <- if (is.null(factor2)) {
      rep(TRUE, nrow(comp))
    } else {
      comp$f1_1 == comp$f1_2 | comp$f2_1 == comp$f2_2
    }
  } else {
    if (!is.function(keep)) stop("`keep` must be NULL or a function.")
    keep_row <- keep(comp)
    if (!is.logical(keep_row) || length(keep_row) != nrow(comp) || anyNA(keep_row)) {
      stop("`keep(comp)` must return one non-missing logical value per comparison.")
    }
  }
  comp <- comp[keep_row, , drop = FALSE]
  rownames(comp) <- NULL
  if (nrow(comp) == 0L) stop("The comparison filter removed every comparison.")

  normality <- data.frame(
    group = groups$group,
    n = groups$n,
    shapiro_p = NA_real_,
    stringsAsFactors = FALSE
  )
  normality$shapiro_p <- vapply(values, function(z) {
    if (length(z) < 3L || length(z) > 5000L || length(unique(z)) < 3L) {
      return(NA_real_)
    }
    tryCatch(stats::shapiro.test(z)$p.value, error = function(e) NA_real_)
  }, numeric(1))

  if (is.function(test)) {
    test_fun <- test
    test_name <- "custom"
  } else {
    test <- match.arg(test)
    if (test == "auto") {
      test <- if (all(!is.na(normality$shapiro_p)) &&
                  all(normality$shapiro_p > normality_alpha)) "t" else "wilcox"
    }
    test_fun <- if (test == "t") stats::t.test else stats::wilcox.test
    test_name <- test
  }

  test_one <- function(i) {
    x <- values[[comp$group1[i]]]
    y <- values[[comp$group2[i]]]
    if (length(x) < min_n || length(y) < min_n) {
      return(list(
        p = NA_real_,
        error = paste0("fewer than ", min_n, " observations in a group")
      ))
    }
    ans <- tryCatch(
      do.call(test_fun, c(list(x = x, y = y), test_args)),
      error = function(e) e
    )
    if (inherits(ans, "error")) {
      list(p = NA_real_, error = conditionMessage(ans))
    } else if (is.null(ans$p.value) || length(ans$p.value) != 1L) {
      list(p = NA_real_, error = "test did not return one `p.value`")
    } else {
      list(p = as.numeric(ans$p.value), error = NA_character_)
    }
  }

  tested <- lapply(seq_len(nrow(comp)), test_one)
  comp$p_value <- vapply(tested, `[[`, numeric(1), "p")
  comp$test_error <- vapply(tested, `[[`, character(1), "error")
  comp$test <- test_name
  comp$padj <- NA_real_
  valid_p <- is.finite(comp$p_value)
  comp$padj[valid_p] <- stats::p.adjust(
    comp$p_value[valid_p], method = p_adjust_method
  )

  all_comparisons <- comp
  comparisons <- comp[!is.na(comp$padj) & comp$padj <= alpha, , drop = FALSE]
  rownames(comparisons) <- NULL

  structure(
    list(
      comparisons = comparisons,
      all_comparisons = all_comparisons,
      groups = groups,
      normality = normality,
      plot_data = plot_data,
      observed_y = d[[response]],
      observed_data = data.frame(.observed_y = d[[response]]),
      method = test_name,
      settings = list(
        x_layout = x_layout,
        dodge_width = dodge_width,
        facet = NULL,
        facet_scales = facet_scales
      )
    ),
    class = "pairwise_comparison_result"
  )
}


#' Position Pairwise Significance Brackets
#'
#' Adds annotations and compact y positions to a filtered pairwise-comparison
#' data frame. The first argument may be the object returned by
#' [make_pairwise_comparisons()] or a comparison data frame accompanied by
#' `groups` and `observed_y`.
#'
#' @param comparisons A `pairwise_comparison_result` returned by
#'   [make_pairwise_comparisons()], or a data frame containing `padj`,
#'   `xmin`, and `xmax` columns.
#' @param groups `NULL` when `comparisons` is a `pairwise_comparison_result`;
#'   otherwise, a data frame containing numeric `x` and `y_max` columns.
#' @param observed_y `NULL` when `comparisons` is a
#'   `pairwise_comparison_result`. Without facets, it may otherwise be a
#'   numeric vector of finite response values. With `facet`, supply a data
#'   frame containing the facet column and either `.observed_y`, `observed_y`,
#'   or exactly one other numeric column; alternatively, supply a named list or
#'   a numeric vector whose names identify facets.
#' @inheritParams compact_pairwise_signif
#'
#' @return The final comparison data frame, including `annotation`, `roof`,
#'   `span`, `y_position`, `group`, and `lane_order`. Without facets, its
#'   recommended two-value y-axis range is stored in `attr(result, "ylim")`.
#'   With facets, that attribute is a named list of panel ranges and
#'   `attr(result, "facet_limits")` is a long data frame suitable for an
#'   invisible `geom_blank()` scale-training layer. Physical layout parameters
#'   are stored in `attr(result, "layout_settings")`.
#'
#' @seealso [make_pairwise_comparisons()], [compact_pairwise_signif()]
#' @export
#' @md
position_pairwise_significance <- function(
    comparisons,
    groups = NULL,
    observed_y = NULL,
    facet = NULL,
    facet_scales = c("free_y", "fixed", "free", "free_x"),
    annotation = c("q", "stars", "none"),
    binarize_significance = FALSE,
    q_digits = 2L,
    textsize = 3.88,
    lineheight = 1.2,
    panel_height_mm = 85,
    data_gap_mm = 2,
    bar_gap_mm = 1.5,
    top_margin_mm = 1.5,
    y_min = NULL) {

  facet_scales_missing <- missing(facet_scales)
  if (inherits(comparisons, "pairwise_comparison_result")) {
    comparison_result <- comparisons
    comparisons <- comparison_result$comparisons
    if (is.null(groups)) groups <- comparison_result$groups
    if (is.null(facet)) facet <- comparison_result$settings$facet
    if (facet_scales_missing &&
        !is.null(comparison_result$settings$facet_scales)) {
      facet_scales <- comparison_result$settings$facet_scales
    }
    if (is.null(observed_y)) {
      observed_y <- if (!is.null(facet)) {
        comparison_result$observed_data
      } else {
        comparison_result$observed_y
      }
    }
  }

  if (!is.data.frame(comparisons)) stop("`comparisons` must be a data frame or comparison result.")
  if (!is.data.frame(groups)) stop("`groups` must be a data frame.")
  comparison_columns <- c("padj", "xmin", "xmax")
  group_columns <- c("x", "y_max")
  if (!all(comparison_columns %in% names(comparisons))) {
    stop("`comparisons` is missing: ",
         paste(setdiff(comparison_columns, names(comparisons)), collapse = ", "))
  }
  if (!all(group_columns %in% names(groups))) {
    stop("`groups` is missing: ",
         paste(setdiff(group_columns, names(groups)), collapse = ", "))
  }
  if (!all(vapply(comparisons[comparison_columns], is.numeric, logical(1)))) {
    stop("`padj`, `xmin`, and `xmax` in `comparisons` must be numeric.")
  }
  if (!all(vapply(groups[group_columns], is.numeric, logical(1)))) {
    stop("`x` and `y_max` in `groups` must be numeric.")
  }
  if (any(!is.finite(as.matrix(comparisons[comparison_columns])))) {
    stop("`padj`, `xmin`, and `xmax` in `comparisons` must be finite.")
  }
  if (any(!is.finite(as.matrix(groups[group_columns])))) {
    stop("`x` and `y_max` in `groups` must be finite.")
  }
  if (any(comparisons$xmin > comparisons$xmax)) {
    stop("Every `xmin` must be less than or equal to its `xmax`.")
  }
  facet_scales <- match.arg(facet_scales)

  if (!is.null(facet)) {
    if (!is.character(facet) || length(facet) != 1L || is.na(facet)) {
      stop("`facet` must be NULL or one column name.")
    }
    if (!facet %in% names(comparisons)) {
      stop("`comparisons` does not contain the facet column `", facet, "`.")
    }
    if (!facet %in% names(groups)) {
      stop("`groups` does not contain the facet column `", facet, "`.")
    }

    observed_by_facet <- if (is.data.frame(observed_y)) {
      if (!facet %in% names(observed_y)) {
        stop("Faceted `observed_y` data does not contain `", facet, "`.")
      }
      y_column <- if (".observed_y" %in% names(observed_y)) {
        ".observed_y"
      } else if ("observed_y" %in% names(observed_y)) {
        "observed_y"
      } else {
        candidates <- setdiff(
          names(observed_y)[vapply(observed_y, is.numeric, logical(1))],
          facet
        )
        if (length(candidates) != 1L) {
          stop("Faceted `observed_y` must contain `.observed_y`, `observed_y`, ",
               "or exactly one non-facet numeric column.")
        }
        candidates
      }
      split(observed_y[[y_column]], as.character(observed_y[[facet]]))
    } else if (is.list(observed_y) && !is.null(names(observed_y))) {
      observed_y
    } else if (is.numeric(observed_y) && !is.null(names(observed_y))) {
      split(unname(observed_y), names(observed_y))
    } else {
      stop("With `facet`, supply `observed_y` as a data frame containing the ",
           "facet and response, a named list, or a numeric vector named by facet.")
    }

    facet_values <- unique(as.character(groups[[facet]]))
    facet_values <- facet_values[!is.na(facet_values)]
    if (!length(facet_values)) stop("No non-missing facet values were found.")
    missing_observed <- facet_values[vapply(
      facet_values,
      function(facet_value) is.null(observed_by_facet[[facet_value]]),
      logical(1)
    )]
    if (length(missing_observed)) {
      stop("No `observed_y` values were supplied for facet(s): ",
           paste(missing_observed, collapse = ", "), ".")
    }

    panel_argument <- function(x, facet_value, argument, allow_null = FALSE) {
      if (is.null(x) && allow_null) return(NULL)
      if (length(x) == 1L) return(unname(x))
      if (!is.null(names(x)) && facet_value %in% names(x)) {
        return(unname(x[[facet_value]]))
      }
      stop("`", argument, "` must have length one or be named by facet.")
    }

    ## Facets with a fixed y scale must be solved against one common axis
    ## span. Temporary, non-overlapping x blocks let the ordinary packing
    ## algorithm lay out every panel at once without treating brackets from
    ## different panels as collisions.
    if (facet_scales %in% c("fixed", "free_x")) {
      common_argument <- function(x, argument, allow_null = FALSE) {
        if (is.null(x) && allow_null) return(NULL)
        values <- unname(x)
        if (!length(values) || anyNA(values) ||
            length(unique(values)) != 1L) {
          stop("With a fixed y scale, `", argument,
               "` must be one value (or identical values for every facet).")
        }
        values[[1L]]
      }

      common_panel_height <- common_argument(
        panel_height_mm, "panel_height_mm"
      )
      common_y_min <- common_argument(y_min, "y_min", allow_null = TRUE)
      temporary_groups <- groups
      temporary_comparisons <- comparisons
      x_values <- c(groups$x, comparisons$xmin, comparisons$xmax)
      x_values <- x_values[is.finite(x_values)]
      x_width <- if (length(x_values)) diff(range(x_values)) else 0
      if (!is.finite(x_width) || x_width <= 0) x_width <- 1
      shifts <- stats::setNames(
        (seq_along(facet_values) - 1) * (2 * x_width + 2),
        facet_values
      )
      group_shift <- unname(shifts[as.character(groups[[facet]])])
      comparison_shift <- unname(
        shifts[as.character(comparisons[[facet]])]
      )
      temporary_groups$x <- temporary_groups$x + group_shift
      temporary_comparisons$xmin <- temporary_comparisons$xmin +
        comparison_shift
      temporary_comparisons$xmax <- temporary_comparisons$xmax +
        comparison_shift

      positioned <- position_pairwise_significance(
        comparisons = temporary_comparisons,
        groups = temporary_groups,
        observed_y = unlist(observed_by_facet[facet_values], use.names = FALSE),
        facet = NULL,
        facet_scales = facet_scales,
        annotation = annotation,
        binarize_significance = binarize_significance,
        q_digits = q_digits,
        textsize = textsize,
        lineheight = lineheight,
        panel_height_mm = common_panel_height,
        data_gap_mm = data_gap_mm,
        bar_gap_mm = bar_gap_mm,
        top_margin_mm = top_margin_mm,
        y_min = common_y_min
      )
      positioned$xmin <- comparisons$xmin
      positioned$xmax <- comparisons$xmax
      positioned$group <- seq_len(nrow(positioned))

      common_ylim <- attr(positioned, "ylim")
      facet_limits <- data.frame(
        y = rep(common_ylim, times = length(facet_values))
      )
      facet_limits[[facet]] <- rep(facet_values, each = 2L)
      attr(positioned, "ylim") <- stats::setNames(
        rep(list(common_ylim), length(facet_values)), facet_values
      )
      attr(positioned, "facet_limits") <- facet_limits
      attr(positioned, "layout_settings") <- list(
        facet = facet,
        facet_scales = facet_scales,
        panel_height_mm = panel_height_mm,
        textsize = textsize,
        binarize_significance = binarize_significance
      )
      return(positioned)
    }

    panel_results <- lapply(facet_values, function(facet_value) {
      comparisons_i <- comparisons[
        as.character(comparisons[[facet]]) == facet_value, , drop = FALSE
      ]
      groups_i <- groups[
        as.character(groups[[facet]]) == facet_value, , drop = FALSE
      ]
      observed_i <- observed_by_facet[[facet_value]]
      if (is.null(observed_i)) {
        stop("No `observed_y` values were supplied for facet `", facet_value, "`.")
      }
      positioned_i <- position_pairwise_significance(
        comparisons = comparisons_i,
        groups = groups_i,
        observed_y = observed_i,
        facet = NULL,
        facet_scales = facet_scales,
        annotation = annotation,
        binarize_significance = binarize_significance,
        q_digits = q_digits,
        textsize = textsize,
        lineheight = lineheight,
        panel_height_mm = panel_argument(
          panel_height_mm, facet_value, "panel_height_mm"
        ),
        data_gap_mm = data_gap_mm,
        bar_gap_mm = bar_gap_mm,
        top_margin_mm = top_margin_mm,
        y_min = panel_argument(y_min, facet_value, "y_min", allow_null = TRUE)
      )
      positioned_i[[facet]] <- rep(facet_value, nrow(positioned_i))
      list(data = positioned_i, ylim = attr(positioned_i, "ylim"))
    })

    positioned <- do.call(rbind, lapply(panel_results, `[[`, "data"))
    rownames(positioned) <- NULL
    positioned$group <- seq_len(nrow(positioned))

    facet_limits <- do.call(rbind, Map(function(result, facet_value) {
      answer <- data.frame(y = result$ylim)
      answer[[facet]] <- rep(facet_value, 2L)
      answer
    }, panel_results, facet_values))
    rownames(facet_limits) <- NULL

    attr(positioned, "ylim") <- stats::setNames(
      lapply(panel_results, `[[`, "ylim"), facet_values
    )
    attr(positioned, "facet_limits") <- facet_limits
    attr(positioned, "layout_settings") <- list(
      facet = facet,
      facet_scales = facet_scales,
      panel_height_mm = panel_height_mm,
      textsize = textsize,
      binarize_significance = binarize_significance
    )
    return(positioned)
  }

  if (!is.numeric(observed_y)) stop("`observed_y` must be numeric.")
  observed_y <- observed_y[is.finite(observed_y)]
  if (!length(observed_y)) stop("`observed_y` has no finite values.")
  if (!is.logical(binarize_significance) ||
      length(binarize_significance) != 1L || is.na(binarize_significance)) {
    stop("`binarize_significance` must be TRUE or FALSE.")
  }
  if (panel_height_mm <= 0 || textsize < 0 || lineheight <= 0 ||
      data_gap_mm < 0 || bar_gap_mm < 0 || top_margin_mm < 0) {
    stop("Physical-size arguments must be non-negative and `panel_height_mm`/`lineheight` positive.")
  }

  if (is.function(annotation)) {
    annotation_fun <- annotation
    annotation_name <- "custom"
  } else {
    annotation_name <- match.arg(annotation)
    annotation_fun <- NULL
  }

  data_min <- min(observed_y)
  data_max <- max(observed_y)
  raw_span <- data_max - data_min
  if (is.null(y_min)) {
    scale_min <- if (raw_span > 0) {
      data_min
    } else {
      data_min - max(abs(data_min), 1) * 0.1
    }
  } else {
    scale_min <- y_min
    if (!is.finite(scale_min) || scale_min > data_min) {
      stop("`y_min` must be finite and no greater than the data minimum.")
    }
  }
  base_span <- data_max - scale_min
  if (!is.finite(base_span) || base_span <= 0) {
    base_span <- max(abs(data_max), 1) * 0.1
  }

  sig <- comparisons
  if (!is.null(annotation_fun)) {
    sig$annotation <- vapply(sig$padj, function(z) {
      ans <- annotation_fun(z)
      if (length(ans) != 1L) {
        stop("The annotation function must return one label per q-value.")
      }
      as.character(ans)
    }, character(1))
  } else if (annotation_name == "q") {
    fmt <- paste0("q = %.", as.integer(q_digits), "g")
    sig$annotation <- sprintf(fmt, sig$padj)
  } else if (annotation_name == "stars") {
    if (binarize_significance) {
      sig$annotation <- rep("*", nrow(sig))
    } else {
      sig$annotation <- ifelse(sig$padj <= 1e-4, "****",
                        ifelse(sig$padj <= 1e-3, "***",
                        ifelse(sig$padj <= 1e-2, "**", "*")))
    }
  } else {
    sig$annotation <- rep("", nrow(sig))
  }

  if (nrow(sig) == 0L) {
    sig$roof <- numeric(0)
    sig$span <- numeric(0)
    sig$y_position <- numeric(0)
    sig$group <- integer(0)
    sig$lane_order <- integer(0)
    attr(sig, "ylim") <- c(scale_min, data_max)
    attr(sig, "layout_settings") <- list(
      panel_height_mm = panel_height_mm,
      textsize = textsize,
      binarize_significance = binarize_significance
    )
    return(sig)
  }

  n_lines <- lengths(strsplit(sig$annotation, "\n", fixed = TRUE))
  annotation_mm <- ifelse(
    nzchar(sig$annotation), textsize * lineheight * n_lines, 0
  )

  sig$roof <- vapply(seq_len(nrow(sig)), function(i) {
    crossed <- groups$x >= sig$xmin[i] & groups$x <= sig$xmax[i]
    if (!any(crossed)) stop("A bracket crosses no group x coordinates.")
    max(groups$y_max[crossed])
  }, numeric(1))
  sig$span <- sig$xmax - sig$xmin

  packing_order <- order(sig$roof, sig$span, sig$xmin, sig$xmax)
  overlaps <- function(i, j) {
    !(sig$xmax[i] < sig$xmin[j] || sig$xmin[i] > sig$xmax[j])
  }

  layout_at_span <- function(axis_span) {
    y <- rep(NA_real_, nrow(sig))
    for (k in seq_along(packing_order)) {
      i <- packing_order[k]
      candidate <- sig$roof[i] + (data_gap_mm / panel_height_mm) * axis_span
      if (k > 1L) {
        previous <- packing_order[seq_len(k - 1L)]
        previous <- previous[vapply(
          previous, function(j) overlaps(i, j), logical(1)
        )]
        if (length(previous)) {
          above_previous <- y[previous] +
            ((annotation_mm[previous] + bar_gap_mm) / panel_height_mm) * axis_span
          candidate <- max(candidate, above_previous)
        }
      }
      y[i] <- candidate
    }
    y
  }

  axis_span <- base_span
  converged <- FALSE
  for (iter in seq_len(200L)) {
    y_now <- layout_at_span(axis_span)
    occupied_top <- max(
      data_max,
      y_now + (annotation_mm / panel_height_mm) * axis_span
    )
    new_span <- occupied_top - scale_min +
      (top_margin_mm / panel_height_mm) * axis_span
    if (!is.finite(new_span) || new_span > base_span * 1e8) break
    if (abs(new_span - axis_span) <= 1e-10 * max(1, axis_span)) {
      axis_span <- new_span
      converged <- TRUE
      break
    }
    axis_span <- new_span
  }
  if (!converged) {
    stop("The requested brackets/text cannot fit in `panel_height_mm`. ",
         "Increase that value (and export the plot at the corresponding height), ",
         "or show fewer comparisons.")
  }

  sig$y_position <- layout_at_span(axis_span)
  sig$group <- seq_len(nrow(sig))
  sig$lane_order <- match(seq_len(nrow(sig)), packing_order)
  attr(sig, "ylim") <- c(scale_min, scale_min + axis_span)
  attr(sig, "layout_settings") <- list(
    panel_height_mm = panel_height_mm,
    textsize = textsize,
    binarize_significance = binarize_significance
  )
  sig
}


#' Test Pairwise Comparisons and Compactly Place Significance Brackets
#'
#' Performs all relevant pairwise tests for a one- or two-factor design,
#' adjusts the resulting p-values, prepares annotations, and uses a
#' collision-aware algorithm to place significance brackets. Brackets begin
#' above the observations they actually cross, and brackets with disjoint
#' horizontal intervals may share the same vertical space.
#'
#' @param data A data frame containing the response and grouping variables.
#' @param response A length-one character string naming the numeric response
#'   column in `data`.
#' @param factor1 A length-one character string naming the first grouping
#'   column in `data`.
#' @param factor2 `NULL` for a one-factor design, or a length-one character
#'   string naming the second grouping column.
#' @param facet `NULL` for an un-faceted analysis, or a length-one character
#'   string naming the facet column. Tests, multiplicity adjustment, and bar
#'   positioning are performed independently within each facet.
#' @param test Either one of `"t"`, `"wilcox"`, or `"auto"`, or a function.
#'   `"t"` uses [stats::t.test()] (Welch's test unless changed through
#'   `test_args`), and `"wilcox"` uses [stats::wilcox.test()]. A custom
#'   function must accept arguments named `x` and `y` and return a list with a
#'   length-one `p.value` component.
#' @param test_args A named list of additional arguments passed to the test
#'   function. For example, `list(var.equal = TRUE)` requests equal-variance
#'   t-tests.
#' @param normality_alpha A length-one numeric value used only when
#'   `test = "auto"`. Shapiro-Wilk tests are applied to all observed cells. A
#'   single test family is then used for every contrast: t-tests only when all
#'   cell-level normality p-values exceed `normality_alpha`, and Wilcoxon tests
#'   otherwise.
#' @param min_n The minimum number of finite observations required in each
#'   member of a comparison. Comparisons failing this requirement are retained
#'   in `all_comparisons` with a missing p-value and an explanatory
#'   `test_error`.
#' @param keep `NULL` or a comparison-filter function. The function receives a
#'   data frame with one row per possible pair and columns `group1`, `group2`,
#'   `f1_1`, `f1_2`, and, for two-factor designs, `f2_1` and `f2_2`. It must
#'   return one non-missing logical value per row. When `NULL`, all one-factor
#'   comparisons are retained; for two factors, comparisons are retained when
#'   factor 1 or factor 2 is equal between the two cells.
#' @param p_adjust_method A method accepted by [stats::p.adjust()]. Adjustment
#'   is performed only after `keep` has removed irrelevant comparisons.
#' @param alpha The adjusted-p-value threshold. Comparisons with
#'   `padj <= alpha` are drawn.
#' @param annotation One of `"q"`, `"stars"`, or `"none"`, or a function that
#'   converts one q-value into one character label. Star thresholds are
#'   `0.05`, `0.01`, `0.001`, and `0.0001`.
#' @param binarize_significance Logical. When `TRUE` and
#'   `annotation = "stars"`, every comparison retained by the adjusted-p-value
#'   filter is labeled with exactly `"*"`, regardless of its significance
#'   level. This does not change testing, FDR adjustment, or the `alpha`
#'   threshold. It has no effect on q-value, blank, or custom annotations.
#' @param q_digits The number of significant digits used when
#'   `annotation = "q"`.
#' @param x_layout For a two-factor design, either `"dodge"` to place factor 2
#'   within factor 1 or `"interaction"` to place every observed cell on a
#'   separate discrete x position. Ignored for a one-factor design.
#' @param facet_scales The scale behavior supplied to ggplot2 faceting:
#'   `"fixed"`, `"free_y"`, `"free_x"`, or `"free"`. With fixed x scales
#'   (`"fixed"` or `"free_y"`), x coordinates use global factor levels. With
#'   free x scales, coordinates are recalculated within each facet. With free
#'   y scales (`"free_y"` or `"free"`), bracket positions and y limits are
#'   solved independently per facet. With fixed y scales, one common y range
#'   is used so physical gaps remain consistent across panels.
#' @param dodge_width Total dodge width used to calculate factor-2 centres.
#'   Use the same value in [ggplot2::position_dodge()] and
#'   [ggplot2::position_jitterdodge()].
#' @param interaction_sep Character string placed between factor values when
#'   constructing labels for `x_layout = "interaction"`.
#' @param textsize Annotation size in millimetres, following the `textsize`
#'   convention used by [ggsignif::geom_signif()].
#' @param lineheight Multiplicative line-height estimate used to reserve
#'   vertical annotation space.
#' @param panel_height_mm Intended height, in millimetres, of the rendered
#'   plotting panel. This is the panel containing the data, not the complete
#'   figure including titles, axes, margins, or legends. It is required to
#'   convert physical text sizes into y-axis data units. For unequal-height
#'   free-y facets, this may be a numeric vector named by facet. Fixed-y facets
#'   require one value (or identical values for every facet).
#' @param data_gap_mm Minimum physical gap, in millimetres, between the local
#'   data maximum and the first bracket above it.
#' @param bar_gap_mm Minimum physical gap, in millimetres, between the top of a
#'   lower bracket annotation and the next overlapping bracket.
#' @param top_margin_mm Minimum physical space, in millimetres, retained above
#'   the highest annotation.
#' @param y_min `NULL` to use the observed response minimum as the lower limit
#'   for layout calculations, or a finite numeric value no greater than the
#'   observed minimum. Supplying the intended plot minimum, commonly zero for
#'   non-negative measurements, improves physical-size calibration. For
#'   free-y facets this may be a numeric vector named by facet; fixed-y facets
#'   require a common value.
#' @param make_layer Logical. If `TRUE`, create a [ggsignif::geom_signif()]
#'   layer. If `FALSE`, calculate tests and positions without requiring
#'   ggplot2 or ggsignif.
#' @param geom_args A named list of arguments passed to
#'   [ggsignif::geom_signif()]. These values override the helper's layer
#'   defaults; for example, `list(tip_length = 0.01, colour = "navy")`.
#'
#' @details
#' The wrapper delegates statistical work to [make_pairwise_comparisons()] and
#' passes its result to [position_pairwise_significance()]. Both helpers are
#' exported, so callers may inspect or modify the filtered comparisons between
#' the two stages. The wrapper then creates the optional ggsignif layer and
#' assembles the backward-compatible result object.
#'
#' Rows with missing grouping values, missing responses, or non-finite
#' responses are excluded from testing and layout calculations. The input
#' factors are tracked separately rather than parsed from a concatenated group
#' name, so factor levels may safely contain underscores or other separators.
#'
#' `test = "auto"` implements a deliberately conservative convenience rule;
#' it does not select a different test for each contrast. In confirmatory work,
#' choosing the test in advance is generally preferable to normality
#' pre-testing.
#'
#' For every retained comparison, the initial bracket roof is the largest
#' observation among all plotted group centres crossed by that bracket. The
#' comparisons are ordered by roof and horizontal span, then greedily assigned
#' to the lowest non-colliding vertical position. Text height and gaps are
#' converted from millimetres by solving against the final y range. If the
#' requested brackets cannot physically fit in `panel_height_mm`, the function
#' stops and asks for a taller panel instead of silently allowing overlap.
#'
#' The layout assumes a continuous linear y scale, unpaired samples, and
#' vertical significance brackets. Without facets, use
#' `coord_cartesian(ylim = result$ylim, expand = FALSE, clip = "off")`. With
#' facets, the returned layer list includes invisible per-panel scale anchors;
#' use `scale_y_continuous(expand = expansion(mult = c(0, 0)))` and do not set
#' one global y limit.
#'
#' @return An object of class `compact_pairwise_signif`, which is a list with
#'   the following components:
#'   \describe{
#'     \item{comparisons}{A data frame containing the significant comparisons,
#'       raw and adjusted p-values, annotations, bracket coordinates, local
#'       roofs, and packing order.}
#'     \item{all_comparisons}{A data frame containing every biologically
#'       retained comparison, including non-significant results and test
#'       errors.}
#'     \item{groups}{One row per observed factor cell, with its plotting
#'       coordinate, sample size, and observed maximum.}
#'     \item{normality}{Cell sample sizes and Shapiro-Wilk p-values used by the
#'       automatic test-selection rule.}
#'     \item{plot_data}{A copy of `data` with standardized `.signif_x`,
#'       `.signif_y`, `.signif_fill`, and `.signif_group` columns.}
#'     \item{layer}{A ggsignif layer for an un-faceted result. For facets, a
#'       list containing a scale-anchor layer and, when needed, the ggsignif
#'       layer. It is `NULL` when `make_layer = FALSE`.}
#'     \item{ylim}{A two-value range without facets, or a named list of ranges
#'       with facets.}
#'     \item{facet_limits}{Panel-specific scale limits in long form, or `NULL`.}
#'     \item{method}{The test family actually used.}
#'     \item{settings}{Selected layout and physical-size settings.}
#'   }
#'
#' @seealso [make_pairwise_comparisons()], [position_pairwise_significance()],
#'   [stats::t.test()], [stats::wilcox.test()], [stats::p.adjust()],
#'   [ggsignif::geom_signif()]
#' @keywords htest
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("ggsignif", quietly = TRUE)) {
#'   dat <- transform(ToothGrowth, dose = factor(dose))
#'
#'   sig <- compact_pairwise_signif(
#'     data = dat,
#'     response = "len",
#'     factor1 = "dose",
#'     test = "t",
#'     annotation = "stars",
#'     binarize_significance = TRUE,
#'     panel_height_mm = 85,
#'     y_min = 0
#'   )
#'
#'   ggplot2::ggplot(
#'     sig$plot_data,
#'     ggplot2::aes(x = .signif_x, y = .signif_y)
#'   ) +
#'     ggplot2::geom_boxplot() +
#'     sig$layer +
#'     ggplot2::coord_cartesian(
#'       ylim = sig$ylim,
#'       expand = FALSE,
#'       clip = "off"
#'     )
#' }
#' @md
compact_pairwise_signif <- function(
    data,
    response,
    factor1,
    factor2 = NULL,
    facet = NULL,
    test = c("t", "wilcox", "auto"),
    test_args = list(),
    normality_alpha = 0.05,
    min_n = 2L,
    keep = NULL,
    p_adjust_method = "fdr",
    alpha = 0.05,
    annotation = c("q", "stars", "none"),
    binarize_significance = FALSE,
    q_digits = 2L,
    x_layout = c("dodge", "interaction"),
    facet_scales = c("free_y", "fixed", "free", "free_x"),
    dodge_width = 0.75,
    interaction_sep = "_",
    textsize = 3.88,
    lineheight = 1.2,
    panel_height_mm = 85,
    data_gap_mm = 2,
    bar_gap_mm = 1.5,
    top_margin_mm = 1.5,
    y_min = NULL,
    make_layer = TRUE,
    geom_args = list()) {

  if (!is.logical(make_layer) || length(make_layer) != 1L || is.na(make_layer)) {
    stop("`make_layer` must be TRUE or FALSE.")
  }
  if (!is.list(geom_args)) stop("`geom_args` must be a list.")

  comparison_result <- make_pairwise_comparisons(
    data = data,
    response = response,
    factor1 = factor1,
    factor2 = factor2,
    facet = facet,
    test = test,
    test_args = test_args,
    normality_alpha = normality_alpha,
    min_n = min_n,
    keep = keep,
    p_adjust_method = p_adjust_method,
    alpha = alpha,
    x_layout = x_layout,
    facet_scales = facet_scales,
    dodge_width = dodge_width,
    interaction_sep = interaction_sep
  )

  positioned_comparisons <- position_pairwise_significance(
    comparisons = comparison_result,
    facet = facet,
    facet_scales = facet_scales,
    annotation = annotation,
    binarize_significance = binarize_significance,
    q_digits = q_digits,
    textsize = textsize,
    lineheight = lineheight,
    panel_height_mm = panel_height_mm,
    data_gap_mm = data_gap_mm,
    bar_gap_mm = bar_gap_mm,
    top_margin_mm = top_margin_mm,
    y_min = y_min
  )
  ylim <- attr(positioned_comparisons, "ylim")
  facet_limits <- attr(positioned_comparisons, "facet_limits")
  layout_settings <- attr(positioned_comparisons, "layout_settings")

  layer <- NULL
  if (make_layer && (!is.null(facet) || nrow(positioned_comparisons) > 0L)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Install `ggplot2`, or call with `make_layer = FALSE`.")
    }
    significance_layer <- NULL
    if (nrow(positioned_comparisons) > 0L) {
      if (!requireNamespace("ggsignif", quietly = TRUE)) {
        stop("Install `ggsignif`, or call with `make_layer = FALSE`.")
      }
      layer_defaults <- list(
        mapping = ggplot2::aes(
          xmin = xmin, xmax = xmax, annotations = annotation,
          y_position = y_position, group = group
        ),
        data = positioned_comparisons,
        manual = TRUE,
        inherit.aes = FALSE,
        margin_top = 0,
        step_increase = 0,
        tip_length = 0,
        textsize = textsize,
        vjust = 0
      )
      significance_layer <- withCallingHandlers(
        do.call(
          ggsignif::geom_signif,
          utils::modifyList(layer_defaults, geom_args)
        ),
        warning = function(w) {
          msg <- conditionMessage(w)
          expected <- grepl("Ignoring unknown aesthetics", msg, fixed = TRUE) &&
            all(vapply(
              c("xmin", "xmax", "annotations", "y_position"),
              grepl, logical(1), x = msg, fixed = TRUE
            ))
          if (expected) invokeRestart("muffleWarning")
        }
      )
    }

    if (is.null(facet)) {
      layer <- significance_layer
    } else {
      scale_anchor_layer <- ggplot2::geom_blank(
        mapping = ggplot2::aes(y = y),
        data = facet_limits,
        inherit.aes = FALSE
      )
      layer <- Filter(Negate(is.null), list(scale_anchor_layer, significance_layer))
    }
  }

  out <- list(
    comparisons = positioned_comparisons,
    all_comparisons = comparison_result$all_comparisons,
    groups = comparison_result$groups,
    normality = comparison_result$normality,
    plot_data = comparison_result$plot_data,
    layer = layer,
    ylim = ylim,
    facet_limits = facet_limits,
    method = comparison_result$method,
    settings = utils::modifyList(comparison_result$settings, layout_settings)
  )
  structure(out, class = "compact_pairwise_signif")
}


#' Print a Compact Pairwise Significance Layout
#'
#' Prints the selected test family and comparison counts for an object returned
#' by [compact_pairwise_signif()].
#'
#' @param x An object of class `compact_pairwise_signif`.
#' @param ... Additional arguments reserved for compatibility with the
#'   [base::print()] generic. They are currently ignored.
#'
#' @return `x`, invisibly.
#' @export
#' @md
print.compact_pairwise_signif <- function(x, ...) {
  cat("Compact pairwise significance layout\n")
  cat("  test:", x$method, "\n")
  cat("  relevant comparisons tested:", nrow(x$all_comparisons), "\n")
  cat("  comparisons retained:", nrow(x$comparisons), "\n")
  invisible(x)
}


# Variables evaluated non-standardly inside ggplot2 aesthetic mappings and in
# the package examples. Declaring them prevents false-positive R CMD check
# notes without adding rlang as a required dependency.
utils::globalVariables(c(
  ".signif_fill", ".signif_x", ".signif_y",
  "annotation", "group", "xmax", "xmin", "y", "y_position"
))
