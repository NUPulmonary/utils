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
    test = c("t", "wilcox", "auto"),
    test_args = list(),
    normality_alpha = 0.05,
    min_n = 2L,
    keep = NULL,
    p_adjust_method = "fdr",
    alpha = 0.05,
    x_layout = c("dodge", "interaction"),
    dodge_width = 0.75,
    interaction_sep = "_") {

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  nm <- c(response, factor1, factor2)
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
  if (is.null(factor2)) x_layout <- "interaction"

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
      method = test_name,
      settings = list(x_layout = x_layout, dodge_width = dodge_width)
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
#'   `pairwise_comparison_result`; otherwise, the finite response values used
#'   to establish the y-axis range.
#' @inheritParams compact_pairwise_signif
#'
#' @return The final comparison data frame, including `annotation`, `roof`,
#'   `span`, `y_position`, `group`, and `lane_order`. Its recommended two-value
#'   y-axis range is stored in `attr(result, "ylim")`, and the physical layout
#'   parameters are stored in `attr(result, "layout_settings")`.
#'
#' @seealso [make_pairwise_comparisons()], [compact_pairwise_signif()]
#' @export
#' @md
position_pairwise_significance <- function(
    comparisons,
    groups = NULL,
    observed_y = NULL,
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

  if (inherits(comparisons, "pairwise_comparison_result")) {
    comparison_result <- comparisons
    comparisons <- comparison_result$comparisons
    if (is.null(groups)) groups <- comparison_result$groups
    if (is.null(observed_y)) observed_y <- comparison_result$observed_y
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
#'   convert physical text sizes into y-axis data units.
#' @param data_gap_mm Minimum physical gap, in millimetres, between the local
#'   data maximum and the first bracket above it.
#' @param bar_gap_mm Minimum physical gap, in millimetres, between the top of a
#'   lower bracket annotation and the next overlapping bracket.
#' @param top_margin_mm Minimum physical space, in millimetres, retained above
#'   the highest annotation.
#' @param y_min `NULL` to use the observed response minimum as the lower limit
#'   for layout calculations, or a finite numeric value no greater than the
#'   observed minimum. Supplying the intended plot minimum, commonly zero for
#'   non-negative measurements, improves physical-size calibration.
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
#' The layout assumes one panel, a continuous linear y scale, unpaired samples,
#' and vertical significance brackets. For exact physical spacing, add the
#' returned layer and use
#' `coord_cartesian(ylim = result$ylim, expand = FALSE, clip = "off")`.
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
#'     \item{layer}{A ggsignif layer, or `NULL` when `make_layer = FALSE` or no
#'       comparison passes `alpha`.}
#'     \item{ylim}{Recommended lower and upper y limits.}
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
    test = test,
    test_args = test_args,
    normality_alpha = normality_alpha,
    min_n = min_n,
    keep = keep,
    p_adjust_method = p_adjust_method,
    alpha = alpha,
    x_layout = x_layout,
    dodge_width = dodge_width,
    interaction_sep = interaction_sep
  )

  positioned_comparisons <- position_pairwise_significance(
    comparisons = comparison_result,
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
  layout_settings <- attr(positioned_comparisons, "layout_settings")

  layer <- NULL
  if (make_layer && nrow(positioned_comparisons) > 0L) {
    if (!requireNamespace("ggplot2", quietly = TRUE) ||
        !requireNamespace("ggsignif", quietly = TRUE)) {
      stop("Install `ggplot2` and `ggsignif`, or call with `make_layer = FALSE`.")
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
    layer <- withCallingHandlers(
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

  out <- list(
    comparisons = positioned_comparisons,
    all_comparisons = comparison_result$all_comparisons,
    groups = comparison_result$groups,
    normality = comparison_result$normality,
    plot_data = comparison_result$plot_data,
    layer = layer,
    ylim = ylim,
    method = comparison_result$method,
    settings = c(comparison_result$settings, layout_settings)
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
  "annotation", "group", "xmax", "xmin", "y_position"
))
