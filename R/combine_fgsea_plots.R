# Run fgsea for one pathway across comparisons and overlay the results.
# Input: a chosen pathway, the full pathway collection, named rank vectors,
#        and plotting/fgsea options.
# Output: a ggplot2 object with curves, gene-hit ticks, ES extrema, and labels.
# Environment: R with fgsea and ggplot2 installed. Comparisons should use the
# same gene universe, ranking direction, and gene-level statistic.

#' Combine fgsea plots for one pathway across comparisons
#'
#' @param pathway Name of the pathway to display in `pathways`, or its character
#'   vector of gene identifiers. A gene vector must match exactly one pathway.
#' @param pathways Named list of all gene sets to test with
#'   [fgsea::fgseaMultilevel()].
#' @param ranks Named list of numeric gene-level ranking vectors. Each vector
#'   must have gene identifiers as names; list names label the comparisons.
#' @param label_size Positive numeric size for the peak labels, passed to
#'   [ggplot2::geom_text()].
#' @param nproc Nonnegative integer number of workers passed to
#'   [fgsea::fgseaMultilevel()]. The default `0` uses fgsea's backend setting.
#' @param min_size Positive integer minimum pathway size passed as `minSize` to
#'   [fgsea::fgseaMultilevel()].
#' @param seed Optional integer seed reset immediately before each comparison's
#'   fgsea run. Use the seed from separately run comparisons to reproduce their
#'   adjusted p-values.
#'
#' @details Runs [fgsea::fgseaMultilevel()] separately for each comparison,
#'   testing the full `pathways` list, then obtains the selected pathway's
#'   running scores and ticks from [fgsea::plotEnrichmentData()]. Adjusted
#'   p-values are calculated over the pathways retained in each comparison.
#'   `fgseaMultilevel()` uses random sampling. If `seed` is `NULL`, the current
#'   random-number state advances across comparisons.
#'
#' @return A ggplot2 object with curves, separate tick rows, dashed lines at
#'   each curve's minimum and maximum ES, and q/NES labels just above the
#'   dashed line with the largest absolute ES for each comparison.
#' @export
combine_fgsea_plots = function(pathway, pathways, ranks, label_size = 3,
                               nproc = 0, min_size = 1, seed = NULL) {
  pathway_names = names(pathways)
  if (!is.list(pathways) || length(pathways) == 0L ||
      is.null(pathway_names) || anyNA(pathway_names) ||
      any(!nzchar(pathway_names)) || anyDuplicated(pathway_names)) {
    stop("pathways must be a nonempty named list with unique names")
  }
  if (!is.character(pathway) || length(pathway) == 0L ||
      anyNA(pathway) || any(!nzchar(pathway))) {
    stop("pathway must be a pathway name or a character vector of genes")
  }
  
  # A named selection is simplest; matching a gene vector preserves prior calls.
  if (length(pathway) == 1L && pathway %in% pathway_names) {
    pathway_name = pathway
  } else {
    matches = pathway_names[vapply(pathways, function(genes) {
      setequal(genes, pathway)
    }, logical(1))]
    if (length(matches) != 1L) {
      stop("pathway genes must match exactly one entry in pathways")
    }
    pathway_name = matches
  }
  selected_pathway = pathways[[pathway_name]]
  
  if (!is.numeric(label_size) || length(label_size) != 1L ||
      !is.finite(label_size) || label_size <= 0) {
    stop("label_size must be one positive, finite number")
  }
  if (!is.numeric(nproc) || length(nproc) != 1L ||
      !is.finite(nproc) || nproc < 0 || nproc != floor(nproc)) {
    stop("nproc must be one nonnegative integer")
  }
  if (!is.numeric(min_size) || length(min_size) != 1L ||
      !is.finite(min_size) || min_size < 1 || min_size != floor(min_size)) {
    stop("min_size must be one positive integer")
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
       seed != floor(seed) || abs(seed) > .Machine$integer.max)) {
    stop("seed must be NULL or one integer within R's integer range")
  }
  
  comparison_names = names(ranks)
  if (!is.list(ranks) || length(ranks) == 0L ||
      is.null(comparison_names) || anyNA(comparison_names) ||
      any(!nzchar(comparison_names)) || anyDuplicated(comparison_names)) {
    stop("ranks must be a nonempty list with unique comparison names")
  }
  
  valid_ranks = vapply(ranks, function(stats) {
    is.numeric(stats) && length(stats) > 1L && all(is.finite(stats)) &&
      !is.null(names(stats)) && !anyNA(names(stats)) &&
      all(nzchar(names(stats))) && !anyDuplicated(names(stats))
  }, logical(1))
  if (!all(valid_ranks)) {
    stop("each rank vector must contain finite numbers with unique gene names")
  }
  
  # Test the full collection, then select the pathway used in the plot.
  analyses = lapply(ranks, function(stats) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    result = fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      nproc = nproc,
      minSize = min_size
    )
    selected_result = result[result$pathway == pathway_name, ]
    if (nrow(selected_result) != 1L) {
      stop("selected pathway was excluded; check min_size and gene overlap")
    }
    enrichment = fgsea::plotEnrichmentData(
      pathway = selected_pathway,
      stats = stats
    )
    list(result = selected_result, enrichment = enrichment)
  })
  
  # Preserve comparison names when assembling the three plot layers.
  curves = do.call(rbind, Map(function(analysis, comparison) {
    data.frame(analysis$enrichment$curve, comparison = comparison)
  }, analyses, comparison_names))
  
  ticks = do.call(rbind, Map(function(analysis, comparison) {
    data.frame(rank = analysis$enrichment$ticks$rank,
               comparison = comparison)
  }, analyses, comparison_names))
  
  extrema = do.call(rbind, Map(function(analysis, comparison) {
    data.frame(
      comparison = comparison,
      yintercept = c(analysis$enrichment$negES,
                     analysis$enrichment$posES)
    )
  }, analyses, comparison_names))
  
  # Match each label to the dashed ES line with the largest magnitude.
  annotation_data = do.call(rbind, Map(function(analysis, comparison) {
    neg_es = analysis$enrichment$negES
    pos_es = analysis$enrichment$posES
    extreme_es = if (abs(neg_es) > abs(pos_es)) neg_es else pos_es
    nes = analysis$result$NES
    is_negative = if (is.na(nes)) extreme_es < 0 else nes < 0
    data.frame(
      line_es = extreme_es,
      is_negative = is_negative,
      comparison = comparison,
      annotation = paste0(
        "q = ", trimws(format(analysis$result$padj,
                              digits = 3, scientific = TRUE)),
        "   NES = ", sprintf("%.2f", nes)
      )
    )
  }, analyses, comparison_names))
  
  # Stagger gene-hit ticks below the curves so comparisons remain distinct.
  score_span = diff(range(curves$ES))
  if (score_span == 0) {
    score_span = 1
  }
  tick_height = 0.05 * score_span
  ticks$y = min(curves$ES) -
    1.5 * match(ticks$comparison, comparison_names) * tick_height
  ticks$yend = ticks$y + 0.6 * tick_height
  
  # Labels sit above their dashed lines, against the left or right plot edge.
  annotation_data$label_x = max(curves$rank) *
    ifelse(annotation_data$is_negative, 0.02, 0.98)
  annotation_data$label_y = annotation_data$line_es + 0.015 * score_span
  annotation_data$label_hjust = ifelse(annotation_data$is_negative, 0, 1)
  
  # Keep the legend in the same order as the supplied rank vectors.
  curves$comparison = factor(curves$comparison, levels = comparison_names)
  ticks$comparison = factor(ticks$comparison, levels = comparison_names)
  extrema$comparison = factor(extrema$comparison, levels = comparison_names)
  annotation_data$comparison = factor(annotation_data$comparison,
                                      levels = comparison_names)
  
  ggplot2::ggplot(curves, ggplot2::aes(x = rank, y = ES,
                                       color = comparison)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_segment(
      data = ticks,
      mapping = ggplot2::aes(x = rank, xend = rank, y = y, yend = yend,
                             color = comparison),
      inherit.aes = FALSE,
      linewidth = 0.3,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(
      data = extrema,
      mapping = ggplot2::aes(yintercept = yintercept, color = comparison),
      linetype = "dashed",
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60") +
    ggplot2::geom_text(
      data = annotation_data,
      mapping = ggplot2::aes(x = label_x, y = label_y, label = annotation,
                             color = comparison,
                             hjust = label_hjust),
      inherit.aes = FALSE,
      size = label_size,
      vjust = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.12, 0.15))
    ) +
    ggplot2::labs(x = "Gene rank", y = "Running enrichment score",
                  color = "Comparison") +
    ggplot2::theme_bw()
}