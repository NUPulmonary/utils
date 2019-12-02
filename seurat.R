library(cowplot)
library(ggplot2)
library(Seurat)
library(SearchTrees)


..OUT_DIR <- NULL

setOutDir <- function(dir) {
    assign("..OUT_DIR", dir, envir = .GlobalEnv)
}

rasterizePlotWithLegend <- function(p) {
    legend <- get_legend(p)
    p <- AugmentPlot(p)
    return(plot_grid(p, legend, rel_widths = c(1, 0.3)))
}

plotWrapper <- function(name, f) {
    if (is.null(..OUT_DIR)) {
        stop("out_dir is missing, cannot save pdfs")
    }
    pdf(paste(..OUT_DIR, paste0(name, ".pdf"), sep = "/"), width = 12, height = 10)
    plot(f())
    dev.off()
}

plotMarkers <- function(obj, markers = NULL, meta = NULL) {
    if (is.null(markers)) {
        markers <- "FABP4, SPP1, CCL2, CCL3, CSF1R, CSF1, CD1C, FCN1, CD3E,
                    EPCAM, CD79A, JCHAIN, MARCO, MRC1, KIT, CHI3L1, MMP9,
                    CLEC4C, GZMB, CCL22, CCR7, IRF4, CLEC9A, XCR1, THBD, CLEC10A,
                    FCGR2B, nCount_RNA, nFeature_RNA"
    }

    for (marker in strsplit(markers, ",\\s*")[[1]]) {
        plotWrapper(marker, function() {
            rasterizePlotWithLegend(plotGene(
                obj,
                marker,
                meta
            ) + ggtitle(marker))
        })
    }
}

plotCondition <- function(obj, meta, min.cutoff = 0, max.cutoff = 1) {
    unique_meta <- unique(unlist(obj[[meta]]))
    if (length(unique_meta) != 2) {
        stop(paste("Field", meta, "must have 2 unique values, but has", length(unique_meta)))
    }
    result <- c()
    data <- FetchData(obj, c("UMAP_1", "UMAP_2", meta))
    for (x in unique_meta) {
        tree <- createTree(data[data[[meta]] == x, ])
        neighbors <- knnLookup(tree, newdat = data[data[[meta]] != x, ], k = 1)
        coords <- data.frame(
            data[data[[meta]] != x, c("UMAP_1", "UMAP_2")],
            data[data[[meta]] == x, ][neighbors, c("UMAP_1", "UMAP_2")]
        )
        dists <- apply(coords, 1, function(row) {
            sqrt((row[3] - row[1])**2 + (row[4] - row[2])**2)
        })
        if (x == unique_meta[1]) {
            dists <- dists * -1
        }
        r_names <- names(result)
        result <- c(
            result,
            dists
        )
        names(result) <- c(
            r_names,
            rownames(coords)
        )
    }
    qq <- quantile(result, probs = seq(min.cutoff, max.cutoff))
    result[result < qq[1]] <- qq[1]
    result[result > qq[2]] <- qq[2]
    data$nndist <- result[rownames(data)]
    data <- data[order(abs(data$nndist)), ]
    result_range <- range(result)
    midpoint <- abs(result_range[1]) / (result_range[2] - result_range[1])
    return(
        ggplot(data, aes(x = UMAP_1, y = UMAP_2, col = nndist)) +
            geom_point(size = 0.5) +
            scale_color_gradientn(
                colors = c("red", "white", "blue"),
                values = c(0, midpoint, 1),
                breaks = c(result_range[1], 0, result_range[2]),
                labels = c(unique_meta[2], "", unique_meta[1])
            ) +
            theme_classic()
    )
}

plotGene <- function(obj, gene, meta) {
  unique_meta <- unique(unlist(obj[[meta]]))
  if (length(unique_meta) != 2) {
    stop(paste("Field", meta, "must have 2 unique values, but has", length(unique_meta)))
  }

  data <- FetchData(obj, c("UMAP_1", "UMAP_2", meta, gene))
  data[data[[meta]] == unique_meta[2], gene] <- data[[gene]][data[[meta]] == unique_meta[2]] * -1
  data <- data[order(abs(data[[gene]])), ]
  result_range <- range(data[[gene]])
  return(
    ggplot(data, aes_string(x = "UMAP_1", y = "UMAP_2", col = gene)) +
      geom_point(size = 1.5) +
      scale_color_gradient2(
        low = "red", high = "blue",
        mid = "lightgray",
        breaks = c(result_range[1], 0, result_range[2]),
        labels = c(
          paste0(unique_meta[2], " max (", sprintf("%.2f", result_range[1] * -1), ")"),
          "",
          paste0(unique_meta[1], " max (", sprintf("%.2f", result_range[2]), ")")
        )
      ) +
      theme_classic()
  )
}

plotClusterComposition <- function(obj, meta, relative = FALSE) {
    data <- FetchData(obj, c("ident", meta))
    plot <- ggplot(data, aes_string(x = "ident", fill = meta))
    if (relative) {
        plot <- plot + geom_bar(position = "fill")
    } else {
        plot <- plot + geom_bar()
    }
    return(
        plot + theme_classic()
    )
}

describe <- function(obj, out_dir, meta, markers = NULL) {
    setOutDir(out_dir)

    DefaultAssay(obj) <- "RNA"
    plotWrapper("clusters", function() { AugmentPlot(DimPlot(obj, label = TRUE, pt.size = 2)) })
    plotMarkers(obj, markers, meta)
    plotWrapper(meta, function() { plotCondition(obj, meta) })
    plotWrapper(
        paste0(meta, "-composition"),
        function() { plotClusterComposition(obj, meta) }
    )
    plotWrapper(
        paste0(meta, "-composition-rel"),
        function() { plotClusterComposition(obj, meta, relative = TRUE) }
    )

    setOutDir(NULL)
}
