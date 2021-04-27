library(edgeR)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(fgsea)
library(stylo)
library(gplots)
library(pheatmap)


..OUT_DIR <- NULL

set_out_dir <- function(dir) {
    assign("..OUT_DIR", dir, envir = .GlobalEnv)
}

plot_wrapper <- function(name, f, width = 12, height = 10) {
    if (is.null(..OUT_DIR)) {
        stop("out_dir is missing, cannot save pdfs")
    }

    pdf(paste(..OUT_DIR, paste0(name, ".pdf"), sep = "/"), width = width, height = height)
    p <- f()
    plot(p)
    dev.off()
    return(p)
}

extract_info_from_name <- function(x, uRNA_pattern, name_pattern, columns, group_pattern) {
    x <- as.character(x)
    if (grepl(uRNA_pattern, x, perl = TRUE)) {
        return(c(rep(NA, length(columns)), "Uni RNA"))
    }
    if (is.null(name_pattern)) {
        match <- regexec(group_pattern, x, perl = TRUE)
        return(c(rep(NA, length(columns)), regmatches(x, match)[[1]][2]))
    }
    match <- regexec(name_pattern, x, perl = TRUE)
    groups <- regmatches(x, match)[[1]]
    groups <- groups[2:length(groups)]
    return(c(groups, paste(groups, collapse = "_")))
}

load_samples <- function(
    samples_file = "samples.txt",
    uRNA_pattern = "uRNA",
    name_pattern = NULL,
    columns = c(),
    group_pattern = "S\\d\\d\\_(.+)\\_R\\d+",
    sample_id_from_name = FALSE
) {
    samples <- read.delim(samples_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1:2]
    colnames(samples) <- c("sample_id", "sample_name")
    if (sample_id_from_name) {
        samples$sample_id <- substr(samples$sample_name, 1, 3)
    }
    samples <- cbind(
        samples,
        t(as.data.frame(lapply(
            samples$sample_name,
            extract_info_from_name,
            uRNA_pattern = uRNA_pattern,
            name_pattern = name_pattern,
            columns = columns,
            group_pattern = group_pattern
        )))
    )
    rownames(samples) <- NULL
    colnames(samples) <- c("sample_id", "sample_name", columns, "group")
    return(samples)
}


load_counts <- function(file, samples = NULL) {
    if (substr(file, nchar(file) - 3, nchar(file)) == ".csv") {
        counts <- read.delim(file, header = TRUE, sep = ",", row.names = 1)
    } else {
        counts <- read.delim(file, header = TRUE, sep = "\t")
    }

    if (colnames(counts)[1] == "Chr") {
        counts <- counts[, -(1:5)]
    }
    if (!is.null(samples)) {
        if (base::all(colnames(counts) %in% samples$sample_name)
            && ncol(counts) == nrow(samples))
        {
            counts <- counts[, samples$sample_name]
        } else {
            colnames(counts) <- samples$sample_name
        }
    }
    return(counts)
}


plot_cor <- function(data) {
    cor.matrix <- lsa::cosine(as.matrix(data))
    corrplot::corrplot(
        cor.matrix,
        method = "color",
        tl.cex = .8,
        tl.col = "black",
        cl.lim = c(0,1),
        tl.srt = 45
    )
}


plot_lib_sizes <- function(data) {
    df <- data.frame(lib.size = sort(colSums(data)))
        ggplot(df, aes(x = reorder(rownames(df), lib.size), y = lib.size)) +
        geom_col() +
        ggtitle("Sorted library sizes") +
        theme(axis.text.x = element_text(angle = 90, vjust = .5)) +
        xlab("")
}


get_anno <- function(counts, organism = "mmusculus") {
    mart <- useEnsembl(biomart = "ensembl",
                       dataset = sprintf("%s_gene_ensembl", organism),
                       mirror = "useast") # options: useast, uswest, asia

    anno <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = rownames(counts),
                  mart = mart)
    return(anno)
}


create_experiment <- function(samples,
                              counts,
                              anno,
                              exclude = c(),
                              group_column = "group",
                              remove_urna = TRUE
) {
    group <- samples[[group_column]]
    keep_idx <- !(samples$sample_id %in% exclude)
    if (remove_urna) {
        keep_idx <- keep_idx & (group != "Uni RNA")
    }
    columns <- as.character(samples$sample_name[keep_idx])
    samples <- samples[samples$sample_name %in% columns, ]
    group <- samples[[group_column]]
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    design <- design[, colnames(design)[colSums(design) > 0]]
    rownames(design) <- samples$sample_name
    return(list(
        samples = samples,
        counts = counts[, columns],
        design = design,
        group_column = group_column,
        anno = anno
    ))
}


add_dge <- function(exp) {
    if (!is.list(exp) && "samples" %in% names(exp)) {
        dge <- DGEList(exp$counts, group = exp$samples[[exp$group_column]])
    } else {
        dge <- DGEList(exp$counts)
    }
    dge$genes <- data.frame(
        external_gene_name = exp$anno[match(rownames(dge), exp$anno$ensembl_gene_id), "external_gene_name"]
    )
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = exp$design)
    exp$dge <- dge
    return(exp)
}


add_pca <- function(exp, top = 500) {
    data <- cpm(exp$dge, log = TRUE)
    gene.vars <- rowSums((data - rowMeans(data))^2)
    keep_idx <- order(gene.vars, decreasing = TRUE)[1:min(top, nrow(data))]

    exp$pca <- prcomp(t(data[keep_idx, ]), center = TRUE, scale = TRUE)
    exp$pca.var <- exp$pca$sdev^2 / sum(exp$pca$sdev^2)
    return(exp)
}


add_kmeans <- function(exp, clusters) {
    set.seed(42)

    contrasts <- c()
    group_names <- colnames(exp$design)
    first_group = group_names[1]
    for (col in group_names[2:length(group_names)]) {
        contrasts <- c(
            contrasts,
            sprintf("%s-%s", first_group, col)
        )
    }
    contrasts <- makeContrasts(contrasts = contrasts, levels = exp$design)

    fit <- glmFit(exp$dge, design = exp$design)
    res <- glmLRT(fit, contrast = contrasts)
    most.de <- topTags(res, n = nrow(res$table))
    sign <- rownames(most.de$table)[most.de$table$FDR < .05]

    sign.data <- cpm(exp$dge, log = TRUE)[sign, ]
    sign.data <- t(scale(t(sign.data)))
    sign.data <- sign.data[!is.na(sign.data[, 1]) & !is.na(rownames(sign.data)), ]
    z <- kmeans(sign.data, clusters, nstart = 500, iter.max = 50)
    sign.data <- sign.data[names(sort(z$cluster)), ]

    gene_table <- NULL
    for (i in sort(unique(z$cluster))) {
        cl <- most.de$table[names(z$cluster[z$cluster == i]),]
        cl$cluster <- i
        cl <- cl[order(cl$FDR), c("external_gene_name", "cluster", "logCPM", "LR", "FDR")]
        if (is.null(gene_table)) {
            gene_table <- cl
        } else {
            gene_table <- rbind(gene_table, cl)
        }
    }

    exp$kmeans_degs <- most.de$table
    exp$kmeans_data <- sign.data
    exp$kmeans_cluster <- z$cluster
    exp$kmeans_dist <- dist.cosine(t(exp$kmeans_data))
    exp$kmeans_genes <- gene_table

    return(exp)
}


plot_mds <- function(exp, dims = 1:2, label = FALSE) {
    data <- plotMDS(exp$dge, plot = FALSE, dim.plot = dims)
    data <- data.frame(
        x = data$x,
        y = data$y,
        group = exp$samples[[exp$group_column]]
    )
    plot <- ggplot(data, aes(x = x, y = y, col = group)) +
        geom_point() +
        ggtitle("MDS plot of all samples") +
        xlab(sprintf("Dim %d", dims[1])) +
        ylab(sprintf("Dim %d", dims[2])) +
        theme_minimal() +
        coord_fixed()

    if (label) {
        plot <- plot +
            geom_label_repel(aes(label = exp$samples$sample_id), show.legend = FALSE)
    }

    return(plot)
}


plot_pca <- function(exp, dims = 1:2, label = FALSE, pt.size = 1) {
    df <- data.frame(exp$pca$x[, dims])
    colnames(df) <- c("x", "y")
    df <- cbind(df, exp$samples[match(rownames(df), exp$samples$sample_name), ])

    plot <- ggplot(df, aes_string("x", "y", col = exp$group_column)) +
        geom_point(size = pt.size) +
        ggtitle("PCA of all samples") +
        xlab(sprintf(
            "PC%d (%.1f%% explained variance)",
            dims[1],
            exp$pca.var[dims[1]] * 100
        )) +
        ylab(sprintf(
            "PC%d (%.1f%% explained variance)",
            dims[2],
            exp$pca.var[dims[2]] * 100
        )) +
        theme_minimal() +
        coord_fixed()

    if (label) {
        plot <- plot +
            geom_label_repel(aes(label = exp$samples$sample_id), show.legend = FALSE)
    }

    return(plot)
}


plot_gene2 <- function(exp, gene, label=TRUE, group_order=NULL) { # TODO: add raw=TRUE for raw counts
    gene.id <- exp$anno$ensembl_gene_id[exp$anno$external_gene_name == gene][1]
    df <- data.frame(group = exp$samples[[exp$group_column]])
    if (!is.null(group_order)) {
        df$group <- factor(df$group, levels = group_order)
    }
    rownames(df) <- exp$samples$sample_id
    df[[gene]] <- exp$dge$counts[gene.id, ]
    plot <- ggplot(df, aes_string(x = "group", y = sprintf("`%s`", gene))) +
        theme_minimal() +
        geom_boxplot() +
        geom_dotplot(binaxis = "y", stackdir = "center")
    if (label) {
        plot <- plot + geom_label_repel(
            aes(label = rownames(df)),
            na.rm = TRUE,
            nudge_x = .3,
            size = 3
        )
    }
    return(plot)
}


plot_gene <- function(data.obj, data.samples, gene.anno, gene, label=TRUE) {
    gene.id <- gene.anno$ensembl_gene_id[gene.anno$external_gene_name == gene][1]
    df <- data.frame(group=data.samples$group)
    rownames(df) <- data.samples$sample_id
    if (class(data.obj) == "data.frame") { # simple typing TODO: redo
        df[[gene]] <- t(data.obj[gene.id, data.samples$sample_name])
    } else {
        df[[gene]] <- data.obj$dge$counts[gene.id, ]
    }
    plot <- ggplot(df, aes_string(x = "group", y = sprintf("`%s`", gene))) +
        theme_minimal() +
        geom_boxplot() +
        geom_dotplot(binaxis = 'y', stackdir = 'center')
    if (label) {
        plot <- plot + geom_label_repel(
            aes(label = rownames(df)),
            na.rm = TRUE,
            nudge_x = .3,
            size = 3
        )
    }
    return(plot)
}


create_ma_plot <- function(genes) {
    genes$significance <- rep("Not significant", length.out = nrow(genes))
    genes$significance[genes$FDR < .05] <- "FDR < 0.05"
    plot <- ggplot(genes, aes(x = logCPM, y = logFC), size = 6) +
        geom_point(aes(colour = significance), alpha = .5) +
        scale_color_manual(
            values = c("#6194BC", "gray39"),
            labels = c(
                paste0('FDR < 0.05 (', sum(genes$FDR < .05), ')'),
                paste0('Not significant (', sum(genes$FDR >= .05), ')')
            )
        ) +
        theme(
            legend.position = c(0.9, 0.9),
            legend.title=element_blank()
        )
    return(plot)
}


create_vol_plot <- function(genes, up, down, up.name = "", highlight.genes) {
    up <- up[up$FDR < .05, ]
    down <- down[down$FDR < .05, ]

    # grob <- grobTree(textGrob(
    #   paste0("Down: ", nrow(down), "\n", "Up: ", nrow(up)),
    #   x=0,
    #   y=0.95,
    #   hjust=-0.1,
    #   vjust = 1,
    #   gp = gpar(fontface = "italic", fontsize = 16)
    # ))

    #genes$gene <- rownames(genes)
    upregs <- genes$FDR < .05 & genes$logFC > 0
    downregs <- genes$FDR < .05 & genes$logFC < 0
    up.legend <- paste0("Up in ", up.name, " (", sum(upregs), ")")
    down.legend <- paste0("Down in ", up.name, " (", sum(downregs), ")")

    genes$significance <- "Not significant"
    genes$significance[upregs] <- up.legend
    genes$significance[downregs] <- down.legend
    genes$significance <- factor(
        genes$significance,
        levels = c(
            "Not significant",
            down.legend,
            up.legend
        )
    )

    if (length(highlight.genes) > 0) {
        top.genes <- genes[genes$external_gene_name %in% highlight.genes, ]
    } else {
        top.genes <- data.frame(rbind(
            genes[rownames(up)[1:10], ],
            genes[rownames(down)[1:10], ]
        ))
    }

    plot <- ggplot(genes, aes(x = logFC, y = -log10(FDR)), size = 6) +
        geom_point(
            aes(color = significance),
            alpha = .6,
            shape = 20,
            size = 3
        ) +
        theme_minimal() +
        scale_color_manual(
            values = c("#ddccbb", "#306090", "#961e2e"),
            labels = c(
                paste0("Not significant (", sum(genes$significance == "Not significant"), ")"),
                down.legend,
                up.legend
            )
        ) +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 10),
            legend.title = element_blank(),
            plot.title = element_text(size = 16)
        ) +
        geom_text_repel(
            data = top.genes,
            aes(label = external_gene_name),
            size = 5,
            box.padding = unit(0.5, "lines"),
            point.padding = unit(0.1, "lines")
        )
    return(plot)
}


compute_de_genes <- function(exp, contrast, highlight.genes, lfc = 0) {
    up.name <- rownames(contrast)[contrast > 0]
    fit <- glmFit(exp$dge, exp$design)
    if (lfc > 0) {
        genes <- glmTreat(fit, contrast = contrast, lfc = lfc)
    } else {
        genes <- glmLRT(fit, contrast = contrast)
    }
    most.de <- topTags(genes, n = nrow(genes$table))
    upreg <- most.de$table[most.de$table$logFC > 0, ]
    upreg <- upreg[order(upreg$logFC, decreasing = TRUE), ]
    downreg <- most.de$table[most.de$table$logFC < 0, ]
    downreg <- downreg[order(downreg$logFC), ]

    return(list(
        de = genes,
        genes = most.de$table,
        up = upreg,
        down = downreg,
        vol.plot = create_vol_plot(
            most.de$table,
            upreg,
            downreg,
            up.name = up.name,
            highlight.genes
        ),
        ma.plot = create_ma_plot(most.de$table)
    ))
}


compute_de_genes2 <- function(exp, contrast, lfc = 0) {
    fit <- glmFit(exp$dge, exp$design)
    if (lfc > 0) {
        genes <- glmTreat(fit, contrast = contrast, lfc = lfc)
    } else {
        genes <- glmLRT(fit, contrast = contrast)
    }
    most.de <- topTags(genes, n = nrow(genes$table))
    upreg <- most.de$table[most.de$table$logFC > 0, ]
    upreg <- upreg[order(upreg$logFC, decreasing = TRUE), ]
    downreg <- most.de$table[most.de$table$logFC < 0, ]
    downreg <- downreg[order(downreg$logFC), ]

    return(list(
        de = genes,
        genes = most.de$table,
        up = upreg,
        down = downreg
    ))
}


add_de <- function(exp, contrast, name, title, lfc = 0) {
    if (is.null(exp$dea)) {
        exp$dea <- list()
    }
    exp$dea[[name]] <- compute_de_genes(exp, contrast, lfc)
    de.obj <- compute_de_genes2(exp, contrast, lfc)
    return(exp)
}


plot_volcano <- function(
    exp,
    name,
    up.name = "",
    highlight.genes = NULL,
    alpha = 0.6,
    pt.size = 3
) {
    genes <- exp$dea[[name]]$genes
    up <- exp$dea[[name]]$up
    down <- exp$dea[[name]]$down

    upregs <- genes$FDR < .05 & genes$logFC > 0
    downregs <- genes$FDR < .05 & genes$logFC < 0
    up.legend <- paste0("Up in ", up.name, " (", sum(upregs), ")")
    down.legend <- paste0("Down in ", up.name, " (", sum(downregs), ")")

    genes$significance <- "Not significant"
    genes$significance[upregs] <- up.legend
    genes$significance[downregs] <- down.legend
    genes$significance <- factor(
        genes$significance,
        levels = c(
            "Not significant",
            down.legend,
            up.legend
        )
    )

    if (is.null(highlight.genes)) {
        top.genes <- data.frame(rbind(
            genes[rownames(up)[1:10], ],
            genes[rownames(down)[1:10], ]
        ))
    } else {
        top.genes <- genes[genes$external_gene_name %in% highlight.genes, ]
    }

    plot <- ggplot(genes, aes(x = logFC, y = -log10(FDR))) +
        geom_point(
            aes(color = significance),
            alpha = alpha,
            shape = 20,
            size = pt.size
        ) +
        theme_minimal() +
        scale_color_manual(
            values = c("#ddccbb", "#306090", "#961e2e"),
            labels = c(
                paste0("Not significant (", sum(genes$significance == "Not significant"), ")"),
                down.legend,
                up.legend
            )
        ) +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 10),
            legend.title = element_blank(),
            plot.title = element_text(size = 16)
        ) +
        geom_text_repel(
            data = top.genes,
            aes(label = external_gene_name),
            size = 5,
            box.padding = unit(1, "lines"),
            point.padding = unit(0, "lines"),
            min.segment.length = 0
        )
    return(plot)
}


run_de <- function(exp, contrast, name, title, highlight.genes = c(), lfc = 0) {
    dir.create(name, showWarnings = FALSE)
    de.obj <- compute_de_genes(exp, contrast, highlight.genes, lfc)

    pdf(paste(name, paste0(name, "-volcano.pdf"), sep = "/"), width = 9, height = 6)
    plot(de.obj$vol.plot + ggtitle(title))
    dev.off()

    write.csv(
        de.obj$genes[order(de.obj$genes$logFC, decreasing = TRUE), ],
        paste(name, paste0(name, "-list-full.csv"), sep = "/")
    )

    fdr.genes <- de.obj$genes[de.obj$genes$FDR < .05, ]
    fdr.genes <- fdr.genes[order(fdr.genes$logFC, decreasing = TRUE), ]
    write.csv(fdr.genes, paste(name, paste0(name, "-list.csv"), sep = "/"))

    write.table(
        na.omit(fdr.genes$external_gene_name[fdr.genes$logFC > 0]),
        paste(name, paste0(name, "-up-genes.txt"), sep = "/"),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

    write.table(
        na.omit(fdr.genes$external_gene_name[fdr.genes$logFC < 0]),
        paste(name, paste0(name, "-down-genes.txt"), sep = "/"),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

    write.table(
        na.omit(de.obj$genes$external_gene_name),
        paste(name, "background.txt", sep = "/"),
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

    # entrez <- mapIds(
    #   org.Mm.eg.db,
    #   keys = rownames(de.obj$de),
    #   keytype = "ENSEMBL",
    #   column = "ENTREZID"
    # )
    # entrez[is.na(entrez)] <- "NA"
    # rownames(de.obj$de) <- make.unique(entrez)
    # go <- goana(de.obj$de, species = "Mm")
    # write.csv(
    #   topGO(go, ontology = "BP", number = 50),
    #   paste(name, paste0(name, "-goana.csv"), sep = "/")
    # )
    return(de.obj)
}


run_enrichment <- function(contrast, pathways = "../../docs/mouse-isr.gmt") {
    isr.paths <- gmtPathways(pathways)
    # isr.paths <- examplePathways
    for (path in names(isr.paths)) {
      isr.paths[[path]] <- unique(isr.paths[[path]])
    }
    #ups <- contrast$genes[contrast$genes$logFC > 0, ]
    #downs <- contrast$genes[contrast$genes$logFC < 0,]
    #all <- list(up = ups, down = downs)

    result <- list()
    #for (dir in names(all)) {
    genes <- contrast$genes
    metric <- genes$LR * sign(genes$logFC)
    #ranks <- genes[order(genes$logFC, decreasing = TRUE), "logFC", drop = FALSE]
    rranks <- metric[order(metric, decreasing = TRUE)]
    #rranks <- unlist(ranks)
    names(rranks) <- rownames(genes)[order(metric, decreasing = TRUE)]
    # names(rranks) <- as.character(gene.anno$entrezgene_id[match(rownames(ranks), gene.anno$ensembl_gene_id)])
    rranks <- rranks[!is.na(names(rranks))]
    enr <- fgsea(
        pathways = isr.paths,
        stats = rranks,
        nperm = 10000
    )
    pval <- enr[enr$pathway == "isr-mouse", "padj"]
    plot <- plotEnrichment(isr.paths$`isr-mouse`, rranks) +
        ggtitle(paste0("ISR from Calico, padj=", pval))
    result[["table"]] <- enr
    result[["plot"]] <- plot
    #}

    return(result)
}


plot_kmeans <- function(exp, groups, colnames = NULL) {
    data <- exp$kmeans_data
    if (!is.null(colnames)) {
        colnames(data) <- colnames
    }

    clusters <- length(unique(exp$kmeans_cluster))
    mycol <- colorpanel(1000, "blue", "white", "red")
    df <- data.frame(group = factor(exp$kmeans_cluster))
    rownames(df) <- names(exp$kmeans_cluster)
    qq <- pheatmap(
        data,
        col = mycol,
        border_color = NA,
        cluster_rows = FALSE,
        annotation_row = df,
        show_rownames = FALSE,
        gaps_row = cumsum(table(sort(exp$kmeans_cluster)))[1:clusters],
        cutree_cols = groups,
        clustering_distance_cols = exp$kmeans_dist,
        clustering_method = "ward.D2"
    )
    return(qq)
}
