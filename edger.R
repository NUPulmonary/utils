library(edgeR)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(fgsea)


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
    group_pattern = "S\\d\\d\\_(.+)\\_R\\d+"
) {
    samples <- read.delim(samples_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[, 1:2]
    colnames(samples) <- c("sample_id", "sample_name")
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
        colnames(counts) <- samples$sample_name
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


plot_mds <- function(data, group, dim.plot = c(1, 2)) {
    data <- plotMDS(data, plot = FALSE, dim.plot = dim.plot)
    data <- data.frame(x = data$x, y = data$y, group = group)
    return(ggplot(data, aes(x = x, y = y, col = group)) + geom_point())
}


get_anno <- function(counts) {
    mart <- useEnsembl(biomart = "ensembl",
                       dataset = "mmusculus_gene_ensembl",
                       mirror = "uswest") # options: useast, uswest, asia

    anno <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id'),
                  filters = 'ensembl_gene_id',
                  values = rownames(counts),
                  mart = mart)
    return(anno)
}

run_pca <- function(data, anno, design = NULL, samples = NULL, group_column = "group") {
    if (!is.null(samples)) {
        dge <- DGEList(data, group = samples[[group_column]])
    } else {
        dge <- DGEList(data)
    }
    dge$genes <- data.frame(
        external_gene_name = anno[match(rownames(dge), anno$ensembl_gene_id), "external_gene_name"]
    )
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = design)

    dge.pca <- prcomp(cpm(dge, log = TRUE), center = TRUE, scale = TRUE)

    df <- data.frame(dge.pca$rotation[, 1:2])
    df$group <- samples[[group_column]][match(rownames(df), samples$sample_name)]

    # TODO redo
    df <- cbind(
        df,
        t(data.frame(lapply(strsplit(as.character(
            rownames(df)
        ), '_'), `[`, 2:5)))
    )
    colnames(df) <- c(colnames(df)[1:3], "Age", "Bleo", "CellType", "Treatment")

    pca.fig <- ggplot(df) +
        geom_point(aes(PC1, PC2, col = group)) +
        ggtitle("PCA of all samples")

    #eig.fig <- fviz_eig(dge.pca, main = "PC variance explained")
    return(list(
        dge = dge,
        mds.fig = plot_mds(dge, samples[[group_column]]),
        pca = dge.pca,
        pca.fig = pca.fig
        #eig.fig = eig.fig
    ))
}


subset_experiment <- function(samples, counts, exclude = c(), group_column = "group") {
    group <- samples[[group_column]]
    columns <- as.character(
        samples$sample_name[group != "Uni RNA" & !(samples$sample_id %in% exclude)]
    )
    samples <- samples[samples$sample_name %in% columns, ]
    group <- samples[[group_column]]
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    design <- design[, colnames(design)[colSums(design) > 0]]
    rownames(design) <- samples$sample_name
    return(list(
        samples = samples,
        counts = counts[, columns],
        design = design
    ))
}


plot_gene <- function(data.obj, data.samples, gene.anno, gene) {
    gene.id <- gene.anno$ensembl_gene_id[gene.anno$external_gene_name == gene]
    df <- data.frame(group=data.samples$group)
    rownames(df) <- data.samples$sample_id
    df[[gene]] <- data.obj$dge$counts[gene.id, ]
    return(ggplot(df, aes_string(x = "group", y = gene)) +
        theme_minimal() +
        geom_boxplot() +
        geom_dotplot(binaxis = 'y', stackdir = 'center') +
        geom_label_repel(
            aes(label = rownames(df)),
            na.rm = TRUE,
            nudge_x = .3,
            size = 3
        )
    )
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


compute_de_genes <- function(dge, design, contrast, highlight.genes) {
    up.name <- rownames(contrast)[contrast > 0]
    fit <- glmFit(dge, design)
    genes <- glmLRT(fit, contrast = contrast)
    # TODO: explore
    # genes <- glmTreat(fit, contrast = contrast)#, lfc = log2(1.5))
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


run_de <- function(obj, design, contrast, name, title, highlight.genes = c()) {
    dir.create(name, showWarnings = FALSE)
    de.obj <- compute_de_genes(obj$dge, design, contrast, highlight.genes)

    pdf(paste(name, paste0(name, "-volcano.pdf"), sep = "/"), width = 9, height = 6)
    plot(de.obj$vol.plot + ggtitle(title))
    dev.off()

    fdr.genes <- de.obj$genes[de.obj$genes$FDR < .05, ]
    fdr.genes <- fdr.genes[order(fdr.genes$logFC, decreasing = TRUE), ]
    write.csv(fdr.genes, paste(name, paste0(name, "-list.csv"), sep = "/"))

    write.table(
        fdr.genes$external_gene_name[fdr.genes$logFC > 0], 
        paste(name, paste0(name, "-up-genes.txt"), sep = "/"), 
        row.names = FALSE, 
        col.names = FALSE, 
        quote = FALSE
    )

    write.table(
        fdr.genes$external_gene_name[fdr.genes$logFC < 0], 
        paste(name, paste0(name, "-down-genes.txt"), sep = "/"), 
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


run_enrichment <- function(contrast) {
    isr.paths <- gmtPathways("../../docs/mouse-isr.gmt")
    # isr.paths <- examplePathways
    for (path in names(isr.paths)) {
      isr.paths[[path]] <- unique(isr.paths[[path]])
    }
    ranks <- contrast$genes[order(contrast$genes$LR), "LR", drop = FALSE]
    rranks <- unlist(ranks)
    names(rranks) <- rownames(ranks)
    # names(rranks) <- as.character(gene.anno$entrezgene_id[match(rownames(ranks), gene.anno$ensembl_gene_id)])
    rranks <- rranks[!is.na(names(rranks))]
    result <- fgsea(
        pathways = isr.paths,
        stats = rranks,
        nperm = 10000
    )
    pval <- result[result$pathway == "isr-mouse-extended", "padj"]
    plot <- plotEnrichment(isr.paths$`isr-mouse-extended`, rranks) + 
        ggtitle(paste0("ISR extended from Calico, padj=", pval))
    return(list(
        table = result,
        plot = plot
    ))
}