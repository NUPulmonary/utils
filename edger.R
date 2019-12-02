library(edgeR)
library(ggplot2)
library(biomaRt)
library(ggrepel)


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

    anno <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values= rownames(counts),
                  mart = mart)
    return(anno)
}

run_pca <- function(data, anno, design = NULL, samples = NULL) {
    if (!is.null(samples)) {
        dge <- DGEList(data, group = samples$group)
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
    df$group <- samples$group[match(rownames(df), samples$sample_name)]
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
        mds.fig = plot_mds(dge, samples$group),
        pca = dge.pca,
        pca.fig = pca.fig
        #eig.fig = eig.fig
    ))
}


subset_experiment <- function(samples, counts, exclude = c()) {
    columns <- as.character(
        samples$sample_name[samples$group != "Uni RNA" & !(samples$sample_id %in% exclude)]
    )
    samples <- samples[samples$sample_name %in% columns, ]
    design <- model.matrix(~0 + samples$group)
    colnames(design) <- levels(samples$group)
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
    rownames(df) <- data.samples$sample_name
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
