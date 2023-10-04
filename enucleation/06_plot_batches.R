
library(tidyverse)
library(SingleCellExperiment)
source("dataset.R")

main = function() {
    plot_batches()
    make_umap_with_ref()
    plot_batches_with_ref()
}

plot_batches = function() {
    sce = load_current()
    sce = sce[, sample.int(ncol(sce))] # randomize order of cells
    sce$condition = ifelse(grepl("1L", sce$sample) | grepl("3L", sce$sample), "CTL", "ENU")
    
    ## H1 level
    reducedDim(sce, "UMAP") = reducedDim(sce,"UMAP_H1")
    p1 = plot_raster_umap_sce(sce, "sample") + scale_color_manual(values = sample_cols()) + ggtitle("Merged")
    sample_p = lapply(sample_names(), function(my_sample) {
        sce$my_label = ifelse(sce$sample == my_sample, my_sample, "other")
        my_cols = setNames(c(sample_cols()[my_sample], gray(0.95, 0.05)), c(my_sample, "other"))
        return(plot_raster_umap_sce(sce, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    })
    p = patchwork::wrap_plots(c(list(p1), sample_p), ncol=3) & theme_void() & guides(col = "none")
    ggsave("fig/batches_h1.pdf", p)
    p = plot_raster_umap_sce(sce, "condition") + theme_void() + scale_color_manual(values = condition_cols())
    ggsave("fig/h1_condition.pdf", p + guides(col="none"))

    ## H2 level
    sce = sce[, sce$h1 == "Glutamatergic"]
    reducedDim(sce, "UMAP") = reducedDim(sce,"UMAP_H2")
    p1 = plot_raster_umap_sce(sce, "sample") + scale_color_manual(values = sample_cols()) + ggtitle("Merged")
    sample_p = lapply(sample_names(), function(my_sample) {
        sce$my_label = ifelse(sce$sample == my_sample, my_sample, "other")
        my_cols = setNames(c(sample_cols()[my_sample], gray(0.95, 0.05)), c(my_sample, "other"))
        return(plot_raster_umap_sce(sce, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    })
    p = patchwork::wrap_plots(c(list(p1), sample_p), ncol=3) & theme_void() & guides(col = "none")
    ggsave("fig/batches_h2.pdf", p)
    p = plot_raster_umap_sce(sce, "condition") + theme_void() + scale_color_manual(values = condition_cols())
    ggsave("fig/h2_condition.pdf", p + guides(col="none"))
    
    ## H3 level
    sce = sce[, sce$h2 == "L2/3 IT"]
    reducedDim(sce, "UMAP") = reducedDim(sce,"UMAP_H3")
    p = plot_raster_umap_sce(sce, "condition") + theme_void() + scale_color_manual(values = condition_cols())
    ggsave("fig/h3_condition.pdf", p + guides(col="none"))
    reducedDim(sce, "UMAP") = reducedDim(sce,"UMAP_H2")
    p = plot_raster_umap_sce(sce, "condition", alpha=0.5) + scale_color_manual(values = condition_cols()) + coord_cartesian(xlim = c(-4,0),c(-9,-5)) + theme_void()
    p
    ggsave("fig/h3_condition_v2.pdf", p + guides(col="none"))
}

load_current = function() {
    barseq = load_barseq()
    sample_cell = strsplit(colnames(barseq), split = ".", fixed = TRUE) # replace by regex for efficiency
    barseq$sample = sapply(sample_cell, "[", 1)

    umap = column_to_rownames(as.data.frame(read_csv("data/umap.csv.gz", show_col_types=FALSE)), "cell_id")
    reducedDim(barseq, "UMAP_H1") = umap[colnames(barseq), ]

    umap = as.matrix(column_to_rownames(as.data.frame(read_csv("analysis/glu/umap.csv.gz", show_col_types=FALSE)), "cell_id"))
    umap_full = matrix(NA, ncol = 2, nrow = ncol(barseq), dimnames = list(colnames(barseq), c("UMAP1", "UMAP2")))
    umap_full[rownames(umap),] = umap
    reducedDim(barseq, "UMAP_H2") = as.data.frame(umap_full)

    umap = as.matrix(column_to_rownames(as.data.frame(read_csv("analysis/subglu/subumap.csv.gz", show_col_types=FALSE)), "cell_id"))
    umap_full = matrix(NA, ncol = 2, nrow = ncol(barseq), dimnames = list(colnames(barseq), c("UMAP1", "UMAP2")))
    umap_full[rownames(umap),] = umap
    reducedDim(barseq, "UMAP_H3") = as.data.frame(umap_full)

    clusters = load_labels()
    barseq$h1 = clusters$h1
    barseq$h2 = clusters$h2
    barseq$h3 = clusters$h3

    return(barseq)
}

plot_raster_umap_sce = function(barseq, color_by, text_by=NULL, point_size=0.1, text_size=3, raster.dpi=150, alpha=1) {
    to_plot = as_tibble(reducedDim(barseq, "UMAP"))
    colnames(to_plot) = c("UMAP_1", "UMAP_2")
    if (color_by %in% rownames(barseq)) {
        to_plot[[color_by]] = logcounts(barseq)[color_by,]
    } else {
        to_plot[[color_by]] = barseq[[color_by]]
    }
    if (!is.null(text_by)) {
        to_plot[[text_by]] = barseq[[text_by]]
    }    
    plot_raster_umap(to_plot, color_by, text_by, point_size, text_size, raster.dpi, alpha)
}

plot_raster_umap = function(to_plot, color_by, text_by=NULL, point_size=0.1, text_size=3, raster.dpi=150, alpha=1) {
    # trick to make dot size as small as possible: https://stackoverflow.com/questions/34638902/point-size-in-ggplot-2-0-0
    result=to_plot %>%
        ggplot(aes(x=UMAP_1, y=UMAP_2, col=.data[[color_by]])) +
        ggrastr::geom_point_rast(size=point_size, shape=16, stroke=0, raster.dpi = raster.dpi, alpha=alpha) +
        theme_classic()
    if (is.numeric(to_plot[[color_by]])) {
        result = result + scale_color_viridis_c()
    }
    if (!is.null(text_by)) {
        to_plot = to_plot %>%
            group_by(.data[[text_by]]) %>%
            summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        result = result +
            geom_text(data=to_plot, aes(label=.data[[text_by]]), size=text_size, color="black")
    }
    return(result)
}

load_reference = function(label_file = "data/ref_labels.csv") {
    barseq = load_ref()
    barseq$sample = "reference"
    clusters = read_csv(label_file)
    barseq = barseq[, as.character(clusters$sample)]
    barseq$h1 = clusters$class
    barseq$h2 = clusters$subclass
    barseq$h3 = clusters$cluster
    colLabels(barseq) = barseq$cluster_id
    return(barseq)
}

make_umap_with_ref = function(n_pca = 30, k_umap = 15) {
    dir.create("analysis/ref", FALSE, TRUE)
    f1 = "analysis/ref/pca_h1.rds"
    f2 = "analysis/ref/pca_h2.rds"
    if (!file.exists(f1) | !file.exists(f2)) {
        current = load_current()
        current = current[, current$sample %in% c("D078_1L", "D078_2L")]
        ref = load_reference()
        sce = MetaNeighbor::mergeSCE(list(current = current, ref = ref))
        rm(current); rm(ref); gc()
        logcounts(sce) = log1p(convert_to_cpm(counts(sce), 10))
        if (!file.exists(f1)) {
            pca = scater::calculatePCA(sce, ncomponents=n_pca, exprs_values="logcounts")
            saveRDS(pca, f1)
            rm(pca)
        }
        if (!file.exists(f2)) {
            sce = sce[, sce$h1 == "Glutamatergic"]
            pca = scater::calculatePCA(sce, ncomponents=n_pca, exprs_values="logcounts")
            saveRDS(pca, f2)
            rm(pca)
        }
        rm(sce); gc()
    }
    filename = "analysis/ref/umap_h1.csv.gz"
    if (!file.exists(filename)) {
        pca = readRDS("analysis/ref/pca_h1.rds")
        umap = scater::calculateUMAP(
            pca, transposed=TRUE, n_neighbors=k_umap,
            external_neighbors=TRUE, BNPARAM = BiocNeighbors::AnnoyParam()
        )
        colnames(umap) = c("UMAP1","UMAP2")
        write_csv(as_tibble(umap, rownames = "cell_id"), filename)
        rm(umap); rm(pca); gc()
    }
    filename = "analysis/ref/umap_h2.csv.gz"
    if (!file.exists(filename)) {
        pca = readRDS("analysis/ref/pca_h2.rds")
        umap = scater::calculateUMAP(
            pca, transposed=TRUE, n_neighbors=k_umap,
            external_neighbors=TRUE, BNPARAM = BiocNeighbors::AnnoyParam()
        )
        colnames(umap) = c("UMAP1","UMAP2")
        write_csv(as_tibble(umap, rownames = "cell_id"), filename)
        rm(umap); rm(pca); gc()
    }
}

plot_batches_with_ref = function() {    
    label = load_labels()
    h1 = setNames(label$h1, label$cell_id)
    h2 = setNames(label$h2, label$cell_id)
    
    ## H1 level
    umap = read_csv("analysis/ref/umap_h1.csv.gz", show_col_types = FALSE)
    umap = rename(umap, UMAP1 = "UMAP_1", UMAP2 = "UMAP_2")
    umap = umap[sample.int(nrow(umap)),]
    umap$sample = substr(umap$cell_id, 1, 7)
    umap$sample = ifelse(umap$sample %in% sample_names(), umap$sample, "REF")
    umap$cell_type = h1[umap$cell_id]
    all_samples = names(sample_cols())[names(sample_cols()) %in% unique(umap$sample)]
    my_sample_cols = sample_cols()
#    my_sample_cols["REF"] = gray(0.5)
    p1 = plot_raster_umap(umap, "sample", text_by = "cell_type") + scale_color_manual(values = my_sample_cols) + ggtitle("Merged")
    sample_p = lapply(all_samples, function(my_sample) {
        umap$my_label = ifelse(umap$sample == my_sample, my_sample, "other")
        my_cols = setNames(c(my_sample_cols[my_sample], gray(0.95, 0.05)), c(my_sample, "other"))
        return(plot_raster_umap(umap, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    })
    p = patchwork::wrap_plots(c(list(p1), sample_p), ncol=2) & theme_void() & guides(col = "none")
    ggsave("fig/batches_ref_h1.pdf", p)

    ## H2 level
    umap = read_csv("analysis/ref/umap_h2.csv.gz", show_col_types = FALSE)
    umap = rename(umap, UMAP1 = "UMAP_1", UMAP2 = "UMAP_2")
    umap = umap[sample.int(nrow(umap)),]
    umap$sample = substr(umap$cell_id, 1, 7)
    umap$sample = ifelse(umap$sample %in% sample_names(), umap$sample, "REF")
    umap$cell_type = h2[umap$cell_id]
    all_samples = names(sample_cols())[names(sample_cols()) %in% unique(umap$sample)]
    my_sample_cols = sample_cols()
#    my_sample_cols["REF"] = gray(0.5)
    p1 = plot_raster_umap(umap, "sample", text_by = "cell_type") + scale_color_manual(values = my_sample_cols) + ggtitle("Merged")
    sample_p = lapply(all_samples, function(my_sample) {
        umap$my_label = ifelse(umap$sample == my_sample, my_sample, "other")
        my_cols = setNames(c(my_sample_cols[my_sample], gray(0.95, 0.05)), c(my_sample, "other"))
        return(plot_raster_umap(umap, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    })
    p = patchwork::wrap_plots(c(list(p1), sample_p), ncol=2) & theme_void() & guides(col = "none")
    ggsave("fig/batches_ref_h2.pdf", p)
}

if (sys.nframe() == 0) {
    main()
}