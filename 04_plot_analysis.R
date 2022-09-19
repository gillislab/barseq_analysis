
library(tidyverse)
library(SingleCellExperiment)
library(rasterpdf)
source("dataset.R")


main = function() {
    #plot_analysis("analysis/whole", "figs/whole")
    #plot_analysis("analysis/glu", "figs/glu")
    #plot_analysis("analysis/gaba", "figs/gaba")
    plot_subglu()
}

plot_analysis = function(input_dir, output_dir, max_cells = 2000000) {
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)    
    barseq = load_dataset(input_dir, max_cells)
    plot_umap(barseq, output_dir)
    plot_gradients(barseq, output_dir)
    plot_marker_scores(barseq, output_dir)
    plot_slice_summary(barseq, output_dir)
    plot_slices(barseq, output_dir)
    plot_cell_types(barseq, output_dir)
}

plot_subglu = function(input_dir = "analysis/subglu", output_dir = "figs/subglu") {
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)
    barseq = load_dataset(input_dir)
    plot_umap(barseq, output_dir,1)
    plot_gradients(barseq, output_dir)

    ct_dir = file.path(output_dir, "cell_types")
    dir.create(ct_dir, showWarnings = FALSE, recursive=TRUE)    
    reducedDim(barseq, "UMAP") = as.data.frame(read_csv(file.path(input_dir, "subumap.csv")))
    subclass_id = read_csv("analysis/glu/cluster.csv")
    subclass_id = subclass_id[match(colnames(barseq), subclass_id$sample),]$label
    subclass_annotation = deframe(read_csv("analysis/glu/cluster_annotation.csv"))
    subclass_name = subclass_annotation[subclass_id]
    for (my_subclass in sort(unique(subclass_name))) {
        filename = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', my_subclass)))
        subdata = barseq
        subdata$label[subclass_name!=my_subclass] = NA
        plot_umap(subdata[, !is.na(subdata$label)], ct_dir, filename=paste0(filename, "_umap.pdf"))
        plot_gradients(subdata[, !is.na(subdata$label)], ct_dir, filename=paste0(filename, "_gradients.pdf"))
        plot_slice_summary(subdata, ct_dir, filename=paste0(filename, "_slice.pdf"))
    }
}

load_dataset = function(input_dir, max_cells=2000000) {
    barseq = load_barseq(normalization_factor = 10)
    barseq$slice = as.factor(barseq$slice)
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    barseq$n_genes = colSums(counts(barseq)>0)
    
    clusters = read_csv(file.path(input_dir, "cluster.csv"))
    cluster_annotation = deframe(read_csv(file.path(input_dir, "cluster_annotation.csv")))
    barseq = barseq[, as.character(clusters$sample)]
    barseq$cluster_id = as.factor(clusters$label)
    colLabels(barseq) = cluster_annotation[barseq$cluster_id]
    reducedDim(barseq, "UMAP") = as.data.frame(read_csv(file.path(input_dir, "umap.csv")))
    
    if (ncol(barseq) > max_cells) {
        set.seed(17)
        cells_to_keep = sample.int(ncol(barseq), max_cells)
        barseq = barseq[, cells_to_keep]
    }
    return(barseq)
}

plot_umap = function(barseq, output_dir, text_size=3, filename="preliminary_analysis.pdf") {
    point_size = 70 / sqrt(ncol(barseq))
    pdf(file.path(output_dir, filename))
    # try alpha = 0.05
    p = plot_raster_umap(barseq, color_by="label", text_by="subclass", text_size=text_size, point_size=point_size) +
        ggtitle("Preliminary clustering") +
        guides(col="none") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "n_genes", point_size=point_size) +
        ggtitle("Number of detected genes") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "Slc17a7", point_size=point_size) +
        ggtitle("Excitatory marker") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "Gad1", point_size=point_size) +
        ggtitle("Inhibitory marker") +
        coord_fixed()
    print(p)
    dev.off()
}

plot_raster_umap = function(barseq, color_by, text_by=NULL, point_size=0.1, text_size=3, raster.dpi=150, alpha=1) {
    to_plot = as_tibble(reducedDim(barseq, "UMAP"))
    colnames(to_plot) = c("UMAP_1", "UMAP_2")
    if (color_by %in% rownames(barseq)) {
        to_plot[[color_by]] = logcounts(barseq)[color_by,]
    } else {
        to_plot[[color_by]] = barseq[[color_by]]
    }
    # trick to make dot size as small as possible: https://stackoverflow.com/questions/34638902/point-size-in-ggplot-2-0-0
    result=to_plot %>%
        ggplot(aes(x=UMAP_1, y=UMAP_2, col=.data[[color_by]])) +
        ggrastr::geom_point_rast(size=point_size, shape=16, stroke=0, raster.dpi = raster.dpi, alpha=alpha) +
        theme_classic()
    if (is.numeric(to_plot[[color_by]])) {
        result = result + scale_color_viridis_c()
    }
    if (!is.null(text_by)) {
        to_plot[[text_by]] = barseq[[text_by]]
        to_plot = to_plot %>%
            group_by(.data[[text_by]]) %>%
            summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        result = result +
            geom_text(data=to_plot, aes(label=.data[[text_by]]), size=text_size, color="black")
    }
    return(result)
}

plot_gradients = function(barseq, output_dir, filename="gradients.pdf") {
    set.seed(17)
    point_size = 70 / sqrt(ncol(barseq))
    barseq$depth = sqrt(barseq$depth_x^2 + barseq$depth_y^2)
    barseq$slice = as.numeric(barseq$slice)
    barseq = barseq[, sample.int(ncol(barseq))]
    pdf(file.path(output_dir, filename))
    p = plot_raster_umap(barseq, color_by = "slice", point_size=point_size) +
        ggtitle("Anterior-Posterior") +
        scale_color_distiller(palette = "RdYlBu") +
        labs(col = "Slice") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "aligned_x", point_size=point_size) +
        ggtitle("Medial-Lateral") +
        scale_color_distiller(palette = "RdYlBu") +
        labs(col = "Medial->Lateral") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "aligned_y", point_size=point_size) +
        ggtitle("Dorso-ventral") +
        scale_color_distiller(palette = "RdYlBu") +
        labs(col = "Ventral->Dorsal") +
        coord_fixed()
    print(p)
    dev.off()
}

plot_marker_scores = function(barseq, output_dir) {
    marker_scores = load_marker_scores(barseq)
    MetaMarkers::plot_marker_scores(
        marker_scores, reducedDim(barseq, "UMAP"),
        normalize_scores=FALSE, rasterize_umap=TRUE, point_size=0.1
    ) +
        theme_classic() +
        ggtitle("Marker score")
    ggsave(file.path(output_dir, "marker_scores.pdf"), width = 15, height=15)
}

load_marker_scores = function(barseq) {
    marker_scores = read_csv("analysis/marker_scores.csv") %>%
        column_to_rownames("sample_id") %>%
        as.matrix()
    marker_scores = t(marker_scores) / colMaxs(marker_scores)
    marker_scores = marker_scores[,as.character(colnames(barseq))]
    return(marker_scores)
}

plot_slice_summary = function(barseq, output_dir, filename="all_slices.pdf") {    
    n_slices = max(as.numeric(unfactor(barseq$slice)))
    pdf(file.path(output_dir, filename))
    for (i in 1:n_slices) {
        print(plot_slice(barseq, i))
    }
    dev.off()
}

plot_slices = function(barseq, output_dir) {
    slice_dir = file.path(output_dir, "slices")
    dir.create(slice_dir, showWarnings = FALSE, recursive=TRUE)
    n_slices = max(as.numeric(unfactor(barseq$slice)))
    unique_labels = levels(as.factor(barseq$label))
    for (i in 1:n_slices) {
        pdf(file.path(slice_dir, paste0("slice_", i, ".pdf")))
        print(plot_slice(barseq, i))
        for (l in unique_labels) {
            print(plot_slice_ct(barseq, i, l))
        }
        dev.off()
    }
}

plot_cell_types = function(barseq, output_dir) {
    ct_dir = file.path(output_dir, "cell_types")
    dir.create(ct_dir, showWarnings = FALSE, recursive=TRUE)
    n_slices = max(as.numeric(unfactor(barseq$slice)))
    for (l in unique(barseq$label)) {
        file_name = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', l)))
        pdf(file.path(ct_dir, paste0(file_name, ".pdf")))
        for (i in 1:n_slices) {
            print(plot_slice_ct(barseq, i, l))
        }
        dev.off()
    }
}

plot_slice = function(barseq, slice_id, raster.dpi=150) {
    as.data.frame(colData(barseq)) %>%
        mutate(label = as.factor(label)) %>%
        filter(slice == slice_id) %>%
        ggplot(aes(x = aligned_x, y = aligned_y, col = label)) +
        ggrastr::geom_point_rast(size=0.5, shape=16, stroke=0, raster.dpi = raster.dpi) +
        theme_classic() +
        labs(x=NULL, y=NULL, col=NULL, title=paste0("Slice ", slice_id)) +
        coord_equal() +
        lims(x = range(barseq$aligned_x), y = range(barseq$aligned_y)) +
        theme(legend.text=element_text(size=6), legend.position = "bottom") +
        scale_color_discrete(na.value="gray90", drop=FALSE)
}

plot_slice_ct = function(barseq, slice_id, ct) {
    all_labels = as.factor(barseq$label)
    my_cols = make_cluster_cols(all_labels)
    full_plot = plot_slice(barseq, slice_id)
    my_cols[names(my_cols) != ct] = "gray90"
    my_plot = full_plot +
        ggtitle(paste0("Slice ", slice_id, ", Cluster ", ct)) +
        scale_color_manual(values = my_cols, drop=FALSE)
    return(my_plot)
}

make_cluster_cols = function(cluster_names) {
    cluster_names = levels(cluster_names)
    result = scales::hue_pal()(length(cluster_names))
    names(result) = cluster_names
    return(result)
}

if (sys.nframe() == 0) {
    main()
}
