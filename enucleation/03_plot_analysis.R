
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

plot_subglu = function(input_dir = "analysis/subglu", output_dir = "fig/subglu") {
    message("H3 level plots...")
    set.seed(17)    
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)

    # load and preformat data
    message("Loading data...")
    barseq = load_dataset(input_dir, max_cells = Inf)
    barseq$slice_number = as.numeric(as.character(barseq$slice))
    barseq$slice_id = paste(barseq$sample, str_pad(barseq$slice, 2, pad="0"), sep="|")
    barseq$original_label = barseq$label
    barseq = barseq[, sample.int(ncol(barseq))] # randomize order for plotting
    subclass_id = deframe(read_csv("analysis/glu/cluster.csv.gz", show_col_types=FALSE))
    barseq$subclass_id = as.factor(subclass_id[colnames(barseq)])
    subclass_annotation = deframe(read_csv("analysis/glu/cluster_annotation.csv", show_col_types=FALSE))
    barseq$subclass_name = subclass_annotation[as.character(barseq$subclass_id)]
    
    # global plots
    if (FALSE) {
    message("Global plots...")
    plot_umap(barseq, output_dir, 1)
    plot_gradients(barseq, output_dir)
    plot_ct_prop(barseq, output_dir)
    }
    
    # make individual plots for each subclass
    if (FALSE) {
    ct_dir = file.path(output_dir, "cell_types")
    dir.create(ct_dir, showWarnings = FALSE, recursive=TRUE)    
    umap = column_to_rownames(as.data.frame(read_csv(file.path(input_dir, "subumap.csv.gz"), show_col_types=FALSE)), "cell_id")
    reducedDim(barseq, "UMAP") = umap[colnames(barseq),]
    message("Cell-type-specific plots...")
    for (my_subclass in sort(unique(barseq$subclass_name))) {
        message(my_subclass)
        filename = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', my_subclass)))
        barseq$label = barseq$original_label
        barseq$label[barseq$subclass_name!=my_subclass] = NA
        subclass_data = barseq[, !is.na(barseq$label)]
        plot_umap(subclass_data, ct_dir, filename=paste0(filename, "_umap.pdf"))
        pdf(file.path(ct_dir, paste0(filename, "_batch.pdf")))
        for (my_sample in sort(unique(subclass_data$sample))) {
            subclass_data$my_label = ifelse(subclass_data$sample == my_sample, my_sample, "other")
            my_cols = setNames(c("red", gray(0.9, 0.1)), c(my_sample, "other"))
            p = plot_raster_umap(subclass_data, "my_label") +
                scale_color_manual(values = my_cols) +
                ggtitle(my_sample) +
                coord_fixed()
            print(p)
        }
        dev.off()
        plot_gradients(subclass_data, ct_dir, filename=paste0(filename, "_gradients.pdf"))
        plot_slice_summary(barseq, ct_dir, filename=paste0(filename, "_slice.pdf"))
    }
    }
    ct_dir = file.path(output_dir, "cell_types_manual_h2_h3")
    dir.create(ct_dir, showWarnings = FALSE, recursive=TRUE)    
    if ("subclass" %in% names(colData(barseq))) {
        for (my_subclass in sort(unique(barseq$subclass))) {
            message(my_subclass)
            filename = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', my_subclass)))
            barseq$label = barseq$original_label
            barseq$label[barseq$subclass!=my_subclass] = NA
            plot_slice_summary(barseq, ct_dir, filename=paste0(filename, "_slice.pdf"))
        }
    } 
}

load_dataset = function(input_dir, max_cells=2000000) {
    barseq = load_barseq()
    sample_cell = strsplit(colnames(barseq), split = ".", fixed = TRUE) # replace by regex for efficiency
    barseq$sample = sapply(sample_cell, "[", 1)
    barseq$slice = as.factor(sapply(strsplit(sapply(sample_cell, "[", 2), split = "_", fixed = TRUE),"[", 1))
    barseq$n_genes = colSums(logcounts(barseq)>0)
    
    # load annotations and UMAP
    clusters = read_csv(file.path(input_dir, "cluster.csv.gz"), show_col_types=FALSE)
    cluster_annotation = read_csv(file.path(input_dir, "cluster_annotation.csv"), show_col_types=FALSE)
    clusters = left_join(clusters, cluster_annotation, by=c("label"="cluster_id"))
    # keep only clusters not annotated as "NA"
    barseq = barseq[, as.character(clusters$sample)]
    barseq$cluster_id = as.factor(clusters$label)
    colLabels(barseq) = clusters$cluster_name
    barseq$subclass = clusters$subclass
    umap = column_to_rownames(as.data.frame(read_csv(file.path(input_dir, "umap.csv.gz"), show_col_types=FALSE)), "cell_id")
    reducedDim(barseq, "UMAP") = umap[colnames(barseq), ]
    
    # load cell positions
    cellpos = read_csv("data/cell_position.csv.gz", show_col_types=FALSE)
    cellpos$cell_id = paste(cellpos$sample, cellpos$cell_id, sep=".")
    cellpos = cellpos[match(colnames(barseq), cellpos$cell_id),]
    barseq$aligned_x = cellpos$pos_x
    barseq$aligned_y = cellpos$pos_y
    
    if (ncol(barseq) > max_cells) {
        cells_to_keep = sample.int(ncol(barseq), max_cells)
        barseq = barseq[, cells_to_keep]
    }
    return(barseq)
}

plot_umap = function(barseq, output_dir, text_size=3, filename="preliminary_analysis.pdf") {
    point_size = 70 / sqrt(ncol(barseq))
    alpha = min(1, 100000/ncol(barseq))
    pdf(file.path(output_dir, filename))
    # try alpha = 0.05
    color_by = ifelse("subclass" %in% names(colData(barseq)), "subclass", "label")
    if ("subclass" %in% names(colData(barseq))) {
        p = plot_raster_umap(barseq, color_by=color_by, text_by="subclass", text_size=text_size, point_size=point_size, alpha=alpha) +
            ggtitle("Preliminary clustering") +
            guides(col="none") +
            coord_fixed()
        print(p)
    }
    p = plot_raster_umap(barseq, color_by=color_by, text_by="label", text_size=text_size, point_size=point_size, alpha=alpha) +
        ggtitle("Preliminary clustering") +
        guides(col="none") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "n_genes", point_size=point_size, alpha=alpha) +
        ggtitle("Number of detected genes") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "Slc17a7", point_size=point_size, alpha=alpha) +
        ggtitle("Excitatory marker") +
        coord_fixed()
    print(p)
    p = plot_raster_umap(barseq, color_by = "Gad1", point_size=point_size, alpha=alpha) +
        ggtitle("Inhibitory marker") +
        coord_fixed()
    print(p)
    stats = tibble(cell_type = barseq$label, n_genes=barseq$n_genes, gaba = logcounts(barseq)["Gad1",], glu = logcounts(barseq)["Slc17a7",])
    p = ggplot(stats, aes(cell_type, n_genes)) +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL, y="# genes detected")
    print(p)
    p = ggplot(stats, aes(cell_type, glu)) +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL, y="Slc17a7 expression")
    print(p)
    p = stats %>%
        group_by(cell_type) %>%
        summarize(glu = mean(glu>0)) %>%
        ggplot(aes(cell_type, glu)) +
        geom_col() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL, y="Fraction cells expressing Slc17a7") +
        geom_hline(yintercept = 0.95, linetype = "dashed")
    print(p)
    p = ggplot(stats, aes(cell_type, gaba)) +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL, y="Gad1 expression")
    print(p)
    p = stats %>%
        group_by(cell_type) %>%
        summarize(gaba = mean(gaba>0)) %>%
        ggplot(aes(cell_type, gaba)) +
        geom_col() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL, y="Fraction cells expressing Gad1") +
        geom_hline(yintercept = 0.95, linetype = "dashed")
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

plot_ct_prop = function(barseq, output_dir) {
    full_stats = as_tibble(colData(barseq)) %>%
        mutate(condition = ifelse(grepl("1L", sample) | grepl("3L", sample), "CTL", "ENU")) %>%
        group_by(label, subclass, condition, sample) %>%
        tally() %>%
        group_by(subclass, condition, sample) %>%
        mutate(f=n/sum(n))
    mean_stats = full_stats %>%
        group_by(label, condition, subclass) %>%
        summarize(me = mean(f))
    p = ggplot(mean_stats, aes(x=label,fill=condition)) +
        geom_col(aes(y=me), position = "dodge") +
        geom_point(data=full_stats, aes(y=f), position = position_dodge(0.8)) +
        facet_wrap(~ subclass, scales = "free") +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(y="Proportion")
    ggsave(file.path(output_dir, "ct_prop.pdf"), width = 20, height = 20)
    full_stats = full_stats %>%
        mutate(sample = substr(sample, 1, 4)) %>%
        group_by(label, sample, subclass) %>%
        summarize(me = f[condition=="ENU"] - f[condition=="CTL"])
    p = ggplot(full_stats, aes(x=label,y=me,fill=subclass)) +
        geom_boxplot(aes(y=me), show.legend = FALSE) +
        geom_jitter(show.legend=FALSE) +
        geom_hline(yintercept = 0, linetype="dashed") +
        facet_wrap(~ subclass, scales = "free") +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(y="Proportion difference within litter (ENO - CTL)")
    ggsave(file.path(output_dir, "ct_prop_paired.pdf"), width = 40/3, height = 40/3)
}
    
    
plot_gradients = function(barseq, output_dir, filename="gradients.pdf", slice_slot="slice_number") {
    set.seed(17)
    point_size = 70 / sqrt(ncol(barseq))
    all_samples = sort(unique(barseq$sample))
    pdf(file.path(output_dir, filename))
    p = lapply(all_samples, function(my_sample) {
        plot_raster_umap(barseq[, barseq$sample == my_sample], color_by = slice_slot, point_size=point_size) +
            scale_color_distiller(palette = "RdYlBu") +
            ggtitle(my_sample) +
            coord_fixed() +
            theme_void() +
            guides(col = "none")
    })
    p = patchwork::wrap_plots(p)+patchwork::plot_annotation(title="Posterior->Anterior (blue->red)")
    print(p)
    p = lapply(all_samples, function(my_sample) {
        plot_raster_umap(barseq, color_by = "aligned_x", point_size=point_size) +
            scale_color_distiller(palette = "RdYlBu") +
            ggtitle(my_sample) +
            coord_fixed() +
            theme_void() +
            guides(col = "none")
    })
    p = patchwork::wrap_plots(p)+patchwork::plot_annotation(title="Medial->Lateral (blue->red)")
    print(p)
    p = lapply(all_samples, function(my_sample) {
        plot_raster_umap(barseq, color_by = "aligned_y", point_size=point_size) +
            scale_color_distiller(palette = "RdYlBu") +
            ggtitle(my_sample) +
            coord_fixed() +
            theme_void() +
            guides(col = "none")
    })
    p = patchwork::wrap_plots(p)+patchwork::plot_annotation(title="Ventral->Dorsal (blue->red)")
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
    marker_scores = read_csv("analysis/marker_scores.csv", show_col_types=FALSE) %>%
        column_to_rownames("sample_id") %>%
        as.matrix()
    marker_scores = t(marker_scores) / colMaxs(marker_scores)
    marker_scores = marker_scores[,as.character(colnames(barseq))]
    return(marker_scores)
}

plot_slice_summary = function(barseq, output_dir, filename="all_slices.pdf") {    
    all_slices = sort(unique(barseq$slice_id))
    pdf(file.path(output_dir, filename))
    for (i in all_slices) {
        print(plot_slice(barseq, i))
    }
    dev.off()
}

plot_slices = function(barseq, output_dir) {
    slice_dir = file.path(output_dir, "slices")
    dir.create(slice_dir, showWarnings = FALSE, recursive=TRUE)
    all_slices = sort(unique(barseq$slice_id))
    unique_labels = levels(as.factor(barseq$label))
    for (i in all_slices) {
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
    all_slices = sort(unique(barseq$slice_id))
    for (l in unique(barseq$label)) {
        file_name = gsub('\\.', '_', make.names(gsub('[\\(\\)]', '', l)))
        pdf(file.path(ct_dir, paste0(file_name, ".pdf")))
        for (i in all_slices) {
            print(plot_slice_ct(barseq, i, l))
        }
        dev.off()
    }
}

plot_slice = function(barseq, slice_to_plot, raster.dpi=150) {
    as.data.frame(colData(barseq)) %>%
        mutate(label = as.factor(label)) %>%
        filter(slice_id == slice_to_plot) %>%
        ggplot(aes(x = aligned_x, y = aligned_y, col = label)) +
        ggrastr::geom_point_rast(size=0.5, shape=16, stroke=0, raster.dpi = raster.dpi) +
        theme_void() +
        labs(x=NULL, y=NULL, col=NULL, title=paste0("Slice ", slice_to_plot)) +
        coord_equal() +
        lims(x = range(barseq$aligned_x), y = range(barseq$aligned_y)) +
        guides(colour = guide_legend(override.aes = list(size=2))) +
        theme(legend.text=element_text(size=8), legend.position = "bottom") +
        scale_color_discrete(na.value="gray90", drop=FALSE) 
}

plot_slice_ct = function(barseq, slice_to_plot, ct) {
    all_labels = as.factor(barseq$label)
    my_cols = make_cluster_cols(all_labels)
    full_plot = plot_slice(barseq, slice_to_plot)
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
