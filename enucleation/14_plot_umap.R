
library(tidyverse)
library(patchwork)
library(SingleCellExperiment)
library(rasterpdf)
source("dataset.R")


main = function() {
    # load and preformat data
    message("Loading data...")
    barseq = load_dataset("analysis/subglu")
    barseq$condition = ifelse(grepl("1L", barseq$sample) | grepl("3L", barseq$sample), "CTL", "ENU")
    
    # plot cell types of interest
    my_data = as_tibble(colData(barseq)) %>%
        select(CCFparentname, CCFname, h2, h3, h3_name, subclass_name, condition, sample) %>%
        cbind(as_tibble(reducedDim(barseq, "UMAP"))) %>%
        filter(CCFparentname %in% ctx_roi()) %>%
        mutate(ccf = CCFparentname, litter = substr(sample, 1, 4))
    my_ccf = c("VISp")
    
    # select background
    pdf("fig/umap_shift/select_background.pdf", width = 40/3)
    for (ct_prefix in c("l23it","l45it","l5it","l6it","l56it","l6b", "l5et", "l6ct","np","car3")) {
        my_ct = list("l23it"="L2/3 IT",
                     "l45it"="L4/5 IT",
                     "l5it"="L5 IT",
                     "l6it"="L6 IT",
                     "l56it"=c("L5 IT","L6 IT"),
                     "l6b"=c("L6b","L6b/CT"),
                     "l5et"="L5 ET",
                     "l6ct"="L6 CT",
                     "np"="NP",
                     "car3"="Car3")[[ct_prefix]]
        shift = read_csv(glue::glue("analysis/shift_comp/{ct_prefix}_h2.csv"), show_col_types = FALSE)
        shift$jsd = shift_jsd(my_data, my_ct)[shift$ccf]
        p1 = ggplot(shift, aes(log_or/log(2), -log(jsd)/log(10), label=ccf)) +
            geom_text() +
            theme_classic(base_size = 16) +
            labs(x="Neighborhood shift (log2(OR))",y="Similarity to UMAP shift (-log10(JSD))")
        p2 = plot_diff_density(my_data, my_ct, my_ccf, nlev = 10) + guides(fill = "none")
        p3 = plot_ctl_density_by_ccf(my_data, my_ct)
        print((p1 | (p2/p3)) + patchwork::plot_layout(widths = c(1,2)))
    }
    dev.off()
    
    # plot shift + selected background regions
    nlev=10
    pdf("fig/umap_shift/composition_shift.pdf", width = 60/3)   
    plot_composition_shift(my_data, c("L2/3 IT"), my_ccf, c("AUDp","RSPagl"), nlev = nlev)
    plot_composition_shift(my_data, c("L4/5 IT"), my_ccf, c("VISa","SSp-tr","VISal","VISl"), nlev = nlev)
    plot_composition_shift(my_data, c("L5 IT"), my_ccf, c("ACAv"), nlev = nlev)
    plot_composition_shift(my_data, c("L6 IT"), my_ccf, c("ACAv"), nlev = nlev)
    plot_composition_shift(my_data, c("L5 IT","L6 IT"), my_ccf, c("ACAv"), nlev = nlev)
    plot_composition_shift(my_data, c("L6b", "L6b/CT"), my_ccf, c("TEa","AUDv"), nlev = nlev)
    plot_composition_shift(my_data, c("L5 ET"), my_ccf, c("VISl","RSPagl","VISpor"), nlev = nlev)
    plot_composition_shift(my_data, c("L6 CT"), my_ccf, c("TEa","RSPagl"), nlev = nlev)
    plot_composition_shift(my_data, c("NP"), my_ccf, c("ACAv","TEa"), nlev = nlev)
    plot_composition_shift(my_data, c("Car3"), my_ccf, c("AUDv"), nlev = nlev)
    dev.off()
    
    # focus on L6b/CT, L6b, CT
    #my_data = as_tibble(colData(barseq)) %>%
    #    select(CCFparentname, CCFname, h2, h3, h3_name, subclass_name, condition, sample) %>%
    #    cbind(as_tibble(reducedDim(barseq, "UMAP_full"))) %>%
    #    filter(CCFparentname %in% ctx_roi())
    #pdf("fig/composition_shift_L6b_ct.pdf", width = 30/3)
    #plot_composition_shift(my_data, c("L6b", "L6b/CT", "L6 CT"), c("VISp","VISpor","MOs"), umap_lims = list(c(-5,0),c(7,15)))
    #dev.off()
    #pdf("fig/composition_shift_by_litter.pdf", width = 30/3)
    #plot_composition_shift_by_litter(barseq, my_data, c("L6b", "L6b/CT", "L6 CT"), c("VISp","VISpor","MOs"))
    #dev.off()
    
    # interpret shifts
    plot_interpreted_shift(my_data, "L2/3 IT")
    ggsave("fig/composition_shift_l23it.pdf")

    plot_interpreted_shift(my_data, "L6 IT")
    ggsave("fig/composition_shift_l6it.pdf")

    plot_interpreted_shift(my_data, "L4/5 IT")
    ggsave("fig/composition_shift_l45it.pdf")
    
    plot_interpreted_shift(my_data, c("L6b", "L6b/CT"))
    ggsave("fig/composition_shift_l6b.pdf")

    plot_interpreted_shift(my_data, c("L6 CT"))
    ggsave("fig/composition_shift_ct.pdf")
    
    # recompute UMAP for L6b/CT
    sce = barseq[, barseq$h2 %in% c("L6b","L6b/CT","L6 CT")]
    pca = readRDS("analysis//glu/pca.rds")
    reducedDim(sce, "PCA") = pca[colnames(sce),]
    #logcounts(sce) = log1p(cpm(sce))
    #sce = scater::runPCA(sce, ncomponents=30, exprs_values="logcounts")
    sce = scater::runUMAP(sce, dimred="PCA", n_neighbors=15)
    l6_data = as_tibble(colData(sce)) %>%
        select(CCFparentname, CCFname, h2, h3, h3_name, subclass_name, condition, sample) %>%
        cbind(as_tibble(reducedDim(sce, "UMAP"))) %>%
        filter(CCFparentname %in% ctx_roi())    
    plot_interpreted_shift(l6_data, c("L6b","L6b/CT","L6 CT"))
    ggsave("fig/composition_shift_L6b_ct.pdf")
}

load_dataset = function(input_dir) {
    barseq = load_glu()
    sample_cell = strsplit(colnames(barseq), split = ".", fixed = TRUE) # replace by regex for efficiency
    barseq$sample = sapply(sample_cell, "[", 1)
    barseq$slice = as.factor(sapply(strsplit(sapply(sample_cell, "[", 2), split = "_", fixed = TRUE),"[", 1))
    barseq$n_genes = colSums(counts(barseq)>0)
    barseq$slice_number = as.numeric(as.character(barseq$slice))
    barseq$slice_id = paste(barseq$sample, str_pad(barseq$slice, 2, pad="0"), sep="|")
    barseq$original_label = barseq$label
    barseq = barseq[, sample.int(ncol(barseq))] # randomize order for plotting
    subclass_id = deframe(read_csv("analysis/glu/cluster.csv.gz", show_col_types=FALSE))
    barseq$subclass_id = as.factor(subclass_id[colnames(barseq)])
    subclass_annotation = deframe(read_csv("analysis/glu/cluster_annotation.csv", show_col_types=FALSE))
    barseq$subclass_name = subclass_annotation[as.character(barseq$subclass_id)]
    umap = column_to_rownames(as.data.frame(read_csv(file.path(input_dir, "subumap.csv.gz"), show_col_types=FALSE)), "cell_id")
    reducedDim(barseq, "UMAP") = umap[colnames(barseq),]
    umap = column_to_rownames(as.data.frame(read_csv("analysis/glu/umap.csv.gz", show_col_types=FALSE)), "cell_id")
    reducedDim(barseq, "UMAP_full") = umap[colnames(barseq),]
    return(barseq)
}

load_glu = function() {
    labels = load_labels()
    result = lapply(set_names(sample_names()), function(my_sample) {
        sce = load_sample(my_sample)
        my_labels = labels[match(colnames(sce), labels$cell_id),]
        colData(sce) = cbind(colData(sce), select(my_labels, -cell_id))
        sce = sce[, sce$h1 == "Glutamatergic" & sce$h2 %in% ctx_subclass()]
        return(sce)
    })
    result = MetaNeighbor::mergeSCE(result)
    return(result)
}

plot_raster_umap = function(barseq, color_by, text_by=NULL, point_size=0.1, text_size=3, raster.dpi=150, alpha=1, umap_name="UMAP") {
    to_plot = as_tibble(reducedDim(barseq, umap_name))
    colnames(to_plot) = c("UMAP1", "UMAP2")
    if (color_by %in% rownames(barseq)) {
        to_plot[[color_by]] = logcounts(barseq)[color_by,]
    } else {
        to_plot[[color_by]] = barseq[[color_by]]
    }
    if (!is.null(text_by)) {
        to_plot[[text_by]] = barseq[[text_by]]
    }    
    plot_raster_umap_(to_plot, color_by, text_by, point_size, text_size, raster.dpi, alpha)
}

plot_raster_umap_ = function(to_plot, color_by, text_by=NULL, point_size=0.1, text_size=3, raster.dpi=150, alpha=1) {
    # trick to make dot size as small as possible: https://stackoverflow.com/questions/34638902/point-size-in-ggplot-2-0-0
    result=to_plot %>%
        ggplot(aes(x=UMAP1, y=UMAP2, col=.data[[color_by]])) +
        ggrastr::geom_point_rast(size=point_size, shape=16, stroke=0, raster.dpi = raster.dpi, alpha=alpha) +
        theme_classic()
    if (is.numeric(to_plot[[color_by]])) {
        result = result + scale_color_viridis_c()
    }
    if (!is.null(text_by)) {
        to_plot = to_plot %>%
            group_by(.data[[text_by]]) %>%
            summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
        result = result +
            geom_text(data=to_plot, aes(label=.data[[text_by]]), size=text_size, color="black")
    }
    return(result)
}

plot_composition_shift = function(my_data, my_h2, my_ccf, background_ccf, nlev=7) {
    p0 = plot_raster_umap_(my_data[my_data$h2 %in% my_h2,], "h3_name", "h3_name", alpha = 0.2) +
        guides(col="none") +
        theme_void() +
        coord_equal()
    background_data = filter(my_data, CCFparentname %in% background_ccf)
    p1 = plot_density_v2(my_data, my_h2, my_ccf, nlev = nlev)
    p2 = plot_ctl_density(my_data, my_h2, nbin = 50, nlev=nlev) + ggtitle("all regions")
    p3 = plot_ctl_density(background_data, my_h2, nbin = 50, nlev=nlev) + ggtitle("background regions (see right)")
    p4 = plot_ctl_density_by_ccf(background_data, my_h2, nbin = 50, nlev = nlev)
    (p0 | (p1 / (p2|p3)) | p4) + patchwork::plot_annotation(my_h2)
}

plot_composition_shift_original = function(my_data, my_ct, my_ccf, umap_lims = NULL, nbin=20, nlev=7, smoothing="kde", n_interp=3) {
    p1 = plot_raster_umap_(my_data[my_data$h2 %in% my_ct,], "h3_name", "h3_name", alpha = 0.2) +
        guides(col="none") +
        theme_void() +
        theme(aspect.ratio = 1)
    p2 = plot_ccf_overlap(my_data, my_ct) + guides(col="none",size="none") + ggtitle("Main repartition in CTL")
    p3 = plot_density_v2(my_data, my_ct, my_ccf, nbin = nbin, nlev = nlev, smoothing=smoothing, n_interp=n_interp)
    p4 = plot_ctl_density(my_data, my_ct, nbin = 2*nbin, nlev = nlev, smoothing=smoothing, n_interp=n_interp/2)
    if (!is.null(umap_lims)) {
        p1 = p1 + coord_equal(xlim = umap_lims[[1]], ylim = umap_lims[[2]])
        p3 = p3 & coord_equal(xlim = umap_lims[[1]], ylim = umap_lims[[2]])
        p4 = p4 + coord_equal(xlim = umap_lims[[1]], ylim = umap_lims[[2]])
    }
    ((p1/p4) + plot_layout(heights=c(1,1)) | p3) + plot_layout(widths = c(1,2))
}

plot_density_legacy = function(my_data, my_ct, my_ccf) {
    my_data %>%
        dplyr::filter(h2 %in% my_ct) %>%
        mutate(ccf = ifelse(CCFparentname %in% my_ccf, CCFparentname, "other")) %>%
        mutate(ccf = factor(ccf, c(my_ccf, "other"))) %>%
        ggplot(aes(x=UMAP1,y=UMAP2)) +
        geom_density_2d_filled(contour_var = "ndensity") +
        facet_grid(ccf~condition) +
        theme_void() +
        coord_equal()
}

plot_ctl_density = function(my_data, my_ct, nbin = 20, nlev=7, lev_type="linear", smoothing="kde", n_interp=3) {
    rdylbu = colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdBu"))(2*nlev-1)
    my_data = filter(my_data, h2 %in% my_ct & condition == "CTL")
    if (smoothing != "kde") {
        ctl_data = my_data %>%
            mutate(umap_x = cut(UMAP1, seq(min(UMAP1), max(UMAP1), length.out = nbin), include.lowest = TRUE),
                   umap_y = cut(UMAP2, seq(min(UMAP2), max(UMAP2), length.out = nbin), include.lowest = TRUE)) %>%
            group_by(umap_x, umap_y) %>%
            tally() %>%
            ungroup() %>%
            mutate(f = n / max(n)) %>%
            ungroup()
    } else {
        ctl_data = with(my_data, smooth_kde2d(UMAP1, UMAP2, n=nbin))
    }
    if (smoothing == "linear") {
        ctl_data = ctl_data %>%
            reframe(smooth_density(umap_x, umap_y, f, n_interp)) %>%
            drop_na()
    }
    p_ctl = ctl_data %>%
        ggplot(aes(x=as.numeric(umap_x),y=as.numeric(umap_y),z=f)) +
        theme_void() +
        coord_equal() +
        scale_fill_manual(values = rdylbu[nlev:(2*nlev-1)])
    if (lev_type == "linear") {
        p_ctl = p_ctl + stat_contour_filled(bins = nlev, color=NA)
    } else {
        p_ctl = p_ctl + stat_contour_filled(breaks=quantile(summarized_data$f, seq(0,1,length.out = nlev+1)))
    }
    p_ctl & guides(fill = "none")
}

plot_density_v2 = function(my_data, my_ct, my_ccf, nbin = 20, nlev=7, lev_type="linear", smoothing="kde", n_interp=3) {
    my_data = my_data %>%
        filter(h2 %in% my_ct & CCFparentname %in% my_ccf) %>%
        mutate(ccf = factor(CCFparentname, my_ccf))
    if (smoothing != "kde") {
        summarized_data = my_data %>%
            mutate(umap_x = cut(UMAP1, seq(min(UMAP1), max(UMAP1), length.out = nbin), include.lowest = TRUE),
                   umap_y = cut(UMAP2, seq(min(UMAP2), max(UMAP2), length.out = nbin), include.lowest = TRUE)) %>%
            group_by(umap_x, umap_y, ccf, condition) %>%
            tally() %>%
            group_by(ccf, condition) %>%
            mutate(f = n / max(n)) %>%
            ungroup() %>%
            complete(nesting(umap_x,umap_y),ccf,condition, fill = list(f=0))
    }  else {
        summarized_data = my_data %>%
            group_by(ccf, condition) %>%
            reframe(smooth_kde2d(UMAP1, UMAP2, n=nbin))
    }
    if (smoothing == "linear") {
        summarized_data = summarized_data %>%
            group_by(ccf, condition) %>%
            reframe(smooth_density(umap_x, umap_y, f, n_interp)) %>%
            drop_na()
    }
    rdylbu = colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdBu"))(2*nlev-1)
    p2a = summarized_data %>%
        filter(ccf %in% my_ccf & condition == "CTL") %>%
        ggplot(aes(x=as.numeric(umap_x),y=as.numeric(umap_y),z=f)) +
        theme_void() +
        coord_equal() +
        facet_grid(ccf~.) +
        scale_fill_manual(values = rdylbu[nlev:(2*nlev-1)])
    p2b = summarized_data %>%
        filter(ccf %in% my_ccf & condition == "ENU") %>%
        ggplot(aes(x=as.numeric(umap_x),y=as.numeric(umap_y),z=f)) +
        theme_void() +
        coord_equal() +
        facet_grid(ccf~.) +
        scale_fill_manual(values = rev(rdylbu[1:nlev]))
    if (lev_type == "linear") {
        p2a = p2a + stat_contour_filled(bins = nlev, color=NA)
        p2b = p2b + stat_contour_filled(bins = nlev, color=NA)
    } else {
        p2a = p2a + stat_contour_filled(breaks=quantile(summarized_data$f, seq(0,1,length.out = nlev+1)))
        p2b = p2b + stat_contour_filled(breaks=quantile(summarized_data$f, seq(0,1,length.out = nlev+1)))
    }
    p_diff = plot_diff_density(my_data, my_ct, my_ccf, nbin, nlev, lev_type, smoothing=smoothing, n_interp=n_interp)
    p_shift = (p2a|p2b|p_diff) #& theme(strip.text = element_blank())
    p_shift & guides(fill = "none")
}

smooth_kde2d = function(x,y,nbin,lims=c(range(x),range(y))) {
    dens2d = MASS::kde2d(x, y, n=nbin,lims=lims)
    result = tibble(umap_x = rep(dens2d$x, nbin), umap_y = rep(dens2d$y, each=nbin), f = c(dens2d$z))
    return(result)
}

smooth_density = function(umap_x, umap_y, f, n_interp=2) {
    umap_x = as.numeric(umap_x)
    umap_y = as.numeric(umap_y)
    target_x = seq(min(umap_x), max(umap_x), 1/n_interp)
    target_y = seq(min(umap_y), max(umap_y), 1/n_interp)
    result = interp::interp(umap_x, umap_y, f, output = "grid", xo = target_x, yo=target_y, method="linear")$z
    rownames(result) = as.character(target_x)
    colnames(result) = as.character(target_y)
    result = as_tibble(result, rownames = "umap_x") %>%
        pivot_longer(-umap_x, values_to = "f", names_to = "umap_y") %>%
        mutate(umap_x = as.numeric(umap_x), umap_y = as.numeric(umap_y))
    return(result)
}

plot_diff_density = function(my_data, my_ct, my_ccf, nbin=50, nlev = 10, lev_type = "linear", smoothing="kde", n_interp=3, relative=FALSE) {
    my_data = my_data %>%
            mutate(ccf = CCFparentname, litter = substr(sample, 1, 4)) %>%
            filter(h2 %in% my_ct & ccf %in% my_ccf) %>%
            mutate(ccf = factor(ccf, my_ccf))
    if (smoothing == "kde") {
        summarized_data = my_data %>%
            group_by(ccf, condition, litter) %>%
            reframe(smooth_kde2d(UMAP1, UMAP2, nbin, lims = c(range(my_data$UMAP1), range(my_data$UMAP2))))
    } else {
        summarized_data = my_data %>%
            mutate(umap_x = cut(UMAP1, seq(min(UMAP1), max(UMAP1), length.out = nbin), include.lowest = TRUE),
                   umap_y = cut(UMAP2, seq(min(UMAP2), max(UMAP2), length.out = nbin), include.lowest = TRUE)) %>%
            group_by(umap_x, umap_y, ccf, condition, litter) %>%
            tally() %>%
            group_by(ccf, condition, litter) %>%
            mutate(f = n / sum(n)) %>%
            ungroup() %>%
            complete(nesting(umap_x,umap_y),ccf,condition,litter,fill = list(f=0))
    }
    summarized_data = summarized_data %>%
        pivot_wider(id_cols=c(umap_x, umap_y, ccf, litter), names_from = "condition", values_from = "f") %>%
        mutate(delta_f = if(relative) {log2((ENU+1e-2)/(CTL+1e-2))} else {(ENU-CTL)} )
    average_data = summarized_data %>%
        group_by(umap_x,umap_y,ccf) %>%
        summarize(delta_f = mean(delta_f))# %>%
        #group_by(ccf) %>%
        #mutate(delta_f = delta_f / max(abs(min(delta_f)), max(delta_f)))
    if (smoothing == "linear") {
        average_data = average_data %>%
            group_by(ccf) %>%
            reframe(smooth_density(umap_x, umap_y, delta_f, n_interp)) %>%
            rename("f" = "delta_f") %>%
            drop_na()
    }
    min_delta = min(average_data$delta_f)
    max_delta = max(average_data$delta_f)
    highest_delta = max(-min_delta, max_delta)
    breaks = seq(-highest_delta, highest_delta, length.out = 2*nlev)
    rdylbu = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdBu"))(2*nlev-1))
    keep_break = breaks >= min_delta & breaks <= max_delta
    keep_break = c(FALSE, keep_break[-length(keep_break)]) | c(keep_break[-1], FALSE)
    keep_col = keep_break[-1] & keep_break[-length(keep_break)]
    breaks = breaks[keep_break]
    rdylbu = rdylbu[keep_col]
    p = average_data %>%
        ggplot(aes(x=as.numeric(umap_x),y=as.numeric(umap_y),z=delta_f)) +
        theme_void() +
        coord_equal() +
        facet_grid(ccf~.) +
        scale_fill_manual(values = rdylbu)
    if (lev_type == "linear") {
        p = p + stat_contour_filled(breaks = breaks, color=NA)
    } else {
        stop("not implemented")
        #p = p + stat_contour_filled(breaks=quantile(average_data$delta_f, seq(0,1,length.out = nlev+1)))
    }
    return(p)
}

plot_ccf_overlap = function(my_data, my_ct, threshold=0.05) {
    sub_data = filter(my_data, h2 %in% my_ct & condition == "CTL")
    my_overlap = jaccard_index(sub_data$h3_name, sub_data$CCFparentname)
    my_overlap = my_overlap[, colMaxs(my_overlap) > threshold]
    p1 = my_overlap %>%
        tidy_overlap() %>%
#        mutate(type = factor(type, rev(cell_type_order()))) %>%
#        mutate(CCF = factor(CCF, c(my_ccf, "other"))) %>%
        ggplot(aes(x=CCF,y=type,size=overlap, color=overlap)) +
        geom_point() +
        theme_classic(base_size = 16) +
#        scale_size_area(max_size = 5) +
        scale_color_gradient(low = "gray90", high = "blue4", limits=c(0,NA)) +
        theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5,size=10),
              axis.text.y = element_text(size=10)) +
        labs(y = NULL, x=NULL)
    p1
}

overlap = function(x, y) {
    keep = !is.na(x) & !is.na(y)
    X = label_matrix(x[keep])
    Y = label_matrix(y[keep])
    overlap = Matrix::crossprod(X,Y)
    result = overlap  / rowSums(overlap)
    return(as.matrix(result))
}

jaccard_index = function(x, y) {
    keep = !is.na(x) & !is.na(y)
    X = label_matrix(x[keep])
    Y = label_matrix(y[keep])
    overlap = Matrix::crossprod(X,Y)
    union = outer(colSums(X), colSums(Y), "+")
    result = overlap/union
    return(as.matrix(result))
}

label_matrix = function(labels) {
    labels = as.factor(labels)
    result = Matrix::sparseMatrix(i=1:length(labels),j=as.numeric(labels),
                                  dimnames = list(names(labels), levels(labels)))
    return(1*result)
}

tidy_overlap = function(x) {
    as_tibble(x, rownames="type") %>%
        pivot_longer(-type, names_to = "CCF", values_to = "overlap")    
}

compute_umap_density = function(my_data, my_ct, nbin=50) {
    result = my_data %>%
        filter(h2 %in% my_ct) %>%
        mutate(ccf = CCFparentname) %>%
        mutate(umap_x = cut(UMAP1, seq(min(UMAP1), max(UMAP1), length.out = nbin), include.lowest = TRUE),
               umap_y = cut(UMAP2, seq(min(UMAP2), max(UMAP2), length.out = nbin), include.lowest = TRUE)) %>%
        group_by(umap_x, umap_y, ccf, condition) %>%
        tally() %>%
        group_by(ccf, condition) %>%
        mutate(f = n / sum(n)) %>%
        ungroup() %>%
        complete(nesting(umap_x,umap_y),ccf,condition, fill = list(f=0))
    return(result)
}


shift_jsd = function(my_data, my_ct) {
    all_ccf = sort(unique(my_data$CCFparentname))
    all_density = compute_umap_density(my_data, my_ct)
    enu_shift = all_density %>%
        filter(ccf == "VISp") %>%
        pivot_wider(id_cols = c(umap_x, umap_y), names_from = condition, values_from = f) %>%
        mutate(diff = ENU-CTL) %>%
        mutate(diff = pmax(diff, 0)) %>%
        select(umap_x, umap_y, diff)
    enu_shift$diff = enu_shift$diff / sum(enu_shift$diff)
    jsd = sapply(set_names(all_ccf), function(my_ccf) {
        ccf_density = all_density %>%
            filter(ccf == my_ccf & condition == "CTL") %>%
            full_join(enu_shift, by = c("umap_x", "umap_y")) %>%
            select(p = f, q = diff)
        ccf_density$m = (ccf_density$p + ccf_density$q)/2
        kl_p = with(ccf_density, sum(p*log(p/m), na.rm = TRUE))
        kl_q = with(ccf_density, sum(q*log(q/m), na.rm = TRUE))
        jsd = (kl_p + kl_q)/2
        return(jsd)
    })
    return(jsd)
}

plot_interpreted_shift = function(my_data, my_ct) {
    all_ccf = sort(unique(my_data$CCFparentname))
    all_density = compute_umap_density(my_data, my_ct)
    jsd = sapply(set_names(all_ccf), function(my_ccf) {
        ccf_density = all_density %>%
            filter(ccf == my_ccf) %>%
            pivot_wider(id_cols = c(umap_x, umap_y), names_from = condition, values_from = f)
        ccf_density$M = (ccf_density$CTL + ccf_density$ENU)/2
        kl_ctl = with(ccf_density, sum(CTL*log(CTL/M), na.rm = TRUE))
        kl_enu = with(ccf_density, sum(ENU*log(ENU/M), na.rm = TRUE))
        jsd = (kl_ctl + kl_enu)/2
        return(jsd)
    })   
    enu_shift = all_density %>%
        filter(ccf == "VISp") %>%
        pivot_wider(id_cols = c(umap_x, umap_y), names_from = condition, values_from = f) %>%
        mutate(diff = ENU-CTL) %>%
        mutate(diff = pmax(diff, 0)) %>%
        select(umap_x, umap_y, diff)
    enu_shift$diff = enu_shift$diff / sum(enu_shift$diff)
    jsd = sapply(set_names(all_ccf), function(my_ccf) {
        ccf_density = all_density %>%
            filter(ccf == my_ccf & condition == "CTL") %>%
            full_join(enu_shift, by = c("umap_x", "umap_y")) %>%
            select(p = f, q = diff)
        ccf_density$m = (ccf_density$p + ccf_density$q)/2
        kl_p = with(ccf_density, sum(p*log(p/m), na.rm = TRUE))
        kl_q = with(ccf_density, sum(q*log(q/m), na.rm = TRUE))
        jsd = (kl_p + kl_q)/2
        return(jsd)
    })
    p = plot_composition_shift(my_data, my_ct, c("VISp", names(sort(jsd))[1:3]))
    p
}

plot_ctl_density_by_ccf = function(my_data, my_ct, nbin = 20, nlev=7) {
    rdylbu = colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdBu"))(2*nlev-1)
    my_data = filter(my_data, h2 %in% my_ct & condition == "CTL")
    ctl_data = my_data %>%
        group_by(CCFparentname) %>%
        reframe(smooth_kde2d(UMAP1, UMAP2, n=nbin)) %>%
        group_by(CCFparentname) %>%
        mutate(f = f / max(f))
    p_ctl = ctl_data %>%
        ggplot(aes(x=as.numeric(umap_x),y=as.numeric(umap_y),z=f)) +
        theme_void() +
        coord_equal() +
        scale_fill_manual(values = rdylbu[nlev:(2*nlev-1)]) +
        facet_wrap(~ CCFparentname) +
        stat_contour_filled(bins = nlev, color=NA)
    p_ctl & guides(fill = "none")
}


if (sys.nframe() == 0) {
    main()
}
