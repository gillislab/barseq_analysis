
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

make_bins = function(x, n_bins) {
    x = cut(x, breaks = seq(min(x,na.rm = TRUE),max(x, na.rm = TRUE),length.out = n_bins))
}

compute_flatmap = function(sce, group_by, supergroup_by=NULL) {
    result = tibble(x=sce$x, y=sce$y, group = sce[[group_by]],
                    supergroup = if(is.null(supergroup_by)) {"all"} else {sce[[supergroup_by]]}) %>%
        filter(!is.na(x) & ! is.na(y)) %>%
        group_by(x,y,group,supergroup) %>%
        tally() %>%
        ungroup() %>%
        complete(nesting(group,supergroup), nesting(x,y), fill=list("n"=0)) %>%
        group_by(x,y,supergroup) %>%
        mutate(n_total = sum(n), f=n/n_total) %>%
        group_by(group, supergroup) %>%
        mutate(f01 = f / max(f)) %>%
        ungroup()
    return(result)
}

plot_flatmap_by_group= function(to_plot, n_threshold=0) {
    p = to_plot %>%
        mutate(f01 = ifelse(n_total >= n_threshold, f01, 0)) %>%
        ggplot(aes(x,y,fill=f01)) +
        geom_raster(show.legend = FALSE) +
        facet_wrap(~ group) +
        theme_void() +
        scale_fill_gradient(low = gray(0.95), high = "darkblue", limits=c(0,NA)) +
        coord_equal()
    return(p)
}
