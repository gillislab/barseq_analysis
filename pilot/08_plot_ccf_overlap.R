
library(tidyverse)
library(Matrix)
source("dataset.R")


main = function() {
    # Load annotated BARseq
    barseq = load_barseq(normalization_factor = 10)
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    barseq$study_name = paste("barseq", barseq$slice, sep="_")
    clusters = read_csv("analysis/labels.csv")
    assertthat::are_equal(clusters$sample, colnames(barseq))
    barseq$cluster = clusters$cluster
    barseq$subclass = clusters$subclass
    barseq$class = clusters$class
    barseq = barseq[, barseq$class == "Glutamatergic"]
    ctx = barseq[, barseq$subclass %in% ctx_subclass()]
    non_ctx = barseq[, !(barseq$subclass %in% ctx_subclass())]
    
    subclass_overlap = tidy_overlap(overlap(ctx$subclass, ctx$CCFparentname))
    p1 = subclass_overlap %>%
        mutate(type = factor(type, rev(ctx_subclass()))) %>%
        ggplot(aes(x=CCF,y=type,size=overlap, color=overlap)) +
        geom_point() +
        theme_bw() +
        scale_size_area(max_size = 5) +
        scale_color_gradient(low = "gray90", high = "blue4", limits=c(0,NA)) +
        theme(axis.text.x = element_text(size=6, angle = 90, hjust=1,vjust=0.5),
              axis.text.y = element_text(size=6)) +
        labs(y = "Cortical subclasses (H2)", x=NULL)
    p1
    ggsave("figs/ccf/ctx_subclass.pdf", p1)

    threshold = 0.05
    subclass_overlap = overlap(non_ctx$cluster, non_ctx$ccf_name)
    subclass_overlap = subclass_overlap[, colMaxs(subclass_overlap) > threshold]
    col_order = order_rows_according_to_cols(t(subclass_overlap), 100)
    p1 = subclass_overlap %>%
        tidy_overlap() %>%
        mutate(type = factor(type, rev(rownames(subclass_overlap)))) %>%
        mutate(CCF = factor(CCF, colnames(subclass_overlap)[col_order])) %>%
        ggplot(aes(x=CCF,y=type,size=overlap, color=overlap)) +
        geom_point() +
        theme_bw() +
        scale_size_area(max_size = 5) +
        scale_color_gradient(low = "gray90", high = "blue4", limits=c(0,NA)) +
        theme(axis.text.x = element_text(size=6, angle = 90, hjust=1,vjust=0.5),
              axis.text.y = element_text(size=6)) +
        labs(y = "Non-cortical types (H3)", x=NULL)
    p1
    ggsave("figs/ccf/nonctx_cluster.pdf", p1)

    threshold = 0.05
    my_overlap = overlap(ctx$cluster, ctx$CCFparentname)
    my_overlap = my_overlap[, colMaxs(my_overlap) > threshold]
#    col_order = order_rows_according_to_cols(t(my_overlap), 100)
    p1 = my_overlap %>%
        tidy_overlap() %>%
        mutate(type = factor(type, rev(cell_type_order()))) %>%
#        mutate(CCF = factor(CCF, colnames(my_overlap)[col_order])) %>%
        mutate(CCF = as.factor(CCF)) %>%
        ggplot(aes(x=CCF,y=type,size=overlap, color=overlap)) +
        geom_point() +
        theme_bw() +
        scale_size_area(max_size = 5) +
        scale_color_gradient(low = "gray90", high = "blue4", limits=c(0,NA)) +
        theme(axis.text.x = element_text(size=6, angle = 90, hjust=1,vjust=0.5),
              axis.text.y = element_text(size=6)) +
        labs(y = "Cortical types (H3)", x=NULL)
    p1
    ggsave("figs/ccf/ctx_cluster.pdf", p1)

    threshold = 0.05
    all_subclasses = ctx_subclass()
    all_subclasses[grepl("RSP", all_subclasses)] = "RSP"
    all_subclasses = unique(all_subclasses)
    all_plots = lapply(all_subclasses, function(s) {
        is_s = !is.na(ctx$subclass) & grepl(s, ctx$subclass)
        my_overlap = overlap(ctx$cluster[is_s], ctx$CCFname[is_s])
        my_overlap = my_overlap[, colMaxs(my_overlap) > threshold]
        col_order = order_rows_according_to_cols(t(my_overlap), 100)
        my_plot = my_overlap %>%
            tidy_overlap() %>%
            mutate(type = factor(type, rev(cell_type_order()))) %>%
            mutate(CCF = factor(CCF, colnames(my_overlap)[col_order])) %>%
            ggplot(aes(x=CCF,y=type,size=overlap, color=overlap)) +
            geom_point() +
            theme_bw() +
            scale_size_area(max_size = 5, limits=c(0,0.5)) +
            scale_color_gradient(low = "gray90", high = "blue4", limits=c(0,0.5)) +
            theme(axis.text.x = element_text(size=6, angle = 90, hjust=1,vjust=0.5),
                  axis.text.y = element_text(size=6)) +
            labs(y = NULL, x=NULL) +
            ggtitle(s)
        return(my_plot)
    })
    pdf("figs/ccf/ctx_cluster_by_subclass.pdf", width = 10, height = 10)
    do.call(cowplot::plot_grid, c(all_plots, ncol = 2, align="hv"))
    dev.off()
    p_no_legend = lapply(all_plots, function(p) { p + guides(col="none", size="none") })
    pdf("figs/ccf/ctx_cluster_by_subclass_no_legend.pdf", width = 10, height = 10)
    do.call(cowplot::plot_grid, c(p_no_legend, ncol = 2, align="hv"))
    dev.off()
}

jaccard_index = function(x, y) {
    keep = !is.na(x) & !is.na(y)
    X = label_matrix(x[keep])
    Y = label_matrix(y[keep])
    overlap = crossprod(X,Y)
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

filter_ccf = function(jaccard_stats, topn=3) {
    to_keep = jaccard_stats %>%
        arrange(type, desc(jaccard_index)) %>%
        group_by(type) %>%
        filter(rank(-jaccard_index) <= topn) %>%
        pull(CCF) %>%
        unique()
    jaccard_stats %>%
        filter(CCF %in% to_keep) %>%
        mutate(CCF = factor(CCF, levels=to_keep)) %>%
        filter(jaccard_index>1e-10)
}

overlap = function(x, y) {
    keep = !is.na(x) & !is.na(y)
    X = label_matrix(x[keep])
    Y = label_matrix(y[keep])
    overlap = crossprod(X,Y)
    result = overlap  / rowSums(overlap)
    return(as.matrix(result))
}

tidy_overlap = function(x) {
    as_tibble(x, rownames="type") %>%
        pivot_longer(-type, names_to = "CCF", values_to = "overlap")    
}

order_rows_according_to_cols = function(M, alpha = 1) {
    M <- M**alpha
    row_score <- colSums(t(M)*seq_len(ncol(M)), na.rm=TRUE)/rowSums(M, na.rm=TRUE)
    return(order(row_score))
}

if (sys.nframe() == 0) {
    main()
}
