# -*- coding: utf-8 -*-
# for the rescaling NMF loadings, see
# https://stats.stackexchange.com/questions/46469/in-non-negative-matrix-factorization-are-the-coefficients-of-features-comparabl
# -> best practice is to scale W using l1 or l2-normalization, then scale H accordingly
# is it necessary here? -> does not change most analyses (unaffected by scaling)

library(tidyverse)
library(SingleCellExperiment)
library(NMF)
source("dataset.R")
source("gene_set_enrichment.R")
#source("21_make_pseudobulks.R")


main = function() {
    # load BARseq data
    subclass_pseudobulks = readRDS("data/subclass_pseudobulks.rds")
    subclass_pseudobulks = subclass_pseudobulks[!rownames(subclass_pseudobulks) %in% c("Slc17a7", "Gad1"),]
    subclass_pseudobulks = subclass_pseudobulks[, !(subclass_pseudobulks$subclass %in% c("Car3", "RSP UL", "RSP DL"))]
    cpm(subclass_pseudobulks) = convert_to_cpm(counts(subclass_pseudobulks),
                                               total_counts = median(colSums(counts(subclass_pseudobulks))))
    pseudobulk_matrix = as_tibble(t(cpm(subclass_pseudobulks))) %>%
        mutate(subclass = subclass_pseudobulks$subclass, bin_id = subclass_pseudobulks$bin_id) %>%
        pivot_longer(c(-subclass, -bin_id), names_to = "gene", values_to = "expr") %>%
        mutate(var = paste(gene,subclass,sep="|")) %>%
        pivot_wider(bin_id, names_from = var, values_from = expr) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=sqrt(colSums(pseudobulk_matrix**2))) # l2-normalization
    
    # load nmf results 
    nmf_tidy = read_csv("analysis/spatial_patterns/spatial_nmf_scaled.csv")[1:11] %>%
        pivot_longer(-bin_id, names_to = "nmf", values_to = "coef") %>%
        mutate(nmf = nmf_name()[nmf]) %>%
        drop_na() %>%
        add_pos()
    nmf_loading = read_csv("analysis/spatial_patterns/spatial_nmf_scaled_loadings.csv") %>%
        mutate(nmf = nmf_name()[nmf]) %>%
        drop_na()
    nmf_result = readRDS("analysis/spatial_patterns/spatial_nmf_scaled.rds")
    w = basis(nmf_result)
    h = coef(nmf_result)
    var_explained = sapply(1:ncol(w), function(i) {
        1 - colSums((scaled_pseudobulk - w[,i,drop=FALSE] %*% h[i,,drop=FALSE])**2)
    })
    colnames(var_explained) = paste0("NMF", 1:ncol(w))
    var_tidy = as_tibble(var_explained, rownames="var") %>%
        pivot_longer(-var, names_to = "nmf", values_to = "var_exp") %>%
        separate(var, sep = "\\|", into = c("gene","subclass")) %>%
        mutate(nmf = nmf_name()[nmf]) %>%
        drop_na()
    total_var_explained = enframe(1 - colSums((scaled_pseudobulk - w %*% h)**2), "")

    # questions:
    # - which CCF regions are preferentially associated with each factor?
    # - which genes are preferentially associated with each factor?
    # - which subclasses are preferentially associated with each factor?
    plot_curated_nmf(nmf_tidy)
    plot_ccf_overlap(nmf_tidy)
    plot_gene_nmf_association(nmf_tidy, nmf_loading, subclass_pseudobulks)
    plot_subclass_nmf_association(nmf_loading)
    plot_cluster_nmf_spatial_association(nmf_tidy)
    plot_cluster_nmf_gene_association(nmf_loading)
}

nmf_name = function() {
    result = c(
        "NMF1"="NMF1_TEa_AUD",
        "NMF7"="NMF7_SSs_AUD",
        "NMF3"="NMF3_ACA_MOs",
        "NMF10"="NMF10_RSP_ACA",
        "NMF4"="NMF4_MOp_SSp",
        "NMF8"="NMF8_SSp",
        "NMF5"="NMF5_VIS"
    )
    return(result)
}

plot_curated_nmf = function(nmf_tidy) {
    nmf_tidy %>%
        mutate(nmf = factor(nmf, nmf_name())) %>%
        group_by(nmf) %>%
        mutate(coef = coef / max(coef)) %>%
        ungroup() %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=coef)) +
        geom_rect() +
        facet_wrap(~ nmf) +
        theme_classic() +
        labs(x="M-L axis",y="A-P axis") +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,NA))
    ggsave("figs/spatial_patterns/gene_association/nmf_curated.pdf")
}

plot_ccf_overlap = function(nmf_tidy) {
    roi = load_roi_info()
    nmf_matrix = nmf_tidy %>%
        pivot_wider(bin_id, names_from = nmf, values_from = coef) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    nmf_matrix = nmf_matrix[, nmf_name()]
    roi_tidy = roi %>%
        filter(bin_id %in% rownames(nmf_matrix)) %>%
        pivot_wider(names_from = roi, values_from = n_cells, values_fill = 0) 
    roi_matrix = as.matrix(column_to_rownames(roi_tidy, "bin_id"))
    roi_matrix = roi_matrix[rownames(nmf_matrix),]

    ccf_overlap = cor(nmf_matrix, roi_matrix, method="s")
    positive_overlap = ccf_overlap
    positive_overlap[positive_overlap<0] = 0
    to_plot = ccf_overlap[, order_rows_according_to_cols(t(positive_overlap))]
    p1 = as_tibble(to_plot, rownames = "nmf") %>%
        pivot_longer(-nmf, names_to = "roi", values_to = "s_cor") %>%
        mutate(nmf = factor(nmf, rev(rownames(to_plot)))) %>%
        mutate(roi = factor(roi, colnames(to_plot))) %>%
        ggplot(aes(x=roi,y=nmf,fill=s_cor)) +
        geom_raster() +
        theme_bw(base_size=10) +
        labs(x=NULL,y=NULL) +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))  +
        scale_fill_gradient2(low = "white", mid = "gray90", high = "darkblue", midpoint = 0)
    p1
    ggsave("figs/spatial_patterns/gene_association/ccf_associations.pdf",p1)
    
    # plot strongest associations
    roi_pos = roi_tidy %>%
        pivot_longer(-bin_id, names_to = "nmf", values_to = "coef") %>%
        add_pos()
    my_threshold = 0.25
    pdf("figs/spatial_patterns/gene_association/ccf_top.pdf")
    for (my_nmf in rownames(ccf_overlap)) {
        p1 = plot_nmf(nmf_tidy, my_nmf)
        top_ccf = round(sort(ccf_overlap[my_nmf,], decreasing = TRUE), 2)
        top_ccf = top_ccf[top_ccf > my_threshold]
        p2 = lapply(names(top_ccf), function(my_ccf) {
            plot_nmf(roi_pos, my_ccf) +
                ggtitle(paste0(my_ccf, " (", top_ccf[my_ccf], ")"))
        })
        do.call(gridExtra::grid.arrange, c(list(p1), p2, list(ncol=2)))
    }
    dev.off()
}

load_roi_info = function() {
    barseq = load_annotated_data()
    barseq$bin_id = barseq$slice*100 + barseq$ml
    barseq = barseq[, !(barseq$subclass %in% c("Car3", "RSP UL", "RSP DL"))]
    result = tibble(bin_id = barseq$bin_id, roi = barseq$CCFparentname) %>%
        drop_na(roi) %>%
        group_by(bin_id, roi) %>%
        tally(name = "n_cells") %>%
        ungroup()
    return(result)
}

add_pos = function(df) {
    pseudobulks = readRDS("data/subclass_pseudobulks.rds")
    bin_to_pos = tibble(bin_id=pseudobulks$bin_id, ap=pseudobulks$ap,
                        ml_min=pseudobulks$ml_min, ml_max=pseudobulks$ml_max) %>%
        distinct()
    result = left_join(df, bin_to_pos)
    return(result)
}

order_rows_according_to_cols = function(M, alpha = 100) {
    base_score = 1-rowSums(M)
    M <- M**alpha
    row_score <- colSums(t(M)*seq_len(ncol(M)), na.rm=TRUE)/rowSums(M, na.rm=TRUE)
    return(order(base_score + row_score))
}

plot_nmf = function(nmf_tidy, my_nmf) {
    nmf_tidy %>%
        filter(nmf == my_nmf) %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=coef)) +
        geom_rect(show.legend=FALSE) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        labs(x="M-L axis",y="A-P axis") +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,NA)) +
        ggtitle(my_nmf)
}


plot_gene_nmf_association = function(nmf_tidy, var_tidy, subclass_pseudobulks) {
    gene_loading = var_tidy %>%
        group_by(gene, nmf) %>%
        summarize(var_exp = mean(var_exp), .groups="drop")
    # compute null var. exp.: permute var. exp. across each NMF factor
    set.seed(17)
    null_var_exp = var_tidy %>%
        group_by(nmf) %>%
        mutate(var_exp = sample(var_exp)) %>%
        group_by(gene, nmf) %>%
        summarize(var_exp = mean(var_exp), .groups="drop")
    # null quantiles -> 11%, 18%, 21%, 25%
    quantile(null_var_exp$var_exp, probs = c(0.5, 0.9, 0.95, 0.99))
    my_threshold = quantile(null_var_exp$var_exp, probs = c(0.99))
    bind_rows(list(observed = gene_loading, null = null_var_exp), .id = "type") %>%
        ggplot(aes(x = var_exp, col = type)) +
        geom_density() +
        theme_bw() +
        geom_vline(xintercept = my_threshold, linetype = "dashed") +
        facet_wrap(~ nmf)
    # keep only genes that have an association > 99% null quantile with some NMF
    gene_loading = gene_loading %>%
        pivot_wider(gene, names_from = "nmf", values_from = "var_exp") %>%
        column_to_rownames("gene") %>%
        as.matrix()
    gene_loading = gene_loading[, nmf_name()]
    keep_gene = rowMaxs(gene_loading) > my_threshold
    to_plot = gene_loading[keep_gene,]
    to_plot = to_plot[order_rows_according_to_cols(to_plot),]
    p1 = as_tibble(to_plot, rownames="gene") %>%
        pivot_longer(-gene, names_to = "nmf", values_to = "var_exp") %>%
        mutate(nmf = factor(nmf, rev(colnames(to_plot)))) %>%
        mutate(gene = factor(gene, rownames(to_plot))) %>%
        mutate(significant = var_exp > my_threshold) %>%
        mutate(var_exp=100*var_exp) %>%
        ggplot(aes(x=gene,y=nmf)) +
        geom_point(aes(size=var_exp, col=significant, fill = var_exp), pch=21, stroke=0.2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))  +
        scale_radius(range = c(0,4)) +
        scale_color_manual(values = c("white", "black")) +
        scale_fill_gradient2(low = "white", mid = "gray90", high = "darkblue", midpoint = 100*my_threshold, limits=c(0, NA)) +
        labs(x=NULL,y=NULL)
    p1
    ggsave("figs/spatial_patterns/gene_association/gene_associations.pdf", p1, height = 2)
    # filter for genes that are associated with a single NMF?
    n_nmf = rowSums(gene_loading > my_threshold)
    p1 = enframe(n_nmf, "gene", "n_nmf") %>%
        ggplot(aes(x = n_nmf)) +
        geom_bar() +
        theme_bw(base_size = 20) +
        labs(x = "Number of NMF factors", "Number of genes")
    p1
    ggsave("figs/spatial_patterns/gene_association/gene_number_associations.pdf", p1)
    keep_gene = n_nmf == 1
    to_plot = gene_loading[keep_gene,]
    to_plot = to_plot[order_rows_according_to_cols(to_plot),]
    p1 = as_tibble(to_plot, rownames="gene") %>%
        pivot_longer(-gene, names_to = "nmf", values_to = "var_exp") %>%
        mutate(nmf = factor(nmf, rev(colnames(to_plot)))) %>%
        mutate(gene = factor(gene, rownames(to_plot))) %>%
        mutate(significant = var_exp > my_threshold) %>%
        mutate(var_exp=100*var_exp) %>%
        ggplot(aes(x=gene,y=nmf)) +
        geom_point(aes(size=var_exp, col=significant, fill = var_exp), pch=21, stroke=0.2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))  +
        scale_radius(range = c(0,4)) +
        scale_color_manual(values = c("white", "black")) +
        scale_fill_gradient2(low = "white", mid = "gray90", high = "darkblue", midpoint = 100*my_threshold, limits=c(0, NA)) +
        labs(x=NULL,y=NULL)
    p1
    ggsave("figs/spatial_patterns/gene_association/gene_single_associations.pdf", p1, height = 2)
    
    # show best gene(s) associated with each NMF -> keep top 3?
    keep_gene = n_nmf == 1 & (rownames(gene_loading) != "Slc30a3")
    to_plot = gene_loading[keep_gene,]
    pdf("figs/spatial_patterns/gene_association/genes_top.pdf")
    for (my_nmf in colnames(to_plot)) {
        top_genes = round(100*sort(to_plot[,my_nmf], decreasing = TRUE)[1:3])
        p1 = plot_nmf(nmf_tidy, my_nmf)
        p2 = plot_gene(subclass_pseudobulks, names(top_genes)[1]) +
            ggtitle(paste0(names(top_genes)[1], " (", top_genes[1], ")"))
        p3 = plot_gene(subclass_pseudobulks, names(top_genes)[2]) +
            ggtitle(paste0(names(top_genes)[2], " (", top_genes[2], ")"))
        p4 = plot_gene(subclass_pseudobulks, names(top_genes)[3]) +
            ggtitle(paste0(names(top_genes)[3], " (", top_genes[3], ")"))
        gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)
    }
    dev.off()
}

plot_gene = function(barseq, g) {
    bin_expr = as_tibble(t(counts(barseq))) %>%
        mutate(ml_bin=barseq$ml_bin, ml_min=barseq$ml_min, ml_max=barseq$ml_max, ap=barseq$ap) %>%
        group_by(ml_bin, ml_min, ml_max, ap) %>%
        summarize_all(sum) %>%
        pivot_longer(c(-ml_bin, -ml_min, -ml_max, -ap), names_to = "gene", values_to = "expr") %>%
        mutate(expr = expr / sum(expr) * 1e4) %>%
        ungroup() %>%
        filter(gene == g) %>%
        select(ml_min, ml_max, ap, expr)
    p1 = bin_expr %>%
        ggplot(aes(fill=expr)) +
        geom_rect(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1), show.legend=FALSE) +
        ggtitle(g) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        labs(x="M-L axis",y="A-P axis") +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,NA)) 
    return(p1)
}

plot_subclass_nmf_association = function(nmf_loading) {
    # potential issue -> factors and loadings are unnormalized. Is this an issue when comparing gene loadings?
    subclass_loading = nmf_loading %>%
        group_by(subclass, nmf) %>%
        summarize(loading = mean(loading), .groups="drop")
    # compute null loadings: permute loadings across each NMF factor
    set.seed(17)
    null_loading = nmf_loading %>%
        group_by(nmf) %>%
        mutate(loading = sample(loading)) %>%
        group_by(subclass, nmf) %>%
        summarize(loading = mean(loading), .groups="drop")
    # null quantiles -> 0.007, 0.008, 0.0083, 0.009
    quantile(null_loading$loading, probs = c(0.5, 0.9, 0.95, 0.99))
    my_threshold = quantile(null_loading$loading, probs = c(0.99))
    p1 = bind_rows(list(observed = subclass_loading, null = null_loading), .id = "type") %>%
        ggplot(aes(x = loading, col = type)) +
        geom_density() +
        theme_bw() +
        geom_vline(xintercept = my_threshold, linetype = "dashed") +
        facet_wrap(~ nmf)
    p1
    ggsave("figs/spatial_patterns/gene_association/subclass_null.pdf")
    # -> almost no significant association (except maybe NMF8?)
    
    # let’s start plotting
    subclass_loading = subclass_loading %>%
        pivot_wider(subclass, names_from = "nmf", values_from = "loading") %>%
        column_to_rownames("subclass") %>%
        as.matrix()
    subclass_loading = subclass_loading[, nmf_name()]
    p1 = as_tibble(subclass_loading, rownames="subclass") %>%
        pivot_longer(-subclass, names_to = "nmf", values_to = "loading") %>%
        mutate(nmf = factor(nmf, rev(colnames(subclass_loading)))) %>%
        mutate(subclass = factor(subclass, ctx_subclass())) %>%
        ggplot(aes(x=subclass,y=nmf,fill=loading)) +
        geom_raster() +
#        geom_text(aes(label = ifelse(loading>my_threshold,"*",NA)), size=10, col="white") +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
        scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,NA))
    p1
    ggsave("figs/spatial_patterns/gene_association/subclass_loading.pdf", p1)
    p1 = as_tibble(subclass_loading, rownames="subclass") %>%
        pivot_longer(-subclass, names_to = "nmf", values_to = "loading") %>%
        mutate(subclass = factor(subclass, ctx_subclass())) %>%
        group_by(subclass) %>%
        summarize(total_loading = sum(loading)) %>%
        ggplot(aes(x=subclass,y=total_loading)) +
        geom_col() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
    p1
    ggsave("figs/spatial_patterns/gene_association/subclass_total.pdf", p1)
    # some subclasses have higher loadings overall
    # - intuition: lower subtype diversity (~ less noise) -> L2/3 IT, L6 IT, CT
    # -> switch to subclass-normalized loadings
    subclass_loading = subclass_loading / rowSums(subclass_loading)
    p1 = as_tibble(subclass_loading, rownames="subclass") %>%
        pivot_longer(-subclass, names_to = "nmf", values_to = "loading") %>%
        mutate(nmf = factor(nmf, rev(colnames(subclass_loading)))) %>%
        mutate(subclass = factor(subclass, ctx_subclass())) %>%
        ggplot(aes(x=subclass,y=nmf,fill=loading)) +
        geom_raster() +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
        scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,NA)) 
    p1
    ggsave("figs/spatial_patterns/gene_association/subclass_normalized_loading.pdf", p1)    
}


plot_cluster_nmf_spatial_association = function(nmf_tidy) {
    cluster = load_cluster_info()
    nmf_matrix = nmf_tidy %>%
        pivot_wider(bin_id, names_from = nmf, values_from = coef) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    nmf_matrix = nmf_matrix[, nmf_name()]
    cluster_tidy = cluster %>%
        filter(bin_id %in% rownames(nmf_matrix)) %>%
        pivot_wider(names_from = cell_type, values_from = n_cells, values_fill = 0) 
    cluster_matrix = as.matrix(column_to_rownames(cluster_tidy, "bin_id"))
    cluster_matrix = cluster_matrix[rownames(nmf_matrix),]
    
    # compute associations using Spearman correlation
    # the magnitude of the Spearman correlation depends on sparsity -> rescale
    overlap = cor(nmf_matrix, cluster_matrix, method="s")
    overlap = t(scale(t(overlap)))
    positive_overlap = overlap
    positive_overlap[positive_overlap<0] = 0
    to_plot = overlap[, order_rows_according_to_cols(t(positive_overlap))]
    write(colnames(to_plot), "figs/spatial_patterns/gene_association/ct_order.txt")
    p1 = as_tibble(to_plot, rownames = "nmf") %>%
        pivot_longer(-nmf, names_to = "cluster", values_to = "s_cor") %>%
        mutate(subclass = cluster_to_subclass(cluster)) %>%
        mutate(nmf = factor(nmf, rev(rownames(to_plot)))) %>%
        mutate(cluster = factor(cluster, colnames(to_plot))) %>%
        mutate(subclass = factor(subclass, ctx_subclass())) %>%
        ggplot(aes(x=cluster,y=nmf)) +
        geom_point(aes(size=s_cor, col = s_cor)) +
        facet_grid(~ subclass, scales = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=6), legend.position = "top") +
        scale_radius(range = c(0,3)) +
#        scale_size_area(max_size = 3) +
        scale_color_gradient2(low = "darkorange", mid = "gray90", high = "darkblue") +
        labs(x=NULL,y=NULL)
    p1
    ggsave("figs/spatial_patterns/gene_association/cluster_associations_space.pdf", width = 13, height = 3)
    
    # plot strongest associations
    cluster_pos = cluster_tidy %>%
        pivot_longer(-bin_id, names_to = "nmf", values_to = "coef") %>%
        add_pos()
    my_threshold = 1
    pdf("figs/spatial_patterns/gene_association/clusters_top.pdf")
    for (my_nmf in rownames(overlap)) {
        p1 = plot_nmf(nmf_tidy, my_nmf)
        top_cluster = round(sort(overlap[my_nmf,], decreasing = TRUE), 2)
        top_cluster = top_cluster[top_cluster > my_threshold]
        p2 = lapply(names(top_cluster), function(my_cluster) {
            plot_nmf(cluster_pos, my_cluster) +
                ggtitle(paste0(my_cluster, " (", top_cluster[my_cluster], ")"))
        })
        do.call(gridExtra::grid.arrange, c(list(p1), p2, list(ncol=3)))
    }
    dev.off()

    # plot all cluster patterns
    cluster_pos = cluster_tidy %>%
        pivot_longer(-bin_id, names_to = "cluster", values_to = "n_cells") %>%
        add_pos() %>%
        mutate(subclass = cluster_to_subclass()[cluster]) %>%
        mutate(cluster = factor(cluster, cell_type_order())) %>%
        mutate(subclass = factor(subclass, ctx_subclass()))
    pdf("figs/spatial_patterns/gene_association/cluster_distribution.pdf")
    for (my_subclass in ctx_subclass()) {
        keep = cluster_pos$subclass == my_subclass
        if (sum(keep) == 0) { next }
        print(plot_clusters(cluster_pos[keep,]))
    }
    dev.off()
}

load_cluster_info = function() {
    barseq = load_annotated_data()
    barseq$bin_id = barseq$slice*100 + barseq$ml
    barseq = barseq[, !(barseq$subclass %in% c("Car3", "RSP UL", "RSP DL"))]
    result = tibble(bin_id = barseq$bin_id, cell_type = barseq$cluster) %>%
        group_by(bin_id, cell_type) %>%
        tally(name = "n_cells") %>%
        ungroup()
    return(result)
}

cluster_to_subclass = function(cluster_name) {
    clusters = read_csv("analysis/labels.csv")
    result = select(clusters, cluster, subclass) %>%
        drop_na() %>%
        distinct() %>%
        deframe()
    return(result[cluster_name])
}

plot_clusters = function(cluster_tidy) {
    cluster_tidy %>%
        group_by(cluster) %>%
        mutate(expr = n_cells / max(n_cells)) %>%
        ungroup() %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=expr)) +
        geom_rect(show.legend=FALSE) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        labs(x="M-L axis",y="A-P axis") +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,1)) +
        facet_wrap(~cluster)
}

plot_cluster_nmf_gene_association = function(tidy_var) {
    # marker statistics of top genes
    ct_order = readLines("figs/spatial_patterns/gene_association/ct_order.txt")
    marker_stats = read_csv("analysis/markers/cluster_markers_by_subclass.csv") %>%
        filter(cell_type %in% ct_order)

    # extract genes significantly associated with NMF components
    gene_loading = var_tidy %>%
        group_by(gene, nmf) %>%
        summarize(var_exp = mean(var_exp), .groups="drop")
    gene_loading = gene_loading %>%
        pivot_wider(gene, names_from = "nmf", values_from = "var_exp") %>%
        column_to_rownames("gene") %>%
        as.matrix()
    gene_loading = gene_loading[, nmf_name()]
    # keep only genes that have an association > 99% null quantile with some NMF
    my_threshold = 0.254
    keep_gene = rowMaxs(gene_loading) > my_threshold & (rownames(gene_loading) != "Slc30a3")
    nmf_genes = gene_loading[keep_gene,]
    n_nmf = rowSums(gene_loading > my_threshold)
    keep_gene = n_nmf == 1 & (rownames(gene_loading) != "Slc30a3")
    nmf_genes2 = gene_loading[keep_gene,]
    nmf_genes = nmf_genes[order_rows_according_to_cols(nmf_genes),]
    nmf_genes2 = nmf_genes2[order_rows_according_to_cols(nmf_genes2),]
    to_plot = cluster_matrix
    my_genes = as_tibble(nmf_genes2, rownames="gene") %>%
        pivot_longer(-gene, names_to = "nmf", values_to = "loading") %>%
        filter(loading > my_threshold) %>%
        group_by(nmf) %>%
        slice_max(loading, n = 100) %>%
        ungroup() %>%
        select(gene, nmf, loading)
    
    # IT types
    p1 = marker_stats %>%
        mutate(cluster = factor(cell_type, ct_order)) %>%
        mutate(subclass = factor(group, ctx_subclass())) %>%
        inner_join(my_genes, by = "gene") %>%
        mutate(nmf = factor(nmf, nmf_name())) %>%
        filter(grepl("IT", subclass)) %>%
        ggplot(aes(x=cluster,y=gene,col=log2(fold_change))) +
        geom_point(aes(size=log2(fold_change_detection), fill = log2(fold_change), col=log_fdr<log(0.05)), pch=21) +
        facet_grid(nmf ~ subclass, scales = "free") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=6)) +
        scale_radius(range = c(0,3)) +
#        scale_size_area(max_size = 3) +
        scale_fill_gradient2(low = "darkorange", mid = "gray90", high = "darkblue") +
        scale_color_manual(values=c(gray(1,0),"black")) +
        labs(x=NULL,y=NULL)
    p1
    ggsave("figs/spatial_patterns/gene_association/cluster_associations_gene_it_example.pdf", height = 9)
    # non-IT types
    p1 = marker_stats %>%
        mutate(cluster = factor(cell_type, ct_order)) %>%
        mutate(subclass = factor(group, ctx_subclass())) %>%
        inner_join(my_genes, by = "gene") %>%
        mutate(nmf = factor(nmf, nmf_name())) %>%
        filter(!grepl("IT", subclass)) %>%
        ggplot(aes(x=cluster,y=gene,col=log2(fold_change))) +
        geom_point(aes(size=log2(fold_change_detection), fill = log2(fold_change), col=log_fdr<log(0.05)), pch=21) +
        facet_grid(nmf ~ subclass, scales = "free") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=6)) +
        scale_radius(range = c(0,3)) +
#        scale_size_area(max_size = 3) +
        scale_fill_gradient2(low = "darkorange", mid = "gray90", high = "darkblue") +
        scale_color_manual(values=c(gray(1,0),"black")) +
        labs(x=NULL,y=NULL)
    p1
    # -> transform into enrichment values
    
    # Hypergeometric test
    all_genes = unique(marker_stats$gene)
    marker_sets = marker_stats %>%
        filter(fold_change > 1.5 & log_fdr < log(0.05)) %>%
        select(cell_type, gene) %>%
        with(., split(.$gene, .$cell_type))
    nmf_sets = my_genes %>%
        with(., split(.$gene, .$nmf))
    ct_enrichment = gene_set_enrichment(marker_sets, nmf_sets, all_genes)
    p1 = ct_enrichment %>%
        mutate(cluster = factor(candidate_set, ct_order)) %>%
        mutate(subclass = factor(cluster_to_subclass()[candidate_set], ctx_subclass())) %>%
        mutate(nmf = factor(reference_set, rev(nmf_name()))) %>%
        mutate(odds_ratio = pmax(pmin(odds_ratio,8), 1/8)) %>%
        ggplot(aes(x=cluster,y=nmf)) +
        geom_point(aes(size=log2(odds_ratio), col = log2(odds_ratio))) +
        facet_grid(~ subclass, scales = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=6), legend.position = "top") +
        scale_radius(range = c(0,3)) +
#        scale_size_area(max_size = 3) +
        scale_color_gradient2(low = "darkorange", mid = "gray90", high = "darkblue") +
        labs(x=NULL,y=NULL)
    p1
    ggsave("figs/spatial_patterns/gene_association/cluster_associations_gene_hyper.pdf", width = 13, height = 3)
    
    # AUROC enrichment (FC or FDR ranking?)
    marker_stats = mutate(marker_stats, mlog_fdr = -log_fdr)
    my_stat = "fold_change" # mlog_fdr or fold_change
    marker_strength = marker_stats %>%
        pivot_wider(gene, names_from = "cell_type", values_from = all_of(my_stat)) %>%
        column_to_rownames("gene") %>%
        as.matrix()
    nmf_sets = my_genes %>%
        with(., split(.$gene, .$nmf)) %>%
        gene_set_to_matrix(rownames(marker_strength))
    ct_enrichment = as_tibble(compute_aurocs(marker_strength, nmf_sets), rownames="nmf") %>%
        pivot_longer(-nmf, names_to = "cell_type", values_to = "auroc")
    p1 = ct_enrichment %>%
        mutate(cluster = factor(cell_type, ct_order)) %>%
        mutate(subclass = factor(cluster_to_subclass()[cell_type], ctx_subclass())) %>%
        mutate(nmf = factor(nmf, rev(nmf_name()))) %>%
        ggplot(aes(x=cluster,y=nmf)) +
        geom_point(aes(size=auroc, col = auroc)) +
        facet_grid(~ subclass, scales = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=6), legend.position = "top") +
        scale_radius(range = c(0,3)) +
#        scale_size_area(max_size = 3) +
        scale_color_gradient2(low = "darkorange", mid = "gray90", high = "darkblue", midpoint = 0.5) +
        labs(x=NULL,y=NULL)
    p1
    ggsave(paste0("figs/spatial_patterns/gene_association/cluster_associations_gene_auroc_", my_stat ,".pdf"), width = 13, height = 3)
    
    # DE score (based on logFC and -logFDR)
    marker_stats = mutate(marker_stats, log2fc = pmin(log2(fold_change),2))
    marker_stats = mutate(marker_stats, mlog_fdr = pmin(-log_fdr, 20))
    my_stat = "log2fc" # mlog_fdr or log2fc
    de_score = marker_stats %>%
        select(cell_type, gene, de_score = .data[[my_stat]]) %>%
        inner_join(my_genes) %>%
        group_by(cell_type, nmf) %>%
        summarize(de_score = mean(de_score), .groups="drop")
    p1 = de_score %>%
        mutate(cluster = factor(cell_type, ct_order)) %>%
        mutate(subclass = factor(cluster_to_subclass()[cell_type], ctx_subclass())) %>%
        mutate(nmf = factor(nmf, rev(nmf_name()))) %>%
        ggplot(aes(x=cluster,y=nmf)) +
        geom_point(aes(size=de_score, col = de_score)) +
        facet_grid(~ subclass, scales = "free_x") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=6), legend.position = "top") +
        scale_radius(range = c(0,3)) +
#        scale_size_area(max_size = 3) +
        scale_color_gradient2(low = "darkorange", mid = "gray90", high = "darkblue") +
        labs(x=NULL,y=NULL)
    p1
    ggsave(paste0("figs/spatial_patterns/gene_association/cluster_associations_gene_de_score_", my_stat ,".pdf"), width = 13, height = 3)
    
    # show one example for a specific subclass? side-by-side NMF / subtype pattern?
    # -> maybe L5 IT because there aren’t many subtypes?
}

if (sys.nframe() == 0) {
    main()
}
