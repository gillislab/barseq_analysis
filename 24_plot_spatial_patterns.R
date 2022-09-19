# -*- coding: utf-8 -*-

library(tidyverse)
library(ica)
library(NMF)
source("dataset.R")


main = function() {
    my_genes = c("Nnat", "Slc24a2", "C1ql3","Lpp", "Ctgf", "Lmo4", "Rcan2","Tenm3")
    plot_nmf_scaled(my_genes)
    plot_nmf(my_genes)
    plot_pca(my_genes)
    plot_ica(my_genes)
} 

plot_nmf_scaled = function(my_genes) {
    # load nmf results 
    my_nmf = readRDS("analysis/spatial_patterns/spatial_nmf_scaled.rds")
    w = basis(my_nmf)
    colnames(w) = paste0("NMF", 1:ncol(w))
    h = coef(my_nmf)
    rownames(h) = paste0("NMF", 1:ncol(w))
    nmf_tidy = read_csv("analysis/spatial_patterns/spatial_nmf_scaled.csv")
    nmf_loading = read_csv("analysis/spatial_patterns/spatial_nmf_scaled_loadings.csv")

    # variance explained?
    pseudobulk_matrix = load_pseudobulk_matrix()
    scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=sqrt(colSums(pseudobulk_matrix**2))) # l2-normalization
    #total_var = sum(colVars(scaled_pseudobulk))
    total_var = sum(scaled_pseudobulk**2)
    n_f = 1:ncol(w)
    var_exp = sapply(n_f, function(i) {
        #residual_var = sum(colVars(scaled_pseudobulk - w[,1:i,drop=FALSE]%*%h[1:i,,drop=FALSE]))
        residual_var = sum((scaled_pseudobulk - w[,1:i,drop=FALSE]%*%h[1:i,,drop=FALSE])**2)
        (total_var - residual_var) / total_var
    })
    tibble(n_f, var_exp) %>%
        ggplot(aes(x=n_f, y=100*var_exp)) +
        geom_col(fill="skyblue3") +
        theme_bw(base_size = 20) +
        labs(x="Number of NMF factors", y="Cumulative variance explained")
    ggsave("figs/spatial_patterns/nmf_scaled/var_exp.pdf")
    # main factors
    my_mfs = paste0("NMF", 1:ncol(w))
    nmf_tidy %>%
        select(bin_id, ap, ml_min, ml_max, all_of(my_mfs)) %>%
        pivot_longer(c(-bin_id,-ap,-ml_min,-ml_max), names_to = "nmf", values_to = "expr") %>%
        group_by(nmf) %>%
        mutate(expr = expr / max(expr)) %>%
        ungroup() %>%
        mutate(nmf = factor(nmf, my_mfs)) %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=expr)) +
        geom_rect() +
        facet_wrap(~ nmf) +
        theme_classic() +
        labs(x="M-L axis",y="A-P axis") +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
              legend.position = c(7/8,1/6)) +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,NA))
    ggsave("figs/spatial_patterns/nmf_scaled/nmf_1_10.pdf")
    nmf_loading %>%
        filter(gene %in% my_genes) %>%
        ggplot(aes(nmf,loading)) +
        geom_boxplot() +
        facet_wrap(~ gene)
    ggsave("figs/spatial_patterns/nmf_scaled/gene_boxplot.pdf")
    gene_pc_association = nmf_loading %>%
        group_by(gene, nmf) %>%
        summarize(t = median(loading), .groups="drop")
    gene_pc_association %>%
        filter(gene %in% my_genes) %>%
        mutate(gene = factor(gene, rev(my_genes))) %>%
        ggplot(aes(x = nmf, y=gene, fill=t)) +
        geom_raster() +
        scale_fill_gradient(low="white", high="darkblue")
    ggsave("figs/spatial_patterns/nmf_scaled/gene_association.pdf")
    association = nmf_loading %>%
        group_by(subclass, nmf) %>%
        summarize(t = sum(loading**2), .groups="drop")
    association %>%
        ggplot(aes(x = nmf, y=subclass, fill=t)) +
        geom_raster() +
        scale_fill_gradient(low="white", high="darkblue")
    ggsave("figs/spatial_patterns/nmf_scaled/ct_association.pdf")
    # how much do the top NMFs explain in each cell type?
    # can we show that these NMFs correspond to "generic"/"cell type unrelated" patterns?
    scaling_weights = colSds(pseudobulk_matrix)
    total_var = as_tibble(t(scaled_pseudobulk)*scaling_weights, rownames="var") %>%
        separate(var, c("gene", "ct"), "\\|") %>%
        pivot_longer(c(-gene,-ct), names_to = "bin_id", values_to = "expr") %>%
        group_by(ct) %>%
        summarize(total_var = sum(expr**2))
    nmf_components = lapply(set_names(1:ncol(w)), function(i) { w[,i,drop=FALSE]%*%h[i,,drop=FALSE]})
    nmf_var = bind_rows(lapply(nmf_components, function(c) {
        as_tibble(t(c)*scaling_weights, rownames="var") %>%
        separate(var, c("gene", "ct"), "\\|") %>%
        pivot_longer(c(-gene,-ct), names_to = "bin_id", values_to = "expr") %>%
        group_by(ct) %>%
        summarize(nmf_var = sum(expr**2), .groups="drop")
    }), .id="nmf")
    var_exp = inner_join(total_var, nmf_var) %>%
        mutate(var_exp = nmf_var / total_var)
    var_exp %>%
        ggplot(aes(x = as.integer(nmf), y= ct, fill=var_exp*100)) +
        geom_raster() +
        scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,NA))
    # these loadings aren’t very interesting, should we re-project the data?
    anova_res = scale(pseudobulk_matrix, scale=FALSE)
    total_var = enframe(colSums(anova_res**2), "var", "total_var") %>%
        separate(var, c("gene", "ct"), "\\|") %>%
        group_by(ct) %>%
        summarize(total_var = sum(total_var), .groups="drop")
    scaled_w = scale(w, center=FALSE, scale=sqrt(colSums(w**2)))
    gplots::heatmap.2(crossprod(scaled_w))
    crossprod(scaled_w)
    loadings = crossprod(scaled_w, anova_res)
    nmf_var = as_tibble(t(loadings), rownames = "var") %>%
        separate(var, c("gene", "ct"), "\\|") %>%
        pivot_longer(c(-gene, -ct), names_to = "nmf", values_to = "loading") %>%
        group_by(ct, nmf) %>%
        summarize(nmf_var = sum(loading**2), .groups="drop")
    var_exp = inner_join(total_var, nmf_var) %>%
        mutate(var_exp = nmf_var / total_var)
    var_exp %>%
        ggplot(aes(x = nmf, y= ct, fill=var_exp*100)) +
        geom_raster() +
        scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,NA))
    ggsave("figs/spatial_patterns/nmf_scaled/ct_association_centered.pdf")
}

load_pseudobulk_matrix = function() {
    pseudobulks = readRDS("data/subclass_pseudobulks.rds")
    pseudobulks = pseudobulks[!rownames(pseudobulks) %in% c("Slc17a7", "Gad1"),]
    cpm(pseudobulks) = convert_to_cpm(counts(pseudobulks), total_counts = median(colSums(counts(pseudobulks))))
    subclass_pseudobulks = as_tibble(t(cpm(pseudobulks))) %>%
        mutate(subclass = pseudobulks$subclass, bin_id = pseudobulks$bin_id) %>%
        pivot_longer(c(-subclass, -bin_id), names_to = "gene", values_to = "expr")
    pseudobulk_matrix = subclass_pseudobulks %>%
        filter(!(subclass %in% c("Car3", "RSP UL", "RSP DL"))) %>%
        mutate(var = paste(gene,subclass,sep="|")) %>%
        pivot_wider(bin_id, names_from = var, values_from = expr) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    return(pseudobulk_matrix)
}

plot_nmf = function(my_genes) {
    # load results
    my_nmf = readRDS("analysis/spatial_patterns/spatial_nmf.rds")
    w = basis(my_nmf)
    colnames(w) = paste0("NMF", 1:ncol(w))
    h = coef(my_nmf)
    rownames(h) = paste0("NMF", 1:ncol(w))
    nmf_tidy = read_csv("analysis/spatial_patterns/spatial_nmf.csv")
    nmf_loading = read_csv("analysis/spatial_patterns/spatial_nmf_loadings.csv")
    
    # plot all factors
    my_mfs = paste0("NMF", 1:ncol(w))
    nmf_tidy %>%
        select(bin_id, ap, ml_min, ml_max, all_of(my_mfs)) %>%
        pivot_longer(c(-bin_id,-ap,-ml_min,-ml_max), names_to = "nmf", values_to = "expr") %>%
        mutate(nmf = factor(nmf, my_mfs)) %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=expr)) +
        geom_rect() +
        facet_wrap(~ nmf) +
        theme_classic() +
        scale_fill_gradient(low = "white", high = "darkblue", limits=c(0,NA))
    ggsave("figs/spatial_patterns/nmf/nmf_1_10.pdf")
    nmf_loading %>%
        filter(gene %in% my_genes) %>%
        ggplot(aes(nmf,loading)) +
        geom_boxplot() +
        facet_wrap(~ gene)
    ggsave("figs/spatial_patterns/nmf/gene_boxplot.pdf")
    gene_pc_association = nmf_loading %>%
        group_by(gene, nmf) %>%
        summarize(t = median(loading), .groups="drop")
    gene_pc_association %>%
        filter(gene %in% my_genes) %>%
        mutate(gene = factor(gene, rev(my_genes))) %>%
        ggplot(aes(x = nmf, y=gene, fill=t)) +
        geom_raster() +
        scale_fill_gradient(low="white", high="darkblue")
    ggsave("figs/spatial_patterns/nmf/gene_association.pdf")
}

plot_pca = function(my_genes) {
    my_pca = readRDS("analysis/spatial_patterns/spatial_pcs.rds")
    pc_tidy = read_csv("analysis/spatial_patterns/spatial_pcs.csv")
    
    # variance explained
    pc_s = summary(my_pca)
    tibble(var_exp=100*pc_s$importance[2,1:20], pc=1:20) %>%
        ggplot(aes(x=pc,y=var_exp)) +
        geom_col(fill="skyblue3") +
        lims(y=c(0,NA)) +
        labs(x="PC", y="Variance explained (%)")
    ggsave("figs/spatial_patterns/pca/pc_weight.pdf")
    tibble(var_exp=100*pc_s$importance[2,]) %>%
        ggplot(aes(x=var_exp)) +
        geom_histogram(fill="skyblue3", bins = 50) +
        labs(x="Variance explained (%)") +
        scale_x_log10() +
        geom_vline(xintercept = 0.7, linetype="dashed")
    ggsave("figs/spatial_patterns/pca/pc_weight_dist.pdf")
    
    # PC1=batch effect, PC2=M and L, PC3=M-L gradient, PC7=A-P gradient
    my_pcs = paste0("PC", 1:9)
    pc_tidy %>%
        select(bin_id, ap, ml_min, ml_max, all_of(my_pcs)) %>%
        pivot_longer(c(-bin_id,-ap,-ml_min,-ml_max), names_to = "pc", values_to = "expr") %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=expr)) +
        geom_rect() +
        facet_wrap(~ pc) +
        theme_classic() +
        scale_fill_gradient2(low = "darkorange", mid = "white", high = "darkblue",
                             limits = c(-4,4), oob = scales::squish)
    ggsave("figs/spatial_patterns/pca/pc_1_9.pdf")

    # loadings?
    # anova_res = my_pca$x %*% t(my_pca$rotation)
    # are weighted loadings a better indication of key patterns?
    pc_loading = as_tibble(sweep(my_pca$rotation[,1:9], 2, my_pca$sdev[1:9], "*"), rownames = "var") %>%
        pivot_longer(-var, names_to = "pc", values_to = "loading") %>%
        separate(var, sep = "\\|", into = c("gene","subclass"))
    # association PC - gene
    gene_loading = pc_loading %>%
        group_by(gene, pc) %>%
        summarize(loading = mean(loading), .groups="drop") %>%
        pivot_wider(c(gene), names_from = pc, values_from = loading)
    gene_loading %>%
        ggplot(aes(PC1, PC2, label=gene)) +
        geom_text()
    ggsave("figs/spatial_patterns/pca/gene_pc_12.pdf")
    gene_loading %>%
        ggplot(aes(PC2, PC3, label=gene)) +
        geom_text()
    ggsave("figs/spatial_patterns/pca/gene_pc_23.pdf")
    pc_loading %>%
        filter(gene %in% my_genes) %>%
        ggplot(aes(pc,loading)) +
        geom_boxplot() +
        facet_wrap(~ gene)
    ggsave("figs/spatial_patterns/pca/gene_pc_boxplot.pdf")
    gene_pc_association = pc_loading %>%
        group_by(gene, pc) %>%
        summarize(t = mean(loading) / sd(loading) * sqrt(n()), .groups="drop")
    gene_pc_association %>%
        filter(gene %in% my_genes) %>%
        mutate(gene = factor(gene, rev(my_genes))) %>%
        ggplot(aes(x = pc, y=gene, fill=t)) +
        geom_raster() +
        scale_fill_gradient2(low="darkorange", mid="white", high="darkblue", midpoint=0)
    ggsave("figs/spatial_patterns/pca/gene_pc_association.pdf")
    # association PC - CT
    ct_loading = pc_loading %>%
        group_by(subclass, pc) %>%
        summarize(loading = mean(loading), .groups="drop") %>%
        pivot_wider(c(subclass), names_from = pc, values_from = loading)
    ct_loading %>%
        ggplot(aes(PC1, PC2, label=subclass)) +
        geom_text()
    ct_loading %>%
        ggplot(aes(PC2, PC3, label=subclass)) +
        geom_text()
    pc_ct_association = pc_loading %>%
        group_by(subclass, pc) %>%
        summarize(t = mean(loading) / sd(loading) * sqrt(n()), .groups="drop")
    pc_ct_association %>%
        ggplot(aes(x = pc, y=subclass, fill=t)) +
        geom_raster() +
        scale_fill_gradient2(low="darkorange", mid="white", high="darkblue", midpoint=0)
    ggsave("figs/spatial_patterns/pca/ct_pc_association.pdf")
    pc_loading %>%
        ggplot(aes(pc,loading)) +
        geom_boxplot() +
        facet_wrap(~ subclass)
    ggsave("figs/spatial_patterns/pca/ct_pc_boxplot.pdf")
} 
    
plot_ica = function(my_genes) {
    my_ica = readRDS("analysis/spatial_patterns/spatial_ics.rds")
    ic_tidy = read_csv("analysis/spatial_patterns/spatial_ics.csv")
    ic_loading = read_csv("analysis/spatial_patterns/spatial_ic_loadings.csv")

    # variance explained
    tibble(var_exp=100*my_ica$vafs[1:20], ic=1:20) %>%
        ggplot(aes(x=ic,y=var_exp)) +
        geom_col(fill="skyblue3") +
        lims(y=c(0,NA)) +
        labs(x="IC", y="Variance explained (%)")
    ggsave("figs/spatial_patterns/ica/ic_weight.pdf")
    tibble(var_exp=100*my_ica$vafs) %>%
        ggplot(aes(x=var_exp)) +
        geom_histogram(fill="skyblue3", bins = 20) +
        labs(x="Variance explained (%)") +
        scale_x_log10() +
        geom_vline(xintercept = 0.7, linetype="dashed")
    ggsave("figs/spatial_patterns/ica/ic_weight_dist.pdf")
    my_ics = paste0("IC", 1:16)
    ic_tidy %>%
        select(bin_id, ap, ml_min, ml_max, all_of(my_ics)) %>%
        pivot_longer(c(-bin_id,-ap,-ml_min,-ml_max), names_to = "ic", values_to = "expr") %>%
        mutate(ic = factor(ic, my_ics)) %>%
        ggplot(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1,fill=expr)) +
        geom_rect() +
        facet_wrap(~ ic) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) +
        scale_fill_gradient2(low = "darkorange", mid = "white", high = "darkblue",
                             limits = c(-6,6), oob = scales::squish)
    ggsave("figs/spatial_patterns/ica/ic_1_16.pdf")
    # ICA is way too sensitive to outliers!
}

if (sys.nframe() == 0) {
    main()
}
