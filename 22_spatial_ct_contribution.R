# -*- coding: utf-8 -*-

library(tidyverse)
library(SingleCellExperiment)
source("dataset.R")


main = function() {
    # gene-wise computation
    # - ct variability: in each cubelet, what fraction of variance is accounted by ct?
    # - spatial variability: in each ct, what fraction of variance is accounted by space?
    # fraction of variance -> ANOVA
    # - ct variability: gene_i ~ ct_i + epsilon -> gene is important for cell typing
    # - spatial variability: gene_i ~ bin_i + epsilon -> gene has a lot of spatial variability
    
    # subclass level
    subclass_pseudobulks = readRDS("data/subclass_pseudobulks.rds")
    subclass_pseudobulks = subclass_pseudobulks[!rownames(subclass_pseudobulks) %in% c("Slc17a7", "Gad1"),]
    subclass_pseudobulks = subclass_pseudobulks[, !(subclass_pseudobulks$subclass %in% c("Car3", "RSP UL", "RSP DL"))]
    cpm(subclass_pseudobulks) = convert_to_cpm(counts(subclass_pseudobulks),
                                               total_counts = median(colSums(counts(subclass_pseudobulks))))
    tidy_pseudobulks = as_tibble(t(cpm(subclass_pseudobulks))) %>%
        mutate(subclass = subclass_pseudobulks$subclass, bin_id = subclass_pseudobulks$bin_id) %>%
        pivot_longer(c(-subclass, -bin_id), names_to = "gene", values_to = "expr")
    subclass_anova = tidy_pseudobulks %>%
        group_by(gene) %>%
        nest() %>%
        summarize(aov = lapply(data, function(df) { bind_rows(list(
                cell_type = one_way_anova(df, "subclass"),
                space = one_way_anova(df, "bin_id")
            ), .id = "factor")})) %>%
        unnest(aov)
    p1 = subclass_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        ggplot(aes(x=cell_type,y=space,label=gene)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw(base_size=20) +
        labs(x="Var. explained (H2 types)",y="Var. explained (space)")
    p1
    ggsave("figs/spatial_patterns/ct_vs_spatial/subclass.pdf", p1)
    
    my_genes = c("Ctgf", "Slc24a2", "Tenm3", "Nnat")
    pdf("figs/spatial_patterns/ct_vs_spatial/subclass_examples.pdf", width=5)
    for (g in my_genes) {
        plot_gene(subclass_pseudobulks, g)
    }
    dev.off()

    # cluster level
    cluster_pseudobulks = readRDS("data/cluster_pseudobulks.rds")
    cluster_pseudobulks = cluster_pseudobulks[!rownames(cluster_pseudobulks) %in% c("Slc17a7", "Gad1"),]
    cluster_pseudobulks = cluster_pseudobulks[, !(cluster_pseudobulks$subclass %in% c("Car3", "RSP UL", "RSP DL"))]
    cpm(cluster_pseudobulks) = convert_to_cpm(counts(cluster_pseudobulks),
                                              total_counts = median(colSums(counts(cluster_pseudobulks))))
    tidy_pseudobulks = as_tibble(t(cpm(cluster_pseudobulks))) %>%
        mutate(subclass = cluster_pseudobulks$subclass, cluster = cluster_pseudobulks$cluster,
               bin_id = cluster_pseudobulks$bin_id) %>%
        pivot_longer(c(-subclass, -cluster, -bin_id), names_to = "gene", values_to = "expr")
    cluster_anova = tidy_pseudobulks %>%
        group_by(gene) %>%
        nest() %>%
        summarize(aov = lapply(data, function(df) { bind_rows(list(
                cell_type = one_way_anova(df, "cluster"),
                subclass = one_way_anova(df, "subclass"),
                space = one_way_anova(df, "bin_id")
            ), .id = "factor")})) %>%
        unnest(aov)
    write_csv(cluster_anova, "analysis/anova/cluster_anova.csv")
    p1 = cluster_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        ggplot(aes(x=cell_type,y=space,label=gene)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw(base_size=20) +
        labs(x="Var. explained (H3 types)",y="Var. explained (space)")
    p1
    ggsave("figs/spatial_patterns/ct_vs_spatial/cluster.pdf", p1)
    p1 = cluster_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        ggplot(aes(x=subclass,y=space,label=gene)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw(base_size=20) +
        labs(x="Var. explained (H2 types)",y="Var. explained (space)")
    p1
    ggsave("figs/spatial_patterns/ct_vs_spatial/subclass_by_cluster.pdf", p1)
    p1 = cluster_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        ggplot(aes(subclass,cell_type,label=gene)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw(base_size=20) +
        labs(x="Var. explained (H2 types)",y="Var. explained (H3 types)")
    ggsave("figs/spatial_patterns/ct_vs_spatial/cluster_delta.pdf")
    p1 = cluster_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        mutate(delta_var = cell_type-subclass) %>%
        ggplot(aes(delta_var,space,label=gene)) +
        geom_text() +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        theme_bw(base_size=20) +
        labs(x="Delta(var. exp.) (H3 vs H2 types)",y="Var. explained (space)")
    p1
    cluster_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        mutate(delta_var = cell_type-subclass) %>%
        with(., cor(.$delta_var, .$space, method="s")) # 0.23
    ggsave("figs/spatial_patterns/ct_vs_spatial/space_vs_cluster_delta.pdf", p1)
    p1 = cluster_anova %>%
        pivot_wider(gene, names_from = "factor", values_from = "var_exp") %>%
        mutate(delta_var = cell_type-subclass) %>%
        filter(delta_var>0.1) %>%
        mutate(gene = fct_reorder(gene, delta_var,.desc = TRUE)) %>%
        ggplot(aes(gene, delta_var,label=round(100*delta_var))) +
        geom_col(show.legend = FALSE, fill="skyblue3") +
#        geom_text(nudge_y = 0.01) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5, size=12)) +
        labs(x=NULL,y="Increase in variance explained") +
        scale_y_continuous(labels = scales::percent)
    p1
    ggsave("figs/spatial_patterns/ct_vs_spatial/cluster_delta_summary.pdf", p1)
}

one_way_anova = function(df, factor_name) {
    result = df %>%
        mutate(Y_mean = mean(expr)) %>%
        group_by(.data[[factor_name]]) %>%
        mutate(Y_fac = mean(expr)) %>%
        ungroup() %>%
        mutate(alpha = (Y_fac - Y_mean)**2,
               eps = (expr - Y_fac)**2) %>%
        summarize(N = n(), K = length(unique(.data[[factor_name]])),
                  SS_fac = sum(alpha), SS_eps = sum(eps),
                  MSS_fac = sum(alpha)/(K-1), MSS_eps = sum(eps)/(N-K),
                  var_exp = SS_fac / (SS_fac + SS_eps), .groups="drop") %>%
        mutate(F = MSS_fac / MSS_eps,
               pval = pf(F, K-1, N-K, lower.tail = FALSE))
    return(result)
}

plot_gene = function(barseq, g) {
    bin_expr = as_tibble(t(counts(barseq))) %>%
        mutate(ml_bin=barseq$ml_bin, ml_min=barseq$ml_min, ml_max=barseq$ml_max, ap=barseq$ap) %>%
        group_by(ml_bin, ml_min, ml_max, ap) %>%
        summarize_all(sum) %>%
        pivot_longer(c(-ml_bin, -ml_min, -ml_max, -ap), names_to = "gene", values_to = "expr") %>%
        mutate(expr = sqrt(expr / sum(expr) * 1e4)) %>%
        ungroup() %>%
        filter(gene == g) %>%
        select(ml_min, ml_max, ap, expr)
    subclass_expr = tibble(ml_min=barseq$ml_min, ml_max=barseq$ml_max, ml_bin=barseq$ml_bin, ap=barseq$ap,
                       expr=sqrt(cpm(barseq)[g,]), label=barseq$subclass)
    p1 = bin_expr %>%
        ggplot(aes(fill=expr)) +
        geom_rect(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1), show.legend=FALSE) +
        ggtitle(g) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,NA)) +
        labs(x="M-L axis",y="A-P axis")
    p2 = subclass_expr %>%
        ggplot(aes(fill=expr)) +
        geom_rect(aes(xmin=ml_min,xmax=ml_max,ymin=ap,ymax=ap+1), show.legend=FALSE) +
        ggtitle(g) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        scale_fill_gradient(low = "gray90", high = "darkblue", limits=c(0,NA)) +
        labs(x="M-L axis",y=NULL) +
        facet_wrap(~label)
    p3 = subclass_expr %>%
        group_by(ap, label) %>%
        summarize(expr = mean(expr), n = n(), .groups="drop") %>%
        ggplot(aes(ap,expr,col=label)) +
        geom_line(aes(group=label), show.legend=FALSE) +
        ggtitle(g) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        labs(x="A-P axis",y="Expression")
    p4 = subclass_expr %>%
        group_by(ml_bin, label) %>%
        summarize(expr = mean(expr), n = n(), .groups="drop") %>%
        ggplot(aes(ml_bin,expr,col=label)) +
        geom_line(aes(group=label), show.legend=FALSE) +
        ggtitle(g) +
        theme_bw() +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
        labs(x="M-L axis",y=NULL)
    gridExtra::grid.arrange(p1,p2,p3,p4)
}

if (sys.nframe() == 0) {
    main()
}
