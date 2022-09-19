# -*- coding: utf-8 -*-

library(tidyverse)
library(SingleCellExperiment)
library(NMF)
source("dataset.R")


main = function() {
    barseq = load_annotated_data()
    pseudobulks = make_pseudobulks(barseq)
    pseudobulks = pseudobulks[!rownames(pseudobulks) %in% c("Slc17a7", "Gad1"),]
    pseudobulks = pseudobulks[!grepl("unused", rownames(pseudobulks)),]
    cpm(pseudobulks) = convert_to_cpm(counts(pseudobulks), total_counts = median(colSums(counts(pseudobulks))))
    cluster_patterns = tibble(ap = barseq$slice, ml = barseq$ml_axis, cluster = barseq$cluster, subclass = barseq$subclass) %>%
        group_by(ap, ml, subclass, cluster) %>%
        tally() %>%
        group_by(ap,ml) %>%
        mutate(cubelet_prop = n / sum(n)) %>%
        group_by(ap,ml,subclass) %>%
        mutate(subclass_prop = n / sum(n)) %>%
        ungroup()

    plot_ct_patterns(cluster_patterns)
    # NMF3 and NMF7 play prominent roles in nearly all cell types -> why is that? are they gradient of each other?
    # first: which genes are driving these patterns?
    nmf_tidy = read_csv("analysis/spatial_nmf_scaled.csv")
    nmf_loading = read_csv("analysis/spatial_nmf_scaled_loadings.csv")
    gene_loading = nmf_loading %>%
        group_by(gene, nmf) %>%
        summarize(loading = mean(loading), .groups="drop") %>%
        pivot_wider(c(gene), names_from = nmf, values_from = loading)
    gene_loading %>%
        ggplot(aes(NMF3, NMF7, label=gene)) +
        geom_text()
    # recompute loadings based on overall expression patterns?
    w = basis(readRDS("analysis/spatial_nmf_scaled.rds"))
    colnames(w) = paste0("NMF", 1:ncol(w))
    scaled_w = scale(w, center=FALSE, scale=sqrt(colSums(w**2)))
    expression_matrix = as_tibble(t(as.matrix(counts(barseq)))) %>%
        mutate(bin_id = barseq$slice*100 + barseq$ml_axis_legacy) %>%
        pivot_longer(-bin_id, names_to = "gene", values_to = "expr") %>%
        group_by(bin_id,gene) %>%
        summarize(expr = mean(expr), .groups="drop") %>%
        pivot_wider(bin_id, names_from = "gene", values_from = "expr") %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    expression_matrix = expression_matrix[rownames(scaled_w),]
    loadings = cor(scaled_w, expression_matrix)
    as_tibble(t(loadings), rownames="gene") %>%
        ggplot(aes(x=NMF3,y=NMF7,label=gene)) +
        geom_text()
    ggsave("figs/coexp/nmf_scaled_gene_37.pdf")
    plot_gene(barseq, "Timp2", "subclass")
    plot_gene(barseq, "Tenm3", "subclass")
    plot_gene(barseq, "Lmo4", "subclass")
    plot_gene(barseq, "Ncald", "subclass")
    # NMF3
    plot_gene(barseq, "Rcan2", "subclass")
    plot_gene(barseq, "Syt2", "subclass")
    plot_gene(barseq, "Lmo4", "subclass")
    plot_gene(barseq, "Slc24a2", "subclass")
    # same thing with cell type proportions?
    cluster_to_subclass = deframe(distinct(tibble(cluster = barseq$cluster, subclass=barseq$subclass)))
    prop_matrix = tibble(bin_id = barseq$slice*100+barseq$ml_axis_legacy,
                         cluster = barseq$cluster, subclass = barseq$subclass) %>%
        group_by(bin_id, subclass, cluster) %>%
        tally() %>%
        group_by(bin_id,subclass) %>%
        mutate(subclass_prop = n / sum(n)) %>%
        ungroup() %>%
        pivot_wider(bin_id, names_from = cluster, values_from = subclass_prop, values_fill = 0) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    prop_matrix = prop_matrix[rownames(scaled_w),]
    loadings = cor(scaled_w, prop_matrix)
    as_tibble(t(loadings), rownames="cluster") %>%
        ggplot(aes(x=NMF3,y=NMF7,label=cluster)) +
        geom_text()
    ggsave("figs/coexp/nmf_scaled_cluster_37.pdf")   
    pdf("figs/ct_patterns/nmf_scaled_cluster_heatmap.pdf")
    for (my_subclass in unique(cluster_to_subclass)) {
        to_plot = as_tibble(t(loadings), rownames="cluster") %>%
            filter(cluster_to_subclass[cluster] %in% my_subclass) %>%
            pivot_longer(-cluster, names_to = "nmf", values_to = "correlation") %>%
            ggplot(aes(x=nmf,y=cluster,fill=correlation)) +
            geom_raster() +
            ggtitle(my_subclass) +
            scale_fill_gradient2(low = "darkorange", mid = "white", high = "darkblue")
        print(to_plot)
    }
    dev.off()
    my_subclass = "L4/5 IT"
    cluster_a = "L4/5 IT ML-A"
    cluster_b = "L4/5 IT UL"
    my_pseudobulks = pseudobulks[, pseudobulks$cluster %in% c(cluster_a, cluster_b)]
    expr_diff = as_tibble(t(log1p(cpm(my_pseudobulks)))) %>%
        mutate(cluster = my_pseudobulks$cluster, bin_id = my_pseudobulks$bin_id) %>%
        pivot_longer(c(-cluster, -bin_id), names_to = "gene", values_to = "expr") %>%
        pivot_wider(c(gene, bin_id), names_from = "cluster", values_from = "expr") %>%
        drop_na() %>%
        mutate(delta_expr = .data[[cluster_a]]-.data[[cluster_b]]) %>%
        group_by(gene) %>%
        summarize(auroc = mean(delta_expr>0) + mean(delta_expr==0)/2, .groups="drop")
    expr_diff %>% arrange(auroc)
    pdf("figs/ct_patterns/L45_AvsUL.pdf")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Lamp5")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Cbln2")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Zmat4")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Nrsn1")
    dev.off()
    cluster_a = "L4/5 IT ML-A"
    cluster_b = "L4/5 IT ML-P"
    my_pseudobulks = pseudobulks[, pseudobulks$cluster %in% c(cluster_a, cluster_b)]
    expr_diff = as_tibble(t(log1p(cpm(my_pseudobulks)))) %>%
        mutate(cluster = my_pseudobulks$cluster, bin_id = my_pseudobulks$bin_id) %>%
        pivot_longer(c(-cluster, -bin_id), names_to = "gene", values_to = "expr") %>%
        pivot_wider(c(gene, bin_id), names_from = "cluster", values_from = "expr") %>%
        drop_na() %>%
        mutate(delta_expr = .data[[cluster_a]]-.data[[cluster_b]]) %>%
        group_by(gene) %>%
        summarize(auroc = mean(delta_expr>0) + mean(delta_expr==0)/2, .groups="drop")
    expr_diff %>% arrange(auroc)
    pdf("figs/ct_patterns/L45_AvsP.pdf")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Nrsn1")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Ncald")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Rcan2")
    plot_gene(barseq[, barseq$subclass %in% my_subclass], "Camk4")
    dev.off()
}

old_pc = function() {
    pc_tidy = read_csv("analysis/spatial_pcs.csv")
    pcs = pc_tidy %>%
        select(-ap, -ml) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    bin_expr = as_tibble(t(cpm(pseudobulks))) %>%
        mutate(bin_id = pseudobulks$bin_id, ct = paste(pseudobulks$subclass,pseudobulks$cluster,sep="|")) %>%
        pivot_longer(c(-bin_id, -ct), names_to = "gene", values_to = "expr") %>%
        mutate(gene_ct = paste(gene, ct, sep="|")) %>%
        pivot_wider(c(bin_id), names_from = gene_ct, values_from = expr, values_fill = 0) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    pc_expr = crossprod(pcs[,1:9], bin_expr[rownames(pcs),]) %>%
        as_tibble(rownames = "pc") %>%
        pivot_longer(-pc, names_to = "gene_ct", values_to = "expr") %>%
        separate("gene_ct", c("gene", "subclass", "cluster"), "\\|")
    my_subclass = "CT"
    subclass_expr = filter(pc_expr, subclass == my_subclass)
    ct_diff = inner_join(subclass_expr, subclass_expr, by=c("pc", "gene", "subclass"), suffix=c("_a","_b")) %>%
        filter(cluster_b > cluster_a) %>%
        mutate(delta_expr = expr_b - expr_a) %>%
        select(-expr_a, -expr_b)
    ct_diff %>%
        mutate(ct = paste(cluster_a, cluster_b, sep="|")) %>%
        group_by(ct, pc) %>%
        summarize(delta_expr = mean(abs(delta_expr)), .groups="drop") %>%
        ggplot(aes(pc, delta_expr)) +
        geom_boxplot() 
#        geom_line(aes(group=ct, col=ct), show.legend = FALSE)
    ct_diff %>%
        mutate(ct = paste(cluster_a, cluster_b, sep="|")) %>%
        group_by(ct, pc) %>%
        summarize(delta_expr = mean(abs(delta_expr)), .groups="drop") %>%
        ggplot(aes(pc, ct, fill=delta_expr)) +
        geom_raster()
    ct_a = "CT CTX A"
    ct_b = "CT CTX P/RSP"
    ct_diff %>%
        filter(cluster_a == ct_a & cluster_b == ct_b) %>%
        group_by(pc) %>%
        summarize(delta_expr = mean(abs(delta_expr)), .groups="drop") %>%
        ggplot(aes(pc, delta_expr)) +
        geom_col()
    pc1_genes = ct_diff %>%
        filter(cluster_a == ct_a & cluster_b == ct_b & pc =="PC1") %>%
        arrange(desc(abs(delta_expr))) %>%
        pull(gene)
    sub_barseq = barseq[, barseq$cluster %in% c(ct_a, ct_b)]
    plot_gene(sub_barseq, pc1_genes[1])
    
    # how much do CT proportions correlate with PCs?
    barseq$bin_id = barseq$slice*100 + barseq$ml_axis
    my_subclass = "CT"
    ct_prop = tibble(subclass = barseq$subclass, cluster = barseq$cluster, bin_id = barseq$bin_id) %>%
        filter(subclass == my_subclass) %>%
        group_by(bin_id, cluster) %>%
        tally() %>%
        ungroup() %>%
        pivot_wider(bin_id, names_from = "cluster", values_from = "n", values_fill =) %>%
        column_to_rownames("bin_id") %>%
        as.matrix()
    ct_prop = ct_prop / rowSums(ct_prop)
    ct_pc_association = crossprod(ct_prop[rownames(pcs),], pcs[,1:9] / 539)
    as_tibble(ct_pc_association, rownames = "cluster") %>%
        pivot_longer(-cluster, names_to = "pc", values_to = "loading") %>%
        ggplot(aes(x=cluster, y=pc, fill=loading)) +
        geom_raster() +
        scale_fill_gradient2(low="darkorange", mid="white", high="darkblue", midpoint = 0)
    ggsave("figs/coexp/ct_example_heatmap.pdf")
    ct_prop_tidy = as_tibble(ct_prop, rownames = "bin_id") %>%
        pivot_longer(-bin_id, names_to = "cell_type", values_to = "prop") %>%
        mutate(ap = round(as.numeric(bin_id)/100),
               ml = as.numeric(bin_id) %% 100)
    my_cts = c("CT CTX A", "CT CTX P/RSP")
    ct_prop_tidy %>%
        filter(cell_type %in% my_cts) %>%
        group_by(ap, ml, cell_type) %>%
        summarize(prop = mean(prop)) %>%
        ggplot(aes(x = ap, y = prop, col=cell_type)) +
        geom_line() +
        facet_wrap(~ml)
    ggsave("figs/coexp/ct_example_pc1.pdf")
    my_cts = c("CT CTX A/ACA", "CT CTX UL")
    ct_prop_tidy %>%
        filter(cell_type %in% my_cts) %>%
        group_by(ap, ml, cell_type) %>%
        summarize(prop = mean(prop)) %>%
        ggplot(aes(x = ml, y = prop, col=cell_type)) +
        geom_line() +
        facet_wrap(~ap)
    ggsave("figs/coexp/ct_example_pc2.pdf")
}

load_annotated_data = function(ccf=FALSE) {
    barseq = load_barseq()
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    clusters = read_csv("analysis/labels.csv")
    assertthat::are_equal(clusters$sample, colnames(barseq))
    barseq$cluster = clusters$cluster
    barseq$subclass = clusters$subclass
    keep_subclass = c("CT", "L2/3 IT", "L4/5 IT", "L5 IT", "L6 IT", "L6b", "NP", "PT")
    barseq = barseq[, barseq$subclass %in% keep_subclass]
    barseq = barseq[, !is.na(barseq$bin_id)]
    if (ccf) {
        barseq$x = barseq$ccf_x
        barseq$y = barseq$ccf_y
        barseq$z = barseq$ccf_z
        barseq = barseq[, !is.na(barseq$ccf_x)]
    } else {
        # slices are 200um away from each other
        barseq$x = barseq$aligned_x
        barseq$y = barseq$aligned_y
        barseq$z = barseq$slice * 200
    }
    ml_axis = tibble(x = barseq$bin_ywarped) %>%
        mutate(new_x = as.numeric(cut(x, quantile(x, probs = seq(0,1,l=21)), include.lowest = TRUE)))
    barseq$ml_axis = ml_axis$new_x
    ml_axis_legacy = tibble(x = barseq$bin_ywarped, slice = barseq$slice) %>%
        group_by(slice) %>%
        mutate(new_x = as.numeric(cut(x, quantile(x, probs = seq(0,1,l=21)), include.lowest = TRUE)))
    barseq$ml_axis_legacy = ml_axis_legacy$new_x
    return(barseq)
}

make_pseudobulks = function(barseq) {
    bin_id = barseq$slice*100 + barseq$ml_axis
    expr = as_tibble(as.matrix(t(counts(barseq)))) %>%
        mutate(subclass = barseq$subclass, cluster = barseq$cluster, bin_id = bin_id) %>%
        pivot_longer(c(-subclass, -cluster, -bin_id), names_to = "gene", values_to = "expr") %>%
        group_by(subclass, cluster, bin_id, gene) %>%
        summarize(expr = sum(expr), .groups="drop")
    expr_matrix = expr %>%
        mutate(pseudobulk_id = paste(bin_id,subclass,cluster,sep="|")) %>%
        pivot_wider(c(gene), names_from = "pseudobulk_id", values_from = "expr") %>%
        column_to_rownames("gene") %>%
        as.matrix()
    result = SingleCellExperiment(list(counts = expr_matrix))
    sample_metadata = strsplit(colnames(result), "|", TRUE)
    result$bin_id = as.integer(sapply(sample_metadata, "[", 1))
    result$subclass = sapply(sample_metadata, "[", 2)
    result$cluster = sapply(sample_metadata, "[", 3)
    result$ap = round(result$bin_id / 100)
    result$ml = result$bin_id %% 100
    return(result)
}

get_ap = function(pseudobulks, bin_id) {
    bin_to_ap = as.data.frame(colData(pseudobulks)) %>%
        select(bin_id, ap) %>%
        distinct() %>%
        deframe()
    return(bin_to_ap[bin_id])
}

get_ml = function(pseudobulks, bin_id) {
    bin_to_ml = as.data.frame(colData(pseudobulks)) %>%
        select(bin_id, ml) %>%
        distinct() %>%
        deframe()
    return(bin_to_ml[bin_id])
}

plot_ct_patterns = function(cluster_patterns) {
    pdf("figs/ct_patterns/cluster_patterns.pdf")
    to_plot = cluster_patterns %>%
        group_by(ml, ap) %>%
        summarize(n = sum(n), .groups="drop") %>%
        ggplot(aes(x=ml, y=ap, fill = n)) +
        geom_raster() +
        ggtitle("Total number of cells") +
        theme_classic() +
        scale_fill_gradient(low = "white", high = "darkblue")
    print(to_plot)
    to_plot = cluster_patterns %>%
        group_by(ml, ap, subclass) %>%
        summarize(n = sum(n), .groups="drop") %>%
        group_by(ml, ap) %>%
        mutate(prop = n / sum(n)) %>%
        ggplot(aes(x=ml, y=ap, fill = prop)) +
        geom_raster() +
        facet_wrap(~ subclass) +
        ggtitle("Total proportion of cells by subclass") +
        theme_classic() +
        scale_fill_gradient(low = "white", high = "darkblue")
    print(to_plot)
    for (my_subclass in sort(unique(cluster_patterns$subclass))) {
        to_plot = cluster_patterns %>%
            filter(subclass == my_subclass) %>%
            ggplot(aes(x=ml, y=ap, fill = subclass_prop)) +
            geom_raster() +
            facet_wrap(~ cluster) +
            theme_classic() +
            scale_fill_gradient(low = "white", high = "darkblue")
        print(to_plot)
    }
    dev.off()
}

plot_gene = function(barseq, g, label="cluster") {
        p1 = tibble(x=barseq$bin_ywarped,y=barseq$slice,expr=counts(barseq)[g,],label=barseq[[label]]) %>%
            mutate(x = cut(x, breaks = quantile(x, seq(0,1,0.02)), include.lowest = TRUE)) %>%
            group_by(x,y) %>%
            summarize(expr = mean(expr), n = n(), .groups="drop") %>%
            ggplot(aes(x,y,fill=sqrt(expr))) +
            geom_tile() +
            ggtitle(g) +
            theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
            labs(x="M-L axis",y="A-P axis")
        p2 = tibble(x=barseq$bin_ywarped,y=barseq$slice,expr=counts(barseq)[g,],label=barseq[[label]]) %>%
            mutate(x = cut(x, breaks = quantile(x, seq(0,1,0.02)), include.lowest = TRUE)) %>%
            group_by(x,y,label) %>%
            summarize(expr = mean(expr), n = n(), .groups="drop") %>%
            ggplot(aes(x,y,fill=sqrt(expr))) +
            geom_tile()  +
            theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
            labs(x="M-L axis",y="A-P axis") +
            facet_wrap(~label)
        p3 = tibble(ap=barseq$slice,ml=barseq$ml_axis, expr=cpm(barseq)[g,],label=barseq[[label]]) %>%
            group_by(ap, label) %>%
            summarize(expr = mean(expr), n = n(), .groups="drop") %>%
            filter(n>20) %>%
            ggplot(aes(ap,sqrt(expr),col=label)) +
            geom_line(aes(group=label)) +
            ggtitle(g)
        p4 = tibble(ap=barseq$slice,ml=barseq$ml_axis, expr=cpm(barseq)[g,],label=barseq[[label]]) %>%
            group_by(ml, label) %>%
            summarize(expr = mean(expr), n = n(), .groups="drop") %>%
            filter(n>20) %>%
            ggplot(aes(ml,sqrt(expr),col=label)) +
            geom_line(aes(group=label)) +
            ggtitle(g)
        gridExtra::grid.arrange(p1,p2,p3,p4)
}

if (sys.nframe() == 0) {
    main()
}
