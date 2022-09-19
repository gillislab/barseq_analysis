# -*- coding: utf-8 -*-

library(tidyverse)
library(SingleCellExperiment)
source("dataset.R")


main = function() {
    barseq = load_annotated_data()
    make_pseudobulks()
}

load_annotated_data = function() {
    barseq = load_barseq()
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    clusters = read_csv("analysis/labels.csv")
    assertthat::are_equal(clusters$sample, colnames(barseq))
    barseq$cluster = clusters$cluster
    barseq$subclass = clusters$subclass
    barseq = barseq[, barseq$subclass %in% ctx_subclass()]
    barseq = barseq[, !is.na(barseq$bin_id)]
    barseq = barseq[!grepl("unused", rownames(barseq)),]
    
    ml_axis = tibble(x = barseq$bin_ywarped, slice = barseq$slice) %>%
        group_by(slice) %>%
        mutate(ml_range = cut(x, quantile(x, probs = seq(0,1,l=21)), include.lowest = TRUE),
               ml_bin = as.numeric(cut(x, quantile(x, probs = seq(0,1,l=21)), include.lowest = TRUE))) %>%
        ungroup() %>%
        mutate(x_tmp = str_sub(ml_range, 2, -2)) %>% 
        separate(x_tmp, c("ml_min", "ml_max"), sep = ",") %>% 
        mutate_at(c("ml_min", "ml_max"), as.double)
    barseq$ml = ml_axis$ml_bin
    barseq$ml_min = ml_axis$ml_min
    barseq$ml_max = ml_axis$ml_max
    return(barseq)
}

make_pseudobulks = function(barseq) {
    bin_id = barseq$slice*100 + barseq$ml
    ml_pos = tibble(bin_id = bin_id, ml = barseq$ml, ml_min = barseq$ml_min, ml_max = barseq$ml_max) %>%
        distinct()
    bin_to_ml = deframe(select(ml_pos, bin_id, ml))
    bin_to_ml_min = deframe(select(ml_pos, bin_id, ml_min))
    bin_to_ml_max = deframe(select(ml_pos, bin_id, ml_max))
    tidy_expr = as_tibble(as.matrix(t(counts(barseq)))) %>%
        mutate(subclass = barseq$subclass, cluster = barseq$cluster, bin_id = bin_id) %>%
        pivot_longer(c(-subclass, -cluster, -bin_id), names_to = "gene", values_to = "expr")
    
    # subclass level
    expr_matrix = tidy_expr %>%
        group_by(subclass, bin_id, gene) %>%
        summarize(expr = sum(expr), .groups="drop") %>%
        mutate(pseudobulk_id = paste(bin_id,subclass,sep="|")) %>%
        pivot_wider(c(gene), names_from = "pseudobulk_id", values_from = "expr") %>%
        column_to_rownames("gene") %>%
        as.matrix()
    result = SingleCellExperiment(list(counts = expr_matrix))
    sample_metadata = strsplit(colnames(result), "|", TRUE)
    result$bin_id = as.integer(sapply(sample_metadata, "[", 1))
    result$subclass = sapply(sample_metadata, "[", 2)
    result$ap = round(result$bin_id / 100)
    result$ml_bin = bin_to_ml[as.character(result$bin_id)]
    result$ml_min = bin_to_ml_min[as.character(result$bin_id)]
    result$ml_max = bin_to_ml_max[as.character(result$bin_id)]
    pseudobulk_id = as.factor(paste(bin_id,barseq$subclass,sep="|"))
    n_cells = tabulate(pseudobulk_id)
    names(n_cells) = levels(pseudobulk_id)
    result$n_cells = n_cells[colnames(result)]
    saveRDS(result, "data/subclass_pseudobulks.rds")

    # cluster level
    expr_matrix = tidy_expr %>%
        group_by(subclass, cluster, bin_id, gene) %>%
        summarize(expr = sum(expr), .groups="drop") %>%
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
    result$ml_bin = bin_to_ml[as.character(result$bin_id)]
    result$ml_min = bin_to_ml_min[as.character(result$bin_id)]
    result$ml_max = bin_to_ml_max[as.character(result$bin_id)]
    pseudobulk_id = as.factor(paste(bin_id,barseq$subclass,barseq$cluster,sep="|"))
    n_cells = tabulate(pseudobulk_id)
    names(n_cells) = levels(pseudobulk_id)
    result$n_cells = n_cells[colnames(result)]
    saveRDS(result, "data/cluster_pseudobulks.rds")
}

if (sys.nframe() == 0) {
    main()
}
