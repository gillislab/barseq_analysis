# -*- coding: utf-8 -*-

library(tidyverse)
library(SingleCellExperiment)
library(ica)
library(NMF)
source("dataset.R")


main = function() {
    pseudobulks = readRDS("data/subclass_pseudobulks.rds")
    pseudobulks = pseudobulks[!rownames(pseudobulks) %in% c("Slc17a7", "Gad1"),] # consider removing Slc30a3
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
    bin_to_pos = as.data.frame(colData(pseudobulks)) %>%
        select(bin_id, ap, ml_min, ml_max) %>%
        distinct()
    
    run_nmf_scaled(pseudobulk_matrix, bin_to_pos)
    run_nmf_scaled_l1(pseudobulk_matrix, bin_to_pos)
    run_nmf(pseudobulk_matrix, bin_to_pos)
    run_pca(pseudobulk_matrix, bin_to_pos)
    run_ica(pseudobulk_matrix, bin_to_pos)
}

run_nmf_scaled = function(pseudobulk_matrix, bin_to_pos) {
    # pseudobulk ~ W.H, where W = spatial patterns (#bin x #factor), H = loadings (#factor x (#subclass.#gene))
    # for scaled NMF, should we choose l1, l2-normalization or variance normalization?
    #scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=colSds(pseudobulk_matrix)) # equal variance
    #scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=colSums(pseudobulk_matrix)) # l1-normalization
    scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=sqrt(colSums(pseudobulk_matrix**2))) # l2-normalization
    my_nmf = nmf(scaled_pseudobulk, 10) # consider setting seed parameter for reproducibility
    saveRDS(my_nmf, "analysis/spatial_patterns/spatial_nmf_scaled.rds")
    
    w = basis(my_nmf)
    colnames(w) = paste0("NMF", 1:ncol(w))
    nmf_tidy = as_tibble(w, rownames = "bin_id") %>%
        mutate(ap = get_ap(bin_to_pos, bin_id),
               ml_min = get_ml_min(bin_to_pos, bin_id),
               ml_max = get_ml_max(bin_to_pos, bin_id))
    write_csv(nmf_tidy, "analysis/spatial_patterns/spatial_nmf_scaled.csv")
    
    h = coef(my_nmf)
    rownames(h) = paste0("NMF", 1:ncol(w))
    nmf_loading = as_tibble(t(h), rownames = "var") %>%
        pivot_longer(-var, names_to = "nmf", values_to = "loading") %>%
        separate(var, sep = "\\|", into = c("gene","subclass"))
    write_csv(nmf_loading, "analysis/spatial_patterns/spatial_nmf_scaled_loadings.csv")
}

run_nmf_scaled_l1 = function(pseudobulk_matrix, bin_to_pos) {
    # pseudobulk ~ W.H, where W = spatial patterns (#bin x #factor), H = loadings (#factor x (#subclass.#gene))
    # for scaled NMF, should we choose l1, l2-normalization or variance normalization?
    #scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=colSds(pseudobulk_matrix)) # equal variance
    scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=colSums(pseudobulk_matrix)) # l1-normalization
    #scaled_pseudobulk = scale(pseudobulk_matrix, center=FALSE, scale=sqrt(colSums(pseudobulk_matrix**2))) # l2-normalization
    my_nmf = nmf(scaled_pseudobulk, 10)
    saveRDS(my_nmf, "analysis/spatial_patterns/spatial_nmf_scaled_l1.rds")
    
    w = basis(my_nmf)
    colnames(w) = paste0("NMF", 1:ncol(w))
    nmf_tidy = as_tibble(w, rownames = "bin_id") %>%
        mutate(ap = get_ap(bin_to_pos, bin_id),
               ml_min = get_ml_min(bin_to_pos, bin_id),
               ml_max = get_ml_max(bin_to_pos, bin_id))
    write_csv(nmf_tidy, "analysis/spatial_patterns/spatial_nmf_scaled_l1.csv")
    
    h = coef(my_nmf)
    rownames(h) = paste0("NMF", 1:ncol(w))
    nmf_loading = as_tibble(t(h), rownames = "var") %>%
        pivot_longer(-var, names_to = "nmf", values_to = "loading") %>%
        separate(var, sep = "\\|", into = c("gene","subclass"))
    write_csv(nmf_loading, "analysis/spatial_patterns/spatial_nmf_scaled_l1_loadings.csv")
}

get_ap = function(bin_to_pos, bin_id) {
    bin_to_ap = deframe(select(bin_to_pos, bin_id, ap))
    return(bin_to_ap[bin_id])
}

get_ml_min = function(bin_to_pos, bin_id) {
    bin_to_ml = deframe(select(bin_to_pos, bin_id, ml_min))
    return(bin_to_ml[bin_id])
}

get_ml_max = function(bin_to_pos, bin_id) {
    bin_to_ml = deframe(select(bin_to_pos, bin_id, ml_max))
    return(bin_to_ml[bin_id])
}

run_nmf = function(pseudobulk_matrix, bin_to_pos) {
    my_nmf = nmf(pseudobulk_matrix, 10)
    saveRDS(my_nmf, "analysis/spatial_patterns/spatial_nmf.rds")
    
    w = basis(my_nmf)
    colnames(w) = paste0("NMF", 1:ncol(w))
    nmf_tidy = as_tibble(w, rownames = "bin_id") %>%
        mutate(ap = get_ap(bin_to_pos, bin_id),
               ml_min = get_ml_min(bin_to_pos, bin_id),
               ml_max = get_ml_max(bin_to_pos, bin_id))
    write_csv(nmf_tidy, "analysis/spatial_patterns/spatial_nmf.csv")
    
    h = coef(my_nmf)
    rownames(h) = paste0("NMF", 1:ncol(w))
    nmf_loading = as_tibble(t(h), rownames = "var") %>%
        pivot_longer(-var, names_to = "nmf", values_to = "loading") %>%
        separate(var, sep = "\\|", into = c("gene","subclass"))
    write_csv(nmf_loading, "analysis/spatial_patterns/spatial_nmf_loadings.csv")
}

run_pca = function(pseudobulk_matrix, bin_to_pos) {
    # PCA on ANOVA residuals -> PCA on scaled pseudobulk matrix
    scaled_res = scale(pseudobulk_matrix)
    my_pca = prcomp(scaled_res)
    saveRDS(my_pca, "analysis/spatial_patterns/spatial_pcs.rds")
    my_pcs = scale(my_pca$x, center=FALSE, scale=my_pca$sdev)
    pc_tidy = as_tibble(my_pcs, rownames = "bin_id") %>%
        mutate(ap = get_ap(bin_to_pos, bin_id),
               ml_min = get_ml_min(bin_to_pos, bin_id),
               ml_max = get_ml_max(bin_to_pos, bin_id))
    write_csv(pc_tidy, "analysis/spatial_patterns/spatial_pcs.csv")
}

run_ica = function(pseudobulk_matrix, bin_to_pos) {
    # X = tcrossprod(S,M) + E, where S are the source signals and M is the mixing matrix
    # M = P%*%R, where R is an orthogonal rotation
    # Y = tcrossprod(X,Q), where Y is the whitened version of X, Q is the pseudoinverse of P (whitening matrix)
    scaled_res = scale(pseudobulk_matrix)
    my_ica = icafast(scaled_res, nc=100)
    saveRDS(my_ica, "analysis/spatial_patterns/spatial_ics.rds")
    
    ic_names = paste0("IC", 1:ncol(my_ica$S))
    ic_tidy = as_tibble(my_ica$S, rownames = "bin_id", .name_repair = ~ic_names) %>%
        mutate(ap = get_ap(bin_to_pos, bin_id),
               ml_min = get_ml_min(bin_to_pos, bin_id),
               ml_max = get_ml_max(bin_to_pos, bin_id))
    write_csv(ic_tidy, "analysis/spatial_patterns/spatial_ics.csv")

    loadings = my_ica$M
    dimnames(loadings) = list(colnames(scaled_res), ic_names)
    ic_loading = as_tibble(loadings[,1:9], rownames = "var") %>%
        pivot_longer(-var, names_to = "ic", values_to = "loading") %>%
        separate(var, sep = "\\|", into = c("gene","subclass"))
    write_csv(ic_loading, "analysis/spatial_patterns/spatial_ic_loadings.csv")
}

if (sys.nframe() == 0) {
    main()
}
