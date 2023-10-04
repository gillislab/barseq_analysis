
library(tidyverse)
library(patchwork)
source("dataset.R")
source("aurocs.R")


main = function() {
    barseq = load_dataset()
    metadata = as_tibble(colData(barseq), rownames = "cell_id") %>%
        select(cell_id, h1, h2, h3=h3_name, ccf=CCFparentname, study_id, litter, condition)
    
    # compute nearest neighbors and cache results, ~1h30
    filename = "analysis/neighbors.rds"
    if (file.exists(filename)) {
        neighbors = readRDS(filename)
    } else {
        set.seed(17)
        system.time({ neighbors = compute_neighbors(barseq) })
        neighbors = neighbors[colnames(barseq),]
        saveRDS(neighbors, filename)
        gc()
    }
    
    # regular kNN classifier -> issues with composition bias, compare to randomized versions?
    set.seed(17)
    all.equal(metadata$cell_id, rownames(neighbors))
    ctl_comp = tibble(cell_id = rownames(neighbors),
                      local_comp = rowMeans(matrix((metadata$condition == "CTL")[neighbors], ncol = ncol(neighbors))))
    result = ctl_comp %>%
        mutate(predicted = ifelse(local_comp + runif(n(), -0.1,0.1), "CTL", "ENU")) %>%
        left_join(metadata) %>%
        group_by(h3,ccf,litter) %>%
        mutate(randomized_comp = sample(local_comp))
    h3_comp = result %>%
        group_by(h3, ccf) %>%
        reframe(enframe(sapply(unique(litter), function(l) { mean(condition[litter != l] == "CTL") }), "litter", "h3_comp")) %>%
        ungroup()
    result = left_join(result, h3_comp)
    
    # overall performance
    compute_aurocs(result$local_comp, result$condition == "CTL") # 0.49
    compute_aurocs(result$randomized_comp, result$condition == "CTL") # 0.49
    
    # by litter
    perf_by_litter = result %>%
        group_by(litter) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  accuracy = mean(predicted == condition),n = n())
    perf_by_litter

    # by H3
    perf_by_h3 = result %>%
        group_by(litter, h3) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  accuracy = mean(predicted == condition),n = n()) %>%
        group_by(h3) %>%
        filter(mean(n)>100) %>%
        ungroup()
    perf_by_h3 %>%
        group_by(h3) %>%
        summarize(auroc = median(auroc)) %>%
        ggplot(aes(x=auroc)) +
        geom_histogram(bins=20) +
        theme_bw(base_size = 16)
    p_h3 = perf_by_h3 %>%
        group_by(h3) %>%
        filter(sum(auroc>0.5) >= 3 & median(auroc) > 0.5) %>%
        ungroup() %>%
        mutate(h3 = fct_reorder(h3, auroc, .desc = TRUE)) %>%
        ggplot(aes(x=h3, y=auroc)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL,y="Enucleation predictability (AUROC)")
    p_h3

    # by CCF
    perf_by_ccf = result %>%
        drop_na() %>%
        group_by(litter, ccf) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  auroc_h3 = c(compute_aurocs(h3_comp, condition=="CTL")),
                  auroc_h3_v2 = c(compute_aurocs(randomized_comp, condition=="CTL")),
                  accuracy = mean(predicted == condition),n = n()) %>%
        group_by(ccf) %>%
        filter(mean(n)>100) %>%
        ungroup()
    perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc = median(auroc)) %>%
        ggplot(aes(x=auroc)) +
        geom_histogram(bins = 20) +
        theme_bw(base_size = 16)
    p_ccf = perf_by_ccf %>%
        group_by(ccf) %>%
        filter(sum(auroc>0.5) >= 3 & median(auroc) > 0.5) %>%
        ungroup() %>%
        mutate(ccf = fct_reorder(ccf, auroc, .desc = TRUE)) %>%
        ggplot(aes(x=ccf, y=auroc)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL,y="Enucleation predictability (AUROC)") +
        lims(y=c(0,1))
    p_ccf
    perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc = median(auroc), auroc_h3 = median(auroc_h3)) %>%
        ggplot(aes(x=auroc,y=auroc_h3,label=ccf)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        coord_equal() +
        lims(x=c(0.35,0.65),y=c(0.35,0.65)) +
        labs(x="Neighbor prediction (AUROC)",y="H3 composition prediction (AUROC)")
    ggsave("fig/enucleation_pred/neighbors_vs_composition.pdf")
    perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc = median(auroc), auroc_h3_v2 = median(auroc_h3_v2)) %>%
        ggplot(aes(x=auroc,y=auroc_h3_v2,label=ccf)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        coord_equal() +
        lims(x=c(0.35,0.65),y=c(0.35,0.65)) +
        labs(x="Neighbor prediction (AUROC)",y="H3 composition prediction (AUROC)")
    ggsave("fig/enucleation_pred/neighbors_vs_composition_v2.pdf")
    perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc_h3 = median(auroc_h3), auroc_h3_v2 = median(auroc_h3_v2)) %>%
        ggplot(aes(x=auroc_h3,y=auroc_h3_v2,label=ccf)) +
        geom_text() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        coord_equal() +
        lims(x=c(0.35,0.65),y=c(0.35,0.65)) +
        labs(x="Composition predictor (AUROC)",y="Randomized composition prediction (AUROC)")
    ggsave("fig/enucleation_pred/composition_vs_composition_v2.pdf")
    perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc = median(auroc), auroc_h3 = median(auroc_h3), auroc_h3_v2 = median(auroc_h3_v2)) %>%
        arrange(desc(auroc)) %>%
        head()
    perf_by_ccf %>%
        group_by(ccf) %>%
        filter(sum(auroc>0.5) >= 3 & median(auroc) > 0.5) %>%
        ungroup() %>%
        mutate(ccf = fct_reorder(ccf, auroc_h3_v2, .desc = TRUE)) %>%
        ggplot(aes(x=ccf, y=auroc_h3_v2)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL,y="Enucleation predictability (AUROC)")
    
    # by H3 x CCF
    perf = result %>%
        drop_na() %>%
        group_by(litter, h3, ccf) %>%
        mutate(null_condition = sample(condition)) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  auroc_h3 = c(compute_aurocs(h3_comp, condition=="CTL")),
                  p_auroc = auroc_p_value(auroc, as.matrix(condition=="CTL")),
                  null_auroc = c(compute_aurocs(local_comp, null_condition=="CTL")),
                  p_null_auroc = auroc_p_value(null_auroc, as.matrix(null_condition=="CTL")),
                  accuracy = mean(predicted == condition),n = n(),
                  min_ctl = min(local_comp),
                  max_ctl = max(local_comp)) %>%
        group_by(h3, ccf) %>%
        filter(mean(n)>100) %>%
        ungroup() %>%
        mutate(fdr_auroc = p.adjust(p_auroc, "fdr"))
    perf %>%
        group_by(h3, ccf) %>%
        summarize(mean_auroc = median(auroc), mean_null_auroc = median(null_auroc), auroc_h3 = median(auroc_h3)) %>%
        pivot_longer(cols = c(mean_auroc, mean_null_auroc, auroc_h3), names_to = "stat", values_to = "auroc") %>%
        ggplot(aes(x=auroc,fill=stat)) +
        geom_histogram(bins = 20, position = "dodge") # null does not work properly -> what’s a good null?
    p_h3_ccf_hist = perf %>%
        group_by(h3,ccf) %>%
        summarize(auroc = median(auroc)) %>%
        ggplot(aes(x=auroc)) +
        geom_histogram(binwidth = 0.01) +
        geom_vline(xintercept = 0.5, linetype="dashed") +
        theme_classic(base_size = 16) +
        labs(x="Enucleation predictability (AUROC)",y="#(H3 x CCF)") +
        lims(x=c(0,1))
    p_h3_ccf_hist
    p_h3_ccf = perf %>%
        mutate(ccf = paste(h3, ccf, sep = ", ")) %>%
        group_by(ccf) %>%
        summarize(auroc = median(auroc)) %>%
        ungroup() %>%
        mutate(ccf = fct_reorder(ccf, auroc, .desc = TRUE)) %>%
        ggplot(aes(x=ccf, y=auroc, label=ccf)) +
        geom_point() +
        ggrepel::geom_text_repel(data=.%>% filter(auroc>0.6), force_pull = 0.01, force=10, size=2, max.overlaps = Inf) +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_blank()) +
        labs(x="H3 x CCF (ordered by performance)",y="Enucleation predictability (AUROC)") +
        lims(y=c(0,1))
    p_h3_ccf
    p_h3_ccf_v2 = perf %>%
        mutate(ccf = paste(h3, ccf, sep = ", ")) %>%
        group_by(ccf) %>%
        filter(median(auroc) > 0.6) %>%
        ungroup() %>%
        mutate(ccf = fct_reorder(ccf, auroc, .desc = TRUE)) %>%
        ggplot(aes(x=ccf, y=auroc)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL,y="Enucleation predictability (AUROC)") +
        lims(y=c(0,1))
    p_h3_ccf_v2
    # additional stats
    perf %>%
        mutate(ccf = paste(h3, ccf, sep = ", ")) %>%
        group_by(ccf) %>%
        filter(median(auroc) > 0.6) %>%
        summarize(auroc = median(auroc), min_ctl = min(min_ctl), max_ctl = max(max_ctl)) %>%
        arrange(desc(auroc))
    
    # what happens if we average over ccf? -> non-compositional performance
    normal_perf = perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc_ccf = mean(auroc))
    noncomp_perf = perf %>%
        group_by(ccf, h3) %>%
        summarize(auroc = mean(auroc), n = sum(n)) %>%
        group_by(ccf) %>%
        summarize(auroc_within = sum(n*auroc) / sum(n))
    comp_perf = perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc_h3 = mean(auroc_h3))
    inner_join(noncomp_perf, comp_perf) %>%
        ggplot(aes(auroc_h3, auroc_within,label=ccf)) +
        geom_text() +
        coord_equal() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        lims(x=c(0.35,0.65),y=c(0.35,0.65)) +
        labs(x="H3 composition predictor (AUROC)",y="Within-H3 neighbor predictor (weighted mean AUROC)")
    ggsave("fig/enucleation_pred/within_h3_vs_composition.pdf")
    inner_join(noncomp_perf, normal_perf) %>%
        ggplot(aes(auroc_ccf, auroc_within,label=ccf)) +
        geom_text() +
        coord_equal() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        lims(x=c(0.35,0.65),y=c(0.35,0.65)) +
        labs(x="Neighbor predictor (AUROC)",y="Within-H3 neighbor predictor (weighted mean AUROC)")
    ggsave("fig/enucleation_pred/within_h3_vs_composition.pdf")
    result %>%
        filter(ccf == "VISp") %>%
        group_by(h3,litter) %>%
        summarize(h3_comp = mean(h3_comp), local_comp = mean(local_comp)) %>%
        ggplot(aes(h3_comp, local_comp)) +
        geom_point() +
        facet_wrap(~litter)

    # final plots
    p = p_ccf | p_h3_ccf_hist | p_h3_ccf_v2
    p
    ggsave("fig/enucleation_pred/summary.pdf",p,width=3*20/3)

    ###########################################################
    ## what happens if we pick neighbors within H3 types?
    ###########################################################
    filename = "analysis/neighbors_by_h3.rds"
    if (file.exists(filename)) {
        neighbors = readRDS(filename)
    } else {
        # ~50 minutes
        system.time({
            neighbors = lapply(unique(barseq$h3), function(my_h3) {
                is_my_h3 = barseq$h3 == my_h3
                idx = compute_neighbors(barseq[, is_my_h3])
                # convert to indices in the full dataset
                idx = matrix(seq_len(ncol(barseq))[is_my_h3][idx], ncol = ncol(idx), dimnames = dimnames(idx))
                return(idx)
            })
        })
        neighbors = do.call(rbind, neighbors)
        neighbors = neighbors[colnames(barseq),]
        saveRDS(neighbors, filename)
        gc()
    }

    # kNN prediction
    set.seed(17)
    all.equal(metadata$cell_id, rownames(neighbors))
    ctl_comp = tibble(cell_id = rownames(neighbors),
                      local_comp = rowMeans(matrix((metadata$condition == "CTL")[neighbors], ncol = ncol(neighbors))))
    result_h3 = ctl_comp %>%
        mutate(predicted = ifelse(local_comp + runif(n(), -0.1,0.1), "CTL", "ENU")) %>%
        left_join(metadata)

    # overall performance
    compute_aurocs(result_h3$local_comp, result_h3$condition == "CTL") # 0.495
    
    # by CCF
    perf_h3_by_ccf = result_h3 %>%
        group_by(litter, ccf) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  accuracy = mean(predicted == condition),n = n()) %>%
        group_by(ccf) %>%
        filter(mean(n)>100) %>%
        ungroup()
    perf_h3_by_ccf %>%
        inner_join(perf_by_ccf, by=c("litter", "ccf", "n"), suffix=c("_h3", "_all")) %>%
        ggplot(aes(x=auroc_all,y=auroc_h3)) +
        geom_point() +
        theme_classic(base_size = 16) # almost identical (expected)
    
    # by H3 x CCF
    perf_stratified = result_h3 %>%
        drop_na() %>%
        group_by(litter, h3, ccf) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  min_ctl = min(local_comp),
                  max_ctl = max(local_comp),
                  n=n()) %>%
        group_by(h3, ccf) %>%
        filter(mean(n)>100) %>%
        ungroup()
    
    ###########################################################
    ## what happens if we pick neighbors within H3xCCF types?
    ###########################################################
    filename = "analysis/neighbors_by_h3_ccf.rds"
    if (file.exists(filename)) {
        neighbors = readRDS(filename)
    } else {
        h3_ccf = paste(barseq$h3, barseq$CCFparentname)
        # ~35min
        system.time({
            neighbors = lapply(unique(h3_ccf), function(my_h3) {
                message(my_h3)
                is_my_h3 = h3_ccf == my_h3
                if (sum(is_my_h3) < 100) { return(NULL) } # only keep combinations containing at least 100 neurons
                idx = compute_neighbors(barseq[, is_my_h3])
                # convert to indices in the full dataset
                idx = matrix(seq_len(ncol(barseq))[is_my_h3][idx], ncol = ncol(idx), dimnames = dimnames(idx))
                return(idx)
            })
        })
        neighbors = do.call(rbind, neighbors)
        saveRDS(neighbors, filename)
        gc()
    }

    # kNN prediction
    set.seed(17)
    ctl_comp = tibble(cell_id = rownames(neighbors),
                      local_comp = rowMeans(matrix((metadata$condition == "CTL")[neighbors], ncol = ncol(neighbors))))
    ctl_comp = drop_na(ctl_comp)
    result_h3ccf = ctl_comp %>%
        mutate(predicted = ifelse(local_comp + runif(n(), -0.1,0.1), "CTL", "ENU")) %>%
        left_join(metadata)

    # overall performance
    compute_aurocs(result_h3ccf$local_comp, result_h3ccf$condition == "CTL") # 0.49
    
    # by CCF
    perf_stratified = result_h3ccf %>%
        group_by(litter, ccf) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")), n = n()) %>%
        group_by(ccf) %>%
        filter(mean(n)>100) %>%
        ungroup()
    perf_stratified %>%
        inner_join(perf_by_ccf, by=c("litter", "ccf"), suffix=c("_stratified", "_all")) %>%
        group_by(ccf) %>%
        summarize(auroc_all = median(auroc_all), auroc_h3 = median(auroc_h3), auroc_stratified = median(auroc_stratified)) %>%
        ggplot(aes(x=auroc_h3,y=auroc_stratified,label=ccf)) +
        geom_text() +
        theme_classic(base_size = 16) # correlated
    
    # by H3 x CCF
    perf_stratified = result_h3ccf %>%
        group_by(litter, h3, ccf) %>%
        summarize(auroc = c(compute_aurocs(local_comp, condition=="CTL")),
                  min_ctl = min(local_comp),
                  max_ctl = max(local_comp),
                  n=n()) %>%
        group_by(h3, ccf) %>%
        filter(mean(n)>100) %>%
        ungroup()
    perf_stratified %>%
        group_by(h3,ccf) %>%
        summarize(auroc = median(auroc)) %>%
        ggplot(aes(x=auroc)) +
        geom_histogram(binwidth = 0.01) +
        geom_vline(xintercept = 0.5, linetype="dashed") +
        theme_classic(base_size = 16) +
        labs(x="Enucleation predictability (AUROC)",y="#(H3 x CCF)") +
        lims(x=c(0,1))
    perf_stratified %>%
        group_by(h3,ccf) %>%
        summarize(auroc = median(auroc)) %>%
        arrange(desc(auroc)) %>%
        head()
    perf_stratified %>%
        mutate(ccf = paste(h3, ccf, sep = ", ")) %>%
        group_by(ccf) %>%
        filter(median(auroc) > 0.65) %>%
        ungroup() %>%
        mutate(ccf = fct_reorder(ccf, auroc, .desc = TRUE)) %>%
        ggplot(aes(x=ccf, y=auroc)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        geom_hline(yintercept = 0.5, linetype = "dashed") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x=NULL,y="Enucleation predictability (AUROC)") +
        lims(y=c(0,1))
    
    ## compare different predictors
    normal_perf = perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc_ccf = mean(auroc))
    noncomp_perf = perf_stratified %>%
        group_by(ccf, h3) %>%
        summarize(auroc = mean(auroc, na.rm = TRUE), n = sum(n)) %>%
        group_by(ccf) %>%
        summarize(auroc_within = sum(n*auroc) / sum(n))
    comp_perf = perf_by_ccf %>%
        group_by(ccf) %>%
        summarize(auroc_h3 = mean(auroc_h3))
    inner_join(noncomp_perf, comp_perf) %>%
        ggplot(aes(auroc_h3, auroc_within,label=ccf)) +
        geom_text() +
        coord_equal() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        lims(x=c(0.3,0.7),y=c(0.3,0.7)) +
        labs(x="H3 composition predictor (AUROC)",y="Within-H3 neighbor predictor (weighted mean AUROC)")
    ggsave("fig/enucleation_pred/composition_vs_stratified.pdf")
    inner_join(noncomp_perf, normal_perf) %>%
        ggplot(aes(auroc_ccf, auroc_within,label=ccf)) +
        geom_text() +
        coord_equal() +
        geom_abline(slope=1, linetype="dashed") +
        theme_classic(base_size = 16) +
        lims(x=c(0.3,0.7),y=c(0.3,0.7)) +
        labs(x="Neighbor predictor (AUROC)",y="Neighbor predictor, stratified by H3 (weighted mean AUROC)")
    ggsave("fig/enucleation_pred/neighbors_vs_stratified.pdf")
    
    # why are the non-compositional and the neighbor predictor so correlated?
    visp_predictors = inner_join(result, result_h3ccf, by=c("cell_id","litter","ccf","h3"), suffix=c("_neighbor","_stratified")) %>%
        filter(ccf == "VISp")
    visp_predictors %>%
        group_by(h3,litter) %>%
        summarize(local_comp_neighbor = mean(local_comp_neighbor), local_comp_stratified = mean(local_comp_stratified)) %>%
        ggplot(aes(local_comp_neighbor, local_comp_stratified)) +
        geom_point() +
        facet_wrap(~litter)
    visp_predictors %>%
        group_by(h3,litter) %>%
        summarize(local_comp_neighbor = mean(local_comp_neighbor), h3_comp = mean(h3_comp)) %>%
        ggplot(aes(local_comp_neighbor, h3_comp)) +
        geom_point() +
        facet_wrap(~litter)
    visp_predictors %>%
        group_by(h3,litter) %>%
        summarize(local_comp_stratified = mean(local_comp_stratified), h3_comp = mean(h3_comp)) %>%
        ggplot(aes(local_comp_stratified, h3_comp)) +
        geom_point() +
        facet_wrap(~litter)
    visp_predictors %>%
        filter(h3 == "L2/3 IT DL") %>%
        ggplot(aes(local_comp_neighbor, local_comp_stratified)) +
        geom_point(alpha=0.15) +
        facet_wrap(~litter)
    visp_predictors %>%
        filter(h3 == "L2/3 IT DL") %>%
        ggplot(aes(local_comp_neighbor, h3_comp)) +
        geom_point(alpha=0.15) +
        facet_wrap(~litter)

    # what happens if we randomize the neighbor predictor within H3xCCF combinations?


    #########################################
    ### compare with differential abundance
    #########################################
    dc = read_csv("analysis/differential_composition.csv", show_col_types = FALSE)
    full_perf = perf %>%
        group_by(h3, ccf) %>%
        summarize(auroc = median(auroc)) %>%
        left_join(dc, by=c("h3"="cell_type","ccf"="ccf1"))
    pcc = full_perf %>%
        drop_na() %>%
        with(., cor(.$Estimate,.$auroc))
    full_perf %>%
        mutate(my_label = paste(h3, ccf)) %>%
        ggplot(aes(x=auroc,y=Estimate)) +
        geom_point() +
        ggrepel::geom_text_repel(data = . %>% filter(auroc>0.6), aes(label=my_label), size=3) +
        labs(x="Enucleation predictability (AUROC)",y="Enrichment in enucleated samples (log2FC)") +
        annotate("text", x=0.4,y=1.5,label=round(pcc,2)) +
        theme_classic(base_size = 16)
    ggsave("fig/enucleation_pred/predictability_vs_composition.pdf")
}

load_dataset = function(input_dir="analysis/subglu") {
    barseq = load_glu()
    sample_cell = strsplit(colnames(barseq), split = ".", fixed = TRUE) # replace by regex for efficiency
    barseq$sample = sapply(sample_cell, "[", 1)    
    barseq$condition = ifelse(grepl("1L", barseq$sample) | grepl("3L", barseq$sample), "CTL", "ENU")
    barseq$litter = substr(barseq$sample, 1, 4)
    barseq = barseq[, barseq$CCFparentname %in% ctx_roi()]
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

compute_neighbors = function(sce, type=c("across","within"), k=100) {
    type = match.arg(type)
    norm_data = normalize_cols(counts(sce))
    results = lapply(set_names(unique(sce$litter)), function(litter) {
        is_litter = sce$litter == litter
        is_ref = if(type == "across") {!is_litter} else {is_litter}
        neighbors = find_neighbors(norm_data[, is_ref], norm_data[, is_litter], k) # returns index of neighbors, in order
        if (ncol(neighbors) < k) {
            tmp = matrix(NA, nrow = nrow(neighbors), ncol=k, dimnames = dimnames(neighbors))
            tmp[,1:ncol(neighbors)] = neighbors
            neighbors = tmp
        }
        cell_id = matrix(seq_len(ncol(sce))[is_ref][neighbors], ncol=k) # convert to indices based on full data (not just ref)
        rownames(cell_id) = colnames(sce)[is_litter]
        return(cell_id)
    })
    results = do.call(rbind, results)
    results = results[colnames(sce),]
    return(results)
}

# Scale matrix such that all colums sum to 0 and have l2-norm of 1
normalize_cols <- function(M, ranked = TRUE) {
  result <- as.matrix(M)
  if (ranked) {
    result <- matrixStats::colRanks(result, ties.method = "average",
                                    preserveShape = TRUE)
  }
  result <- scale_cols(result)
  dimnames(result) <- dimnames(M)
  return(result)
}

# This is not exactly equivalent to scale (by a factor of sqrt(nrow(M-1)))
# and is a little faster.
# The point is that the dot product of vectors scaled this way is the
# correlation coefficient (which is not true using conventional scaling)
scale_cols <- function(M) {
    cm <- colMeans(M)
    cnorm <- 1 / sqrt(colSums(M**2) - nrow(M) * cm**2)
    matrixStats::t_tx_OP_y(matrixStats::t_tx_OP_y(M, cm, "-"), cnorm, "*")
}

find_neighbors = function(ref, target, k) {
    # BNPARAM -> KmknnParam() [DEFAULT], VptreeParam(), AnnoyParam() [APPROX.]
    # approximate version is slower here (probably because of low dim.)
    result=BiocNeighbors::queryKNN(t(ref),t(target),k, BNPARAM=BiocNeighbors::AnnoyParam())
    return(result$index)
}

fisher_method = function(pval) {
    pchisq(-2*sum(log(pval)), df = 2*length(pval), lower.tail = FALSE, log.p = TRUE)
}

my_fdr = function(p_values, log.p = TRUE) {
    i = length(p_values):1L
    o = order(p_values, decreasing = TRUE)
    n = length(p_values)
    result = rep(0, n)
    if (log.p) {
        result[o] = pmin(0, cummin(log(n) - log(i) + p_values[o]))
    } else {
        result[o] = pmin(1, cummin(n/i * p_values[o]))
    }
    return(result)
}

if (sys.nframe() == 0) {
    main()
}