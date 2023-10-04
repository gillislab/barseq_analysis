
library(tidyverse)
library(patchwork)
source("dataset.R")
source("aurocs.R")
source("visualization.R")


main = function() {
    barseq = load_dataset()

    average_ctl_correlation(barseq, "L2/3 IT", "L2/3 IT ML_2", "VISp", "_l23it")
    local_composition(barseq, "L2/3 IT", "VISp", "_l23it")

    average_ctl_correlation(barseq, "L4/5 IT", "L4/5 IT ML-P", "VISp", "_l45it")
    local_composition(barseq, "L4/5 IT", "VISp", "_l45it")
    
    if (!file.exists("results/ctl_comp_across.csv.gz")) {
        system.time({ ctl_comp = ctl_composition(barseq) }) # ~1h30
        write_csv(ctl_comp, "analysis/ctl_comp_across.csv.gz")
    }
    if (!file.exists("results/ctl_comp_within.csv.gz")) {
        system.time({ ctl_comp = ctl_composition(barseq, type = "within") }) # ~ 1h
        write_csv(ctl_comp, "analysis/ctl_comp_within.csv.gz")
    }
    ctl_comp = read_csv("analysis/ctl_comp_across.csv.gz")
    plot_ctl_comp(barseq, ctl_comp)
    
    neighborhood_shift(barseq, "L2/3 IT", "VISp", "l23it_")
    neighborhood_shift(barseq, "L4/5 IT", "VISp", "l45it_")
    neighborhood_shift(barseq, c("L5 IT"), "VISp", "l5it_")
    neighborhood_shift(barseq, c("L6 IT"), "VISp", "l6it_")
    neighborhood_shift(barseq, c("L5 IT","L6 IT"), "VISp", "l56it_")
    neighborhood_shift(barseq, c("L6b", "L6b/CT"), "VISp", "l6b_")
    neighborhood_shift(barseq, "L6 CT", "VISp", "l6ct_")
    neighborhood_shift(barseq, "L5 ET", "VISp", "l5et_")
    neighborhood_shift(barseq, "NP", "VISp", "np_")
    neighborhood_shift(barseq, "Car3", "VISp", "car3_")
    all_h2 = unique(barseq$h2)
    neighborhood_shift(barseq, all_h2, "VISp", "all_")
}

average_ctl_correlation = function(barseq, my_ct, subtype, my_ccf, suffix="", fig_dir = "fig/shift_comp") {
    ref = barseq[, barseq$condition == "CTL" & barseq$h2 == my_ct]
    query = barseq[, barseq$h2 == my_ct]

    ## average correlation to CTL
    # compute average brain region expression (in correlation space)(metaneighbor code)
    query_score =  lapply(set_names(unique(query$litter)), function(query_litter) {
        keep_cell = ref$litter != query_litter
        ref_profile = MetaNeighbor::trainModel(rownames(ref), ref[,keep_cell], "counts", ref$h3_name[keep_cell], ref$CCFparentname[keep_cell])
        n_cells = unlist(ref_profile[1,])
        ref_profile = as.matrix(ref_profile[-1,])
        ref_profile = t(t(ref_profile)/n_cells)
        ref_profile = ref_profile[, n_cells>=100]
        # preprocess query data -> correlation space
        keep_cell = query$litter == query_litter
        query_profile = normalize_cols(counts(query[, keep_cell]))
        # average correlation
        query_score = crossprod(ref_profile, query_profile)
        ct_region = strsplit(rownames(query_score), "|", fixed = TRUE)
        query_score = as_tibble(query_score) %>%
            mutate(ref_ct = map_chr(ct_region, 1), ref_ccf = map_chr(ct_region, 2)) %>%
            pivot_longer(c(-ref_ct, -ref_ccf), names_to = "cell_id", values_to = "corr")
        query_metadata = tibble(cell_id = colnames(query), query_ct = query$h3_name, query_ccf = query$CCFparentname, query_condition = query$condition, query_sample = query$sample)
        query_score = inner_join(query_score, query_metadata)
        query_score = filter(query_score, ref_ct == query_ct & ref_ccf == query_ccf) %>%
            select(sample = query_sample, ccf = ref_ccf, ct = ref_ct, corr, condition = query_condition)
        gc()
        return(query_score)
    })
    query_score = bind_rows(query_score, .id = "litter")
    
    ## plot results
    # show composition shift
    p = tibble(ccf = query$CCFparentname, ct = query$h3_name, sample=query$sample, litter = query$litter) %>%
        filter(ccf == my_ccf) %>%
        dplyr::count(ccf, ct, sample, litter) %>%
        group_by(sample) %>%
        mutate(f = n / sum(n)) %>%
        ggplot(aes(x=litter,y=f,fill=sample)) +
        facet_wrap(~ ct) +
        scale_fill_manual(values = sample_cols()) +
        geom_col(show.legend = FALSE, position = "dodge") +
        theme_classic(base_size = 16) +
        labs(x=NULL, y="% neurons") +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        scale_y_continuous(labels = scales::percent)
    ggsave(file.path(fig_dir, glue::glue("shift_comp{suffix}.pdf")), p)
    
    # show shift in average correlation across H3 types
    cor_shift = query_score %>%
        filter(ccf == my_ccf) %>%
        mutate(is_enu = condition == "ENU") %>%
        group_by(ct) %>%
        summarize(auroc = c(compute_aurocs(corr, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu))),
                  delta_cor = mean(corr[is_enu]) - mean(corr[!is_enu])) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    p1 = query_score %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=condition,y=corr)) +
        facet_wrap(~ ct) +
        geom_violin(aes(fill=condition), show.legend = FALSE) +
        geom_boxplot(width=0.25) +
        geom_text(data = cor_shift, aes(x=1.5, y=0.6, label=signif(fdr,1))) +
        theme_classic(base_size = 16) +
        ggtitle(my_ccf) +
        labs(x = NULL, y="Average correlation to CTL neurons")
    cor_shift = query_score %>%
        filter(ccf == my_ccf) %>%
        mutate(is_enu = condition == "ENU") %>%
        group_by(ct, litter) %>%
        summarize(auroc = c(compute_aurocs(corr, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu))),
                  delta_cor = mean(corr[is_enu]) - mean(corr[!is_enu])) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    p2 = query_score %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=litter,y=corr)) +
        facet_wrap(~ ct) +
        geom_violin(aes(fill=sample), show.legend = FALSE) +
        geom_boxplot(aes(group=sample), width=0.25, position = position_dodge(0.9)) +
        geom_text(data = cor_shift, aes(x=litter, y=0.2, label=signif(fdr,1))) +
        theme_classic(base_size = 16) +
        ggtitle(my_ccf) +
        scale_fill_manual(values = sample_cols()) +
        labs(x = NULL, y="Average correlation to CTL neurons")
    p = p1 | p2
    ggsave(file.path(fig_dir, glue::glue("shift_corr{suffix}.pdf")), p, width = 40/3)

    # show shift in correlation across ccf
    p1 = query_score %>%
        filter(ct == subtype) %>%
        ggplot(aes(x=condition,y=corr)) +
        facet_wrap(~ ccf) +
        geom_violin(aes(fill=condition)) +
        geom_boxplot(width=0.25) +
        theme_classic(base_size = 16) +
        ggtitle(subtype)
    cor_shift_ct = query_score %>%
        filter(ct == subtype) %>%
        mutate(is_enu = condition == "ENU") %>%
        group_by(ccf) %>%
        summarize(auroc = c(compute_aurocs(corr, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu))),
                  delta_cor = mean(corr[is_enu]) - mean(corr[!is_enu])) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    p2 = ggplot(cor_shift_ct, aes(x=2*(auroc-0.5), y=-log10(fdr), label=ccf)) +
        geom_point() +
        ggrepel::geom_text_repel(data = . %>% filter(fdr<=0.01), max.overlaps = Inf) +
        geom_hline(yintercept = -log10(0.01), linetype="dashed") +
        theme_classic(base_size = 16) +
        ggtitle(subtype) +
        labs(x="Correlation shift (normalized Wilcoxon)")
    p = p1 | p2
    ggsave(file.path(fig_dir, glue::glue("shift_across_ccf{suffix}.pdf")), p, width = 40/3)
    
    # simple visualizations -> UMAP + flatmap
    #subquery = query[, query$CCFparentname == my_ccf]
    #p2 = plot_raster_umap(subquery[, subquery$condition == "ENU"], color_by = "h3_name", text_by = "h3_name", point_size = 5, alpha = 0.01) + ggtitle("ENU") + guides(col="none")
    #p1 = plot_raster_umap(subquery[, subquery$condition == "CTL"], color_by = "h3_name", text_by = "h3_name", point_size = 5, alpha = 0.01) + ggtitle("CTL") + guides(col="none")
    #p1
    #p2
    flatmap_enu = compute_flatmap(query[, query$condition == "ENU"], "h3_name")
    flatmap_ctl = compute_flatmap(query[, query$condition == "CTL"], "h3_name")
    p1 = plot_flatmap_by_group(flatmap_ctl, n_threshold = 20) + ggtitle("CTL")
    p2 = plot_flatmap_by_group(flatmap_enu, n_threshold = 20) + ggtitle("ENU")
    flatmap_diff = inner_join(flatmap_enu, flatmap_ctl, by = c("group","supergroup","x","y"), suffix = c("_enu", "_ctl"))
    flatmap_diff$f01 = flatmap_diff$f01_enu - flatmap_diff$f01_ctl
    flatmap_diff$n_total = pmin(flatmap_diff$n_total_ctl, flatmap_diff$n_total_enu)
    p_diff = plot_flatmap_by_group(flatmap_diff, n_threshold = 20) +
        scale_fill_gradient2(low=RColorBrewer::brewer.pal("RdBu", n = 11)[11],
                             high=RColorBrewer::brewer.pal("RdBu", n = 11)[1])
    p = p1|p2|(p_diff + ggtitle("ENU-CTL"))
    ggsave(file.path(fig_dir, glue::glue("shift_flatmap{suffix}.pdf")), p, width = 60/3)
}

local_composition = function(barseq, my_h2, my_ccf, suffix="", fig_dir = "fig/shift_comp") {
    ## assess mixity of single-neuron neighborhoods
    # - metric 1: ratio of ENU vs CTL (local imbalance)
    # - metric 2: auroc of ENU vs CTL (local mixity)
    sce = barseq[, barseq$h2 == my_h2]
    is_ccf = sce$CCFparentname == my_ccf
    norm_data = normalize_cols(counts(sce))
    # find top 100 neighbors for each cell in the CT + CCF of interest-> stratify by litter
    # ~ 4 minutes for L2/3 IT
    results = lapply(set_names(unique(sce$litter)), function(litter) {
        is_litter = sce$litter == litter
        neighbors = find_neighbors(norm_data[, !is_litter], norm_data[, is_litter & is_ccf], 100)
        is_ctl = matrix(sce$condition[!is_litter][neighbors], nrow(neighbors)) == "CTL"
        ccf = matrix(sce$CCFparentname[!is_litter][neighbors], nrow(neighbors))
        result = tibble(
            cell_id = colnames(sce)[is_litter & is_ccf],
            local_visp = rowMeans(ccf == "VISp"),
            local_comp = rowMeans(is_ctl),
            local_auc = compute_aurocs(1:100, t(is_ctl)*1),
            local_pval = auroc_p_value(local_auc, t(is_ctl)*1, two_tailed = TRUE)
        )
    })
    results = bind_rows(results, .id = "litter")
    cell_idx = match(results$cell_id, colnames(sce))
    results$condition = sce$condition[cell_idx]
    results$sample = sce$sample[cell_idx]
    results$ct = sce$h3_name[cell_idx]
    results$ccf = sce$CCFparentname[cell_idx]

    sce$permuted_condition = tibble(litter = sce$litter, condition = sce$condition, ccf=sce$CCFparentname, h3 = sce$h3) %>%
        group_by(litter, ccf, h3) %>%
        mutate(condition = sample(condition)) %>%
        pull(condition)
    ctl = lapply(set_names(unique(sce$litter)), function(litter) {
        is_litter = sce$litter == litter
        neighbors = find_neighbors(norm_data[, !is_litter], norm_data[, is_litter & is_ccf], 100)
        is_ctl = matrix(sce$permuted_condition[!is_litter][neighbors], nrow(neighbors)) == "CTL"
        ccf = matrix(sce$CCFparentname[!is_litter][neighbors], nrow(neighbors))
        result = tibble(
            cell_id = colnames(sce)[is_litter & is_ccf],
            local_visp = rowMeans(ccf == "VISp"),
            local_comp = rowMeans(is_ctl),
            local_auc = compute_aurocs(1:100, t(is_ctl)*1),
            local_pval = auroc_p_value(local_auc, t(is_ctl)*1, two_tailed = TRUE)
        )
    })
    ctl = bind_rows(ctl, .id = "litter")
    cell_idx = match(ctl$cell_id, colnames(sce))
    ctl$condition = sce$condition[cell_idx]
    ctl$sample = sce$sample[cell_idx]
    ctl$ct = sce$h3_name[cell_idx]
    ctl$ccf = sce$CCFparentname[cell_idx]
    ctl %>% group_by(litter) %>% summarize(mean(local_comp))
    
    # what fraction of neighbors are in VISp
    results %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=condition, y=local_visp)) +
        geom_violin(show.legend = FALSE) +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        labs(x=NULL,y="VISp identity (fraction neighbors in VISp)")
    
    shift_pval = results %>%
        mutate(is_enu = condition == "ENU") %>%
        filter(ccf == my_ccf) %>%
        group_by(ct, litter) %>%
        summarize(auroc = c(compute_aurocs(local_comp, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu)))) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    average_comp = mean(sce$condition == "CTL")
    p1 = results %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=litter, y=local_comp)) +
        geom_violin(aes(fill=sample), show.legend = FALSE) +
        geom_boxplot(aes(group=sample), width=0.25, position = position_dodge(0.9)) +
        geom_text(data=shift_pval, aes(y=0.9,label=signif(fdr,1))) +
        geom_hline(yintercept = average_comp, linetype="dashed") +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        scale_fill_manual(values = sample_cols()) +
        labs(x=NULL,y="Local composition (fraction neighbors in CTL)")
    shift_pval = results %>%
        mutate(is_enu = condition == "ENU") %>%
        filter(ccf == my_ccf) %>%
        group_by(ct) %>%
        summarize(auroc = c(compute_aurocs(local_comp, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu)))) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    p2 = results %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=condition, y=local_comp)) +
        geom_violin(aes(fill=condition), show.legend = FALSE) +
        geom_boxplot(width=0.25) +
        geom_text(data=shift_pval, aes(x=1.5, y=0.9,label=signif(fdr,1))) +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        labs(x=NULL,y="Local composition (fraction neighbors in CTL)")
    p = p2 | p1
    ggsave(file.path(fig_dir, glue::glue("shift_local_comp{suffix}.pdf")), p, width=40/3)

    shift_pval = ctl %>%
        mutate(is_enu = condition == "ENU") %>%
        filter(ccf == my_ccf) %>%
        group_by(ct, litter) %>%
        summarize(auroc = c(compute_aurocs(local_comp, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu)))) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    p1 = ctl %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=litter, y=local_comp)) +
        geom_violin(aes(fill=sample), show.legend = FALSE) +
        geom_boxplot(aes(group=sample), width=0.25, position = position_dodge(0.9)) +
        geom_text(data=shift_pval, aes(y=0.9,label=signif(fdr,1))) +
        geom_hline(yintercept = average_comp, linetype="dashed") +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        scale_fill_manual(values = sample_cols()) +
        labs(x=NULL,y="Local composition (fraction neighbors in CTL)")
    shift_pval = ctl %>%
        mutate(is_enu = condition == "ENU") %>%
        filter(ccf == my_ccf) %>%
        group_by(ct) %>%
        summarize(auroc = c(compute_aurocs(local_comp, is_enu)),
                  pval = c(auroc_p_value(auroc, as.matrix(is_enu)))) %>%
        mutate(fdr = p.adjust(pval, "fdr"))
    p2 = ctl %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=condition, y=local_comp)) +
        geom_violin(aes(fill=condition), show.legend = FALSE) +
        geom_boxplot(width=0.25) +
        geom_text(data=shift_pval, aes(x=1.5, y=0.9,label=signif(fdr,1))) +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        labs(x=NULL,y="Local composition (fraction neighbors in CTL)")
    p = p2 | p1
    ggsave(file.path(fig_dir, glue::glue("shift_local_comp_permuted{suffix}.pdf")), p, width=40/3)

    p1=results %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=litter, y=2*local_auc-1)) +
        geom_violin(aes(fill=sample), show.legend = FALSE) +
        geom_boxplot(aes(group=sample), width=0.25, position = position_dodge(0.9)) +
        geom_hline(yintercept = 0, linetype="dashed") +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        scale_fill_manual(values = sample_cols()) +
        labs(x=NULL,y="Mixity index (normalized AUROC)")
    p1
    p2=results %>%
        mutate(local_fdr = p.adjust(local_pval, "fdr")) %>%
        filter(ccf == my_ccf) %>%
        ggplot(aes(x=litter, y=local_fdr)) +
        geom_violin(aes(fill=sample), show.legend = FALSE) +
        geom_boxplot(aes(group=sample), width=0.25, position = position_dodge(0.9)) +
        geom_hline(yintercept = 0.05, linetype="dashed") +
        facet_wrap(~ ct) +
        theme_classic(base_size = 16) +
        scale_fill_manual(values = sample_cols()) +
        labs(x=NULL,y="Mixity index (FDR)")
    p=p1|p2
    ggsave(file.path(fig_dir, glue::glue("shift_local_mixity{suffix}.pdf")), p, width = 40/3)
           
    all_ct = lapply(unique(results$ct), function(my_h3) {
        my_cor = results %>%
            filter(ct == my_h3 & condition == "CTL") %>%
            with(., cor(local_comp, local_visp))
        p1=results %>%
            filter(ct == my_h3 & condition == "CTL") %>%
            ggplot(aes(x=local_visp,y=local_comp)) +
            geom_jitter(alpha=0.1) +
            theme_classic(base_size = 16) +
            ggtitle(paste(my_h3, round(my_cor,2))) +
            labs(x = "VISp identity\n(fraction neighbors in VISp)", y = "Local composition (fraction neighbors in CTL)")
        p2=results %>%
            filter(ct == my_h3) %>%
            mutate(local_visp = cut(local_visp, breaks = quantile(local_visp, probs=seq(0,1,by=0.2)), include.lowest = TRUE)) %>%
            dplyr::count(local_visp, condition) %>%
            group_by(condition) %>%
            mutate(f = n / sum(n)) %>%
            ggplot(aes(x=local_visp, fill=condition)) +
            geom_col(aes(y = f), position="dodge") +
            theme_classic(base_size = 16) +
            theme(legend.position = "top", axis.text.x = element_text(angle=45, hjust=1)) +
            labs(x = "VISp identity\n(fraction neighbors in VISp)", y = "Density")
        p1 | p2
    })
    pdf(file.path(fig_dir, glue::glue("shift_detail{suffix}.pdf")))
    print(all_ct)
    dev.off()
}

ctl_composition = function(sce, type=c("across","within"), k=100) {
    type = match.arg(type)
    norm_data = normalize_cols(counts(sce))
    results = lapply(set_names(unique(sce$litter)), function(litter) {
        is_litter = sce$litter == litter
        is_ref = if(type == "across") {!is_litter} else {is_litter}
        neighbors = find_neighbors(norm_data[, is_ref], norm_data[, is_litter], k)
        is_ctl = matrix(sce$condition[is_ref][neighbors], nrow(neighbors)) == "CTL"
        result = tibble(
            cell_id = colnames(sce)[is_litter],
            local_comp = rowMeans(is_ctl)
        )
    })
    results = bind_rows(results)
    return(results)
}

plot_ctl_comp = function(barseq, ctl_comp) {
    metadata = as_tibble(colData(barseq), rownames = "cell_id") %>%
        select(cell_id, ccf = CCFparentname, h1, h2, h3 = h3_name, study_id, condition)
    to_plot = inner_join(ctl_comp, metadata, by="cell_id")
    
    p1_data = to_plot %>%
        group_by(h3) %>%
        summarize(p_ctl = mean(condition == "CTL"), n = n()) %>%
        filter(n >= 100)
    p1 = p1_data %>%
        ggplot(aes(x = 100*p_ctl)) +
        geom_histogram(bins = 15, fill="skyblue2", color="black") +
        theme_classic(base_size = 16) +
        labs(x="%CTL neurons",y="# H3 types")
    p1
    range(p1_data$p_ctl); mean(p1_data$p_ctl) # 0.45, 0.57; 0.5

    p2_data = to_plot %>%
        group_by(h3, ccf) %>%
        summarize(p_ctl = mean(condition == "CTL"), n = n()) %>%
        filter(n >= 100)
    p2 = p2_data %>%
        ggplot(aes(x = 100*p_ctl)) +
        geom_histogram(bins = 20, fill="skyblue2", color="black") +
        theme_classic(base_size = 16) +
        labs(x="%CTL neurons",y="# (H3 types x CCF region)")
    p2_data %>% arrange(p_ctl)
    range(p2_data$p_ctl); mean(p2_data$p_ctl) # 0.25 (VISp), 0.8 (VISpm); 0.51
    p2b = p2_data %>%
        ggplot(aes(x = p_ctl)) +
        geom_histogram(bins = 20, fill="skyblue2") +
        theme_classic(base_size = 16) +
        labs(x="%CTL neurons",y="# (H3 types x CCF region)") +
        facet_wrap(~ ccf)
    p2b
    p2c = p2_data %>%
        ggplot(aes(x = h3, y = p_ctl)) +
        geom_boxplot() +
        theme_classic(base_size = 12) +
        labs(y="%CTL neurons",x=NULL) +
        theme(axis.text.x = element_text(angle=45,hjust=1))
    p2c

    p3 = to_plot %>%
        ggplot(aes(x = 100*local_comp)) +
        geom_histogram(fill="skyblue2", color="black") +
        theme_classic(base_size = 16) +
        labs(x="%CTL neurons in neighbors",y="# neurons")
    p3
    quantile(to_plot$local_comp, probs = c(0,0.0005,0.005,0.025,0.5,0.975,0.995,0.9995,1)) 
    # 100% interval: 14% - 91%
    # 99.9% interval: 27% - 80%
    # 99% interval: 32% - 74%
    # 95%: 37% - 69%
    p3_data = to_plot %>%
        group_by(h3, ccf) %>%
        summarize(local_comp = sort(local_comp)[10])
    p3_data = to_plot %>%
        group_by(h3, ccf) %>%
        summarize(local_comp = quantile(local_comp, 0.01))
    arrange(p3_data, local_comp) %>% head()
    p3b = p3_data %>%
        ggplot(aes(x = 100*local_comp)) +
        geom_histogram(fill="skyblue2", color="black") +
        theme_classic(base_size = 16) +
        labs(x="%CTL neurons in neighbors",y="# (H3 types x CCF region)")
    p3b
    range(p3_data$local_comp) # 33% - 65%
    
    ggsave("fig/h3_comp/h3_comp_pooled.pdf", p1)
    ggsave("fig/h3_comp/h3_comp_ccf.pdf", p2)
    ggsave("fig/h3_comp/h3_comp_ccf_by_ccf.pdf", p2b)
    ggsave("fig/h3_comp/neighborhood_comp.pdf", p3)
    ggsave("fig/h3_comp/neighborhood_comp_by_h3.pdf", p3)
}

load_dataset = function(input_dir="analysis/subglu") {
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
    
    barseq$x = fct_rev(make_bins(barseq$CCF_streamlines_x, 100))
    barseq$y = fct_rev(make_bins(barseq$CCF_streamlines_y, 100))
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

plot_density_umap = function(sce, color_by, nbin=50, nlev = 10, smoothing=FALSE, n_interp=3, umap_name = "UMAP") {
    to_plot = as_tibble(reducedDim(sce, umap_name))
    colnames(to_plot) = c("UMAP1", "UMAP2")
    if (color_by %in% rownames(sce)) {
        to_plot[[color_by]] = logcounts(sce)[color_by,]
    } else {
        to_plot[[color_by]] = sce[[color_by]]
    }
    to_plot = to_plot %>%
        mutate(umap_x = cut(UMAP1, seq(min(UMAP1), max(UMAP1), length.out = nbin), include.lowest = TRUE),
               umap_y = cut(UMAP2, seq(min(UMAP2), max(UMAP2), length.out = nbin), include.lowest = TRUE)) %>%
        group_by(umap_x, umap_y) %>%
        summarize("{color_by}" := mean(.data[[color_by]]), .groups = "drop")
    if (smoothing) {
        to_plot = to_plot %>%
            reframe(smooth_density(umap_x, umap_y, .data[[color_by]], n_interp)) %>%
            rename("f" = color_by) %>%
            drop_na()
    }
    m = min(to_plot[[color_by]])
    M = max(to_plot[[color_by]])
    absM = max(-m, M)
    breaks = seq(-absM, absM, length.out = 2*nlev)
    my_cols = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 11, "RdBu"))(2*nlev-1))
    keep_break = breaks >= m & breaks <= M
    keep_break = c(FALSE, keep_break[-length(keep_break)]) | c(keep_break[-1], FALSE)
    keep_col = keep_break[-1] & keep_break[-length(keep_break)]
    breaks = breaks[keep_break]
    my_cols = my_cols[keep_col]
    p = to_plot %>%
        ggplot(aes(x=as.numeric(umap_x),y=as.numeric(umap_y),z=.data[[color_by]])) +
        stat_contour_filled(breaks = breaks) +
        scale_fill_manual(values = my_cols) +
        theme_void() +
        coord_equal()
    return(p)
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

neighborhood_shift = function(barseq, my_h2, my_ccf, prefix="", output_dir="analysis/shift_comp", fig_dir="fig/shift_comp/") {
    neighbor_file = file.path(output_dir, glue::glue("{prefix}neighbors.rds"))
    if (!file.exists(neighbor_file)) {
        sce = barseq[, barseq$h2 %in% my_h2]
        is_ccf = sce$CCFparentname == my_ccf
        is_ctl = sce$condition == "CTL"
        norm_data = normalize_cols(counts(sce))
        # find top 100 neighbors for each cell in the CT + CCF of interest-> stratify by litter
        # ~ 6 minutes for L2/3 IT and k=1000
        k=100
        system.time({
        results = lapply(set_names(unique(sce$litter)), function(litter) {
            is_litter = sce$litter == litter
            is_ref = !is_litter & is_ctl
            is_query = is_litter & is_ccf
            neighbors = find_neighbors(norm_data[, is_ref], norm_data[, is_query], k)
            cell_id_ = colnames(sce)[is_query]
            ccf = matrix(sce$CCFparentname[is_ref][neighbors], nrow(neighbors))
            result = tibble(cell_id = rep(cell_id_, k), rank = rep(1:k, each=length(cell_id_)), ccf = as.vector(ccf))
        })
        })
        results = bind_rows(results, .id = "litter")
        cell_idx = match(results$cell_id, colnames(sce))
        results$condition = sce$condition[cell_idx]
        results$h3 = sce$h3_name[cell_idx]
        saveRDS(results, neighbor_file)
    }
    results = readRDS(neighbor_file)
        
    # shift in neighborhood across CCF areas?
    result_file = file.path(output_dir, glue::glue("{prefix}h2_by_litter.csv"))
    if (!file.exists(result_file)) {
        top_k = 20
        results_h2_by_litter = results %>%
            filter(rank <= top_k) %>%
            dplyr::count(litter, ccf, condition) %>%
            pivot_wider(id_cols = everything(), names_from = "condition", values_from = "n", values_fill = 0) %>%
            group_by(litter) %>%
            mutate(a=ENU,b=sum(ENU)-ENU,c=CTL,d=sum(CTL)-CTL) %>%
            mutate(f_ENU = ENU/sum(ENU),
                   f_CTL = CTL/sum(CTL),
                   delta = f_ENU-f_CTL,
                   log_or = log(a+1) - log(b+1) - log(c+1) + log(d+1),
                   log_or_se = sqrt(1/a + 1/b + 1/c + 1/d),
                   pval_enu = phyper(a - 0.5, a+b, c+d, a+c, lower.tail = FALSE),
                   pval_ctl = phyper(a + 0.5, a+b, c+d, a+c, lower.tail = TRUE),
                   pval = pmin(2*pmin(pval_enu, pval_ctl), 1),
                   fdr = p.adjust(pval, "fdr")) %>%
            select(-a, -b, -c, -d)
        write_csv(results_h2_by_litter, result_file)
    }
    results_h2_by_litter = read_csv(result_file, show_col_types = FALSE)
    
    result_file = file.path(output_dir, glue::glue("{prefix}h2.csv"))
    if (!file.exists(result_file)) {
        results_h2 = results_h2_by_litter %>%
            group_by(ccf) %>%
            summarize(log_or = mean(log_or),
                      delta = mean(delta),
                      pval_enu = fisher_method(pval_enu),
                      pval_ctl = fisher_method(pval_ctl),
                      log_pval = pmin(2*pmin(pval_enu, pval_ctl), 1),
                      log_fdr = my_fdr(log_pval), .groups="drop")
        write_csv(results_h2, result_file)
    }
    results_h2 = read_csv(result_file, show_col_types = FALSE)

    ## make plots
    p1 = results_h2 %>%
        mutate(ccf = fct_reorder(ccf, log_fdr*sign(log_or))) %>%
        ggplot(aes(x=ccf,y=-log_fdr/log(10)*sign(log_or))) +
        geom_col(aes(fill=log_or), position = "dodge") +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
        labs(x=NULL,y="Neighborhood enrichment (-log10(FDR))") +
        scale_fill_gradient2()
    p4 = results_h2_by_litter %>%
        ungroup() %>%
        mutate(ccf = fct_reorder(ccf, log10(fdr)*sign(log_or))) %>%
        ggplot(aes(x=ccf,y=-log10(fdr)*sign(log_or))) +
        geom_col(aes(fill=log_or), position = "dodge") +
        facet_grid(litter ~ .) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "top") +
        labs(x=NULL,y="Neighborhood enrichment (-log10(FDR))") +
        scale_fill_gradient2()
    p2 = results_h2 %>%
        mutate(ccf = fct_reorder(ccf, -log_or)) %>%
        ggplot(aes(x=ccf,y=log_or)) +
        geom_col() +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "top") +
        labs(x=NULL,y="Neighborhood enrichment (log2 odds ratio)")
    p3 = results_h2 %>%
        ungroup() %>%
        mutate(ccf_neighbor = fct_reorder(ccf, -(delta))) %>%
        ggplot(aes(x=ccf_neighbor,y=delta)) +
        geom_col() +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
        labs(x=NULL,y="Neighborhood increase (fraction)")
    pdf(file.path(fig_dir, glue::glue("{prefix}neighborhood.pdf")))
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    dev.off()

    ## by H3 types?
    if (FALSE) {
        results_h3_by_litter = results %>%
            filter(rank <= top_k) %>%
            dplyr::count(litter, h3, ccf, condition) %>%
            pivot_wider(id_cols = everything(), names_from = "condition", values_from = "n", values_fill = 0) %>%
            group_by(litter, h3) %>%
            mutate(a=ENU,b=sum(ENU)-ENU,c=CTL,d=sum(CTL)-CTL) %>%
            mutate(f_ENU = ENU/sum(ENU),
                   f_CTL = CTL/sum(CTL),
                   delta = f_ENU-f_CTL,
                   log_or = log(a) - log(b) - log(c) + log(d),
                   log_or_se = sqrt(1/a + 1/b + 1/c + 1/d),
                   pval_enu = phyper(a - 0.5, a+b, c+d, a+c, lower.tail = FALSE),
                   pval_ctl = phyper(a + 0.5, a+b, c+d, a+c, lower.tail = TRUE),
                   pval = pmin(2*pmin(pval_enu, pval_ctl), 1),
                   fdr = p.adjust(pval, "fdr")) %>%
            select(-a, -b, -c, -d)
        results_h3 = results_h3_by_litter %>%
            group_by(h3, ccf) %>%
            summarize(log_or = mean(log_or),
                      delta = mean(delta),
                      pval_enu = fisher_method(pval_enu),
                      pval_ctl = fisher_method(pval_ctl),
                      log_pval = pmin(2*pmin(pval_enu, pval_ctl), 1),
                      log_fdr = my_fdr(log_pval), .groups="drop")
        results_h3_by_litter %>%
            select(litter, h3, ccf, ctl = f_CTL, enu = f_ENU) %>%
            pivot_longer(c(ctl,enu), names_to = "condition", values_to = "f") %>%
            ggplot(aes(ccf,f)) +
            geom_col(aes(fill = condition), position = "dodge") +
            facet_grid(h3 ~ .) +
            theme(axis.text.x = element_text(angle=90,vjust = 0.5, hjust = 1))
        results_h3 %>%
            mutate(ccf = fct_reorder(ccf, log_fdr*sign(log_or))) %>%
            ggplot(aes(x=ccf,y=-log_fdr/log(10)*sign(log_or))) +
            geom_col(aes(fill=log_or), position = "dodge") +
            facet_grid(h3 ~ .) +
            theme_classic(base_size = 16) +
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
            labs(x=NULL,y="Neighborhood enrichment (-log10(FDR))") +
            scale_fill_gradient2()
        results_h3 %>%
            mutate(ccf = fct_reorder(ccf, -log_or)) %>%
            ggplot(aes(x=ccf,y=log_or)) +
            geom_col() +
            facet_grid(h3 ~ .) +
            theme_classic(base_size = 16) +
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "top") +
            labs(x=NULL,y="Neighborhood enrichment (log2 odds ratio)")
        results_h3 %>%
            mutate(ccf_neighbor = fct_reorder(ccf, -(delta))) %>%
            ggplot(aes(x=ccf_neighbor,y=delta)) +
            geom_col() +
            facet_grid(h3 ~ .) +
            theme_classic(base_size = 16) +
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
            labs(x=NULL,y="Neighborhood increase (fraction)")
    }
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