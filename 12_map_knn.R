
library(tidyverse)
library(patchwork)
source("dataset.R")


main = function() {
    # load barseq
    barseq = load_barseq()
    barseq$study_name = paste("barseq", barseq$slice, sep="_")
    clusters = read_csv("analysis/labels.csv")
    assertthat::are_equal(clusters$sample, colnames(barseq))
    barseq$cluster = clusters$cluster
    barseq$subclass = clusters$subclass
    barseq$class = clusters$class
    barseq = barseq[, barseq$class == "Glutamatergic"]
    table(barseq$CCFparentname)
        
    # load train dataset
    label_name = "cluster_label"
    dataset = load_full_yao()
    dataset = dataset[, dataset$class_label == "Glutamatergic"]
    dataset$label = dataset[[label_name]]
    table(dataset$joint_region_label)

    #predict_cell_type(barseq, dataset)
    plot_predictions(barseq, dataset)
}

predict_cell_type = function(barseq, dataset) {        
    # common set + normalize data + rescale genes (optional)
    blacklist = c("Slc17a7", "Gad1")
    gene_set = intersect(rownames(barseq), rownames(dataset))
    gene_set = gene_set[!gene_set %in% blacklist]
    f_gene = rowMeans(counts(barseq)[gene_set,]) / rowMeans(counts(dataset)[gene_set,])
    barseq_data = normalize_cols(counts(barseq)[gene_set,])
    train_data = normalize_cols(counts(dataset)[gene_set,]*f_gene)
    #barseq_data = normalize_cols(log1p(convert_to_cpm(counts(barseq)[gene_set,], 1000)), ranked = FALSE)
    #train_data = normalize_cols(log1p(convert_to_cpm(counts(dataset)[gene_set,], 1000)), ranked = FALSE)
    #barseq_data = log1p(convert_to_cpm(counts(barseq)[gene_set,], 1000))
    #train_data = log1p(convert_to_cpm(counts(dataset)[gene_set,], 1000))
    
    # make predictions
    neighbors = find_neighbors(train_data, barseq_data, 15) # MNN instead?, ~5 min
    pred = predict_cell_type(neighbors, dataset$cluster_label, colnames(barseq), barseq$cluster)
    #hist(pred$confidence)
    write_csv(pred, "analysis/knn_mapping/knn_predicted_clusters.csv")
}

plot_predictions = function(barseq, dataset) {
    cluster_to_subclass = setNames(barseq$subclass, barseq$cluster)
    cluster_to_subclass = cluster_to_subclass[unique(names(cluster_to_subclass))]
    yao_to_subclass = setNames(dataset$subclass_label, dataset$cluster_label)
    yao_to_subclass = yao_to_subclass[unique(names(yao_to_subclass))]

    pred = read_csv("analysis/knn_mapping/knn_predicted_clusters.csv")
    pred$n_reads = colSums(counts(barseq))
    pred$n_genes = colSums(counts(barseq)>0)
    pred_summary = pred %>%
#        mutate(predicted_ct = sample(predicted_ct)) %>% # used to compute null jaccard indices
        group_by(cell_type, predicted_ct) %>%
        tally() %>%
        ungroup()
    ctx_pred_matrix = pred_summary %>%
        filter(cluster_to_subclass[cell_type] %in% ctx_subclass()) %>%
        pivot_wider(c(cell_type), names_from = predicted_ct, values_from = n, values_fill = 0) %>%
        column_to_rownames("cell_type") %>%
        as.matrix()
    nonctx_pred_matrix = pred_summary %>%
        filter(!(cluster_to_subclass[cell_type] %in% ctx_subclass())) %>%
        pivot_wider(c(cell_type), names_from = predicted_ct, values_from = n, values_fill = 0) %>%
        column_to_rownames("cell_type") %>%
        as.matrix()
        
    # ROI metadata
    ref_metadata = tibble(cluster = dataset$cluster_label, subclass = dataset$subclass_label, roi = dataset$joint_region_label) %>%
        group_by(cluster, subclass, roi) %>%
        tally() %>%
        ungroup() %>%
        compute_roi_enrichment() %>%
        mutate(fdr = p.adjust(exp(log_pval), "fdr")) %>%
        group_by(subclass) %>%
        mutate(f_cells = n / sum(n)) %>%
        ungroup()
    barseq_metadata = tibble(cluster = barseq$cluster, subclass = barseq$subclass, roi = barseq$CCFparentname) %>%
        drop_na() %>%
        mutate(roi = simplify_roi(roi)) %>%
        group_by(cluster, subclass, roi) %>%
        tally() %>%
        ungroup() %>%
        compute_roi_enrichment() %>%
        mutate(fdr = p.adjust(exp(log_pval), "fdr")) %>%
        group_by(subclass) %>%
        mutate(f_cells = n / sum(n)) %>%
        ungroup()
    
    # overlap with reference (what fraction of target cells predicted as each ref type)?
    # summary_matrix = ctx_pred_matrix / rowSums(ctx_pred_matrix)
    # Jaccard coefficient
    summary_matrix = ctx_pred_matrix / (outer(rowSums(ctx_pred_matrix), colSums(ctx_pred_matrix), "+") - ctx_pred_matrix)
    summary_matrix = summary_matrix[sort(rownames(summary_matrix)),]
    nonctx_summary_matrix = nonctx_pred_matrix / (outer(rowSums(nonctx_pred_matrix), colSums(nonctx_pred_matrix), "+") - nonctx_pred_matrix)
    nonctx_summary_matrix = nonctx_summary_matrix[sort(rownames(nonctx_summary_matrix)),]
    # hist(summary_matrix, breaks = 50) # null jaccard never reaches 0.04
    # quantile(summary_matrix, probs = c(0.5, 0.9, 0.95, 0.99, 0.999, 1)) # null = 0.001, 0.005, 0.008, 0.015, 0.025, 0.04

    # summary plot for ctx types (only consider reference type with a hit > threshold)
    threshold = 0.05
    to_plot = summary_matrix
    to_plot = to_plot[, colSums(to_plot>threshold)>0, drop=FALSE]
    to_plot = to_plot[cell_type_order(),]
    col_order = order_rows_according_to_cols(t(to_plot), 100)
    p1 = as_tibble(to_plot, rownames = "cluster") %>%
        pivot_longer(-cluster, names_to = "predicted_ct", values_to = "f") %>%
        mutate(cluster = factor(cluster, cell_type_order())) %>%
        mutate(cluster = fct_rev(cluster)) %>%
        mutate(predicted_ct = factor(predicted_ct, colnames(to_plot)[col_order])) %>%
        ggplot(aes(predicted_ct, cluster, size=ifelse(f==0,NA,f), col=f)) +
        geom_point() +
        scale_color_gradient(low = "gray90", high = "blue4") +
        scale_size_area(max_size=3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size=6)) +
        theme(axis.text.y = element_text(size=6), legend.position = "bottom") +
        labs(x=NULL,y=NULL)
    ggsave("figs/knn_mapping/yao_mapping_summary.pdf", p1)
    
    # summary plot for non-ctx types
    to_plot = nonctx_summary_matrix
    to_plot = to_plot[, colSums(to_plot>threshold)>0, drop=FALSE]
    col_order = order_rows_according_to_cols(t(to_plot), 100)
    p1 = as_tibble(to_plot, rownames = "cluster") %>%
        pivot_longer(-cluster, names_to = "predicted_ct", values_to = "f") %>%
        mutate(cluster = as.factor(cluster)) %>%
        mutate(cluster = fct_rev(cluster)) %>%
        mutate(predicted_ct = factor(predicted_ct, colnames(to_plot)[col_order])) %>%
        ggplot(aes(predicted_ct, cluster, size=ifelse(f==0,NA,f), col=f)) +
        geom_point() +
        scale_color_gradient(low = "gray90", high = "blue4") +
        scale_size_area(max_size=3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size=6)) +
        theme(axis.text.y = element_text(size=6), legend.position = "bottom") +
        labs(x=NULL,y=NULL)
    ggsave("figs/knn_mapping/yao_mapping_summary_nonctx.pdf", p1)
    
    # detailed plots per CTX subclass
    max_jaccard = max(summary_matrix)
    max_cells = max(max(ref_metadata$f_cells), max(barseq_metadata$f_cells))
    my_txt_size=11
    all_plots = lapply(ctx_subclass(), function(my_subclass) {
        to_plot = summary_matrix[cluster_to_subclass[rownames(summary_matrix)] == my_subclass,,drop=FALSE]
        to_plot = to_plot[, colSums(to_plot>threshold)>0, drop=FALSE]
        #row_order = na.omit(rownames(to_plot)[match(cell_type_order(), rownames(to_plot))])
        #to_plot = to_plot[row_order,]
        #col_order = order_rows_according_to_cols(t(to_plot), 100)
        p1 = as_tibble(to_plot, rownames = "cluster") %>%
            pivot_longer(-cluster, names_to = "predicted_ct", values_to = "f") %>%
#            mutate(cluster = factor(cluster, rownames(to_plot))) %>%
            mutate(cluster = factor(cluster, cell_type_order())) %>%
            mutate(cluster = fct_rev(cluster)) %>%
#            mutate(predicted_ct = factor(predicted_ct, colnames(to_plot)[col_order])) %>%
            mutate(predicted_ct = as.factor(predicted_ct)) %>%
            ggplot(aes(predicted_ct, cluster)) +
            geom_point(aes(fill=f, col=f>0.05, size=ifelse(f==0,NA,f)), pch=21, show.legend=FALSE) +
            scale_fill_gradient(low = "gray90", high = "blue4", limits=c(0, max_jaccard)) +
            scale_color_manual(values=c(gray(1,0), "black")) +
            scale_size_area(max_size=5, limits=c(0,max_jaccard)) +
            theme_bw(base_size = my_txt_size) +
            theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
            labs(x=NULL,y=NULL) +
            ggtitle(my_subclass)
        p2 = barseq_metadata %>%
            filter(subclass == my_subclass) %>%
            filter(roi %in% roi_to_plot()) %>%
            mutate(roi = factor(roi, roi_to_plot())) %>%
            mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
            mutate(n = pmin(n,max_cells)) %>%
            complete(cluster, roi) %>%
            mutate(cluster = factor(cluster, cell_type_order())) %>%
            mutate(cluster = fct_rev(cluster)) %>%
            ggplot(aes(x=roi, y=cluster)) +
            geom_point(aes(size=f_cells, fill = log_OR, col=fdr<0.05), pch=21, show.legend=FALSE) +
            theme_bw(base_size = my_txt_size) +
            theme(axis.text.x = element_blank()) +
            scale_size_area(max_size = 5) +
            scale_fill_gradient2(low = "darkblue", mid = "gray90", high = "darkorange2") +
            scale_color_manual(values=c(gray(1,0),"black")) +
            labs(x = NULL, y = NULL)
        p3 = ref_metadata %>%
            filter(subclass %in% subclass_mapping()[[my_subclass]]) %>%
            mutate(cluster = as.factor(cluster)) %>%
            mutate(cluster = fct_rev(cluster)) %>%
            filter(roi %in% roi_to_plot()) %>%
            mutate(roi = factor(roi, roi_to_plot())) %>%
            mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
            mutate(n = pmin(n,max_cells)) %>%
            complete(cluster, roi) %>%
            ggplot(aes(x=roi, y=cluster)) +
            geom_point(aes(size=f_cells, fill = log_OR, col=fdr<0.05), pch=21, show.legend=FALSE) +
            theme_bw(base_size = my_txt_size) +
            theme(axis.text.x = element_text(angle = 45, hjust=1)) +
            scale_size_area(max_size = 5) +
            scale_fill_gradient2(low = "darkblue", mid = "gray90", high = "darkorange2") +
            scale_color_manual(values=c(gray(1,0),"black")) +
            labs(x = NULL, y = NULL)
        p = p1 | (p2/p3)
        return(p)
    })
    p = wrap_plots(all_plots, ncol=4)
    ggsave("figs/knn_mapping/yao_mapping.pdf", p, height = 3*20/3, width=4*20/3)
    
    # focus on ET types
    my_subclass = "PT"
    to_plot = summary_matrix[cluster_to_subclass[rownames(summary_matrix)] == my_subclass,,drop=FALSE] %>%
        as_tibble(rownames = "barseq_type") %>%
        pivot_longer(-barseq_type, names_to = "ref_type", values_to = "overlap") %>%
        filter(overlap > 0.1) %>%
        group_by(ref_type) %>%
        slice_max(overlap, n=1) %>%
        group_by(barseq_type) %>%
#        slice_max(overlap, n=1) %>%
        summarize(type = c(barseq_type[1], ref_type),
                  dataset = c("BARseq", rep("yao", n())),
                  all_types = paste(type, collapse="\n"), .groups="drop") %>%
        inner_join(rbind(barseq_metadata, ref_metadata), by=c("type"="cluster"))
    my_roi = roi_to_plot()
    my_roi = my_roi[!my_roi %in% c("ENT","PIR")]
    p1 = to_plot %>%
        mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
        filter(roi %in% my_roi) %>%
        mutate(roi = factor(roi, my_roi)) %>%
        group_by(all_types) %>%
        complete(roi, type, fill = list(n=0,fdr=1,log_OR=0)) %>%
        ungroup() %>%
        mutate(type = factor(type, c(cell_type_order(),sort(colnames(summary_matrix))))) %>%
        ggplot(aes(x=roi,y=log_OR)) +
        geom_col(aes(group=type, fill = dataset), position=position_dodge(0.8), width=0.8, show.legend = FALSE) +
        geom_text(aes(label=ifelse(fdr<0.05, "*", ""), group=type), position = position_dodge(0.8)) +
        facet_wrap(~ all_types, ncol=3) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
        scale_y_continuous(breaks = c(-3,0,3), limits=c(-3.5,3.5)) +
        scale_fill_manual(values=c("darkorange2","gray80")) +
        labs(x = "Brain area", y = "Regional enrichment (log odds ratio)")
    p1
    ggsave("figs/knn_mapping/yao_mapping_et_types.pdf")
    # same plot as a dot plot (easier to read)
    ct_order = to_plot %>%
        select(barseq_type, type) %>%
        distinct() %>%
        with(., split(.$type, .$barseq_type))
    ct_order = unlist(ct_order[cell_type_order()])
    p1 = to_plot %>%
#        mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
        filter(roi %in% my_roi) %>%
        mutate(roi = factor(roi, my_roi)) %>%
        group_by(all_types) %>%
        complete(roi, type, fill = list(n=0,fdr=1,log_OR=0)) %>%
        ungroup() %>%
        mutate(type = factor(type, rev(ct_order))) %>%
        ggplot(aes(x=roi,y=type)) +
        geom_point(aes(size=pmin(-log10(fdr),3), fill = log_OR, col=fdr<0.05), pch=21) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_radius(range = c(0,6)) +
#        scale_fill_manual(values=c("darkorange2","gray80")) +
        scale_fill_gradient2(low = "darkblue", mid = "gray90", high = "darkorange2") +
        scale_color_manual(values=c(gray(1,0),"black")) +
        labs(x = "Brain area", y = NULL)
    p1
    ggsave("figs/knn_mapping/yao_mapping_et_types_v2.pdf")
    p1 = to_plot %>%
#        mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
        filter(roi %in% my_roi) %>%
        mutate(roi = factor(roi, my_roi)) %>%
        group_by(all_types) %>%
        complete(roi, type, fill = list(n=0,fdr=1,log_OR=0)) %>%
        ungroup() %>%
        mutate(type = factor(type, rev(ct_order))) %>%
        ggplot(aes(x=roi,y=type)) +
        geom_point(aes(size=jaccard, fill = log_OR, col=fdr<0.05), pch=21) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_size(trans="log10", range = c(0,6)) +
#        scale_fill_manual(values=c("darkorange2","gray80")) +
        scale_fill_gradient2(low = "darkblue", mid = "gray90", high = "darkorange2") +
        scale_color_manual(values=c(gray(1,0),"black")) +
        labs(x = "Brain area", y = NULL)
    p1
    ggsave("figs/knn_mapping/yao_mapping_et_types_v3.pdf")
    p1 = to_plot %>%
#        mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
        filter(roi %in% my_roi) %>%
        mutate(roi = factor(roi, my_roi)) %>%
        group_by(all_types) %>%
        complete(roi, type, fill = list(n=0,fdr=1,log_OR=0)) %>%
        ungroup() %>%
        mutate(type = factor(type, rev(ct_order))) %>%
        ggplot(aes(x=roi,y=type)) +
        geom_point(aes(size=jaccard, fill = log_OR, col=fdr<0.05), pch=21) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_size(range = c(0,8)) +
        scale_fill_gradient2(low = "darkblue", mid = "gray90", high = "darkorange2") +
        scale_color_manual(values=c(gray(1,0),"black")) +
        labs(x = "Brain area", y = NULL)
    p1
    ggsave("figs/knn_mapping/yao_mapping_et_types_v4.pdf")
    p1 = to_plot %>%
#        mutate(log_OR = pmin(pmax(log_OR, -3), 3)) %>%
        filter(roi %in% my_roi) %>%
        mutate(roi = factor(roi, my_roi)) %>%
        group_by(all_types) %>%
        complete(roi, type, fill = list(n=0,fdr=1,log_OR=0)) %>%
        ungroup() %>%
        mutate(type = factor(type, rev(ct_order))) %>%
        ggplot(aes(x=roi,y=type)) +
        geom_point(aes(size=f_cells, fill = log_OR, col=fdr<0.05), pch=21) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust=1)) +
        scale_size_area(max_size = 8) +
        scale_fill_gradient2(low = "darkblue", mid = "gray90", high = "darkorange2") +
        scale_color_manual(values=c(gray(1,0),"black")) +
        labs(x = "Brain area", y = NULL)
    p1
    ggsave("figs/knn_mapping/yao_mapping_et_types_v5.pdf")
}

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

predict_cell_type = function(neighbors, ref_label, query_id, query_label) {
    labels = matrix(ref_label[neighbors], nrow(neighbors))
    pred = as_tibble(labels) %>%
        mutate(sample = query_id) %>%
        pivot_longer(-sample, values_to = "prediction", names_to = NULL) %>%
        group_by(sample, prediction) %>%
        tally() %>%
        group_by(sample) %>%
        summarize(predicted_ct = prediction[which.max(n)],
                  confidence = n[which.max(n)] / sum(n), .groups = "drop") %>%
        inner_join(tibble(sample = query_id, cell_type = query_label)) %>%
        select(sample, cell_type, predicted_ct, confidence)
    return(pred)
}

order_rows_according_to_cols = function(M, alpha = 1) {
    M <- M**alpha
    row_score <- colSums(t(M)*seq_len(ncol(M)), na.rm=TRUE)/rowSums(M, na.rm=TRUE)
    return(order(row_score))
}

simplify_roi = function(roi) {
    result = roi
    result[startsWith(roi, "ACA")] = "ACA"
    result[startsWith(roi, "AI")] = "AI"
    result[startsWith(roi, "AUD")] = "AUD"
    result[startsWith(roi, "ENT")] = "ENT"
    result[startsWith(roi, "RSP")] = "RSP"
    result[startsWith(roi, "SSp")] = "SSp"
    result[startsWith(roi, "VIS")] = "VIS"
    result[roi %in% "CA1"] = "HIP"
    result[roi %in% c("MOs", "FRP")] = "MOs-FRP"
    result[roi %in% c("PL", "ILA") | startsWith(roi, "ORB")] = "PL-ILA-ORB"
    result[roi %in% c("SSs","GU","VISC","AIp")] = "SSs-GU-VISC-AIp"
    result[roi %in% c("TEa","PERI","ECT")] = "TEa-PERI-ECT"
    return(result)
}

roi_to_plot = function() {
    c("ENT", "PIR", "TEa-PERI-ECT", "RSP", "ACA", "AI", "PL-ILA-ORB", "MOs-FRP", "MOp", "SSs-GU-VISC-AIp", "SSp", "AUD", "PTLp", "VIS")
}

compute_roi_enrichment = function(roi_data) {
    roi_total = roi_data %>%
        group_by(roi) %>%
        summarize(n = sum(n), .groups = "drop") %>%
        deframe()
    cluster_total = roi_data %>%
        group_by(cluster) %>%
        summarize(n = sum(n), .groups = "drop") %>%
        deframe()
    N = sum(roi_data$n)
    # Each test is reduced to a 2x2 contingency table:
    #         ROI_j  ROI_-j
    # CT_i  |   n  |  Ci-n |
    # CT_-i | Rj-n |N-Ci-Rj|
    #  -> OR = n*(N-Ci-Rj) / [(Ci-n)*(Rj-n)]
    #  -> white drawn = n, total drawn = Ci, #white= Rj, #black = N-Rj 
    #  -> Jaccard Index = n / (Ci+Rj-n)
    result = roi_data %>%
        mutate(Rj = roi_total[roi], Ci = cluster_total[cluster]) %>%
        # phyper() computes the proportion of scenarios that are *strictly* worse,
        # the -0.5 correction makes sure that the equally worse scenario is included
        mutate(log_pval = phyper(n-0.5, Rj, N-Rj, Ci, lower.tail=FALSE, log.p = TRUE),
               log_OR = log(n) + log(N-Ci-Rj) - log(Ci-n) - log(Rj-n),
               jaccard = n / (Ci+Rj-n)) %>%
        select(-Rj, -Ci)
    return(result)
}

#reannotate_barseq = function(cluster) {
#    new_label = c(
#        "BLA_L6M_1" = "L5_IT-like_BLA",
#        "BLA_L6M_2" = "CT-like_PIRI",
#        "BLA_L6M_3" = "CT_L6b-like_M",
#        "BLA_L6M_4" = "L6_IT-like_BLA",
#        "BLA_L6M_5" = "L6_IT-like_M_1",
#        "BLA_L6M_6" = "L6_IT-like_M_2",
#        "BLA_L6M_7" = "L6b-like_M"
#    )
#    is_changing = cluster %in% names(new_label)
#    result = cluster 
#    result[is_changing] = new_label[result[is_changing]]
#    return(result)
#}

subclass_mapping = function() {
    list("L2/3 IT"="L2/3 IT CTX", "L4/5 IT"="L4/5 IT CTX",
         "L5 IT"="L5 IT CTX", "RSP UL"="L2/3 IT PPP",
         "RSP DL"=c("L4 RSP-ACA","L5 IT CTX"),
         "L6 IT"="L6 IT CTX","PT"="L5 PT CTX",
         "NP"=c("NP SUB","NP PPP","L5 NP CTX"),
         "CT"="L6 CT CTX", "L6b"=c("L6b/CT ENT","L6b CTX"),
         "Car3"="Car3")
}

if (sys.nframe() == 0) {
    main()
}
