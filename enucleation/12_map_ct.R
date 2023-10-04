
library(tidyverse)
library(hdf5r)
source("dataset.R")


main = function() {
    # load datasets
    current = load_current(); gc()
    colLabels(current) = current$h3_name
    
    # vs original
    ref = load_reference(); gc()
    pred = predict_cell_type(current, ref)
    write_csv(pred, "analysis/predicted_h3_ref.csv")
    p = plot_predictions(current, ref, pred, ref_ct_order = cell_type_order_ref())
    ggsave("fig/mapping_ref.pdf", p)
    
    # vs cheng
    cheng = load_cheng22(); gc()
    cheng$h3 = cheng$Type # Type = final ("harmonized" cell type), cluster = initial type
    cheng$h2 = cheng$Subclass # same story Subclass vs subclass
    colLabels(cheng) = cheng$h3
    # using barseq as reference
    pred = predict_cell_type(cheng, current)
    write_csv(pred, "analysis/predicted_h3_cheng22.csv")
    p = plot_predictions(cheng, current, pred)
    ggsave("fig/mapping_cheng22.pdf", p)    
    # visp only, barseq as reference
    colLabels(cheng) = cheng$Type
    pred = predict_cell_type(cheng, visp)
    write_csv(pred, "analysis/predicted_h3_cheng22_visp.csv")
    p = plot_predictions(cheng, visp, pred)
    ggsave("fig/mapping_cheng22_visp.pdf", p)    
}

load_current = function() {
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

load_reference = function(label_file = "../clustering/data/ref_labels.csv") {
    barseq = load_ref()
    barseq$sample = "reference"
    clusters = read_csv(label_file)
    clusters = filter(clusters, subclass %in% ctx_subclass_ref())
    barseq = barseq[, as.character(clusters$sample)]
    barseq$h2 = as.factor(clusters$subclass)
    barseq$h3 = as.factor(clusters$cluster)
    colLabels(barseq) = barseq$h3
    return(barseq)
}

load_cheng22 = function() {
    h5_file = H5File$new("data/cheng22_p28.h5ad", mode="r")
    genes = readDataSet(h5_file[["var/_index"]])
    samples = readDataSet(h5_file[["obs/_index"]])
    counts = Matrix::sparseMatrix(
        i=readDataSet(h5_file[["raw/X/indices"]])+1,
        p=readDataSet(h5_file[["raw/X/indptr"]]),
        x=readDataSet(h5_file[["raw/X/data"]]),
        dims=c(length(genes), length(samples)),
        dimnames = list(genes, samples)
    )
    metadata = read.csv("data/cheng22_p28_metadata.csv", row.names = 1)
    sce = SingleCellExperiment(assays = list(counts = counts), colData = metadata)
    return(sce)
}

predict_cell_type = function(query, ref) {        
    # common set + normalize data + rescale genes (optional)
    blacklist = c("Slc17a7", "Gad1", "unused-1","unused-2","unused-3","unused-4","unused-5")
    gene_set = intersect(rownames(query), rownames(ref))
    gene_set = gene_set[!gene_set %in% blacklist]
    f_gene = rowMeans(counts(query)[gene_set,]) / rowMeans(counts(ref)[gene_set,])
    query_data = normalize_cols(counts(query)[gene_set,])
    ref_data = normalize_cols(counts(ref)[gene_set,]*f_gene)
    #barseq_data = normalize_cols(log1p(convert_to_cpm(counts(barseq)[gene_set,], 1000)), ranked = FALSE)
    #train_data = normalize_cols(log1p(convert_to_cpm(counts(dataset)[gene_set,], 1000)), ranked = FALSE)
    #barseq_data = log1p(convert_to_cpm(counts(barseq)[gene_set,], 1000))
    #train_data = log1p(convert_to_cpm(counts(dataset)[gene_set,], 1000))
    
    # make predictions
    neighbors = find_neighbors(ref_data, query_data, 15) # MNN instead?, ~5 min
    pred = transfer_labels(neighbors, colLabels(ref), colnames(query), colLabels(query))
    #hist(pred$confidence)
    return(pred)
}

plot_predictions = function(query, ref, pred, ref_ct_order = NULL) {
    pred$n_reads = colSums(counts(query))
    pred$n_genes = colSums(counts(query)>0)
    pred_summary = pred %>%
#        mutate(predicted_ct = sample(predicted_ct)) %>% # used to compute null jaccard indices
        group_by(cell_type, predicted_ct) %>%
        tally() %>%
        ungroup()
    pred_matrix = pred_summary %>%
        pivot_wider(id_cols=c(cell_type), names_from = predicted_ct, values_from = n, values_fill = 0) %>%
        column_to_rownames("cell_type") %>%
        as.matrix()
        
    # overlap with reference (what fraction of target cells predicted as each ref type)?
    # summary_matrix = ctx_pred_matrix / rowSums(ctx_pred_matrix)
    # Jaccard coefficient
    summary_matrix = pred_matrix / (outer(rowSums(pred_matrix), colSums(pred_matrix), "+") - pred_matrix)
    summary_matrix = summary_matrix[sort(rownames(summary_matrix)),]
    # hist(summary_matrix, breaks = 50) # null jaccard never reaches 0.04
    # quantile(summary_matrix, probs = c(0.5, 0.9, 0.95, 0.99, 0.999, 1)) # null = 0.001, 0.005, 0.008, 0.015, 0.025, 0.04

    # summary plot for ctx types (only consider reference type with a hit > threshold)
    threshold = 0.05
    to_plot = t(summary_matrix)
    to_plot = to_plot[, colSums(to_plot>threshold)>0, drop=FALSE]
    if (!is.null(ref_ct_order)) { to_plot = to_plot[ref_ct_order,] }
    col_order = order_rows_according_to_cols(t(to_plot), 100)
    p1 = as_tibble(to_plot, rownames = "cluster") %>%
        pivot_longer(-cluster, names_to = "predicted_ct", values_to = "f") %>%
        mutate(cluster = factor(cluster, rownames(to_plot))) %>%
        mutate(predicted_ct = factor(predicted_ct, colnames(to_plot)[rev(col_order)])) %>%
        ggplot(aes(cluster, predicted_ct, size=ifelse(f==0,NA,f), col=f)) +
        geom_point() +
        scale_color_gradient(low = "gray90", high = "blue4") +
        scale_size_area(max_size=3) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, size=6)) +
        theme(axis.text.y = element_text(size=6), legend.position = "bottom") +
        labs(x=NULL,y=NULL)
    return(p1)
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

transfer_labels = function(neighbors, ref_label, query_id, query_label) {
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
    result = roi_data %>%
        mutate(Rj = roi_total[roi], Ci = cluster_total[cluster]) %>%
        # phyper() computes the proportion of scenarios that are *strictly* worse,
        # the -0.5 correction makes sure that the equally worse scenario is included
        mutate(log_pval = phyper(n-0.5, Rj, N-Rj, Ci, lower.tail=FALSE, log.p = TRUE),
               log_OR = log(n) + log(N-Ci-Rj) - log(Ci-n) - log(Rj-n)) %>%
        select(-Rj, -Ci)
    return(result)
}

if (sys.nframe() == 0) {
    main()
}
