
library(tidyverse)
source("dataset.R")


main = function() {
    # reference dataset
    dataset = load_yao21("scC", "MOp")
    dataset = dataset[, dataset$class_label == "Glutamatergic"]
    dataset = dataset[, grepl("CTX", dataset$subclass_label) | grepl("Car3", dataset$subclass_label)]
    n_cells = tabulate(as.factor(dataset$cluster_label))
    ct_to_keep = levels(as.factor(dataset$cluster_label))[n_cells>=50]
    dataset = dataset[, dataset$cluster_label %in% ct_to_keep]
    logcounts(dataset) = log1p(convert_to_cpm(counts(dataset), 10000))

    # BARseq panel
    barseq = load_barseq()
    clusters = read_csv("analysis/labels.csv")
    assertthat::are_equal(clusters$sample, colnames(barseq))
    barseq$cluster = clusters$cluster
    barseq$subclass = clusters$subclass
    barseq = barseq[, clusters$class == "Glutamatergic"]
    blacklist = c("Slc17a7", "Gad1")
    barseq_panel = intersect(rownames(barseq), rownames(dataset))
    barseq_panel = barseq_panel[!barseq_panel %in% blacklist]

    #assess_predictability(barseq, dataset, barseq_panel)
    plot_results(barseq, dataset, barseq_panel)
}

assess_predictability = function(barseq, dataset, barseq_panel) {    
    is_mop = !is.na(barseq$CCFparentname) & barseq$CCFparentname == "MOp"

    # HVG space
    dec = scran::modelGeneVar(dataset)
    hvg = scran::getTopHVGs(dec, n=2000)
    my_data = logcounts(dataset)[hvg,]
    neighbors = find_neighbors(my_data, my_data, 16) # MNN instead?, ~8 minutes
    pred_subclass = predict_cell_type(neighbors[,-1], dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_hvg_subclass.csv")
    pred_cluster = predict_cell_type(neighbors[,-1], dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_hvg_cluster.csv")

    # PCA space
    dataset = scater::runPCA(dataset, ncomponents=30, exprs_values="logcounts", subset_row=hvg)
    my_data = t(reducedDim(dataset, "PCA"))
    neighbors = find_neighbors(my_data, my_data, 16) # MNN instead?
    pred_subclass = predict_cell_type(neighbors[,-1], dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_pca_subclass.csv")
    pred_cluster = predict_cell_type(neighbors[,-1], dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_pca_cluster.csv")

    # barseq space
    my_data = logcounts(dataset)[barseq_panel,]
    neighbors = find_neighbors(my_data, my_data, 16) # MNN instead?
    pred_subclass = predict_cell_type(neighbors[,-1], dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_barseq_panel_subclass.csv")
    pred_cluster = predict_cell_type(neighbors[,-1], dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_barseq_panel_cluster.csv")

    # barseq PCA space
    dataset = scater::runPCA(dataset, ncomponents=30, exprs_values="logcounts", subset_row=barseq_panel)
    my_data = t(reducedDim(dataset, "PCA"))
    neighbors = find_neighbors(my_data, my_data, 16) # MNN instead?
    pred_subclass = predict_cell_type(neighbors[,-1], dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_barseq_pca_subclass.csv")
    pred_cluster = predict_cell_type(neighbors[,-1], dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_barseq_pca_cluster.csv")
    
    # barseq downsampled uniformly -> keep ~ 1 read out of 7
    #  - mean(counts(dataset)[barseq_panel,]) -> ~3.7 UMIs
    #  - mean(counts(barseq)[barseq_panel, is_mop]) -> 0.52 UMIs
    #  - mean(counts(barseq)[barseq_panel,]) / mean(counts(dataset)[barseq_panel,]) -> 0.14 
    set.seed(17)
    query_data = downsample_data(counts(dataset)[barseq_panel,], 1/7)
    query_data = normalize_cols(log1p(convert_to_cpm(query_data, 10)))
    ref_data = normalize_cols(log1p(convert_to_cpm(counts(dataset)[barseq_panel,], 10)))
    neighbors = find_neighbors(ref_data, query_data, 16)
    neighbors = remove_self(neighbors)
    pred_subclass = predict_cell_type(neighbors, dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_downsampling_unif_subclass.csv")
    pred_cluster = predict_cell_type(neighbors, dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_downsampling_unif_cluster.csv")
    
    # barseq downsampled by gene
    gene_counts = rowMeans(counts(barseq)[barseq_panel, is_mop])
    set.seed(17)
    query_data = downsample_data_by_gene(counts(dataset)[barseq_panel,], gene_counts)
    query_data = normalize_cols(query_data)
    #query_data = normalize_cols(log1p(convert_to_cpm(query_data, 1000)), ranked=FALSE)
    ref_data = normalize_cols(counts(dataset)[barseq_panel,])
    #ref_data = normalize_cols(logcounts(dataset)[barseq_panel,], ranked=FALSE)
    neighbors = find_neighbors(ref_data, query_data, 16)
    neighbors = remove_self(neighbors)
    pred_subclass = predict_cell_type(neighbors, dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_downsampling_genewise_subclass.csv")
    pred_cluster = predict_cell_type(neighbors, dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_downsampling_genewise_cluster.csv")
    # what if both query and ref are downsampled? -> rescale ref genes
    f_gene = target_gene_counts / rowMeans(count_matrix)
    ref_data = normalize_cols(counts(dataset)[barseq_panel,]*f_gene)
    neighbors = find_neighbors(ref_data, query_data, 16)
    neighbors = remove_self(neighbors)
    pred_subclass = predict_cell_type(neighbors, dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_downsampling_ref_subclass.csv")
    pred_cluster = predict_cell_type(neighbors, dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_downsampling_ref_cluster.csv")
    # what if we embed vectors in reference PCA space?
    # my_ref = prcomp(t(ref_data)) -> full PCA
    # all.equal(ref_data, t(my_ref$x %*% t(my_ref$rotation)) + my_ref$center) -> project into PC space
    # all.equal(my_ref$x, t(ref_data - rowMeans(ref_data)) %*% my_ref$rotation) -> project into original space
    ref_data = normalize_cols(counts(dataset)[barseq_panel,])
    pca = scater::calculatePCA(ref_data, ncomponents=30)
    rot_pca = attributes(pca)$rotation[rownames(ref_data),]
    ref_pca = t(t(ref_data - rowMeans(ref_data)) %*% rot_pca)
    # all.equal(ref_pca, pca[,])
    query_pca = t(t(query_data - rowMeans(query_data, na.rm = TRUE)) %*% rot_pca)
    neighbors = find_neighbors(ref_pca, query_pca, 16)
    neighbors = remove_self(neighbors)
    pred_subclass = predict_cell_type(neighbors, dataset$subclass_label, colnames(dataset), dataset$subclass_label)
    write_csv(pred_subclass, "analysis/knn_mapping/knn_downsampling_pca_subclass.csv")
    pred_cluster = predict_cell_type(neighbors, dataset$cluster_label, colnames(dataset), dataset$cluster_label)
    write_csv(pred_cluster, "analysis/knn_mapping/knn_downsampling_pca_cluster.csv")
}

plot_results = function(barseq, dataset, barseq_panel) {
    is_mop = !is.na(barseq$CCFparentname) & barseq$CCFparentname == "MOp"

    # summary plots
    # "hvg",  "barseq_pca", "downsampled_unif", "downsampled_gene", "downsampled_qc", 
    gene_set_order = c("hvg_pca"="Highly variable genes\n(2000 genes, 30 PCs)",
                       "barseq_panel"="BARseq panel\n(107 genes)",
                       "downsampled_ref"="BARseq panel\n(matched sensitivity)",
                       "permuted_labels"="Permuted labels")
    gene_counts = rowMeans(counts(barseq)[barseq_panel, is_mop])
    set.seed(17)
    downsampled_data = downsample_data_by_gene(counts(dataset)[barseq_panel,], gene_counts)
    passes_qc = (colSums(downsampled_data) >= 20) & (colSums(downsampled_data>0)>=5)
    set.seed(17)
    permuted_labels = tibble(cell_type = dataset$subclass_label, predicted_ct = sample(cell_type))
    subclass_results = bind_rows(list(
        hvg = compute_accuracy(read_csv("analysis/knn_mapping/knn_hvg_subclass.csv")),
        hvg_pca = compute_accuracy(read_csv("analysis/knn_mapping/knn_pca_subclass.csv")),
        barseq_panel = compute_accuracy(read_csv("analysis/knn_mapping/knn_barseq_panel_subclass.csv")),
        barseq_pca = compute_accuracy(read_csv("analysis/knn_mapping/knn_barseq_pca_subclass.csv")),
        downsampled_unif = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_unif_subclass.csv")),
        downsampled_gene = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_subclass.csv")),
        downsampled_ref = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_ref_subclass.csv")),
        downsampled_pca = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_pca_subclass.csv")),
        downsampled_qc = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_subclass.csv")[passes_qc,]),
        permuted_labels = compute_accuracy(permuted_labels)
    ), .id = "gene_set")
    p1 = subclass_results %>%
        mutate(gene_set = factor(gene_set_order[gene_set], gene_set_order)) %>%
        drop_na(gene_set) %>%
        ggplot(aes(x=gene_set,y=accuracy, fill = gene_set)) +
        geom_boxplot(show.legend = FALSE) +
        theme_bw(base_size = 20) +
        scale_fill_brewer(palette="Pastel1") +
        theme(axis.text.x = element_text(angle=30, hjust = 1)) +
        labs(x=NULL, y="Cell type predictability (H2 accuracy)")
    p1
    ggsave("figs/knn_mapping/knn_subclass_predictability.pdf", p1)

    set.seed(17)
    permuted_labels = tibble(cell_type = dataset$cluster_label, predicted_ct = sample(cell_type))
    cluster_results = bind_rows(list(
        hvg = compute_accuracy(read_csv("analysis/knn_mapping/knn_hvg_cluster.csv")),
        hvg_pca = compute_accuracy(read_csv("analysis/knn_mapping/knn_pca_cluster.csv")),
        barseq_panel = compute_accuracy(read_csv("analysis/knn_mapping/knn_barseq_panel_cluster.csv")),
        barseq_pca = compute_accuracy(read_csv("analysis/knn_mapping/knn_barseq_pca_cluster.csv")),
        downsampled_unif = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_unif_cluster.csv")),
        downsampled_gene = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_cluster.csv")),
        downsampled_ref = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_ref_cluster.csv")),
        downsampled_pca = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_pca_cluster.csv")),
        downsampled_qc = compute_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_cluster.csv")[passes_qc,]),
        permuted_labels = compute_accuracy(permuted_labels)
    ), .id = "gene_set")
    p1 = cluster_results %>%
        mutate(gene_set = factor(gene_set_order[gene_set], gene_set_order)) %>%
        drop_na(gene_set) %>%
        ggplot(aes(x=gene_set,y=accuracy, fill = gene_set)) +
        geom_boxplot(show.legend = FALSE) +
        theme_bw(base_size = 20) +
        scale_fill_brewer(palette = "Pastel1") +
        theme(axis.text.x = element_text(angle=30, hjust = 1)) +
        labs(x=NULL, y="Cell type predictability (H3 accuracy)")
    ggsave("figs/knn_mapping/knn_cluster_predictability.pdf", p1)

    # micro accuracy
    set.seed(17)
    permuted_labels = tibble(cell_type = dataset$subclass_label, predicted_ct = sample(cell_type))
    subclass_results = bind_rows(list(
        hvg = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_hvg_subclass.csv")),
        hvg_pca = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_pca_subclass.csv")),
        barseq_panel = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_barseq_panel_subclass.csv")),
        barseq_pca = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_barseq_pca_subclass.csv")),
        downsampled_unif = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_unif_subclass.csv")),
        downsampled_gene = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_subclass.csv")),
        downsampled_ref = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_ref_subclass.csv")),
        downsampled_pca = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_pca_subclass.csv")),
        downsampled_qc = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_subclass.csv")[passes_qc,]),
        permuted_labels = compute_micro_accuracy(permuted_labels)
    ), .id = "gene_set")
    p1 = subclass_results %>%
        mutate(gene_set = factor(gene_set_order[gene_set], gene_set_order)) %>%
        drop_na(gene_set) %>%
        ggplot(aes(x=gene_set,y=accuracy, fill = gene_set)) +
        geom_col(show.legend = FALSE) +
        geom_text(aes(label = round(accuracy,2)), size=6, nudge_y=-0.04) +
        theme_bw(base_size = 20) +
        scale_fill_brewer(palette = "Pastel1") +
        theme(axis.text.x = element_text(angle=30, hjust = 1)) +
        labs(x=NULL, y="Cell type predictability (H2 accuracy)") +
        lims(y=c(0,1))
    p1
    ggsave("figs/knn_mapping/knn_subclass_predictability_micro.pdf", p1)

    set.seed(17)
    permuted_labels = tibble(cell_type = dataset$cluster_label, predicted_ct = sample(cell_type))
    cluster_results = bind_rows(list(
        hvg = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_hvg_cluster.csv")),
        hvg_pca = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_pca_cluster.csv")),
        barseq_panel = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_barseq_panel_cluster.csv")),
        barseq_pca = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_barseq_pca_cluster.csv")),
        downsampled_unif = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_unif_cluster.csv")),
        downsampled_gene = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_cluster.csv")),
        downsampled_ref = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_ref_cluster.csv")),
        downsampled_pca = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_pca_cluster.csv")),
        downsampled_qc = compute_micro_accuracy(read_csv("analysis/knn_mapping/knn_downsampling_genewise_cluster.csv")[passes_qc,]),
        permuted_labels = compute_micro_accuracy(permuted_labels)
    ), .id = "gene_set")
    p1 = cluster_results %>%
        mutate(gene_set = factor(gene_set_order[gene_set], gene_set_order)) %>%
        drop_na(gene_set) %>%
        ggplot(aes(x=gene_set,y=accuracy, fill = gene_set)) +
        geom_col(show.legend = FALSE) +
        geom_text(aes(label = round(accuracy,2)), size=6, nudge_y=-0.04) +
        theme_bw(base_size = 20) +
        scale_fill_brewer(palette = "Pastel1") +
        theme(axis.text.x = element_text(angle=30, hjust = 1)) +
        labs(x=NULL, y="Cell type predictability (H3 accuracy)") +
        lims(y=c(0,1))
    p1
    ggsave("figs/knn_mapping/knn_cluster_predictability_micro.pdf", p1)
    
    # expected Jaccard index distribution?
    j_index = bind_rows(list(
        hvg_pca = compute_jaccard(read_csv("analysis/knn_mapping/knn_pca_cluster.csv")),
        barseq_panel = compute_jaccard(read_csv("analysis/knn_mapping/knn_barseq_panel_cluster.csv")),
        downsampled_ref = compute_jaccard(read_csv("analysis/knn_mapping/knn_downsampling_ref_cluster.csv")),
        permuted_labels = compute_jaccard(permuted_labels)
    ), .id = "gene_set")
    p1 = j_index %>%
        mutate(gene_set = factor(gene_set_order[gene_set], gene_set_order)) %>%
        drop_na(gene_set) %>%
        ggplot(aes(x=gene_set,y=jaccard, fill = gene_set)) +
        geom_boxplot(show.legend = FALSE) +
        theme_bw(base_size = 20) +
        scale_fill_brewer(palette = "Pastel1") +
        theme(axis.text.x = element_text(angle=30, hjust = 1)) +
        labs(x=NULL, y="H3 prediction self-overlap (Jaccard index)") +
        geom_hline(yintercept = 0.05, linetype="dashed") +
        lims(y=c(0,1))
    p1
    ggsave("figs/knn_mapping/knn_cluster_predictability_jaccard.pdf", p1)    
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

compute_accuracy = function(pred) {
    pred %>%
        group_by(cell_type) %>%
        summarize(accuracy = mean(predicted_ct == cell_type), n_cells = n(), .groups = "drop")
}

compute_micro_accuracy = function(pred) {
    pred %>%
        summarize(accuracy = mean(predicted_ct == cell_type))
}

downsample_data = function(count_matrix, f) {
    to_keep = rbinom(n=length(count_matrix@x), size = count_matrix@x, prob = f)
    result = count_matrix
    result@x = as.numeric(to_keep)
    return(result)
}

downsample_data_by_gene = function(count_matrix, target_gene_counts) {
    f_gene = target_gene_counts / rowMeans(count_matrix)
    f = f_gene[count_matrix@i+1]
    result = count_matrix
    result@x = result@x * floor(f)
    to_keep = rbinom(n=length(count_matrix@x), size = count_matrix@x, prob = f%%1)
    result@x = result@x + as.numeric(to_keep)
    return(result)
}

remove_self = function(neighbors) {
    last_index = ncol(neighbors)
    tneighbors = t(neighbors)
    result = neighbors
    for (i in 1:nrow(neighbors)) {
        self_index = which(tneighbors[,i]==i)
        result[i, self_index] = neighbors[i, last_index]
    }
    result = result[, -last_index]
    return(result)
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

compute_jaccard = function(pred) {
    confusion_matrix = pred %>%
        mutate(cell_type = as.factor(cell_type)) %>%
        mutate(predicted_ct = factor(predicted_ct, levels(cell_type))) %>%
        group_by(cell_type, predicted_ct) %>%
        tally() %>%
        ungroup() %>%
        pivot_wider(c(cell_type), names_from = predicted_ct, values_from = n, values_fill = 0, names_expand = TRUE) %>%
        column_to_rownames("cell_type") %>%
        as.matrix()
    j_index = confusion_matrix / (outer(rowSums(confusion_matrix), colSums(confusion_matrix), "+") - confusion_matrix)
    j_index = j_index[,rownames(j_index)]
    result = diag(j_index)
    names(result) = rownames(j_index)
    return(enframe(result, "cell_type", "jaccard"))
}

if (sys.nframe() == 0) {
    main()
}
