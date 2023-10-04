# https://bioconductor.org/books/release/OSCA/clustering.html

library(tidyverse)
library(SingleCellExperiment)
source("dataset.R")


main = function() {
    barseq = load_barseq()
    subcluster_glu(barseq)
}

subcluster_glu = function(barseq) {
    set.seed(17)
    output_dir = file.path("analysis", "subglu")
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)
    n_pca=30
    k_umap=15
    k_cluster=15
    
    file.copy("analysis/glu/umap.csv.gz", output_dir)
    barseq = filter_barseq(barseq, "analysis/glu")
    labels = sapply(strsplit(barseq$label, split = " ", fixed=TRUE), "[", 1)
    clusters = list()
    cluster_annotation = list()
    umap=list()
    for (l in sort(unique(labels))) {
        message("Computing H3 clusters for ", l, " (", sum(labels==l), "cells).")
        subdata = barseq[, labels==l]
        gc()
        message("PCA: ", paste(round(system.time({
        subdata = scater::runPCA(subdata, ncomponents=n_pca, exprs_values="logcounts")
        })), collapse = " "))
        gc()
        message("UMAP: ", paste(round(system.time({
        subdata = scater::runUMAP(subdata, dimred="PCA", n_neighbors=k_umap)
        })), collapse = " "))
        gc()
        message("Clustering: ", paste(round(system.time({
        result = cluster_data(subdata, k_cluster, "jaccard")    
        })), collapse = " "))
        gc()
        result$label = paste0(l, "_", result$label)
        clusters[[l]] = result
        cluster_annotation[[l]] = tibble(cluster_id=sort(unique(result$label)), cluster_name=cluster_id)
        umap[[l]] = as.data.frame(reducedDim(subdata, "UMAP")) %>%
            mutate(cell_id = colnames(subdata))
    }
    clusters = bind_rows(clusters)
    clusters = clusters[match(colnames(barseq), clusters$sample),]
    cluster_annotation = bind_rows(cluster_annotation)
    umap = bind_rows(umap)
    umap = umap[match(colnames(barseq), cell_id$cell_id),]
    colnames(umap) = c("UMAP1","UMAP2","cell_id")
    write_csv(clusters, file.path(output_dir, "cluster.csv.gz")) 
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
    write_csv(umap, file.path(output_dir, "subumap.csv.gz"))
}

filter_barseq = function(barseq, input_dir, min_cluster_size=10) {
    clusters = read_csv(file.path(input_dir, "cluster.csv.gz"))
    cluster_annotation = deframe(read_csv(file.path(input_dir, "cluster_annotation.csv")))
    # keep only clusters not annotated as "NA"
    is_not_na = names(cluster_annotation)[!is.na(cluster_annotation)]
    clusters = filter(clusters, label %in% is_not_na)
    barseq = barseq[,as.character(clusters$cell_id)]
    barseq$label = cluster_annotation[as.character(clusters$label)]
    return(barseq)
}

analyze_barseq = function(barseq, output_name, n_pca=30, k_umap=100, k_cluster=50) {
    output_dir = file.path("analysis", make.names(output_name))
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)

    barseq = scater::runPCA(barseq, ncomponents=n_pca, exprs_values="logcounts")
    barseq = scater::runUMAP(barseq, dimred="PCA", n_neighbors=k_umap)
    write_csv(as.data.frame(reducedDim(barseq, "UMAP")), file.path(output_dir, "umap.csv"))
    
    clusters = cluster_data(barseq, k_cluster, "jaccard")
    write_csv(clusters, file.path(output_dir, "cluster.csv")) 
    
    # create default cluster annotation
    labels = as.factor(clusters$label)
    cluster_annotation = tibble(cluster_id = levels(labels)) %>%
        mutate(cluster_name = paste0(output_name, "_", cluster_id))
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
}

cluster_data = function(barseq, k=10, type="rank") {
    # type -> "rank" (preservation of neighbor rank),
    # "number" (number of shared neighbors), "jaccard"
    snn = scran::buildSNNGraph(barseq, k=k, use.dimred = "PCA", type=type,
                               BNPARAM = BiocNeighbors::AnnoyParam(),
                               BPPARAM = BiocParallel::MulticoreParam(16))
    gc()
    #label = leiden::leiden(snn)
    #label = igraph::cluster_walktrap(g)$membership
    label = igraph::cluster_louvain(snn)$membership
    result = tibble(sample = colnames(barseq), label)
    # only for walktrap
    #table(igraph::cut_at(community.walktrap, n=5))
    return(result)
}

if (sys.nframe() == 0) {
    main()
}
