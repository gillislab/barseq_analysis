# https://bioconductor.org/books/release/OSCA/clustering.html

library(tidyverse)
library(SingleCellExperiment)
source("dataset.R")


main = function() {
    barseq = load_barseq(normalization_factor = 10)
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    
    #analyze_whole(barseq)
    #analyze_glu(barseq)
    analyze_gaba(barseq)
    #subcluster_glu(barseq)
}

analyze_whole = function(barseq) {
    analyze_barseq(barseq, "whole", k_umap=15, k_cluster=15)
}

analyze_glu = function(barseq) {
    barseq = filter_barseq(barseq, "analysis/whole")
    analyze_barseq(barseq[, startsWith(barseq$label, "GLU")], "glu", k_umap=15, k_cluster=15)
}

analyze_gaba = function(barseq) {
    barseq = filter_barseq(barseq, "analysis/whole")
    analyze_barseq(barseq[, startsWith(barseq$label, "GABA")], "gaba", k_umap=15, k_cluster=15)
}

subcluster_glu = function(barseq) {
    set.seed(17)
    output_dir = file.path("analysis", "subglu")
    dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)
    n_pca=30
    k_umap=15
    k_cluster=15
    
    file.copy("analysis/glu/umap.csv", output_dir)
    barseq = filter_barseq(barseq, "analysis/glu")
    labels = sapply(strsplit(barseq$label, split = " ", fixed=TRUE), "[", 1)
    clusters = list()
    cluster_annotation = list()
    umap=list()
    for (l in sort(unique(labels))) {
        subdata = barseq[, labels==l]
        subdata = scater::runPCA(subdata, ncomponents=n_pca, exprs_values="logcounts")    
        subdata = scater::runUMAP(subdata, dimred="PCA", n_neighbors=k_umap)
        result = cluster_data(subdata, k_cluster, "jaccard")    
        result$label = paste0(l, "_", result$label)
        clusters[[l]] = result
        cluster_annotation[[l]] = tibble(cluster_id=sort(unique(result$label)), cluster_name=cluster_id)
        umap[[l]] = as.data.frame(reducedDim(subdata, "UMAP")) %>%
            mutate(sample = colnames(subdata))
    }
    clusters = bind_rows(clusters)
    clusters = clusters[match(colnames(barseq), clusters$sample),]
    cluster_annotation = bind_rows(cluster_annotation)
    umap = bind_rows(umap)
    umap = umap[match(colnames(barseq), umap$sample),1:2]
    write_csv(clusters, file.path(output_dir, "cluster.csv")) 
    write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
    write_csv(umap, file.path(output_dir, "subumap.csv"))
    
    # TODO: compute markers
}

filter_barseq = function(barseq, input_dir) {
    clusters = read_csv(file.path(input_dir, "cluster.csv"))
    cluster_annotation = deframe(read_csv(file.path(input_dir, "cluster_annotation.csv")))
    barseq = barseq[,as.character(clusters$sample)]
    barseq$label = cluster_annotation[clusters$label]
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
    #label = leiden::leiden(snn)
    #label = igraph::cluster_walktrap(g)$membership
    label = igraph::cluster_louvain(snn)$membership
    result = tibble(sample = colnames(barseq), label)
    # only for walktrap
    #table(igraph::cut_at(community.walktrap, n=5))
    return(result)
}

subcluster_barseq = function(barseq, n_pca=30, k_cluster=15) {
    return(clusters)
}

if (sys.nframe() == 0) {
    main()
}
