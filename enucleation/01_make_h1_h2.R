# https://satijalab.org/seurat/articles/seurat5_sketch_analysis.html
# https://satijalab.org/seurat/articles/parsebio_sketch_integration

library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
source("dataset.R")
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")


main = function() {
    recompute = FALSE
    if (recompute) {
        combine_all()
        gc()
        make_h1() # ~5h with most efficient steps (25 minutes clustering)
        h1_labels = read_csv("data/clusters_leiden.csv.gz")
        coi = h1_labels$cell_id[h1_labels$label %in% c(4,5,6,7,9)] # glutamatergic clusters
        make_h2(coi, "analysis/glu") # ~2h30 total (10 minute clustering)
    }
}

combine_all = function() {
    all_counts = lapply(sample_names(), load_sample_counts)
    all_counts = do.call(cbind, all_counts)
    gc()
    
    # SCE format
    sce = SingleCellExperiment(assays = list(counts=all_counts, logcounts = log1p(convert_to_cpm(all_counts, 10)))) # 7M cells
    saveRDS(sce, "data/sce.rds")
    rm(sce)
    gc()
    
    # gather metadata from individual datasets
    my_samples = sample_names()
    h1_labels = select(bind_rows(lapply(my_samples, load_sample_annotation, annotation_level = "whole")), sample, h1=cluster_name)
    h2_labels = rbind(
        select(bind_rows(lapply(my_samples, load_sample_annotation, annotation_level = "glu")), sample, h2=cluster_name),
        select(bind_rows(lapply(my_samples, load_sample_annotation, annotation_level = "gaba")), sample, h2=cluster)
    )
    h3_labels = select(bind_rows(lapply(my_samples, load_sample_annotation, annotation_level = "subglu")), sample, h3=cluster)
    labels = left_join(h1_labels, h2_labels) %>% left_join(h3_labels)
    write_csv(labels, "data/annotation.csv.gz")
    
    # gather cell positions from individual datasets
    cellpos = lapply(set_names(sample_names()), function(s) {
        result = as_tibble(colData(load_sample(s)), rownames = "cell_id")
    })
    cellpos = bind_rows(cellpos, .id = "sample")
    cellpos = select(cellpos, -angle, -depth_x, -depth_y)
    write_csv(cellpos, "data/cell_position.csv.gz")
}

make_h1 = function(n_pca = 30, k_cluster = 15, k_umap = 15, resolution=0.2) {
    ## Bioconductor analysis
    if (!file.exists("data/pca.rds")) {
        sce = readRDS("data/sce.rds")
        # ~ 2min30 for 9.1M cells
        pca = scater::calculatePCA(sce, ncomponents=n_pca, exprs_values="logcounts")
        saveRDS(pca, "data/pca.rds")
        rm(pca); rm(sce); gc()
    }
    if (!file.exists("data/umap.csv.gz")) {
        pca = readRDS("data/pca.rds")
        # uses uwot::umap internally -> check differences with Seurat, use uwot-learn directly?
        #  - main difference with seurat "euclidean" vs "cosine"?
        #  - annoy is bad at finding true nearest neighbors -> https://jlmelville.github.io/uwot/pycompare
        # multiple options -> external_neighbor + annoy, multiple threads etc.
        #  - option 1 = single thread + annoy -> ~3h, up to 35 GB (swapping)
        umap = scater::calculateUMAP(
            pca, transposed=TRUE, n_neighbors=k_umap,
            external_neighbors=TRUE, BNPARAM = BiocNeighbors::AnnoyParam()
        )
        #  - option 2 = default parameters -> XXX
        #umap = scater::calculateUMAP(pca, transposed=TRUE, n_neighbors=k_umap)
        #  - option 3 = 16 threads + "native" NN search -> XXX
        #umap = scater::calculateUMAP(pca, transposed=TRUE, n_neighbors=k_umap, n_threads=16)
        colnames(umap) = c("UMAP1","UMAP2")
        write_csv(as_tibble(umap, rownames = "cell_id"), "data/umap.csv.gz")
        rm(umap); rm(pca); gc()
    }
        
    ## Seurat (classic / sketching / downsampling)
    # classic
    if (!file.exists("data/snn_srt.rds")) {
        pca = readRDS("data/pca.rds")[,1:10]
        gc()
        # ~ 1h30, up to 80GB for 9.1M cells
        snn = FindNeighbors(pca, k.param = k_cluster) # /!\ note: 10 PCs by default
        saveRDS(snn, "data/snn_srt.rds")
        rm(pca); rm(snn); gc()
    }
    if (!file.exists("data/clusters_leiden.csv.gz")) {
        snn = as(readRDS("data/snn_srt.rds")$snn, "dgCMatrix")
        snn = igraph::graph_from_adjacency_matrix(snn, weighted=TRUE, mode="undirected")
        gc()
        # ~ 26minutes, up to 50GB for 9.1M cells
        message("Leiden igraph: ", paste(round(system.time({
        clusters = igraph::cluster_leiden(snn, objective_function='modularity', resolution = resolution)
        })), collapse = " "))
        result = tibble(cell_id = clusters$name, label = clusters$membership)
        write_csv(result, "data/clusters_leiden.csv.gz")
        rm(clusters); rm(result); rm(snn); gc()
    }
}

make_h2 = function(coi, output_dir, n_pca = 30, k_cluster = 15, k_umap = 15, resolution=0.8) {
    message("Running H2 level analysis for ", length(coi), " cells.")
    dir.create(output_dir, showWarnings=FALSE)
    if (!file.exists(file.path(output_dir, "pca.rds"))) {
        sce = readRDS("data/sce.rds")[, coi]; gc()
        # 1min30 for 4M cells (cluster 1 core)
        message("PCA: ", paste(round(system.time({
            pca = scater::calculatePCA(sce, ncomponents=n_pca, exprs_values="logcounts")
        })), collapse=" "))
        saveRDS(pca, file.path(output_dir, "pca.rds"))
        rm(pca); rm(sce); gc()
    }
    if (!file.exists(file.path(output_dir, "umap.csv.gz"))) {
        pca = readRDS(file.path(output_dir, "pca.rds")); gc()
        # 1h40 for 4M cells (cluster 1 core)
        message("UMAP1: ", paste(round(system.time({
            umap = scater::calculateUMAP(
                pca, transposed=TRUE, n_neighbors=k_umap,
                external_neighbors=TRUE, BNPARAM = BiocNeighbors::AnnoyParam()
            )
        })), collapse=" "))
        colnames(umap) = c("UMAP1","UMAP2")
        write_csv(as_tibble(umap, rownames = "cell_id"), file.path(output_dir, "umap.csv.gz"))
        rm(umap); gc()
        rm(pca); gc()
    }
    if (!file.exists(file.path(output_dir, "snn_srt.rds"))) {
        pca = readRDS(file.path(output_dir, "pca.rds")); gc()
        # ~ 40min for 4M cells (10PCs, cluster 1 core)
        # ~ 50min for 4M cells (30PCs, laptop)
        message("SNN: ", paste(round(system.time({
            snn = FindNeighbors(pca, k.param = k_cluster)
        })), collapse=" "))
        saveRDS(snn, file.path(output_dir, "snn_srt.rds"))
        rm(pca); rm(snn); gc()
    }
    if (!file.exists(file.path(output_dir, "cluster.csv.gz"))) {
        snn = as(readRDS(file.path(output_dir, "snn_srt.rds"))$snn, "dgCMatrix")
        snn = igraph::graph_from_adjacency_matrix(snn, weighted=TRUE, mode="undirected")
        gc()
        # ~2h for 4M cells (30PCs, -1 iteratons, laptop)
        message("Leiden: ", paste(round(system.time({
            clusters = igraph::cluster_leiden(snn, objective_function='modularity', resolution_parameter = resolution, n_iterations = -1)
        })), collapse=" "))
        result = tibble(cell_id = clusters$name, label = clusters$membership)
        write_csv(result, file.path(output_dir, "cluster.csv.gz"))
        rm(clusters); rm(result); rm(snn); gc()
    }
}
 
if (sys.nframe() == 0) {
    main()
}