# https://satijalab.org/seurat/articles/seurat5_sketch_analysis.html
# https://satijalab.org/seurat/articles/parsebio_sketch_integration
# requires Seurat v5 and BPCells package
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# remotes::install_github("bnprks/BPCells", quiet = TRUE)

library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
#library(BPCells)
source("dataset.R")
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")


main = function() {
    recompute = FALSE
    if (recompute) {
        combine_all()
        gc()
        analyze_whole() # ~5h with most efficient steps (25 minutes clustering)
        annotate_h1()
        h1_labels = read_csv("data/clusters_leiden.csv.gz")
        coi = h1_labels$cell_id[h1_labels$label %in% c(4,5,6,7,9)] # glutamatergic clusters
        make_h2(coi, "analysis/glu") # ~2h30 total (10 minute clustering)
        annotate_h2()
    }
    h2_labels = read_csv("analysis/glu/clusters_louvain.csv.gz")
    ## this is a test to estimate clustering speed and choose clustering strategy at h3 level
    ## see 02_make_h3.R for actual clustering
    #coi = h2_labels$cell_id[h2_labels$label == 2]
    #make_h3(coi, "analysis/l45it") 
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

    # BPCell format (for Seurat analysis) -> sketching not ideal, clustering runs better on all cells
    #write_matrix_dir(mat=all_counts, dir="data/bpcell_counts")
    #rm(all_counts)
    #gc()
    
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

analyse_whole = function(n_pca = 30, k_cluster = 15, k_umap = 15, resolution=0.2) {
    # we’ll experiment 3 strategies
    #  - Bioconductor analysis on all cells
    #  - Seurat analysis on sketched subset
    #  - Seurat analysis on simple downsampled subset
    
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
    #if (!file.exists("data/snn_bioc.rds")) {
    #    pca = readRDS("data/pca.rds")[, 1:10]
    #    # pc=30: fails after 1h40 due to long vectors..., up to 100~GB (annoy - 1 core - laptop/cluster)
    #    # pc=10: also fails after 1h35 (killed)
    #    system.time({
    #    snn = bluster::makeSNNGraph(pca, k=k_cluster, type="jaccard", BNPARAM = BiocNeighbors::AnnoyParam())
    #    })
    #    gc()
    #    saveRDS(snn, "data/snn_bioc.rds")
    #    rm(pca); rm(snn); gc()
    #}
        
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
    if (!file.exists("data/clusters_louvain.csv.gz")) {
        snn = readRDS("data/snn_srt.rds")$snn
        gc()
        # ~ 5h40, up to 20GB for 9.1M cells
        system.time({
        clusters = FindClusters(snn, resolution = resolution)
        })
        result = tibble(cell_id = colnames(snn), label = clusters$XXX)
        write_csv(result, "data/clusters_louvain.csv.gz")
        rm(clusters); rm(result); rm(snn); gc()
    }
    if (!file.exists("data/clusters_leidenalg.csv.gz")) {
        snn = readRDS("data/snn_srt.rds")$snn
        gc()
        # ~ XX, up to 75GB (heavy swapping on laptop) for 9.1M cells
        message("Leiden srt: ", paste(round(system.time({
        clusters = FindClusters(snn, resolution = resolution, algorithm=4, method="igraph")
        })), collapse = " "))
        result = tibble(cell_id = colnames(snn), label = clusters$membership)
        write_csv(result, "data/clusters_leidenalg.csv.gz")
        rm(result); rm(snn); gc()
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
    
    # sketching
    # -> sketching not ideal, clustering runs better on all cells
    # it’s important to use BPCells at this step, but may not be necessary afterwards
    #srt = SketchData(srt, ncells = 1000000, method = "LeverageScore", sketched.assay = "sketch")
    #LayerData(srt, "scale.data") = as.matrix(LayerData(srt, "data")) # do *not* scale the data -> actually, should we scale?
    #srt = ScaleData(srt)
    # ~ 1 minute for 200K cells
    #srt = FindNeighbors(srt, k.param = k_cluster)
    # ~ 2 minutes for 200K cells
    #srt = FindClusters(srt, resolution=0.2)
    #DimPlot(srt, alpha=0.5)
    #FeaturePlot(srt, "sketch_Gad1", label=TRUE) & NoLegend()
    #FeaturePlot(srt, "sketch_Slc17a7", label=TRUE) & NoLegend()
    # ~ 10 minutes for 7M cells
    #system.time({
    #srt = ProjectData(srt, assay="RNA", full.reduction = "pca.full", sketched.assay = "sketch",
    #                  sketched.reduction = "pca", umap.model = "umap", dims = 1:n_pca,
    #                  refdata = list(cluster_full = "seurat_clusters"))
    #})
    #saveRDS(srt, "data/srt_whole.rds")
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
        # 1h55 for 4M cells (cluster 1 core)
        #message("UMAP2: ", paste(round(system.time({
        #    umap2 = scater::calculateUMAP(pca, transposed=TRUE, n_neighbors=k_umap)
        #})), collapse=" "))
        #colnames(umap) = c("UMAP1","UMAP2")
        #write_csv(as_tibble(umap, rownames = "cell_id"), file.path(output_dir, "umap2.csv.gz"))
        #rm(umap2);
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
    if (!file.exists(file.path(output_dir, "clusters_leiden.csv.gz"))) {
        snn = as(readRDS(file.path(output_dir, "snn_srt.rds"))$snn, "dgCMatrix")
        snn = igraph::graph_from_adjacency_matrix(snn, weighted=TRUE, mode="undirected")
        gc()
        # ~2h for 4M cells (30PCs, -1 iteratons, laptop)
        message("Leiden: ", paste(round(system.time({
            clusters = igraph::cluster_leiden(snn, objective_function='modularity', resolution_parameter = resolution, n_iterations = -1)
        })), collapse=" "))
        result = tibble(cell_id = clusters$name, label = clusters$membership)
        write_csv(result, file.path(output_dir, "clusters_leiden.csv.gz"))
        rm(clusters); rm(result); rm(snn); gc()
    }
    if (!file.exists(file.path(output_dir,"clusters_louvain.csv.gz"))) {
        snn = readRDS(file.path(output_dir, "snn_srt.rds"))$snn
        gc()
        # leiden interrupted after ~50min, up to 40GB (heavy swapping on laptop) for 9.1M cells
        # louvain ~4h
        message("Louvain: ", paste(round(system.time({
        clusters = FindClusters(snn, resolution = resolution)
        })), collapse = " "))
        result = tibble(cell_id = rownames(clusters), label = clusters[,1])
        write_csv(result, file.path(output_dir,"clusters_louvain.csv.gz"))
        rm(clusters); rm(result); rm(snn); gc()
    }
    if (!file.exists(file.path(output_dir, "clusters_sketch.csv.gz"))) {
        srt = CreateSeuratObject(counts = open_matrix_dir("data/bpcell_counts/"))
        srt = srt[,coi]
        srt = NormalizeData(srt, scale.factor = 10)
        # sketching: ~16min for 4M cells
        n_sketch = 200000
        system.time({ srt = SketchData(srt, ncells = n_sketch, method = "LeverageScore", sketched.assay = "sketch") })
        LayerData(srt, "scale.data") = as.matrix(LayerData(srt, "data")) # do *not* scale the data -> actually, should we scale?
        #srt = ScaleData(srt)
        srt = RunPCA(srt, features = rownames(srt), npcs = n_pca)
        # ~ 1 minute for 200K cells
        system.time({ srt = FindNeighbors(srt, k.param = k_cluster) })
        # ~ 2 minutes for 200K cells
        system.time({ srt = FindClusters(srt, resolution=resolution) })
        # ~ 30 minutes for 4M cells
        system.time({
            srt = ProjectData(srt, assay="RNA", full.reduction = "pca.full", sketched.assay = "sketch",
                              sketched.reduction = "pca", dims = 1:n_pca, refdata = list(cluster_full = "seurat_clusters"))
        })
        result = tibble(cell_id = colnames(srt), label = srt$cluster_full)
        write_csv(result, file.path(output_dir,"clusters_sketch.csv.gz"))
        rm(srt); rm(result); gc()        
    }
}

make_h3 = function(coi, output_dir, n_pca = 30, k_cluster = 15, k_umap = 15, resolution=1) {
    message("Running H3 level analysis for ", length(coi), " cells.")
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
        message("UMAP: ", paste(round(system.time({
            umap = scater::calculateUMAP(pca, transposed=TRUE, n_neighbors=k_umap)
        })), collapse=" "))
        colnames(umap) = c("UMAP1","UMAP2")
        write_csv(as_tibble(umap, rownames = "cell_id"), file.path(output_dir, "umap.csv.gz"))
        rm(umap); rm(pca); gc()
    }
    if (!file.exists(file.path(output_dir, "snn.rds"))) {
        pca = readRDS(file.path(output_dir, "pca.rds")); gc()
        # ~ 40min for 4M cells (10PCs, cluster 1 core)
        message("SNN: ", paste(round(system.time({
            snn = bluster::makeSNNGraph(pca, k=k_cluster, type="jaccard", BNPARAM = BiocNeighbors::AnnoyParam())
        })), collapse=" "))
        igraph::V(snn)$name = rownames(pca)
        saveRDS(snn, file.path(output_dir, "snn.rds"))
        rm(pca); rm(snn); gc()
    }
    if (!file.exists(file.path(output_dir, "clusters_leiden.csv.gz"))) {
        snn = readRDS(file.path(output_dir, "snn.rds"))
        gc()
        # ~10min for 4M cells (10PCs, 2 iterations, cluster 1 core)
        message("Leiden: ", paste(round(system.time({
            clusters = igraph::cluster_leiden(snn, objective_function='modularity', resolution_parameter = resolution, n_iterations = -1)
        })), collapse=" "))
        result = tibble(cell_id = clusters$name, label = clusters$membership)
        write_csv(result, file.path(output_dir, "clusters_leiden.csv.gz"))
        rm(clusters); rm(result); gc()
        message("Louvain: ", paste(round(system.time({
            clusters = igraph::cluster_louvain(snn, resolution = resolution)
        })), collapse=" "))
        result = tibble(cell_id = clusters$name, label = clusters$membership)
        write_csv(result, file.path(output_dir, "clusters_louvain.csv.gz"))
        rm(clusters); rm(result); rm(snn); gc()
    }
}

annotate_h1 = function() {
    sce = readRDS("data/sce.rds")
    sce$sample = sapply(strsplit(colnames(sce), ".", TRUE), "[", 1)
    sce$label = as.factor(read_csv("data/clusters_leiden.csv.gz")$label)
    reducedDim(sce, "UMAP") = as.matrix(column_to_rownames(read_csv("data/umap.csv.gz"), "cell_id"))
    sce = sce[, sample.int(ncol(sce))]
    gc()

    pdf("fig/h1_types.pdf")
    print(plot_raster_umap(sce, "label", "label", alpha=0.1))
    print(plot_raster_umap(sce, "Slc17a7", "label", alpha=0.1))
    print(plot_raster_umap(sce, "Gad1", "label", alpha=0.1))
    for (my_sample in sort(unique(sce$sample))) {
        sce$my_label = ifelse(sce$sample == my_sample, my_sample, "other")
        my_cols = setNames(c("red", gray(0.9, 0.1)), c(my_sample, "other"))
        print(plot_raster_umap(sce, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    }
    dev.off()
    # -> 4,5,6,7,9 * GLU
    # -> 2 GABA
}

annotate_h2 = function() {
    sce = readRDS("data/sce.rds")
    labels = read_csv("analysis/glu/clusters_leiden.csv.gz")
    sce = sce[, labels$cell_id]
    sce$label = as.factor(labels$label)
    sce$sample = sapply(strsplit(colnames(sce), ".", TRUE), "[", 1)
    reducedDim(sce, "UMAP") = as.matrix(column_to_rownames(read_csv("analysis/glu/umap.csv.gz"), "cell_id"))
    sce = sce[, sample.int(ncol(sce))]
    gc()
    
    cluster_size = table(sce$label)
    to_keep = names(cluster_size)[cluster_size > 10]
    sce = sce[, sce$label %in% to_keep,]
    gc()

    manual_labels = read_csv("data/annotation.csv.gz") %>%
        filter(startsWith(h1, "GLU")) %>%
        select(sample, h2) %>%
        deframe()
    sce$manual_label = manual_labels[colnames(sce)]
    sce$louvain_label = as.factor(deframe(read_csv("analysis/glu/clusters_louvain.csv.gz"))[colnames(sce)])
    sce$sketch_label = as.factor(deframe(read_csv("analysis/glu/clusters_sketch.csv.gz"))[colnames(sce)])
    
    # map leiden labels to manual labels
    cluster_annotation = as_tibble(colData(sce)[,c("label", "manual_label")]) %>%
        group_by(label, manual_label) %>%
        tally() %>%
        group_by(label) %>%
        mutate(f = n/sum(n)) %>%
        slice_max(f, n=1) %>%
        select(cluster_id = label, cluster_name = manual_label)
    write_csv(cluster_annotation, "analysis/glu/cluster_annotation.csv")

    pdf("fig/h2_types.pdf")
    print(plot_raster_umap(sce, "label", "label", alpha=0.01) + ggtitle("Leiden") + NoLegend())
    print(plot_raster_umap(sce, "manual_label", "manual_label", alpha=0.01)  + ggtitle("Manual annotation") + NoLegend())
    print(plot_raster_umap(sce, "louvain_label", "louvain_label", alpha=0.01) + ggtitle("Louvain") + NoLegend())
    print(plot_raster_umap(sce, "sketch_label", "sketch_label", alpha=0.01) + ggtitle("Sketch 200k, Louvain") + NoLegend())
    for (my_sample in sort(unique(sce$sample))) {
        sce$my_label = ifelse(sce$sample == my_sample, my_sample, "other")
        my_cols = setNames(c("red", gray(0.9, 0.1)), c(my_sample, "other"))
        print(plot_raster_umap(sce, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    }
    dev.off()
}

annotate_h3 = function() {
    sce = readRDS("data/sce.rds")
    labels = read_csv("analysis/l45it/clusters_leiden.csv.gz")
    sce = sce[, labels$cell_id]
    sce$label = as.factor(labels$label)
    sce$sample = sapply(strsplit(colnames(sce), ".", TRUE), "[", 1)
    reducedDim(sce, "UMAP") = as.matrix(column_to_rownames(read_csv("analysis/l45it/umap.csv.gz"), "cell_id"))
    sce = sce[, sample.int(ncol(sce))]
    gc()

    manual_labels = read_csv("data/annotation.csv.gz") %>%
        filter(sample %in% colnames(sce)) %>%
        select(sample, h3) %>%
        deframe()
    sce$manual_label = manual_labels[colnames(sce)]
    sce$louvain_label = as.factor(deframe(read_csv("analysis/l45it/clusters_louvain.csv.gz"))[colnames(sce)])

    pdf("fig/h3_types.pdf")
    print(plot_raster_umap(sce, "label", "label", alpha=0.5) + ggtitle("Leiden") + NoLegend())
    print(plot_raster_umap(sce, "manual_label", "manual_label", alpha=0.5)  + ggtitle("Manual annotation") + NoLegend())
    print(plot_raster_umap(sce, "louvain_label", "louvain_label", alpha=0.5) + ggtitle("Louvain") + NoLegend())
    for (my_sample in sort(unique(sce$sample))) {
        sce$my_label = ifelse(sce$sample == my_sample, my_sample, "other")
        my_cols = setNames(c("red", gray(0.9, 0.1)), c(my_sample, "other"))
        print(plot_raster_umap(sce, "my_label") + scale_color_manual(values = my_cols) + ggtitle(my_sample))
    }
    dev.off()
    
    labels = read_csv("analysis//l45it/clusters_louvain.csv.gz")
    labels$sample = sapply(strsplit(labels$cell_id, ".", TRUE), "[", 1)
    labels$condition = ifelse(grepl("1L", labels$sample) | grepl("3L", labels$sample), "ctl", "enu")
    labels %>%
        mutate(label = as.factor(label)) %>%
        group_by(sample, condition, label) %>%
        tally() %>%
        group_by(sample, condition) %>%
        mutate(f = n / sum(n)) %>%
        ggplot(aes(x=sample,y=f,fill=condition)) +
        geom_col(position = "dodge") +
        facet_wrap(~ label) +
        RotatedAxis()
    ggsave("fig/h3_prop.pdf")

    labels = read_csv("analysis/glu/clusters_louvain.csv.gz")
    labels$sample = sapply(strsplit(labels$cell_id, ".", TRUE), "[", 1)
    labels$condition = ifelse(grepl("1L", labels$sample) | grepl("3L", labels$sample), "ctl", "enu")
    labels %>%
        group_by(sample, condition, label) %>%
        tally() %>%
        group_by(sample, condition) %>%
        mutate(f = n / sum(n)) %>%
        ungroup() %>%
        filter(label<10) %>%
        ggplot(aes(x=sample,y=f,fill=condition)) +
        geom_col(position = "dodge") +
        facet_wrap(~ label) +
        RotatedAxis()
    ggsave("fig/h2_prop.pdf")
}


plot_raster_umap = function(barseq, color_by, text_by=NULL, point_size=0.1, text_size=3, raster.dpi=150, alpha=1) {
    to_plot = as_tibble(reducedDim(barseq, "UMAP"))
    colnames(to_plot) = c("UMAP_1", "UMAP_2")
    if (color_by %in% rownames(barseq)) {
        to_plot[[color_by]] = logcounts(barseq)[color_by,]
    } else {
        to_plot[[color_by]] = barseq[[color_by]]
    }
    # trick to make dot size as small as possible: https://stackoverflow.com/questions/34638902/point-size-in-ggplot-2-0-0
    result=to_plot %>%
        ggplot(aes(x=UMAP_1, y=UMAP_2, col=.data[[color_by]])) +
        ggrastr::geom_point_rast(size=point_size, shape=16, stroke=0, raster.dpi = raster.dpi, alpha=alpha) +
        theme_classic()
    if (is.numeric(to_plot[[color_by]])) {
        result = result + scale_color_viridis_c()
    }
    if (!is.null(text_by)) {
        to_plot[[text_by]] = barseq[[text_by]]
        to_plot = to_plot %>%
            group_by(.data[[text_by]]) %>%
            summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))
        result = result +
            geom_text(data=to_plot, aes(label=.data[[text_by]]), size=text_size, color="black")
    }
    return(result)
}

 
if (sys.nframe() == 0) {
    main()
}