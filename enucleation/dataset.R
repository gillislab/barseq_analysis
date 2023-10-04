
library(tidyverse)
library(SingleCellExperiment)

BARSEQ_DIR = "~/data/brain/barseq_ctx"


sample_names = function() {
    c("D076_1L","D076_4L","D077_1L", "D077_2L","D078_1L","D078_2L","D079_3L","D079_4L")
}

load_barseq = function(slice = NULL) {
    sce = readRDS("data/sce.rds")
    if (!is.null(slice)) {
        sce = sce[, sce$slice %in% slice]
    }
    return(sce)
}

convert_gene_names = function(genes) {
    to_convert = c("Tafa1" = "Fam19a1", "Tafa2" = "Fam19a2", "Ccn2"="Ctgf")
    needs_conversion = genes %in% names(to_convert)
    genes[needs_conversion] = to_convert[genes[needs_conversion]]
    return(genes)
}

convert_to_cpm = function(M, total_counts = 1000000) {
    normalization_factor = Matrix::colSums(M) / total_counts
    if (is(M, "dgCMatrix")) {
        M@x = M@x / rep.int(normalization_factor, diff(M@p))
        return(M)
    } else {
        return(scale(M, center = FALSE, scale = normalization_factor))
    }
}

load_sample = function(sample_name, normalization_factor = 10, min_genes = 5, min_counts = 20) {
    sce = readRDS(file.path(BARSEQ_DIR, "230401", "v2", glue::glue("{sample_name}.rds")))
    rownames(sce) = convert_gene_names(rownames(sce))
    sce = sce[, colSums(counts(sce)) >= min_counts & colSums(counts(sce)>0) >= min_genes]
    colnames(sce) = paste0(sample_name, ".", sce$slice, "_", colnames(sce))
    cpm(sce) = convert_to_cpm(counts(sce), normalization_factor)
    return(sce)
}

load_sample_counts = function(sample_name, min_genes=5, min_counts=20) {
    sce = readRDS(file.path(BARSEQ_DIR, "230401/v1", glue::glue("{sample_name}.rds")))
    result = counts(sce)
    rownames(result) = convert_gene_names(rownames(result))
    colnames(result) = glue::glue("{sample_name}.{sce$slice}_{colnames(result)}")
    result = result[, colSums(result) >= min_counts & colSums(result>0) >= min_genes]
    return(result)
}

load_sample_annotation = function(sample_name, annotation_level) {
    metadata_dir = file.path(BARSEQ_DIR, "230401/v1", sample_name, "analysis", annotation_level)
    result = read_csv(file.path(metadata_dir, "cluster.csv"))%>%
        left_join(read_csv(file.path(metadata_dir, "cluster_annotation.csv")), by = c("label"="cluster_id")) %>%
        select(-label) %>%
        mutate(sample = glue::glue("{sample_name}.{sample}"))
    return(result)
}

load_labels = function() {
    read_csv("analysis/type_labels.csv.gz", show_col_types = FALSE)
}


load_ref = function(filename="barseq_210630.rds", slice = NULL, normalization_factor = 10, min_genes = 5, min_counts = 20, metadata_filename="barseq_210630_metadata.csv") {
    sce = readRDS(file.path(BARSEQ_DIR, filename))
    rownames(sce) = convert_gene_names(rownames(sce))
    sce = sce[, colSums(counts(sce)) >= min_counts & colSums(counts(sce)>0) >= min_genes]
    colnames(sce) = paste0(sce$slice, "_", colnames(sce))
    cpm(sce) = convert_to_cpm(counts(sce), normalization_factor)
    if (!is.null(slice)) {
        sce = sce[, sce$slice %in% slice]
    }
    metadata = read_csv(file.path(BARSEQ_DIR, metadata_filename)) %>%
        mutate(sample_name = paste0(slice, "_", id)) %>%
        select(-slice, -id)
    m = match(colnames(sce), metadata$sample_name)
    colData(sce) = cbind(colData(sce), metadata[m,])
    return(sce)
}

ctx_subclass = function() {
    c("L2/3 IT", "L4/5 IT", "L5 IT", "RSP/ACA UL", "RSP/ACA ML", "L6 IT", "L5 ET", "NP", "L6 CT", "L6b/CT", "L6b", "Car3")
}

ctx_subclass_ref = function() {
    c("L2/3 IT", "L4/5 IT", "L5 IT", "RSP UL", "RSP DL", "L6 IT", "PT", "NP", "CT", "L6b", "Car3")
}

cell_type_order_ref = function() {
    c("L2/3 IT UL-P/ENT", "L2/3 IT M-L/PIR", "L2/3 IT ML-A", "L2/3 IT M", "L2/3 IT ML-P", "L2/3 IT DL", "L2/3 IT Unclear/LQ",
      "L4/5 IT UL", "L4/5 IT ML-A", "L4/5 IT ML-P", "L4/5 IT M", "L4/5 IT M-L", "L4/5 IT P-L/LA", "L4/5 IT DL",
      "L5 IT UL", "L5 IT M", "L5 IT M-L/ENT DL/EPv", "L5 IT DL",
      "RSPd UL", "RSPv UL", "RSP L4", "RSP L5",
      "L6 IT UL", "L6 IT ML-A", "L6 IT L", "L6 IT ML-P", "L6 IT DL",
      "PT AUD", "PT P RSP/IT-like?",  "PT CTX Unclear",
      "PT CTX M/RSPd A", "PT CTX ML", "PT CTX MOp/s DL", "PT CTX A-M/L", "PT CTX P-L/RSPd-agl/TEa",
      "PT CTX UL", "PT CTX P", "PT RSPv/PPP",
      "NP CTX L5", "NP CTX L6/M/L", "NP CTX L/RSP/SUBv-sp", "NP RSP", "NP SUBd-sp", "NP PPP",
      "CT CTX UL", "CT CTX A", "CT CTX A/ACA", "CT RSP", "CT CTX P/RSP", "CT Unclear/LQ", "CT CTX A/L-V",
      "L6b CTX UL", "L6b CTX A/EPd?", "L6b CTX A L-V/EPd", "L6b CTX P",
      "Car3 CTX UL", "Car3 CLA/EPd", "Car3 EPd/CLA", "Car3 CTX P-L", "Car3 CTX DL")
}

ctx_roi = function() {
    read_delim("data/Major_cortical_areas.csv", show_col_types = FALSE)$Areas
}

read_panel = function(panel_name) {
    read_csv(file.path("panel", paste0(panel_name, ".csv"))) %>%
        select(cell_type, gene)
}

subclass_markers = function() {
    marker_panel = read_panel("final_panel") %>%
        select(cell_type, gene) %>%
        filter(!(cell_type %in% c("IT", "IT RSP", "IT TPE-ENT", "IT/PT", "PT"))) %>%
        mutate(group = "all")
    return(marker_panel)
}

read_clusters = function() {
    read_csv("analysis/cluster.csv")
}

sample_cols = function() {
    c(c("REF"=gray(0.8)), setNames(RColorBrewer::brewer.pal("Paired", n=8), sample_names()))
}

condition_cols = function() {
    c("CTL"="black","ENU"="#EA1E8C")
}