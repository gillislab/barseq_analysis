
library(tidyverse)
library(SingleCellExperiment)
source("~/archive/projects/common/brain_datasets.R", chdir=TRUE)
source("~/archive/data//brain/yao21//yao21.R")

BARSEQ_DIR = "~/archive/data/brain/barseq_ctx"


dataset_names = function() {
    c("barseq", "tasic18_v1", "tasic18_alm", biccn_datasets())
}

load_my_dataset = function(dataset_name) {
    if (startsWith(dataset_name, "tasic18")) {
        dataset = load_dataset(dataset_name)
        dataset$subclass = get_tasic_subclass(dataset)
    } else if (dataset_name == "barseq") {
        dataset = load_barseq("barseq_210424")
        dataset$subclass = "unknown"
    } else {
        dataset = load_biccn_dataset(dataset_name)
    }
    return(dataset)
}

get_tasic_subclass = function(dataset) {
    result = dataset$subclass
    result[dataset$cluster == "L6 IT VISp Car3"] = "L6 IT Car3"
    result[dataset$subclass == "NP"] = "L5/6 NP"
    result[dataset$subclass == "L5 PT"] = "L5 ET"
    return(result)
}

load_barseq = function(filename="barseq_210630.rds", slice = NULL, normalization_factor = 10, min_genes = 5, min_counts = 20, metadata_filename="barseq_210630_metadata.csv") {
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

ctx_subclass = function() {
    c("L2/3 IT", "L4/5 IT", "L5 IT", "RSP UL", "RSP DL", "L6 IT", "PT", "NP", "CT", "L6b", "Car3")
}

cell_type_order = function() {
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

previous_panel = function() {
    my_genes = unique(read_panel("test_panel")$gene)
    my_genes = c(my_genes, "Slc17a7", "Gad1", "Slc30a3")
    return(my_genes)
}

read_clusters = function() {
    read_csv("analysis/cluster.csv")
}

load_full_yao = function() {
    yao_filename = "data/full_yao.rds"
    if (file.exists(yao_filename)) {
        result = readRDS(yao_filename)
    } else {
        result = make_full_yao()
        saveRDS(result, yao_filename)
    }
    return(result)
}

make_full_yao = function() {
    barseq = load_barseq()
    gene_set = rownames(barseq)
    gene_set = gene_set[!startsWith(gene_set, "unused")]
    
    all_datasets = data.frame(technology = yao21_technologies()) %>%
        group_by(technology) %>%
        summarize(brain_region = yao21_brain_regions(technology)) %>%
        transpose()  
    
    result = lapply(set_names(all_datasets), function(d) {
        print(paste(d$technology, d$brain_region))
        dataset = load_yao21(d$technology, d$brain_region)
        dataset = dataset[gene_set,]
        dataset$study_name=dataset$study_id
        return(dataset)
    })
    result = MetaNeighbor::mergeSCE(result)
    return(result)
}