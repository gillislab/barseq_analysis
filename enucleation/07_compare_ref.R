
library(tidyverse)
source("dataset.R")


main = function() {
    current = load_current()
    current$n_genes = colSums(logcounts(current)>0)
    ref = load_reference()
    ref$n_genes = colSums(counts(ref)>0)
    ref$sample = "REF"
    ref$condition = "REF"
    logcounts(ref) = log1p(cpm(ref))
    all = MetaNeighbor::mergeSCE(sce_list = list(current = current, ref = ref))
    rm(ref); rm(current); gc()
    all$n_reads = colSums(counts(all))

    as_tibble(colData(all)) %>%
        ggplot(aes(x=n_genes)) +
        geom_histogram(binwidth = 1) +
        facet_grid(sample ~ .) +
        theme_classic() +
        geom_vline(xintercept = 28, linetype="dashed") +
        geom_vline(xintercept = 37, linetype="dashed") +
        labs(x = "# genes detected per cell")
    ggsave("fig/ref_n_genes.pdf")

    # median # genes 28 vs 37
    median(all$n_genes[all$sample=="REF"]) # 28
    median(all$n_genes[all$sample!="REF"]) # 37
    as_tibble(colData(all)) %>%
        mutate(sample = factor(sample, names(sample_cols()))) %>%
        ggplot(aes(x=sample, y=n_genes, fill = sample)) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_manual(values = sample_cols()) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) +
        labs(y = "# genes detected per cell", x = NULL)
    ggsave("fig/ref_n_genes_bp.pdf")

    # median # reads
    median(all$n_reads[all$sample=="REF"]) # 57
    median(all$n_reads[all$sample!="REF"]) # 88
    as_tibble(colData(all)) %>%
        mutate(sample = factor(sample, names(sample_cols()))) %>%
        ggplot(aes(x=sample, y=n_reads, fill = sample)) +
        geom_boxplot(show.legend = FALSE) +
        scale_fill_manual(values = sample_cols()) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45, hjust = 1)) +
        labs(y = "# reads per cell", x = NULL) +
        scale_y_log10()
    ggsave("fig/ref_n_reads_bp.pdf")

    dr = lapply(set_names(unique(all$sample)), function(my_sample) {
        enframe(rowMeans(logcounts(all)[, all$sample == my_sample] > 0), "gene", "dr")
    })
    dr = bind_rows(dr, .id = "sample")
    dr %>%
        mutate(sample = ifelse(sample == "original", "original", "current")) %>%
        pivot_wider(id_cols = "gene", names_from = "sample", values_from = "dr", values_fn = mean) %>%
        ggplot(aes(x=original, y=current, label=gene)) +
        geom_text() +
        geom_abline(slope = 1, linetype="dashed") +
        theme_classic() +
        labs(x = "Detection rate in original dataset", y="Detection rate in current dataset")
    ggsave("fig/ref_dr.pdf")
    dr %>%
        mutate(sample = ifelse(sample == "original", "original", "current")) %>%
        pivot_wider(id_cols = "gene", names_from = "sample", values_from = "dr", values_fn = mean) %>%
        mutate(dr_fc = current / original, dr_diff = current - original) %>%
        ggplot(aes(x=original, y=dr_diff, label=gene)) +
        geom_text() +
        geom_hline(yintercept = 0, linetype="dashed") +
        theme_classic() +
        labs(x = "Detection rate in original dataset", y="Detection rate difference (current - original)")
    ggsave("fig/ref_dr_diff.pdf")
    dr %>%
        mutate(sample = ifelse(sample == "original", "original", "current")) %>%
        pivot_wider(id_cols = "gene", names_from = "sample", values_from = "dr", values_fn = mean) %>%
        mutate(dr_fc = current / original, dr_diff = current - original) %>%
        filter(dr_diff > 0.1) %>%
        ggplot(aes(x=original, y=dr_diff, label=gene)) +
        ggrepel::geom_text_repel() +
        theme_classic() +
        labs(x = "Detection rate in original dataset", y="Detection rate difference (current - original)")
    ggsave("fig/ref_dr_diff_zoom.pdf")
    dr %>%
        mutate(sample = ifelse(sample == "original", "original", "current")) %>%
        pivot_wider(id_cols = "gene", names_from = "sample", values_from = "dr", values_fn = mean) %>%
        mutate(dr_fc = current / original, dr_diff = current - original) %>%
        ggplot(aes(x=original, y=dr_fc, label=gene)) +
        geom_text() +
        geom_hline(yintercept = 1, linetype="dashed") +
        theme_classic() +
        scale_y_log10() +
        labs(x = "Detection rate in original dataset", y="Detection rate FC (current / original)")
    ggsave("fig/ref_dr_fc.pdf")
}

load_current = function(label_file = "analysis/type_labels.csv.gz") {
    barseq = load_barseq()
    sample_cell = strsplit(colnames(barseq), split = ".", fixed = TRUE) # replace by regex for efficiency
    barseq$sample = sapply(sample_cell, "[", 1)
    barseq$slice = as.factor(sapply(strsplit(sapply(sample_cell, "[", 2), split = "_", fixed = TRUE),"[", 1))
    clusters = read_csv(label_file, show_col_types = FALSE)
    barseq$h1 = clusters$h1
    barseq$h2 = clusters$h2
    barseq$h3 = clusters$h3
    return(barseq)
}

load_reference = function(label_file = "data/ref_labels.csv") {
    barseq = load_ref()
    barseq$sample = "reference"
    clusters = read_csv(label_file)
    clusters = filter(clusters, class == "Glutamatergic")
    barseq = barseq[, as.character(clusters$sample)]
    barseq$subclass = as.factor(clusters$subclass)
    barseq$cluster_id = as.factor(clusters$cluster)
    colLabels(barseq) = barseq$cluster_id
    return(barseq)
}

if (sys.nframe() == 0) {
    main()
}
