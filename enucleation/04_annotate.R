
library(tidyverse)
source("dataset.R")


main = function() {
    ## mapping to original dataset with MetaNeighbor
    recompute_mn = FALSE
    if (recompute_mn) {
        run_mn()
    }
    
    ## h1 level
    compare_h1()
    
    ## finer levels
    all_hits = read_csv("data/mn_cluster.csv.gz", show_col_type = FALSE) %>%
        pivot_longer(-cluster_id, names_to = "ct2", values_to = "auroc") %>%
        separate(cluster_id, c("study_1", "ct_1"), "\\|") %>%
        separate(ct2, c("study_2", "ct_2"), "\\|")
    all_hits = filter(all_hits, (study_1 != study_2) & (study_1 != "reference"))
    all_hits = mutate(all_hits, hit_type = ifelse(study_2 == "reference", "ref_hit", "internal_hit"))
    top_hits = all_hits %>%
        group_by(ct_1, ct_2, hit_type) %>%
        summarize(med_auroc = median(auroc), iqr_auroc = iqr(auroc)) %>%
        group_by(ct_1, hit_type) %>%
        slice_max(med_auroc, n = 3) %>%
        group_by(ct_1) %>%
        mutate(hit_rank = n() - rank(med_auroc)) %>%
        ungroup()
    detailed_top_hits = inner_join(all_hits, top_hits) %>%
        mutate(subclass = str_extract(ct_1, "(.+)_\\d+$", group=1))
    pdf("fig/mn_cluster_matching.pdf")
    for (my_subclass in sort(unique(detailed_top_hits$subclass))) {
        p = detailed_top_hits %>%
            filter(subclass == my_subclass) %>%
            mutate(ct_2 = fct_reorder(paste0(ct_2, "|", ct_1), hit_rank, last)) %>%
            ggplot(aes(x=ct_2,y=auroc, fill = hit_type)) +
            geom_boxplot() +
            facet_wrap(~ ct_1, scales = "free") + 
            scale_x_discrete(labels = function(x) { gsub("\\|.+$", "", x) }) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust=1))
        print(p)
    }
    dev.off()
    best_hit = top_hits %>%
        group_by(ct_1, hit_type) %>%
        mutate(hit_rank = n() - rank(med_auroc)) %>%
        mutate(second_hit_auroc = med_auroc[hit_rank==1]) %>%
        ungroup() %>%
        filter(hit_rank == 0 & hit_type == "ref_hit") %>%
        mutate(subclass = str_extract(ct_1, "(.+)_\\d+$", group=1)) %>%
        select(cluster_id = ct_1, cluster_name = ct_2, subclass) %>%
        mutate(notes = "")
    write_csv(best_hit, "analysis/subglu/mn_hits.csv")
        
    ## also export manual labels
    input_dir = "analysis/subglu/"
    clusters = read_csv(file.path(input_dir, "cluster.csv.gz"))
    manual_labels = read_csv("data/annotation.csv.gz")
    manual_annotation = inner_join(clusters, manual_labels)
    top_manual_annotation = manual_annotation %>%
        group_by(label, h3) %>%
        tally() %>%
        group_by(label) %>%
        mutate(f = n / sum(n)) %>%
        slice_max(f, n=3)
    write_csv(top_manual_annotation, "analysis/subglu/manual_hits.csv")
}

run_mn = function() {
    sce_list = load_all_data()
    mn_model = train_model(sce_list)
    auroc = run_mn_(sce_list, mn_model, one_vs_best=FALSE)
    write_csv(as_tibble(auroc, rownames = "cluster_id"), "data/mn_cluster.csv.gz")
    auroc = run_mn_(sce_list, mn_model, one_vs_best=TRUE)
    write_csv(as_tibble(auroc, rownames = "cluster_id"), "data/mn_cluster_1vbest.csv.gz")
}

train_model = function(sce_list) {
    my_hvg = rownames(sce_list[[1]])
    my_hvg = my_hvg[!(grepl("unused", my_hvg))]
    cluster_model = lapply(sce_list, function(sce) {
        MetaNeighbor::trainModel(my_hvg, sce, study_id = sce$sample, cell_type = sce$cluster_id)
    })
    cluster_model = do.call(cbind, cluster_model)
    gc()
    return(cluster_model)
}

run_mn_ = function(sce_list, mn_model, one_vs_best) {
    auroc = lapply(sce_list, function(sce) {
        result = MetaNeighbor::MetaNeighborUS(trained_model=mn_model, dat=sce, study_id = sce$sample, cell_type = sce$cluster_id, one_vs_best = one_vs_best)
        gc()
        return(result)
    })
    auroc = do.call(rbind, auroc)
    auroc = auroc[colnames(auroc),]
    gc()
    return(auroc)
}

load_all_data = function() {
    result = load_current()
    result$ref = load_reference()
    return(result)
}

load_current = function(input_dir = "analysis/subglu") {
    barseq = load_barseq()
    counts(barseq) = logcounts(barseq)
    sample_cell = strsplit(colnames(barseq), split = ".", fixed = TRUE) # replace by regex for efficiency
    barseq$sample = sapply(sample_cell, "[", 1)
    barseq$slice = as.factor(sapply(strsplit(sapply(sample_cell, "[", 2), split = "_", fixed = TRUE),"[", 1))

    clusters = read_csv(file.path(input_dir, "cluster.csv.gz"))
    cluster_annotation = deframe(read_csv(file.path(input_dir, "cluster_annotation.csv")))
    # keep only clusters not annotated as "NA"
    is_not_na = names(cluster_annotation)[!is.na(cluster_annotation)]
    clusters = filter(clusters, label %in% is_not_na)
    barseq = barseq[, as.character(clusters$sample)]
    barseq$cluster_id = as.factor(clusters$label)
    colLabels(barseq) = cluster_annotation[barseq$cluster_id]
    result = lapply(sample_names(), function(s) { barseq[,barseq$sample==s] })
    return(result)
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

order_rows_according_to_cols = function(M, alpha = 100) {
    M <- M**alpha
    row_score <- colSums(t(M)*seq_len(ncol(M)), na.rm=TRUE)/rowSums(M, na.rm=TRUE)
    return(order(row_score))
}

compare_h1 = function() {
    if (!file.exists("data/mn_h1.csv.gz")) {
        # current data
        current = load_barseq()
        current$sample = "current"
        clusters = read_csv("data/clusters_leiden.csv.gz")
        current$cluster_id = as.factor(clusters$label)
        colLabels(current) = current$cluster_id
        # reference data
        barseq = load_ref()
        barseq$sample = "reference"
        clusters = read_csv("data/ref_labels.csv")
        barseq$cluster_id = as.factor(clusters$class)
        colLabels(barseq) = barseq$cluster_id
        # run MN
        sce_list = list(current = current, ref = barseq)
        mn_model = train_model(sce_list)
        auroc = run_mn_(sce_list, mn_model, one_vs_best=FALSE)    
        write_csv(as_tibble(auroc, rownames = "cluster_id"), "data/mn_h1.csv.gz")
    }
    auroc = as.matrix(column_to_rownames(read_csv("data/mn_h1.csv.gz"), "cluster_id"))
    pdf("fig/mn_h1_matching.pdf")
    MetaNeighbor::plotHeatmap(auroc)
    dev.off()
}

if (sys.nframe() == 0) {
    main()
}