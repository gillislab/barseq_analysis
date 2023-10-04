library(tidyverse)
library(patchwork)
source("dataset.R")


main = function() {
    export_labels()
    plot_composition()
}

export_labels = function() {
    # h1
    result = read_csv("data/clusters_leiden.csv.gz")
    result$h1 = "Other"
    result$h1[result$label %in% c(5,6,7,9)] = "Glutamatergic"
    result$h1[result$label %in% c(2)] = "GABAergic"
    result = select(result, -label)
    
    # gather H3 labels
    h3 = read_csv("analysis/subglu/manual_annotation.csv")
    # change cluster id to reflect manual subclass assignment
    h3 = h3 %>%
        group_by(subclass) %>%
        mutate(h3 = paste0(subclass, "_", 1:n())) %>%
        select(cluster_id, h2 = subclass, h3, h3_name = cluster_name)
    # simplify name if there is a single member in subclass
    h2 = as.factor(h3$h2)
    is_singleton = h3$h2 %in% levels(h2)[tabulate(h2)==1]
    h3$h3[is_singleton] = h3$h2[is_singleton]
    # assign labels to individual cells
    h3_clusters = read_csv("analysis/subglu/cluster.csv.gz")
    names(h3_clusters) = c("cell_id", "cluster_id")
    h3_clusters = left_join(h3_clusters, h3)
    h3_clusters = select(h3_clusters, -cluster_id)
    
    # assemble h1 and h3 labels
    result = result %>% left_join(h3_clusters)
    head(result %>% filter(h1 == "Glutamatergic"))
    is_unclear = result$h1 == "Glutamatergic" & is.na(result$h3)
    result$h2[is_unclear] = "Unclear"
    result$h3[is_unclear] = "Unclear_0"
    result$h3_name[is_unclear] = "Unclear/unclustered"
    
    write_csv(result, "analysis/type_labels.csv")
}
 
plot_composition = function() {
    result = read_csv("analysis//type_labels.csv") %>%
        mutate(sample = substr(cell_id,1,7)) %>%
        mutate(condition = ifelse(substr(cell_id,6,6) %in% c("1","3"), "CTL", "ENU"))
    result$litter = substr(result$sample, 1, 4)
    p1 = result %>%
        ggplot(aes(x=h1, fill = sample)) +
        geom_bar(position = "dodge", show.legend = FALSE) +
        facet_grid(litter ~ .) +
        scale_fill_manual(values=sample_cols()) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=10)) +
        labs(x=NULL,y=NULL) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 2))
    p2 = result %>%
        filter(h1 == "Glutamatergic" & h2 != "Unclear") %>%
        group_by(sample, condition, litter, h2) %>%
        tally() %>%
        group_by(sample, condition, litter) %>%
        mutate(f = n / sum(n)) %>%
        ggplot(aes(x=h2, fill = sample)) +
        geom_col(aes(y=f), position = "dodge", show.legend = FALSE) +
        facet_grid(litter ~ .) +
        scale_fill_manual(values=sample_cols()) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=10)) +
        labs(x=NULL,y=NULL) +
        scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 2))
    p3 = result %>%
        filter(h2 %in% ctx_subclass()) %>%
        mutate(h2 = factor(h2, ctx_subclass())) %>%
        mutate(h3 = fct_reorder(h3_name, as.numeric(h2))) %>%
        group_by(sample, condition, litter, h3) %>%
        tally() %>%
        group_by(sample, condition, litter) %>%
        mutate(f = n / sum(n)) %>%
        ggplot(aes(x=h3, fill = sample)) +
        geom_col(aes(y=f), position = "dodge", show.legend = FALSE) +
        facet_grid(litter ~ .) +
        scale_fill_manual(values=sample_cols()) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=45,hjust=1,size=10)) +
        labs(x=NULL,y=NULL) +
        scale_y_continuous(labels = scales::percent, breaks = scales::pretty_breaks(n = 2))
    p = ((p1|p2)+plot_layout(widths = c(1,4)))/p3
    ggsave("fig/sample_composition.pdf", p, height = 25/3, width = 40/3)
}

if (sys.nframe() == 0) {
    main()
}
