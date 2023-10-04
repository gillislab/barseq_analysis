
library(tidyverse)
library(MetaMarkers)
source("dataset.R")


main = function() {
    # Load annotated BARseq (glutamatergic cells)
    barseq = load_barseq(normalization_factor = 10)
    logcounts(barseq) = log1p(cpm(barseq))/log(2)
    barseq$study_name = paste("barseq", barseq$slice, sep="_")
    clusters = read_csv("analysis/labels.csv")
    assertthat::are_equal(clusters$sample, colnames(barseq))
    barseq$cluster = clusters$cluster
    barseq$subclass = clusters$subclass
    barseq$class = clusters$class
    barseq = barseq[, barseq$class == "Glutamatergic" & barseq$subclass %in% ctx_subclass()]
    
    output_dir = "analysis/markers"
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    #compute_markers_(barseq, output_dir)
    plot_markers(barseq, output_dir)
}

compute_markers_ = function(barseq, output_dir) {
    subclass_markers = compute_markers(cpm(barseq), barseq$subclass)
    write_csv(subclass_markers, file.path(output_dir, "subclass_markers.csv"))

    cluster_markers = compute_markers(cpm(barseq)*1000, barseq$cluster)
    write_csv(cluster_markers, file.path(output_dir, "cluster_markers.csv"))

    cluster_markers = compute_markers(cpm(barseq)*1000, barseq$cluster, barseq$subclass)
    write_csv(cluster_markers, file.path(output_dir, "cluster_markers_by_subclass.csv"))
}

plot_markers = function(barseq, output_dir) {
    subclass_markers = read_csv(file.path(output_dir, "subclass_markers.csv"))
    subclass_markers %>%
        group_by(cell_type) %>%
        filter(rank(-auroc)<6)
    my_genes = c(
        "Slc30a3", # IT (T18)
        "Cux2", # L2/3-L4/5 IT (Y21b)
        "Rasgrf2", # L2/3 IT, panel Y21a
        "Rorb", # L4/5 IT (Y21a)
        "Etv1", # L5 IT, NP-PT-L5 IT (Y21b)
        "Scnn1a", # RSP -> Tshz2 + Slc30a3 + Scnn1a (Y21b), Coro6 (novel marker)
        "C1ql3", # L6 IT, panel Y21a
        "Fezf2", # non-IT marker, L5 IT (Y21a)
        "Rab3c", # PT -> panel Y21a; see also Hcn1 (panel Y21a)
        "Tle4", # CT-NP-L6b marker (Y21b)
        "Tshz2", # NP (T18, Y21b)
        "Foxp2", # CT (T18, Y21a, Y21b)
        "Ctgf", # L6b (Y21b); see also Cplx3 (Y21b)
        "Synpr" # Car3 -> panel Y21a
    )
    # single-cell, dot plot -> fraction expressing, min-max normalized expression
    marker_expr = as_tibble(as.matrix(t(logcounts(barseq)[my_genes,]))) %>%
        mutate(subclass = barseq$subclass, cell_type = barseq$cluster) %>%
        pivot_longer(c(-subclass, -cell_type), names_to = "gene", values_to = "expr") %>%
        group_by(gene, cell_type, subclass) %>%
        summarize(dr = mean(expr>0), expr = mean(expr), .groups="drop") %>%
        group_by(gene) %>%
        mutate(scaled_expr = (expr - min(expr)) / max(expr)) %>%
        ungroup()
    p1 = marker_expr %>%
        mutate(cell_type = factor(cell_type, cell_type_order())) %>%
        mutate(gene = factor(gene, rev(my_genes))) %>%
        ggplot(aes(x = cell_type, y = gene)) +
        geom_point(aes(size = dr, col = scaled_expr)) +
        theme_bw() +
        scale_radius(range = c(0,3)) +
        scale_color_gradient2(low = "white", mid = "gray80", high = "darkorange2", midpoint = 0.25, limits=c(0,1)) +
        theme(legend.position = "top",
              axis.text.x = element_text(angle=90, hjust=1,vjust=0.5),
              axis.text.y = element_text(size=12)) +
        labs(x=NULL, y=NULL)
    ggsave("figs/marker_expression.pdf", p1)
}

if (sys.nframe() == 0) {
    main()
}
