
library(tidyverse)
library(Matrix)
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
    
    expr = counts(barseq)
    labels = label_matrix(barseq$cluster)
    ct_expr = as.matrix(expr %*% labels)
    ct_dist = as.dist(1-cor(ct_expr, method="s"))
    my_hclust = hclust(ct_dist, method="average")
    ct_order = my_hclust$labels
    ct_weight = as.numeric(factor(ct_order, cell_type_order()))
    my_dend = as.dendrogram(my_hclust)
    my_dend2 = reorder(my_dend, ct_weight, agglo.FUN = "mean")
    pdf("figs/hierarchy_ctx.pdf")
    plot(my_dend2)
    dev.off()
}

label_matrix = function(labels) {
    labels = as.factor(labels)
    result = Matrix::sparseMatrix(i=1:length(labels),j=as.numeric(labels),
                                  dimnames = list(names(labels), levels(labels)))
    return(1*result)
}

if (sys.nframe() == 0) {
    main()
}
