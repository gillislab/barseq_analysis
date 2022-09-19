
library(tidyverse)
library(SingleCellExperiment)
library(MetaMarkers)
source("dataset.R")


main = function() {
    output_dir = "analysis"
    dir.create(output_dir, showWarnings = FALSE)
    
    marker_panel = subclass_markers()
    barseq = load_barseq(normalization_factor = 10)
    logcounts(barseq) = log1p(cpm(barseq))/log(2)

    marker_scores = score_cells(logcounts(barseq), marker_panel)
    rownames(marker_scores) = get_cell_type(rownames(marker_scores))
    marker_enrichment = compute_marker_enrichment(marker_scores)
    
    write_csv(as_tibble(t(marker_scores), rownames="sample_id"), file.path(output_dir, "marker_scores.csv"))
    write_csv(as_tibble(t(marker_enrichment), rownames="sample_id"), file.path(output_dir, "marker_enrichment.csv"))
}


if (sys.nframe() == 0) {
    main()
}
