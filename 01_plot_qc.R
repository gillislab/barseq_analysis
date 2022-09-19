
library(tidyverse)
source("dataset.R")


main = function() {
    # current QC criteria
    min_genes = 5
    min_counts = 20

    output_dir = file.path("figs", "QC")
    dir.create(output_dir,recursive = TRUE, showWarnings = FALSE)
    barseq = load_barseq(min_genes = 0, min_counts = 0)
    qc = data.frame(
        sample = colnames(barseq),
        n_genes = colSums(counts(barseq) > 0),
        n_counts = colSums(counts(barseq))
    )
    length(unique(qc$sample)) # 2.26M cells
    mean(qc$n_counts==0) # 17% of cells with no counts
    sum(qc$n_counts>0) # 1.87M cells
    qc = filter(qc, n_counts>0)

    qc %>%
        group_by(n_genes) %>%
        tally(name = "n_cells") %>%
        ungroup() %>%
        ggplot(aes(x = n_genes, y = n_cells)) +
        geom_ribbon(aes(ymin=0,ymax=n_cells), fill = "lightblue2") +
        theme_bw(base_size = 20) +
        geom_vline(xintercept = min_genes, linetype="dashed") +
        labs(x="Number of genes", y="Number of cells")
    ggsave(file.path(output_dir, "n_genes.pdf"))
    
    qc %>%
        mutate(n_counts = ifelse(n_counts>200,200,n_counts)) %>%
        group_by(n_counts) %>%
        tally(name = "n_cells") %>%
        ungroup() %>%
        ggplot(aes(x = n_counts, y = n_cells)) +
        geom_ribbon(aes(ymin=0,ymax=n_cells), fill = "lightblue2") +
        theme_bw(base_size = 20) +
        geom_vline(xintercept = min_counts, linetype="dashed") +
        labs(x="Number of unique reads", y="Number of cells")
    ggsave(file.path(output_dir, "n_counts.pdf"))
    
    passes_qc = qc$n_genes >= min_genes & qc$n_counts >= min_counts
    mean(qc[passes_qc, "n_genes"]) # average of 27 genes / cell
    mean(qc[passes_qc, "n_counts"]) # average of 60 counts / cell
    f = round(100*mean(passes_qc))
    n = round(sum(passes_qc)/1e6,2)
    qc %>%
        ggplot(aes(x = n_genes, y = n_counts)) +
        ggrastr::geom_jitter_rast(alpha=0.01) +
        geom_density_2d() +
        theme_classic(base_size=20) +
        scale_x_log10() +
        scale_y_log10() +
        geom_hline(yintercept = min_counts, linetype = "dashed") +
        geom_vline(xintercept = min_genes, linetype = "dashed") +
        ggtitle(paste0(n, "M (", f, "%) cells pass QC criteria"))
    ggsave(file.path(output_dir, "counts_vs_genes.pdf"))
}

if (sys.nframe() == 0) {
    main()
}
