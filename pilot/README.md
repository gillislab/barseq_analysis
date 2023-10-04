Analysis of BARseq data
===

Overview
---

This collection of R scripts was used to analyze a BARseq dataset of 40 slices covering a mouse cortex hemisphere. The analysis contains 3 major steps:
 - Preprocessing and clustering (scripts 01 to 09): low-quality neurons are filtered out, data are normalized and clustered, marker expression and localization is assessed to label cell types.
 - Comparison to a reference single-cell RNA-seq dataset (scripts 11 to 12): a nearest neighbor approach is used to map BARseq cell types to a reference scRNAseq dataset covering the hippocampus and isocortex (Yao et al., Cell, 2021).
 - Extraction of spatial patterns from the BARseq dataset (scripts 21 to 26): single-cell data are converted to pseudo-bulk format (one expression vector per cell type and spatial bin), expression is broken down into variability attributable to cell types vs pure spatial patterns, spatial patterns are extracted using NMF, NMF patterns are associated with genes and cell types.
 
The directory contains two helper scripts:
 - dataset.R: gives all the other scripts access to data (BARseq, Yao) and metadata (marker genes and ordering of cell types for plots).
 - gene_set_enrichment.R: hypergeometric test.
  
Step 1: preprocessing and clustering
---

Important! Scripts 03-05 of the pipeline are run iteratively to obtain cell types:
 - During the first iteration, preliminary clusters and marker expression are used to subdivide cells as "excitatory", "inhibitory" and "other" (H1 types).
 - During the second iteration, "excitatory" and "inhibitory" types are clustered separately (H2 types).
 - During the third iteration, each excitatory H2 type is clustered separately into even finer subtypes (H3 types).
The second iteration can only be performed after the first iteration (scripts 03-05) has been completed, the third iteration only after the second iteration has been completed. In total, scripts 03-05 will be run three times each (03->04->05 -> 03->04->05 -> 03->04->05). In the main() function, be careful to uncomment the lines corresponding to the current level of iteration and comment all other lines.
 
In order, run:
 - 01_plot_qc.R: generates QC plot and determine appropriate thresholds for filtering out low quality cells (based on the number of reads and number of detected genes).
 - 02_compute_marker_scores.R: computes the aggregate marker expression of previously defined marker sets corresponding to each major excitatory subclass (e.g., L2/3 IT, PT, L6b).
 - Iterative clustering (run 3 times, commenting and uncommenting appropriate lines, see above):
    - 03_analyze_barseq.R: clusters data using standard single-cell workflow (PCA -> kNN -> Louvain).
    - 04_plot_analysis.R: makes a series of graphs, including UMAP with cell types, UMAP with marker expression, localization of types along physical space. These visualization were designed to facilitate the assignement of cell types.
    - 05_annotate_clusters.R: enables the user to replace computer-generated cluster names with human-readable and biologically motivated names. We recommend re-running script 04 after this step to obtain UMAP visualizations with human-readable cell types.
 - 06_aggregate_labels.R: once the 3 rounds of clustering are done, collects all the labels and generates a single-cell table associating each single-cell with its H1, H2, and H3 labels.
 - 07_compute_marker_expression.R: for each H2 and H3 cell type, identifies and plots markers.
 - 08_plot_ccf_overlap.R: plots overlap of H2 and H3 cell types with standardized CCF regions.
 - 09_make_hierarchy.R: performs hierarchical clustering of H2 types.
 
 
 Step 2: mapping to reference data
 ---
 
 In order, run:
  - 11_assess_predictability.R: uses a surrogate scRNAseq dataset *instead of* BARseq data to evaluate how the restricted gene panel and the change of sensitivity from scRNAseq to BARseq affects our ability to assign cell types at the single-cell level.
  - 12_map_knn.R: maps BARseq cells to reference cell types (using the same mapping procedure as 11), then summarizes and plots result by evaluating the overlap between BARseq types and scRNAseq-predicted types.
  

Step 3: extraction of spatial patterns
---

In order, run:
 - 21_make_pseudobulks.R: aggregates the data at the pseudobulk level, i.e., sums the expression of single cells that belong to the same cell type and spatial bin, creating one gene expression vector for each cell type and spatial bin. Two levels of aggregation (H2 and H3 types).
 - 22_spatial_ct_contribution.R: performs ANOVA on the pseudobulk data to disentangle variability explained by types vs variability explained by spatial variability alone. Used to find genes that are mostly explained by variability across types (marker genes) and genes mostly explained by variability across spatial bins (spatial gradients/patterns).
 - 23_find_spatial_patterns.R: extracts spatial patterns using different strategies (PCA/NMF/ICA, in the end NMF was preferred). The pseudobulk matrices are reordered such that the spatial bins become features (genes and cell types become variables), then standard matrix factorizations are applied, generating a matrix of spatial patterns on one side, a loading matrix for each gene and cell type on the other side.
 - 24_plot_spatial_patterns.R: creates summary plots showing the extracted spatial patterns and variance explained by these patterns.
 - 25_pc_gene_association.R: identifies and plots genes that are strongly associated with each NMF pattern (genes that have similar spatial expression patterns as the NMF pattern).
 - 26_pc_ct_association.R: identifies and plots cell types that are strongly associated with each NMF pattern (cell types that have similar spatial distributions as the NMF pattern). 