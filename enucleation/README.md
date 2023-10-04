Analysis of BARseq data
===

Overview
---

This collection of R scripts was used to analyze a BARseq dataset combining 8 mouse brains from an enucleation experiment (8 mice from 4 litters, 1 control and 1 enucleated mouse per litter). The analysis contains 2 major steps:
 - Preprocessing and clustering (scripts 01 to 07): low-quality neurons are filtered out, data are normalized and clustered, marker expression and localization is assessed to label cell types, quick comparison of annotations and QC metrics with pilot experiment.
 - Secondary analyses (scripts 11 to 15): differential analyses between control and enucleated samples, formal mapping of cell types with pilot experiment, assessment and visualization of within-H3 expression shifts.
 
The directory contains three helper scripts:
 - dataset.R: gives all the other scripts access to data (BARseq, Yao) and metadata (marker genes and ordering of cell types for plots).
 - aurocs.R: functions used to perform ROC tests.
 - visualization.R: functions to generate UMAP plots and flatmaps.
  
Step 1: preprocessing and clustering
---
 
In order, run:
 - 01_make_h1_h2.R: combine data from the 8 mice, apply Leiden clustering to the whole dataset (H1 level), then re-apply Leiden clustering to glutamatergic cells (H2 level).
 - 02_make_h3.R: apply Louvain clustering to each H2 type (H3 level).
 - 03_plot_analysis.R: makes a series of graphs, including UMAP with cell types, UMAP with marker expression, localization of types along physical space. These visualization were designed to facilitate the assignement of cell types.
 - 04_annotate.R: map cell types to pilot experiment and manual annotations based on localization to facilitate annotation.
 - 05_make_labels.R: gather all annotations (H1, H2, H3 level) into a single table and plot composition of each sample.
 - 06_plot_batches.R: visualize batch effects, gauge mixing of samples in UMAP space and across types.
 - 07_compare_ref.R: compare basic QC metrics (# genes, # reads, per-gene sensitivity) with pilot experiment.
 
 
 Step 2: secondary analyses
 ---
 
 In order, run:
  - 11_map_ct.R: use kNN predictor to compare annotations between enucleation experiment, pilot experiment, and published visual deprivation experiment (Cheng et al., 2022).
  - 12_differential_composition.R: apply ANOVA model to identify H3 x CCF combinations that are significantly enriched or depleted in enucleated samples.
  - 13_predict_condition.R: use kNN predictor to gauge local mixity of neurons from control and enucleated samples.
  - 14_plot_umap.R: visualize local (within-H3) composition shifts.