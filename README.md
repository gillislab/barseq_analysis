Analysis of BARseq data
===


This repository contains collections of scripts that were used to analyze two BARseq experiments
 - pilot: Clustering, annotation and spatial analysis of a pilot experiment with 40 slices from a single mouse brain.
 - enucleation: Clustering and annotation of a dataset with 8 mouse brains (4 control, 4 enucleated).

All the scripts are R scripts designed and run on R version 4.2 on Ubuntu 22.04 (16 cores, 64GB of RAM) using the command line (`Rscript <name_of_script>.R`).
To fully reproduce the analyses, scripts take approximately 2 hours to complete for the pilot experiment, 10 hours for the enucleation experiment.

The scripts require the following packages, which can be installed from CRAN or Bioconductor:
tidyverse, hd5fr, rasterpdf, ggrastr, ggrepel, patckwork, Matrix, matrixStats, ica, NMF,
scater, scran, igraph, Seurat, SingleCellExperiment, MetaNeighbor, BiocNeighbors, BiocParallel.
They also require the MetaMarkers package, which can be installed from Github (`devtools::install_github("gillislab/MetaMarkers")`)
The installation process is usually fast, but can take up to 15 minutes starting from an empty R installation.

Please see READMEs inside each directory for more details about the analyses.
