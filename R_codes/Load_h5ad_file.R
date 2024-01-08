#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: h5ad file
#
# OUTPUT: QC and cell type annotation
library(ggspavis)
library(ggplot2)
library(STutility)
library(EBImage) 
library(SpotClean)
library(S4Vectors)
library(SpatialExperiment)
library(STdeconvolve)
library(Matrix)
library(rhdf5)
library(zellkonverter)
library(scater)
library(patchwork)
library(scran)
###
Sys.setenv(VROOM_CONNECTION_SIZE = 1048576)  # Setting buffer size to 1 MB
min_cells <- 3  # Minimum number of cells a gene must be expressed in
min_features = 200  # Minimum number of features a cell must have

# Read the UMI counts
umi_counts_hvCOC <- readr::read_delim("rds_files/heart_organoides/GSE200277_hiPSC_UMIcounts.txt.gz", delim = "\t")
umi_counts_hvCOC <- column_to_rownames(umi_counts_hvCOC, var = "AAACCTGCACAGCCCA-1_hvCOC_01")

umi_counts_hvCOC_filtered <- umi_counts_hvCOC[rowSums(umi_counts_hvCOC > 0) >= min_cells, ]

cell_feature_counts_hvCOC = rowSums(umi_counts_hvCOC_filtered > 0)  # Count the number of features per cell
umi_counts_hvCOC_filtered = umi_counts_hvCOC_filtered[cell_feature_counts_hvCOC >= min_features, ]
#all(sapply(umi_counts_hvCOC_filtered, is.numeric))
non_numeric_columns <- sapply(umi_counts_hvCOC_filtered, function(column) !all(is.numeric(column)))
names(umi_counts_hvCOC_filtered)[non_numeric_columns]
umi_counts_hvCOC_filtered <- umi_counts_hvCOC_filtered[, !non_numeric_columns]
#all(sapply(umi_counts_hvCOC_filtered, is.numeric))


umi_counts_hESC <- readr::read_delim("rds_files/heart_organoides/GSE200277_hESC_UMIcounts.txt.gz", delim = "\t")
umi_counts_hESC <- column_to_rownames(umi_counts_hESC, var = "AAACGAACAAGGCTTT-1_CS_01")

umi_counts_hESC_filtered <- umi_counts_hESC[rowSums(umi_counts_hESC > 0) >= min_cells, ]

cell_feature_counts_hESC = rowSums(umi_counts_hESC_filtered > 0)  # Count the number of features per cell
umi_counts_hESC_filtered = umi_counts_hESC_filtered[cell_feature_counts_hESC >= min_features, ]
#all(sapply(umi_counts_hESC_filtered, is.numeric))
non_numeric_columns <- sapply(umi_counts_hESC_filtered, function(column) !all(is.numeric(column)))
names(umi_counts_hESC_filtered)[non_numeric_columns]
umi_counts_hESC_filtered <- umi_counts_hESC_filtered[, !non_numeric_columns]
#all(sapply(umi_counts_hESC_filtered, is.numeric))


org_hear_se_hvCOC <- CreateSeuratObject(counts = umi_counts_hvCOC_filtered)
saveRDS(org_hear_se_hvCOC, "rds_files/heart_organoides/seurat_hvCOC.rds")
org_hear_se_hESC <- CreateSeuratObject(counts = umi_counts_hESC_filtered)
saveRDS(org_hear_se_hESC, "rds_files/heart_organoides/seurat_hESC.rds")

