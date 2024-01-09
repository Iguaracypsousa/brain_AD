setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: h5ad file
#
# OUTPUT: QC and cell type annotation
library(Seurat)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(scater)
###Load data
scRNA_brain <- readH5AD("AD_Brain_analysis/427_ROSMAP_Data/rds_files/course_data_analysis/RNA.h5ad")

#######
gene_names <- mapIds(org.Hs.eg.db, keys=rownames(scRNA_brain), keytype="SYMBOL", columns="SYMBOL",column="SYMBOL")
###


##### get a few gene names
grep("^MT-",rowData(scRNA_brain)$SYMBOL,value = T)
grep("^RP[LS]",rowData(scRNA_brain)$SYMBOL,value = T)
grep("ATP8",rowData(scRNA_brain)$SYMBOL,value = T)
#######
columns(EnsDb.Hsapiens.v86)
ensdb_genes <- genes(EnsDb.Hsapiens.v86)
MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_name
is_mito <- rownames(scRNA_brain) %in% MT_names
table(is_mito)
#####Basic QC
assayNames(scRNA_brain)
assays(scRNA_brain)$counts <- assays(scRNA_brain)$X
#
scRNA_brain_cell <- perCellQCMetrics(scRNA_brain,subsets=list(Mito=is_mito))
scRNA_brain_feature <- perFeatureQCMetrics(scRNA_brain)
head(scRNA_brain_cell)

##
plotHighestExprs(scRNA_brain, exprs_values = "counts", 
                 feature_names_to_plot = "SYMBOL")
# Only if scRNA_brain_cell contains relevant data that should be in a SingleCellExperiment object
new_sce_object <- SingleCellExperiment(assays = list(counts = assay(scRNA_brain)), colData = scRNA_brain_cell)

scRNA_brain_cell <- addPerCellQC(scRNA_brain_cell, subsets=list(Mito=is_mito))
scRNA_brain_cell <- addPerFeatureQC(scRNA_brain_cell)
####
scATAC_brain <- readRDS("AD_Brain_analysis/427_ROSMAP_Data/rds_files/course_data_analysis/PeakMatrix.TSS6.cleaned.rds")




Heart_MI <- readRDS("AD_Brain_analysis/427_ROSMAP_Data/rds_files/course_data_analysis/RNA.h5ad")
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

