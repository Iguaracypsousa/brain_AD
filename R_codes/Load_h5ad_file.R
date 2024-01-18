reticulate::use_condaenv("/hps/software/users/petsalaki/users/iguaracy/miniconda3/envs/my_r_env/", required = TRUE)
setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: h5ad file
#
# OUTPUT: QC and cell type annotation
suppressPackageStartupMessages({
  library("zellkonverter")
  library("ggplot2")
  library("reticulate")
  library("SingleCellExperiment")
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library("EnsDb.Hsapiens.v86")
  library("scater")
  library("Matrix")
  library("Seurat")
  library("DoubletFinder")
  library("dplyr")
})

sc <- import("scanpy")
#Converting from Python to R
###Load data
scRNA_brain <- sc$read_h5ad("427_ROSMAP_Data/rds_files/course_data_analysis/RNA.h5ad")

##downsample by 5% by cell type
# Using reticulate to run Python code
cell_types <- scRNA_brain$obs$major.celltype

# Get unique cell types and their counts
unique_cell_types <- unique(cell_types)
cell_type_counts <- table(cell_types)

# Calculate the number of cells to sample for each cell type (5%)
sample_sizes <- round(cell_type_counts * 0.05)

set.seed(123)  # For reproducibility
sampled_indices <- unlist(lapply(unique_cell_types, function(cell_type) {
  indices <- which(cell_types == cell_type)
  sample(indices, sample_sizes[cell_type])
}))

py_run_string("import anndata")
scRNA_brain_sampled <- scRNA_brain[sampled_indices, ]
#########

#################### extracting counts

py_run_string("
import pandas as pd
# Example: Exporting a DataFrame to CSV
# Adjust the command according to the data you need

dense_matrix = scRNA_brain_sampled.X.toarray()

# Extract gene names (for row labels)
genes = scRNA_brain_sampled.var_names.to_list()

df = pd.DataFrame(dense_matrix.T, index=genes)
df.to_csv('427_ROSMAP_Data/rds_files/course_data_analysis/counts_downsample_RNA.csv', index=False)
")
counts <- read.csv('427_ROSMAP_Data/rds_files/course_data_analysis/counts_downsample_RNA.csv')
# Assuming 'your_matrix' is your matrix
column_sums <- colSums(counts)


colnames(counts) <- py$scRNA_brain_sampled$obs_names$to_list()
rownames(counts) <- py$scRNA_brain_sampled$var_names$to_list()


# Assuming 'your_matrix' is your matrix
column_sums <- colSums(counts)

# Extract the counts matrix
counts_matrix <- scRNA_brain$X
counts_matrix <- py_to_r(scRNA_brain$X)
# If you need to convert it to a dense matrix in R, you can do so using reticulate's py_to_r function
# (Be cautious with large datasets as this can be memory-intensive)
counts_matrix_dense <- reticulate::py_to_r(counts_matrix)
# Get the count data
counts <- assay(scRNA_brain, "X")  # replace "X" with the name of your assay if different
# Check if there are any non-zero values in the matrix
any(counts != 0)

# Extract cell-level metadata
cell_metadata <- as.data.frame(colData(scRNA_brain))


#####
# Create Seurat object using counts
scRNA_brain_seurat <- CreateSeuratObject(counts = counts)

# Add cell-level metadata
scRNA_brain_seurat@meta.data <- cell_metadata
# If you also want to add gene-level metadata
# (Replace 'features' with the correct slot if necessary)
rownames(scRNA_brain_seurat)
##checking metadata 
scRNA_brain_seurat@meta.data %>% 
  View()

###downsample to get 5% of each cell type in this dataset
# Calculate 5% of the cells for each cell type
celltype_counts <- table(scRNA_brain_seurat$major.celltype)
downsample_counts <- round(celltype_counts * 0.05)

# Perform the downsampling
downsampled_cells <- lapply(names(downsample_counts), function(cell_type) {
  cells_to_sample <- which(scRNA_brain_seurat$major.celltype == cell_type)
  sample(cells_to_sample, downsample_counts[cell_type])
})

# Combine all the cells
downsampled_cells <- unlist(downsampled_cells)

# Subset the Seurat object to keep only the downsampled cells
scRNA_brain_seurat_downsampled <- subset(scRNA_brain_seurat, cells = downsampled_cells)
scRNA_brain_seurat_downsampled
# downsampled_seurat now contains 5% of each cell type
##save data
saveRDS(scRNA_brain_seurat_downsampled, file = "EBI_course_2024/rds_files/scRNA_brain_seurat_downsampled.rds")

#######
###QC 
###QC 

scRNA_brain_seurat <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat_downsampled.rds")
DefaultAssay(scRNA_brain_seurat) <- "RNA"
#Let’s make violin plots of the selected metadata features. 
jpeg("EBI_course_2024/figs/scRNA_brain_seurat_VlnPlot.jpeg",width = 10, height = 7, units = 'in', res=300)
VlnPlot(scRNA_brain_seurat, features = c("n_genes","n_counts","pct_mito","pct_ribo"), ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
dev.off()

##Let’s plot some of the metadata features against each other and see how they correlate. 
##The number above each plot is a Pearson correlation coefficient.

jpeg("EBI_course_2024/figs/scRNA_brain_seurat_FeatureScatter_nCounts_RNA_pct_mito.jpeg",width = 6, height = 5, units = 'in', res=300)
FeatureScatter(scRNA_brain_seurat, feature1 = "n_counts", feature2 = "pct_mito")
dev.off()

jpeg("EBI_course_2024/figs/scRNA_brain_seurat_FeatureScatter_nCounts_RNA_nFeaturess_RNA.jpeg",width = 6, height = 5, units = 'in', res=300)
FeatureScatter(scRNA_brain_seurat, feature1 = "n_counts", feature2 = "n_genes")
dev.off()

jpeg("EBI_course_2024/figs/scRNA_brain_seurat_FeatureScatter_nCounts_RNA_pct_ribo.jpeg",width = 6, height = 5, units = 'in', res=300)
FeatureScatter(scRNA_brain_seurat, feature1 = "n_counts", feature2 = "pct_ribo")
dev.off()

jpeg("EBI_course_2024/figs/scRNA_brain_seurat_FeatureScatter_pct_ribo_pct_mito.jpeg",width = 6, height = 5, units = 'in', res=300)
FeatureScatter(scRNA_brain_seurat, feature1 = "pct_ribo", feature2 = "pct_mito")
dev.off()
#clean memory
rm(C)
gc()
##subset scRNA_brain_seurat
subset(
  scRNA_brain_seurat,
  n_genes>750 & 
    n_genes < 8000 & 
    pct_mito < 10 
  #percent.Largest.Gene < 10
) -> scRNA_brain_seurat

scRNA_brain_seurat
saveRDS(scRNA_brain_seurat, file = "EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")

###cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
NormalizeData(scRNA_brain_seurat, normalization.method = "LogNormalize", scale.factor = 10000) -> scRNA_brain_seurat
scRNA_brain_seurat <- CellCycleScoring(scRNA_brain_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
jpeg("EBI_course_2024/figs/scRNA_brain_seurat_VlnPlot_S.Score_G2M.Score.jpeg",width = 10, height = 7, units = 'in', res=300)
VlnPlot(scRNA_brain_seurat, features = c("S.Score", "G2M.Score"),
        ncol = 4, pt.size = 0.1)
dev.off()
#######
matrix_counts <- scRNA_brain_seurat@assays$RNA$counts
#Check if there are any NA values in the matrix
anyNA(matrix_counts)
# Attempt to convert the matrix to numeric and check for NAs
anyNA(as.numeric(matrix_counts))
##
# Check if there are any non-zero values in the matrix
any(matrix_counts != 0)
# Get cell type information
cell_types <- scRNA_brain_seurat@meta.data$major.celltype

# Sum the counts per cell type
counts_per_cell_type <- vapply(unique(cell_types), function(ct) {
  sum(matrix_counts[, cell_types == ct])
}, numeric(1))

# Create a named vector
names(counts_per_cell_type) <- unique(cell_types)

counts_per_cell_type

scRNA_brain_seurat <- FindVariableFeatures(scRNA_brain_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA_brain_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scRNA_brain_seurat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
jpeg("EBI_course_2024/figs/scRNA_brain_seurat_VariableFeaturePlot.jpeg",width = 6, height = 5, units = 'in', res=300)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
dev.off
gc()
#######important

scRNA_brain_seurat = ScaleData(scRNA_brain_seurat, vars.to.regress = c("n_genes", "pct_mito","S.Score", "G2M.Score"),
                     verbose = F)
#clean memory
#rm(data.filt)

###PCA

scRNA_brain_seurat <- RunPCA(scRNA_brain_seurat, features = VariableFeatures(object = scRNA_brain_seurat))
jpeg("EBI_course_2024/figs/scRNA_brain_seurat_ElbowPlot.jpeg",width = 8, height = 7, units = 'in', res=300)
ElbowPlot(scRNA_brain_seurat, ndims = 50)
dev.off()
scRNA_brain_seurat <- RunUMAP(scRNA_brain_seurat, dims = 1:15, verbose = F)
gc()
saveRDS(scRNA_brain_seurat, file = "EBI_course_2024/rds_files/scRNA_brain_seurat.SYMBOL_doublets.rds")


##checking metadata for doublets
scRNA_brain_seurat@meta.data %>% 
  View()
#scrublets were alredy removed before sharing the data. score cutoff value of 0.3 of was applied to remove doublets.
##run this part in the cluster
#Can run parameter optimization with paramSweep, but skip for now.

# sweep.res <- paramSweep_v3(data.filt) sweep.stats <-
# summarizeSweep(sweep.res, GT = FALSE) bcmvn <- find.pK(sweep.stats)
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cellscells.


#nExp <- round(ncol(scRNA_brain_seurat) * 0.04)  # expect 4% doublets
#scRNA_brain_seurat <- doubletFinder_v3(scRNA_brain_seurat, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
# name of the DF prediction can change, so extract the correct column name.

#load new rds file
scRNA_brain_seurat <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat.DB.rds")
DF.name = colnames(scRNA_brain_seurat@meta.data)[grepl("DF.classification", colnames(scRNA_brain_seurat@meta.data))]

jpeg("EBI_course_2024/figs/scRNA_brain_seurat_cowplot_doublets.jpeg",width = 12, height = 7, units = 'in', res=300)
cowplot::plot_grid(ncol = 2, DimPlot(scRNA_brain_seurat, group.by = "region") + NoAxes(),
                   DimPlot(scRNA_brain_seurat, group.by = DF.name) + NoAxes())
dev.off()

jpeg("EBI_course_2024/figs/scRNA_brain_seurat_VlnPlot_doublets.jpeg",width = 8, height = 7, units = 'in', res=300)
VlnPlot(scRNA_brain_seurat, features = "n_genes", group.by = DF.name, pt.size = 0.1)
dev.off()

dim(scRNA_brain_seurat)##14523 110780
scRNA_brain_seurat = scRNA_brain_seurat[, scRNA_brain_seurat@meta.data[, DF.name] == "Singlet"]
dim(scRNA_brain_seurat)##14523 106349






# ######for slides
# jpeg("EBI_course_2024/figs/scRNA_brain_seurat_UMAP_scRNA_brain_seurat_bf_reg.jpeg",width = 7, height = 5, units = 'in', res=300)
# DimPlot(scRNA_brain_seurat, reduction = "umap", group.by = "cell_source",pt.size = 0.01,raster=FALSE)
# dev.off()
# 
# jpeg("EBI_course_2024/figs/scRNA_brain_seurat_UMAP_scRNA_brain_seurat_af_reg.jpeg",width = 7, height = 5, units = 'in', res=300)
# DimPlot(scRNA_brain_seurat, reduction = "UMAP", group.by = "cell_source",pt.size = 0.01,raster=FALSE)
# dev.off()

##save data
saveRDS(scRNA_brain_seurat, file = "EBI_course_2024/EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")
