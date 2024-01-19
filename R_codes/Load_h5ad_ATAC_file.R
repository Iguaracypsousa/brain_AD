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
  library("Signac")
  library("biovizBase")
})

sc <- import("scanpy")
#Converting from Python to R
###Load data
# Read the h5ad file
scRNA_ATAC_brain <- readRDS("427_ROSMAP_Data/rds_files/course_data_analysis/PeakMatrix.TSS6.cleaned.rds")

# Check cell type information
metadata <- as.data.frame(colData(scRNA_ATAC_brain))
table(colData(scRNA_ATAC_brain)$Celltype1)

##downsample by 5% by cell type
set.seed(123) # For reproducibility
downsampled_cells <- lapply(split(colnames(scRNA_ATAC_brain), colData(scRNA_ATAC_brain)$Celltype1), function(cells) {
  sample(cells, size = ceiling(length(cells) * 0.05))
})

# Flatten the list to get a vector of cell names
downsampled_cells <- unlist(downsampled_cells)

scRNA_ATAC_brain_downsampled_data <- scRNA_ATAC_brain[, downsampled_cells]
# Check the dimensions of the downsampled data
dim(scRNA_ATAC_brain_downsampled_data)

# Check the new distribution of cell types
table(colData(scRNA_ATAC_brain_downsampled_data)$Celltype4)
#########


###Convert RangedSummarizedExperiment to seurat object 


# Extract the counts matrix again
counts <- assay(scRNA_ATAC_brain_downsampled_data, "PeakMatrix")

# Create a new Seurat object
scRNA_ATAC_brain_seurat <- CreateSeuratObject(counts = counts)

# Add metadata (if needed)
metadata <- as.data.frame(colData(scRNA_ATAC_brain_downsampled_data))
rownames(metadata) <- colnames(scRNA_ATAC_brain_seurat)
scRNA_ATAC_brain_seurat@meta.data <- metadata

# Extract genomic ranges and create a Chromatin Assay
granges <- rowRanges(scRNA_ATAC_brain_downsampled_data)
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
genome(annotation) <- "hg38"
chrom_assay <- CreateChromatinAssay(counts = counts, ranges = granges, genome = 'hg38',annotation = annotation)
scRNA_ATAC_brain_seurat[['peaks']] <- chrom_assay

########
##save data
saveRDS(scRNA_ATAC_brain_seurat, file = "EBI_course_2024/rds_files/scRNA_ATAC_brain_downsampled_data.rds")

######QC amd processing

##
scRNA_ATAC_brain_seurat <- readRDS("EBI_course_2024/rds_files/scRNA_ATAC_brain_downsampled_data.rds")

# compute nucleosome signal score per cell
scRNA_ATAC_brain_seurat <- NucleosomeSignal(object = scRNA_ATAC_brain_seurat)
######Transcriptional start site (TSS) enrichment score.  Poor ATAC-seq experiments typically will have a low TSS enrichment score
jpeg("EBI_course_2024/figs/scRNA_ATAC_brain_TSSEnrichment.jpeg",width = 8, height = 7, units = 'in', res=300)
DensityScatter(scRNA_ATAC_brain_seurat, x = 'nCount_peaks', y = 'TSSEnrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
##distribution of each QC metric separately using a violin plot:
jpeg("EBI_course_2024/figs/scRNA_ATAC_brain_QC_metrics.jpeg",width = 8, height = 7, units = 'in', res=300)
VlnPlot(
  object = scRNA_ATAC_brain_seurat,
  features = c('nCount_peaks', 'TSSEnrichment', 'BlacklistRatio'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()
#Finally we remove cells that are outliers for these QC metrics.
scRNA_ATAC_brain_seurat <- subset(
  x = scRNA_ATAC_brain_seurat,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    BlacklistRatio < 0.05 &
    TSSEnrichment > 3
)
scRNA_ATAC_brain_seurat
######Normalization and linear dimensional reduction
scRNA_ATAC_brain_seurat <- FindTopFeatures(scRNA_ATAC_brain_seurat, min.cutoff = 5)
scRNA_ATAC_brain_seurat <- RunTFIDF(scRNA_ATAC_brain_seurat)
scRNA_ATAC_brain_seurat <- RunSVD(scRNA_ATAC_brain_seurat)

##
DepthCor(scRNA_ATAC_brain_seurat)

#Non-linear dimension reduction and clustering
scRNA_ATAC_brain_seurat <- RunUMAP(object = scRNA_ATAC_brain_seurat, reduction = 'lsi', dims = 2:30)
DimPlot(object = scRNA_ATAC_brain_seurat, label = TRUE, group.by = "Celltype1") + NoLegend()

##save data
saveRDS(scRNA_ATAC_brain_seurat, file = "EBI_course_2024/rds_files/scRNA_ATAC_brain_seurat_QC.rds")

