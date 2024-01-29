setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scATACseq data
# 
# INPUT: h5ad file
#
# OUTPUT: QC and cell type annotation
suppressPackageStartupMessages({
  library("Seurat")
  library("Matrix")
  library("EnsDb.Hsapiens.v86")
  library("AnnotationDbi")
  library("Signac")
  library("dplyr")
  library("SummarizedExperiment")
})

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
scRNA_ATAC_brain_seurat[['ATAC']] <- chrom_assay

########
##save data
saveRDS(scRNA_ATAC_brain_seurat, file = "EBI_course_2024/rds_files/scRNA_ATAC_brain_downsampled_data.rds")


######QC and processing
scRNA_ATAC_brain_seurat <- readRDS("EBI_course_2024/rds_files/scRNA_ATAC_brain_downsampled_data.rds")

##
scRNA_ATAC_brain_seurat <- readRDS("EBI_course_2024/rds_files/scRNA_ATAC_brain_downsampled_data.rds")
DefaultAssay(scRNA_ATAC_brain_seurat) <- 'ATAC'
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
    nCount_peaks < 40000 &
    BlacklistRatio < 0.05 &
    TSSEnrichment > 3
)
scRNA_ATAC_brain_seurat
######Normalization and linear dimensional reduction
scRNA_ATAC_brain_seurat <- FindTopFeatures(scRNA_ATAC_brain_seurat, min.cutoff = 5)
scRNA_ATAC_brain_seurat <- RunTFIDF(scRNA_ATAC_brain_seurat)
scRNA_ATAC_brain_seurat <- RunSVD(scRNA_ATAC_brain_seurat)


#Non-linear dimension reduction and clustering
scRNA_ATAC_brain_seurat <- RunUMAP(object = scRNA_ATAC_brain_seurat, reduction = 'lsi', dims = 2:30)
DimPlot(object = scRNA_ATAC_brain_seurat, label = TRUE, group.by = "Celltype1") + NoLegend()


#Add the Gene score to the ATAC-seq object 
metadata_scRNA_ATAC_brain_seurat <- scRNA_ATAC_brain_seurat@meta.data
scRNA_ATAC_brain_geneScore <- readRDS("427_ROSMAP_Data/rds_files/course_data_analysis/GeneScoreMatrix.TSS6.cleaned.rds")
##get Gene score 
metadata_geneScore <- as.data.frame(colData(scRNA_ATAC_brain_geneScore))
# Extract gene names or identifiers from rowData
gene_names <- rowData(scRNA_ATAC_brain_geneScore)$name


# Check the first few names to ensure they are what you expect
head(gene_names)


# Find common row names with scRNA_ATAC_brain_seurat
DefaultAssay(scRNA_ATAC_brain_seurat) <- 'RNA'
common_rows <- intersect(rownames(metadata_scRNA_ATAC_brain_seurat), rownames(metadata_geneScore))

# Subset both data frames to keep only rows with common row names
metadata_geneScore <- metadata_geneScore[common_rows, ]
##select fragments
metadata_geneScore

# Extract the counts matrix again
GeneScoreMatrix <- assay(scRNA_ATAC_brain_geneScore, "GeneScoreMatrix")
# Ensure the length of gene_names matches the number of rows in filtered_matrix
if (length(gene_names) == nrow(GeneScoreMatrix)) {
  rownames(GeneScoreMatrix) <- gene_names
} else {
  stop("The number of gene names does not match the number of rows in the filtered matrix.")
}
# Assuming `metadata_geneScore` holds relevant identifiers in its rows or a specific column
# If the metadata contains cell identifiers and you want to filter cells:
common_cells <- intersect(colnames(GeneScoreMatrix), rownames(metadata_geneScore))  # Adjust as necessary
# Filter columns
filtered_matrix <- GeneScoreMatrix[, common_cells]
rownames(filtered_matrix)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
scRNA_ATAC_brain_seurat[['ACTIVITY']] <- CreateAssayObject(counts = filtered_matrix)

##save data
saveRDS(scRNA_ATAC_brain_seurat, file = "EBI_course_2024/rds_files/scRNA_ATAC_brain_seurat_QC.rds")
