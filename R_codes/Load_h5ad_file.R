reticulate::use_condaenv("/hps/software/users/petsalaki/users/iguaracy/miniconda3/envs/my_r_env/", required = TRUE)
setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: h5ad file
#
# OUTPUT: QC and cell type annotation
suppressPackageStartupMessages({
  library("reticulate")
  library("ggplot2")
  library("SingleCellExperiment")
  library("scater")
  library("Seurat")
  library("Matrix")
})
sc <- import("scanpy")
#Converting from Python to R
###Load data
scRNA_brain <- sc$read_h5ad("AD_Brain_analysis/427_ROSMAP_Data/rds_files/course_data_analysis/RNA.h5ad")
#
scRNA_brain
head(scRNA_brain$obs)
head(scRNA_brain$var)
scRNA_brain$X[1:5, 1:5]
#######
##
ggplot(scRNA_brain$obs, aes(x = n_counts, y = n_genes, colour = major.celltype)) +
  geom_point()
#Creating a Seurat object from AnnData
exprs <- t(scRNA_brain$X)
colnames(exprs) <- scRNA_Brain$obs_names$to_list()
rownames(exprs) <- scRNA_Brain$var_names$to_list()
# Create the Seurat object
seurat <- CreateSeuratObject(exprs)
# Set the expression assay
seurat <- SetAssayData(seurat, "data", exprs)
# Add observation metscRNA_Brain
seurat <- AddMetscRNA_Brain(seurat, scRNA_Brain$obs)
# Add fetaure metscRNA_Brain
seurat[["RNA"]][["n_cells"]] <- scRNA_Brain$var["n_cells"]
# Add embedding
embedding <- scRNA_Brain$obsm["X_umap"]
rownames(embedding) <- scRNA_Brain$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")



saveRDS(org_hear_se_hESC, "rds_files/heart_organoides/seurat_hESC.rds")

