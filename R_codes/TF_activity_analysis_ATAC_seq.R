setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq and scATACseq  data
# 
# INPUT: rds file
#
# OUTPUT: TF activity
suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("EnsDb.Hsapiens.v86")
  library("JASPAR2020")
  library("TFBSTools")
  library("patchwork")
})
###
# Load your scATAC-seq data
scATAC_brain <- readRDS("EBI_course_2024/rds_files/scRNA_ATAC_brain_ACTIVITY_norm.rds")
# Replace with your data loading code

# Identify peaks (accessible regions)
peaks <- scATAC_brain@assays$ATAC

# Annotate peaks with TF motifs
motif <- getMatrixSet(JASPAR2020, opts)
motif_enrichment <- FindMotifs(object = atac_data, features = peaks, motif = motif)

# Calculate TF activity
# This is a simplistic representation. In practice, you'll need to consider the context of your specific data and experiment.
tf_activity <- rowSums(motif_enrichment)