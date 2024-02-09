setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: rds file
#
# OUTPUT: DEGs
suppressPackageStartupMessages({
  library("Seurat")
  library("tidyverse")
  library("magrittr")
  library("dplyr")
})
# ##load table

scRNA_brain <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")
Idents(scRNA_brain) <- scRNA_brain$major.celltype
scRNA_brain@meta.data$cell_type <- scRNA_brain$major.celltype
# Assuming your Seurat object is named 'scRNA_brain'
# Split the data by 'Pathology'
data_split <- SplitObject(scRNA_brain, split.by = "Pathology")

# Initialize a list to store results
markers_list <- list()

# Loop over the split data
for (i in names(data_split)) {
  print(paste("Processing", i))
  
  # Get cell types
  cell_types <- unique(data_split[[i]]@meta.data$cell_type)
  
  # Loop through each cell type
  for (cell_type in cell_types) {
    print(paste("Analyzing cell type", cell_type))
    
    # Check the number of cells in this cell type
    num_cells <- sum(data_split[[i]]@meta.data$cell_type == cell_type)
    
    if (num_cells < 3) {
      print(paste("Skipping cell type", cell_type, "- fewer than 3 cells"))
      next  # Skip to the next cell type
    }
    
    # Perform differential expression analysis
    markers <- FindMarkers(data_split[[i]], ident.1 = cell_type,only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
    
    # Store the results
    markers_list[[paste(i, cell_type, sep = "_")]] <- markers
    
    # Optional: Save results to files
    write.csv(markers, paste0("EBI_course_2024/Table_DEGs/Brain_DEGs_", i, "_", cell_type, ".csv"))
  }
}

# markers_list now contains all the differential expression results