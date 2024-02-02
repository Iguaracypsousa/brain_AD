setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: rds file
#
# OUTPUT: Cell-cell communication results
suppressPackageStartupMessages({
  library("Seurat")
  library("tidyverse")
  library("magrittr")
  library("liana")
  library("dplyr")
})


####load data
scATAC_brain <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")

Idents(scATAC_brain) <- scATAC_brain$major.celltype
# Split the data by 'Pathology'
data_split <- SplitObject(scATAC_brain, split.by = "Pathology")

# Initialize a list to store LIANA results
liana_results <- list()

# Run LIANA for each subset
for (subset_name in names(data_split)) {
  print(paste("Processing", subset_name))
  
  # Assuming that the LIANA wrap function works directly on a Seurat object
  liana_subset_result <- liana_wrap(data_split[[subset_name]])
  
  # Store the result
  liana_results[[subset_name]] <- liana_subset_result
  
}

# Now, liana_results contains all the LIANA results for each subset

# Process and save results for each subset
for (subset_name in names(liana_results)) {
  print(paste("Aggregating and filtering results for", subset_name))
  
  # Aggregate results
  aggregated_liana <- liana_results[[subset_name]] %>%
    liana_aggregate()
  
  # Filter results
  filtered_liana <- aggregated_liana %>%
    filter(aggregate_rank <= 0.05)  # Adjust this threshold as needed
  
  # Save the filtered results to a CSV file
  csv_filename <- paste0("EBI_course_2024/Cell_cell_results/LIANA_scRNA_Brain_", subset_name, ".csv")
  write.csv(filtered_liana, csv_filename, row.names = FALSE)
  
  print(paste("Results saved for", subset_name, "in", csv_filename))
}

