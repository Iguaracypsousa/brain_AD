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
  library("ggplot2")
  library("patchwork")
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
# Initialize an empty list for storing aggregated results for all pathologies
aggregated_liana_all <- list()
# Process and save results for each subset
for (subset_name in names(liana_results)) {
  print(paste("Aggregating and filtering results for", subset_name))
  
  # Assuming liana_aggregate() and the filtering are correctly applied on your liana results
  # Aggregate results
  aggregated_liana <- liana_results[[subset_name]] %>% 
    liana_aggregate()
  
  # Filter results
  filtered_liana <- aggregated_liana %>%
    dplyr::filter(aggregate_rank <= 0.05)  # Ensure dplyr is specified to avoid method conflicts
  
  # Store the aggregated (and optionally filtered) results in the list
  aggregated_liana_all[[subset_name]] <- filtered_liana
  
  # Save the filtered results to a CSV file
  csv_filename <- paste0("EBI_course_2024/Cell_cell_results/LIANA_scRNA_Brain_", subset_name, ".csv")
  write.csv(filtered_liana, csv_filename, row.names = FALSE)
  
  print(paste("Results saved for", subset_name, "in", csv_filename))
}
#### Check some cell-cell interactions by pathology
major.celltype <- levels(scATAC_brain$major.celltype)
major.celltype

# Specify the cell type of interest
cell_type_of_interest <- "Mic"  # Replace with the cell type(s) of you're interest

##cell sending the signal (Ligand)
###Non AD
# Adjust the function call based on how 'liana_dotplot' is implemented and expects its arguments - sender
dotplot_Mic_source_nonAD <- liana_dotplot(
  aggregated_liana_all$nonAD,
  source_groups = c(cell_type_of_interest),  # Adjust this based on how your cell types are categorized (source or target)
  target_groups = NULL,  # Adjust target groups as needed or based on your specific analysis requirements
  ntop = 20  # Assuming you want the top 20 interactions
)
dotplot_Mic_source_nonAD
##Early AD
dotplot_Mic_source_earlyAD <- liana_dotplot(
  aggregated_liana_all$earlyAD,
  source_groups = c(cell_type_of_interest),  # Adjust this based on how your cell types are categorized (source or target)
  target_groups = NULL,  # Adjust target groups as needed or based on your specific analysis requirements
  ntop = 20  # Assuming you want the top 20 interactions
)
dotplot_Mic_source_earlyAD

##Late AD
dotplot_Mic_source_lateAD <- liana_dotplot(
  aggregated_liana_all$lateAD,
  source_groups = c(cell_type_of_interest),  # Adjust this based on how your cell types are categorized (source or target)
  target_groups = NULL,  # Adjust target groups as needed or based on your specific analysis requirements
  ntop = 20  # Assuming you want the top 20 interactions
)
dotplot_Mic_source_lateAD

###
jpeg("EBI_course_2024/figs/Microglia_source_cell_type_nonAD_EarlyAD_LateAD.jpeg",width = 32, height = 12,units = "in", res=200)


(dotplot_Mic_source_nonAD | dotplot_Mic_source_earlyAD | dotplot_Mic_source_lateAD) +
  plot_annotation(title = "LIANA Dot Plots: Non AD vs Early AD vs Late AD")

dev.off()

##cell receiving the signal (Receptor)
##
###Non AD
# Adjust the function call based on how 'liana_dotplot' is implemented and expects its arguments - sender
dotplot_Mic_target_nonAD <- liana_dotplot(
  aggregated_liana_all$nonAD,
  source_groups = NULL,  # Adjust this based on how your cell types are categorized (source or target)
  target_groups = c(cell_type_of_interest),  # Adjust target groups as needed or based on your specific analysis requirements
  ntop = 20  # Assuming you want the top 20 interactions
)
dotplot_Mic_target_nonAD
##Early AD
dotplot_Mic_target_earlyAD <- liana_dotplot(
  aggregated_liana_all$earlyAD,
  source_groups = NULL,  # Adjust this based on how your cell types are categorized (source or target)
  target_groups = c(cell_type_of_interest),  # Adjust target groups as needed or based on your specific analysis requirements
  ntop = 20  # Assuming you want the top 20 interactions
)
dotplot_Mic_target_earlyAD

##Late AD
dotplot_Mic_target_lateAD <- liana_dotplot(
  aggregated_liana_all$lateAD,
  source_groups = NULL,  # Adjust this based on how your cell types are categorized (source or target)
  target_groups = c(cell_type_of_interest),  # Adjust target groups as needed or based on your specific analysis requirements
  ntop = 20  # Assuming you want the top 20 interactions
)
dotplot_Mic_target_lateAD

###
jpeg("EBI_course_2024/figs/Microglia_target_cell_type_nonAD_EarlyAD_LateAD.jpeg",width = 32, height = 12,units = "in", res=200)


(dotplot_Mic_target_nonAD | dotplot_Mic_target_earlyAD | dotplot_Mic_target_lateAD) +
  plot_annotation(title = "LIANA Dot Plots: Non AD vs Early AD vs Late AD")

dev.off()
