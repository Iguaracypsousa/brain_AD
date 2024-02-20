setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: CCC results, TF activity and DEGs
#
# OUTPUT: Seed nodes for scphuEGO - receptors, TFs and DEGs
suppressPackageStartupMessages({
  library("tidyverse")
  library("magrittr")
  library("dplyr")
  library("enrichR")
})
#######
##Here we will perform enrichment analysis on the receptors and TFs for each cell type by pathology
##Then we will only keep the DEGs downstream to enriched pathways for receptors and TFs
#The goal is to decrease the noise when extracting the networks signalling cell type specific, instead of using all DEGs for each cell cluster
## scphuego uses the three layers of signal propagation: receptors, TFs and DEGs
##Load BP pathways and Genes
BP_Genes_GO_BP <- read_csv("EBI_course_2024/BP_Pathways/BP_All_genes.csv")
##choose the dataset for enrichment
dbs <- listEnrichrDbs()

#
# Define the groups and cell types

scRNA_brain <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")
cell_types <- levels(scRNA_brain$major.celltype)
groups <- levels(scRNA_brain$Pathology)
groups
###get the differentially expressed receptors by cell type
# Define the function to load tables and process data
load_tables <- function(cell_type, degs_file, ccc_file) {
  # Load merged CCC results and ensure it has the correct columns
  CCC_Brain_results <- read_csv(ccc_file, show_col_types = FALSE) %>%
    select(source,target,ligand.complex,receptor.complex,aggregate_rank)
  
  # Load DEGs and handle missing column names
  DEGs <- read_csv(degs_file, show_col_types = FALSE) %>%
    rename_with(~ ifelse(. == "...1", "gene", .), everything())
  
  # Get unique receptors and ligands from CCC_Brain_results
  unique_receptors <- distinct(CCC_Brain_results, receptor.complex)
 
  
  # Filter the DEGs for those that are receptors in CCC_Brain_results
  DEG_receptors <- DEGs %>%
    filter(gene %in% unique_receptors$receptor.complex)
  
  
  # Save the filtered data
  write_csv(DEG_receptors, paste0("EBI_course_2024/Table_receptors/", group, "_", cell_type, "_DEG_receptors.csv"))
}

# Loop over each group
for (group in groups) {
  ccc_file <- paste0("EBI_course_2024/Cell_cell_results/LIANA_scRNA_Brain_", group, ".csv")
  
  # Process each cell type for the group
  for (cell_type in cell_types) {
    degs_file <- paste0("EBI_course_2024/Table_DEGs/Brain_DEGs_", group, "_", cell_type, ".csv")
    
    # Check if the DEGs file exists before processing
    if (file.exists(degs_file)) {
      load_tables(cell_type, degs_file, ccc_file)
    } else {
      message(paste("DEGs file for cell type", cell_type, "and group", group, "does not exist. Skipping..."))
    }
  }
}

# Define the function to get enriched genes Fro receptors and TFs by pathology
get_enriched_genes <- function(receptor_file, tf_RNA_file,tf_ATAC_file,cell_type, output_file) {
  
  # Read receptor data
  receptors <- read_csv(receptor_file) %>%
    select(gene, avg_log2FC)
  colnames(receptors) <- c("SYMBOL","logFC")
  
  # Read TF data
  tfs_RNA <- read_csv(tf_RNA_file) %>%
    filter(condition == cell_type) %>%
    select(source, score) %>%
    arrange(desc(score)) %>%
    slice_head(n = 10)
  colnames(tfs_RNA) <- c("SYMBOL", "score")
  
  tfs_ATAC <- read_csv(tf_ATAC_file) %>%
    filter(cell_type == !!cell_type) %>%
    select(motif.name, fold.enrichment) %>%
    arrange(desc(fold.enrichment)) %>%
    slice_head(n = 10)
  
  
  # Combine receptor and TF data
  rec_tf <- bind_rows(receptors, tfs_RNA)
  rec_tf <- bind_rows(rec_tf, tfs_ATAC)
  
  # Perform enrichment analysis
  enrich_results <- enrichr(rec_tf$SYMBOL, database = "GO_Biological_Process_2023")
  
  # Filter for significant results
  significant_results <- enrich_results$GO_Biological_Process_2023 %>%
    filter(Adjusted.P.value < 0.05) %>%
    pull(Term) # Extracting just the Term column
  go_terms <- gsub(".*\\(GO:(\\d+)\\).*", "GO:\\1", significant_results)
  
  
  # Get genes associated with significant pathways
  significant_genes <- BP_Genes_GO_BP %>%
    filter(GO_Terms %in% go_terms)
  
  # Write significant genes and pathways to file
  write_csv(significant_genes, output_file)
  
}

# Loop over each group and cell type
for (group in groups) {
  for (cell_type in cell_types) {
    # Construct file paths for receptors and TF activities
    receptor_file <- paste0("EBI_course_2024/Table_receptors/", group,"_", cell_type, "_DEG_receptors.csv")
    tf_RNA_file <- paste0("EBI_course_2024/TF_activity/scRNA_brain_", group, "_collectTRI.csv")
    tf_ATAC_file <- paste0("EBI_course_2024/TF_activity/scATAC_brain_", group, ".csv")
    
    # Define output files
    output_file <- paste0("EBI_course_2024/BP_Pathways/Recp_TFs_", cell_type, "_", group, "_enriched.csv")
    
    # Check if receptor and TF activity files exist
    if (file.exists(receptor_file) && file.exists(tf_RNA_file) && file.exists(tf_ATAC_file)) {
      # Run the function
      get_enriched_genes(receptor_file, tf_RNA_file,tf_ATAC_file, cell_type, output_file)
    } else {
      message(paste("Files for group", group, "and cell type", cell_type, "do not exist. Skipping..."))
    }
  }
}
###########Now get the TF, Receptors and DEGs present on those pathways to generate the networks
process_cell_type_group <- function(cell_type, group) {
  # Construct file paths for DEGs, receptors, and TF activities
  degs_file <- paste0("EBI_course_2024/Table_DEGs/Brain_DEGs_", group, "_", cell_type, ".csv")
  genes_paths_file <- paste0("EBI_course_2024/BP_Pathways/Recp_TFs_", cell_type, "_", group, "_enriched.csv")
  receptors_file <- paste0("EBI_course_2024/Table_receptors/", group,"_", cell_type, "_DEG_receptors.csv")
  tf_RNA_file <- paste0("EBI_course_2024/TF_activity/scRNA_brain_", group, "_collectTRI.csv")
  tf_ATAC_file <- paste0("EBI_course_2024/TF_activity/scATAC_brain_", group, ".csv")
  
  # Check if the necessary files exist
  if (file.exists(degs_file) && file.exists(genes_paths_file) && file.exists(receptors_file) && file.exists(tf_RNA_file) && file.exists(tf_ATAC_file)) {
    # Load necessary data
    uniprot_to_gene <- read.table("EBI_course_2024/Network_files/uniprot_to_gene.tab", header = F, sep = "\t") %>%
      rename(Protein = V1, SYMBOL = V2) %>%
      distinct()  # Ensure no duplicate entries
    UNIPROT_ID <- read_tsv(file = 'EBI_course_2024/Network_files/uniprot_ID.tsv') %>%
      rename(Protein = colnames(.)[1]) %>%
      distinct()  # Ensure no duplicate entries
    uniprot_to_gene <- right_join(uniprot_to_gene, UNIPROT_ID, by = "Protein")
    
    # Ensure uniprot_to_gene has unique SYMBOL entries
    uniprot_to_gene <- uniprot_to_gene %>%
      group_by(SYMBOL) %>%
      summarise(Protein = first(Protein))  # Resolve many-to-many relationships
    
    # Load DEGs
    degs_data <- read_csv(degs_file) %>%
      select(SYMBOL = 1, RWRS = avg_log2FC)  # Update column names as per your file
    
    # Load genes from enriched pathways
    genes_paths_data <- read_csv(genes_paths_file) %>%
      select(SYMBOL = SYMBOL)  # Update column names as per your file
    degs_filtered <- inner_join(degs_data, genes_paths_data, by = "SYMBOL")
    
    # Load receptors
    receptors_data <- read_csv(receptors_file) %>%
      select(SYMBOL = gene, RWRS = avg_log2FC)  # Update column names as per your file
    degs_filtered <- degs_filtered %>%
      filter(!SYMBOL %in% receptors_data$SYMBOL)
    
    # Load TF data
    tfs_RNA <- read_csv(tf_RNA_file) %>%
      filter(condition == cell_type) %>%
      select(source, score) %>%
      arrange(desc(score)) %>%
      slice_head(n = 10)
    colnames(tfs_RNA) <- c("SYMBOL", "RWRS")
    
    tfs_ATAC <- read_csv(tf_ATAC_file) %>%
      filter(cell_type == !!cell_type) %>%
      select(motif.name, fold.enrichment) %>%
      arrange(desc(fold.enrichment)) %>%
      slice_head(n = 10)
    colnames(tfs_ATAC) <- c("SYMBOL", "RWRS")
    tfs_data <- bind_rows(tfs_RNA, tfs_ATAC)
    
    # Merge tables
    merged_data <- full_join(degs_filtered, receptors_data) %>%
      full_join(tfs_data) %>%
      inner_join(uniprot_to_gene, by = "SYMBOL")%>%
      distinct(SYMBOL, .keep_all = TRUE)
    
    # Save the data
    output_file_csv <- paste0("EBI_course_2024/seed_nodes_scphuEGO/", cell_type, "_", group, "_seed_nodes_phuego.csv")
    output_file_txt <- paste0("EBI_course_2024/seed_nodes_scphuEGO/", cell_type, "_", group, "_seed_nodes_phuego.txt")
    write_csv(merged_data, output_file_csv)
    merged_data_2 <- merged_data[, c("Protein", "RWRS")]
    write.table(merged_data_2, output_file_txt, col.names = F, sep = "\t", row.names = F, quote = F)
  } else {
    message("Files for group ", group, " and cell type ", cell_type, " do not exist. Skipping...")
  }
}


####
# Loop over each group and cell type
for (group in groups) {
  for (cell_type in cell_types) {
    process_cell_type_group(cell_type, group)
  }
}

