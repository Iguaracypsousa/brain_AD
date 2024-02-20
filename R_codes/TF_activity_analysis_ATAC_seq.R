setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
Sys.setenv(LD_LIBRARY_PATH = paste("/hps/software/users/petsalaki/users/iguaracy/miniconda3/envs/my_r_env/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
dyn.load("/hps/software/users/petsalaki/users/iguaracy/miniconda3/envs/my_r_env/lib/libgsl.so.25")

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
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("JASPAR2020")
  library("TFBSTools")
  library("patchwork")
  library("motifmatchr")
  library("ggseqlogo")
})
###
# Load your scATAC-seq data
scATAC_brain <- readRDS("EBI_course_2024/rds_files/scATAC_brain_ACTIVITY_norm.rds")
Idents(scATAC_brain) <- scATAC_brain$Celltype1
active.ident <- scATAC_brain@active.ident
####RNA cell annotation
new.cluster.ids <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc", "Vas")
current.cluster.ids <- levels(scATAC_brain)

scATAC_brain@active.ident <- plyr::mapvalues(x = scATAC_brain@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(scATAC_brain, pt.size = 0.001,label = T, raster=FALSE,label.size = 4)

cell_types <- levels(scATAC_brain@active.ident)
# Assuming 'scATAC_brain' is your Seurat object with a ChromatinAssay named 'ATAC'
# First, set the default assay if not already set
DefaultAssay(scATAC_brain) <- 'ATAC'

# Split the data by 'Pathology'
data_split <- SplitObject(scATAC_brain, split.by = "Pathology")


# Plot the UMAP
p1 <- DimPlot(scATAC_brain, label = TRUE, pt.size = 0.1) + NoLegend()
p1



# Create a PWM object from JASPAR
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE) # 9606 is the NCBI Taxonomy ID for Homo sapiens
)

# Convert PFMatrixList to PWMList
pwm <- lapply(pfm, function(x) toPWM(x, type = "prob"))

# Iterate over each subset
data_split_motifs <- lapply(data_split, function(subset) {
  # Add motifs to each subset
  AddMotifs(
    object = subset,
    pfm,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
})
# Initialize an empty list to store results
da_peaks_results <- list()

# Loop through each cell type
for(cell_type in cell_types) {
  da_peaks_results[[cell_type]] <- lapply(data_split_motifs, function(subset) {
    # Check if the cell type has enough cells
    if(sum(Idents(subset) == cell_type) >= 3) {
      return(FindMarkers(
        object = subset,
        ident.1 = cell_type,
        ident.2 = NULL, # Compare to all other cells
        only.pos = TRUE,
        test.use = 'LR',
        min.pct = 0.25,
        latent.vars = 'nCount_peaks'
      ))
    } else {
      # Return NULL or an informative message if the cell group is too small
      return(NULL)
    }
  })
}
# get top differentially accessible peaks

# After collecting DA peaks, now identify and store top DA peaks based on p-value threshold
top_da_peaks_results <- list()

# Loop through each cell type and condition to filter top DA peaks
for(cell_type in names(da_peaks_results)) {
  top_da_peaks_results[[cell_type]] <- lapply(da_peaks_results[[cell_type]], function(markers) {
    if(!is.null(markers)) {
      # Filter for peaks with p_val < 0.005
      top_peaks <- markers[markers$p_val < 0.005, ]
      return(rownames(top_peaks))
    } else {
      return(NULL)
    }
  })
}

######
# Assuming you have a list of Seurat objects named data_split_motifs with conditions as names
# And a list named top_da_peaks_results with the top differentially accessible peaks

enriched_motifs_results <- list()

# Loop through each condition in your data_split_motifs list
for (condition in names(data_split_motifs)) {
  # Loop through each cell type
  for (cell_type in names(top_da_peaks_results)) {
    # Check if there are any DA peaks for the current cell type and condition
    if (!is.null(top_da_peaks_results[[cell_type]][[condition]])) {
      # Retrieve the DA peaks for the current cell type and condition
      da_peaks <- top_da_peaks_results[[cell_type]][[condition]]
      
      # Ensure there are enough peaks to analyze
      if (length(da_peaks) > 0) {
        # Perform motif enrichment analysis using the FindMotifs function from Seurat
        enriched_motifs_results[[cell_type]][[condition]] <- FindMotifs(
          object = data_split_motifs[[condition]], 
          features = da_peaks
        )
      } else {
        # Handle the case where there are no significant DA peaks
        cat("No significant DA peaks for", cell_type, "in", condition, "\n")
        enriched_motifs_results[[cell_type]][[condition]] <- NULL
      }
    } else {
      # Handle the case where DA peaks data is missing for the cell type and condition
      cat("DA peaks data missing for", cell_type, "in", condition, "\n")
      enriched_motifs_results[[cell_type]][[condition]] <- NULL
    }
  }
}

# After the loop, enriched_motifs_results will contain the results of FindMotifs for all cell types and conditions

# Assuming the enriched_motifs_results for each condition is a dataframe
# Initialize an empty list for each condition to store the motif results
enriched_motifs_dfs <- list(
  nonAD = list(),
  earlyAD = list(),
  lateAD = list()
)

# The rest of your code goes here...
# Once you've run your FindMotifs analysis and have the enriched_motifs_results filled...

# Loop through the results and combine them into data frames for each condition
for (cell_type in cell_types) {
  for (condition in names(enriched_motifs_dfs)) {
    # Check if there are results for the current cell type and condition
    if (!is.null(enriched_motifs_results[[cell_type]]) && 
        !is.null(enriched_motifs_results[[cell_type]][[condition]])) {
      # Assume that the result can be coerced into a data frame
      motifs_df <- as.data.frame(enriched_motifs_results[[cell_type]][[condition]])
      
      # Add a column for the cell type to the data frame
      motifs_df$cell_type <- cell_type
      
      # Append the data frame to the appropriate list in enriched_motifs_dfs
      enriched_motifs_dfs[[condition]][[cell_type]] <- motifs_df
    }
  }
}

# Combine the lists of data frames into single data frames for each condition
final_dfs <- lapply(enriched_motifs_dfs, function(condition_list) {
  # Use bind_rows from dplyr to combine all data frames in the list
  bind_rows(condition_list)
})

# Name the list elements for clarity
names(final_dfs) <- c("nonAD", "earlyAD", "lateAD")
nonAD_TFs <- final_dfs$nonAD %>% filter(p.adjust < 0.05)

earlyAD_TFs <- final_dfs$earlyAD %>% filter(p.adjust < 0.05)
lateAD_TFs <- final_dfs$lateAD %>% filter(p.adjust < 0.05)

write_csv(nonAD_TFs, "EBI_course_2024/TF_activity/scATAC_brain_nonAD.csv")
write_csv(earlyAD_TFs, "EBI_course_2024/TF_activity/scATAC_brain_earlyAD.csv")
write_csv(lateAD_TFs, "EBI_course_2024/TF_activity/scATAC_brain_lateAD.csv")

###checking the motifs
#choose a cell type and condition / here Microglia and nonAD
MotifPlot(
  object = data_split_motifs$nonAD,
  motifs = head(rownames(enriched_motifs_results$Mic$nonAD))
)




