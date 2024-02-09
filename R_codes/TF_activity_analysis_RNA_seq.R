setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq data
# 
# INPUT: rds file
#
# OUTPUT: TF activity
suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("dorothea")
  library("readr")
  library("purrr")
  library("decoupleR")
 })


# Load the Seurat object
scRNA_brain <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")


######Identify differential expressed genes across conditions
# Calculate the average expression by cell type
avg_exp <- AverageExpression(object = scRNA_brain, assays = "RNA",slot = "counts", group.by = c("major.celltype", "Pathology"))
# Extract the expression matrix
exp_matrix <- avg_exp$RNA

# Apply log1p transformation
log_avg_exp <- log1p(exp_matrix)

# Check the first few entries of the log-transformed data
print(head(log_avg_exp))

# Calculate the sum of expression values for each gene
gene_sums <- rowSums(log_avg_exp)
# Plot a histogram of the row sums
hist(gene_sums, breaks = 100, main = "Distribution of Gene Expression Sums")

# Set a minimum expression threshold
min_exp_threshold <- 1

# Filter out the low expressed genes
filtered_log_avg_exp <- log_avg_exp[gene_sums > min_exp_threshold, ]
filt_gene_sums <- rowSums(filtered_log_avg_exp)
hist(filt_gene_sums, breaks = 100, main = "Distribution of Gene Expression Sums")

####Load dorothea regulons
Regulon_file <- read_csv("EBI_course_2024/TF_activity/Regulon_file_get_collectri.csv")

# Run wmean
scRNA_brain_acts <- run_wmean(mat=filtered_log_avg_exp, net=Regulon_file, .source='source', .target='target',
                                 .mor='mor', times = 1000, minsize = 25)

scRNA_brain_acts <- scRNA_brain_acts %>%
  filter(statistic == 'norm_wmean', p_value < 0.05)
scRNA_brain_acts

# Filter rows for early_AD
early_AD <- scRNA_brain_acts %>%
  filter(grepl("earlyAD", condition))

# Filter rows for non_AD
non_AD <- scRNA_brain_acts %>%
  filter(grepl("nonAD", condition))

# Filter rows for lateAD
lateAD <- scRNA_brain_acts %>%
  filter(grepl("lateAD", condition))

write_csv(early_AD, "EBI_course_2024/TF_activity/scRNA_brain_earlyAD_collectTRI.csv")
write_csv(non_AD, "EBI_course_2024/TF_activity/scRNA_brain_nonAD_collectTRI.csv")
write_csv(lateAD, "EBI_course_2024/TF_activity/scRNA_brain_lateAD_collectTRI.csv")
