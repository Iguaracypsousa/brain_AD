setwd("/nfs/research/petsalaki/users/iguaracy/scRNA_seq_analysis/R/AD_Brain_analysis/")
#!/hps/software/users/petsalaki/users/tkoutsandreas/R-4.1.2/bin/Rscript

# AUTHOR: Iguaracy Sousa, Postdoc EMBL-EBI
#
# Human adult brain scRNAseq and scATACseq  data
# 
# INPUT: rds file
#
# OUTPUT: integrated data
suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("EnsDb.Hsapiens.v86")
  library("ggplot2")
  library("cowplot")
  library("RColorBrewer")
})
###Seurat pipeline
###Load files 
scRNA_brain <- readRDS("EBI_course_2024/rds_files/scRNA_brain_seurat_QC.rds")
scATAC_brain <- readRDS("EBI_course_2024/rds_files/scRNA_ATAC_brain_seurat_QC.rds")
####satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
p1 <- DimPlot(scRNA_brain, group.by = "major.celltype", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(scATAC_brain,group.by = "Celltype1", label = TRUE) + NoLegend() + ggtitle("ATAC")
##
jpeg("EBI_course_2024/figs/scRNA_ATAC_brain_cell_types_.jpeg",width = 11, height = 7, units = 'in', res=300)
p1 + p2
dev.off()


##Identifying anchors between scRNA-seq and scATAC-seq datasets
# normalize gene activities
DefaultAssay(scATAC_brain) <- "ACTIVITY"
scATAC_brain <- NormalizeData(scATAC_brain)
scATAC_brain <- ScaleData(scATAC_brain, features = rownames(scATAC_brain))

transfer.anchors <- FindTransferAnchors(reference = scRNA_brain, query = scATAC_brain, features = VariableFeatures(object = scRNA_brain),
                                        reference.assay = "RNA", query.assay = "ACTIVITY",reduction = "cca")

#Annotate scATAC-seq cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = scRNA_brain$major.celltype,
                                     weight.reduction = scATAC_brain[["lsi"]], dims = 2:30)

scATAC_brain <- AddMetaData(scATAC_brain, metadata = celltype.predictions)

###

scATAC_brain$annotation_correct <- scATAC_brain$predicted.id == scATAC_brain$Celltype1
p1 <- DimPlot(scATAC_brain, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(scATAC_brain, group.by = "Celltype1", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2

####
predictions <- table(scATAC_brain$Celltype1, scATAC_brain$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(scATAC_brain$Celltype1 == scATAC_brain$predicted.id))
incorrect <- length(which(scATAC_brain$Celltype1 != scATAC_brain$predicted.id))
data <- FetchData(scATAC_brain, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2

######Merge process
# first add dataset-identifying metadata
scATAC_brain$dataset <- "ATAC"
scRNA_brain$dataset <- "RNA"
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(scRNA_brain)
refdata <- GetAssayData(scRNA_brain, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = scATAC_brain[["lsi"]],
                           dims = 2:30)
scATAC_brain[["RNA"]] <- imputation

merged_RNA_ATAC <- merge(x = scRNA_brain, y = scATAC_brain)

# process the combined dataset
# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
merged_RNA_ATAC <- ScaleData(merged_RNA_ATAC, features = genes.use, do.scale = FALSE)
merged_RNA_ATAC <- RunPCA(merged_RNA_ATAC, features = genes.use, verbose = FALSE)
merged_RNA_ATAC <- RunUMAP(merged_RNA_ATAC, dims = 1:30)

###
# Assuming you have your predictions in a format that can be easily added to the Seurat object
# Add the predicted cell types to the scATAC_brain object
scATAC_brain$predicted_celltype <- scATAC_brain$predicted.id # your predicted cell type vector

# Since you've already merged the datasets, you'll need to add the predicted cell type to the merged object
# First, we need to find out which cells in the merged object come from the ATAC dataset
# Assuming that 'dataset' metadata is preserved in the merged object

# Get indices of ATAC cells in the merged object
atac_cells <- WhichCells(merged_RNA_ATAC, expression = dataset == "ATAC")

# Add the predictions to the merged object
merged_RNA_ATAC$major.celltype[atac_cells] <- scATAC_brain$predicted_celltype
Idents(scATAC_brain) <- scATAC_brain$predicted_celltype
# Now, you can visualize or analyze the merged object with the added predictions
# For example, replotting UMAP with the new cell type annotations
DimPlot(merged_RNA_ATAC, group.by = "predicted_celltype")

DimPlot(merged_RNA_ATAC, group.by = c("dataset", "major.celltype"))


merged_RNA_ATAC_metadata <- merged_RNA_ATAC@meta.data
######Gene for microgolia
jpeg("figs/FeaturePlot_ECs_All_markers_and_VWF.jpeg",width = 13, height = 9, units = 'in', res=150)
FeaturePlot(merged_RNA_ATAC,c("SPI1","IRF8","CEBPA"),min.cutoff = 0,raster=FALSE)
dev.off()



### ploting Cells marekrs
DefaultAssay(scRNA_brain) <- "RNA"

SPI1_RNA <- FeaturePlot(scRNA_brain,"SPI1",min.cutoff = 0,raster=FALSE) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("SPI1 - Microglia - RNA")


DefaultAssay(scATAC_brain) <- "ACTIVITY"

SPI1_ATAC <- FeaturePlot(scATAC_brain,"SPI1",min.cutoff = 0,raster=FALSE) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("SPI1 - Microglia - ATAC")



jpeg("EBI_course_2024/figs/UMAP_Microglia_Marker_RNA_ATAC.jpeg",width = 10, height =7, units = 'in', res=200)
plot_grid(SPI1_RNA,SPI1_ATAC)
dev.off()
