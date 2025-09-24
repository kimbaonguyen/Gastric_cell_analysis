#########Doublet Detection#########
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)

library(DoubletFinder)
library(Seurat)
library(remotes)

ls("package:DoubletFinder")

setwd("F:/Kim_Nguyen/Secondary_Reclustering_of_Myeloid_cells")                                                              
save.image(file = "F:/Kim_Nguyen/Secondary_Reclustering_of_Myeloid_cells/Reclustering_Myeloid_cells_06032025_Doublet_Detection.RData")
load("F:/Kim_Nguyen/Secondary_Reclustering_of_Myeloid_cells/Reclustering_Myeloid_cells_06032025_Doublet_Detection.RData")

## Load the .rds file
Reclustering_Doublet_Detection <- readRDS("Reclustering_Myeloid_cells_06022025_v3.rds")
gastric_single_cell <- readRDS("F:/Kim_Nguyen/all_object_QC_UMAP (1).rds")

## Add "Annotation of cell type" column 
DimPlot(gastric_single_cell, reduction = "umap")
gastric_single_cell@meta.data$Annotation_of_cell_types = "NA"
DimPlot(gastric_single_cell, reduction = "umap", group.by = c("Annotation_of_cell_types", "SCT_snn_res.0.2"))  

gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "0"), "Annotation_of_cell_types"] = "ILC2"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "1"), "Annotation_of_cell_types"] = "Macrophage"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "2"), "Annotation_of_cell_types"] = "PreDC2s.DC2s.ProDC3s.DC3s"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "3"), "Annotation_of_cell_types"] = "CD4 or CD8"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "4"), "Annotation_of_cell_types"] = "gamma delta T cells"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "5"), "Annotation_of_cell_types"] = "Monocyte"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "6"), "Annotation_of_cell_types"] = "NK cell"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "7"), "Annotation_of_cell_types"] = "B cells"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "8"), "Annotation_of_cell_types"] = "Treg (Regulatory T) cell"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "9"), "Annotation_of_cell_types"] = "DC1"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "10"), "Annotation_of_cell_types"] = "CCR7_and_DCs"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "11"), "Annotation_of_cell_types"] = "Mast_cells"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "12"), "Annotation_of_cell_types"] = "Plasma"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "13"), "Annotation_of_cell_types"] = "endothelial cells and stromal cells"
gastric_single_cell@meta.data[which(gastric_single_cell@meta.data$SCT_snn_res.0.2 == "14"), "Annotation_of_cell_types"] = "Mast_cells_2"

###################Re-clustering analysis on Myeloid Cells #################
Idents(gastric_single_cell) <- gastric_single_cell@meta.data$Annotation_of_cell_types
Macrophage.Monocyte.DCs <- subset(gastric_single_cell, idents = c("Macrophage", "Monocyte","PreDC2s.DC2s.ProDC3s.DC3s", "DC1", "CCR7_and_DCs" ))

##############Preprocess Seurat Object#################
DefaultAssay(Macrophage.Monocyte.DCs) <- "RNA"
Macrophage.Monocyte.DCs <- NormalizeData(Macrophage.Monocyte.DCs, normalization.method = "LogNormalize", scale.factor = 10000)
Macrophage.Monocyte.DCs <- FindVariableFeatures(Macrophage.Monocyte.DCs)
Macrophage.Monocyte.DCs <- ScaleData(Macrophage.Monocyte.DCs)
Macrophage.Monocyte.DCs <- RunPCA(Macrophage.Monocyte.DCs)  ###DoubletFinder works on PCA space
Macrophage.Monocyte.DCs <- FindNeighbors(Macrophage.Monocyte.DCs, dims = 1:10) 
Macrophage.Monocyte.DCs <- FindClusters(Macrophage.Monocyte.DCs, resolution = 0.5)
Macrophage.Monocyte.DCs <- RunUMAP(Macrophage.Monocyte.DCs, dims = 1:10)

#Macrophage.Monocyte.DCs <- NormalizeData(Macrophage.Monocyte.DCs, normalization.method = Macrophage.Monocyte.DCs@commands$NormalizaData.RNA@params$normalization.method, scale.factor = Macrophage.Monocyte.DCs@commands$NormalizeData.SCT@params$scale.factor, margin = Macrophage.Monocyte.DCs@commands$NormalizeData.SCT@params$margin)
#SCTransform(Macrophage.Monocyte.DCs)
#Macrophage.Monocyte.DCs <- NormalizeData(Macrophage.Monocyte.DCs, normalization.method = Macrophage.Monocyte.DCs@commands$NormalizaData.SCT@params$normalization.method, scale.factor = Macrophage.Monocyte.DCs@commands$NormalizeData.RNA@params$scale.factor, margin = Macrophage.Monocyte.DCs@commands$NormalizeData.RNA@params$margin)
#Macrophage.Monocyte.DCs <- RunPCA(Macrophage.Monocyte.DCs, assay = "SCT")  ###DoubletFinder works on PCA space
#Macrophage.Monocyte.DCs <- FindNeighbors(Macrophage.Monocyte.DCs, dims = 1:10, assay = "SCT") 
#Macrophage.Monocyte.DCs <- FindClusters(Macrophage.Monocyte.DCs, resolution = 0.5, assay = "SCT")
#Macrophage.Monocyte.DCs <- RunUMAP(Macrophage.Monocyte.DCs, dims = 1:10, assay = "SCT")

options(future.globals.maxSize = 8 * 1024^3) #8GB added due to error of future functions

                                ##################Estimate number of Doublets################
#######Find optimal pK (parameter sweep)######
sweep.res.list.Myeloid.cells <- paramSweep(Macrophage.Monocyte.DCs, PCs = 1:10, sct = FALSE) 
#sct = T means using SCTransform-based normalization 
#sct = F means  standard normalization

#Summarize the results for classification of each tested pK value
sweep.stats.Myeloid.cells <- summarizeSweep(sweep.res.list.Myeloid.cells, GT = FALSE)

#Find the optimal pK value with a table of pK values and their corresponding BCmetrics (higher = better doublet discrimination)
bcmvn.Myeloid.cells <- find.pK(sweep.stats.Myeloid.cells)#pK is a turning parameter that controls how artificial doublets are simulated
best.pK <- as.numeric(as.character(bcmvn.Myeloid.cells[which.max(bcmvn.Myeloid.cells$BCmetric), "pK"]))
#>NULL, meaning sweep results were generated correctly and contain data
print(bcmvn.Myeloid.cells)

#Homotypic Doublet Proportion Estimate
homo.doublet.est <- modelHomotypic(Macrophage.Monocyte.DCs@meta.data$RNA_snn_res.0.5)
nExp_poi.Myeloid.cells <- round(0.07 * nrow(Macrophage.Monocyte.DCs@meta.data)) 
# 7% of total cells
### Common rule of thumb is 5-10% of total cells
nExp_poi.adj.Myeloid.cells <- round(nExp_poi.Myeloid.cells*(1-homo.doublet.est)) 

#Run DoubletFinder with classification stringencies
Macrophage.Monocyte.DCs.doubletFinder <- doubletFinder(Macrophage.Monocyte.DCs, PCs = 1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi.Myeloid.cells, reuse.pANN = NULL, sct = FALSE) #1
colnames(Macrophage.Monocyte.DCs.doubletFinder@meta.data)
Macrophage.Monocyte.DCs.doubletFinder <- doubletFinder(Macrophage.Monocyte.DCs.doubletFinder, PCs = 1:10, pN = 0.25, pK = best.pK, nExp = nExp_poi.adj.Myeloid.cells, reuse.pANN = "pANN_0.25_0.08_1113", sct = FALSE) #2

#Extract doublets from 1st classification column
doublets_1113 <- rownames(Macrophage.Monocyte.DCs.doubletFinder@meta.data[Macrophage.Monocyte.DCs.doubletFinder@meta.data$DF.classifications_0.25_0.08_1113 == "Doublet",])

#Extract doublets from 2nd classification column
doublets_983 <- rownames(Macrophage.Monocyte.DCs.doubletFinder@meta.data[Macrophage.Monocyte.DCs.doubletFinder@meta.data$classifications_0.25_0.08_983 == "Doublet",])

#Visualize doublets on UMAP
Idents(Macrophage.Monocyte.DCs.doubletFinder) <- "DF.classifications_0.25_0.08_983"
DimPlot(Macrophage.Monocyte.DCs.doubletFinder, reduction = "umap", group.by = c("DF.classifications_0.25_0.08_983", "Annotation_of_cell_types"))
DimPlot(Macrophage.Monocyte.DCs.doubletFinder, reduction = "umap", group.by = c("DF.classifications_0.25_0.08_983"))

Idents(Macrophage.Monocyte.DCs.doubletFinder) <- "DF.classifications_0.25_0.08_1113"
DimPlot(Macrophage.Monocyte.DCs.doubletFinder, reduction = "umap", group.by = c("DF.classifications_0.25_0.08_1113", "Annotation_of_cell_types"))
