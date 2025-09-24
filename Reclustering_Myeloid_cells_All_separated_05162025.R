##### Monocytes #####
Monocyte <- subset(gastric_single_cell, idents = "Monocyte") ##extracting only cells with identity "Monocyte"
     table(Idents(gastric_single_cell)) ## Check labels
     Monocyte <- NormalizeData(Macrophage) ## optional
     Monocyte <- FindVariableFeatures(Monocyte) ##Find variable genes ## optional
     Monocyte <- ScaleData(Monocyte)
     
Monocyte <- RunPCA(Monocyte)
     DimPlot(Monocyte, reduction = "pca") ## Viewing PCA Plot
     
Monocyte <- FindNeighbors(Monocyte, dims = 1:10) ## Find neighbors
     
Monocyte <- FindClusters(Monocyte, resolution = 0.2) ##Find clusters
     
Monocyte <- RunUMAP(Monocyte, dims = 1:10)
     DimPlot(Monocyte, reduction = "umap", label = TRUE)

##### Pre-DC2s, DC2s, proDC3s, DC3s #####
PreDC2s.DC2s.ProDC3s.DC3s <- subset(gastric_single_cell, idents = "PreDC2s.DC2s.ProDC3s.DC3s")     
     table(Idents(gastric_single_cell)) ## Check labels
     PreDC2s.DC2s.ProDC3s.DC3s <- NormalizeData(PreDC2s.DC2s.ProDC3s.DC3s) ## optional
     PreDC2s.DC2s.ProDC3s.DC3s <- FindVariableFeatures(PreDC2s.DC2s.ProDC3s.DC3s) ##Find variable genes ## optional
     PreDC2s.DC2s.ProDC3s.DC3s <- ScaleData(PreDC2s.DC2s.ProDC3s.DC3s)
     
PreDC2s.DC2s.ProDC3s.DC3s <- RunPCA(PreDC2s.DC2s.ProDC3s.DC3s)
     DimPlot(PreDC2s.DC2s.ProDC3s.DC3s, reduction = "pca") ## Viewing PCA Plot
     
PreDC2s.DC2s.ProDC3s.DC3s <- FindNeighbors(PreDC2s.DC2s.ProDC3s.DC3s, dims = 1:10) ## Find neighbors
     
PreDC2s.DC2s.ProDC3s.DC3s <- FindClusters(PreDC2s.DC2s.ProDC3s.DC3s, resolution = 0.2) ##Find clusters

PreDC2s.DC2s.ProDC3s.DC3s <- RunUMAP(PreDC2s.DC2s.ProDC3s.DC3s, dims = 1:10)
     DimPlot(PreDC2s.DC2s.ProDC3s.DC3s, reduction = "umap", label = TRUE)

##### DC1 #####    
DC1 <- subset(gastric_single_cell, idents = "DC1")     
     table(Idents(gastric_single_cell)) ## Check labels
     DC1 <- NormalizeData(DC1) ## optional
     DC1 <- FindVariableFeatures(DC1) ##Find variable genes ## optional
     DC1 <- ScaleData(DC1)
     
DC1 <- RunPCA(DC1)
     DimPlot(DC1, reduction = "pca") ## Viewing PCA Plot
     
DC1 <- FindNeighbors(DC1, dims = 1:10) ## Find neighbors
     
DC1 <- FindClusters(DC1, resolution = 0.2) ##Find clusters
     
DC1 <- RunUMAP(DC1, dims = 1:10)
     DimPlot(DC1, reduction = "umap", label = TRUE)
     
##### CCR7 and DCs ##### 
CCR7_and_DCs <- subset(gastric_single_cell, idents = "CCR7_and_DCs")    
     table(Idents(gastric_single_cell))     
     CCR7_and_DCs <- NormalizeData(CCR7_and_DCs)
     CCR7_and_DCs <- FindVariableFeatures(CCR7_and_DCs)
     CCR7_and_DCs <- ScaleData(CCR7_and_DCs)
     
CCR7_and_DCs <- RunPCA(CCR7_and_DCs)
     DimPlot(CCR7_and_DCs, reduction = "pca")
     
CCR7_and_DCs <- FindNeighbors(CCR7_and_DCs, dims = 1:10)

CCR7_and_DCs <- FindClusters(CCR7_and_DCs, resolution = 0.2)
     
CCR7_and_DCs <- RunUMAP(CCR7_and_DCs, dims = 1:10)
     DimPlot(CCR7_and_DCs, reduction = "umap", label = TRUE)   
     
##### Mast cells ##### 
Mast_cells <- subset(gastric_single_cell, idents = "Mast_cells")
     Mast_cells <- NormalizeData(Mast_cells)
     Mast_cells <- FindVariableFeatures(Mast_cells)
     Mast_cells <- ScaleData(Mast_cells)
     
Mast_cells <- RunPCA(Mast_cells)
     DimPlot(Mast_cells, reduction = "pca")
    
Mast_cells <- FindNeighbors(Mast_cells, dims = 1:10)

Mast_cells <- FindClusters(Mast_cells, resolution = 0.2)

Mast_cells <- RunUMAP(Mast_cells, dims = 1:10)     
     DimPlot(Mast_cells, reduction = "umap", label = TRUE)
     
##### Mast cells-2 #####   
Mast_cells_2 <- subset(gastric_single_cell, idents = "Mast_cells_2")
     Mast_cells_2 <- NormalizeData(Mast_cells_2)
     Mast_cells_2 <- FindVariableFeatures(Mast_cells_2)
     Mast_cells_2 <- ScaleData(Mast_cells_2)
     
Mast_cells_2 <- RunPCA(Mast_cells_2)
     DimPlot(Mast_cells_2, reduction = "pca")
     
Mast_cells_2 <- FindNeighbors(Mast_cells_2, dims = 1:10)

Mast_cells_2 <- FindClusters(Mast_cells_2, resolution = 0.2)

Mast_cells_2 <- RunUMAP(Mast_cells_2, dims = 1:10)
     DimPlot(Mast_cells_2, reduction = "umap", label = TRUE)
     