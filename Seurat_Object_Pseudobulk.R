library(dplyr)
library(Seurat)
library(patchwork)

## Set the path to the .rds file
file_path <- "F:/Kim_Nguyen/all_object_QC_UMAP (1).rds"
## Load the .rds file
gastric_single_cell <- readRDS(file_path)
## Check the data contents
metadata <- gastric_single_cell@meta.data

# Assessing Meta Data for Grouping
gastric_sc_metadata = gastric_single_cell@meta.data

gastric_single_cell$celltype = Idents(gastric_single_cell)
# gastric_single_cell$sample = gastric_single_cell$orig.ident

# Aggregate counts by both cell type and sample
pseudobulk_gastric_counts = AggregateExpression(gastric_single_cell,
                                                group.by = c("Annotation_of_cell_types","sample"),
                                                assays = "RNA",
                                                slot = "counts")


# Define mapping of g0, g1, g2... to actual cell types
metadata$Annotation_of_cell_types <- recode(metadata$Annotation_of_cell_types,
                             "g0" = "ILC2",
                             "g1" = "Macrophage",
                             "g2" = "Pre-DC2s, DC2s, proDC3s, DC3s",
                             "g3" = "CD4/CD8 T cells",
                             "g4" = "Gamma T cells",
                             "g5" = "Monocytes",
                             "g6" = "NK cells",
                             "g7" = "B cells", 
                             "g8" = "Treg cell",
                             "g9" = "DC1",
                             "g10" = "CCR7, DCs",
                             "g11" = "Mast cells", 
                             "g12" = "Plasma",
                             "g13" = "endothelical cells, stomal cells",
                             "g14" = "Mast cells-2")


pseudobulk_gastric_matrix = pseudobulk_gastric_counts$RNA

################### Prepare Pseudobulk Matric for DE Analysis #################

########## Convert Wide Matrix to Long Format  ################

# Convert to dataframe format
pseudobulk_gastric_df = as.data.frame(as.matrix(pseudobulk_gastric_matrix))
View(pseudobulk_gastric_df)
pseudobulk_gastric_df$gene = rownames(pseudobulk_gastric_df)

write.csv(pseudobulk_gastric_matrix, file = "pseudobulk_matrix.csv", row.names = TRUE)


# Reshape: Separate "Annotation_of_cell_types_sample" into separate columns
pseudobulk_gastric_df_long = pseudobulk_gastric_df %>%
  tidyr::pivot_longer(cols = -gene, names_to = "Annotation_of_cell_types_sample", values_to = "counts") %>%
  tidyr::separate(Annotation_of_cell_types_sample, into = c("Annotation_of_cell_types", "sample"), sep = "_")


################### Perform DE Analysis with DESeq2 ###################
# Create sample metadata dataframe
coldata = unique(pseudobulk_gastric_df_long[, c("sample", "Annotation_of_cell_types")])
View(coldata)
coldata$Annotation_of_cell_types_sample = paste(coldata$Annotation_of_cell_types, coldata$sample, sep = "_")
coldata$condition <- ifelse(grepl("MNU", coldata$sample), "treated", "untreated")
coldata$genotype <- ifelse(grepl("KO", coldata$sample), "KO", "WT")
View(coldata)
rownames(coldata) <- coldata$Annotation_of_cell_types_sample
View(coldata)
write.csv(coldata, file = "sample_metadata.csv", row.names = TRUE)

setwd("F:/Kim_Nguyen/")


# Convert pseudobulk matrix into DESeq2 format
dds <- DESeqDataSetFromMatrix(countData = pseudobulk_gastric_matrix, 
                              colData = coldata, 
                              design = ~ genotype + condition)
dds

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Get the results
res <- results(dds)

# View the summary of results (default: log2 fold changes for condition/genotype comparisons)
summary(res)

# MA Plot
plotMA(res, main = "DESeq2 MA Plot", ylim = c(-5, 5))


write.csv(pseudobulk_gastric_matrix, file = "pseudobulk_matrix.csv", row.names = TRUE)

