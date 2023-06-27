###### Building unsupervised reference #########


############ Weinkopf Dataset - Skin Leish ###############
library(Seurat)
library(Matrix)
library(ggplot2)
library(sctransform)
library(dplyr)
library(patchwork)
library(SingleR)
library(stringr)
library(ggpubr)
dir_data <- paste0(getwd(), "/data")
dir_results<- paste0(getwd(), "/results")
dir_figures<- paste0(getwd(), "/figures")

exp <- c("GSM5509913_high.W2.Inf.filtered_counts", "GSM5509914_high.W5.Naive.filtered_counts", "GSM5509915_high.W6.Naive.filtered_counts",  "GSM5509916_mid.W1.Inf.filtered_counts", "GSM5509917_mid.W2.Inf.filtered_counts", "GSM5509918_mid.W5.Naive.filtered_counts","GSM5509919_mid.W6.Naive.filtered_counts")
exp_cell_id <- c("high.W2", "high.W5", "high.W6", "mid1.W1", "mid2.W2", "mid1.W5", "mid2.W6")
exp_cell_id_complete <- c("high.W1","high.W2", "high.W5", "high.W6", "mid1.W1", "mid2.W2", "mid1.W5", "mid2.W6")

gse_p <- Read10X(paste0(dir_data, "/GSE181720_RAW/","GSM5509912_high.W1.Inf.filtered_counts"))
merged.list <- CreateSeuratObject(gse_p, project = "high.W1")
for(i in 1:length(exp)) {
    gse <- Read10X(paste0(dir_data, "/GSE181720_RAW/",exp[i]))
    seurat <- CreateSeuratObject(gse, project = exp_cell_id[i])
    merged.list <- c(merged.list, seurat)
}
seurat_merge <- merge(merged.list[[1]], c(merged.list[[2]], merged.list[[3]], merged.list[[4]], merged.list[[5]], merged.list[[6]], merged.list[[7]], merged.list[[8]]), add.cell.ids = exp_cell_id_complete, project = "leish_samples")
table(seurat_merge$orig.ident)

### SAVING THE RESULTS ###
saveRDS(seurat_merge, file = paste0(dir_results, "/leish_samples-dataset.rds"))
rm(gse, gse_p, seurat, exp, exp_cell_id, i, exp_cell_id_complete)

### READING THE RESULTS ###
getwd()
list.files()
seurat_merge <- readRDS(paste0(dir_results, "/leish_samples-dataset.rds"))
dim(seurat_merge)
table(seurat_merge$orig.ident)

#Load metadata
metadata <- read.delim(paste0(dir_data, "/GSE181720_combined.qc/metadata.tsv"), sep="", header=TRUE, quote="\"")
table(metadata$orig.ident)

#identify idents on seurat data
seurat_merge$nCount_RNA <- NULL
seurat_merge$nFeature_RNA <- NULL
rownames(metadata) <- paste0(rownames(metadata), "-1")

#join the two metadata tables
seurat_new <- AddMetaData(seurat_merge, metadata)

#Remove NA values from object
table(!is.na(seurat_new@meta.data$orig.ident))
seurat_data <- subset(seurat_new, subset = percent.mt != is.na(percent.mt))
immune_cells <- c("B cells", "Basophils", "Cytotoxic T cells", "Dendritic Cells 1", "Dendritic Cells 2", "Dermal Dendritic Cells", "DN T-lymphocytes", "Helper T Cells", "Inflammatory Monocytes", "Macrophage", "Mast cells", "Neutrophils", "NK cells", "Proliferating T Cells", "Resident Macrophage", "Schwann cells", "Tregs", "γδ-T cells")
weinkopff_immune <- subset(seurat_data, subset = cell_type %in% immune_cells)
dim(weinkopff_immune)
table(weinkopff_immune$cell_type)

#Checking the quality control values
VlnPlot(weinkopff_immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


############ Cheng Dataset - Skin Allergy ###############

##Loading the counts provided by the author
cheng_dt <- paste0(dir_data, "/Immune.rds")
cheng_immune <- readRDS(cheng_dt)

cheng_immune <- RenameIdents(cheng_immune,'cDC1' = "Dendritic Cells 1", cDC2 = "Dedritic Cells 2", "M/B" = "Mast cells", "Neu" = "Neutrophils", "NK" = "NK cells", "dγdT" = "γδ-T cells", "M/MdM" = "Mono/MonoDerivedMac", "LC" = "Langerhans Cells")
head(cheng_immune)

old_names <- c("cDC1", "cDC2", "Mac", "M/B", "Neu", "NK", "dγdT", "M/MdM", "LC")
new_names <- c("Dendritic Cells 1", "Dendritic Cells 2", "Macrophage", "Mast cells", "Neutrophils", "NK cells", "γδ-T cells", "M/MdM", "Langerhans cells")
for(i in 1:length(old_names)) {
  cheng_immune$ident <- str_replace(cheng_immune$ident, paste0("^",old_names[i],"$"), new_names[i])
}
head(cheng_immune)
table(cheng_immune$ident)

##### Map reference ######

skin.list <- list(weinkopff_immune, cheng_immune)

for (i in 1:length(skin.list)) {
  skin.list[[i]] <- NormalizeData(skin.list[[i]], verbose = FALSE)
  skin.list[[i]] <- FindVariableFeatures(skin.list[[i]], selection.method = "vst", nfeatures = 2000,
                                             verbose = FALSE)
}

skin.anchors <- FindIntegrationAnchors(object.list = skin.list, dims = 1:30)
skin.integrated <- IntegrateData(anchorset = skin.anchors, dims = 1:30)   

# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(skin.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
skin.integrated <- ScaleData(skin.integrated, verbose = FALSE)
skin.integrated <- RunPCA(skin.integrated, npcs = 30, verbose = FALSE)
skin.integrated <- RunUMAP(skin.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

### SAVING THE RESULTS ###
saveRDS(skin.integrated, file = paste0(dir_results, "/skin-integrated-dataset.rds"))

### READING THE RESULTS ###
getwd()
list.files()
seurat_integrated <- readRDS(paste0(dir_results, "/skin-integrated-dataset.rds"))
dim(seurat_integrated)

p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend()
p1 + p2
skin.integrated = seurat_integrated


##### TESTANDO REFERENCE MAP COM CLUSTER T CELL ####
load(paste0(dir_data,"/lm.combined.anno.Robj"))

immune_group <- c("CD4 and CD8 T Cells", "ILC2 and ILC3 cells", "NK NKT and γδ T cells", "Proliferating T cell", "T regs", "Eosinophils", "Basophils", "Neutrophil", "Macrophages", "Monocyte-derived tissue macrophages ", "Monocytes")
immune_data <- subset(lm.combined.anno, subset = fine_celltypes %in% immune_group)
immune_data$overall_group <- factor(x = immune_data$overall_group, levels = c("naive", "short_term_exposure", "long_term_exposure"))
dim(immune_data)
table(immune_data$fine_celltypes)

immune_umap <- RunUMAP(immune_data, dims = 1:30, verbose = FALSE)
immune_neigh <- FindNeighbors(immune_umap, dims = 1:30, verbose = FALSE)
immune_cluster <- FindClusters(immune_neigh, verbose = FALSE, resolution = 1, graph.name = "integrated_nn")

skin.anchors <- FindTransferAnchors(reference = skin.integrated, query = immune_cluster,
                                        dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = skin.anchors, refdata = skin.integrated$cell_type,
                            dims = 1:30)
immune_cluster <- AddMetaData(immune_cluster, metadata = predictions)

table(immune_cluster$predicted.id)

library(ggplot2)
library(cowplot)
library(patchwork)

skin.integrated <- RunUMAP(skin.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
immune_cluster <- MapQuery(anchorset = skin.anchors, reference = skin.integrated, query = immune_cluster,
                           refdata = list(celltype = "cell_type"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(skin.integrated, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(immune_cluster, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

DimPlot(immune_cluster, reduction = "ref.umap", group.by = "predicted.celltype", label = T, split.by = "overall_group", ncol = 2)


VlnPlot(immune_cluster, c("Cd16", "Ncam1"), group.by = "predicted.id")
