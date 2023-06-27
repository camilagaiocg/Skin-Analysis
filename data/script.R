remotes::install_github("mojaveazure/seurat-disk")

###Loading packages###
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(ggplot2)
library(sctransform)
library(dplyr)
library(patchwork)
library(SingleR)
library(stringr)
library(ggpubr)

###Getting the address
getwd()
#dir.create
dir_data <- paste0(getwd(), "/data")
dir_results<- paste0(getwd(), "/results")
dir_figures<- paste0(getwd(), "/figures")

#Loading data
load(paste0(dir_data,"/lm.combined.anno.Robj"))

#Checking the metadata
dim(lm.combined.anno)
head(lm.combined.anno@meta.data)
table(lm.combined.anno@meta.data$overall_group)
VlnPlot(lm.combined.anno, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Reorganize the object to plot the UMAP
group_name <- "umap_skin"
object_skin <- lm.combined.anno

#DimPlot separated by group
pdf(file=paste0(dir_figures, "/",group_name,"_separated",".pdf"), width = 30, height = 10)
groups <- c("naive", "short_term_exposure", "long_term_exposure")
naive_umap <- create_umaps("naive")
st_umap <- create_umaps("short_term_exposure")
lt_umap <- create_umaps("long_term_exposure")

ggarrange(naive_umap, st_umap, lt_umap + rremove("x.text"),
          ncol = 3, nrow = 1)
dev.off()

plotCellTypeProps(x = object_skin, clusters = NULL, sample = NULL)

##########################################################################################################################################
######################################### RECLUSTERING SUBSETS OF CELLS ##################################################################
##########################################################################################################################################


table(lm.combined.anno@meta.data$fine_celltypes)

##### SUBSET T CELLS ####
t_cells_group <- c("CD4 and CD8 T Cells", "ILC2 and ILC3 cells", "NK NKT and γδ T cells", "Proliferating T cell", "T regs")
t_cells_data <- subset(lm.combined.anno, subset = fine_celltypes %in% t_cells_group)
t_cells_data$overall_group <- factor(x = t_cells_data$overall_group, levels = c("naive", "short_term_exposure", "long_term_exposure"))
dim(t_cells_data)
table(t_cells_data$fine_celltypes)

t_cells_umap <- RunUMAP(t_cells_data, dims = 1:30, verbose = FALSE) 
t_cells_neigh <- FindNeighbors(t_cells_umap, dims = 1:30, verbose = FALSE) 
t_cells_cluster <- FindClusters(t_cells_neigh, verbose = FALSE, resolution = 1, graph.name = "integrated_nn")

DimPlot(t_cells_cluster,
        reduction = "umap",
        label.size = 5,
        pt.size = 0.5,
        label = T,
        split.by = 'overall_group',
        ncol = 2,
)
##### SUBSET EO/BASO/NO ####
granulocytes_group <- c("Eosinophils", "Basophils", "Neutrophil")
granulocytes_data <- subset(lm.combined.anno, subset = fine_celltypes %in% granulocytes_group)
granulocytes_data$overall_group <- factor(x = granulocytes_data$overall_group, levels = c("naive", "short_term_exposure", "long_term_exposure"))
dim(granulocytes_data)
table(granulocytes_data$fine_celltypes)

granulocytes_umap <- RunUMAP(granulocytes_data, dims = 1:30, verbose = FALSE) 
granulocytes_neigh <- FindNeighbors(granulocytes_umap, dims = 1:30, verbose = FALSE) 
granulocytes_cluster <- FindClusters(granulocytes_neigh, verbose = FALSE, resolution = 1, graph.name = "integrated_snn")

DimPlot(granulocytes_cluster,
        reduction = "umap",
        label.size = 5,
        pt.size = 0.5,
        label = T,
        split.by = 'overall_group',
        ncol = 2
)
##### MONO/MAC ####
mono_mac_group <- c("Macrophages", "Monocyte-derived tissue macrophages ", "Monocytes")
mono_mac_data <- subset(lm.combined.anno, subset = fine_celltypes %in% mono_mac_group)
mono_mac_data$overall_group <- factor(x = mono_mac_data$overall_group, levels = c("naive", "short_term_exposure", "long_term_exposure"))
dim(mono_mac_data)
table(mono_mac_data$fine_celltypes)

mono_mac_umap <- RunUMAP(mono_mac_data, dims = 1:30, verbose = FALSE) 
mono_mac_neigh <- FindNeighbors(mono_mac_umap, dims = 1:30, verbose = FALSE) 
mono_mac_cluster <- FindClusters(mono_mac_neigh, verbose = FALSE, resolution = 1, graph.name = "integrated_snn")

DimPlot(mono_mac_cluster,
        reduction = "umap",
        label.size = 5,
        pt.size = 0.5,
        label = T,
        split.by = 'overall_group',
        ncol = 2
)

##########################################################################################################################################
################################################# AUXILIAR FUNCTIONS #####################################################################
##########################################################################################################################################

create_umaps = function(i) {
  if(i == "naive")
    label <- "Naive"
  else if(i == "short_term_exposure")
    label <- "ST"
  else 
    label <- "LT"
  sub <- subset(object_skin, subset = overall_group == i)
  umap <- DimPlot(sub,
         reduction = "umap",
         label.size = 3,
         pt.size = 0.5,
         label = T,
         label.box = T,
         label.color = "white",
         ncol = 2,
         repel = T,
  ) + NoLegend()  + labs(title = label) + CenterTitle()
  return(umap)
} 
