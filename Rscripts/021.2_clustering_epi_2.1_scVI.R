# from 030_cluster_str.R
analysis_step <- "021.2_clustering_epi_2.1_scVI"

# load packages ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)
library(SeuratWrappers)
library(reticulate)
options(future.globals.maxSize = 1e10)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path_root <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path_root, res_path))

# Helper function ----
cluster_harm = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- IntegrateLayers(
    object = seu, method = HarmonyIntegration,
    orig.reduction = "pca", 
    new.reduction = "harmony")
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
  return(seu)
}

recluster_harm = function(seu_obj, npcs, res){
  seu[["RNA"]]$scale.data <- NULL
  seu[["RNA"]]$data <- NULL
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- IntegrateLayers(
    object = seu, method = HarmonyIntegration,
    orig.reduction = "pca", 
    new.reduction = "harmony")
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
  return(seu)
}

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Load data ----
seu <- readRDS("RDSfiles/seu_020_epi.RDS")
# load("RDSfiles/epi_type_names.Rdata")

# ITER1: cluster w/ harmony ----
plot_path <- file.path(plot_path_root, "iter1")
fs::dir_create(c(plot_path))
seu <- subset(seu, subset = celltype %in% c("Pariet.", "EEC", "Gob."), invert = T)
# seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident) 
seu <- recluster_harm(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

# Check markers by FindAllMarkers
seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  # dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 1000) %>%
  ungroup() -> markers
# openxlsx2::write_xlsx(markers, file.path(res_path, "markers.xlsx"))

### c21, c26: stromal contamination
### plasma cell contamination is difficult to resolve... > leave as it is
### >>> decided to remove c21, c26

# ITER2: remove probable contamination, then re-cluster ----
plot_path <- file.path(plot_path_root, "iter2")
fs::dir_create(c(plot_path))
seu <- subset(seu, idents = c(21,26), invert = TRUE)
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
seu <- recluster_harm(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

DimPlot(seu, group.by = "celltype") + NoAxes()
ggsave("tissue_type.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

### not well clustered by marker gene expression...
### even better before with Priet., EEC & Gob.
### >>> try scVI integration









# Add cellgroup annotation ----
ArtEC_names <- colnames(seu)[seu$seurat_clusters %in% c(7,9)]
VenEC_names <- colnames(seu)[seu$seurat_clusters %in% c(1,5,10,14)]
CapEC_names <- colnames(seu)[seu$seurat_clusters %in% c(0,2:4,6,8,12,13,15,16)]
ProlifEC_names <- colnames(seu)[seu$seurat_clusters %in% c(11)]
save(ArtEC_names, VenEC_names, CapEC_names, ProlifEC_names, file = file.path("RDSfiles", "endo_type_names.Rdata"))

seu$celltype2 <- ""
seu$celltype2[ArtEC_names] <- "ArtEC"
seu$celltype2[VenEC_names] <- "VenEC"
seu$celltype2[CapEC_names] <- "CapEC"
seu$celltype2[ProlifEC_names] <- "ProlifEC"
seu$celltype2 <- factor(seu$celltype2, levels = c("ArtEC", "VenEC", "CapEC", "ProlifEC"))
DimPlot(seu, group.by = "celltype2", cols = "polychrome") & NoAxes()
ggsave("celltype2.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Dot plot
Idents(seu) <- "celltype2"
markers <- FindAllMarkers(seu, only.pos = TRUE)
openxlsx2::write_xlsx(markers, file.path(res_path, paste0("markers_endo_type.xlsx")))
features <- c("SEMA3G", "GJA4", "ACKR1", "SELP", "RGCC", "CA4", "MKI67", "TOP2A")
DotPlot(seu, group.by = "celltype2", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5.5, height = 4, units = "in", dpi = 150)

# save RDS
saveRDS(seu, file = file.path("RDSfiles", "seu_032_endo.RDS"))

# additional plots ----
# to compare tumor vs non-tumor etc.
# exclude tissue_type peri_NT (only 2 cells)
seu <- subset(seu, subset = tissue_type == "peri_NT", invert = T)
seu$tissue_type <- droplevels(seu$tissue_type)
seu$tissue_type <- fct_recode(seu$tissue_type, "peritoneal"="peri_T", "non_tum"="pri_NT", "primary"="pri_T")
seu$tissue_type <- fct_relevel(seu$tissue_type, "primary")

DimPlot(seu, group.by = "tissue_type") + NoAxes()
ggsave("tissue_type.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

add_feat <- "ENG"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            max.cutoff = "q80", min.cutoff = "q0"
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

DotPlot(seu, features = "LRG1", group.by = "Laurens")
ggsave("LRG1_dotplot.png", path = plot_path, width = 4, height = 4, units = "in", dpi = 150)

# VlnPlot(seu, features = "LRG1", group.by = "tissue_type", pt.size = 1) & NoLegend()
# ggsave("LRG1_vln.png", path = plot_path, width = 3, height = 3, units = "in", dpi = 150)
