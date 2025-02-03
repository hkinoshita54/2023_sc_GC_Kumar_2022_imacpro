# 
analysis_step <- "030_clustering_str"

# load packages ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path_root <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

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
seu_all <- readRDS("RDSfiles/seu_010_raw.RDS")
load("RDSfiles/cellgroup_names.Rdata")

# ITER1: cluster w/ harmony ----
plot_path <- file.path(plot_path_root, "iter1")
fs::dir_create(c(plot_path))
seu <- seu_all[, str_names]
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident) 
seu <- cluster_harm(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

# ITER2: remove probable heterotypic multiplets, then re-cluster ----
plot_path <- file.path(plot_path_root, "iter2")
fs::dir_create(c(plot_path))
seu <- subset(seu, idents = c(10,11,17,21), invert = TRUE)
seu <- recluster_harm(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

# ITER3: remove probable heterotypic multiplets, then re-cluster ----
plot_path <- file.path(plot_path_root, "iter3")
fs::dir_create(c(plot_path))
seu <- subset(seu, idents = c(13,22), invert = TRUE)
seu <- recluster_harm(seu, npcs = 50, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 4.5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

# Add cellgroup annotation ----
fib_names <- colnames(seu)[seu$seurat_clusters %in% c(0,2,5,6,7,9,13,17,18,21)]
endo_names <- colnames(seu)[seu$seurat_clusters %in% c(1,3,8,10,11,16,19)]
myo_names <- colnames(seu)[seu$seurat_clusters %in% c(4,12,14,15,22)]
glia_names <- colnames(seu)[seu$seurat_clusters %in% c(20)]
save(fib_names, endo_names, myo_names, glia_names, file = file.path("RDSfiles", "str_type_names.Rdata"))

seu$celltype <- ""
seu$celltype[fib_names] <- "Fib."
seu$celltype[endo_names] <- "Endo."
seu$celltype[myo_names] <- "Myo."
seu$celltype[glia_names] <- "Glia"
seu$celltype <- factor(seu$celltype, levels = c("Fib.", "Endo.", "Myo.", "Glia"))
DimPlot(seu, group.by = "celltype") & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# save RDS
saveRDS(seu, file = file.path("RDSfiles", "seu_030_str.RDS"))

# additional plots ----
add_feat <- "ENG"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            max.cutoff = "q80", min.cutoff = "q0"
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

