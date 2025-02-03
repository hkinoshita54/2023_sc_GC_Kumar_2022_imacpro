# 
analysis_step <- "020_clustering_epi"

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
seu <- seu_all[, epi_names]
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
seu <- subset(seu, idents = c(1), invert = TRUE)
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

DimPlot(seu, group.by = "tissue_type") + NoAxes()
ggsave("tissue_type.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# Add cellgroup annotation ----
prog_names <- colnames(seu)[seu$seurat_clusters %in% c(7,10,12)]
pit_names <- colnames(seu)[seu$seurat_clusters %in% c(1,3,11)]
neck_names <- colnames(seu)[seu$seurat_clusters %in% c(5,14)]
chief_names <- colnames(seu)[seu$seurat_clusters %in% c(4,13,15,19,27)]
pariet_names <- colnames(seu)[seu$seurat_clusters %in% c(18)]
eec_names <- colnames(seu)[seu$seurat_clusters %in% c(6,23)]
gob_names <- colnames(seu)[seu$seurat_clusters %in% c(24)]
epi_nos_names <- colnames(seu)[seu$seurat_clusters %in% c(0,2,8,9,16,17,20,21,22,25,26)]
save(prog_names, pit_names, neck_names, chief_names, pariet_names, eec_names, gob_names, epi_nos_names, 
     file = file.path("RDSfiles", "epi_type_names.Rdata"))

seu$celltype <- ""
seu$celltype[prog_names] <- "Prog."
seu$celltype[pit_names] <- "Pit"
seu$celltype[neck_names] <- "Neck"
seu$celltype[chief_names] <- "Chief"
seu$celltype[pariet_names] <- "Pariet."
seu$celltype[eec_names] <- "EEC"
seu$celltype[gob_names] <- "Gob."
seu$celltype[epi_nos_names] <- "EpiNOS"
seu$celltype <- factor(seu$celltype, levels = c("Prog.", "Pit", "Neck", "Chief", "Pariet.", "EEC", "Gob.", "EpiNOS"))
DimPlot(seu, group.by = "celltype", cols = "polychrome") & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# save RDS
saveRDS(seu, file = file.path("RDSfiles", "seu_020_epi.RDS"))

# additional plots ----
add_feat <- "MSLN"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            max.cutoff = "q80", min.cutoff = "q0"
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

