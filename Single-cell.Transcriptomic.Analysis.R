# ==================================================================================================
# scRNA-seq workflow (GSE122960):
#   1) Environment setup + package loading
#   2) Untar raw archive and read 10X HDF5 (.h5) files
#   3) Build Seurat objects per sample and merge into one object
#   4) Add QC metrics, visualize QC, and filter low-quality cells
#   5) Normalize, HVG selection, scaling, PCA, Harmony integration
#   6) Graph-based clustering at multiple resolutions and UMAP/t-SNE embedding
#   7) SingleR-based coarse annotation using HumanPrimaryCellAtlas (pruned labels)
#   8) Marker detection (FindAllMarkers) based on SingleR identities
#   9) Visualization: cell count barplot, cell proportion barplot, violin plot, dot plot
#
# Inputs:
#   - GSE122960_RAW.tar (contains multiple *.h5 files)
# Outputs (examples):
#   - obj.Rdata
#   - markers.RData
#   - sce.all.annotated.RData
# ==================================================================================================
rm(list = ls())
options(stringsAsFactors = FALSE, check.names = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(dplyr)
  library(future)
  library(clustree)
  library(cowplot)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(hdf5r)
  library(ggsci)
})
untar("GSE122960_RAW.tar", exdir = "GSE122960_RAW")
raw_dir <- "GSE122960_RAW/"
samples <- list.files(raw_dir, pattern = "\\.h5$", full.names = FALSE)
if (length(samples) == 0) stop("No .h5 files found under: ", raw_dir)
sample_info <- data.frame(
  filename = samples,
  sample_type = ifelse(grepl("IPF", samples), "IPF", "Donor"),
  sample_id = gsub(".*_(IPF|Donor)_(\\d+).*", "\\1_\\2", samples),
  stringsAsFactors = FALSE
)
sce_list <- lapply(samples, function(fname) {
  fpath <- gsub(" ", "", file.path(raw_dir, fname))
  tmp <- Seurat::Read10X_h5(fpath)
  if (is.list(tmp)) {
    gene_key <- grep("Gene", names(tmp), value = TRUE)[1]
    if (is.na(gene_key)) gene_key <- names(tmp)[1]
    counts <- tmp[[gene_key]]
  } else {
    counts <- tmp
  }
  CreateSeuratObject(
    counts = counts,
    project = sample_info$sample_id[match(fname, sample_info$filename)],
    min.cells = 5,
    min.features = 300
  )
})
sce_all <- merge(
  x = sce_list[[1]],
  y = sce_list[-1],
  add.cell.ids = sample_info$sample_id
)
if (packageVersion("Seurat") >= "5.0.0") {
  sce_all <- JoinLayers(sce_all)
  counts_layer <- LayerData(sce_all, assay = "RNA", layer = "counts")
} else {
  counts_layer <- GetAssayData(sce_all, slot = "counts")
}
print(paste0("Total cells (pre-QC): ", ncol(sce_all)))
print(table(sce_all$orig.ident))
print(head(sce_all@meta.data, 3))
sce_all[["percent.mt"]] <- PercentageFeatureSet(sce_all, pattern = "^MT-")
sce_all[["percent.rp"]] <- PercentageFeatureSet(sce_all, pattern = "^RP[SL]")
sce_all[["percent.hb"]] <- PercentageFeatureSet(sce_all, pattern = "^HB[^(P)]")
print(head(sce_all@meta.data, 3))
print(dim(sce_all))
VlnPlot(
  sce_all,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp", "percent.hb"),
  ncol = 3,
  pt.size = 0,
  group.by = "orig.ident"
)
sce_all <- subset(
  sce_all,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5000 &
    nCount_RNA < 20000 &
    percent.mt < 20 &
    percent.hb < 1
)
print(dim(sce_all))
print(paste0("Total cells (post-QC): ", ncol(sce_all)))
p_sc1 <- FeatureScatter(sce_all, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
p_sc2 <- FeatureScatter(sce_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
print(p_sc1 + p_sc2)
obj_file <- "obj.Rdata"
if (!file.exists(obj_file)) {
  sce_all <- sce_all %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(features = rownames(.)) %>%
    RunPCA(pc.genes = VariableFeatures(.)) %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:20, reduction = "harmony") %>%
    FindClusters(graph.name = "RNA_snn", resolution = c(0.05, 0.1, 0.2, 0.5)) %>%
    RunUMAP(dims = 1:20, reduction = "harmony") %>%
    RunTSNE(dims = 1:20, reduction = "harmony")
  save(sce_all, file = obj_file)
}
load(obj_file)
ElbowPlot(sce_all)
UMAPPlot(sce_all, label = TRUE)
TSNEPlot(sce_all, label = TRUE)
suppressPackageStartupMessages({
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(SingleCellExperiment)
})
sce_all <- JoinLayers(sce_all, assay = "RNA")
DefaultAssay(sce_all) <- "RNA"
print(head(colnames(sce_all@meta.data), 20))
print(DimPlot(sce_all, reduction = "tsne", group.by = "RNA_snn_res.0.05", label = TRUE) + ggtitle("res.0.05"))
print(DimPlot(sce_all, reduction = "tsne", group.by = "RNA_snn_res.0.1",  label = TRUE) + ggtitle("res.0.1"))
print(DimPlot(sce_all, reduction = "tsne", group.by = "RNA_snn_res.0.2",  label = TRUE) + ggtitle("res.0.2"))
print(DimPlot(sce_all, reduction = "tsne", group.by = "RNA_snn_res.0.2",  label = TRUE) + ggtitle("Single cell clustering (res.0.2)"))
keep_types <- c(
  "T_cells", "NK_cell", "Macrophage", "Monocyte", "Epithelial_cells",
  "Endothelial_cells", "DC", "B_cell", "Neutrophils", "Fibroblasts"
)
bpp <- MulticoreParam(workers = 4)
sc_sce <- as.SingleCellExperiment(sce_all, assay = "RNA")
full_ref <- HumanPrimaryCellAtlasData()
ref_idx <- which(full_ref$label.main %in% keep_types)
ref_sub <- full_ref[, ref_idx]
singleR_pred <- SingleR(
  test = sc_sce,
  ref = ref_sub,
  labels = ref_sub$label.main,
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  BPPARAM = bpp
)
sce_all$SR_sub_pruned <- singleR_pred$pruned.labels
sce_all$SR_sub_pruned2 <- ifelse(
  sce_all$SR_sub_pruned %in% keep_types,
  sce_all$SR_sub_pruned,
  "Others"
)
sce_all$SR_sub_pruned2 <- factor(sce_all$SR_sub_pruned2, levels = c(keep_types, "Others"))
print(DimPlot(sce_all, reduction = "umap", group.by = "SR_sub_pruned2", label = TRUE) + ggtitle("Single cell clustering and annotation"))
print(DimPlot(sce_all, reduction = "tsne", group.by = "SR_sub_pruned2", label = TRUE) + ggtitle("Single cell clustering and annotation"))
print(table(sce_all$SR_sub_pruned2))
save(sce_all, file = "sce.all.annotated.RData")
scRNA <- sce_all
scRNA$SR_sub_pruned2 <- factor(scRNA$SR_sub_pruned2)
Idents(scRNA) <- "SR_sub_pruned2"
markers <- FindAllMarkers(
  object = scRNA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
allmarkers <- markers
save(allmarkers, file = "markers.RData")
dat_count <- as.data.frame(table(Idents(scRNA)))
dat_count$label <- paste(dat_count$Var1, dat_count$Freq, sep = ":")
suppressPackageStartupMessages({
  library(paletteer)
  library(forcats)
})
dat_count$Var1 <- fct_reorder(dat_count$Var1, dat_count$Freq)
p_count <- ggplot(dat_count, aes(x = Freq, y = Var1, fill = Var1)) +
  scale_fill_paletteer_d("ggsci::category20_d3") +
  geom_bar(stat = "identity") +
  geom_text(aes(x = 0, label = label), hjust = 0) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")
  )
print(p_count)
p_count_axis <- ggplot(dat_count, aes(x = Freq, y = Var1, fill = Var1)) +
  scale_fill_paletteer_d("ggsci::category20_d3") +
  geom_bar(stat = "identity") +
  geom_text(aes(x = 0, label = label), hjust = 0) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x.top = element_line(colour = "white"),
    axis.line.y.right = element_line(colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 1, 0.5, 1), "cm")
  )
print(p_count_axis)
scRNA$SR_sub_pruned2 <- factor(
  scRNA$SR_sub_pruned2,
  levels = c(
    "T_cells", "NK_cell", "Macrophage", "Monocyte", "Epithelial_cells",
    "Endothelial_cells", "DC", "B_cell", "Neutrophils", "Fibroblasts",
    "Others"
  )
)
p_prop <- ggplot(scRNA@meta.data, aes(x = orig.ident, fill = SR_sub_pruned2)) +
  geom_bar(position = "fill", alpha = 0.9, width = 0.5) +
  scale_fill_paletteer_d("ggsci::category20_d3") +
  theme_classic() +
  coord_fixed(ratio = 4) +
  coord_flip() +
  labs(x = "Sample (orig.ident)", y = "Cell Proportion", fill = "Cell Type")
print(p_prop)
violin_genes <- c(
  "CD3D", "CD4", "CD8A",
  "NCAM1", "FCGR3A", "NKG7",
  "CD68", "ADGRE1", "CD163",
  "CD14", "CX3CR1", "CCR2",
  "EPCAM", "KRT18", "KRT19",
  "PECAM1", "VWF", "CDH5",
  "CD1C", "CLEC9A",
  "CD19", "MS4A1", "CD79A",
  "CSF3R", "S100A8", "CXCL3",
  "PDGFRA", "DCN", "COL1A1"
)
dup_genes <- violin_genes[duplicated(violin_genes)]
print(dup_genes)
scRNA[["RNA"]] <- JoinLayers(scRNA[["RNA"]])
expr_mat <- as.matrix(GetAssayData(object = scRNA, assay = "RNA", layer = "data"))
vln_df <- expr_mat %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(violin_genes)) %>%
  rownames_to_column("cell") %>%
  mutate(cluster = scRNA$SR_sub_pruned2) %>%
  pivot_longer(cols = -c(cell, cluster), names_to = "gene", values_to = "expr") %>%
  mutate(gene = factor(gene, levels = violin_genes))
base_cols <- paletteer_d("ggsci::category20_d3")
my_cols <- colorRampPalette(base_cols)(length(unique(vln_df$cluster)))
p_vln <- ggplot(vln_df, aes(x = expr, y = cluster, fill = cluster)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = my_cols) +
  facet_grid(. ~ gene, scales = "free_y", switch = "x") +
  scale_x_continuous(expand = c(0, 0), position = "top") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x.top = element_blank(),
    axis.text.x.top = element_text(hjust = 1, size = 7),
    axis.text.y = element_text(color = "black", size = 14),
    axis.ticks.y = element_line(color = "black"),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0, "cm"),
    strip.text.y = element_text(angle = 0, size = 14, hjust = 0),
    strip.background.y = element_blank()
  )
print(p_vln)
load("markers.RData")
Idents(scRNA) <- "SR_sub_pruned2"
dot_genes <- c(
  "CD3D", "CD4", "CD8A",
  "NCAM1", "FCGR3A", "NKG7",
  "CD68", "ADGRE1", "CD163",
  "CD14", "CX3CR1", "CCR2",
  "EPCAM", "KRT18", "KRT19",
  "PECAM1", "VWF", "CDH5",
  "CD1C", "CLEC9A",
  "CD19", "MS4A1", "CD79A",
  "CSF3R", "S100A8", "CXCL3",
  "PDGFRA", "DCN", "COL1A1"
)
p_dot <- DotPlot(
  object = scRNA,
  features = dot_genes,
  cols = "PRGn",
  group.by = "SR_sub_pruned2"
) +
  RotatedAxis() +
  labs(
    x = "Genes",
    y = "Cell Type",
    color = "Avg Expression",
    size = "Percent Expressing"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_dot)