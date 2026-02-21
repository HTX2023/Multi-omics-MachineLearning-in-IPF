# ==================================================================================================
# Complete WGCNA workflow (teaching-style, end-to-end)
#
# Prerequisites (objects must exist in the environment):
#   - exp      : expression matrix/data.frame (genes x samples) OR numeric matrix (genes x samples)
#   - clinical : sample metadata data.frame with rownames = sample IDs and column "Group"
#               where Group contains at least: "IPF" and "Normal"
#
# Main steps:
#   1) Prepare expression matrix (samples x genes) using top MAD genes
#   2) Remove low-quality samples/genes and detect outliers by sample clustering
#   3) Build trait matrix and align traits to expression samples
#   4) Choose soft-thresholding power
#   5) Construct network and identify modules (blockwiseModules)
#   6) Module–trait correlation heatmap
#   7) GS vs MM scatterplot for a chosen module and chosen trait
#   8) TOM heatmap on a random subset of genes
#   9) Eigengene network plots
#  10) Export module gene lists to Excel
#
# Outputs (files):
#   - 2.sampleClustering.png
#   - 4.Soft_threshold.png
#   - 5.DendroAndColors.png
#   - genes.RData
#   - 6.labeledHeatmap.png
#   - 7.MM_GS_scatterplot.png
#   - 8.Sub_netheatmap.png
#   - 9.Eigengene_dendrogram.png
#   - 10.Eigengene_heatmap.png
#   - WGCNA_module_genes.xlsx
# ==================================================================================================
rm(list = ls())
options(stringsAsFactors = FALSE, check.names = FALSE)
suppressPackageStartupMessages({
  library(WGCNA)
  library(tinyarray)
  library(stringr)
  library(gplots)
  library(openxlsx)
})
enableWGCNAThreads()
stopifnot(exists("exp"), exists("clinical"))
if (!"Group" %in% colnames(clinical)) stop("clinical must contain a column named 'Group'.")
exp_num <- as.matrix(exp)
storage.mode(exp_num) <- "numeric"
if (anyNA(exp_num)) warning("Expression matrix contains NA values; WGCNA will attempt to filter bad genes/samples.")
n_keep <- 5000
if (nrow(exp_num) < n_keep) n_keep <- nrow(exp_num)
mad_rank <- order(apply(exp_num, 1, mad, na.rm = TRUE), decreasing = TRUE)
datExpr0 <- t(exp_num[mad_rank[seq_len(n_keep)], , drop = FALSE])
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
sampleTree <- hclust(dist(datExpr0), method = "average")
png("2.sampleClustering.png", width = 5000, height = 2000, res = 300)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
dev.off()
cutHeight <- 130
clust <- cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
keepSamples <- (clust != 0)
datExpr0 <- datExpr0[keepSamples, , drop = FALSE]
datExpr <- datExpr0
traitData <- data.frame(
  IPF = as.numeric(clinical$Group == "IPF"),
  CTL = as.numeric(clinical$Group == "Normal"),
  stringsAsFactors = FALSE
)
rownames(traitData) <- rownames(clinical)
datTraits <- traitData
sampleNames <- rownames(datExpr)
datTraits <- datTraits[sampleNames, , drop = FALSE]
traitColors <- numbers2colors(datTraits, signed = FALSE)
sampleTree2 <- hclust(dist(datExpr), method = "average")
png("3.sampleDendro_traitHeatmap.png", width = 5000, height = 2200, res = 300)
plotDendroAndColors(
  sampleTree2,
  traitColors,
  groupLabels = names(datTraits),
  main = "Sample dendrogram and trait heatmap"
)
dev.off()
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 12
cex1 <- 0.9
png("4.Soft_threshold.png", width = 3000, height = 1500, res = 300)
par(mfrow = c(1, 2))
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = cex1, col = "red")
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "testTOM",
  verbose = 3
)
print(table(net$colors))
mergedColors <- labels2colors(net$colors)
png("5.DendroAndColors.png", width = 2000, height = 1200, res = 300)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
gm <- data.frame(net.colors = net$colors, color = moduleColors, stringsAsFactors = FALSE)
genes <- split(rownames(gm), gm$color)
save(genes, file = "genes.RData")
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
png("6.labeledHeatmap.png", width = 2000, height = 2000, res = 300)
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module–Trait Relationships"
)
dev.off()
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)
i_trait <- 1
module_to_plot <- "tan"
instrait <- datTraits[, i_trait]
geneTraitSignificance <- as.data.frame(cor(datExpr, instrait, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste0("GS.", names(datTraits)[i_trait])
names(GSPvalue) <- paste0("p.GS.", names(datTraits)[i_trait])
png("7.MM_GS_scatterplot.png", width = 2000, height = 2000, res = 300)
column <- match(module_to_plot, modNames)
moduleGenesFlag <- (moduleColors == module_to_plot)
verboseScatterplot(
  abs(geneModuleMembership[moduleGenesFlag, column]),
  abs(geneTraitSignificance[moduleGenesFlag, 1]),
  xlab = paste("Module Membership in", module_to_plot, "module"),
  ylab = "Gene significance",
  main = "Module membership vs. gene significance",
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = module_to_plot
)
dev.off()
nSelect <- 4000
if (nGenes < nSelect) nSelect <- nGenes
set.seed(10)
dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = softPower)
select_idx <- sample(nGenes, size = nSelect)
selectTOM <- dissTOM[select_idx, select_idx]
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select_idx]
myheatcol <- colorpanel(250, "red", "orange", "lemonchiffon")
png("8.Sub_netheatmap.png", width = 2000, height = 2000, res = 300)
plotDiss <- selectTOM^7
diag(plotDiss) <- NA
TOMplot(
  plotDiss,
  selectTree,
  selectColors,
  col = myheatcol,
  main = "Network heatmap plot (selected genes)"
)
dev.off()
MET <- orderMEs(cbind(MEs, instrait))
png("9.Eigengene_dendrogram.png", width = 2000, height = 2000, res = 300)
plotEigengeneNetworks(
  MET,
  "Eigengene dendrogram",
  marDendro = c(0, 4, 2, 0),
  plotHeatmaps = FALSE
)
dev.off()
png("10.Eigengene_heatmap.png", width = 2000, height = 2000, res = 300)
plotEigengeneNetworks(
  MET,
  "Eigengene adjacency heatmap",
  plotDendrograms = FALSE,
  marHeatmap = c(4, 5, 2, 2),
  xLabelsAngle = 90
)
dev.off()
moduleGenes <- lapply(names(genes), function(mod) data.frame(Gene = genes[[mod]], stringsAsFactors = FALSE))
names(moduleGenes) <- names(genes)
maxLen <- max(vapply(moduleGenes, nrow, numeric(1)))
geneTable <- as.data.frame(matrix(NA_character_, nrow = maxLen, ncol = length(moduleGenes)))
colnames(geneTable) <- names(moduleGenes)
for (mod in names(moduleGenes)) {
  geneList <- moduleGenes[[mod]]$Gene
  geneTable[[mod]] <- c(geneList, rep(NA_character_, maxLen - length(geneList)))
}
write.xlsx(geneTable, file = "WGCNA_module_genes.xlsx", rowNames = FALSE)
cat("Module gene lists have been saved to 'WGCNA_module_genes.xlsx'\n")