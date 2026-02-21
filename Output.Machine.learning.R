# ==================================================================================================
# Machine-learning pipeline (binary classification) with core-gene intersection, model training,
# AUC benchmarking, heatmap visualization, cohort-wise ROC curves, and gene-wise ROC curves.
# Requirements:
#   1) pipline.Machine.Learning.R must define: scaleData(), RunML(), ExtractVar(), CalPredictScore(),
#      PredictClass(), RunEval()
#   2) refer.ML.R must define: SimpleHeatmap()
# Inputs:
#   - train_merge_data.csv / test_merge_data.csv (samples as rows; features as columns; last column = Group)
#   - 113_ML_methods.txt (column "x": model combinations, e.g., "Stepglm[both]+NaiveBayes")
# Outputs:
#   - gene_intersection_upset.pdf, core_genes.txt
#   - model.MLmodel.rds, model.riskMatrix.txt, model.classMatrix.txt, model.genes.txt, model.AUCmatrix.txt
#   - ML_AUC_heatmap.pdf, ROC.<Cohort>.pdf, ROC.genes.pdf
# ==================================================================================================
source("pipline.Machine.Learning.R")
options(stringsAsFactors = FALSE, check.names = FALSE)
suppressPackageStartupMessages({
  library(UpSetR)
  library(pROC)
  library(glmnet)
  library(RColorBrewer)
})
set.seed(123)
# =======================================
# 1) Load train/test data and split labels
# =======================================
train_df <- read.table("train_merge_data.csv", header = TRUE, row.names = 1, sep = ",")
test_df  <- read.table("test_merge_data.csv", header = TRUE, row.names = 1, sep = ",")
train_features_df <- train_df[, -ncol(train_df), drop = FALSE]
train_labels_df   <- train_df[,  ncol(train_df), drop = FALSE]
test_features_df  <- test_df[,  -ncol(test_df),  drop = FALSE]
test_labels_df    <- test_df[,   ncol(test_df),  drop = FALSE]
test_labels_df$Cohort <- gsub("(.+)\\_(.+)\\_(.+)", "\\1", rownames(test_df))
# ===============================================================
# 2) Harmonize features across cohorts and standardize expression
# ===============================================================
common_features <- intersect(colnames(train_features_df), colnames(test_features_df))
x_train <- as.matrix(train_df[, common_features, drop = FALSE])
x_test  <- as.matrix(test_df[,  common_features, drop = FALSE])
x_train <- scaleData(x_train, centerFlags = TRUE, scaleFlags = TRUE)
x_test  <- scaleData(x_test, cohort = test_labels_df$Cohort, centerFlags = TRUE, scaleFlags = TRUE)
# ============================================
# 3) Variable pre-selection (multiple methods)
# ============================================
methods_tbl <- read.table("113_ML_methods.txt", header = TRUE, sep = "\t", check.names = FALSE)
methods_all <- methods_tbl$x
class_var <- "Group"
pretrain_methods <- c("Lasso", "glmBoost", "RF", "Stepglm[both]", "Stepglm[backward]")
pretrain_vars <- list()
time_list <- list()
for (m in pretrain_methods) {
  t_used <- system.time({
    pretrain_vars[[m]] <- RunML(method = m, Train_set = x_train, Train_label = train_labels_df, mode = "Variable", classVar = class_var)
  })
  time_list[[m]] <- t_used[["elapsed"]]
  cat(sprintf("Method [%s] took %.2f seconds\n", m, t_used[["elapsed"]]))
}
pretrain_vars[["simple"]] <- colnames(x_train)
# =====================================
# 4) UpSet plot + core-gene intersection
# =====================================
gene_lists <- list(
  Lasso     = pretrain_vars$Lasso,
  glmBoost  = pretrain_vars$glmBoost,
  RF        = pretrain_vars$RF,
  Step_both = pretrain_vars$`Stepglm[both]`,
  Step_back = pretrain_vars$`Stepglm[backward]`
)
upset_data <- fromList(gene_lists)
set_colors <- c("#DF0A1F", "#1C5BA7", "#019E73", "#ED621B", "#E477C1")
pdf("gene_intersection_upset.pdf", width = 10, height = 6)
upset(upset_data,
      sets = names(gene_lists),
      order.by = "freq",
      text.scale = 1.2,
      matrix.color = set_colors,
      mainbar.y.label = "Gene number intersected",
      sets.x.label = "Gene number selected")
dev.off()
core_genes <- Reduce(intersect, gene_lists)
if (length(core_genes) == 0) {
  warning("Core-gene intersection is empty. Falling back to all common features.")
  core_genes <- colnames(x_train)
}
write.table(core_genes, "core_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
variable_core <- core_genes
min_selected_var <- length(core_genes)
# ===========================================
# 5) Train models (using core genes as inputs)
# ===========================================
models <- list()
x_train_bk <- x_train
for (method_name in methods_all) {
  cat(match(method_name, methods_all), ":", method_name, "\n")
  parts <- strsplit(method_name, "\\+")[[1]]
  if (length(parts) == 1) parts <- c("simple", parts)
  x_train_sub <- x_train_bk[, variable_core, drop = FALSE]
  fit <- RunML(method = parts[2], Train_set = x_train_sub, Train_label = train_labels_df, mode = "Model", classVar = class_var)
  if (length(ExtractVar(fit)) < min_selected_var) {
    fit <- NULL
  }
  models[[method_name]] <- fit
}
x_train <- x_train_bk
rm(x_train_bk)
models <- models[!vapply(models, is.null, logical(1))]
saveRDS(models, "model.MLmodel.rds")
# ===========================================================
# 6) Predict risk scores and classes on (train + test) samples
# ===========================================================
models <- readRDS("model.MLmodel.rds")
valid_methods <- names(models)
combined_x <- rbind.data.frame(x_train[, variable_core, drop = FALSE], x_test[, variable_core, drop = FALSE])
rs_list <- list()
for (m in valid_methods) {
  rs_list[[m]] <- CalPredictScore(fit = models[[m]], new_data = combined_x)
}
risk_mat <- as.data.frame(t(do.call(rbind, rs_list)))
risk_tab <- cbind(id = rownames(risk_mat), risk_mat)
write.table(risk_tab, "model.riskMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
class_list <- list()
for (m in valid_methods) {
  class_list[[m]] <- PredictClass(fit = models[[m]], new_data = combined_x)
}
class_mat <- as.data.frame(t(do.call(rbind, class_list)))
class_tab <- cbind(id = rownames(class_mat), class_mat)
write.table(class_tab, "model.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# ==========================================
# 7) Export selected features for each model
# ==========================================
fea_list <- list()
for (m in valid_methods) {
  fea_list[[m]] <- ExtractVar(models[[m]])
}
fea_df <- do.call(rbind, lapply(names(models), function(nm) {
  data.frame(features = ExtractVar(models[[nm]]), algorithm = nm, stringsAsFactors = FALSE)
}))
write.table(fea_df, file = "model.genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# ===========================================
# 8) Evaluate model AUC (train vs each cohort)
# ===========================================
auc_list <- list()
for (m in valid_methods) {
  auc_list[[m]] <- RunEval(fit = models[[m]],
                           Test_set = x_test[, variable_core, drop = FALSE],
                           Test_label = test_labels_df,
                           Train_set = x_train[, variable_core, drop = FALSE],
                           Train_label = train_labels_df,
                           Train_name = "Train",
                           cohortVar = "Cohort",
                           classVar = class_var)
}
auc_mat <- do.call(rbind, auc_list)
auc_tab <- cbind(Method = rownames(auc_mat), auc_mat)
write.table(auc_tab, "model.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# ============================
# 9) Plot AUC heatmap (ranked)
# ============================
auc_mat2 <- read.table("model.AUCmatrix.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
avg_auc <- apply(auc_mat2, 1, mean)
avg_auc <- sort(avg_auc, decreasing = TRUE)
auc_mat2 <- auc_mat2[names(avg_auc), , drop = FALSE]
avg_auc_fmt <- as.numeric(format(avg_auc, digits = 3, nsmall = 3))
if (ncol(auc_mat2) <= 12) {
  cohort_col <- brewer.pal(n = ncol(auc_mat2), name = "Set3")
} else {
  cohort_col <- colorRampPalette(brewer.pal(12, "Set3"))(ncol(auc_mat2))
}
names(cohort_col) <- colnames(auc_mat2)
source("refer.ML.R")
cellwidth <- 1
cellheight <- 0.5
hm <- SimpleHeatmap(Cindex_mat = auc_mat2,
                    avg_Cindex = avg_auc_fmt,
                    CohortCol = cohort_col,
                    barCol = "#c51b8a",
                    cellwidth = cellwidth,
                    cellheight = cellheight,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE)
pdf(file = "ML_AUC_heatmap.pdf", width = cellwidth * ncol(auc_mat2) + 6, height = cellheight * nrow(auc_mat2) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
# ===========================================================
# 10) Cohort-wise ROC curves for a selected best/target method
# ===========================================================
rs_file <- "model.riskMatrix.txt"
method_to_plot <- "Stepglm[both]+NaiveBayes"
risk_rt <- read.table(rs_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
cohort_id <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(risk_rt))
cohort_id <- gsub("(.*)\\.(.*)", "\\1", cohort_id)
risk_rt$Cohort <- cohort_id
for (coh in unique(risk_rt$Cohort)) {
  rt_sub <- risk_rt[risk_rt$Cohort == coh, , drop = FALSE]
  y_name <- gsub(".*_(.*)$", "\\1", rownames(rt_sub))
  y_bin <- ifelse(y_name == "Normal", 0, 1)
  if (length(unique(y_bin)) != 2) {
    cat("Error: y must contain exactly two unique values. Check data for Cohort:", coh, "\n")
    next
  }
  if (!method_to_plot %in% colnames(rt_sub)) {
    cat("Error: method column not found:", method_to_plot, "in Cohort:", coh, "\n")
    next
  }
  roc1 <- roc(y_bin, as.numeric(rt_sub[, method_to_plot]))
  ci1 <- ci.auc(roc1, method = "bootstrap")
  ci_vec <- as.numeric(ci1)
  pdf(file = paste0("ROC.", coh, ".pdf"), width = 5, height = 4.75)
  plot(roc1, print.auc = TRUE, col = "#756bb1", legacy.axes = TRUE, main = coh)
  text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ci_vec[1]), "-", sprintf("%.03f", ci_vec[3])), col = "#756bb1")
  dev.off()
}
# ==========================================================
# 11) Gene-wise ROC curves in training set (overlayed curves)
#     Note: This block uses the same label definition as the
#     training label column whenever available (more robust).
# ==========================================================
gene_file <- "model.genes.txt"
gene_rt <- read.table(gene_file, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
genes_to_plot <- gene_rt$features[gene_rt$algorithm == method_to_plot]
genes_to_plot <- unique(genes_to_plot)
if (length(genes_to_plot) == 0) {
  stop("No genes found for the specified algorithm in model.genes.txt: ", method_to_plot)
}
expr_train <- t(as.matrix(train_df[, common_features, drop = FALSE]))
if (!class_var %in% colnames(train_labels_df)) {
  stop("Training label column not found: ", class_var)
}
y_train <- ifelse(train_labels_df[[class_var]] == "Normal", 0, 1)
if (length(unique(y_train)) != 2) {
  stop("Training labels must contain exactly two classes after binarization (Normal=0, Other=1).")
}
bio_col <- rainbow(length(genes_to_plot), s = 0.9, v = 0.9)
auc_text <- character(0)
pdf(file = "ROC.genes.pdf", width = 9, height = 9)
for (k in seq_along(genes_to_plot)) {
  g <- genes_to_plot[k]
  if (!g %in% rownames(expr_train)) {
    warning("Gene not found in expression matrix: ", g)
    next
  }
  roc_g <- roc(y_train, as.numeric(expr_train[g, ]))
  if (k == 1) {
    plot(roc_g, print.auc = FALSE, col = bio_col[k], legacy.axes = TRUE, main = "", lwd = 3)
  } else {
    plot(roc_g, print.auc = FALSE, col = bio_col[k], legacy.axes = TRUE, main = "", lwd = 3, add = TRUE)
  }
  auc_text <- c(auc_text, paste0(g, ", AUC=", sprintf("%.3f", as.numeric(roc_g$auc))))
}
legend("bottomright", auc_text, lwd = 3, bty = "n", cex = 0.8, col = bio_col[seq_along(auc_text)], inset = c(0.05, 0))
dev.off()