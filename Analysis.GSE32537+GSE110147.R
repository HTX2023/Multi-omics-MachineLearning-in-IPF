# ==================================================================================================
# GEO expression preprocessing + phenotype parsing + platform annotation (Entrez -> Symbol) +
# differential expression analysis (tinyarray::get_deg_all) + exporting results.
#
# Inputs (local files):
#   - exp_data.csv               : expression matrix (rows = probes/ids, cols = samples)
#   - phenotype_info.csv         : phenotype/metadata table with at least column "title"
#   - GSE151052_family.soft.gz   : SOFT family file containing platform annotation tables
#
# Outputs:
#   - GSE151052_preprocessed.RData : saved objects (exp, Group, clinical) for downstream analysis
#   - deg_filtered_output.xlsx     : DEGs (up/down only) exported to Excel
#
# Notes:
#   - This script assumes that "CTL" and "IPF" are distinguishable from pd$title.
#   - Platform mapping uses org.Hs.eg.db and expects gpl$ORF to contain Entrez IDs.
#   - If duplicated gene symbols occur after mapping, the script keeps the first occurrence.
#   - If missing symbols occur after mapping, the script drops those rows.
# ==================================================================================================
rm(list = ls())
options(timeout = 100000)
options(stringsAsFactors = FALSE, check.names = FALSE)
suppressPackageStartupMessages({
  library(GEOquery)
  library(tinyarray)
  library(stringr)
  library(org.Hs.eg.db)
  library(writexl)
})
exp <- read.table("exp_data.csv", header = TRUE, row.names = 1, sep = ",")
pd  <- read.table("phenotype_info.csv", header = TRUE, row.names = 1, sep = ",")
if (!"title" %in% colnames(pd)) stop("Column 'title' not found in phenotype_info.csv.")
is_ctl <- str_detect(pd$title, "CTL")
is_ipf <- str_detect(pd$title, "IPF")
pd <- pd[is_ctl | is_ipf, , drop = FALSE]
if (nrow(pd) == 0) stop("No samples matched 'CTL' or 'IPF' in pd$title.")
if (!identical(rownames(pd), colnames(exp))) {
  common_samples <- intersect(rownames(pd), colnames(exp))
  if (length(common_samples) == 0) stop("No overlapping sample IDs between pd and exp.")
  exp <- exp[, common_samples, drop = FALSE]
  pd  <- pd[common_samples, , drop = FALSE]
}
is_ctl2 <- str_detect(pd$title, "CTL")
Group <- ifelse(is_ctl2, "Normal", "IPF")
Group <- factor(Group, levels = c("Normal", "IPF"))
print(table(Group))
print(head(data.frame(title = pd$title, Group = Group), 10))
soft_obj <- getGEO(filename = "GSE151052_family.soft.gz")
if (length(soft_obj@gpls) == 0) stop("No GPL objects found in SOFT file: GSE151052_family.soft.gz")
gpl_tbl <- soft_obj@gpls[[1]]@dataTable@table
if (!all(c("ID", "ORF") %in% colnames(gpl_tbl))) {
  stop("Expected columns 'ID' and 'ORF' not found in GPL table extracted from SOFT file.")
}
gpl_tbl$SYMBOL_FROM_DB <- mapIds(
  x = org.Hs.eg.db,
  keys = as.character(gpl_tbl$ORF),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)
mapped_symbols <- gpl_tbl$SYMBOL_FROM_DB[match(rownames(exp), gpl_tbl$ID)]
has_missing <- is.na(mapped_symbols) | mapped_symbols == ""
if (any(has_missing)) {
  warning("Missing gene symbols detected after mapping; corresponding rows will be removed.")
  exp <- exp[!has_missing, , drop = FALSE]
  mapped_symbols <- mapped_symbols[!has_missing]
}
is_dup <- duplicated(mapped_symbols)
if (any(is_dup)) {
  warning("Duplicated gene symbols detected after mapping; keeping the first occurrence.")
  exp <- exp[!is_dup, , drop = FALSE]
  mapped_symbols <- mapped_symbols[!is_dup]
}
rownames(exp) <- mapped_symbols
ids_combined <- data.frame(
  probe_id = rownames(exp),
  symbol   = rownames(exp),
  stringsAsFactors = FALSE
)
dcp <- get_deg_all(
  exp,
  Group,
  ids = ids_combined,
  entriz = FALSE,
  adjust = TRUE,
  logFC_cutoff = 0.58,
  cluster_cols = TRUE
)
print(head(dcp$deg))
print(table(dcp$deg$change))
print(dcp$plots)
clinical <- data.frame(Group = Group, stringsAsFactors = FALSE)
rownames(clinical) <- rownames(pd)
save(exp, Group, clinical, file = "GSE151052_preprocessed.RData")
deg_filtered <- dcp$deg[dcp$deg$change %in% c("up", "down"), , drop = FALSE]
print(head(deg_filtered))
write_xlsx(deg_filtered, "deg_filtered_output.xlsx")