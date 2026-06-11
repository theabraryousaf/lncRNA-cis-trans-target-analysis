suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript 07_trans_targets_pearson.R expression.tsv annotation.tsv outdir PCC_CUTOFF FDR_CUTOFF")
}

expr_file <- args[1]
ann_file <- args[2]
outdir <- args[3]
PCC_CUTOFF <- as.numeric(args[4])
FDR_CUTOFF <- as.numeric(args[5])

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

expr <- fread(expr_file)
ann <- fread(ann_file)

setnames(expr, names(expr)[1], "gene_id")
expr[, gene_id := sub("\\..*$", "", as.character(gene_id))]
ann[, gene_id_base := sub("\\..*$", "", as.character(gene_id_base))]

sample_cols <- setdiff(names(expr), "gene_id")
if (length(sample_cols) < 3) {
  stop("ERROR: At least 3 samples are required for Pearson correlation.")
}

expr[, (sample_cols) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = sample_cols]
expr <- expr[, lapply(.SD, mean, na.rm = TRUE), by = gene_id, .SDcols = sample_cols]

lnc_ids <- ann[gene_type == "lncRNA", unique(gene_id_base)]
pcg_ids <- ann[gene_type == "protein_coding", unique(gene_id_base)]

common_lnc <- intersect(lnc_ids, expr$gene_id)
common_pcg <- intersect(pcg_ids, expr$gene_id)

cat("lncRNAs used:", length(common_lnc), "\n")
cat("protein-coding genes used:", length(common_pcg), "\n")

if (length(common_lnc) == 0) stop("ERROR: No lncRNAs found in expression matrix.")
if (length(common_pcg) == 0) stop("ERROR: No protein-coding genes found in expression matrix.")

lnc_mat <- as.matrix(expr[match(common_lnc, gene_id), ..sample_cols])
pcg_mat <- as.matrix(expr[match(common_pcg, gene_id), ..sample_cols])
rownames(lnc_mat) <- common_lnc
rownames(pcg_mat) <- common_pcg
storage.mode(lnc_mat) <- "numeric"
storage.mode(pcg_mat) <- "numeric"

keep_variable <- function(mat) {
  rowSums(!is.na(mat)) >= 3 & apply(mat, 1, function(z) sd(z, na.rm = TRUE) > 0)
}

lnc_mat <- lnc_mat[keep_variable(lnc_mat), , drop = FALSE]
pcg_mat <- pcg_mat[keep_variable(pcg_mat), , drop = FALSE]

cat("variable lncRNAs used:", nrow(lnc_mat), "\n")
cat("variable protein-coding genes used:", nrow(pcg_mat), "\n")

rmat <- cor(t(lnc_mat), t(pcg_mat), method = "pearson", use = "pairwise.complete.obs")
nmat <- (!is.na(lnc_mat)) %*% t(!is.na(pcg_mat))
df <- nmat - 2

rclip <- pmin(pmax(rmat, -0.999999999999), 0.999999999999)
tmat <- rclip * sqrt(df / pmax(1 - rclip^2, .Machine$double.eps))
pmat <- 2 * pt(-abs(tmat), df = df)
pmat[df < 1 | is.na(rmat)] <- NA_real_

res <- as.data.table(as.data.frame(as.table(rmat), stringsAsFactors = FALSE))
setnames(res, c("lncRNA_gene_id", "target_gene_id", "PCC"))
res[, n_samples := as.vector(nmat)]
res[, pvalue := as.vector(pmat)]
res[, FDR := p.adjust(pvalue, method = "BH")]
res[, abs_PCC := abs(PCC)]
res[, correlation_direction := fifelse(PCC > 0, "positive", fifelse(PCC < 0, "negative", "zero"))]

lnc_info <- ann[gene_type == "lncRNA", .(
  lncRNA_gene_id = gene_id_base,
  lncRNA_symbol = gene_symbol
)]

target_info <- ann[gene_type == "protein_coding", .(
  target_gene_id = gene_id_base,
  target_symbol = gene_symbol
)]

res <- merge(res, lnc_info, by = "lncRNA_gene_id", all.x = TRUE)
res <- merge(res, target_info, by = "target_gene_id", all.x = TRUE)

res[is.na(lncRNA_symbol) | lncRNA_symbol == "", lncRNA_symbol := lncRNA_gene_id]
res[is.na(target_symbol) | target_symbol == "", target_symbol := target_gene_id]

setorder(res, FDR, -abs_PCC)

all_file <- file.path(outdir, "trans_targets_Pearson_all.tsv")
selected_file <- file.path(outdir, "trans_targets_Pearson_selected.tsv")
priority_file <- file.path(outdir, "trans_lncRNA_priority.tsv")

fwrite(res, all_file, sep = "\t")
fwrite(res, sub("\\.tsv$", ".csv", all_file))

selected <- res[!is.na(PCC) & abs_PCC >= PCC_CUTOFF & FDR <= FDR_CUTOFF]
setorder(selected, FDR, -abs_PCC)
fwrite(selected, selected_file, sep = "\t")
fwrite(selected, sub("\\.tsv$", ".csv", selected_file))

make_gene_lists <- function(dt, direction, out_prefix) {
  if (direction != "all") {
    dt <- dt[correlation_direction == direction]
  }
  genes <- unique(dt[, .(target_gene_id, target_symbol)])
  setorder(genes, target_symbol)
  fwrite(genes, file.path(outdir, paste0(out_prefix, ".tsv")), sep = "\t")
  fwrite(genes, file.path(outdir, paste0(out_prefix, ".csv")))
  writeLines(genes$target_symbol, file.path(outdir, paste0(out_prefix, ".txt")))
}

make_gene_lists(selected, "all", "all_trans_target_genes")
make_gene_lists(selected, "positive", "all_trans_positive_target_genes")
make_gene_lists(selected, "negative", "all_trans_negative_target_genes")

if (nrow(selected) == 0) {
  priority <- data.table(
    lncRNA_gene_id = character(),
    lncRNA_symbol = character(),
    trans_target_count = integer(),
    positive_trans_target_count = integer(),
    negative_trans_target_count = integer(),
    max_abs_PCC = numeric(),
    trans_targets = character()
  )
} else {
  priority <- selected[, .(
    trans_target_count = uniqueN(target_gene_id),
    positive_trans_target_count = uniqueN(target_gene_id[correlation_direction == "positive"]),
    negative_trans_target_count = uniqueN(target_gene_id[correlation_direction == "negative"]),
    max_abs_PCC = max(abs_PCC, na.rm = TRUE),
    trans_targets = paste(sort(unique(target_symbol)), collapse = ";")
  ), by = .(lncRNA_gene_id, lncRNA_symbol)]
  setorder(priority, -trans_target_count, -max_abs_PCC, lncRNA_symbol)
}

fwrite(priority, priority_file, sep = "\t")
fwrite(priority, sub("\\.tsv$", ".csv", priority_file))

cat("All lncRNA-PCG trans pairs tested:", nrow(res), "\n")
cat("Selected trans pairs:", nrow(selected), "\n")
cat("Selected positive trans pairs:", nrow(selected[correlation_direction == "positive"]), "\n")
cat("Selected negative trans pairs:", nrow(selected[correlation_direction == "negative"]), "\n")
cat("lncRNAs with selected trans targets:", nrow(priority), "\n")
cat("PCC cutoff:", PCC_CUTOFF, "\n")
cat("FDR cutoff:", FDR_CUTOFF, "\n")
