# ============================================================
# GO:BP enrichment, secondary vs primary tumour.
# Input: edgeR-style DE table, contrast = secondary vs primary.
# Convention: logFC > 0  =>  HIGHER in secondary.
# Extends the earlier angiogenesis-only script (Lrg1 check +
# HALLMARK/GOBP_ANGIOGENESIS) to the full GO:BP collection, so we
# can see which biological-process terms top the ranked list in
# either direction, not just angiogenesis.
# ============================================================

library(data.table)
library(fgsea)
library(msigdbr)

# -----------------------------
# 1. LOAD DE TABLE
# -----------------------------
de <- fread("../data/secondary_vs_primary.csv", header = TRUE)
setnames(de, c("gene", "logFC", "logCPM", "stat", "PValue", "FDR", "direction", "label"))
de <- de[, .(gene, logFC, logCPM, stat, PValue, FDR, direction)]
de <- de[!is.na(logFC) & !is.na(gene) & is.finite(logFC)]
de <- de[!duplicated(gene)]
message(sprintf("Genes in table: %d  (%d up in secondary, %d down)",
                nrow(de), sum(de$logFC > 0), sum(de$logFC < 0)))

# -----------------------------
# 2. TARGETED CHECK: Lrg1
# -----------------------------
lrg1 <- de[grepl("^Lrg1$", gene, ignore.case = TRUE)]
if (nrow(lrg1) == 0) {
  message("Lrg1 not found in table (filtered out, or different symbol).")
} else {
  print(lrg1[, .(gene, logFC, PValue, FDR, direction)])
  message(if (lrg1$logFC > 0)
    sprintf("Lrg1 is HIGHER in secondary (logFC=%.2f, FDR=%.3g)", lrg1$logFC, lrg1$FDR)
    else
      sprintf("Lrg1 is LOWER in secondary (logFC=%.2f, FDR=%.3g)", lrg1$logFC, lrg1$FDR))
}

# -----------------------------
# 3. FULL GO:BP COLLECTION (mouse symbols)
# -----------------------------
# msigdbr >= 10 uses msigdbr(species, collection="C5", subcollection="GO:BP");
# older versions use category/subcategory. Try the new API first, fall back.
m <- tryCatch(
  as.data.table(msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")),
  error = function(e) as.data.table(msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP"))
)
pathways <- split(m$gene_symbol, m$gs_name)
pathways <- lapply(pathways, unique)
message(sprintf("GO:BP gene sets loaded: %d", length(pathways)))

# -----------------------------
# 4. GSEA ON THE FULL RANKED LIST
#    (needs ALL genes, not just the significant ones)
# -----------------------------
ranks <- de$logFC
names(ranks) <- de$gene
ranks <- sort(ranks, decreasing = TRUE)

set.seed(1)
fg <- fgsea(pathways = pathways, stats = ranks, eps = 0, minSize = 10, maxSize = 500)
fg[, leadingEdge := sapply(leadingEdge, function(x) paste(head(x, 10), collapse = ","))]
setorder(fg, padj)

fwrite(fg, "../results/gobp_gsea_full.csv")

# -----------------------------
# 5. TOP TERMS IN EACH DIRECTION
# -----------------------------
cat("\n=== Top 20 GO:BP terms HIGHER in secondary (NES > 0), by padj ===\n")
print(fg[NES > 0][order(padj)][1:20, .(pathway, NES, pval, padj, size, leadingEdge)])

cat("\n=== Top 20 GO:BP terms HIGHER in primary (NES < 0), by padj ===\n")
print(fg[NES < 0][order(padj)][1:20, .(pathway, NES, pval, padj, size, leadingEdge)])

n_sig <- sum(fg$padj < 0.05, na.rm = TRUE)
cat(sprintf("\n%d / %d tested GO:BP terms significant at padj < 0.05\n", n_sig, nrow(fg)))

# -----------------------------
# 6. WHERE DOES ANGIOGENESIS RANK AMONG ALL GO:BP TERMS?
# -----------------------------
angio_terms <- grep("ANGIOGEN|VASCULAR|VASCULATURE|ENDOTHELI|BLOOD_VESSEL",
                     fg$pathway, value = TRUE, ignore.case = TRUE)
cat(sprintf("\n=== %d angiogenesis/vascular-related GO:BP terms tested ===\n", length(angio_terms)))
print(fg[pathway %in% angio_terms][order(padj), .(pathway, NES, pval, padj, size, leadingEdge)])
