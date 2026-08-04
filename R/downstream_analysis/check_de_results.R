# ── Pull specific genes out of the pseudobulk DESeq2 results ────────────────

de <- read.csv("DE_LPT_MET_vs_Primary.csv", row.names = 1)

genes_of_interest <- c("Trem1", "Pglyrp1", "Lum", "Thbd", "Kcnj8",
                        "Col3a1", "Pdgfa", "App", "Lrg1")

result <- de[de$gene %in% genes_of_interest, c("gene", "avg_log2FC", "p_val", "p_val_adj")]
result <- result[order(result$p_val_adj), ]

# Flag which genes weren't found at all (e.g. filtered out for low expression)
missing <- setdiff(genes_of_interest, result$gene)
if (length(missing) > 0) {
  message("Not present in DE results (likely filtered for low expression): ",
          paste(missing, collapse = ", "))
}

print(result, row.names = FALSE)
write.csv(result, "vascular_genes_DESeq2_results.csv", row.names = FALSE)
