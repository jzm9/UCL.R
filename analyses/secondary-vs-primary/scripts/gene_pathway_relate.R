library(data.table)

de <- fread("../data/secondary_vs_primary.csv", header = TRUE)
setnames(de, c("gene", "logFC", "logCPM", "stat", "PValue", "FDR", "direction", "label"))
de <- de[, .(gene, logFC, FDR)]
de <- de[!is.na(logFC) & !is.na(gene) & is.finite(logFC)]
de <- de[!duplicated(gene)]
setorder(de, -logFC)
de[, rank_secondary := .I]
setorder(de, logFC)
de[, rank_primary := .I]
setkey(de, gene)

le <- fread("../results/leading_edges_top_pathways.csv")

up_sec  <- le[group == "TOP_UP_IN_SECONDARY"]
up_prim <- le[group == "TOP_UP_IN_PRIMARY"]

build_matrix <- function(pathway_tbl, rank_col, top_n_genes = 50) {
  top_genes <- de[order(if (rank_col == "rank_secondary") rank_secondary else rank_primary)][1:top_n_genes, gene]
  mat <- matrix(0L, nrow = length(top_genes), ncol = nrow(pathway_tbl),
                dimnames = list(top_genes, pathway_tbl$pathway))
  for (i in seq_len(nrow(pathway_tbl))) {
    genes_in_le <- trimws(strsplit(pathway_tbl$leadingEdge[i], ",")[[1]])
    mat[rownames(mat) %in% genes_in_le, i] <- 1L
  }
  mat
}

mat_sec  <- build_matrix(up_sec, "rank_secondary")
mat_prim <- build_matrix(up_prim, "rank_primary")

# long-format membership tables (for CSV / gene-level summary)
summarize_side <- function(mat, de_rank_col, side_label) {
  genes <- rownames(mat)
  n_pathways <- rowSums(mat)
  pathway_list <- apply(mat, 1, function(r) paste(colnames(mat)[r == 1], collapse = "; "))
  d <- de[gene %in% genes][order(if (side_label == "secondary") rank_secondary else rank_primary)]
  d[, n_top_pathways_leading_edge := n_pathways[gene]]
  d[, pathways := pathway_list[gene]]
  d[, side := side_label]
  d
}

sec_summary  <- summarize_side(mat_sec, "rank_secondary", "secondary")
prim_summary <- summarize_side(mat_prim, "rank_primary", "primary")

fwrite(sec_summary,  "../results/top_secondary_genes_vs_pathways.csv")
fwrite(prim_summary, "../results/top_primary_genes_vs_pathways.csv")

# also dump matrices as CSV (genes as rows) for the heatmap build step
fwrite(as.data.table(mat_sec, keep.rownames = "gene"),  "../results/matrix_secondary_genes_x_pathways.csv")
fwrite(as.data.table(mat_prim, keep.rownames = "gene"), "../results/matrix_primary_genes_x_pathways.csv")

cat("Secondary: top", nrow(mat_sec), "genes x", ncol(mat_sec), "pathways\n")
cat("Primary:  top", nrow(mat_prim), "genes x", ncol(mat_prim), "pathways\n")
cat("\nGenes in >=5 of the top 20 secondary pathways' leading edges:\n")
print(sec_summary[n_top_pathways_leading_edge >= 5, .(gene, logFC, FDR, rank_secondary, n_top_pathways_leading_edge)])
cat("\nGenes in >=5 of the top 20 primary pathways' leading edges:\n")
print(prim_summary[n_top_pathways_leading_edge >= 5, .(gene, logFC, FDR, rank_primary, n_top_pathways_leading_edge)])
