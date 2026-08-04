# ── Full TME composition + expression analysis ───────────────────────────────
# Uses seurat_all_cells_integrated.rds (tumour + full TME + LPT_WT), the object
# saved right before CopyKAT-based tumour filtering — so cell-type proportions
# here reflect the real tissue composition, not the leaky ~14%-contaminated
# "tumour-only" subset.

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

merged <- readRDS("seurat_all_cells_integrated.rds")
DefaultAssay(merged) <- "SCT"

# ── 1. Orientation: what's in this object? ───────────────────────────────────
table(merged$condition)
table(merged$is_tumour)
table(merged$seurat_clusters)

DimPlot(merged, group.by = "seurat_clusters", label = TRUE) + ggtitle("All cells: clusters")
DimPlot(merged, group.by = "condition") + ggtitle("All cells: condition")
DimPlot(merged, group.by = "is_tumour") + ggtitle("All cells: CopyKAT tumour call")

# ── 2. Marker genes per cluster (identify TME cell types properly this time) ─
Idents(merged) <- "seurat_clusters"
merged <- PrepSCTFindMarkers(merged)   # needed before FindMarkers/FindAllMarkers on SCT
all_markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, "all_cells_cluster_markers_full.csv", row.names = FALSE)

top10 <- all_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
write.csv(top10, "all_cells_cluster_markers_top10.csv", row.names = FALSE)
print(top10, n = 300)

# ── 3. Cluster composition: patient (batch check) and condition ─────────────
cat("\n--- Cluster x Patient ---\n")
print(table(merged$seurat_clusters, merged$patient))
cat("\n--- Cluster x Patient, row-normalised ---\n")
print(round(prop.table(table(merged$seurat_clusters, merged$patient), margin = 1), 2))

cat("\n--- Cluster x Condition ---\n")
print(table(merged$seurat_clusters, merged$condition))
cat("\n--- Cluster x Condition, row-normalised ---\n")
print(round(prop.table(table(merged$seurat_clusters, merged$condition), margin = 1), 2))

# ── 4. Formal test: does cluster proportion differ MET vs Primary? ──────────
# Paired by patient (LPT_WT samples aren't patient-matched to MET/Primary, so
# excluded from this specific test — shown separately below for reference).
props <- merged@meta.data %>%
  filter(condition %in% c("LPT_MET", "Primary")) %>%
  count(patient, condition, seurat_clusters) %>%
  group_by(patient, condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_wide <- props %>%
  select(patient, condition, seurat_clusters, prop) %>%
  pivot_wider(names_from = condition, values_from = prop, values_fill = 0)

cluster_prop_test <- prop_wide %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_MET     = mean(LPT_MET, na.rm = TRUE),
    mean_Primary = mean(Primary, na.rm = TRUE),
    p_value      = tryCatch(wilcox.test(LPT_MET, Primary, paired = TRUE)$p.value,
                             error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  arrange(p_value)

print(cluster_prop_test, n = 30)
write.csv(cluster_prop_test, "cluster_proportion_test_MET_vs_Primary.csv", row.names = FALSE)

# ── 5. Reference: composition including LPT_WT (non-tumour baseline, unpaired)
cat("\n--- Cluster proportions by condition, all 3 groups (unpaired reference) ---\n")
props_all <- merged@meta.data %>%
  count(condition, seurat_clusters) %>%
  group_by(condition) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup() %>%
  pivot_wider(id_cols = seurat_clusters, names_from = condition, values_from = prop, values_fill = 0)
print(props_all, n = 30)

# ── 6. Stacked bar plot of composition per patient/condition ────────────────
ggplot(props, aes(x = interaction(patient, condition), y = prop, fill = seurat_clusters)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Patient / Condition", y = "Proportion of cells", fill = "Cluster",
       title = "Cell composition by patient and condition")
ggsave("tme_composition_barplot.pdf", width = 10, height = 6)
