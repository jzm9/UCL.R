# ── Manual cluster annotation (all-cells object) ─────────────────────────
# Labels derived from FindAllMarkers() output in tme_composition_analysis.R.
# Apply this after loading seurat_all_cells_integrated.rds.

library(Seurat)
library(dplyr)
library(ggplot2)

merged <- readRDS("seurat_all_cells_integrated.rds")

cluster_labels <- c(
  "0"  = "Tumour",
  "1"  = "Tumour",
  "2"  = "Smooth muscle",
  "3"  = "Tumour (S/G2M-phase)",
  "4"  = "Tumour (S/G2M-phase)",
  "5"  = "Pericytes",
  "6"  = "Tumour (S/G2M-phase)",
  "7"  = "Oligodendrocyte-like/glia",
  "8"  = "BBB endothelial",
  "9"  = "Schwann cells",
  "10" = "Fibroblasts",
  "11" = "Inflammatory macrophages/monocytes",
  "12" = "Microglia/macrophages",
  "13" = "Schwann cells",
  "14" = "Sex-linked/technical",
  "15" = "Tumour (stress/IEG state)",
  "16" = "Astrocytes",
  "17" = "Neutrophils",
  "18" = "T cells / NK cells",
  "19" = "Endothelial (Pecam1/Sox17-high)",
  "20" = "Homeostatic microglia",
  "21" = "Endothelial (Kdr/Mmrn2-high)",
  "22" = "Low-quality"
)

merged$cell_type <- factor(unname(cluster_labels[as.character(merged$seurat_clusters)]),
                            levels = unique(cluster_labels))

# Broad tumour vs TME flag, useful as a sanity check against CopyKAT's is_tumour
merged$compartment <- ifelse(grepl("^Tumour", merged$cell_type), "Tumour", "TME")

pdf("umap_all_cells_labelled.pdf", width = 10, height = 7)
print(DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE) +
        ggtitle("All cells: annotated clusters"))
print(DimPlot(merged, group.by = "compartment") +
        ggtitle("All cells: tumour vs TME"))
dev.off()

saveRDS(merged, "seurat_all_cells_integrated_labelled.rds")

# ── Composition by labelled cell type, MET vs Primary (paired by patient) ──
props <- merged@meta.data %>%
  filter(condition %in% c("LPT_MET", "Primary")) %>%
  count(patient, condition, cell_type) %>%
  group_by(patient, condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_wide <- props %>%
  select(patient, condition, cell_type, prop) %>%
  tidyr::pivot_wider(names_from = condition, values_from = prop, values_fill = 0)

cell_type_prop_test <- prop_wide %>%
  group_by(cell_type) %>%
  summarise(
    mean_MET     = mean(LPT_MET, na.rm = TRUE),
    mean_Primary = mean(Primary, na.rm = TRUE),
    p_value      = tryCatch(wilcox.test(LPT_MET, Primary, paired = TRUE)$p.value,
                             error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  arrange(p_value)

print(cell_type_prop_test, n = 30)
write.csv(cell_type_prop_test, "cell_type_proportion_test_MET_vs_Primary.csv", row.names = FALSE)

ggplot(props, aes(x = interaction(patient, condition), y = prop, fill = cell_type)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Patient / Condition", y = "Proportion of cells", fill = "Cell type",
       title = "Cell composition by patient and condition (labelled)")
ggsave("cell_type_composition_barplot.pdf", width = 11, height = 6)
