# =============================================================================
#  ADAPT-NK TARGET SCREENING PIPELINE
#  Script 04 — Cluster Annotation
#
#  Input : checkpoint/S02_characterization_export.rds
#          results/S03_results09_AUC_matrix.csv
#          results/S03_results10_DEA_regulon.csv
#  Output: results/S04_results11_DEA_regulon_sig.csv
#          results/S04_results12_source_by_annotated_cluster.csv
#          plots/S04_plot07_DEA_regulon_volcano_plot.pdf
#          plots/S04_plot08_UMAP_annotated.pdf
#          checkpoint/S04_annotation.rds
#          session/S04_session_info.txt
# =============================================================================

# =============================================================================
#  IDENTITY MAPPING
# =============================================================================
#
#   Mature Cytotoxic NK   = c("FGFBP2","GZMB","PRF1","FCGR3A","CX3CR1","SPON2")
#   Tissue-resident NK    = c("CD69","ITGA1","ITGAE","ZNF683","CXCR6")
#   Helper NK             = c("GZMK","SELL","IL7R","XCL1","TCF7")
#   Activated Stressed NK = c("FOS","JUN","JUNB","EGR1","DUSP1","HSPA1A")
#   Proliferating NK      = c("MKI67","TOP2A","STMN1","TYMS","UBE2C")
# 
#  Each panel is a standard, literature-consistent NK subset signature
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

set.seed(42)

COL_BOTH <- "#7B4F9E"
COL_UP   <- "#B2182B"
COL_DOWN <- "#2166AC"
COL_NS   <- "grey80"

cluster_labels <- c(
  "1" = "Mature Cytotoxic NK",
  "2" = "Tissue-resident NK",
  "3" = "Helper NK",
  "4" = "Activated Stressed NK",
  "5" = "Proliferating NK"
)

cluster_colours <- c(
  "Mature Cytotoxic NK"   = "#00BF7D",
  "Tissue-resident NK"    = "#A3A500",
  "Helper NK"             = "#F8766D",
  "Activated Stressed NK" = "#00B0F6",
  "Proliferating NK"      = "#E76BF3"
)

# ── Load checkpoint ───────────────────────────────────────────────────────────
NKG2C_gated_obj <- readRDS("checkpoint/S02_characterization_export.rds")$NKG2C_gated_obj

# =============================================================================
#  SECTION 1: Merge pySCENIC Regulons
# =============================================================================
# Merge pySCENIC results for further analysis

auc_matrix <- read.csv("results/S03_results09_AUC_matrix.csv", row.names = 1,
                       stringsAsFactors = FALSE, check.names = FALSE)
auc_matrix <- auc_matrix[colnames(NKG2C_gated_obj), ]
for (regulon_name in colnames(auc_matrix)) {
  NKG2C_gated_obj[[regulon_name]] <- auc_matrix[[regulon_name]]
}

# =============================================================================
#  SECTION 2: Regulon Volcano Plot Visualization
# =============================================================================
# Visualize regulon DEA on cluster of interest vs other clusters via volcano plot

DEA_regulon <- read.csv("results/S03_results10_DEA_regulon.csv",
                        stringsAsFactors = FALSE)

cluster_id_regulon <- sort(unique(DEA_regulon$NK_cluster))
make_regulon_colour_scale <- function() {
  ggplot2::scale_colour_manual(
    name   = "Regulon category",
    breaks = c(COL_BOTH, COL_UP, COL_DOWN, COL_NS),
    values = c(COL_BOTH, COL_UP, COL_DOWN, COL_NS),
    limits = c(COL_BOTH, COL_UP, COL_DOWN, COL_NS),
    labels = c("Top 5 significant and extreme Log2FC",
               "Top 5 significant or highest Log2FC",
               "Top 5 significant or lowest Log2FC",
               "Not significant"),
    drop   = FALSE,
    guide  = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 1, shape = 16, stroke = 1)
    )
  )
}

volcano_panels_regulon <- lapply(cluster_id_regulon, function(cl_reg) {

  df_DEA_padj_regulon <- DEA_regulon %>%
    dplyr::filter(NK_cluster == cl_reg, padj < 0.05, avg_log2FC > 0.25) %>%
    dplyr::arrange(padj, dplyr::desc(avg_log2FC)) %>%
    dplyr::slice_min(padj, n = 5) %>%
    dplyr::pull(regulon)

  df_DEA_log2FC_regulon <- DEA_regulon %>%
    dplyr::filter(NK_cluster == cl_reg, padj < 0.05, avg_log2FC > 0.25) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC), padj) %>%
    dplyr::slice_max(avg_log2FC, n = 5) %>%
    dplyr::pull(regulon)

  df_DEA_down_padj_regulon <- DEA_regulon %>%
    dplyr::filter(NK_cluster == cl_reg, padj < 0.05, avg_log2FC < -0.25) %>%
    dplyr::arrange(padj, avg_log2FC) %>%
    dplyr::slice_min(padj, n = 5) %>%
    dplyr::pull(regulon)

  df_DEA_down_log2FC_regulon <- DEA_regulon %>%
    dplyr::filter(NK_cluster == cl_reg, padj < 0.05, avg_log2FC < -0.25) %>%
    dplyr::arrange(avg_log2FC, padj) %>%
    dplyr::slice_min(avg_log2FC, n = 5) %>%
    dplyr::pull(regulon)

  df_DEA_overlap_regulon <- intersect(df_DEA_padj_regulon, df_DEA_log2FC_regulon)
  df_DEA_down_overlap_regulon <- intersect(df_DEA_down_padj_regulon,
                                           df_DEA_down_log2FC_regulon)
  DEA_annotated_regulon <- data.frame(regulon = unique(c(df_DEA_padj_regulon,
                                                         df_DEA_log2FC_regulon,
                                                         df_DEA_down_padj_regulon,
                                                         df_DEA_down_log2FC_regulon)),
                                      stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      DEA_color_regulon = dplyr::case_when(
        regulon %in% c(df_DEA_overlap_regulon, df_DEA_down_overlap_regulon)  ~ COL_BOTH,
        regulon %in% c(df_DEA_padj_regulon, df_DEA_log2FC_regulon)           ~ COL_UP,
        regulon %in% c(df_DEA_down_padj_regulon, df_DEA_down_log2FC_regulon) ~ COL_DOWN,
        TRUE ~ "grey80"))


  df_volcano_all_regulon <- DEA_regulon %>%
    dplyr::filter(NK_cluster == cl_reg) %>%
    dplyr::distinct(regulon, .keep_all = TRUE) %>%     
    dplyr::left_join(DEA_annotated_regulon, by = "regulon") %>%
    dplyr::mutate(neglog10p = -log10(padj + 1e-300),
                  significance_regulon = dplyr::case_when(
                    padj < 0.05 & avg_log2FC >=  0.25 ~ COL_UP,
                    padj < 0.05 & avg_log2FC <= -0.25 ~ COL_DOWN,
                    TRUE                              ~ COL_NS),
                  point_col_regulon = dplyr::coalesce(DEA_color_regulon,
                                                      significance_regulon),
                  alpha_val_regulon = dplyr::if_else(!is.na(DEA_color_regulon),
                                                     0.85, 0.4))
  df_annotated_regulon <- dplyr::filter(df_volcano_all_regulon,
                                        !is.na(DEA_color_regulon))
  
  ggplot2::ggplot(df_volcano_all_regulon,
                  ggplot2::aes(x = avg_log2FC, y = neglog10p)) +
    ggplot2::geom_point(
      ggplot2::aes(colour = point_col_regulon, alpha = alpha_val_regulon),
      size = 1.5, stroke = 0, shape  = 16
    ) +
    make_regulon_colour_scale() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25))
    ) +
    ggplot2::scale_alpha_identity() +
    
    ggplot2::geom_vline(xintercept = c(-0.25, 0.25),
                        linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        linetype = "dashed", colour = "grey50") +
    
    ggplot2::geom_point(
      data  = df_annotated_regulon,
      ggplot2::aes(colour = point_col_regulon),
      shape  = 16,
      size   = 3,
      stroke = 1
    ) +
    ggrepel::geom_text_repel(
      data               = df_annotated_regulon,
      ggplot2::aes(
        label    = regulon,
        colour   = point_col_regulon,
        fontface = ifelse(point_col_regulon == COL_BOTH, "bold", "plain")
      ),
      size               = 2.5,
      max.overlaps       = Inf,
      show.legend        = FALSE,
      segment.size       = 0.3,
      segment.alpha      = 0.6,
      box.padding        = 0.6,
      seed               = 42,
      force              = ifelse(cl_reg == "1", 40, 8),
      force_pull         = ifelse(cl_reg == "1", 0.02, 0.1)
    ) +
    ggplot2::labs(
      title = paste0("Cluster ", cl_reg, " vs Other Clusters"),
      x     = "Log2FC",
      y     = expression(-log[10]~padj)
    ) +    
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
})

names(volcano_panels_regulon) <- as.character(cluster_id_regulon)

legend_source_df <- data.frame(
  x = 1:4, y = 1:4,
  point_col_regulon = c(COL_BOTH, COL_UP, COL_DOWN, COL_NS)
)
legend_plot <- ggplot2::ggplot(
  legend_source_df,
  ggplot2::aes(x, y, colour = point_col_regulon)
) +
  ggplot2::geom_point(size = 3) +
  make_regulon_colour_scale() +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position   = "right",
    legend.background = ggplot2::element_blank(),
    legend.key        = ggplot2::element_blank(),
    plot.background   = ggplot2::element_blank()
  )

regulon_legend_grob <- cowplot::get_legend(legend_plot)

combined_regulon_panels <- patchwork::wrap_plots(
  volcano_panels_regulon[as.character(cluster_id_regulon)],
  nrow = 1
)

plot_regulon_volcano <- (
  combined_regulon_panels | patchwork::wrap_elements(regulon_legend_grob)
) +
  patchwork::plot_layout(widths = c(rep(1, length(cluster_id_regulon)), 0.6)) +
  patchwork::plot_annotation(
    title = "Differential Regulon Activity Across NK Cluster",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  )

pdf("plots/S04_plot07_DEA_regulon_volcano_plot.pdf",
    width  = 5 * length(cluster_id_regulon) + 2,
    height = 7)
print(plot_regulon_volcano)
dev.off()

# =============================================================================
#  SECTION 3: Significant Regulons Filtering
# =============================================================================
# Filter on significant regulons

DEA_regulon_sig <- DEA_regulon %>%
  dplyr::filter(padj < 0.05, abs(avg_log2FC) >= 1) %>%
  dplyr::arrange(NK_cluster, padj, dplyr::desc(avg_log2FC))

write.csv(DEA_regulon_sig, "results/S04_results11_DEA_regulon_sig.csv",
          row.names = FALSE)

# =============================================================================
#  SECTION 4: Cluster Identity Assignment
# =============================================================================
# Assign cluster identity after literature confirmation

NKG2C_gated_obj$NK_cluster_annotated <- factor(
  unname(cluster_labels[as.character(NKG2C_gated_obj$NK_cluster)]),
  levels = names(cluster_colours))

# Calculate number of cells of each source per cluster

cells_source_per_cluster <- NKG2C_gated_obj@meta.data %>%
  dplyr::count(NK_cluster_annotated, source, name = "absolute_count") %>%
  dplyr::group_by(NK_cluster_annotated) %>%
  dplyr::mutate(pct = round(absolute_count / sum(absolute_count) * 100, 1)) %>%
  dplyr::arrange(NK_cluster_annotated, dplyr::desc(pct))

write.csv(cells_source_per_cluster,
          "results/S04_results12_source_by_annotated_cluster.csv",
          row.names = FALSE)

pdf("plots/S04_plot08_UMAP_annotated.pdf", width = 8, height = 6)
print(
  DimPlot(NKG2C_gated_obj,
          reduction = "umap",
          group.by  = "NK_cluster_annotated",
          cols      = cluster_colours) +
    ggtitle("NK Clusters Annotated"))
dev.off()

# =============================================================================
#  Save Checkpoint
# =============================================================================

saveRDS(NKG2C_gated_obj, "checkpoint/S04_annotation.rds")
writeLines(capture.output(sessionInfo()), "session/S04_session_info.txt")
