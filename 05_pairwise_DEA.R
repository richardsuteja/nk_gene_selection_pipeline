# =============================================================================
#  ADAPT-NK TARGET SCREENING PIPELINE
#  Script 05 — Pairwise Gene & Regulon DEA [VALIDATION]
#
#  Input : checkpoint/S04_annotation.rds
#          results/S03_results09_AUC_matrix.csv
#  Output: results/S05_results13_DEA_pairwise_dysfunctional_vs_effector.csv
#          results/S05_results14_DEA_pairwise_dysfunctional_vs_effector_sig.csv
#          results/S05_results15_regulon_DEA_pairwise_dysfunctional_vs_effector.csv
#          results/S05_results16_regulon_DEA_pairwise_dysfunctional_vs_effector_sig.csv
#          plots/S05_plot09_DEA_pairwise_dysfunctional_vs_effector_volcano_plot.pdf
#          plots/S05_plot10_regulon_DEA_pairwise_dysfunctional_vs_effector_volcano_plot.pdf
#          checkpoint/S05_pairwise_DEA.rds
#          session/S05_session_info.txt
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
set.seed(42)

# ── Load checkpoint ───────────────────────────────────────────────────────────
NKG2C_gated_obj <- readRDS("checkpoint/S04_annotation.rds")

NKG2C_gene_pairwise <- subset(NKG2C_gated_obj,
                         NK_cluster %in% c("1", "4"))
NKG2C_gene_pairwise$NK_cluster <- droplevels(factor(NKG2C_gene_pairwise$NK_cluster))
Idents(NKG2C_gene_pairwise) <- "NK_cluster"
DefaultAssay(NKG2C_gene_pairwise) <- "RNA"

# =============================================================================
#  SECTION 1: Pairwise Gene DEA on Dysfunctional vs Effector Cluster
# =============================================================================
# Pairwise Gene DEA on Selected Annotated NK Clusters

pairwise_DEA <- FindMarkers(
  NKG2C_gene_pairwise,
  ident.1         = "4",
  ident.2         = "1",
  test.use        = "wilcox",
  min.pct         = 0.10,
  logfc.threshold = 0,
  verbose         = FALSE) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC >=  1.00 ~ "Upregulated in Activated Stressed NK",
      p_val_adj < 0.05 & avg_log2FC <= -1.00 ~ "Upregulated in Mature Cytotoxic NK",
      TRUE ~ "NS"))

pairwise_DEA_sig <- dplyr::filter(pairwise_DEA, sig != "NS")

write.csv(pairwise_DEA,
          "results/S05_results13_DEA_pairwise_dysfunctional_vs_effector.csv",
          row.names = FALSE)
write.csv(pairwise_DEA_sig,
          "results/S05_results14_DEA_pairwise_dysfunctional_vs_effector_sig.csv",
          row.names = FALSE)

# Visualize via Volcano Plot

COL_BOTH <- "#7B4F9E"
COL_UP   <- "#B2182B"
COL_DOWN <- "#2166AC"
COL_NS   <- "grey80"

df_pairwise_DEA_up_padj <- pairwise_DEA %>%
  dplyr::filter(p_val_adj  <  0.05, avg_log2FC > 0) %>%
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_min(p_val_adj, n = 5) %>%
  dplyr::pull(gene)

df_pairwise_DEA_up_log2FC <- pairwise_DEA %>%
  dplyr::filter(p_val_adj  <  0.05, avg_log2FC > 0) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC), p_val_adj) %>%
  dplyr::slice_max(avg_log2FC, n = 5) %>%
  dplyr::pull(gene)

df_pairwise_DEA_up_overlap <- intersect(df_pairwise_DEA_up_padj,
                                        df_pairwise_DEA_up_log2FC)

df_pairwise_DEA_down_padj <- pairwise_DEA %>%
  dplyr::filter(p_val_adj  <  0.05, avg_log2FC < 0) %>%
  dplyr::arrange(p_val_adj, avg_log2FC) %>%
  dplyr::slice_min(p_val_adj, n = 5) %>%
  dplyr::pull(gene)

df_pairwise_DEA_down_log2FC <- pairwise_DEA %>%
  dplyr::filter(p_val_adj  <  0.05, avg_log2FC < 0) %>%
  dplyr::arrange(avg_log2FC, p_val_adj) %>%
  dplyr::slice_min(avg_log2FC, n = 5) %>%
  dplyr::pull(gene)

df_pairwise_DEA_down_overlap <- intersect(df_pairwise_DEA_down_padj,
                                          df_pairwise_DEA_down_log2FC)

pairwise_DEA_annotated_genes <- data.frame(gene = unique(c(df_pairwise_DEA_up_padj,
                                                           df_pairwise_DEA_up_log2FC,
                                                           df_pairwise_DEA_down_padj,
                                                           df_pairwise_DEA_down_log2FC)),
                                           stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    pairwise_DEA_color = dplyr::case_when(
      gene %in% c(df_pairwise_DEA_up_overlap, df_pairwise_DEA_down_overlap) ~ COL_BOTH,
      gene %in% c(df_pairwise_DEA_up_padj, df_pairwise_DEA_up_log2FC)       ~ COL_UP,
      gene %in% c(df_pairwise_DEA_down_padj, df_pairwise_DEA_down_log2FC)   ~ COL_DOWN,
      TRUE ~ "grey80"))

df_volcano_pairwise_DEA <- pairwise_DEA %>%
  dplyr::left_join(
    pairwise_DEA_annotated_genes[, c("gene", "pairwise_DEA_color")],
    by = "gene") %>%
  dplyr::mutate(neglog10p = -log10(p_val_adj + 1e-300),
    significance = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC >=  1.00 ~ COL_UP,
      p_val_adj < 0.05 & avg_log2FC <= -1.00 ~ COL_DOWN,
      TRUE ~ COL_NS),
    point_col = dplyr::coalesce(pairwise_DEA_color, significance),
    alpha_val = dplyr::if_else(!is.na(pairwise_DEA_color), 0.85, 0.4)
  )
df_annotated <- dplyr::filter(df_volcano_pairwise_DEA, !is.na(pairwise_DEA_color))

plot_volcano_pairwise_DEA <- ggplot2::ggplot(df_volcano_pairwise_DEA,
                ggplot2::aes(x = avg_log2FC, y = neglog10p)) +
  ggplot2::geom_point(ggplot2::aes(colour = point_col,alpha  = alpha_val),
    size   = 1.5,
    stroke = 0,
    shape  = 16
  ) +
  ggplot2::scale_colour_identity(
    guide  = ggplot2::guide_legend(
      title        = "Gene category",
      override.aes = list(size = 3)
    ),
    breaks = c(COL_BOTH, COL_UP, COL_DOWN, COL_NS),
    labels = c("Top 5 significant and extreme Log2FC",
               "Top 5 significant or highest Log2FC",
               "Top 5 significant or lowest Log2FC",
               "Not significant")
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::geom_vline(xintercept = c(-1.00, 1.00),
                      linetype = "dashed", colour = "grey50") +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(
    data   = df_annotated,
    ggplot2::aes(colour = point_col),
    shape  = 16,
    size   = 3,
    stroke = 1
  ) +
  ggrepel::geom_text_repel(
    data = df_annotated,
    ggplot2::aes(
      label    = gene,
      colour   = point_col,
      fontface = ifelse(point_col == COL_BOTH, "bold", "plain")
    ),
    size               = 2.5,
    max.overlaps       = Inf,
    show.legend        = FALSE,
    segment.size       = 0.3,
    segment.alpha      = 0.6,
    box.padding        = 0.6,
    seed               = 42,
    force              = 8
  ) +
  ggplot2::labs(
    title = "Differential Gene Expression on Dysfunctional vs Effector NK Cluster",
    x     = "Log2FC",
    y     = expression(-log[10]~padj)
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold")
  )

pdf("plots/S05_plot09_DEA_pairwise_dysfunctional_vs_effector_volcano_plot.pdf",
    width = 9, height = 7)
  print(plot_volcano_pairwise_DEA)
dev.off()

# =============================================================================
#  SECTION 2: Pairwise Regulon DEA on Dysfunctional vs Effector Cluster
# =============================================================================
# Pairwise Regulon DEA on Selected Annotated NK Clusters

auc_matrix <- read.csv("results/S03_results09_AUC_matrix.csv", row.names = 1,
                       stringsAsFactors = FALSE, check.names = FALSE)
pairwise_auc_matrix <- auc_matrix[colnames(NKG2C_gene_pairwise), ]
cluster_id          <- as.character(NKG2C_gene_pairwise$NK_cluster)

pairwise_regulon_DEA <- lapply(colnames(pairwise_auc_matrix), function(reg) {
  auc_reg_c4 <- pairwise_auc_matrix[cluster_id == "4", reg]
  auc_reg_c1 <- pairwise_auc_matrix[cluster_id == "1", reg]
  test      <- wilcox.test(auc_reg_c4, auc_reg_c1)
  data.frame(
    regulon    = reg,
    avg_log2FC = log2((mean(auc_reg_c4) + 1e-9) / (mean(auc_reg_c1) + 1e-9)),
    pval       = test$p.value,
    stringsAsFactors = FALSE)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    gene      = sub("\\(\\+\\)$", "", regulon),
    p_val_adj = p.adjust(pval, method = "fdr"),
    sig       = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC >=  1.00 ~ "Upregulated in Activated Stressed NK",
      p_val_adj < 0.05 & avg_log2FC <= -1.00 ~ "Upregulated in Mature Cytotoxic NK",
      TRUE ~ "NS")) %>%
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC))

pairwise_regulon_DEA_sig <- dplyr::filter(pairwise_regulon_DEA, sig != "NS")

write.csv(pairwise_regulon_DEA,
          "results/S05_results15_regulon_DEA_pairwise_dysfunctional_vs_effector.csv",
          row.names = FALSE)
write.csv(pairwise_regulon_DEA_sig,
          "results/S05_results16_regulon_DEA_pairwise_dysfunctional_vs_effector_sig.csv",
          row.names = FALSE)

# Visualize pairwise regulon DEA via Volcano Plot

df_regulon_up <- pairwise_regulon_DEA %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 1.00) %>%
  dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene)

df_regulon_down <- pairwise_regulon_DEA %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC < -1.00) %>%
  dplyr::arrange(p_val_adj, avg_log2FC) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene)

regulon_DEA_annotated_genes <- data.frame(
  gene = unique(c(df_regulon_up, df_regulon_down)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    regulon_DEA_color = dplyr::case_when(
      gene %in% df_regulon_up   ~ COL_UP,
      gene %in% df_regulon_down ~ COL_DOWN,
      TRUE ~ "grey80"))

df_volcano_pairwise_regulon_DEA <- pairwise_regulon_DEA %>%
  dplyr::left_join(
    regulon_DEA_annotated_genes[, c("gene", "regulon_DEA_color")],
    by = "gene") %>%
  dplyr::mutate(
    neglog10p    = -log10(p_val_adj + 1e-300),
    significance = dplyr::case_when(
      p_val_adj < 0.05 & avg_log2FC >=  1.00 ~ COL_UP,
      p_val_adj < 0.05 & avg_log2FC <= -1.00 ~ COL_DOWN,
      TRUE ~ COL_NS),
    point_col = dplyr::coalesce(regulon_DEA_color, significance),
    alpha_val = dplyr::if_else(!is.na(regulon_DEA_color), 0.85, 0.4)
  )
df_regulon_annotated <- dplyr::filter(df_volcano_pairwise_regulon_DEA, !is.na(regulon_DEA_color))

plot_volcano_pairwise_regulon_DEA <- ggplot2::ggplot(df_volcano_pairwise_regulon_DEA,
                                                  ggplot2::aes(x = avg_log2FC, y = neglog10p)) +
  ggplot2::geom_point(ggplot2::aes(colour = point_col, alpha = alpha_val),
                      size = 1.5, stroke = 0, shape = 16) +
  ggplot2::scale_colour_identity(
    guide  = ggplot2::guide_legend(title = "Regulon category",
                                   override.aes = list(size = 3)),
    breaks = c(COL_UP, COL_DOWN, COL_NS),
    labels = c("Top 10 most significant, up in Activated Stressed NK",
               "Top 10 most significant, up in Mature Cytotoxic NK",
               "Not significant")
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::geom_vline(xintercept = c(-1.00, 1.00),
                      linetype = "dashed", colour = "grey50") +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(
    data   = df_regulon_annotated,
    ggplot2::aes(colour = point_col),
    shape  = 16,
    size   = 3,
    stroke = 1
  ) +
  ggrepel::geom_text_repel(
    data = df_regulon_annotated,
    ggplot2::aes(
      label    = gene,
      colour   = point_col,
      fontface = "plain"
    ),
    size               = 2.5,
    max.overlaps       = Inf,
    show.legend        = FALSE,
    segment.size       = 0.3,
    segment.alpha      = 0.6,
    box.padding        = 0.6,
    seed               = 42,
    force              = 3
  ) +
  ggplot2::labs(
    title = "Differential Regulon Expression on Dysfunctional vs Effector NK Cluster",
    x     = "Regulon avg Log2FC (AUCell activity)",
    y     = expression(-log[10]~p[adj])
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

pdf("plots/S05_plot10_regulon_DEA_pairwise_dysfunctional_vs_effector_volcano_plot.pdf",
    width = 9, height = 7)
print(plot_volcano_pairwise_regulon_DEA)
dev.off()

# =============================================================================
#  Save Checkpoint
# =============================================================================
saveRDS(
  list(NKG2C_gated_obj     = NKG2C_gated_obj,
       NKG2C_gene_pairwise = NKG2C_gene_pairwise,
       cluster_id          = cluster_id,
       pairwise_DEA        = pairwise_DEA),
  "checkpoint/S05_pairwise_DEA.rds"
)
writeLines(capture.output(sessionInfo()), "session/S05_session_info.txt")
