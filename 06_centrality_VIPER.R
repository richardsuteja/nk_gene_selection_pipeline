# =============================================================================
#  ADAPT-NK TARGET SCREENING PIPELINE
#  Script 06 — GRN Centrality + VIPER Activity Inference + amiRNA Target Ranking
#
#  Input : checkpoint/S05_pairwise_DEA.rds
#          checkpoint/S03_adjacencies.csv
#  Output: results/S06_results17_VIPER_result.csv
#          results/S06_results18_amiRNA_knockdown_candidates_VIPER.csv
#          plots/S06_plot11_VIPER_volcano_dysfunctional_vs_effector.pdf
#          plots/S06_plot12_VIPER_candidates_scatter.pdf
#          checkpoint/S06_NKG2C_for_ORACLE.h5ad
#          checkpoint/S06_centrality_VIPER.rds
#          session/S06_session_info.txt
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(viper)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(matrixStats)
library(anndata)
set.seed(42)

# ── Load checkpoint ───────────────────────────────────────────────────────────
checkpoint_S05      <- readRDS("checkpoint/S05_pairwise_DEA.rds")
NKG2C_gated_obj     <- checkpoint_S05$NKG2C_gated_obj
NKG2C_gene_pairwise <- checkpoint_S05$NKG2C_gene_pairwise
cluster_id          <- checkpoint_S05$cluster_id
pairwise_DEA        <- checkpoint_S05$pairwise_DEA
rm(checkpoint_S05)

COL_UP   <- "#B2182B"
COL_DOWN <- "#2166AC"
COL_NS   <- "grey80"

# =============================================================================
#  SECTION 1: Build Directed Weighted Graph Based on Filtered pySCENIC GRN
# =============================================================================

adjacencies_pyscenic <- read.csv("checkpoint/S03_adjacencies.csv",
                               stringsAsFactors = FALSE) %>%
  dplyr::group_by(target) %>%
  dplyr::slice_max(importance, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()
tf_adjacencies_universe <- unique(adjacencies_pyscenic$TF)

weighted_graph <- igraph::graph_from_data_frame(
  d        = adjacencies_pyscenic[, c("TF", "target", "importance")],
  directed = TRUE)
# graph_from_data_frame already attaches the third column as edge attribute
# `importance`; invert it so a strong GRNBoost2 edge becomes a short path.
distance_to_weight <- 1 / igraph::E(weighted_graph)$importance

df_weighted_betweenness <- data.frame(
  gene        = igraph::V(weighted_graph)$name,
  betweenness = igraph::betweenness(weighted_graph, directed = TRUE,
                                    weights = distance_to_weight,
                                    normalized = TRUE),
  stringsAsFactors = FALSE)
df_weighted_betweenness$TF <- df_weighted_betweenness$gene %in% tf_adjacencies_universe

# =============================================================================
#  SECTION 2: TF Activity Inference via VIPER
# =============================================================================
# Build VIPER network

pyscenic_3col_tmp <- tempfile(fileext = ".tsv")
write.table(
  adjacencies_pyscenic[, c("TF", "target", "importance")],
  file = pyscenic_3col_tmp,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
expr_matrix <- as.matrix(Seurat::GetAssayData(NKG2C_gene_pairwise, layer = "data"))

regulon_VIPER <- aracne2regulon(
  afile  = pyscenic_3col_tmp,
  eset   = expr_matrix,
  format = "3col")
unlink(pyscenic_3col_tmp)

# Welch's t-test signature (C4 vs C1)

compute_signature <- function(expr_matrix, group, group1 = "4", group2 = "1") {
  x1 <- expr_matrix[, group == group1]
  x2 <- expr_matrix[, group == group2]

  m1 <- rowMeans(x1); m2 <- rowMeans(x2)
  v1 <- rowVars(x1);  v2 <- rowVars(x2)
  n1 <- ncol(x1);     n2 <- ncol(x2)

  se    <- sqrt(v1 / n1 + v2 / n2)
  tstat <- (m1 - m2) / se
  tstat[!is.finite(tstat)] <- 0
  names(tstat) <- rownames(expr_matrix)
  tstat
}

signature <- compute_signature(expr_matrix, cluster_id, group1 = "4", group2 = "1")

mrs <- msviper(signature, regulon_VIPER, verbose = TRUE)

VIPER_res <- data.frame(
  gene = names(mrs$es$nes),
  nes  = mrs$es$nes,
  pval = mrs$es$p.value,
  stringsAsFactors = FALSE
) %>% dplyr::arrange(pval)

write.csv(VIPER_res, "results/S06_results17_VIPER_result.csv", row.names = FALSE)

# Visualize via volcano plot

df_VIPER_up <- VIPER_res %>%
  dplyr::filter(pval < 0.05, nes > 0) %>%
  dplyr::arrange(pval, dplyr::desc(nes)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene)

df_VIPER_down <- VIPER_res %>%
  dplyr::filter(pval < 0.05, nes < 0) %>%
  dplyr::arrange(pval, nes) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene)

VIPER_annotated_genes <- data.frame(
  gene = unique(c(df_VIPER_up, df_VIPER_down)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    VIPER_color = dplyr::case_when(
      gene %in% df_VIPER_up   ~ COL_UP,
      gene %in% df_VIPER_down ~ COL_DOWN,
      TRUE ~ "grey80"))

df_volcano_VIPER <- VIPER_res %>%
  dplyr::left_join(VIPER_annotated_genes[, c("gene", "VIPER_color")], by = "gene") %>%
  dplyr::mutate(
    neglog10p    = -log10(pval + 1e-300),
    significance = dplyr::case_when(
      pval < 0.05 & nes > 0 ~ COL_UP,
      pval < 0.05 & nes < 0 ~ COL_DOWN,
      TRUE ~ COL_NS),
    point_col = dplyr::coalesce(VIPER_color, significance),
    alpha_val = dplyr::if_else(!is.na(VIPER_color), 0.85, 0.4)
  )
df_VIPER_annotated <- dplyr::filter(df_volcano_VIPER, !is.na(VIPER_color))
plot_volcano_VIPER <- ggplot2::ggplot(df_volcano_VIPER,
                ggplot2::aes(x = nes, y = neglog10p)) +
  ggplot2::geom_point(ggplot2::aes(colour = point_col, alpha = alpha_val),
    size = 1.5, stroke = 0, shape = 16) +
  ggplot2::scale_colour_identity(
    guide  = ggplot2::guide_legend(title = "VIPER TF category",
                                   override.aes = list(size = 3)),
    breaks = c(COL_UP, COL_DOWN, COL_NS),
    labels = c("Top 10 significant, up in Activated Stressed NK",
               "Top 10 significant, up in Mature Cytotoxic NK",
               "Not significant")
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::geom_hline(yintercept = -log10(0.05),
                      linetype = "dashed", colour = "grey50") +
  ggplot2::geom_point(
    data = df_VIPER_annotated,
    ggplot2::aes(colour = point_col),
    shape = 16, size = 3, stroke = 1
  ) +
  ggrepel::geom_text_repel(
    data = df_VIPER_annotated,
    ggplot2::aes(label = gene, colour = point_col),
    size = 2.5, max.overlaps = Inf, show.legend = FALSE,
    segment.size = 0.3, segment.alpha = 0.6, box.padding = 0.6,
    seed = 42, force = 3
  ) +
  ggplot2::labs(
    title = "VIPER TF Activity on Dysfunctional vs Effector NK Cluster",
    x     = "VIPER NES (C4 Exhausted / C1 Cytotoxic)",
    y     = expression(-log[10]~italic(p))
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

pdf("plots/S06_plot11_VIPER_volcano_dysfunctional_vs_effector.pdf",
    width = 9, height = 7)
print(plot_volcano_VIPER)
dev.off()

# Filter for candidate genes

candidates_amiRNA_VIPER <- df_weighted_betweenness %>%
  dplyr::left_join(VIPER_res, by = "gene") %>%
  dplyr::filter(TF, pval < 0.05, nes > 0) %>%
  dplyr::arrange(dplyr::desc(betweenness)) %>%
  dplyr::mutate(rank = dplyr::row_number())

write.csv(candidates_amiRNA_VIPER,
          "results/S06_results18_amiRNA_knockdown_candidates_VIPER.csv",
          row.names = FALSE)

# =============================================================================
#  SECTION 3: Post-VIPER Candidate Scatter
# =============================================================================
# Betweenness centrality vs pairwise log2FC; y-axis is centrality, not
# -log10(p), so this is a scatter rather than a volcano plot.

df_scatter_VIPER <- pairwise_DEA %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::select(gene, avg_log2FC, p_val_adj) %>%
  dplyr::inner_join(df_weighted_betweenness, by = "gene") %>%
  dplyr::mutate(
    is_VIPER_candidate = gene %in% candidates_amiRNA_VIPER$gene,
    point_col = dplyr::if_else(avg_log2FC > 0, COL_UP, COL_DOWN),
    alpha_val = dplyr::if_else(is_VIPER_candidate, 0.85, 0.4)
  )

df_scatter_VIPER_labels <- dplyr::filter(df_scatter_VIPER, is_VIPER_candidate)

plot_VIPER_scatter <- ggplot2::ggplot(df_scatter_VIPER,
                ggplot2::aes(x = avg_log2FC, y = betweenness)) +
  ggplot2::geom_point(
    ggplot2::aes(colour = point_col, alpha = alpha_val),
    size = 1.5, stroke = 0, shape = 16
  ) +
  ggplot2::scale_colour_identity(
    guide  = ggplot2::guide_legend(title = "",
                                   override.aes = list(size = 3, alpha = 1)),
    breaks = c(COL_UP, COL_DOWN),
    labels = c("Up in Activated Stressed NK",
               "Up in Mature Cytotoxic NK")
  ) +
  ggplot2::scale_alpha_identity() +
  ggplot2::geom_point(
    data = df_scatter_VIPER_labels,
    ggplot2::aes(colour = point_col),
    shape = 16, size = 3, stroke = 1,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = df_scatter_VIPER_labels,
    ggplot2::aes(label = gene, colour = point_col),
    size = 2.5, max.overlaps = Inf, show.legend = FALSE,
    segment.size = 0.3, segment.alpha = 0.6, box.padding = 0.6,
    seed = 42, force = 25
  ) +
  ggplot2::labs(
    title    = "VIPER-Selected amiRNA Candidates",
    x        = "Average Log2FC (C4 Exhausted / C1 Cytotoxic)",
    y        = "Weighted Betweenness Centrality"
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

pdf("plots/S06_plot12_VIPER_candidates_scatter.pdf",
    width = 9, height = 7)
print(plot_VIPER_scatter)
dev.off()

# =============================================================================
#  SECTION 4: Export h5ad for ORACLE (Python)
# =============================================================================
# Export for ORACLE single, pairwise, and fully-combined KO simulation

NKG2C_for_oracle <- NKG2C_gated_obj

pyscenic_cols <- grep("^[A-Z].*\\(\\+\\)$", colnames(NKG2C_for_oracle@meta.data),
                    value = TRUE)
hypoxia_col <- grep("HYPOX", colnames(NKG2C_for_oracle@meta.data),
                    value = TRUE, ignore.case = TRUE)[1]
obs_cols <- c("NK_cluster", "NK_cluster_annotated", hypoxia_col,
              head(pyscenic_cols, 50))

metadata_out <- NKG2C_for_oracle@meta.data[, obs_cols, drop = FALSE] %>%
  dplyr::mutate(NK_cluster = as.character(NK_cluster))

embeddings_out <- list(X_umap = Embeddings(NKG2C_for_oracle, "umap"))
if ("scvi" %in% Reductions(NKG2C_for_oracle)) {
  embeddings_out$X_scvi <- Embeddings(NKG2C_for_oracle, "scvi")
}

adata_out <- AnnData(
  X    = t(as.matrix(
    GetAssayData(NKG2C_for_oracle, assay = "RNA", layer = "counts")
  )),
  obs  = metadata_out,
  obsm = embeddings_out
)
adata_out$write_h5ad("checkpoint/S06_NKG2C_for_ORACLE.h5ad")

# =============================================================================
#  Save Checkpoint
# =============================================================================

saveRDS(
  list(NKG2C_gated_obj         = NKG2C_gated_obj,
       NKG2C_gene_pairwise     = NKG2C_gene_pairwise,
       weighted_graph          = weighted_graph,
       regulon_VIPER           = regulon_VIPER,
       mrs                     = mrs,
       candidates_amiRNA_VIPER = candidates_amiRNA_VIPER,
       df_scatter_VIPER        = df_scatter_VIPER),
  "checkpoint/S06_centrality_VIPER.rds")
writeLines(capture.output(sessionInfo()), "session/S06_session_info.txt")
