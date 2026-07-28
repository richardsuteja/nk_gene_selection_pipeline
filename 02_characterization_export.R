# =============================================================================
#  ADAPT-NK TARGET SCREENING PIPELINE
#  Script 02 — Cluster Characterization & pySCENIC Export
#
#  Input : checkpoint/S01_gating_clustering.rds
#  Output: results/S02_results04_DEA_general.csv
#          results/S02_results05_DEA_general_sig.csv
#          results/S02_results06_GSEA_hallmark_sig.csv
#          results/S02_results07_GSEA_immune_sig.csv
#          results/S02_results08_AUCell_Signature.csv
#          plots/S02_plot04_DEA_volcano_plot.pdf
#          plots/S02_plot05_GSEA_NES_heatmap.pdf
#          plots/S02_plot06_AUC_heatmap.pdf
#          checkpoint/S02_NKG2C_gated_for_pyscenic.h5ad
#          checkpoint/S02_characterization_export.rds
#          session/S02_session_info.txt
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(anndata)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

set.seed(42)

# ── Load checkpoint and General Sanity Check ──────────────────────────────────
NKG2C_gated_obj <- readRDS("checkpoint/S01_gating_clustering.rds")
message("Loaded: ", ncol(NKG2C_gated_obj), " cells, ",
        nrow(NKG2C_gated_obj), " genes")
DefaultAssay(NKG2C_gated_obj) <- "RNA"
Idents(NKG2C_gated_obj) <- "NK_cluster"

# ── Gene sets for GSEA ────────────────────────────────────────────────────────
gene_set_hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
gene_set_immunesig<- msigdbr(species = "Homo sapiens", category = "C7",
                             subcategory = "IMMUNESIGDB") %>%
  dplyr::filter(grepl("NK|CD8", gs_name, ignore.case = TRUE)) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# =============================================================================
#  SECTION 1: Differential Expression Analysis
# =============================================================================
# Perform DEA on cluster of interest vs other clusters

DEA_general <- FindAllMarkers(NKG2C_gated_obj,
                              only.pos        = FALSE,
                              test.use        = "wilcox",
                              min.pct         = 0.10,
                              logfc.threshold = 0,
                              verbose         = FALSE
                              ) %>%
  dplyr::arrange(cluster, p_val_adj, dplyr::desc(avg_log2FC))
DEA_general_sig <- DEA_general %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) >= 1) %>%
  dplyr::arrange(cluster, p_val_adj, dplyr::desc(avg_log2FC))

write.csv(DEA_general, "results/S02_results04_DEA_general.csv",
          row.names = FALSE)
write.csv(DEA_general_sig, "results/S02_results05_DEA_general_sig.csv",
          row.names = FALSE)

# Select 5 top genes each by p-value and Log2FC at each end, then
# visualize via volcano plot with selected annotations

COL_BOTH <- "#7B4F9E"
COL_UP   <- "#B2182B"
COL_DOWN <- "#2166AC"
COL_NS   <- "grey80"

cluster_id <- sort(unique(DEA_general$cluster))
volcano_panels <- lapply(cluster_id, function(cl) {
  
  df_DEA_up_padj <- DEA_general %>%
    dplyr::filter(cluster == cl, p_val_adj  <  0.05, avg_log2FC > 0.25) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(avg_log2FC)) %>%
    dplyr::slice_min(p_val_adj, n = 5) %>%
    dplyr::pull(gene)
    
  df_DEA_up_log2FC <- DEA_general %>%
    dplyr::filter(cluster == cl, p_val_adj  <  0.05, avg_log2FC > 0.25) %>%
    dplyr::arrange(dplyr::desc(avg_log2FC), p_val_adj) %>%
    dplyr::slice_max(avg_log2FC, n = 5) %>%
    dplyr::pull(gene)
    
  df_DEA_up_overlap <- intersect(df_DEA_up_padj, df_DEA_up_log2FC)
  
  df_DEA_down_padj <- DEA_general %>%
    dplyr::filter(cluster == cl, p_val_adj  <  0.05, avg_log2FC < -0.25) %>%
    dplyr::arrange(p_val_adj, avg_log2FC) %>%
    dplyr::slice_min(p_val_adj, n = 5) %>%
    dplyr::pull(gene)
    
  df_DEA_down_log2FC <- DEA_general %>%
    dplyr::filter(cluster == cl, p_val_adj  <  0.05, avg_log2FC < -0.25) %>%
    dplyr::arrange(avg_log2FC, p_val_adj) %>%
    dplyr::slice_min(avg_log2FC, n = 5) %>%
    dplyr::pull(gene)
    
  df_DEA_down_overlap <- intersect(df_DEA_down_padj, df_DEA_down_log2FC)
  
  DEA_annotated_genes <- data.frame(gene = unique(c(df_DEA_up_padj,
                                                    df_DEA_up_log2FC,
                                                    df_DEA_down_padj,
                                                    df_DEA_down_log2FC)),
                                    stringsAsFactors = FALSE) %>%
    
    dplyr::mutate(
      DEA_color = dplyr::case_when(
        gene %in% c(df_DEA_up_overlap, df_DEA_down_overlap) ~ COL_BOTH,
        gene %in% c(df_DEA_up_padj,    df_DEA_up_log2FC)    ~ COL_UP,
        gene %in% c(df_DEA_down_padj,  df_DEA_down_log2FC)  ~ COL_DOWN,
        TRUE ~ "grey80"))
  
  df_volcano_all <- DEA_general %>%
    dplyr::filter(cluster == cl) %>%
    dplyr::left_join(DEA_annotated_genes, by = "gene") %>%
    dplyr::mutate(neglog10p = pmin(-log10(p_val_adj + 1e-300),
                       ifelse(cl == "5", 150, 300)),
                  significance = dplyr::case_when(
                    p_val_adj < 0.05 & avg_log2FC >=  0.25 ~ COL_UP,
                    p_val_adj < 0.05 & avg_log2FC <= -0.25 ~ COL_DOWN,
                    TRUE                                   ~ COL_NS),
                  point_col = dplyr::coalesce(DEA_color, significance),
                  alpha_val = dplyr::if_else(!is.na(DEA_color), 0.85, 0.4))
  df_annotated <- dplyr::filter(df_volcano_all, !is.na(DEA_color))
  
  ggplot2::ggplot(df_volcano_all,
                  ggplot2::aes(x = avg_log2FC, y = neglog10p)) +
  ggplot2::geom_point(
    ggplot2::aes(colour = point_col, alpha  = alpha_val),
    size   = 1.5, stroke = 0, shape  = 16
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
    
    ggplot2::geom_vline(xintercept = c(-0.25, 0.25),
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
      force              = ifelse(cl == "5", 40, 8),
      force_pull         = ifelse(cl == "5", 0.1, 1)
    ) +

    ggplot2::labs(
      title = paste0("Cluster ", cl, " vs Other Clusters"),
      x     = "Log2FC",
      y     = expression(-log[10]~padj)
    ) +    
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold")
    )
})

names(volcano_panels) <- cluster_id
plot_volcano_collage <- patchwork::wrap_plots(
  volcano_panels[c("1", "2", "3", "4", "5")],
  nrow   = 1,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Differential Expression Across NK Clusters",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )
  )

pdf("plots/S02_plot04_DEA_volcano_plot.pdf",
    width  = 5 * length(cluster_id),
    height = 7)
print(plot_volcano_collage)
dev.off()

# =============================================================================
#  SECTION 2: Gene Set Enrichment Analysis
# =============================================================================
# Perform GSEA on every NK clusters

gsea_list_hallmark <- list()
for (grp in sort(unique(NKG2C_gated_obj$NK_cluster))) {
  
  gsea_rank_hallmark <- DEA_general %>%
    dplyr::filter(cluster == grp) %>%
    dplyr::mutate(
      neglog10p = -log10(p_val + 1e-300),
      # p_val == 0 (floating-point underflow, common for very strongly DE
      # genes at this sample size) collapses many genes to the same capped
      # neglog10p -- break that tie by |avg_log2FC| rank instead of leaving
      # it to fgsea's internal (arbitrary) tie handling. percent_rank() is
      # bounded to [0,1], so this can only reorder genes already tied at the
      # cap; it never crosses into where genuinely different p-values sit.
      tie_break = dplyr::if_else(p_val == 0,
                                 dplyr::percent_rank(abs(avg_log2FC)),
                                 0),
      rank_stat = sign(avg_log2FC) * (neglog10p + tie_break)
    ) %>%
    dplyr::arrange(dplyr::desc(rank_stat)) %>%
    { setNames(.$rank_stat, .$gene) }
  
  gsea_list_hallmark[[grp]] <- fgsea(
    pathways    = gene_set_hallmark,
    stats       = gsea_rank_hallmark,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 1e-10,
    nproc       = 8
  ) %>%
    dplyr::mutate(NK_cluster = grp)
}

gsea_hallmark <- dplyr::bind_rows(gsea_list_hallmark)
gsea_hallmark_sig <- gsea_hallmark %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(NK_cluster, padj)

gsea_list_immune <- list()
for (grp in sort(unique(NKG2C_gated_obj$NK_cluster))) {
  
  gsea_rank_immune <- DEA_general %>%
    dplyr::filter(cluster == grp) %>%
    dplyr::mutate(
      neglog10p = -log10(p_val + 1e-300),
      tie_break = dplyr::if_else(p_val == 0,
                                 dplyr::percent_rank(abs(avg_log2FC)),
                                 0),
      rank_stat = sign(avg_log2FC) * (neglog10p + tie_break)
    ) %>%
    dplyr::arrange(dplyr::desc(rank_stat)) %>%
    { setNames(.$rank_stat, .$gene) }
  
  gsea_list_immune[[grp]] <- fgsea(
    pathways    = gene_set_immunesig,
    stats       = gsea_rank_immune,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 1e-10,
    nproc       = 8
  ) %>%
    dplyr::mutate(NK_cluster = grp)
}

gsea_immune <- dplyr::bind_rows(gsea_list_immune)
gsea_immune_sig <- gsea_immune %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(NK_cluster, padj)

write.csv(gsea_hallmark_sig %>%
    dplyr::mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")),
  "results/S02_results06_GSEA_hallmark_sig.csv",
  row.names = FALSE)
write.csv(gsea_immune_sig %>%
            dplyr::mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")),
          "results/S02_results07_GSEA_immune_sig.csv",
          row.names = FALSE)

# Select 5 top and bottom enrichment pathways and clean their labels

gsea_annotated_hallmark <- gsea_hallmark_sig %>%
  dplyr::mutate(direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
  dplyr::group_by(NK_cluster, direction) %>%
  dplyr::slice_min(padj, n = 5) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    pathway_clean = pathway %>%
      gsub("^HALLMARK_",   "H: ",  .) %>%
      gsub("^GSE[0-9]+_", "",     .) %>%
      gsub("_", " ",              .),
    pathway_short = ifelse(
      nchar(pathway_clean) > 55,
      paste0(substr(pathway_clean, 1, 52), "..."),
      pathway_clean))%>%
  dplyr::arrange(direction, dplyr::desc(NES)) %>%
  dplyr::mutate(pathway_short = factor(pathway_short,
                                       levels = unique(pathway_short)))

gsea_annotated_immune <- gsea_immune_sig %>%
  dplyr::mutate(direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>%
  dplyr::group_by(NK_cluster, direction) %>%
  dplyr::slice_min(padj, n = 5) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    pathway_clean = pathway %>%
      gsub("^HALLMARK_",   "H: ",  .) %>%
      gsub("^GSE[0-9]+_", "",     .) %>%
      gsub("_", " ",              .),
    pathway_short = ifelse(
      nchar(pathway_clean) > 55,
      paste0(substr(pathway_clean, 1, 52), "..."),
      pathway_clean
    )
  ) %>%
  dplyr::arrange(direction, dplyr::desc(NES)) %>%
  dplyr::mutate(pathway_short = factor(pathway_short,
                                       levels = unique(pathway_short)))

# Visualize via NES heat map

gsea_heatmap <- function(annotated_df, full_gsea_df, plot_title) {
  gsea_cluster_vector <- annotated_df %>%
    dplyr::mutate(pathway_short = as.character(pathway_short)) %>%
    dplyr::group_by(pathway, pathway_short) %>%
    dplyr::slice_min(as.integer(as.character(NK_cluster)), n = 1,
                     with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(pathway, pathway_short, NK_cluster, NES)
  
  path_order <- gsea_cluster_vector %>%
    dplyr::arrange(as.integer(as.character(NK_cluster)), dplyr::desc(NES)) %>%
    dplyr::pull(pathway_short)
  
  pathway_cluster <- full_gsea_df %>%
    dplyr::filter(pathway %in% gsea_cluster_vector$pathway) %>%
    dplyr::left_join(
      dplyr::select(gsea_cluster_vector, pathway, pathway_short),
      by = "pathway") %>%
    dplyr::distinct(NK_cluster, pathway, .keep_all = TRUE) %>%
    dplyr::mutate(
      NK_cluster    = factor(NK_cluster,
        levels = as.character(sort(as.integer(unique(NK_cluster))))),
      pathway_short = factor(pathway_short, levels = rev(path_order)))
  
  ggplot2::ggplot(pathway_cluster,
                  ggplot2::aes(x = NK_cluster, y = pathway_short)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = NES),
      colour    = "white",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_gradient2(
      low      = COL_DOWN,
      mid      = "white",
      high     = COL_UP,
      midpoint = 0,
      limits   = c(-2.5, 2.5),
      oob      = scales::squish,
      name     = "NES"
    ) +
    ggplot2::geom_text(
      data   = dplyr::filter(pathway_cluster, padj < 0.05),
      ggplot2::aes(label = "*"),
      colour = "white",
      size   = 5,
      vjust  = 0.75
    ) +
    ggplot2::labs(title = plot_title, x = "NK Cluster", y = NULL) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      axis.text.y      = ggplot2::element_text(size = 7),
      axis.text.x      = ggplot2::element_text(size = 9),
      axis.ticks       = ggplot2::element_blank(),
      axis.line        = element_blank()
    )
}

plot_gsea_hallmark  <- gsea_heatmap(gsea_annotated_hallmark, gsea_hallmark,
                                   "Hallmark Gene Sets")
plot_gsea_immunesig <- gsea_heatmap(gsea_annotated_immune, gsea_immune,
                                   "ImmuneSigDB NK/CD8+ Gene Sets")

plot_nes_collage <- patchwork::wrap_plots(
  list(plot_gsea_hallmark, plot_gsea_immunesig),
  nrow   = 1,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title    = "GSEA NES Heatmap Across NK Clusters",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 18,
                                                              face = "bold")))

pdf("plots/S02_plot05_GSEA_NES_heatmap.pdf",
    width  = 16,
    height = 4 + max(
      dplyr::n_distinct(gsea_annotated_hallmark$pathway_short),
      dplyr::n_distinct(gsea_annotated_immune$pathway_short)
    ) * 0.4)
print(plot_nes_collage)
dev.off()

# =============================================================================
#  SECTION 3: Visualize Pre-Calculated AUCell Signature
# =============================================================================
# Summarizing pre-calculated AUCell scores attached from Netskar et al. (2024)

auc_cell <- NKG2C_gated_obj@meta.data %>%
  dplyr::select(NK_cluster, dplyr::ends_with("_auc")) %>%
  dplyr::mutate(NK_cluster = as.character(NK_cluster))

auc_sig_tests <- do.call(rbind, lapply(setdiff(colnames(auc_cell),"NK_cluster"),
                                       function(sig) {
  do.call(rbind, lapply(sort(unique(NKG2C_gated_obj$NK_cluster)), function(cl) {
    data.frame(
      signature_raw = sig,
      NK_cluster    = cl,
      p_value       = suppressWarnings(
        wilcox.test(auc_cell[[sig]][auc_cell$NK_cluster == cl],
                    auc_cell[[sig]][auc_cell$NK_cluster != cl])$p.value))}))})) %>%
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

auc_name_clean <- data.frame(
  signature_raw = setdiff(colnames(auc_cell), "NK_cluster"),
  signature     = setdiff(colnames(auc_cell), "NK_cluster") %>%
    gsub("_auc$",     "",       .) %>%
    gsub("^HALLMARK_","H: ",    .) %>%
    gsub("^KEGG_",    "KEGG: ", .) %>%
    gsub("^GOBP_",    "GO: ",   .) %>%
    gsub("_",         " ",      .) %>%
    toupper(),
  stringsAsFactors = FALSE
)

mean_auc <- auc_cell %>%
  dplyr::group_by(NK_cluster) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 ~ mean(.x))) %>%
  tibble::column_to_rownames("NK_cluster")

write.csv(mean_auc, "results/S02_results08_AUCell_Signature.csv")

auc_z  <- scale(as.matrix(mean_auc))

auc_long <- as.data.frame(t(auc_z)) %>%
  tibble::rownames_to_column("signature_raw") %>%
  dplyr::mutate(
    signature = auc_name_clean$signature[
      match(signature_raw, auc_name_clean$signature_raw)
    ]
  ) %>%
  dplyr::select(-signature_raw) %>%
  tidyr::pivot_longer(-signature,
                      names_to  = "NK_cluster",
                      values_to = "z_score") %>%
  dplyr::left_join(
    auc_sig_tests %>%
      dplyr::mutate(NK_cluster = as.character(NK_cluster)) %>%
      dplyr::left_join(auc_name_clean, by = "signature_raw") %>%
      dplyr::select(signature, NK_cluster, padj),
    by = c("signature", "NK_cluster")
  ) %>%
  dplyr::mutate(
    NK_cluster = factor(NK_cluster,
                        levels = as.character(sort(as.integer(unique(NK_cluster)))))
  )

sig_order <- auc_long %>%
  dplyr::group_by(signature) %>%
  dplyr::slice_max(z_score, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(as.integer(as.character(NK_cluster)), dplyr::desc(z_score)) %>%
  dplyr::pull(signature)

auc_long <- auc_long %>%
  dplyr::mutate(signature = factor(signature, levels = rev(sig_order)))

plot_auc_heatmap <- ggplot2::ggplot(auc_long,
  ggplot2::aes(x = NK_cluster, y = signature)
) +
  ggplot2::geom_tile(ggplot2::aes(fill = z_score)) +
  ggplot2::scale_fill_gradient2(
    low      = COL_DOWN,
    mid      = "white",
    high     = COL_UP,
    midpoint = 0,
    limits   = c(-2, 2),
    oob      = scales::squish,
    name     = "AUC\n(z-score)"
  ) +
  ggplot2::geom_text(
    data   = dplyr::filter(auc_long, padj < 0.05),
    ggplot2::aes(label = "*"),
    colour = "white",
    size   = 5,
    vjust  = 0.75
  ) +
  ggplot2::labs(
    title = "AUCell Signatures Across NK Clusters",
    x     = "NK Cluster",
    y     = NULL
  ) +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    plot.title  = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 7),
    axis.ticks  = ggplot2::element_blank(),
    axis.line   = ggplot2::element_blank()
  )

pdf("plots/S02_plot06_AUC_heatmap.pdf",
    width  = 8,
    height = 3 + dplyr::n_distinct(auc_long$signature) * 0.3)
print(plot_auc_heatmap)
dev.off()

# =============================================================================
#  SECTION 4: Exporting h5ad for pySCENIC (Python)
# =============================================================================

adata_out <- AnnData(
  X = t(as.matrix(
    GetAssayData(NKG2C_gated_obj, assay = "RNA", layer = "counts"))),
  obs = NKG2C_gated_obj@meta.data)
adata_out$write_h5ad("checkpoint/S02_NKG2C_gated_for_pyscenic.h5ad")

# =============================================================================
#  Save Checkpoint
# =============================================================================

saveRDS(
  list(NKG2C_gated_obj = NKG2C_gated_obj),
  "checkpoint/S02_characterization_export.rds")
writeLines(capture.output(sessionInfo()), "session/S02_session_info.txt")