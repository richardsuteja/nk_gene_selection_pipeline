# =============================================================================
#  ADAPT-NK TARGET SCREENING PIPELINE
#  Script 01 — Gating, Batch Correction & Leiden Clustering
#
#  Input : adata_all_nk_after_mapping.h5ad    (from Netskar et al., 2024)
#  Output: results/S01_results01_cells_per_source.csv
#          results/S01_results02_total_cells_by_cluster.csv
#          results/S01_results03_table_source_by_cluster.csv
#          plots/S01_plot01_UMAP_no_annotation.pdf
#          plots/S01_plot02_UMAP_by_source.pdf
#          plots/S01_plot03_UMAP_by_clusters.pdf
#          sanity/S01_sanity01_dropout_per_minimum_cells.csv
#          sanity/S01_sanity02_dropout_rate_minimum_cells.csv
#          sanity/S01_sanity03_cumulative_variance_per_pc.csv
#          sanity/S01_sanity04_resolution_k.param_sweep.csv
#          checkpoint/S01_gating_clustering.rds
#          session/S01_session_info.txt
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(reticulate)
use_condaenv("scvi", required=TRUE)
library(Seurat)
library(anndata)
library(cluster)
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(42)
options(future.globals.maxSize     = 8 * 1024^3)

dir.create("plots",      showWarnings = FALSE)
dir.create("results",    showWarnings = FALSE)
dir.create("sanity",     showWarnings = FALSE)
dir.create("checkpoint", showWarnings = FALSE)
dir.create("session",    showWarnings = FALSE)

# =============================================================================
#  SECTION 1: Load h5ad Raw Data
# =============================================================================
# Convert h5ad file to Seurat Object

h5ad <- read_h5ad("adata_all_nk_after_mapping.h5ad")
counts_mat <- t(h5ad$X)
colnames(counts_mat) <- h5ad$obs_names
rownames(counts_mat) <- h5ad$var_names
seurat_obj <- CreateSeuratObject(
  counts = counts_mat,
  min.cells = 3,
  min.features = 200,
  meta.data = h5ad$obs)
seurat_obj <- NormalizeData(seurat_obj)

# =============================================================================
#  SECTION 2: KLRC2 and TiNK Gating
# =============================================================================
# Gate on positive KLRC2 expression and restrict to TiNK

seurat_obj$site_group <- ifelse(grepl("_normal$", seurat_obj$source),
                                "TrNK",
                                "Not TrNK")
seurat_obj$sample_id <- paste(seurat_obj$patient,
                              seurat_obj$source, sep = "-")
KLRC2_expr <- FetchData(seurat_obj,
                        vars = c("KLRC2", "source", "site_group"),
                        layer = "data")
KLRC2_gated_cells <- rownames(KLRC2_expr)[KLRC2_expr$KLRC2 > 0 &
                                          KLRC2_expr$source != "PBMC" &
                                          KLRC2_expr$site_group != "TrNK"]
NKG2C_gated_obj <- subset(seurat_obj,
                          cells = KLRC2_gated_cells)
NKG2C_gated_obj@meta.data %>%
  dplyr::count(source, name = "n_cells") %>%
  write.csv("results/S01_results01_cells_per_source.csv", row.names = FALSE)

# =============================================================================
#  SECTION 3: Sample-level Quality Filter
# =============================================================================
# Sanity check to determine minimum number of cells to remove noise

qf_cell_count <- NKG2C_gated_obj@meta.data %>%
  dplyr::count(sample_id, name = "n_cells")
qf_meta_data <- NKG2C_gated_obj@meta.data %>%
  dplyr::distinct(sample_id, source)
qf_thresholds <- seq(2, 30, by = 2)
qf_titration_list <- vector("list", length(qf_thresholds))

for (i in seq_along(qf_thresholds)) {
  t <- qf_thresholds[i]
  qf_titration_list[[i]] <- qf_cell_count %>%
    dplyr::left_join(qf_meta_data, by = "sample_id") %>%
    dplyr::mutate(passes = n_cells >= t) %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(
      kept               = sum(passes),
      dropped            = sum(!passes),
      percentage_dropped = round(100 * sum(!passes) / dplyr::n(), 1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(threshold = t)
}
  qf_titration_long <- dplyr::bind_rows(qf_titration_list)
  qf_titration_wide <- qf_titration_long %>%
    dplyr::select(source, threshold, percentage_dropped) %>%
    tidyr::pivot_wider(
      names_from   = threshold,
      values_from  = percentage_dropped,
      names_prefix = "n"
    ) %>%
    dplyr::arrange(dplyr::desc(n2))
  qf_rate_change <- qf_titration_long %>%
    dplyr::group_by(threshold) %>%
    dplyr::summarise(
      mean_pct_dropped = round(mean(percentage_dropped), 1),
      .groups = "drop"
    ) %>%
    dplyr::arrange(threshold) %>%
    dplyr::mutate(
      change_from_prev = round(mean_pct_dropped - dplyr::lag(mean_pct_dropped), 1)
    )
write.csv(qf_titration_wide,"sanity/S01_sanity01_dropout_per_minimum_cells.csv",
          row.names = FALSE)
write.csv(qf_rate_change,"sanity/S01_sanity02_dropout_rate_minimum_cells.csv",
          row.names = FALSE)

# Filter at a threshold of 10 owing to large (~2x) drop in Sarcoma at 12

qf_passed_samples <- qf_cell_count %>%
  dplyr::filter(n_cells >= 10) %>%
  dplyr::pull(sample_id)
qf_passed_cells <- colnames(NKG2C_gated_obj)[
  NKG2C_gated_obj@meta.data$sample_id %in% qf_passed_samples]
NKG2C_gated_obj <- subset(NKG2C_gated_obj, cells = qf_passed_cells)

# =============================================================================
#  SECTION 4: scVI Batch Correction
# =============================================================================
# Identify 2000 highly variable genes then sanity estimate using PCA

NKG2C_gated_obj <- FindVariableFeatures(NKG2C_gated_obj,
                                         nfeatures = 2000,
                                         selection.method = "vst")
scvi_hvgs <- VariableFeatures(NKG2C_gated_obj)
raw_counts <- GetAssayData(NKG2C_gated_obj, assay = "RNA",
                           layer = "counts")[scvi_hvgs, ]

NKG2C_gated_obj <- ScaleData(NKG2C_gated_obj,
                             features      = scvi_hvgs,
                             verbose       = FALSE)
NKG2C_gated_obj <- RunPCA(NKG2C_gated_obj,
                          features  = scvi_hvgs,
                          npcs      = 50,
                          seed.use  = 42,
                          verbose   = FALSE)

pca_scree <- NKG2C_gated_obj[["pca"]]@stdev %>%
  {.^2} %>%
  {100 * . / sum(.)} %>%
  cumsum() %>%
  round(1) %>%
  {data.frame(PC = seq_along(.), cum_pct = .)}
write.csv(pca_scree,"sanity/S01_sanity03_cumulative_variance_per_pc.csv",
          row.names = FALSE)

# Set final n_latent at 45: a conservative default informed by a PCA-based
# rough estimate of effective dimension (90% variance at PC 38)

ad <- import("anndata")
scvi <- import("scvi")
scvi_object <- ad$AnnData(X = t(as.matrix(raw_counts)))
scvi_object$obs <- NKG2C_gated_obj@meta.data %>%
  dplyr::mutate(dplyr::across(where(is.factor), droplevels))
scvi$model$SCVI$setup_anndata(scvi_object, batch_key = "patient")
scvi$settings$seed <- as.integer(42)
vae <- scvi$model$SCVI(scvi_object, n_latent = as.integer(45))

# Train scVI model

vae$train(max_epochs = as.integer(400),
          batch_size = as.integer(512),
          early_stopping = TRUE,
          early_stopping_patience = as.integer(45))

# Extract scVI matrix

latent_mat <- vae$get_latent_representation()
rownames(latent_mat) <- colnames(NKG2C_gated_obj)
colnames(latent_mat) <- paste0("scvi_", seq_len(45))

NKG2C_gated_obj[["scvi"]] <- CreateDimReducObject (embeddings = latent_mat,
                                                   key = "scvi_",
                                                   assay = DefaultAssay(NKG2C_gated_obj))

# =============================================================================
#  SECTION 5: UMAP Visualization
# =============================================================================
# Visual confirmation of UMAP plot after scVI batch correction

set.seed(42)
NKG2C_gated_obj <- RunUMAP(NKG2C_gated_obj,
                            dims        = 1:45,
                            reduction   = "scvi",
                            n.neighbors = 15,
                            min.dist    = 0.05,
                            seed.use    = 42)

pdf("plots/S01_plot01_UMAP_no_annotation.pdf", width = 7, height = 6)
  print(DimPlot(NKG2C_gated_obj,
                reduction = "umap") +
          NoLegend())
dev.off()

pdf("plots/S01_plot02_UMAP_by_source.pdf", width = 7, height = 6)
  print(DimPlot(NKG2C_gated_obj,
                reduction = "umap",
                group.by = "source") +
    ggtitle("NK Cells by Source"))
dev.off()

# =============================================================================
#  SECTION 6: Leiden Clustering
# =============================================================================
# Sanity check sweep calculation of silhouette score on resolution and k.param

set.seed(42)
k_sweep <- seq(5, 50, by = 5)
res_sweep <- seq(0.1, 2.0, by = 0.1)

n_cells <- ncol(NKG2C_gated_obj)
idx_k   <- sample(n_cells, min(2000, n_cells))
dist_k  <- dist(latent_mat[idx_k, ])

sweep_results <- lapply(k_sweep, function(k) {
  message(sprintf("Sweeping k = %d...", k))
  sweep_obj <- FindNeighbors(NKG2C_gated_obj,
                           reduction = "scvi",
                           dims      = 1:45, 
                           k.param   = k,
                           verbose   = FALSE)
  res_list <- lapply(res_sweep, function(res) {
    sweep_obj <- FindClusters(sweep_obj,
                            resolution  = res,
                            algorithm   = 4,
                            random.seed = 42,
                            verbose     = FALSE)
    
    n_cluster_sweep <- nlevels(Idents(sweep_obj))
    clusters_int <- as.integer(Idents(sweep_obj))[idx_k]
    average_sweep_sil <- if (n_cluster_sweep > 1) {
      mean(cluster::silhouette(clusters_int, dist_k)[, 3])
    } else {
      NA_real_
    }

    data.frame(
      k_param        = k,
      resolution     = res,
      n_clusters     = n_cluster_sweep,
      avg_silhouette = round(average_sweep_sil, 4)
    )
  })
  dplyr::bind_rows(res_list)
})

grid_sweep_results  <- dplyr::bind_rows(sweep_results) %>%
  dplyr::arrange(desc(avg_silhouette))

write.csv(grid_sweep_results, "sanity/S01_sanity04_resolution_k.param_sweep.csv",
          row.names = FALSE)

# Building KNN from k.param 40 and perform Leiden Clustering on resolution 0.8
# The levels were selected by manual review of full sweep while balancing silhouette
# quality against per-source representation

NKG2C_gated_obj <- FindNeighbors(NKG2C_gated_obj,
                                  reduction = "scvi",
                                  dims      = 1:45,
                                  k.param   = 40)
NKG2C_gated_obj <- FindClusters(NKG2C_gated_obj,
                                 resolution  = 0.8,
                                 algorithm   = 4,
                                 random.seed = 42)
NKG2C_gated_obj$NK_cluster <- as.character(Idents(NKG2C_gated_obj))

cells_per_cluster <- NKG2C_gated_obj@meta.data %>%
  dplyr::count(NK_cluster, name = "total_cells")

write.csv(cells_per_cluster,
          "results/S01_results02_total_cells_by_cluster.csv",
          row.names = FALSE)

cells_source_per_cluster <- NKG2C_gated_obj@meta.data %>%
  dplyr::count(NK_cluster, source, name = "absolute_count") %>%
  dplyr::group_by(NK_cluster) %>%
  dplyr::mutate(pct = round(absolute_count / sum(absolute_count) * 100, 1)) %>%
  dplyr::arrange(NK_cluster, dplyr::desc(pct))

write.csv(cells_source_per_cluster, 
          "results/S01_results03_table_source_by_cluster.csv", 
          row.names = FALSE)

pdf("plots/S01_plot03_UMAP_by_clusters.pdf", width = 7, height = 6)
  print(DimPlot(NKG2C_gated_obj,
                reduction = "umap",
                group.by  = "NK_cluster",
                cols      = c("1" = "#E69F00",
                              "2" = "#56B4E9",
                              "3" = "#009E73",
                              "4" = "#D55E00",
                              "5" = "#999999"),
                label = TRUE,
                repel = TRUE,
                label.size = 5))
dev.off()

# =============================================================================
#  Save Checkpoint
# =============================================================================

saveRDS(NKG2C_gated_obj, "checkpoint/S01_gating_clustering.rds")
writeLines(capture.output(sessionInfo()), "session/S01_session_info.txt")