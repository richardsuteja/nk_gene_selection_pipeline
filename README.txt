ADAPT-NK TARGET SCREENING PIPELINE
====================================

PIPELINE FLOW
====================================

[VALIDATION] steps contextualize the target-selection. The final
NR4A1/ATF3 targets do not depend on these steps' own analytical conclusions

  Script 01 - Gating & Clustering
  (KLRC2/TiNK gating, scVI batch correction, Leiden clustering)
        |
        v
  --------------------------
  |                        |
  v                        v
  Script 02                Script 03
  (DEA + GSEA +            (pySCENIC: GRNBoost2 + cisTarget
   precalculated            + AUCell regulon activity)
   AUCell visualization)
  |                        |
  --------------------------
        |
        v
  Script 04 - Annotation
  (merge pySCENIC regulons, assign cluster identities)
        |
        v
  Script 05 - Pairwise DEA                            [VALIDATION]
  (gene + regulon DE, C4 vs C1)
        |
        v
  Script 06 - Centrality + VIPER
  (GRN betweenness centrality + VIPER TF activity inference)
        |
        v
  19 candidate TFs
        |
        v
  Literature Validation
  (systematic in vivo NK / CD8+ KO-KD search per candidate;
   see Validation.docx)
        |
        v
  FINAL TARGETS: NR4A1 + ATF3
  (the only two candidates with direct causal in vivo NK evidence)
        |
        v
  Script 07 - Oracle KO Simulation                     [VALIDATION]
  (in silico edit simulation: SMAD4, NR4A1, ATF3, and the combined
  SMAD4+NR4A1+ATF3 construct)


SCRIPT SUMMARY
====================================

01_gating_clustering.R
    Input:  adata_all_nk_after_mapping.h5ad
    Output: checkpoint/S01_gating_clustering.rds
            (+ results/, plots/, sanity/ -- see script header)
    Loads the NK atlas, gates on KLRC2+/TiNK, batch-corrects with scVI,
    clusters with Leiden.

02_characterization_export.R
    Input:  checkpoint/S01_gating_clustering.rds
    Output: checkpoint/S02_characterization_export.rds
            checkpoint/S02_NKG2C_gated_for_pyscenic.h5ad
            (+ results/, plots/ -- see script header)
    Runs one-vs-rest DEA and GSEA per cluster, summarizes precalculated
    AUCell signatures, exports an h5ad for pySCENIC. 

03_pyscenic.ipynb (Python)
    Input:  checkpoint/S02_NKG2C_gated_for_pyscenic.h5ad
    Output: checkpoint/S03_adjacencies.csv
            checkpoint/S03_regulons.csv
            results/S03_results09_AUC_matrix.csv
            results/S03_results10_DEA_regulon.csv
            checkpoint/S03_NKG2C_gated_after_pyscenic.h5ad
    GRNBoost2 (curated TF list, Lambert et al. 2018) infers the GRN;
    cisTarget motif-validates it into regulons; AUCell scores per-cell
    regulon activity; one-vs-rest DE on that activity per cluster.

04_annotation.R
    Input:  checkpoint/S02_characterization_export.rds
            results/S03_results09_AUC_matrix.csv
            results/S03_results10_DEA_regulon.csv
    Output: results/S04_results11_DEA_regulon_sig.csv
            results/S04_results12_source_by_annotated_cluster.csv
            plots/S04_plot07_DEA_regulon_volcano_plot.pdf
            plots/S04_plot08_UMAP_annotated.pdf
            checkpoint/S04_annotation.rds
    Merges pySCENIC regulon scores into the Seurat object, assigns
    human-readable cluster identities (NK_cluster_annotated).

05_pairwise_DEA.R  [VALIDATION]
    Input:  checkpoint/S04_annotation.rds
            results/S03_results09_AUC_matrix.csv
    Output: results/S05_results13_DEA_pairwise_dysfunctional_vs_effector.csv
            results/S05_results14_DEA_pairwise_dysfunctional_vs_effector_sig.csv
            results/S05_results15_regulon_DEA_pairwise_dysfunctional_vs_effector.csv
            results/S05_results16_regulon_DEA_pairwise_dysfunctional_vs_effector_sig.csv
            plots/S05_plot09_DEA_pairwise_dysfunctional_vs_effector_volcano_plot.pdf
            plots/S05_plot10_regulon_DEA_pairwise_dysfunctional_vs_effector_volcano_plot.pdf
            checkpoint/S05_pairwise_DEA.rds
    Gene- and regulon-level DE, C4 (Activated Stressed) vs C1 (Mature
    Cytotoxic) specifically.

06_centrality_VIPER.R
    Input:  checkpoint/S05_pairwise_DEA.rds
            checkpoint/S03_adjacencies.csv
    Output: results/S06_results17_VIPER_result.csv
            results/S06_results18_amiRNA_knockdown_candidates_VIPER.csv
            plots/S06_plot11_VIPER_volcano_dysfunctional_vs_effector.pdf
            plots/S06_plot12_VIPER_candidates_scatter.pdf
            checkpoint/S06_NKG2C_for_ORACLE.h5ad
            checkpoint/S06_centrality_VIPER.rds
    Builds the GRN as a directed graph, computes betweenness centrality,
    infers TF activity with VIPER (regulon enrichment on the C4-vs-C1
    signature, using the C1/C4 subset and cluster labels from Script 05's
    checkpoint).

Literature Validation
    Systematic search (see Validation.docx) for in vivo NK or CD8+ T
    cell KO/KD evidence per candidate. Only NR4A1 and ATF3 have direct
    causal in vivo NK evidence; the rest are underpowered or contradicted by 
    literature

07_oracle.ipynb (Python)  [VALIDATION]
    Input:  checkpoint/S06_NKG2C_for_ORACLE.h5ad
            checkpoint/S03_adjacencies.csv
    Output: results/S07_results19_ORACLE_KO_results.csv
            results/S07_results20_synergy_decomposition.csv
            plots/S07_plot13_single_pairwise_synergy_heatmap.pdf
            plots/S07_plot14_triple_synergy_waterfall.pdf
            checkpoint/S07_NKG2C_after_ORACLE.h5ad
            results/S07_results21_ORACLE_KO_results_hypoxia.csv
            results/S07_results22_synergy_decomposition_hypoxia.csv
            results/S07_results23_uniform_vs_hypoxia_comparison.csv
            plots/S07_plot15_single_pairwise_synergy_heatmap_hypoxia.pdf
            plots/S07_plot16_triple_synergy_waterfall_hypoxia.pdf
    In silico corroboration of the already-selected targets. CellOracle
    simulates every single, pairwise and triple combination of SMAD4, NR4A1
    and ATF3 under two scenarios ("ideal" = expression driven to zero;
    "realistic" = 50% of the baseline mean retained)

    Two separate passes:
    (1) Main flow, uniform KO applied to every cell, no hypoxia-gating.
    (2) Hypoxia-regressed section, reruns the same conditions with NR4A1/ATF3
        perturbed on a continuous per-cell basis: linearly interpolated
        between each cell's own baseline and the scenario target, weighted
        by the attached hypoxia AUCell score

ENVIRONMENT SETUP
====================================

  setup_scvi_env.bat      -> env "scvi"      -- Script 01
  setup_pyscenic_env.bat  -> env "pyscenic"  -- Script 03
  setup_oracle_env.bat    -> env "oracle"    -- Script 07


REFERENCE FILES (manual download required)
====================================

  DatabaseExtract_v_1.01.csv
      Lambert et al. 2018, "The Human Transcription Factors", Cell 172(4):650-665.
      http://humantfs.ccbr.utoronto.ca/download.php

  motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
      cisTarget motif-to-TF annotations (hg38, v10nr_clust).
      https://resources.aertslab.org/cistarget/motif2tf/

  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
      cisTarget ranking database, hg38, 10kb up/downstream of TSS.
      https://resources.aertslab.org/cistarget/databases/

Input data:
  adata_all_nk_after_mapping.h5ad
      Netskar et al. 2024, "Pan-cancer profiling of tumor-infiltrating natural
      killer cells through transcriptional reference mapping",
      Nature Immunology 25:1445-1459. doi:10.1038/s41590-024-01884-z
      https://www.nature.com/articles/s41590-024-01884-z

Files:
    checkpoint/  intermediate objects passed between scripts
    results/     final tables (CSV) -- numbered results01-20 consecutively
                 across the whole pipeline
    plots/       figures (PDF) -- numbered plot01-14 consecutively
    sanity/      threshold-titration and parameter-sweep tables (Script 01)
    session/     sessionInfo() / pip freeze logs per script
