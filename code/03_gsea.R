# -----------------------------------------------------------------------------------------------------------------------#
# Script: 03_gsea.R
# Descripción: RNA-seq downstream analysis template - Gene Set Enrichement Analysis
# Autor: Andrés Arroyo Barea
# Fecha: 2026-11-23
# -----------------------------------------------------------------------------------------------------------------------#

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(cli)
  #library(rna.seq.utils)
})

# ------------------------------------------------------------------------------#
# ------------- 1. Validate inputs from Stage 02 and configuration -------------#
# ------------------------------------------------------------------------------
required_objs <- c("cfg", "de_res_processed")
missing_objs <- setdiff(required_objs, ls())
if (length(missing_objs) > 0) {
  cli::cli_abort(
    "Missing required objects: {.val {missing_objs}}. Run earlier pipeline stages first."
  )
}

species <- cfg$gsea$species %||% "Homo sapiens"


# ------------------------------------------------------------------------------#
# -------------------2. Load and standardize MSigDB gene sets ------------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Loading MSigDB gene sets for {.val {species}}...")

msig_tbl <- msigdbr::msigdbr(species = species)

# Map each pathway set label to its MSigDB filter rule
pathway_filters <- list(
  HALLMARK     = \(df) dplyr::filter(df, gs_collection == "H"),
  GOBP         = \(df) dplyr::filter(df, gs_subcollection == "GO:BP"),
  GOMF         = \(df) dplyr::filter(df, gs_subcollection == "GO:MF"),
  GOCC         = \(df) dplyr::filter(df, gs_subcollection == "GO:CC"),
  REACTOME     = \(df) dplyr::filter(df, gs_subcollection == "CP:REACTOME"),
  KEGG         = \(df) dplyr::filter(df, gs_subcollection == "CP:KEGG_LEGACY"),
  KEGG_MEDICUS = \(df) dplyr::filter(df, gs_subcollection == "CP:KEGG_MEDICUS"),
  WP           = \(df) dplyr::filter(df, gs_subcollection == "CP:WIKIPATHWAYS")
)

pathways_list <- purrr::map(
  pathway_filters,
  \(filter_fn) {
    filter_fn(msig_tbl) %>%
      dplyr::transmute(
        term = as.character(gs_name),
        gene = as.character(db_ncbi_gene)
      ) %>%
      dplyr::filter(!is.na(term), nzchar(term), !is.na(gene), nzchar(gene)) %>%
      dplyr::distinct()
  }
)

pathway_summary <- purrr::imap_dfr(
  pathways_list,
  \(tbl, name) {
    tibble(
      pathway_set = name,
      n_terms = dplyr::n_distinct(tbl$term),
      n_genes = dplyr::n_distinct(tbl$gene)
    )
  }
)
print(pathway_summary)


# ------------------------------------------------------------------------------#
# ------------------------3. Build gene rankings per contrast ------------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Building ranked gene lists...")

rankings <- make_gsea_rankings(
  de_res_list = de_res_processed,
  OrgDb = org.Hs.eg.db,
  gene_type = cfg$gsea$gene_type %||% "ENSEMBL",
  gene_type_col = cfg$gsea$gene_type_col %||% "ENSEMBL_ID",
  rank_col = cfg$gsea$rank_col %||% "stat"
)


# ------------------------------------------------------------------------------#
# ------------------------4. Run GSEA across contrasts x pathway sets ----------#
# ------------------------------------------------------------------------------
cli::cli_inform(
  "Running GSEA across {length(rankings)} contrast(s) and {length(pathways_list)} pathway set(s)..."
)

gsea_grid <- tidyr::expand_grid(
  comparison = names(rankings),
  pathway_set = names(pathways_list)
) %>%
  dplyr::mutate(
    analysis_id = paste(comparison, pathway_set, sep = "__"),
    rank = purrr::map(comparison, \(comp) rankings[[comp]]),
    pathways = purrr::map(pathway_set, \(pset) pathways_list[[pset]]),
    result = purrr::pmap(
      list(rank, pathways),
      \(rk, pw) {
        run_gsea(
          rank = rk,
          pathways = pw,
          min_gs_size = cfg$gsea$min_gs_size %||% 15,
          max_gs_size = cfg$gsea$max_gs_size %||% 500,
          pvalue_cutoff = cfg$gsea$pvalue_cutoff %||% 1,
          p_adjust_method = cfg$gsea$p_adjust_method %||% "BH",
          seed = isTRUE(cfg$gsea$seed %||% TRUE)
        )
      }
    ),
    result_processed = purrr::pmap(
      list(result, comparison, pathway_set),
      \(res, comp, pset) {
        process_GSEA_results(
          gsea_result = res,
          comparison = comp,
          pathway_set = pset,
          OrgDb = org.Hs.eg.db,
          nes_threshold = cfg$gsea$nes_threshold %||% 1,
          padj_threshold = cfg$gsea$padj_threshold %||% 0.25,
          set_gene_symbol = isTRUE(cfg$gsea$set_gene_symbol %||% TRUE)
        )
      }
    )
  ) %>%
  dplyr::relocate(analysis_id, comparison, pathway_set)

# --- 5. Save outputs ----------------------------------------------------------
out_dir <- cfg$paths$intermediate_dir %||% "results/intermediate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_path <- file.path(out_dir, "03_gsea_grid.rds")
saveRDS(gsea_grid, file = rds_path)
cli::cli_alert_success("Saved GSEA grid to {.path {rds_path}}")