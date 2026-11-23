# -----------------------------------------------------------------------------------------------------------------------#
# Script: 05_visualizations.R
# Descripción: RNA-seq downstream analysis template - Gene Set Enrichement Analysis
# Autor: Andrés Arroyo Barea
# Fecha: 2026-11-23
# -----------------------------------------------------------------------------------------------------------------------#

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(cli)
  #library(rna.seq.utils)
})

# ------------------------------------------------------------------------------#
# -------- 1. Validate inputs from previous stages and configuration -----------#
# ------------------------------------------------------------------------------
required_objs <- c("cfg", "de_res_processed")
missing_objs <- setdiff(required_objs, ls())
if (length(missing_objs) > 0) {
  cli::cli_abort(
    "Missing required objects: {.val {missing_objs}}. Run earlier stages first."
  )
}

has_dds <- exists("dds")
has_gsea <- exists("gsea_grid")

fig_base <- cfg$paths$figures_dir %||% "results/figures"
fig_de_dir <- file.path(fig_base, "de")
fig_expr_dir <- file.path(fig_base, "expression")
fig_gsea_dir <- file.path(fig_base, "gsea")

purrr::walk(
  c(fig_de_dir, fig_expr_dir, fig_gsea_dir),
  \(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)
)

# ------------------------------------------------------------------------------#
# -------------- 2. Differential Expression plots per contrast -----------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Generating DE plots across {length(de_res_processed)} contrast(s)...")

ma_plots <- purrr::imap(
  de_res_processed,
  \(df, contrast_name) {
    make_MA_plot(
      data = df,
      title = contrast_name,
      y_limit = cfg$plots$ma$y_limit %||% 6
    )
  }
)

lfc_plots <- purrr::imap(
  de_res_processed,
  \(df, contrast_name) {
    make_LFC_plot(
      data = df,
      title = contrast_name
    )
  }
)

pval_plots <- purrr::imap(
  de_res_processed,
  \(df, contrast_name) {
    make_PVAL_plot(
      data = df,
      title = contrast_name
    )
  }
)

volcano_plots <- purrr::imap(
  de_res_processed,
  \(df, contrast_name) {
    make_VOLCANO_plot(
      data = df,
      title = contrast_name,
      x_limit = cfg$plots$volcano$x_limit %||% 5,
      y_limit = cfg$plots$volcano$y_limit %||% 15,
      lfc_threshold = cfg$de$lfc_threshold %||% 0.58,
      padj_threshold = cfg$de$alpha %||% 0.05,
      label_n = cfg$plots$volcano$label_n %||% 10
    )
  }
)

# Export DE plots per contrast
purrr::iwalk(ma_plots, \(p, nm) {
  ggplot2::ggsave(file.path(fig_de_dir, paste0("ma_", nm, ".png")), plot = p, width = 7, height = 5, dpi = 300)
})
purrr::iwalk(lfc_plots, \(p, nm) {
  ggplot2::ggsave(file.path(fig_de_dir, paste0("lfc_", nm, ".png")), plot = p, width = 7, height = 5, dpi = 300)
})
purrr::iwalk(pval_plots, \(p, nm) {
  ggplot2::ggsave(file.path(fig_de_dir, paste0("pval_", nm, ".png")), plot = p, width = 7, height = 5, dpi = 300)
})
purrr::iwalk(volcano_plots, \(p, nm) {
  ggplot2::ggsave(file.path(fig_de_dir, paste0("volcano_", nm, ".png")), plot = p, width = 7, height = 6, dpi = 300)
})


# ------------------------------------------------------------------------------#
# ------------- 3. Expression panels across groups -----------------------------#
# ------------------------------------------------------------------------------
if (has_dds) {
  cli::cli_inform("Generating gene expression panels...")
  
  palette_cfg <- cfg$plots$expression$palette
  palette_vec <- if (!is.null(palette_cfg) && length(palette_cfg) > 0) {
    unlist(palette_cfg)
  } else {
    NULL
  }
  
  expression_plots <- make_all_EXPRESSION_panels(
    dds = dds,
    de_results = de_res_processed,
    comparisons = cfg$contrasts,
    n = cfg$plots$expression$n_genes %||% 12,
    panel_size = cfg$plots$expression$panel_size %||% 6,
    show_boxplot = isTRUE(cfg$plots$expression$show_boxplot %||% FALSE),
    sample_metadata_cols = cfg$plots$expression$metadata_cols %||% (cfg$design$condition_col %||% "genotype"),
    colors = palette_vec
  )
} else {
  cli::cli_inform("Skipping expression panels ({.val dds} not available).")
}

# ------------------------------------------------------------------------------#
# -------------------- 4. GSEA top pathways & barplots -------------------------#
# ------------------------------------------------------------------------------
if (has_gsea) {
  cli::cli_inform("Selecting top GSEA pathways and building barplots...")
  
  top_n_pathways <- cfg$plots$gsea$top_n %||% 5
  
  gsea_plots <- gsea_grid %>%
    dplyr::mutate(
      top_pathways = purrr::map(
        result_processed,
        \(res_tbl) {
          if (is.null(res_tbl) || nrow(res_tbl) == 0) return(NULL)
          select_top_pathways(res_tbl, n = top_n_pathways)
        }
      ),
      barplot = purrr::pmap(
        list(top_pathways, comparison, pathway_set),
        \(top_tbl, comp, pset) {
          if (is.null(top_tbl) || nrow(top_tbl) == 0) return(NULL)
          barplot_GSEA(top_tbl, comp, pset)
        }
      )
    )
} else {
  cli::cli_inform("Skipping GSEA plots ({.val gsea_grid} not available).")
}

cli::cli_alert_success("Stage 04 visualizations completed.")