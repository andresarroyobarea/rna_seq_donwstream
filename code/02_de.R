# -----------------------------------------------------------------------------------------------------------------------#
# Script: 02_de.R
# Descripción: RNA-seq downstream analysis template - Differential expression analysis ((DESeq2 + shrinkage + formatting)
# Autor: Andrés Arroyo Barea
# Fecha: 2026-11-23
# -----------------------------------------------------------------------------------------------------------------------#
suppressPackageStartupMessages({
  library(DESeq2)
  library(purrr)
  library(dplyr)
  library(cli)
  #library(rna.seq.utils) # Paquete rna-seq-utils.
})

# --------------------------------------------------------#
# --- 1. Validate objects and params from 00_setup.R -----#
# --------------------------------------------------------
required_objs <- c("cfg", "dds", "gene_info")
missing_objs <- setdiff(required_objs, ls())
if (length(missing_objs) > 0) {
  cli::cli_abort(
    "Missing required objects in environment: {.val {missing_objs}}. Run 00_setup first."
  )
}

if (!inherits(dds, "DESeqDataSet")) {
  cli::cli_abort("{.arg dds} must be a valid {.cls DESeqDataSet} object.")
}

contrasts_list <- cfg$contrasts
if (length(contrasts_list) == 0) {
  cli::cli_abort("No contrasts defined under {.field contrasts} in project config.")
}

contrasts_list <- purrr::imap(contrasts_list, \(vec, nm) {
  if (is.null(nm) || !nzchar(nm)) {
    cli::cli_abort("All contrast entries must have a non-empty name.")
  }
  if (!is.character(vec) || length(vec) != 3) {
    cli::cli_abort(
      "Contrast {.val {nm}} must be a character vector of length 3: c(factor, numerator, denominator)."
    )
  }
  vec
})

# ------------------------------------------------------------------------------#
# ----------------------- 2. Fitting DESeq2 model ------------------------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Fitting DESeq2 model...")

dds <- DESeq2::DESeq(
  object = dds,
  betaPrior = isTRUE(cfg$de$beta_prior),
  parallel = isTRUE(cfg$de$parallel)
)


# ------------------------------------------------------------------------------#
# -------------------- 3. Extract Deseq2 results by contrast -------------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Extracting DESeq2 results for {length(contrasts_list)} contrast(s)...")

de_res_raw <- purrr::imap(
  contrasts_list,
  \(contrast_vec, contrast_name) {
    cli::cli_bullets(c("*" = "Contrast: {.val {contrast_name}}"))
    DESeq2::results(
      object = dds,
      contrast = contrast_vec,
      alpha = cfg$de$alpha,
      lfcThreshold = cfg$de$lfc_threshold,
      altHypothesis = cfg$de$alt_hypothesis,
      pAdjustMethod = cfg$de$p_adjust_method,
      independentFiltering = cfg$de$independent_filtering,
      parallel = isTRUE(cfg$de$parallel)
    )
  }
)

# ------------------------------------------------------------------------------#
# ------------------------- 4. Apply LFC shrinkage -----------------------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Applying LFC shrinkage ({.val {cfg$de$shrink_type}})...")

de_res_shrunk <- purrr::imap(
  de_res_raw,
  \(res_obj, contrast_name) {
    DESeq2::lfcShrink(
      dds = dds,
      contrast = contrasts_list[[contrast_name]],
      res = res_obj,
      type = cfg$de$shrink_type,
      quiet = TRUE
    )
  }
)

# Check list before combine.
stopifnot(
  length(de_res_raw) == length(de_res_shrunk),
  identical(names(de_res_raw), names(de_res_shrunk))
)

# ------------------------------------------------------------------------------#
# --- 4. Processing and standardize tablas mediante `rna.seq.utils` ------------#
# ------------------------------------------------------------------------------
cli::cli_inform("Processing DE tables with {.pkg rna.seq.utils}...")

de_res_processed <- purrr::imap(
  de_res_raw,
  \(res_obj, contrast_name) {
    df <- process_DE_results(
      res = res_obj,
      res_shrink = de_res_shrunk[[contrast_name]],
      gene_metadata = gene_info
    )
    dplyr::mutate(df, contrast = contrast_name, .before = 1)
  }
)

# ------------------------------------------------------------------------------#
# ----------------- 5. Save intermediate and final results ---------------------#
# ------------------------------------------------------------------------------
dir.create(cfg$paths$intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$paths$tables_dir, recursive = TRUE, showWarnings = FALSE)

de_bundle <- list(
  raw = de_res_raw,
  shrunk = de_res_shrunk,
  processed = de_res_processed,
  parameters = cfg$de,
  contrasts = contrasts_list,
  timestamp = Sys.time()
)

rds_out <- file.path(cfg$paths$intermediate_dir, "02_de_results.rds")
saveRDS(de_bundle, file = rds_out)
cli::cli_alert_success("Saved DE bundle to {.path {rds_out}}")

# Export standardized table by por contrast in TSV format.
purrr::iwalk(
  de_res_processed,
  \(tbl, contrast_name) {
    out_file <- file.path(cfg$paths$tables_dir, paste0("de_", contrast_name, ".tsv"))
    utils::write.table(tbl, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
)
cli::cli_alert_success("Exported {length(de_res_processed)} DE table(s) to {.path {cfg$paths$tables_dir}}")