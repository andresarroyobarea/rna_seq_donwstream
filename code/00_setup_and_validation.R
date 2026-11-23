# -----------------------------------------------------------------------------------------------------------------------#
# Script: 00_setup_and_validation.R
# Descripción: RNA-seq downstream analysis template - Gene Set Enrichement Analysis
# Autor: Andrés Arroyo Barea
# Fecha: 2026-11-23
# -----------------------------------------------------------------------------------------------------------------------#

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(cli)
})

# -----------------------------------------------------------------------------#
# --------------- 1. Locate and load project YAML configuration ---------------#
# -----------------------------------------------------------------------------
config_path <- Sys.getenv(
  "RNASEQ_CONFIG_FILE",
  unset = "config/project_config_template.yml"
)

if (!file.exists(config_path)) {
  fallback_path <- "config/project_config_template.yml"
  if (file.exists(fallback_path)) {
    cli::cli_warn(
      "Config file {.path {config_path}} not found; falling back to {.path {fallback_path}}."
    )
    config_path <- fallback_path
  } else {
    cli::cli_abort("No configuration YAML found at {.path {config_path}}.")
  }
}

cli::cli_inform("Loading configuration from {.path {config_path}}...")
cfg <- yaml::read_yaml(config_path)

required_sections <- c("paths", "design", "contrasts", "de")
missing_sections <- setdiff(required_sections, names(cfg))
if (length(missing_sections) > 0) {
  cli::cli_abort(
    "Configuration is missing required section(s): {.val {missing_sections}}."
  )
}


# -----------------------------------------------------------------------------#
# --------------------- 2. Validate input file paths --------------------------#
# -----------------------------------------------------------------------------
input_keys <- c("counts_file", "metadata_file", "gene_annotation_file")
for (k in input_keys) {
  fpath <- cfg$paths[[k]]
  if (is.null(fpath) || !nzchar(fpath) || !file.exists(fpath)) {
    cli::cli_abort(
      "Missing or invalid path for {.field paths.{k}}: {.path {fpath %||% '<empty>'}}."
    )
  }
}

# -----------------------------------------------------------------------------#
# ------------- 3. Load metadata, counts, gene annotation ---------------------#
# -----------------------------------------------------------------------------
cli::cli_inform("Reading sample metadata, count matrix, and gene annotation...")

meta_df <- utils::read.delim(
  cfg$paths$metadata_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts_raw <- utils::read.delim(
  cfg$paths$counts_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_info <- utils::read.delim(
  cfg$paths$gene_annotation_file,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rownames(meta_df) <- meta_df[[sample_id_col]]

count_mat <- as.matrix(counts_raw[, setdiff(colnames(counts_raw), gene_id_col), drop = FALSE])
rownames(count_mat) <- counts_raw[[gene_id_col]]
storage.mode(count_mat) <- "integer"


# -----------------------------------------------------------------------------#
# ---------- 4. Apply explicit sample exclusions if declared in YAML ----------#
# -----------------------------------------------------------------------------
excluded <- cfg$design$excluded_samples %||% character(0)
if (length(excluded) > 0) {
  cli::cli_inform("Excluding {length(excluded)} sample(s) declared in config: {.val {excluded}}")
  meta_df <- meta_df[!rownames(meta_df) %in% excluded, , drop = FALSE]
  count_mat <- count_mat[, !colnames(count_mat) %in% excluded, drop = FALSE]
}

if (nrow(meta_df) < 2 || ncol(count_mat) < 2) {
  cli::cli_abort("Too few samples remaining after exclusions ({nrow(meta_df)}).")
}


# -----------------------------------------------------------------------------#
# ----- 5. Validate sample alignment, count integrity, and design factors -----#
# -----------------------------------------------------------------------------
common_samples <- intersect(rownames(meta_df), colnames(count_mat))
if (length(common_samples) != nrow(meta_df) || length(common_samples) != ncol(count_mat)) {
  cli::cli_abort(
    "Sample mismatch between metadata ({nrow(meta_df)}) and count matrix ({ncol(count_mat)})."
  )
}

# Align count matrix columns to metadata row order
count_mat <- count_mat[, rownames(meta_df), drop = FALSE]

if (anyNA(count_mat) || any(count_mat < 0L)) {
  cli::cli_abort("Count matrix contains NA or negative values.")
}

cond_col <- cfg$design$condition_col
if (is.null(cond_col) || !cond_col %in% colnames(meta_df)) {
  cli::cli_abort("Design condition column {.val {cond_col}} not found in metadata.")
}

meta_df[[cond_col]] <- factor(meta_df[[cond_col]])

ref_level <- cfg$design$reference_level
if (!is.null(ref_level) && nzchar(ref_level)) {
  if (!ref_level %in% levels(meta_df[[cond_col]])) {
    cli::cli_abort(
      "Reference level {.val {ref_level}} not found in {.field {cond_col}} ({.val {levels(meta_df[[cond_col]])}})."
    )
  }
  meta_df[[cond_col]] <- stats::relevel(meta_df[[cond_col]], ref = ref_level)
}

# Validate contrasts against metadata columns and factor levels
for (nm in names(cfg$contrasts)) {
  ctr <- cfg$contrasts[[nm]]
  if (length(ctr) != 3) {
    cli::cli_abort("Contrast {.val {nm}} must have 3 elements: c(column, numerator, denominator).")
  }
  col_name <- ctr[[1]]
  num_lvl  <- ctr[[2]]
  den_lvl  <- ctr[[3]]
  
  if (!col_name %in% colnames(meta_df)) {
    cli::cli_abort("Contrast {.val {nm}} references missing metadata column {.val {col_name}}.")
  }
  valid_lvls <- unique(as.character(meta_df[[col_name]]))
  missing_lvls <- setdiff(c(num_lvl, den_lvl), valid_lvls)
  if (length(missing_lvls) > 0) {
    cli::cli_abort(
      "Contrast {.val {nm}} references level(s) not present in {.field {col_name}}: {.val {missing_lvls}}."
    )
  }
}


# -----------------------------------------------------------------------------#
# ----------- 6. Build initial DESeqDataSet and persist setup bundle ----------#
# -----------------------------------------------------------------------------
design_formula <- stats::as.formula(cfg$design$formula %||% paste("~", cond_col))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_df,
  design = design_formula
)

inter_dir <- cfg$paths$intermediate_dir %||% "results/intermediate"
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)

setup_bundle <- list(
  cfg = cfg,
  dds = dds,
  gene_info = gene_info,
  n_samples = ncol(dds),
  n_genes = nrow(dds),
  excluded_samples = excluded,
  created_at = Sys.time()
)

bundle_path <- file.path(inter_dir, "00_setup_bundle.rds")
saveRDS(setup_bundle, file = bundle_path)

cli::cli_alert_success(
  "Setup validated ({ncol(dds)} samples x {nrow(dds)} genes). Saved {.path {bundle_path}}."
)
