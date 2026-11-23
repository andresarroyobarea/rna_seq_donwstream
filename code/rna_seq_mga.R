# --------------------------------------------------------------------------------------------------------#
# Script: rna_seq_mga.R
# Descripción: 3'-RNAseq analysis for PGA1 - TP53 AND MGA KO.
# Autor: Andrés Arroyo Barea
# Fecha: 2026-07-10
# --------------------------------------------------------------------------------------------------------#

# ---- Libraries ---- #
pacman::p_load(DESeq2, edgeR, limma, dplyr, tidyverse, openxlsx, 
               reshape2, variancePartition, PCAtools, glmpca, magrittr, RColorBrewer, pheatmap, patchwork,
               PoiClaClu, ggbeeswarm, EnhancedVolcano, ggpubr,
               clusterProfiler, msigdbr, org.Hs.eg.db, enrichplot, ReactomePA,
               purrr, here)

# Load functions

files <- list.files(
  here("code/shared"),
  pattern = "\\.R$",
  full.names = TRUE
)
purrr::walk(
  files,
  source
)

#---------------------------------------------------------------------------------------------#
# -------------------------- 1. DATA PROCESSING ----------------------------------------------#
#---------------------------------------------------------------------------------------------#
data <- read.csv("./results/upstream/counts.tsv", sep = "\t", header = T, skip = 1, check.names = F, row.names = 1)

# Improve names in colnames
colnames(data) <- gsub("_Aligned.sortedByCoord.out.bam$", "", gsub("^.*/", "", colnames(data)))
colnames(data) <- gsub("-", "_", colnames(data))

# Extract ENSEMBL-Gene id correspondencies
gene_info <- data %>% dplyr::select(Chr:gene_name)

# Subset samples.
data <- data %>% dplyr::select(PGA1_Double_KO_S12:PGA1_WT_S4)

#---------------------------------------------------------------------------------------------#
# -------------------------- 2. METADATA PROCESSING ------------------------------------------#
#---------------------------------------------------------------------------------------------#
metadata <- read.xlsx("./metadata/rna_seq_mga_metadata_2026_07_10.xlsx", detectDates = T)
metadata <- metadata %>% column_to_rownames(., 'sample_id')

metadata$genotype <- factor(metadata$genotype, levels = c("WT", "TP53_KO", "MGA_KO", "Double_KO"))
metadata$harvest_batch <- factor(metadata$harvest_batch, levels = c("H1", "H2", "H3"))


glimpse(metadata)

# Turn variables into factors
#metadata <- metadata %>%
#  dplyr::mutate(across(c(cell_status, treatment, replicate, replicate_global, processing_batch), as.factor)) %>%
#  dplyr::mutate(
#    dna_ext_batch = as.factor(case_when(
#      dna_extraction_date == "2026-11-06" ~ 1,
#      dna_extraction_date == "2026-02-20" ~ 2,
#      TRUE ~ NA
#    )), .after = dna_extraction_date
#  ) %>%
#  dplyr::mutate(
#    group = interaction(cell_status, replicate, sep = "_")
#  )


# Reorder data columns according to metadata columns.
data <- data %>% dplyr::select(rownames(metadata))

# Check order
all(colnames(data) == rownames(metadata))

#-------------------------------------------------------------------------------------------------#
# -------------------------- 3. REMOVE DUPICLATED GENES ------------------------------------------#
#-------------------------------------------------------------------------------------------------#

# Seleccionamos todos los ENSEMBL ID QUE MAPEAN MÁS DE UNA VEZ AL MISMO GEN.
ens_dup_gene_sym <- data %>% 
  merge(., gene_info, by = 0) %>%
  rename(ensembl_id = "Row.names") %>%
  dplyr::select(-(Chr:Length)) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  group_by(gene_name) %>%
  arrange(gene_name) %>%
  dplyr::filter(n() > 1) %>%
  pull(ensembl_id)

# Seleccionamos dentro de los ENSEMBL anteriores (482 gene IDs), cuales dejamos en nuestros datos.
# CRITERIO DE SELECCION: Aquellos ENSEMBL ID con mayor mediana de expresión. En caso de empate,
# EL MAYOR ENSEMBL ID.
ens_dup_gene_sym_to_select <- data %>% 
  merge(., gene_info, by = 0) %>%
  rename(ensembl_id = "Row.names") %>%
  dplyr::select(-(Chr:Length)) %>%
  dplyr::filter(!is.na(gene_name)) %>%
  group_by(gene_name) %>%
  arrange(gene_name) %>%
  dplyr::filter(n() > 1) %>%
  ungroup() %>%  
  dplyr::mutate(
    gene_median = rowMedians(
      as.matrix(
        dplyr::select(., PGA1_WT_S1:PGA1_Double_KO_S14) %>% 
          mutate(across(everything(), as.numeric))
      ), 
      na.rm = TRUE
    )
  ) %>% 
  group_by(gene_name) %>%
  # Ordena primero por mediana descendente y, en caso de empate, por ENSEMBL ID ascendente
  arrange(desc(gene_median), ensembl_id) %>%  
  slice(1) %>%  
  ungroup() %>% 
  pull(ensembl_id)


# To select ENSEMBL ID to remove.
ens_dup_gene_sym_to_remove <- setdiff(ens_dup_gene_sym, ens_dup_gene_sym_to_select)

# Remove duplicated ENSEMBL ID.
data <- data %>% filter(!rownames(.) %in% ens_dup_gene_sym_to_remove)

# Los eliminamos de gene_info
gene_info <- gene_info %>% filter(!rownames(.) %in% ens_dup_gene_sym_to_remove)

rm(ens_dup_gene_sym, ens_dup_gene_sym_to_remove, ens_dup_gene_sym_to_select)


#-----------------------------------------------------------------------------------------------#
# -------------------------- 4.1 LIBRARY SIZE AND GENE DETECTION ---------------------------------#
#-------------------------------------------------------------------------------------------------#

# Los clones se seleccionar en base a que estuviesen en homocigosidad. que fuese delección no multiple de 3
# para evitar codon de stop

# Para TP53 y la partte KO De TP53 en el doble se ve que NO había proteína en WG -> OK
# MGA comprobado por sanger y por secuenciacion pero no por proteína.

# MGA intentar ver si a nivel de transcrito hasta que punto (Exon se transcribe)

# Pellet y extracción de RNA

colSums(data) %>%
  enframe(name = "sample", value = "library_size") %>%
  tidyr::extract(
    sample,
    into = c("cell", "genotype", "rep"),
    regex = "^(PGA1)_(.+)_(S\\d+)$"
  ) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "TP53_KO", "MGA_KO", "Double_KO")),
    sample = factor(
      paste(genotype, rep, sep = "_"),
      levels = c(
        "WT_S1", "WT_S2", "WT_S4",
        "TP53_KO_S3", "TP53_KO_S9", "TP53_KO_S10",
        "MGA_KO_S3", "MGA_KO_S19", "MGA_KO_S20",
        "Double_KO_S5", "Double_KO_S12", "Double_KO_S14"
      )
  )) %>%
  ggplot(aes(sample, library_size, fill = genotype)) +
  geom_col(width = 0.8) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "black")
  ) +
  labs(
    x = NULL,
    y = "Library size",
    fill = "Genotype"
  ) +
  scale_fill_manual(values = c(
    WT = "#4C72B0",
    TP53_KO = "#DD8452",
    MGA_KO = "green4",
    Double_KO = "gold3"
  ))

colSums(data)
colSums(data > 0) %>%
  enframe(name = "sample", value = "genes_detected") %>%
  tidyr::extract(
    sample,
    into = c("cell", "genotype", "rep"),
    regex = "^(PGA1)_(.+)_(S\\d+)$"
  ) %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "TP53_KO", "MGA_KO", "Double_KO")),
    sample = factor(
      paste(genotype, rep, sep = "_"),
      levels = c(
        "WT_S1", "WT_S2", "WT_S4",
        "TP53_KO_S3", "TP53_KO_S9", "TP53_KO_S10",
        "MGA_KO_S3", "MGA_KO_S19", "MGA_KO_S20",
        "Double_KO_S5", "Double_KO_S12", "Double_KO_S14"
      )
    )) %>%
  ggplot(aes(sample, genes_detected, fill = genotype)) +
  geom_col(width = 0.8) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "black")
  ) +
  labs(
    x = NULL,
    y = "Genes detected (N reads > 10)",
    fill = "Genotype"
  ) +
  #scale_y_continuous(breaks = seq(0, 15000, by = 3000)) +
  scale_fill_manual(values = c(
    WT = "#4C72B0",
    TP53_KO = "#DD8452",
    MGA_KO = "green4",
    Double_KO = "gold3"
  ))


#---------------------------------------------------------------------------------------------#
# -------------------------- 5. GENE PREFILTERING --------------------------------------------#
#---------------------------------------------------------------------------------------------#
# filterByExpr function was used to select only genes that have sufficiently large counts to be retained in a statistical analysis.
# Genes before a threshold of CPM counts in a k numer of samples were removed. This work is performed by the own function.
# Model matrix

# Sensitivity analysis - NO MGA19
# filterByExpr() retuvo 13.093 genes en el análisis completo y 12.955 tras excluir MGA19, con 12.955 genes compartidos (98,95%) 
# y 138 (1%) genes que únicamente pasan el filtrado en presencia de MGA19.

# Dada la pequeña proporcion de genes que se pierden del pre-filtrado si no se tiene en cuenta la muestra MGA19, se decidió
# trabajar con los 13093 genes del dataset completo en el análisis de sensibilidad para poder hacer una comparación más justa 
# entre todos los genes.

# Global
global_mat <- model.matrix(~0 + genotype, data = metadata)

# Apply filterbyExpress
keep_global <- filterByExpr(data, global_mat)

# 13093 (21.2%) out of 61617 genes were retained
sum(keep_global)

# Filter count table and gene id table according to remaining genes.

# Global 
data <- data[keep_global, ]
gene_info_global <- gene_info %>% dplyr::filter(rownames(.) %in% rownames(data))


rm(global_mat)


#---------------------------------------------------------------------------------------------#
# -------------------------- 4. CREATE DESEQ OBJECT ------------------------------------------#
#---------------------------------------------------------------------------------------------#

# GLOBAL
dds_global <- DESeqDataSetFromMatrix(
  countData = data,
  colData = metadata,
  design = ~ genotype,
  rowData = gene_info_global
  )

# dds sin MGA-S19
dds_global_no_S19 <- dds_global[, colnames(dds_global) != "PGA1_MGA_KO_S19"]

#---------------------------------------------------------------------------------------------#
# -------------------------- 5. EXPLORATORY DATA ANALYSIS ------------------------------------#
#---------------------------------------------------------------------------------------------#

# Vst normalization for each subset.
vst_global <- vst(dds_global, blind = F)

# -------------------------- 5.3.2 TP53 AND MGA ---------------------------####
# TP53: ENSG00000141510
# MGA: ENSG00000174197

dds_global <- DESeq(dds_global, betaPrior = F)

genes <- c("ENSG00000141510", "ENSG00000174197")
           
# TP53 targets
tp53_targets <- c("ENSG00000124762", "ENSG00000135679", "ENSG00000105327")

# Loss of MGA repression mediated by an atypical polycomb complex promotes tumor progression and invasiveness
mga_targets_paper1 <- c("ENSG00000066923", "ENSG00000078487")

# MGA, L3MBTL2 and E2F6 determine genomic binding of the non-canonical Polycomb repressive complex PRC1.6
mga_targets_paper2 <-  c("ENSG00000178403", "ENSG00000161860", "ENSG00000151164")


df <- bind_rows(
  lapply(tp53_targets, function(g) {
    plotCounts(
      dds_global,
      gene = g,
      intgroup = "genotype",
      returnData = TRUE
    ) %>%
      rownames_to_column("Sample") %>%
      mutate(emsembl_id = g)
  })
) %>% dplyr::mutate(clone = str_extract(Sample, "S\\d+$")) %>%
  left_join(
    gene_info_global %>%
      rownames_to_column("emsembl_id") %>%
      dplyr::select(emsembl_id, gene_name),
    by = "emsembl_id"
  )

ggplot(df, aes(genotype, count, fill = genotype, col = genotype)) +
  geom_point(
    size = 5,
    position = position_jitter(width = 0.08)
  ) +
  geom_text_repel(
    aes(label = clone),
    point.padding = 0.3,   # separación respecto al punto
    box.padding = 0.4,     # separación entre etiquetas
    max.overlaps = Inf
  ) +
  facet_wrap(~gene_name, scales = "free_y") +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Normalized counts"
  )  + theme(
    legend.position = "none",
    title = element_text(size = 18),
    axis.text = element_text(size = 9 ),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 18, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 18, colour = "black", face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, colour = "black", face = "plain", margin = margin(t = 10)),
    strip.text = element_text(size = 18)
  ) +
  scale_color_manual(values = c(
    WT = "#4C72B0",
    TP53_KO = "#DD8452",
    MGA_KO = "green4",
    Double_KO = "gold3"
  ))


# -------------------------- 5.4 PCA------------------------------------------------####

#--------------------- 5.1 PCA GLOBAL ------------------------- #

# Select the top 2000 variables genes for PCA.
vars <- apply(assay(vst_global), 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:1000]

# PCA tools
pca_analysis <- pca(assay(vst_global) %>% as.data.frame() %>% filter(rownames(.) %in% top_genes),
                    metadata = metadata, removeVar = 0.0, center = T, scale = F)


# Screeplot - As expected, large cell lines differences.
screeplot(pca_analysis, axisLabSize = 18, titleLabSize = 22, title = "Screeplot PCA")

# Full data was used here
pca_global <- biplot(pca_analysis,
       x = "PC1", y = "PC3",
       colby = 'genotype',
       shape = "harvest_batch",
       colkey = c(
         "WT" = "#4C72B0",
         "TP53_KO" = "#DD8452",
         "MGA_KO" = "green4",
         "Double_KO" = "gold3"
       ),
       colLegendTitle = '',
       #lab = NULL,
      lab = pca_analysis$metadata$clone,
      labSize = 5,
       pointSize = 5,
       hline = 0, vline = 0, 
       legendPosition = 'right', 
       legendLabSize = 14, 
       legendIconSize = 5.0,
       #title = "PCA all samples - Filtered genes",
       #encircle = T,
       encircleFill = F,
       showLoadings = F,
       encircleLineSize = 4, 
       encircleAlpha = 1/5,
       encircleLineCol = "black",
       sizeLoadingsNames = 2)


pca_global + theme_bw() +
  theme(
    legend.position = "right",
    title = element_text(size = 18),
    legend.text = element_text(size = 15),
    axis.text = element_text(size = 18 ),
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 18, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 18, colour = "black", face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, colour = "black", face = "plain", margin = margin(t = 10)),
    strip.text = element_text(size = 18)
  )


# Additional check
# ---- GLM PCA ----
gpca <- glmpca(counts(dds_global), L=2)
gpca.dat <- gpca$factors
gpca.dat$genotype <- dds_global$genotype
gpca.dat$clone <- dds_global$clone
gpca.dat$harvest_batch <- dds_global$harvest_batch

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = genotype, shape = harvest_batch)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

colData(dds_global)



#  Pairs plot
pairsplot(pca_analysis,
          colby = 'genotype',
          shape = "harvest_batch",
          #colkey = c(
          #   "CART" =  "#B39DDB",
          #   "MM" = "#FF8000"
          #),
          #lab = NULL,
          lab = pca_analysis$metadata$clone,
          labSize = 5,
          pointSize = 5,
          hline = 0, vline = 0, 
          legendPosition = 'right', 
          legendLabSize = 14, 
          legendIconSize = 5.0,
          #title = "PCA all samples - Filtered genes",
          #encircle = T,
          encircleFill = F,
          showLoadings = F,
          encircleLineSize = 4, 
          encircleAlpha = 1/5,
          encircleLineCol = "black",
          sizeLoadingsNames = 2)


#  ---- VariancePartition ---
register(
  MulticoreParam(4)
)

# CART
form <- ~ harvest_batch + genotype

varPart_CART <- fitExtractVarPartModel(
  as.data.frame(assay(vst_global)),
  form,
  as.data.frame(colData(vst_global)),
  BPPARAM = MulticoreParam(4)
)

vp_CART <- sortCols(varPart_CART)
plotVarPart(vp_CART)
colData(vst_mm)

plotPercentBars(
  vp_CART[order(vp_CART$treatment,
                decreasing = TRUE)[1:10],]
)

# MM
varPart_MM <- fitExtractVarPartModel(
  as.data.frame(assay(vst_mm)),
  form,
  as.data.frame(colData(vst_mm)),
  BPPARAM = MulticoreParam(4)
)

vp_MM <- sortCols(varPart_MM)
plotVarPart(vp_MM)

plotPercentBars(
  vp_MM[order(vp_MM$treatment,
                decreasing = TRUE)[1:10],]
)

PC_scores <- plotPCA(
  vst_global,
  intgroup=c("harvest_batch","genotype"),
  returnData=F,
  ntop = 2000
)
PC_scores %>%
  group_by(replicate) %>%
  summarise(
    delta_PC1 = PC1[treatment=="UroA"] - PC1[treatment=="DMSO"],
    delta_PC2 = PC2[treatment=="UroA"] - PC2[treatment=="DMSO"]
  )

#---------------------------- 6.CLUSTERING ANALYSIS -------------------------------------####

# ------------ Color format --------------.

# Col dataframe
ann_col_df_global <- data.frame(genotype = metadata$genotype, harvest_batch = metadata$harvest_batch) %>% set_rownames(., value = rownames(metadata))

ann_col_df_global <- ann_col_df_global %>% dplyr::slice(-c(2, 8))

# Create a named vector of colors for the gruops.
grupo_colors_global <- list(genotype = c("WT" = "#4C72B0",
                            "TP53_KO" = "#DD8452",
                            "MGA_KO" = "green4",
                            "Double_KO"= "gold3"),
                            harvest_batch = c("H1" =  "skyblue",
                                              "H2" = "salmon",
                                              "H3" = "grey70"))

# Create a list of annotation row colors
ann_row_col_global <- list(genotype = grupo_colors_global)

# Color palette
heat.colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# ----------- Select Top 2000 most variable genes ---

# Calcular varianza de cada gen
gene_vars_global <- apply(assay(vst_global), 1, var)
# Seleccionar los genes más variables
top_var_genes_global <- head(order(gene_vars_global, decreasing = TRUE), 2000)


# ----- Euclidean distances on top 2000 -------- #
# Global
vst_global <- vst_global[, !colnames(vst_global) %in% c("PGA1_MGA_KO_S19", "PGA1_WT_S2")]
vst_global_top2000 <- assay(vst_global)[top_var_genes_global, ]
global_sample_dist <-dist(t(vst_global_top2000))
global_sample_dist_matrix <- as.matrix(global_sample_dist)


# ------------- heatmap --------------- #
rownames(global_sample_dist_matrix) <- sub("^PGA1_", "", rownames(global_sample_dist_matrix))
colnames(global_sample_dist_matrix) <- sub("^PGA1_", "", colnames(global_sample_dist_matrix))
rownames(ann_col_df_global) <- sub("^PGA1_", "", rownames(ann_col_df_global))

# GLOBAL
pheatmap(global_sample_dist_matrix,
         clustering_distance_cols = global_sample_dist,
         clustering_distance_rows = global_sample_dist,
         annotation_col = ann_col_df_global,
         annotation_colors = grupo_colors_global,
         color = heat.colors,
         fontsize = 15)


colnames(data)

rm(top_var_genes_global, top_var_genes_mec1, top_var_genes_pga,
   grupo_colors_global, grupo_colors_mec1, grupo_colors_pga, gene_vars_global,
   gene_vars_mec1, gene_vars_pga, cor_mat_global, cor_mat_mec1, cor_mat_pga)



#### 6. Differential gene expression ####
dds_global <- DESeq(dds_global, betaPrior = F)

# Sensitivity analysis
dds_global_no_S19 <- DESeq(dds_global_no_S19, betaPrior = F)

# ---- Dispersion plot ----- #
plotDispEsts(dds_global)
plotDispEsts(dds_global_no_S19)


# ---- Normalization assesment - Raw counts vst normalized ---
raw_counts_melt <- counts(dds_global, normalized = F) %>% 
  data.frame() %>% melt() %>%
  rename(sample = variable, expression = value) %>% 
  mutate(expression = log2(expression+1))

norm_counts_melt <- counts(dds_global, normalized = TRUE) %>% 
  data.frame() %>% melt() %>%
  rename(sample = variable, expression = value) %>% 
  mutate(expression = log2(expression+1))

ggplot(raw_counts_melt, aes(x = sample, y = expression)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Non-normalized read counts distribution") +
  ylab("Raw gene expression (log2)") +
  xlab("") +
  theme(
    legend.position = "none",
    title = element_text(size = 18),
    axis.text = element_text(size = 9 ),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.text.x = element_text(size = 15, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 17, colour = "black", face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = 17, colour = "black", face = "plain", margin = margin(t = 10)),
    strip.text = element_text(size = 15)
  ) 

ggplot(norm_counts_melt, aes(x = sample, y = expression)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Normalized read counts distribution") +
  ylab("Normalized gene expression (log2)") +
  xlab("") +
  theme(
    legend.position = "none",
    title = element_text(size = 18),
    axis.text = element_text(size = 9 ),
    axis.text.y = element_text(size = 15, colour = "black"),
    axis.text.x = element_text(size = 15, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 17, colour = "black", face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = 17, colour = "black", face = "plain", margin = margin(t = 10)),
    strip.text = element_text(size = 15)
  ) 

#### DEGs ####

de_comparisons <- list(
  #TP53_KO_vs_WT = c("genotype", "TP53_KO", "WT"),
  MGA_KO_vs_WT = c("genotype", "MGA_KO", "WT"),
  #Double_KO_vs_WT = c("genotype", "Double_KO", "WT"),
  #Double_KO_vs_TP53_KO = c("genotype", "Double_KO", "TP53_KO")
  Double_KO_vs_MGA_KO = c("genotype", "Double_KO", "MGA_KO")
)

# Extract results for all comparisons
de_res_list <- purrr::map(
  de_comparisons,
  \(contrast) {
    results(
      dds_global_no_S19,
      independentFiltering = TRUE,
      alpha = 0.05,
      lfcThreshold = 0.58,
      altHypothesis = "greaterAbs",
      pAdjustMethod = "BH",
      contrast = contrast,
      parallel = TRUE
    )
  }
)

# LFC shrinked
de_res_list_shrink <- purrr::map(
  de_res_list,
  \(res) lfcShrink(
    dds_global,
    res = res,
    type = "ashr"
  )
)

lapply(de_res_list_shrink, head)


# Checks
stopifnot(
  length(de_res_list) == length(de_res_list_shrink),
  identical(
    names(de_res_list),
    names(de_res_list_shrink)
  )
)


# Processing
de_res_list_processed <- purrr::map2(
  de_res_list,
  de_res_list_shrink,
  \(res, res_shrink) {
    
    process_DE_results(
      res = res,
      res_shrink = res_shrink,
      gene_metadata = gene_info
    )
  }
)

lapply(de_res_list_processed_NO_MGA19, head)
lapply(de_res_list_processed_MGA19, head)

# de_res_list_processed_MGA19 <- de_res_list_processed
# de_res_list_processed_NO_MGA19 <- de_res_list_processed

# class(de_res_list_processed)
#library(writexl)


#write_xlsx(de_res_list_processed, path = "./results/downstream/dge/RNA_seq_results_MGA_TP53_KOs_2026_09_18.xlsx")


# --------- MA PLOT ---------------- # 
ma_plots <- purrr::imap(
  de_res_list_processed,
  \(data, contrast) {
    make_MA_plot(
      data = data,
      title = contrast,
      y_limit = 6
    )
  }
)

ma_plots[[1]]

# ----- log2FC distribution ----- #
lfc_plots <- purrr::imap(
  de_res_list_processed,
  \(data, contrast) {
    make_LFC_plot(
      data = data,
      title = contrast
      )
  }
)


# ----- P-value distribution ----- #
pval_plots <- purrr::imap(
  de_res_list_processed,
  \(data, contrast) {
    make_PVAL_plot(
      data = data,
      title = contrast
    )
  }
)

pval_plots[[5]]


# ---- Volcano plot ---- #
volcano_plots <- purrr::imap(
  de_res_list_processed,
  \(data, contrast) {
    make_VOLCANO_plot(
      data = data,
      title = contrast,
      x_limit = 5,
      y_limit = 15,
      lfc_threshold = 0.58,
      padj_threshold = 0.05,
      label_n = 10
    )
  }
)


# ---- Dot plots ---- #
expression_plots <- make_all_EXPRESSION_panels(
  dds = dds_global,
  de_results = de_res_list_processed,
  comparisons = de_comparisons,
  n = 12,
  panel_size = 6,
  show_boxplot = FALSE,
  sample_metadata_cols = "genotype",
  colors = c(
    "WT" = "steelblue",
    "TP53_KO" = "orange3",
    "MGA_KO" = "forestgreen",
    "Double_KO" = "firebrick"
  )
)



# write.xlsx(cart_uroa_vs_dmso_res, "./results/downstream/differential_expression/DGE_results_CART_UroA_2026_07_06.xlsx")
#write.xlsx(mm_uroa_vs_dmso_res, "./results/downstream/differential_expression/DGE_results_MM_JJN3_UroA_2026_07_06.xlsx")

# Sensitivity analysis
# de_res_list_processed_MGA19 <- de_res_list_processed
# de_res_list_processed_NO_MGA19 <- de_res_list_processed

lapply(de_res_list_processed_MGA19, head)

lapply(de_res_list_processed_NO_MGA19, head)


# Rename results variables for analysis without MGA19
de_res_list_processed_NO_MGA19 <- purrr::map(
  de_res_list_processed_NO_MGA19,
  ~ .x %>% dplyr::select(-gene_name, -Chr) %>%
    dplyr::rename_with(
      ~ paste0(.x, "_NO_MGA19"),
      .cols = baseMean:last_col()
    )
)

# Master dataframe for sensitivity analysis.
master_sens_analysis <- purrr::map(
  names(de_res_list_processed_MGA19),
  \(comparison) {
    
    df_with <- de_res_list_processed_MGA19[[comparison]]
    df_without <- de_res_list_processed_NO_MGA19[[comparison]]
    
    df_with %>%
      dplyr::left_join(
        df_without,
        by = "ENSEMBL_ID"
      )
  }
)
names(master_sens_analysis) <- names(de_res_list_processed_MGA19)

master_sens_analysis <- lapply(master_sens_analysis, function(df) {
  df %>%
    mutate(
      
      degs_sens = case_when(
        differentially_expressed == "YES" & differentially_expressed_NO_MGA19 == "YES" ~ "Commom DEG",
        differentially_expressed == "YES" & differentially_expressed_NO_MGA19 == "NO" ~ "DEG only MGA19",
          differentially_expressed == "NO" & differentially_expressed_NO_MGA19 == "YES" ~ "DEG only no MGA19",
        TRUE ~ "NO"
        ),
      
      log2FC_diff = log2FoldChange - log2FoldChange_NO_MGA19,
      
      regulation_change = case_when(
        DE_category == "Down" & DE_category_NO_MGA19 == "Up" ~ TRUE,
        DE_category == "Up" & DE_category_NO_MGA19 == "Down" ~ TRUE,
        TRUE ~ FALSE
      )
    )
})


# Overlapping
deg_with <- master_sens_analysis$Double_KO_vs_MGA_KO %>%
  filter(differentially_expressed == "YES") %>%
  pull(ENSEMBL_ID)

deg_without <- master_sens_analysis$Double_KO_vs_MGA_KO %>%
  filter(differentially_expressed_NO_MGA19 == "YES") %>%
  pull(ENSEMBL_ID)


s <- master_sens_analysis$MGA_KO_vs_WT %>%
  filter(differentially_expressed == "YES" | differentially_expressed_NO_MGA19 == "YES")


ggVennDiagram(
  list(
    `MGA19-Yes` = deg_with,
    `MGA19-No` = deg_without
  ),
  #label = "count",
  label_alpha = 0,
  set_color = "black",
  set_size = 7,
  label_size = 7,
  relative_width = 0
  ) +
  scale_fill_gradient(low = "#F0F0F0", high = "#56B4E9") +
  theme_void() +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) 


# Correlation
ggplot(master_sens_analysis$MGA_KO_vs_WT, aes(x = log2FoldChange, y = log2FoldChange_NO_MGA19)) +
  geom_point(
    size = 3,
    shape = 21,
    fill = "#3C8DBC",
    color = "white",
    stroke = 0.7,
    alpha = 0.9
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#D55E00",
    fill = "#D55E00",
    alpha = 0.15,
    linewidth = 1
  ) +
  stat_cor(
    method = "pearson",
    #aes(label = paste(..rr.label.., ..p.label.., sep = ", ")),
    #label.x.npc = "left",
    #label.y.npc = "top",
    size = 5
  ) +
  #labs(
   # x = "Variable X",
  #  y = "Variable Y"
  #) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

summary(master_sens_analysis$MGA_KO_vs_WT$log2FC_diff)
summary(master_sens_analysis$Double_KO_vs_MGA_KO$log2FC_diff)

ggplot(master_sens_analysis$Double_KO_vs_MGA_KO, aes(x = log2FC_diff)) +
  geom_density(
    col = "darkblue",
    size = 1
  ) +
  theme_bw(base_size = 20)


library(writexl)

getwd()

write_xlsx(
  list(
    MGA_KO_vs_WT = master_sens_analysis$MGA_KO_vs_WT,
    Double_KO_vs_MGA_KO = master_sens_analysis$Double_KO_vs_MGA_KO
  ),
  "RNA_seq_MGA_TP53_KOs_sensitivity_analysis_2027_08_28.xlsx"
)





# -----------------------------
# GSEA with clusterProfiler
# -----------------------------

# -----------------------------
# 1. Load MSigDB databases
# -----------------------------

# HALLMARKS
hallmarks_data <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_collection == "H") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol)

# kegg_medicus
kegg_medicus <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_subcollection == "CP:KEGG_MEDICUS") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol)

# kegg_legacy
kegg_legacy <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_subcollection == "CP:KEGG_LEGACY") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol)

# GO_BP
go_bp <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_subcollection == "GO:BP") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol, db_version)

# GO_BP
go_mf <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_subcollection == "GO:MF") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol, db_version)

# GO_CC
go_cc <- msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(gs_subcollection == "GO:CC") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol, db_version)

# REACTOME
reactome_bp <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::filter(gs_subcollection == "CP:REACTOME") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol, db_version)


# REACTOME
wikipathways <- msigdbr(species = "Homo sapiens") %>% 
  dplyr::filter(gs_subcollection == "CP:WIKIPATHWAYS") %>%
  dplyr::select(gs_name, db_ncbi_gene, db_gene_symbol, db_version)

# List of pathways to evaluate
pathways_list <- list(
  HALLMARK = hallmarks_data,
  GOBP = go_bp,
  GOMF = go_mf,
  GOCC = go_cc,
  REACTOME = reactome_bp,
  KEGG = kegg_legacy,
  KEGG_MEDICUS = kegg_medicus,
  WP = wikipathways
)

# Standarization
pathways_list <- purrr::map(
  pathways_list,
  \(x) {
    x %>%
      dplyr::select(
        term = gs_name,
        gene = db_ncbi_gene
      ) %>%
      dplyr::distinct()
  }
)

stopifnot(
  all(
    purrr::map_lgl(
      pathways_list,
      \(x) all(c("term", "gene") %in% colnames(x))
    )
  )
)

purrr::map(
  pathways_list,
  \(x) {
    x %>%
      dplyr::summarise(
        n_terms = dplyr::n_distinct(term),
        n_genes = dplyr::n_distinct(gene),
        n_na_genes = sum(is.na(gene))
      )
  }
)

# -----------------------------
# 2. Prepare GSEA ranking data
# -----------------------------

# Hacerlo con de_res_list_processed_MGA19 Y de_res_list_processed_NO_MGA19
# Rehacer ambos si es necesario.

# Get rankings by stat
rankings <- make_gsea_rankings(
  de_res_list = de_res_list_processed_MGA19,
  OrgDb = org.Hs.eg.db,
  gene_type = "ENSEMBL",
  gene_type_col = "ENSEMBL_ID",
  rank_col = "stat"
  )

# -----------------------------------------------------
# 3. Create GSEA grid for all comparisons x pathways
# -----------------------------------------------------

# GSEA comparison grid
gsea_grid <- tidyr::expand_grid(
  comparison = names(rankings),
  pathway_set = names(pathways_list)
) %>%
  dplyr::mutate(
    rank = purrr::map(
      comparison,
      \(x) rankings[[x]]
    ),
    pathways = purrr::map(
      pathway_set,
      \(x) pathways_list[[x]]
    )
  ) %>%
  dplyr::mutate(
    analysis_id = paste(
      comparison,
      pathway_set,
      sep = "__"
    )
  )  %>%
  dplyr::relocate(
    analysis_id,
    comparison,
    pathway_set
  )


# Add results in gsea_grid too.
gsea_grid <- gsea_grid %>%
  dplyr::mutate(
    result = purrr::pmap(
      list(rank, pathways),
      \(rank, pathways) {
        
        run_gsea(
          rank, 
          pathways,
          min_gs_size = 15,
          max_gs_size = 500,
          pvalue_cutoff = 1,
          p_adjust_method = "BH",
          seed = TRUE
        )
      }
    )
  )





# ------------------------------------------------------------------------------
# GSEA RES processed
# ------------------------------------------------------------------------------

gsea_grid <- gsea_grid %>%
  dplyr::mutate(
    result_processed = purrr::pmap(
      list(comparison, pathway_set),
      \(comparison, pathways) {
        
        process_GSEA_results(
          gsea_grid = gsea_grid,
          comparison = comparison,
          pathway_set = pathways,
          OrgDb = org.Hs.eg.db,
          nes_threshold = 1,
          padj_threshold = 0.25,
          set_gene_symbol = TRUE
        )
      }
    )
  )


# gsea_grid_MGA19 <- gsea_grid
#gsea_grid_NO_MGA19 <- gsea_grid



# ------------------------------------------------------------------------------
# GSEA TOP PATHWAYS TO PLOT
# ------------------------------------------------------------------------------

gsea_grid <- gsea_grid %>%
  mutate(
    top_pathways = purrr::map(
      result_processed,
      select_top_pathways,
      n = 5
    )
  )


# ------------------------------------------------------------------------------
# GSEA PLOTS
# ------------------------------------------------------------------------------
gsea_plots <- gsea_grid %>%
  dplyr::select(
    analysis_id,
    comparison,
    pathway_set,
    top_pathways
  ) %>%
  dplyr::mutate(
    barplot = purrr::pmap(
      list(
        top_pathways,
        comparison,
        pathway_set
      ),
      barplot_GSEA
    )
  )

gsea_plots$barplot[31]

# ------------------------------------------------------------------------------
# EXPORT GSEA RES processed tables
# ------------------------------------------------------------------------------
purrr::walk(
  names(de_comparisons),
  ~ export_GSEA_results(
    gsea_grid = gsea_grid,
    comparison = .x,
    output_dir = here("results/downstream/gsea")
  )
)


# ---------- SENSITIVITY ANALYSIS --------- #


gsea_sens_analysis <- purrr::map(
  seq_along(gsea_grid_MGA19$analysis_id),
  \(i) {
    
    df_with <- gsea_grid_MGA19$result_processed[[i]] #%>%
      #dplyr::select(
       # pathway_name,
      #  enrichmentScore,
       # NES,
      #  p.adjust,
       # differentially_enriched,
      #  enrichment,
      #  DENR_category
      #)
    
    df_without <- gsea_grid_NO_MGA19$result_processed[[i]] %>%
    #dplyr::select(
        #pathway_name,
    #enrichmentScore,
    #NES,
    #p.adjust,
    #differentially_enriched,
    #enrichment,
    #DENR_category
    #  )%>%
      rename_with(~ paste0(.x, "_NO_MGA19"), .cols = -c(pathway_name, -setSize))
    
    df_with %>%
      dplyr::left_join(
        df_without,
        by = "pathway_name"
      ) %>%
      dplyr::mutate(
        
        # Significance concordance
        gsea_sensitivity = case_when(
          differentially_enriched == "YES" & differentially_enriched_NO_MGA19 == "YES" ~ "Common Enriched",
          differentially_enriched == "YES" & differentially_enriched_NO_MGA19 == "NO"  ~ "Enriched only MGA19",
          differentially_enriched == "NO"  & differentially_enriched_NO_MGA19 == "YES" ~ "Enriched only no MGA19",
          TRUE ~ "NO"
        ),
        
        # Change in NES
        NES_diff = NES - NES_NO_MGA19,
        

        # Direction
        direction_change = case_when(
          NES > 0 & NES_NO_MGA19 < 0 ~ TRUE,
          NES < 0 & NES_NO_MGA19 > 0 ~ TRUE,
          TRUE ~ FALSE
        )
      )
  }
)

write_xlsx(gsea_sens_analysis, path="./results/downstream/gsea/GSEA_sens_analysis.xlsx")


# Summary

gsea_overlap_summary <- map_dfr(
  gsea_sens_analysis,
  \(df) {
    #tab <- table(df$gsea_sensitivity)
    
    tibble(
      with_MGA19 = sum(df$differentially_enriched == "YES"),
      without_MGA19 = sum(df$differentially_enriched_NO_MGA19 == "YES"),
      #common = tab["Common Enriched"],
      common = sum(df$gsea_sensitivity == "Common Enriched"),
      #overlap_pct = 100 * tab["Common Enriched"] / with_MGA19
      overlap_pct = 100 * common / with_MGA19
    )
  },
  .id = "analysis_id"
)


# Correlation plots

plot_correlation <- function(df, title) {
  
  df <- df %>% dplyr::filter(gsea_sensitivity == "Common Enriched")
  
  ggplot(df, 
         aes(x = enrichmentScore, y = enrichmentScore_NO_MGA19)) +
    
    geom_point(
      size = 3,
      shape = 21,
      fill = "#3C8DBC",
      color = "white",
      stroke = 0.7,
      alpha = 0.9
    ) +
    
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "#D55E00",
      fill = "#D55E00",
      alpha = 0.15,
      linewidth = 1
    ) +
    
    stat_cor(
      method = "spearman",
      size = 7
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    
    labs(
      title = title,
      x = "NES with MGA19",
      y = "NES without MGA19"
    ) +
    
    theme_bw(base_size = 26) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
}

gsea_sensitivity_plots <- imap(
  gsea_sens_analysis,
  ~ plot_correlation(.x, .y)
)


gsea_sensitivity_plots[[11]]

# NES Direction
table(gsea_sens_analysis$MGA_KO_vs_WT__KEGG$direction_change, gsea_sens_analysis$MGA_KO_vs_WT__KEGG$gsea_sensitivity)
                       
gsea_overlap_summary <- map_dfr(
  gsea_sens_analysis,
  \(df) {
    #tab <- table(df$gsea_sensitivity)
    
    tibble(
      change = sum(df$direction_change),
      no_change = sum(!df$direction_change),
      common_total = sum(df$gsea_sensitivity == "Common Enriched"),
      common_change = sum(df$gsea_sensitivity == "Common Enriched"& df$direction_change),
      common_no_change = sum(df$gsea_sensitivity == "Common Enriched"& !df$direction_change)
    )
  },
  .id = "analysis_id"
)




plot_nes_diff <- function(df, title) {
  
  ggplot(df, aes(x = NES_diff)) +
    geom_density(
      col = "darkblue",
      size = 1
    ) +
    theme_bw(base_size = 20)  +
    
    labs(
      title = title,
      x = "NES DIFF (MGA19 - NO MGA19)",
    ) +
    
    theme_bw(base_size = 26) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
}



nes_diff_plots <- imap(
  gsea_sens_analysis,
  ~ plot_nes_diff(.x, .y)
)

nes_diff_plots[[16]]
