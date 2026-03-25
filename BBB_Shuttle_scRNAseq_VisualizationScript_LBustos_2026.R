# ============================================================================ #
# Figure 3 and Supplementary Figures - scRNA-seq Visualization Pipeline
# LBustos - Updated 2026 :)
# ============================================================================ #

# ============================================================================ #
# --------------------------------- Setup ------------------------------------ #
# ============================================================================ #
rm(list = ls())

project_dir <- (NULL)  # <-- set this to your local project folder or wherever the necessry data is.

# Necessary Data for this script:
#   1. SCTint_srat_25-1b23_QC_filtered.rds or version that contains UMAP emmbeddings and metadata
#   2. DESeq2 comparisions_Hprt - per Cell Type.csv
#   3. Pseudobulked expression data_Hprt - per Cell Type.csv
#   4. Pseudobulk Normalized Counts_All Genes - per Cell Type.csv

if (is.null(project_dir)) {
  stop("Please set 'project_dir' to your local project folder before running.")
}

setwd(project_dir)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(purrr)
  library(readr)
  library(scales)
  library(grid)
  library(tibble)
})

# ----------------------------------------------------------------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~ Input files and Output folders ~~~~~~~~~~~~~~~~~~~~ #
# ----------------------------------------------------------------------------- #

#Input files - please make sure the file names are correct
seurat_file <- file.path(project_dir, "SCTint_srat_25-1b23_QC_filtered.rds")
per_cell_de_file <- file.path(project_dir, "Per Cell Data", "DESeq2 comparisions_Hprt - per Cell Type.csv")
per_cell_expr_file <- file.path(project_dir, "Per Cell Data", "Pseudobulked expression data_Hprt - per Cell Type.csv")
pseudobulk_counts_file <- file.path(project_dir, "Per Cell Data", "Pseudobulk Normalized Counts_All Genes - per Cell Type.csv")

# Output directories and folder
out_dir <- file.path(project_dir, "outputs")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")

fig_umap_dir <- file.path(fig_dir, "UMAP")
fig_marker_dir <- file.path(fig_dir, "Markers")
fig_bar_dir <- file.path(fig_dir, "Barplots")
fig_pb_dir <- file.path(fig_dir, "Pseudobulk")

tab_summary_dir <- file.path(tab_dir, "Summary")
tab_pb_dir <- file.path(tab_dir, "Pseudobulk")

dirs_to_make <- c(
  out_dir, fig_dir, tab_dir,
  fig_umap_dir, fig_marker_dir, fig_bar_dir, fig_pb_dir,
  tab_summary_dir, tab_pb_dir
)

invisible(lapply(dirs_to_make, dir.create, recursive = TRUE, showWarnings = FALSE))

# ----------------------------------------------------------------------------- #
# ~~~~~~~~~~~~~~~~~~~~~~~~~ User-set options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ----------------------------------------------------------------------------- #
target_gene <- "Tfrc"   # <- change this to the gene you want for pseudobulk panels
zlim <- 2.5
pb_y_break_step <- 2
prox_locked_y_break_step <- 0.2


# ============================================================================ #
# ~~~~~~~~~~~~~~~ Global labels, Color palettes & Helper functions ~~~~~~~~~~~ #
# ============================================================================ #

# Global Labels
celltype_map <- c(
  "Excitatory neurons" = "Excitatory Neurons",
  "Inhibitory neurons" = "Inhibitory Neurons",
  "Endothelial" = "Endothelial Cells",
  "Endothelial cells" = "Endothelial Cells"
)

cell_levels <- c(
  "Excitatory Neurons",
  "Inhibitory Neurons",
  "Astrocytes",
  "Oligodendrocytes",
  "OPCs",
  "Microglia",
  "Endothelial Cells",
  "Pericytes"
)

# Global Palettes :)
treatment_levels_umap <- c("IV-PBS", "IV-siRNA", "ICV-PBS", "ICV-siRNA")
treat_cols_umap <- c(
  "IV-PBS"    = "#E45F3C",
  "IV-siRNA"  = "#E9943A",
  "ICV-PBS"   = "#79123D",
  "ICV-siRNA" = "#ED1F7F"
)

trt_levels_bar <- c("PBS_IV", "IV", "PBS_ICVb", "ICVb")
short_xlabs_bar <- c("PBS", "siRNA", "PBS", "siRNA")

trt_cols_bar <- c(
  "PBS_ICVb" = "#79123D",
  "PBS_IV"   = "#c06000",
  "ICVb"     = "#ED1F7F",
  "IV"       = "#F57F1F"
)

trt_cols_pb <- c(
  "PBS_ICVb" = "#0F6B99",
  "PBS_IV"   = "#6B990F",
  "ICVb"     = "#7EC3E5",
  "IV"       = "#C3E57E"
)

# Cell-type Heatmap Marker Table

markers <- c(
  "Slc17a7","Slc17a6", "Neurod6", # Excitatory Neuron Cell Annotation
  "Gls","Grin1", "Rbfox3", "Nefl", # Excitatory Neuron Extra Markers
  "Gad1","Gad2","Slc32a1", # Inhibitatory Neuron Cell Annotation
  "Slc6a1", # Inhibitatory Extra Marker
  "Aqp4","Gfap","Aldh1l1","Slc1a3", # Astrocyte Cell Annotation
  "Slc1a2","S100b","Sox9", # Astrocyte Extra Markers
  "Mbp","Mog","Plp1", "Mag", # Oligodenrocyte Cell Annotation
  "Cnp","Galc","Ugt8","Fa2h","Sox10", # Oligodenrocyte Extra Markers
  "Pdgfra","Cspg4","Olig2", # OPC Cell Annotation
  "Olig1", # OPC Extra Marker
  "P2ry12","Tmem119","Cx3cr1", # Mircoglia Cell Annotation
  "Aif1","C1qa",  # Mircoglia Extra Markers
  "Pecam1","Cdh5","Cldn5", # Endothelial Cell Cell Annotation
  "Kdr","Vwf", # Endothelial Cell Extra Markers
  "Pdgfrb","Rgs5","Vtn", # Pericyte Cell Annotation
  "Notch3","Kcnj8","Abcc9" # Pericyte Extra Markers
)

# Helper Functions :)

save_plot_multi <- function(plot, filename_base, folder, width, height, dpi = 600, transparent = FALSE) {
  bg_val <- if (transparent) "transparent" else "white"
  
  ggsave(
    filename = file.path(folder, paste0(filename_base, ".pdf")),
    plot = plot,
    width = width, height = height,
    useDingbats = FALSE,
    bg = bg_val
  )
  
  ggsave(
    filename = file.path(folder, paste0(filename_base, ".png")),
    plot = plot,
    width = width, height = height,
    dpi = dpi,
    bg = bg_val
  )
}

save_plot_tiff <- function(plot, filename_base, folder, width, height, dpi = 600, transparent = FALSE) {
  ggsave(
    filename = file.path(folder, paste0(filename_base, ".tiff")),
    plot = plot,
    width = width, height = height,
    dpi = dpi,
    compression = "lzw",
    bg = if (transparent) "transparent" else "white"
  )
}

nice_up <- function(x, step = 0.1) {
  if (!is.finite(x)) return(0)
  ceiling(x / step) * step
}

p_to_stars <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

clean_celltypes <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- dplyr::recode(x, !!!celltype_map, .default = x)
  factor(x, levels = c(cell_levels, setdiff(sort(unique(x)), cell_levels)))
}

clean_treatment_group_umap <- function(x) {
  x <- stringr::str_squish(as.character(x))
  out <- dplyr::case_when(
    stringr::str_detect(x, "^ICV") & stringr::str_detect(x, "PBS") ~ "ICV-PBS",
    stringr::str_detect(x, "^ICV") & stringr::str_detect(x, "Hprt") ~ "ICV-siRNA",
    stringr::str_detect(x, "^IV")  & stringr::str_detect(x, "PBS") ~ "IV-PBS",
    stringr::str_detect(x, "^IV")  & stringr::str_detect(x, "Hprt") ~ "IV-siRNA",
    TRUE ~ NA_character_
  )
  factor(out, levels = treatment_levels_umap)
}

clean_treatment_group_bar <- function(x) {
  x <- stringr::str_squish(as.character(x))
  out <- dplyr::case_when(
    stringr::str_detect(x, "^IV")  & stringr::str_detect(x, "PBS")  ~ "PBS_IV",
    stringr::str_detect(x, "^IV")  & stringr::str_detect(x, "Hprt") ~ "IV",
    stringr::str_detect(x, "^ICV") & stringr::str_detect(x, "PBS")  ~ "PBS_ICVb",
    stringr::str_detect(x, "^ICV") & stringr::str_detect(x, "Hprt") ~ "ICVb",
    TRUE ~ NA_character_
  )
  factor(out, levels = trt_levels_bar)
}

# ============================================================================ #
# --------- Load Seurat object and clean-up Metadata for nice labels --------- #
# ============================================================================ #

if (!file.exists(seurat_file)) stop("Seurat object not found: ", seurat_file)

obj <- readRDS(seurat_file)
stopifnot(identical(rownames(obj@meta.data), colnames(obj)))

obj$celltype_broad <- clean_celltypes(obj$celltype_broad)
obj$treatment_group <- clean_treatment_group_umap(obj$treatment)

# ================================================================================= #
# --- Figure 3B - scRNA-seq Dataset broad cell type UMAP (labels not included) ---- #
# ================================================================================= #

# NOTE: labels in Figure 3B were added in Affinity Photo based off of Legend generated here

Idents(obj) <- "celltype_broad"

p_umap <- DimPlot(
  obj,
  reduction = "umap",
  label = FALSE,
  repel = TRUE,
  pt.size = 0.05,
  alpha = 0.8
) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

print(p_umap)

save_plot_multi(
  p_umap,
  "UMAP_celltype_broad",
  fig_umap_dir,
  width = 7,
  height = 6,
  dpi = 1200
)

save_plot_tiff(
  p_umap,
  "UMAP_celltype_broad",
  fig_umap_dir,
  width = 7,
  height = 6,
  dpi = 1200
)

# ============================================================================ #
# ---------------- Supp Figure: Cell-type Marker Heatmap --------------------- #
# ============================================================================ #

# NOTE: additional cell-type markers were added to this heatmap than what were used in the cell annotation
#       those additional markers are indicated above in the "markers" table in the Global labels section 

markers_present <- markers[markers %in% rownames(obj)]

DefaultAssay(obj) <- "RNA"

avg <- AverageExpression(
  obj,
  features = markers_present,
  group.by = "celltype_broad",
  assays = DefaultAssay(obj),
  slot = "data",
  verbose = FALSE
)[[1]]

cell_order <- levels(obj$celltype_broad)
cell_order <- cell_order[cell_order %in% colnames(avg)]
avg <- avg[, cell_order, drop = FALSE]

avg_z <- t(scale(t(avg)))
avg_z[is.na(avg_z)] <- 0

heat_df <- as.data.frame(avg_z) |>
  tibble::rownames_to_column("gene") |>
  tidyr::pivot_longer(-gene, names_to = "celltype", values_to = "z") |>
  dplyr::mutate(
    celltype = factor(celltype, levels = cell_order),
    gene = factor(gene, levels = markers_present)
  )

p_heat <- ggplot(heat_df, aes(x = gene, y = celltype, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#497E00",
    mid = "#FAF7F8",
    high = "#C10781",
    midpoint = 0,
    limits = c(-zlim, zlim),
    oob = scales::squish,
    name = "Z-score"
  ) +
  scale_y_discrete(limits = rev(levels(heat_df$celltype))) +
  labs(
    x = NULL,
    y = NULL,
    title = "Cell-Type Annotation Markers"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.box.margin = margin(l = -5),
    legend.margin = margin(0, 0, 0, 0)
  )

print(p_heat)

save_plot_tiff(
  p_heat,
  "Heatmap_Cell_Annotation_Markers",
  fig_marker_dir,
  width = 7,
  height = 2.5,
  dpi = 600,
  transparent = TRUE
)

# ============================================================================ #
# --------------- Supp Figure - UMAP by Treatment Group ---------------------- #
# ============================================================================ #

# NOTE: You can change colors above in the "treat_cols_umap" table above in the global pallette section ;)
#       Not used in the manuscript, but useful!

Idents(obj) <- "treatment_group"

p_treat_overlay <- DimPlot(
  obj,
  reduction = "umap",
  pt.size = 0.2,
  alpha = 0.5
) +
  scale_color_manual(values = treat_cols_umap) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  ggtitle("UMAP by Treatment Group")

p_treat_split <- DimPlot(
  obj,
  reduction = "umap",
  split.by = "treatment_group",
  pt.size = 0.2
) +
  scale_color_manual(values = treat_cols_umap) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

print(p_treat_overlay)
print(p_treat_split)

save_plot_multi(
  p_treat_overlay,
  "UMAP_treatment_overlay",
  fig_umap_dir,
  width = 6,
  height = 6
)

save_plot_multi(
  p_treat_split,
  "UMAP_treatment_split",
  fig_umap_dir,
  width = 18,
  height = 6
)

# ============================================================================ #
# --------------- Supp Fig Cell and Treatment summary tables ----------------- #
# ============================================================================ #

# NOTE: Information that was used to plot in GraphPad for the supplemental figure
#       showing relative proportion of annotated cell types across treatment groups.

meta <- as.data.frame(obj@meta.data)

cells_per_treatment <- meta %>%
  dplyr::count(treatment_group) %>%
  dplyr::arrange(treatment_group)

cells_per_celltype_treatment_pct <- meta %>%
  dplyr::count(treatment_group, celltype_broad) %>%
  dplyr::group_by(treatment_group) %>%
  dplyr::mutate(percent = round(100 * n / sum(n), 2)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(treatment_group, celltype_broad)

readr::write_csv(
  cells_per_treatment,
  file.path(tab_summary_dir, "Cells_per_Treatment_Group.csv")
)

readr::write_csv(
  cells_per_celltype_treatment_pct,
  file.path(tab_summary_dir, "Cells_per_Celltype_per_Treatment_Group_with_Percent.csv")
)

# ============================================================================ #
# -- Figure 3D,E: Gene-of-interest (Hprt1 ONLY) Barplots from per-cell data -- #
# ============================================================================ #

#NOTE: Output is barplots ONLY for Hprt1 expression for each cell type across groups.

if (!file.exists(per_cell_de_file)) stop("Missing file: ", per_cell_de_file)
if (!file.exists(per_cell_expr_file)) stop("Missing file: ", per_cell_expr_file)

res_df <- readr::read_csv(per_cell_de_file, show_col_types = FALSE)
expr_df <- readr::read_csv(per_cell_expr_file, show_col_types = FALSE)

expr_df <- expr_df %>%
  dplyr::mutate(
    cell_type = clean_celltypes(cell_type),
    tg = clean_treatment_group_bar(treatment),
    value = as.numeric(expression)
  )

res_df <- res_df %>%
  dplyr::mutate(
    cell_type = clean_celltypes(cell_type)
  )

sum_df <- expr_df %>%
  dplyr::group_by(cell_type, tg) %>%
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    tg = factor(trt_levels_bar, levels = trt_levels_bar),
    cell_type,
    fill = list(mean = NA_real_, sd = NA_real_)
  )

ann_df <- res_df %>%
  dplyr::mutate(
    comparison = stringr::str_squish(comparison),
    cell_type = as.character(cell_type),
    p = as.numeric(p_val_adj),
    p.label = p_to_stars(p),
    group1 = dplyr::case_when(
      stringr::str_detect(comparison, "^IV") ~ "IV",
      stringr::str_detect(comparison, "^ICV") ~ "ICVb",
      TRUE ~ NA_character_
    ),
    group2 = dplyr::case_when(
      stringr::str_detect(comparison, "^IV") ~ "PBS_IV",
      stringr::str_detect(comparison, "^ICV") ~ "PBS_ICVb",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group1), !is.na(group2)) %>%
  dplyr::select(cell_type, group1, group2, p, p.label)

global_y_upper <- {
  mx <- max(c(sum_df$mean + sum_df$sd, expr_df$value), na.rm = TRUE)
  nice_up(mx * 1.35, step = 0.1)
}

locked_breaks <- seq(0, global_y_upper, by = prox_locked_y_break_step)

plot_bar_one_ct <- function(ct,
                            y_mode = c("free", "locked"),
                            y_upper_locked = global_y_upper,
                            save_file = NULL) {
  y_mode <- match.arg(y_mode)
  
  dsum <- sum_df %>% dplyr::filter(cell_type == ct)
  dpts <- expr_df %>% dplyr::filter(cell_type == ct)
  
  ymax_obs <- max(c(dsum$mean + dsum$sd, dpts$value), na.rm = TRUE)
  if (!is.finite(ymax_obs)) ymax_obs <- 0
  
  y_upper <- if (y_mode == "locked") y_upper_locked else ymax_obs * 1.35
  y_breaks <- if (y_mode == "locked") locked_breaks else pretty(c(0, y_upper), n = 4)
  
  ann <- ann_df %>%
    dplyr::filter(cell_type == ct) %>%
    dplyr::distinct(group1, group2, .keep_all = TRUE)
  
  if (nrow(ann) > 0) {
    base <- max(ymax_obs, y_upper * 0.75)
    gap <- 0.12 * y_upper
    ann$y.position <- base + gap * seq_len(nrow(ann))
    ann$y.position <- pmin(ann$y.position, y_upper * 0.98)
  }
  
  p <- ggplot(dsum, aes(x = tg, y = mean, fill = tg)) +
    geom_col(width = 0.65, colour = "black", linewidth = 0.75, na.rm = TRUE) +
    geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = mean + sd), width = 0.22, na.rm = TRUE) +
    geom_jitter(
      data = dpts,
      aes(x = tg, y = value),
      width = 0.12, height = 0, size = 2.5, alpha = 0.55,
      inherit.aes = FALSE, na.rm = TRUE
    ) +
    scale_x_discrete(drop = FALSE, limits = trt_levels_bar, labels = short_xlabs_bar) +
    scale_fill_manual(values = trt_cols_bar, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    facet_wrap(~cell_type, ncol = 1) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "#CCCCCC", color = "black", linewidth = 1.0),
      strip.text = element_text(face = "bold", color = "black", size = 14),
      axis.text.x = element_text(face = "bold", color = "black", size = 12),
      axis.text.y = element_text(face = "bold", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      panel.spacing = grid::unit(0, "lines"),
      plot.margin = margin(6, 10, 10, 6, "pt")
    ) +
    coord_cartesian(clip = "off") +
    geom_vline(xintercept = 2.5, linetype = "dotted", linewidth = 0.8, color = "black")
  
  if (nrow(ann) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      ann,
      label = "p.label",
      y.position = "y.position",
      tip.length = 0.012,
      hide.ns = TRUE,
      size = 6,
      bracket.size = 0.6
    )
  }
  
  if (!is.null(save_file)) {
    ggsave(
      save_file, p,
      width = 3, height = 4, units = "in",
      dpi = 600, compression = "lzw", device = "tiff"
    )
  }
  
  p
}

dir.create(file.path(fig_bar_dir, "Barplots_Hprt_free"), showWarnings = FALSE, recursive = TRUE)
# FREE = Y axis per cell type is picked based on best for that data set

dir.create(file.path(fig_bar_dir, "Barplots_Hprt_locked"), showWarnings = FALSE, recursive = TRUE)
# LOCKED = Y axis is the same across all cell types

cell_types_bar <- levels(expr_df$cell_type)
cell_types_bar <- cell_types_bar[!is.na(cell_types_bar)]

invisible(lapply(cell_types_bar, function(ct) {
  plot_bar_one_ct(
    ct,
    y_mode = "free",
    save_file = file.path(fig_bar_dir, "Barplots_Hprt_free", paste0(gsub("\\s+", "_", ct), "_free.tiff"))
  )
}))

invisible(lapply(cell_types_bar, function(ct) {
  plot_bar_one_ct(
    ct,
    y_mode = "locked",
    save_file = file.path(fig_bar_dir, "Barplots_Hprt_locked", paste0(gsub("\\s+", "_", ct), "_locked.tiff"))
  )
}))

# ============================================================================ #
# ------------- Supp Figure - Any gene Pseudobulk gene panels (Tfrc) --------- #
# ============================================================================ #

# Note: No statisitical comparisons are represented, these is purely plotting
#       log2 expression data of a specific gene for every cell type across 
#       across treatment groups.
#       This will include barplots for individual cell types and csv files with data.
#       You can change colors of bars in the "trt_cols_pb" table above in the global pallette section.

if (!file.exists(pseudobulk_counts_file)) stop("Missing file: ", pseudobulk_counts_file)

pseudobulk_counts_df <- readr::read_csv(pseudobulk_counts_file, show_col_types = FALSE)
stopifnot("gene" %in% names(pseudobulk_counts_df))

pseudobulk_long_df <- pseudobulk_counts_df %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "colname",
    values_to = "value"
  ) %>%
  dplyr::mutate(colname = stringr::str_squish(colname)) %>%
  tidyr::extract(
    col = colname,
    into = c("sample_id", "cell_type", "treatment_raw"),
    regex = "^([^_]+)_(.+)_([^_]+)$",
    remove = TRUE
  ) %>%
  dplyr::mutate(
    cell_type = clean_celltypes(cell_type),
    treatment_raw = stringr::str_squish(treatment_raw),
    value = as.numeric(value),
    tg = dplyr::case_when(
      stringr::str_detect(treatment_raw, "^IV")  & stringr::str_detect(treatment_raw, "PBS")  ~ "PBS_IV",
      stringr::str_detect(treatment_raw, "^IV")  & !stringr::str_detect(treatment_raw, "PBS") ~ "IV",
      stringr::str_detect(treatment_raw, "^ICV") & stringr::str_detect(treatment_raw, "PBS")  ~ "PBS_ICVb",
      stringr::str_detect(treatment_raw, "^ICV") & !stringr::str_detect(treatment_raw, "PBS") ~ "ICVb",
      TRUE ~ NA_character_
    ),
    tg = factor(tg, levels = trt_levels_bar)
  ) %>%
  dplyr::filter(!is.na(tg), !is.na(cell_type), !is.na(sample_id))

get_gene_y_upper_locked <- function(gene_symbol, step = pb_y_break_step) {
  d <- pseudobulk_long_df %>% dplyr::filter(gene == gene_symbol)
  if (nrow(d) == 0) return(0)
  
  sum_all <- d %>%
    dplyr::group_by(cell_type, tg) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  mx <- max(c(sum_all$mean + sum_all$sd, d$value), na.rm = TRUE)
  if (!is.finite(mx)) mx <- 0
  
  nice_up(mx * 1.35, step = step)
}

plot_one_panel_pb <- function(gene_symbol,
                              ct,
                              y_mode = c("free", "locked"),
                              y_upper_locked = NULL,
                              save_file = NULL) {
  y_mode <- match.arg(y_mode)
  
  d <- pseudobulk_long_df %>%
    dplyr::filter(gene == gene_symbol, cell_type == ct) %>%
    dplyr::mutate(cell_type = droplevels(cell_type))
  
  if (nrow(d) == 0) return(NULL)
  
  sum_df_pb <- d %>%
    dplyr::group_by(cell_type, tg) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      tg = factor(trt_levels_bar, levels = trt_levels_bar),
      fill = list(mean = NA_real_, sd = NA_real_)
    ) %>%
    dplyr::mutate(cell_type = ct)
  
  ymax_obs <- max(c(sum_df_pb$mean + sum_df_pb$sd, d$value), na.rm = TRUE)
  if (!is.finite(ymax_obs)) ymax_obs <- 0
  
  if (y_mode == "locked") {
    if (is.null(y_upper_locked)) stop("y_upper_locked is required when y_mode = 'locked'")
    y_upper <- y_upper_locked
    y_breaks <- seq(0, y_upper, by = pb_y_break_step)
  } else {
    y_upper <- ymax_obs * 1.35
    y_breaks <- pretty(c(0, y_upper), n = 4)
  }
  
  p <- ggplot(sum_df_pb, aes(x = tg, y = mean, fill = tg)) +
    geom_col(width = 0.65, colour = "black", linewidth = 0.75, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
      width = 0.22, na.rm = TRUE
    ) +
    geom_jitter(
      data = d,
      aes(x = tg, y = value),
      width = 0.12, height = 0, size = 2.5, alpha = 0.55,
      inherit.aes = FALSE, na.rm = TRUE
    ) +
    scale_x_discrete(drop = FALSE, limits = trt_levels_bar, labels = short_xlabs_bar) +
    scale_fill_manual(values = trt_cols_pb, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    facet_wrap(~cell_type, ncol = 1) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "#CCCCCC", color = "black", linewidth = 1.0),
      strip.text = element_text(face = "bold", color = "black", size = 14),
      axis.text.x = element_text(face = "bold", color = "black", size = 12),
      axis.text.y = element_text(face = "bold", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      panel.spacing = grid::unit(0, "lines"),
      plot.margin = margin(6, 10, 10, 6, "pt")
    ) +
    coord_cartesian(clip = "off") +
    geom_vline(xintercept = 2.5, linetype = "dotted", linewidth = 0.8, color = "black")
  
  if (!is.null(save_file)) {
    ggsave(
      save_file, p,
      width = 3, height = 4, units = "in",
      dpi = 600, compression = "lzw", device = "tiff"
    )
  }
  
  p
}

save_gene_panels <- function(gene_symbol,
                             out_dir,
                             y_mode = c("free", "locked")) {
  y_mode <- match.arg(y_mode)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  y_upper_locked <- NULL
  if (y_mode == "locked") {
    y_upper_locked <- get_gene_y_upper_locked(gene_symbol, step = pb_y_break_step)
  }
  
  cts <- levels(pseudobulk_long_df$cell_type)
  cts <- cts[!is.na(cts)]
  
  invisible(lapply(cts, function(ct) {
    safe_ct <- gsub("[^A-Za-z0-9]+", "_", ct)
    f <- file.path(out_dir, paste0(gene_symbol, "_", safe_ct, ".tiff"))
    
    plot_one_panel_pb(
      gene_symbol = gene_symbol,
      ct = ct,
      y_mode = y_mode,
      y_upper_locked = y_upper_locked,
      save_file = f
    )
  }))
}

export_gene_summary_csv <- function(gene_symbol, out_file) {
  d <- pseudobulk_long_df %>% dplyr::filter(gene == gene_symbol)
  if (nrow(d) == 0) stop("Gene not found: ", gene_symbol)
  
  sum_df_gene <- d %>%
    dplyr::group_by(gene, cell_type, tg) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      n = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(cell_type, tg)
  
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(sum_df_gene, out_file)
  sum_df_gene
}

save_gene_panels(
  target_gene,
  out_dir = file.path(fig_pb_dir, paste0("GenePanels_", target_gene, "_free")),
  y_mode = "free"
)

save_gene_panels(
  target_gene,
  out_dir = file.path(fig_pb_dir, paste0("GenePanels_", target_gene, "_locked")),
  y_mode = "locked"
)

export_gene_summary_csv(
  target_gene,
  out_file = file.path(tab_pb_dir, paste0(target_gene, "_summary_values.csv"))
)

# ----------------------------------- #
##### Session info for GitHub :) ###### 
# ---------------------------------- #

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(as.character(Sys.time()), file.path(out_dir, "run_timestamp.txt"))

cat("Figure 3 scRNA-seq visualization pipeline complete.\n")
