#' @title Embryo Dataset Analysis
#'
#' @description
#' 
#' Central File where the analysis of a scRNA-seq data from an Embryo dataset is analyzed 
#' 
#' @author Pedro Granjo
#' @date 23-03-2026
#' 
#' 


#Now the working directory will be the folder that this RScript is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("~/GitHub/Inferring_aneuploidy/R/Functions_Processing_InferCNV.R")
source("~/GitHub/Inferring_aneuploidy/R/Score_system.R")
source("~/GitHub/Inferring_aneuploidy/R/GSVA.R")
source("~/GitHub/Inferring_aneuploidy/R/GSEA.R")


############################################################################### -
############### Initial Processing of InferCNV Calls ##########################
############################################################################### -

ref_dirs <- c("ref1", "ref2", "ref3")
base_dir <- "C:/Users/pmgra/Documents/VUB/InferCNV/TE_cells_Petroupoulous_02172026"
infer_objs <- discover_infercnv_runs(base_dir,ref_dirs, pattern = "^run\\.final")
#infer_objs <- discover_infercnv_runs(base_dir,ref_dirs, pattern = "^17_HMM_.*\\.infercnv_obj$")

infer_objs_1 <- lapply(infer_objs, function(x){
  return(list(x@expr.data, x@gene_order))
})

infer_objs <- load_and_prepare_infercnv_reference(infer_objs_1)


final_data <- run_fast_cnv_pipeline(infer_objs,max_gap = 80000,
                                    min_reciprocal_overlap = 0.75)


############################################################################### -
############################### Distribution of Gain and Loss #################
############################################################################### -

supported_events <- final_data[["cnvs_supported"]]


# Extract embryo (everything except the last dot and cell number)
supported_events$embryo <- sub("^(.+)\\.[0-9]+$", "\\1", df$cell_name)

# Extract stage (first part before the first dot)
supported_events$stage <- sub("^([^\\.]+)\\..*$", "\\1", df$cell_name)


## OVerall CNVs Length
plot_all_cnv_distributions(supported_events,c(5, 25, 50))


################################################################################# -
############################# Chromossome Arm Info ##############################
################################################################################ -

library(readr)
load("C:/Users/pmgra/Documents/VUB/InferCNV/chromossome_arms.RData")

#Merge cnv_overlap info with the information that we know about the chromossomes
cnv_total <- add_chromosome_info(supported_events,
                                     chromosome_arms,
                                     chr_col = "chr",
                                     start_col = "start",
                                     end_col = "end") 




plot_df <- cnv_total %>%
  dplyr::select(
    whole_chromosome_gain,
    whole_chromosome_loss,
    p_arm_gain,
    p_arm_loss,
    q_arm_gain,
    q_arm_loss
  )

# Convert to long format
plot_long <- plot_df %>%
  pivot_longer(cols = everything(),
               names_to = "event",
               values_to = "percentage") %>%
  filter(!is.na(percentage)) %>%
  mutate(
    level = case_when(
      grepl("whole", event) ~ "Whole chromosome",
      grepl("^p_", event)   ~ "p arm",
      grepl("^q_", event)   ~ "q arm"
    ),
    type = case_when(
      grepl("gain", event) ~ "Gain",
      grepl("loss", event) ~ "Loss"
    )
  )
########################################################################## -
############################# Statistics #################################
########################################################################## -

#install.packages("patchwork")
library(patchwork)


# Create the three plots to see distributions of CNV % based on chromossome
p1 <- make_density_plot("Whole chromosome", plot_long)
p2 <- make_density_plot("p arm", plot_long)
p3 <- make_density_plot("q arm",plot_long)



final_plot <- p1 | p2 | p3

#ggsave("cnv_density_plots_gap_80k_overlap_075_mean_sd_ed.png",
#       final_plot,
#       width = 17,
#       height = 8,
#       dpi = 300)


#########################################################################-
################## Filtered and SCoring of CNV segments ##################
#########################################################################-

#Save RDS object
res <- readRDS("Petropolous_Karyotyping.rds")

chromosome_arms <- res[["chromosome_arms"]]
seurat_obj <- res[["seura_obj"]]
cnv_total <- res[["cnv_table"]]


info <- data.frame(cell_name = rownames(seurat_obj@meta.data), cell_type = seurat_obj@meta.data$predicted_celltype_singler)

cnv_total <- cnv_total %>%
  dplyr::left_join(info, by = c("cell_name"))


all_events <- cnv_total %>%
  mutate(
    event_id = row_number(),
    ds_cell = paste(cell_type, cell_name, sep = "|"),
    event_width = end - start + 1
  )


#Analysis of CNV segment frequency based on the embryo
clustered_events <- run_cnv_locus_analysis(all_events, by = c("embryo"), min_recip = 0.80, cluster_mode = "complete", sample_col = "embryo")

#Number of cells detected per embryo
n <- lapply(unique(cnv_total$embryo), function(embryo){
  data.frame(cell_type = embryo, n_total_cells = sum(grepl(embryo,colnames(seurat_obj))))
})
celltype_sizes <- do.call(rbind,n)


min_thr <- resolve_celltype_thresholds(
  celltype_sizes,
  method = "fraction",
  fraction = 0.05,
  min_threshold = 2,
  max_threshold = Inf
)

high_thr <- min_thr
colnames(high_thr)[colnames(high_thr) == "min_cells_keep"] <- "high_min_cells"
colnames(high_thr)[colnames(high_thr) == "cell_type"] <- "embryo"
colnames(min_thr)[colnames(min_thr) == "cell_type"] <- "embryo"


res_scores <- score_cnv_clusters(
  clustered_events$cnv_locus_summary,
  clustered_events$clustered_events,
  min_threshold_df = min_thr,
  by_union = "embryo",
  high_threshold_df = high_thr,
  threshold_group_col = "embryo",
  min_length_mb = 25
)


#Info with two list elements stage - how many embryos compared to seurat and embryo level - how many cells per embryo 
summary_info <- make_embryo_stage_summary(cnv_filtered, seurat_obj)

 
########################################################################### -
##############################Chromossome Ploting##########################
############################################################################ -

library(dplyr)

cnv_filtered<- res_scores %>%
  filter(confidence == "high")


# 1. Prepare genome (once per hg38)
genome_structure <- prepare_genome_structure(chromosome_arms)


# 2. Map to genome
cnv_mapped <- map_cnv_to_genome(
  cnv_filtered,
  genome_structure,
  threshold = 75
)

cnv_mapped <- cnv_mapped %>%
  arrange(embryo, chr) %>%
  mutate(plot_idx = row_number(), cell_type = factor(cell_type), embryo = factor(embryo), chr = factor(chr))

#Separation of CNV segments based on the embryo - Potential Embryo-based patterns are displayed
dataset_boundaries <- cnv_mapped %>%
  group_by(embryo) %>%
  summarise(max_idx = max(plot_idx), .groups = "drop")

boundary_lines <- head(dataset_boundaries$max_idx + 0.5, -1)


# 5. Plot
ideogram_plot <- plot_ideogram(genome_structure,
                               dim(cnv_filtered)[1],
                               arm_colors = c("p"="#4DBBD5",
                                              "cen"="black",
                                              "q"="#E64B35")) 

heatmap_plot <- plot_cnv_heatmap(
  cnv_mapped,
  genome_structure$chromosome_lengths,
  boundary_lines = boundary_lines,
  show_x_labels = T
)


#ggsave(filename="heat_FWO.jpeg",plot = combined_plot, width = 22, height = 10 ,units = "cm",dpi = 300)
#install.packages("gridExtra")
library(gridExtra)


combined_plot <- assemble_heatmap_with_ideogram(heatmap_plot,
                                                ideogram_plot,
                                                cnv_mapped)

legend_cnv <- cowplot::get_legend(
  heatmap_plot + theme(legend.position = "right")
)

legend_arms <- cowplot::get_legend(
  ideogram_plot + theme(legend.position = "right")
)


legend_combined <- plot_grid(
  legend_cnv,
  legend_arms,
  ncol = 1,
  align = "v"
)

final_grob <- arrangeGrob(
  combined_plot,
  legend_combined,
  ncol = 2,
  widths = c(0.9, 0.1) #for the heatmap to have more space than the legends which makes sense
)

library(grid)
grid.newpage()
grid.draw(
  final_grob 
)


#ggsave(
#  filename = "C:/Users/pmgra/Documents/VUB/InferCNV/Petropoulous__2016/CNV_final_figure_Petropoulos2016.pdf",
#  plot = final_grob,
#  width = 14,
#  height = 10,
#  dpi = 300
#)

######################################################################################### -
################################## GSEA Analysis ########################################
######################################################################################### -

# TE is the only cell type ideal for a pseudo-bulk analysis

cnv_filtered <- cnv_filtered[cnv_filtered$cell_type =="TE",, drop = F]

seu <- label_karyotype_from_aneu_table(seurat_obj, cnv_filtered,cell_col = "cell_name")
seu <- subset(seu, predicted_celltype_singler == "TE")


pb_lineage <- make_pseudobulk(
  seu,
  group_vars = c("predicted_celltype_singler", "karyotype"),
  K = 4,
  min_cells_per_group = 25
)

res <- run_edger_by_lineage(pb_counts = pb_lineage$pb_counts, pb_meta = pb_lineage$pb_meta, cell_type_col ="predicted_celltype_singler")


hallmark_sets <- get_hallmark_sets()


fg <- run_fgsea_for_lineage(
      res$de_by_lineage[[1]],
      hallmark_sets,
      minSize = 10,
      maxSize = 500)




sv_plot <- plot_fgsea_bar(fg$TE,
               padj_cutoff = 0.05,
               n_max = 10,
               title = "TE — Hallmark enrichment",
               show_labels = TRUE) 

#ggsave(plot= sv_plot, filename = "C:/Users/pmgra/Documents/VUB/InferCNV/Petropoulous__2016/enrichmentplot.png",
#       width = 21,
#       height = 10,
#       dpi = 300)


######################################################################################### -
################################## GSVA Analysis ########################################
######################################################################################### -





test <- readRDS("C:/Users/pmgra/Downloads/run.final.infercnv_obj")

library(Seurat)
genes_annot <- test@gene_order
raw_data <- GetAssayData(seu, layer="data")

##1) Prepare gene annotation (rownames -> gene_id)

# gene_annot must have: chr, start, stop (hg38)
genes_df <- genes_annot %>%
  mutate(gene_id = rownames(genes_annot)) %>%
  filter(gene_id %in% rownames(raw_data)) %>%
  transmute(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(stop),
    gene_id = as.character(gene_id)
  )

genes_df <- genes_df %>% filter(!is.na(chr), !is.na(start), !is.na(end), end >= start)

genes_gr <- GRanges(
  seqnames = genes_df$chr,
  ranges   = IRanges(start = genes_df$start, end = genes_df$end),
  gene_id  = genes_df$gene_id
)
gene_len <- width(genes_gr)

names(gene_len) <- mcols(genes_gr)$gene_id

## ---- 2) Prepare CNV segments ----
# cnv_filtered must include: cell_name, chr, cnv_state, start, end
segs_df <- cnv_filtered %>%
  transmute(
    cell = as.character(cell_name),
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end),
    cnv_state = cnv_state
  ) %>%
  filter(!is.na(chr), !is.na(start), !is.na(end), end >= start)

segs_gr <- GRanges(
  seqnames = segs_df$chr,
  ranges   = IRanges(start = segs_df$start, end = segs_df$end),
  cell     = segs_df$cell,
  cnv_state= segs_df$cnv_state
)

## 3) Overlaps + unilateral gene coverage filter (>= 80%)
hits <- findOverlaps(genes_gr, segs_gr, ignore.strand = TRUE)
qh <- queryHits(hits)    # gene indices
sh <- subjectHits(hits)  # segment indices

# overlap bp (vectorized)
ov_bp <- width(pintersect(genes_gr[qh], segs_gr[sh]))

# gene coverage fraction = overlap_bp / gene_length (unilateral)
cov_frac <- ov_bp / gene_len[mcols(genes_gr)$gene_id[qh]]

keep <- cov_frac >= 0.80 & ov_bp > 0
qh <- qh[keep]; sh <- sh[keep]; ov_bp <- ov_bp[keep]

# Long table of qualifying overlaps
hit_df <- tibble(
  gene_id = mcols(genes_gr)$gene_id[qh],
  cell    = mcols(segs_gr)$cell[sh],
  state   = mcols(segs_gr)$cnv_state[sh],
  ov_bp   = ov_bp
)

## 4) Resolve genes overlapping multiple segments in same cell
# Rule: among qualifying segments, choose the segment with max ov_bp (winner-take-most)
gene_cell_call <- hit_df %>%
  group_by(gene_id, cell) %>%
  summarise(
    state_bp = state[which.max(ov_bp)],
    ov_bp_max = max(ov_bp),
    .groups = "drop"
  )

# Convert cnv_state to directional call (-1/0/+1). If you want to keep amplitudes, skip sign().
gene_cell_call <- gene_cell_call %>%
  mutate(call = case_when(state_bp == "gain" ~ 1,
                          state_bp =="loss" ~ -1 )) %>%
  filter(call != 0)
gene_cell_call <- gene_cell_call %>%
  filter(gene_id %in% rownames(raw_data))


## Per-cell burdens
burden <- gene_cell_call %>%
  group_by(cell) %>%
  summarise(
    B_total = n_distinct(gene_id),
    B_gain  = n_distinct(gene_id[call > 0]),
    B_loss  = n_distinct(gene_id[call < 0]),
    .groups = "drop"
  )




gsva_res <- scgsva_pipeline(seu, gene_sets,assay = "SCT")

hallmark_cols <- colnames(gsva_res@meta.data)[
  grepl("^HALLMARK_", colnames(gsva_res@meta.data))
]

btotal <- as.vector(burden$B_total)
names(btotal) <- gsub("cell_", "", burden$cell)

cell_type <- "TE"



lm_res <- run_hallmark_lm_by_celltype(
  gsva_meta = gsva_res@meta.data,
  hallmark_cols = hallmark_cols,
  btotal = btotal,
  cell_type = cell_type,
  min_cells = 10,
  r2_threshold = 0.10
)



######################################################################################### -
################################## UMAP Plot ############################################
######################################################################################### -



plot_umap_focus_TE_Pre_with_aneuploid_ring <- function(seurat_obj,
                                                       reduction = "umap",
                                                       ploidy_col = "ploidy_status",
                                                       celltype_col = "predicted_celltype_singler",
                                                       focus_types = c("TE", "Prelineage"),
                                                       base_size = 0.7,
                                                       base_alpha = 0.45,
                                                       focus_size = 0.9,
                                                       focus_alpha = 0.85,
                                                       ring_size = 2.0,
                                                       ring_stroke = 0.7,
                                                       ring_color = "#8B0000",
                                                       other_grey = "darkgrey",
                                                       title = "UMAP: TE/Prelineage highlighted; aneuploid ring") {
  
  stopifnot(reduction %in% names(seurat_obj@reductions))
  
  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  md  <- seurat_obj@meta.data
  
  if (!ploidy_col %in% colnames(md)) stop("ploidy_col not found: ", ploidy_col)
  if (!celltype_col %in% colnames(md)) stop("celltype_col not found: ", celltype_col)
  
  df <- cbind(as.data.frame(emb), md[, c(ploidy_col, celltype_col), drop = FALSE])
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  df$ploidy   <- df[[ploidy_col]]
  df$celltype <- as.character(df[[celltype_col]])
  df <- df[!is.na(df$ploidy) & !is.na(df$celltype), , drop = FALSE]
  
  df$is_focus <- df$celltype %in% focus_types
  
  df$legend_group <- ifelse(df$celltype %in% focus_types,
                            df$celltype,
                            "Other lineages")
  
  
  bg    <- df[!df$is_focus, , drop = FALSE]
  focus <- df[df$is_focus,  , drop = FALSE]
  aneu  <- focus[focus$ploidy == "Aneuploid", , drop = FALSE]
  
  # set explicit colors for focus groups that are actually present
  present_focus <- intersect(focus_types, unique(focus$celltype))
  focus_colors <- c(
    "TE" = "#1f78b4",
    "Prelineage" = "#33a02c"
  )
  focus_colors <- focus_colors[present_focus]
  
  
  
  
  ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2)) +
    
    # All cells colored by legend group
    ggplot2::geom_point(
      ggplot2::aes(color = legend_group),
      size = base_size,
      alpha = base_alpha
    ) +
    
    # Focus cells drawn again on top (slightly stronger)
    ggplot2::geom_point(
      data = focus,
      ggplot2::aes(color = legend_group),
      size = focus_size,
      alpha = focus_alpha
    ) +
    
    # Aneuploid rings
    ggplot2::geom_point(
      data = aneu,
      shape = 21,
      fill = NA,
      color = ring_color,
      size = ring_size,
      stroke = ring_stroke
    ) +
    
    ggplot2::scale_color_manual(
      values = c(
        "Other lineages" = other_grey,
        "TE" = "#1f78b4",
        "Prelineage" = "#33a02c"
      )
    ) +
    theme_minimal()
}
seurat.obj <- label_karyotype_from_aneu_table(seurat_obj, cnv_filtered, cell_col = "cell_name")

table(seurat_obj$predicted_celltype_singler)

seurat.obj <- add_ploidy_label(seurat.obj, ploidy_col = "karyotype", out_col = "ploidy_status")




p <- plot_umap_focus_TE_Pre_with_aneuploid_ring(
  seurat.obj,
  focus_types = c("TE"),
  ring_color = "#7A0000",   # slightly softer dark red
  ring_stroke = 0.8
)



#ggsave(p, filename = "Umap.pdf",  width = 8,  # smaller width
#       height = 6)




