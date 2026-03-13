
#Now the working directory will be the folder that this RScript is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#I'm assuming that run file cleaned Claudia and the other document Functions Processing and the other one Petropolous are all
# in the same folder

source("Functions_Processing_InferCNV.R")
source("Score_system.R")
source("GSVA.R")
source("GSEA.R")
###############################################################################
############################### Distribution of Gain and Loss #################
###############################################################################
res <- readRDS("Petropolous_Karyotyping.rds")

load("C:/Users/pmgra/Documents/VUB/InferCNV/TE_cells_Petroupoulous_02172026/TE_infercnv_all.RData")
supported_events <- FWO_Eduard_hescp_dataD0_2[["cnv_support_tbl"]]
cnv_filtered <- res[["cnv_table_chromossome_filtered"]]
chromosome_arms <- res[["chromosome_arms"]]
seurat_obj <- res[["seura_obj"]]
cnv_total <- res[["cnv_table"]]


info <- data.frame(cell_name = rownames(seurat_obj@meta.data), cell_type = seurat_obj@meta.data$predicted_celltype_singler)
info$cell_name


cnv_total <- cnv_total %>%
  dplyr::left_join(info, by = c("cell_name_new" = "cell_name"))


supported_events <- final_data[["cnv_support_tbl"]]
#supported_events <- res[["supported_events"]]

#columns like cnv_length, cnv_length_mb to the table
supported_events <- add_additional_columns(supported_events)


## OVerall CNVs Length
plot_all_cnv_distributions(supported_events,c(5, 25, 50))


#################################################################################
############################# Chromossome Arm Info ##############################
################################################################################

library(readr)
load("C:/Users/pmgra/Documents/VUB/InferCNV/chromossome_arms.RData")

#Merge cnv_overlap info with the information that we know about the chromossomes
cnv_arm_class <- add_chromosome_info(supported_events,
                                     chromosome_arms,
                                     chr_col = "chr",
                                     start_col = "start",
                                     end_col = "end") 



#Add stats about whole chromossome and arm into the df total
cnv_total <- calculate_cnv_arm_percentages(
  cnv_arm_class,
  chromosome_arms
)




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

table(plot_long[plot_long$percentage > 80,"event"])


##########################################################################
############################# Statistics #################################
##########################################################################


# Create the three plots to see distributions of CNV % based on chromossome
p1 <- make_density_plot("Whole chromosome", plot_long)
p2 <- make_density_plot("p arm", plot_long)
p3 <- make_density_plot("q arm",plot_long)

# Stack vertically
library(patchwork)

final_plot <- p1 | p2 | p3

ggsave("cnv_density_plots_gap_80k_overlap_075_mean_sd_ed.png",
       final_plot,
       width = 17,
       height = 8,
       dpi = 300)


#########################################################################
################## Subseting based on Chromossomal arm ##################
#########################################################################



all_events <- cnv_total %>%
  mutate(
    event_id = row_number(),
    ds_cell = paste(cell_type, cell_name, sep = "|"),
    event_width = end - start + 1
  )



clustered_events <- run_cnv_locus_analysis(all_events, by = c("embryo"), min_recip = 0.80, cluster_mode = "complete", sample_col = "embryo")

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




table(res_scores$confidence)

#Info with two list elements stage - how many embryos compared to seurat and embryo level - how many cells per embryo 
summary_info <- make_embryo_stage_summary(cnv_filtered, seurat_obj)


###########################################################################
##############################Chromossome Ploting##########################
############################################################################

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

library(gridExtra)

combined_plot_raw <- assemble_heatmap_with_ideogram(heatmap_plot_raw,
                                                    ideogram_plot,
                                                    cnv_mapped)


combined_plot <- assemble_heatmap_with_ideogram(heatmap_plot,
                                                ideogram_plot,
                                                cnv_mapped)

heatmaps_combined <- plot_grid(
  combined_plot_raw,
  combined_plot,
  ncol = 1,
  align = "v",
  labels = c("A", "B")
)


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


ggsave(
  filename = "C:/Users/pmgra/Documents/VUB/InferCNV/Petropoulous__2016/CNV_final_figure_Petropoulos2016.pdf",
  plot = final_grob,
  width = 14,
  height = 10,
  dpi = 300
)


cnv_filtered <- cnv_filtered[cnv_filtered$cell_type =="TE",, drop = F]

seu <- label_karyotype_from_aneu_table(seurat_obj, cnv_filtered,cell_col = "cell_name_new")


table(seu$karyotype, useNA = "ifany")

seu <- subset(seu, predicted_celltype_singler == "TE")

pb_lineage <- make_pseudobulk(
  seu,
  group_vars = c("predicted_celltype_singler", "karyotype"),
  K = 4,
  min_cells_per_group = 25
)

res <- run_edger_by_lineage(pb_counts = pb_lineage$pb_counts, pb_meta = pb_lineage$pb_meta, cell_type_col ="predicted_celltype_singler")


hallmark_sets <- get_hallmark_sets()

fgsea_by_lineage <- lapply(
  names(res$de_by_lineage),
  function(lin) {
    
    de_tbl <- res$de_by_lineage[[lin]]
    
    if (nrow(de_tbl) == 0) return(NULL)
    
    fg <- run_fgsea_for_lineage(
      de_tbl,
      hallmark_sets,
      minSize = 10,
      maxSize = 500
    )
    
    fg$lineage <- lin
    fg
  }
)

names(fgsea_by_lineage) <- names(res$de_by_lineage)



sv_plot <- plot_fgsea_bar(fgsea_by_lineage$TE,
               padj_cutoff = 0.05,
               n_max = 10,
               title = "TE — Hallmark enrichment",
               show_labels = TRUE) 

ggsave(plot= sv_plot, filename = "C:/Users/pmgra/Documents/VUB/InferCNV/Petropoulous__2016/enrichmentplot.png",
       width = 21,
       height = 10,
       dpi = 300)


###GSVA





test <- readRDS("C:/Users/pmgra/Downloads/run.final.infercnv_obj")

library(Seurat)
genes_annot <- test@gene_order


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(Matrix)
})

test <- readRDS("C:/Users/pmgra/Downloads/run.final.infercnv_obj")

library(Seurat)
genes_annot <- test@gene_order
raw_data <- GetAssayData(seu, layer="data")

## ---- 1) Prepare gene annotation (rownames -> gene_id) ----

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

# Defensive checks
stopifnot(all(c("chr","start","end","gene_id") %in% names(genes_df)))
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

## ---- 3) Overlaps + unilateral gene coverage filter (>= 80%) ----
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

## ---- 4) Resolve genes overlapping multiple segments in same cell ----
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


## ---- 5) Per-cell burdens ----
burden <- gene_cell_call %>%
  group_by(cell) %>%
  summarise(
    B_total = n_distinct(gene_id),
    B_gain  = n_distinct(gene_id[call > 0]),
    B_loss  = n_distinct(gene_id[call < 0]),
    .groups = "drop"
  )

hist(burden$B_total)

library(msigdbr)


get_hallmark_sets <- function(species = "Homo sapiens") {
  msigdbr(species = species, category = "H") %>%
    as.data.table() %>%
    split(by = "gs_name", keep.by = FALSE) %>%
    lapply(function(dt) unique(dt$gene_symbol))
}

gene_sets <- get_hallmark_sets()
gsva_res <- scgsva_pipeline(seu, gene_sets,assay = "SCT")

hallmark_cols <- colnames(gsva_res@meta.data)[
  grepl("^HALLMARK_", colnames(gsva_res@meta.data))
]

btotal <- as.vector(burden$B_total)
names(btotal) <- gsub("cell_", "", burden$cell)

cell_type <- "TE"



run_hallmark_lm_by_celltype <- function(gsva_meta,
                                        hallmark_cols,
                                        btotal,
                                        cell_type,
                                        min_cells = 10,
                                        r2_threshold = 0.10,
                                        padj_method = "BH") {
  
  results <- vector("list", length(hallmark_cols))
  names(results) <- hallmark_cols
  
  for (hm in hallmark_cols) {
    
    assess <- data.frame(
      cell = names(btotal),
      pathway_score = gsva_meta[names(btotal), hm, drop = TRUE],
      btotal = btotal,
      cell_type = cell_type,
      stringsAsFactors = FALSE
    )
    
    assess <- assess[complete.cases(assess), ]
    
    models_by_celltype <- lapply(split(assess, assess$cell_type), function(df) {
      if (nrow(df) < min_cells) return(NULL)
      
      # optional: skip if no variation in predictor or response
      if (length(unique(df$btotal)) < 2) return(NULL)
      if (length(unique(df$pathway_score)) < 2) return(NULL)
      
      lm(pathway_score ~ btotal, data = df)
    })
    
    res_hm <- lapply(names(models_by_celltype), function(ct) {
      mod <- models_by_celltype[[ct]]
      if (is.null(mod)) return(NULL)
      
      sm <- summary(mod)
      coef_tab <- sm$coefficients
      
      if (!"btotal" %in% rownames(coef_tab)) return(NULL)
      
      data.frame(
        hallmark = hm,
        cell_type = ct,
        n_cells = nrow(model.frame(mod)),
        beta_btotal = coef_tab["btotal", "Estimate"],
        se_btotal = coef_tab["btotal", "Std. Error"],
        t_btotal = coef_tab["btotal", "t value"],
        p_btotal = coef_tab["btotal", "Pr(>|t|)"],
        r_squared = sm$r.squared,
        adj_r_squared = sm$adj.r.squared,
        stringsAsFactors = FALSE
      )
    })
    
    res_hm <- Filter(Negate(is.null), res_hm)
    
    if (length(res_hm) > 0) {
      res_hm <- do.call(rbind, res_hm)
      res_hm$padj_btotal <- p.adjust(res_hm$p_btotal, method = padj_method)
      res_hm$pass_sig <- res_hm$padj_btotal < 0.05
      res_hm$pass_r2 <- res_hm$r_squared >= r2_threshold
      res_hm$pass_both <- res_hm$pass_sig & res_hm$pass_r2
    } else {
      res_hm <- NULL
    }
    
    results[[hm]] <- res_hm
  }
  
  results <- Filter(Negate(is.null), results)
  
  all_results <- if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    data.frame()
  }
  
  list(
    by_hallmark = results,
    all_results = all_results,
    significant_and_good_r2 = subset(all_results, padj_btotal < 0.05 & r_squared >= r2_threshold)
  )
}


lm_res <- run_hallmark_lm_by_celltype(
  gsva_meta = gsva_res@meta.data,
  hallmark_cols = hallmark_cols,
  btotal = btotal,
  cell_type = cell_type,
  min_cells = 10,
  r2_threshold = 0.10
)


df <- data.frame(cell = names(btotal),
                 B_total = btotal,
                 HALLMARK_P53_PATHWAY = gsva_res@meta.data[names(btotal), "HALLMARK_P53_PATHWAY", drop = TRUE],
                 infercnv_group = "TE",
                 stringsAsFactors = FALSE
                 )

df_te <- subset(df, infercnv_group == "TE")
df_te <- df_te[complete.cases(df_te[, c("B_total", "HALLMARK_P53_PATHWAY")]), ]

mod_te <- lm(HALLMARK_P53_PATHWAY ~ B_total, data = df_te)
sm <- summary(mod_te)

label <- paste0(
  "R² = ", round(sm$r.squared, 3),
  "\n p = ", signif(sm$coefficients["B_total", "Pr(>|t|)"], 3)
)
ggplot(df_te, aes(x = B_total, y = HALLMARK_P53_PATHWAY)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = Inf, y = Inf, label = label,
           hjust = 1.1, vjust = 1.5)  +
  theme_bw() +
  labs(
    title = "Trophectoderm",
    x = "CNV burden (B_total)",
    y = "HALLMARK_P53_PATHWAY"
  )




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
seurat.obj <- label_karyotype_from_aneu_table(seurat_obj, cnv_filtered, cell_col = "cell_name_new")

table(seurat_obj$predicted_celltype_singler)

seurat.obj <- add_ploidy_label(seurat.obj, ploidy_col = "karyotype", out_col = "ploidy_status")




p <- plot_umap_focus_TE_Pre_with_aneuploid_ring(
  seurat.obj,
  focus_types = c("TE"),
  ring_color = "#7A0000",   # slightly softer dark red
  ring_stroke = 0.8
)



ggsave(p, filename = "Umap.pdf",  width = 8,  # smaller width
       height = 6)




