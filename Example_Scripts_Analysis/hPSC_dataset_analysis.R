#' @title Embryo Dataset Analysis
#'
#' @description
#' 
#' Central File where the analysis of a scRNA-seq data from an Embryo dataset is analyzed 
#' 
#' @author Pedro Granjo
#' @date 13-03-2026
#' 
#' 



#Now the working directory will be the folder that this RScript is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#I'm assuming that run file cleaned Claudia and the other document Functions Processing and the other one Petropolous are all
# in the same folder
#' @title Installation of missing packages
#'
#' @description
#'  Installs required packages that are not currently installed 
#' 
#' @param pkgs packages that you need for your analysis
#' @param installer type of installation, if it is from Biocondutor e.g(BiocManager::install) or cran
#' 
#' 
install_if_missing <- function(pkgs, installer) {
  
  missing <- pkgs[!pkgs %in% rownames(installed.packages())]
  
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    installer(missing)
  }
}

source("~/GitHub/Inferring_aneuploidy/R/Functions_Processing_InferCNV.R")
source("~/GitHub/Inferring_aneuploidy/R/Score_system.R")
source("~/GitHub/Inferring_aneuploidy/R/GSVA.R")
source("~/GitHub/Inferring_aneuploidy/R/GSEA.R")


############################################################################################################-
############################# Across and Within Cell type Integration Approach##############################
############################################################################################################-

test_run <- run_full_cnv_pipeline(
  precomputed = "C:/Users/pmgra/Documents/VUB/InferCNV/Edouard/new_data/",
  start_from = "block2",
  save_intermediate = T,
  outdir            = "C:/Users/pmgra/Documents/VUB/Experimental_code/test_output/",
  counts_mx         = counts_mx,
  metadata          = metadata,
  cell_type_col     = "cell_type",
  gene_order_file   = "~/VUB/InferCNV/InferCNV_RScripts/hg38_gencode_v27.txt",
  mode              = "within",
  chr_exclude       = c("MT", "Y"),
  min_max_counts    = c(100, 1e6),
  n_splits_within   = 3,
  base_outdir       = "C:/Users/pmgra/Documents/VUB/InferCNV/Edouard/new_data/",
  cutoff            = 0.1,
  cluster_by_groups = TRUE,
  HMM               = FALSE,
  denoise           = TRUE,
  analysis_mode     = "subclusters",
  window_length     = 140,
  no_plot           = TRUE,
  resume_if_exists  = TRUE,
  base_dir                              = NULL,
  modes                                 = c("within"),
  tool                                  = "infercnv",
  pattern                               = "^run\\.final",
  max_gap                               = 100000,
  min_overlap_consistent_calls          = 0.75,
  min_overlap_multiple_nodes            = 0.6,
  filter_seq_mb_init                    = 5,
  filter_seq_mb_equiv                   = 7,
  min_references                        = 2,
  overlap_method_equiv_cnv_call_merge   = "reciprocal",
  overlap_method_equiv_cnv_after_filter = "reciprocal",
  parallel                              = FALSE,
  cores                                 = 1L,
  clique_mode_consistent                = "connected",
  removed_log_return                    = FALSE,
  chromosome_arms   = chromosome_arms,
  group_cols        = "cell_type",
  cell_col          = "cell_name",
  chr_col           = "chr",
  start_col         = "start",
  end_col           = "end",
  by                          = c("cell_type", "cnv_state"),
  sample_col                  = "cell_type",
  overlap_method              = "reciprocal",
  min_overlap                 = 0.8,
  boundaries_mb               = c(25, 10),
  base_fraction               = 0.05,
  step                        = 0.05,
  min_cap_threshold           = 5L,
  max_cap_threshold           = 25L,
  total_chromosome_permission = 65
) 



test <- process_tool_cnv_runs(
    base_dir = "C:/Users/pmgra/Documents/VUB/InferCNV/Edouard/new_data/",
    modes                                = c("within"),
    tool                                 = "infercnv",
    pattern                              = "^run\\.final",
    max_gap                              = 100000,
    min_overlap_consistent_calls         = 0.75,
    min_overlap_multiple_nodes           = 0.6,
    filter_seq_mb_init                  = 5,
    filter_seq_mb_equiv                  = 7,
    min_references                       = 2,
    overlap_method_equiv_cnv_call_merge  = "reciprocal",
    overlap_method_equiv_cnv_after_filter = "reciprocal",
    parallel                             = FALSE,
    cores                                = 1L,
    clique_mode_consistent               = "connected",
    removed_log_return                   = FALSE,
    metadata = metadata
)




save(test, file = "C:/Users/pmgra/Documents/VUB/FWO/2026_03_FWO_Claudia/D7_hPSC/test.RData")
load("C:/Users/pmgra/Documents/VUB/FWO/2026_03_FWO_Claudia/D7_hPSC/test.RData")

hepsc <- test

load("C:/Users/pmgra/Documents/VUB/InferCNV/chromossome_arms.RData")

lapply(hepsc[["within"]], function(x){find("cnvs_supported_overlaped",names(x))})

f_df <- lapply(hepsc[["within"]], function(x) {
  x[["cnvs_supported_overlaped"]]
}) %>% bind_rows()


cnv_total_within_df <- add_chromosome_info(f_df,
                                 chromosome_arms,
                                 chr_col = "chr",
                                 start_col = "start",
                                 end_col = "end") 

clustered_events_within <- run_cnv_locus_analysis(cnv_total_within_df, by = c("cell_type"), overlap_method = "reciprocal", min_ovelap = 0.8,removed_log_retur = T, sample_col = "cell_type", cell_col = "cell_name" )



celltype_sizes <- metadata %>%
  group_by(cell_type) %>%
  summarise(
    n_total_cells = n())




res_scores_within <- score_cnv_clusters(
  summary_df = clustered_events_within$cnv_locus_summary,
  clustered_events = clustered_events_within$clustered_events,
  cell_sizes = celltype_sizes,
  by_union             = "cell_type",
  boundaries_mb        = c(25, 10),
  base_fraction        = 0.05,
  step                 = 0.05,
  fractions            = NULL,
  threshold_method     = "auto",
  threshold_mode       = "fractions",
  min_cap_threshold    = 5L,
  max_cap_threshold    = 25L,
  max_tiers            = 2L,
  total_chromosome_permission = 65,
  round_fun = ceiling
)

cnv_filtered <- res_scores_within |>
  dplyr::filter(tier == 1)


metadata_sub <- metadata[!(metadata$cell_name %in% cnv_filtered$cell_name), ]
table(metadata_sub$cell_type)

########################################################################### -
##############################Chromossome Ploting##########################
############################################################################ -

# 1. Prepare genome (once per hg38)
genome_structure <- prepare_genome_structure(chromosome_arms)

# 2. Map to genome
cnv_mapped <- map_cnv_to_genome(
  cnv_filtered,
  genome_structure,
  threshold = 80,
  arrange_df_cols = c("cell_type")
)


heatmap_plot <- prepare_cnv_plot(
    cnv_mapped,
    genome_structure,
    grouping_cols = c("cell_type"),
    state_colors  = c(
      "gain" = "#E64B35",
      "loss" = "royalblue4"))
    
    
plot_cnv_karyotype(
  heatmap_plot,
    genome_structure,
    ideogram_ratio  = 0.05,
    arm_colors      = c(
      "p"   = "paleturquoise4",
      "cen" = "black",
      "q"   = "red4"
    ),
    show_legend     = TRUE)


######################################################################################### -
################################## GSEA Analysis ########################################
######################################################################################### -

seu <- label_karyotype_from_aneu_table(seurat_obj_d7_sub, cnv_filtered,cell_col = "cell_name")
raw_data <- GetAssayData(seu, layer="data")
table(seu$karyotype, useNA = "ifany")


pb_lineage <- make_pseudobulk(
  seu,
  group_vars = c("cell_type", "karyotype"),
  K = 5,
  min_cells_per_group = 25
)

res <- run_edger_by_lineage(pb_counts = pb_lineage$pb_counts, pb_meta = pb_lineage$pb_meta)


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



p <- plot_fgsea_bar(fgsea_by_lineage$Amnion,
               padj_cutoff = 0.05,
               n_max = 10,
               title = "Amnion — Hallmark enrichment (fgsea)",
               show_labels = TRUE) 

ggsave("enrichment_plot_nonneuroectoderm.png",
       p,
       width = 20,
       height = 10,
       dpi = 300)


######################################################################################### -
################################## GSVA Analysis ########################################
######################################################################################### -



test <- readRDS("C:/Users/pmgra/Downloads/run.final.infercnv_obj")

library(Seurat)
genes_annot <- test@gene_order


## ---- 1) Prepare gene annotation (rownames -> gene_id) ----

# gene_annot must have: chr, start, stop (hg38)
genes_df <- genes_annot %>%
  dplyr::mutate(gene_id = rownames(genes_annot)) %>%
  dplyr::filter(gene_id %in% rownames(raw_data)) %>%
  dplyr::transmute(
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



gene_sets <- get_hallmark_sets()
gsva_res <- scgsva_pipeline(seurat_obj_d7_sub, gene_sets,assay = "SCT")

hallmark_cols <- colnames(gsva_res@meta.data)[
  grepl("^HALLMARK_", colnames(gsva_res@meta.data))
]

btotal <- as.vector(burden$B_total)
names(btotal) <- burden$cell

cell_type <- seurat_obj_d7_sub@meta.data[names(btotal), "cell_type", drop = TRUE]


lm_res <- run_hallmark_lm_by_celltype(
  gsva_meta = gsva_res@meta.data,
  hallmark_cols = hallmark_cols,
  btotal = btotal,
  cell_type = cell_type,
  min_cells = 10,
  r2_threshold = 0.10
)


models_by_celltype <- lapply(split(assess, assess$cell_type), function(df) {
  if (nrow(df) < 3) return(NULL)
  lm(p53 ~ btotal, data = df)
})


results_by_celltype <- lapply(names(models_by_celltype), function(ct) {
  mod <- models_by_celltype[[ct]]
  if (is.null(mod)) return(NULL)
  
  sm <- summary(mod)
  coef_tab <- sm$coefficients
  
  data.frame(
    cell_type = ct,
    n_cells = nrow(model.frame(mod)),
    beta_btotal = coef_tab["btotal", "Estimate"],
    se_btotal = coef_tab["btotal", "Std. Error"],
    t_btotal = coef_tab["btotal", "t value"],
    p_btotal = coef_tab["btotal", "Pr(>|t|)"],
    r_squared = sm$r.squared
  )
})

results_by_celltype <- do.call(rbind, results_by_celltype)
results_by_celltype



hallmarks <- grep("HALLMARK", names(gsva_res@meta.data), value = TRUE)

results <- lapply(hallmarks, function(h) {
  assess <- data.frame(
    
    pathway = gsva_res@meta.data[names(btotal), h],
    burden = btotal
  )
  
  assess <- na.omit(assess)
  
  model <- lm(pathway ~ burden, data = assess)
  
  data.frame(
    pathway = h,
    beta = coef(model)[2],
    pvalue = summary(model)$coefficients[2,4],
    r2 = summary(model)$r.squared
  )
})
results <- do.call(rbind, results)



############################################################################################ -
################################## Genome Imbalance ########################################
############################################################################################ -

#need to load GSVA R file for functions

genes_gr <- GRanges(
  seqnames = genes_df$chr,
  ranges   = IRanges(start = genes_df$start, end = genes_df$end),
  gene_id  = genes_df$gene_id
)

genes_gr_red <- GenomicRanges::reduce(genes_gr, ignore.strand = TRUE)
expressed_gene_space_bp <- sum(width(genes_gr_red))

genome_size_bp <- genes_df %>%
  group_by(chr) %>%
  summarise(chr_len = max(end))%>%
  summarise(genome_size_bp = sum(chr_len)) %>%
  pull(genome_size_bp)

genome_burden_expr <- compute_percent_genome_imbalanced(cnv_filtered, genome_size_bp = expressed_gene_space_bp)
genome_burden <- compute_percent_genome_imbalanced(cnv_filtered, genome_size_bp = genome_size_bp)

expression_segments <- genome_burden_expr %>%
  left_join(cnv_filtered %>%select(cell_name, cell_type), by = c("cell" = "cell_name"))

genome_segments <- genome_burden %>%
  left_join(cnv_filtered %>%select(cell_name, cell_type), by = c("cell" = "cell_name"))


head(genome_burden_expr)



library(ggplot2)
cols <- c(
  "#1b9e77",
  "#d95f02",
  "#7570b3",
  "#e7298a"
)

g <- ggplot(expression_segments, aes(x = cell_type, y = pct_genome_imbalanced, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = cell_type), width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "",
    y = "% whole expressed genes imbalanced"
  )

p <- ggplot(genome_segments, aes(x = cell_type, y = pct_genome_imbalanced, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = cell_type), width = 0.15, alpha = 0.5, size = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "",
    y = "% genome imbalanced"
  )      

#ggsave("Edo_whole_genome_imbalance.pdf",
#       p,
#       width = 10,
#       height = 8,
#       dpi = 400)

#ggsave("Edo_genome_expressed_imbalance.pdf",
#       g ,
#       width = 10,
#       height = 8,
#       dpi = 400)
