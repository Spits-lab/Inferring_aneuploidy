
base_dir <- "C:/Users/pmgra/Documents/VUB/FWO_Claudia/D7_hPSC/final_data/"


main_folder <- c("across", "within")


r <- setNames(
  lapply(main_folder, function(folder){
    
  r_data <- list.files(paste0(base_dir, folder),full.names = T)
  rnames <- list.files(paste0(base_dir, folder))
  setNames(lapply(r_data, function(data_dir){
    load(data_dir) 
    
    if(folder == "across"){
      supported_events <- filter_cnv_events(
        final_data[["cnvs_filtered"]],
        min_references = 2
      )
    } else if(folder == "within"){
      supported_events <- final_data[["cnv_support_tbl"]]
    }
  }), rnames)
}), main_folder)


# ---------- 2) Standardize tables ----------
std_events <- function(tbl, dataset, mode){
  tbl %>%
    mutate(dataset = dataset, mode = mode) %>%
    select(dataset, mode, cell_name, chr, start, end, cnv_state, everything())
}

get_dataset <- function(x) sub("_.*", "", x)

across_tbl <- bind_rows(
  lapply(names(r$across), \(ds)
         std_events(r$across[[ds]], get_dataset(ds), "across")
  )
)

within_tbl <- bind_rows(
  lapply(names(r$within), \(ds)
         std_events(r$within[[ds]], get_dataset(ds), "within")
  )
)
req_cols <- c("dataset","cell_name","chr","start","end","cnv_state","mode")
stopifnot(all(req_cols %in% colnames(across_tbl)))
stopifnot(all(req_cols %in% colnames(within_tbl)))

# ---------- 3) Pre-filter to shared cells (and shared datasets) ----------
# This removes events from cells that do not exist in both modes.
within_tbl <- within_tbl %>% mutate(ds_cell = paste(dataset, cell_name, sep="|"))
across_tbl <- across_tbl %>% mutate(ds_cell = paste(dataset, cell_name, sep="|"))

shared_ds_cell <- intersect(unique(within_tbl$ds_cell), unique(across_tbl$ds_cell))

within_f <- within_tbl %>% filter(ds_cell %in% shared_ds_cell)
across_f <- across_tbl %>% filter(ds_cell %in% shared_ds_cell)



library(GenomicRanges)
library(dplyr)

# 1) Convert tables -> GRanges, keeping a row id
to_gr <- function(df, prefix){
  stopifnot(all(c("chr","start","end","cell_name","cnv_state") %in% colnames(df)))
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start, end = df$end),
    row_id   = seq_len(nrow(df)),
    cell_name = df$cell_name,
    cnv_state = df$cnv_state,
    prefix    = prefix
  )
}
make_grp <- function(df){
  paste(df$dataset, df$cell_name, df$chr, df$cnv_state, sep="|")
}

within_f <- within_f %>% mutate(grp = make_grp(.), within_row = row_number())
across_f <- across_f %>% mutate(grp = make_grp(.), across_row = row_number())

common_grps <- intersect(unique(within_f$grp), unique(across_f$grp))
within_f <- within_f %>% filter(grp %in% common_grps)
across_f <- across_f %>% filter(grp %in% common_grps)

within_idx <- split(seq_len(nrow(within_f)), within_f$grp)
across_idx <- split(seq_len(nrow(across_f)), across_f$grp)


# ---------- 5) Overlap within each group + reciprocal overlap filter ----------
min_recip <- 0.75

pairs_list <- vector("list", length(common_grps))
names(pairs_list) <- common_grps

for(g in common_grps){
  wi <- within_idx[[g]]
  ai <- across_idx[[g]]
  if(length(wi) == 0 || length(ai) == 0) next
  
  gr_w <- GRanges(
    seqnames = within_f$chr[wi],
    ranges   = IRanges(within_f$start[wi], within_f$end[wi]),
    within_row = within_f$within_row[wi]
  )
  
  gr_a <- GRanges(
    seqnames = across_f$chr[ai],
    ranges   = IRanges(across_f$start[ai], across_f$end[ai]),
    across_row = across_f$across_row[ai]
  )
  
  hits <- findOverlaps(gr_w, gr_a, ignore.strand = TRUE)
  if(length(hits) == 0) next
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  # Vectorized reciprocal overlap
  inter <- width(pintersect(gr_w[q], gr_a[s]))
  ok <- (inter / width(gr_w[q]) >= min_recip) & (inter / width(gr_a[s]) >= min_recip)
  
  if(any(ok)){
    pairs_list[[g]] <- tibble(
      within_row = mcols(gr_w)$within_row[q[ok]],
      across_row = mcols(gr_a)$across_row[s[ok]]
    )
  }
}

pairs_common <- bind_rows(pairs_list)

# Join back full rows
pairs_common <- pairs_common %>%
  left_join(within_f %>% select(within_row, everything()), by = "within_row") %>%
  left_join(across_f %>% select(across_row, everything()), by = "across_row",
            suffix = c("_within", "_across"))

# (Optional extra safety; should already be true due to grp)
pairs_common <- pairs_common %>%
  filter(
    dataset_within == dataset_across,
    cell_name_within == cell_name_across,
    chr_within == chr_across,
    cnv_state_within == cnv_state_across
  )

# ---------- 6) Define within-only and across-only ----------
within_only <- within_f %>%
  anti_join(pairs_common %>% select(within_row) %>% distinct(), by = "within_row") %>%
  select(-within_row)

across_only <- across_f %>%
  anti_join(pairs_common %>% select(across_row) %>% distinct(), by = "across_row") %>%
  select(-across_row)

# ---------- 7) Summaries ----------
summary_counts <- list(
  n_across_total = nrow(across_tbl),
  n_within_total = nrow(within_tbl),
  n_across_filtered = nrow(across_f),
  n_within_filtered = nrow(within_f),
  n_common_pairs = nrow(pairs_common),
  n_within_only = nrow(within_only),
  n_across_only = nrow(across_only)
)

summary_by_dataset <- pairs_common %>%
  count(dataset_within, cnv_state_within, chr_within, name = "n_common") %>%
  arrange(desc(n_common))

# ---------- 8) Save final object ----------
final_obj <- list(
  params = list(min_recip = min_recip, min_references_across = 2),
  across_tbl = across_tbl %>% select(-ds_cell),
  within_tbl = within_tbl %>% select(-ds_cell),
  across_filtered = across_f %>% select(-ds_cell, -grp),
  within_filtered = within_f %>% select(-ds_cell, -grp),
  pairs_common = pairs_common %>% select(-ds_cell_within, -ds_cell_across),
  within_only = within_only %>% select(-ds_cell, -grp),
  across_only = across_only %>% select(-ds_cell, -grp),
  summary_counts = summary_counts,
  summary_by_dataset = summary_by_dataset
)

out_file <- file.path(base_dir, sprintf("integrated_across_within_common_minRecip%.2f.rds", min_recip))
saveRDS(final_obj, out_file)

message("Saved RDS: ", out_file)
message("Summary: ", paste(names(summary_counts), unlist(summary_counts), sep="=", collapse=" | "))





pairs_common <- integrated_across_within_common_minRecip0.75[["pairs_common"]][,grep("_within",colnames(integrated_across_within_common_minRecip0.75[["pairs_common"]])),,drop=F]
colnames(pairs_common) <- gsub("_within", "", colnames(pairs_common))

final_df <- bind_rows(integrated_across_within_common_minRecip0.75$within_only, integrated_across_within_common_minRecip0.75$across_only, pairs_common)

#columns like cnv_length, cnv_length_mb to the table
supported_events <- add_additional_columns(final_df)


## OVerall CNVs Length
plot_all_cnv_distributions(supported_events,c(5, 25, 50))


#################################################################################
############################# Chromossome Arm Info ##############################
################################################################################

library(readr)
load("C:/Users/pmgra/Documents/VUB/InferCNV/chromossome_arms.RData")

#Merge cnv_overlap info with the information that we know about the chromossomes

integrated_res <- readRDS("~/VUB/FWO_Claudia/integrated_across_within_common_minRecip0.75.rds")


pairs_common_clean <- integrated_res$pairs_common %>%
  select(ends_with("_within"),references_across) %>%
  rename_with(~ sub("_within$", "", .x)) %>%
  mutate(mode = "both") %>%
  select(-grp)
  


supported_events <- bind_rows(integrated_res$across_only,
                              integrated_res$within_only,
pairs_common_clean)


supported_events <- add_additional_columns(supported_events)



cnv_arm_class <- add_chromosome_info(supported_events,
                                     chromosome_arms,
                                     chr_col = "chr",
                                     start_col = "start",
                                     end_col = "end") 



#Add stats about whole chromossome and arm into the df total
cnv_arm_class <- cnv_arm_class %>%
  filter(chr %in% chromosome_arms$chr)

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




##########################################################################
############################# Statistics #################################
##########################################################################


# Create the three plots to see distributions of CNV % based on chromossome
p1 <- make_density_plot("Whole chromosome", plot_long, threshold = 80)
p2 <- make_density_plot("p arm", plot_long, threshold = 80 )
p3 <- make_density_plot("q arm",plot_long, threshold = 80)

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
    ds_cell = paste(dataset, cell_name, sep = "|"),
    event_width = end - start + 1
  )



clustered_events <- run_cnv_locus_analysis(all_events, min_recip = 0.8, cluster_mode = "complete",sample_col = "dataset")

celltype_sizes <- celltype_sizes <- c(
  Undifferentiated = 437,
  Neuroectoderm = 1099,
  `Non neural ectoderm` = 412,
  Amnion = 102
)



min_thr <- resolve_celltype_thresholds(
  celltype_sizes,
  method = "fraction",
  fraction = 0.05,
  min_threshold = 3,
  max_threshold = 25
  )



high_thr <- min_thr
colnames(high_thr)[colnames(high_thr) == "min_cells_keep"] <- "high_min_cells"

colnames(high_thr)[colnames(high_thr) == "cell_type"] <- "dataset"
colnames(min_thr)[colnames(min_thr) == "cell_type"] <- "dataset"

res_scores <- score_cnv_clusters(
  clustered_events$cnv_locus_summary,
  clustered_events$clustered_events,
  min_threshold_df = min_thr,
  high_threshold_df = high_thr,
  threshold_group_col = "dataset",
  min_length_mb = 25
)

cnv_filtered <- res_scores %>% 
  filter(confidence == "high")


length(unique(cnv_filtered$cell_name[cnv_filtered$dataset == "Neuroectoderm"]))

###########################################################################
##############################Chromossome Ploting##########################
############################################################################




library(dplyr)



cnv_filtered <- cnv_filtered %>% arrange(dataset,chr)

# 1. Prepare genome (once per hg38)
genome_structure <- prepare_genome_structure(chromosome_arms)


# 2. Map to genome
cnv_mapped <- map_cnv_to_genome(
  cnv_filtered,
  #cnv_filtered,
  genome_structure,
  threshold = 80
)



# 4. Compute embryo boundaries (once)

cnv_mapped <- cnv_mapped %>%
  arrange(dataset,chr) %>%   # must match your plot order
  mutate(plot_idx = row_number())

cnv_mapped <- cnv_mapped %>%
  mutate(dataset = factor(dataset))

cnv_mapped$dataset
dataset_boundaries <- cnv_mapped %>%
  group_by(dataset) %>%
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


get_hallmark_sets <- function(species = "Homo sapiens") {
  msigdbr(species = species, category = "H") %>%
    as.data.table() %>%
    split(by = "gs_name", keep.by = FALSE) %>%
    lapply(function(dt) unique(dt$gene_symbol))
}

gene_sets <- get_hallmark_sets()
gsva_res <- scgsva_pipeline(seurat_obj_d7_sub, gene_sets,assay = "SCT")

hallmark_cols <- colnames(gsva_res@meta.data)[
  grepl("^HALLMARK_", colnames(gsva_res@meta.data))
]

btotal <- as.vector(burden$B_total)
names(btotal) <- burden$cell

cell_type <- seurat_obj_d7_sub@meta.data[names(btotal), "cell_type", drop = TRUE]



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


compute_percent_genome_imbalanced <- function(cnv_df, genome_size_bp = 3.1e9) {
  
  required_cols <- c("cell_name", "chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(cnv_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  segs_df <- cnv_df %>%
    transmute(
      cell = as.character(cell_name),
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    filter(!is.na(cell), !is.na(chr), !is.na(start), !is.na(end), end >= start)
  
  split_segs <- split(segs_df, segs_df$cell)
  
  out <- lapply(names(split_segs), function(cl) {
    x <- split_segs[[cl]]
    
    gr <- GRanges(
      seqnames = x$chr,
      ranges = IRanges(start = x$start, end = x$end)
    )
    
    gr_red <- reduce(gr, ignore.strand = TRUE)
    bp_imbalanced <- sum(width(gr_red))
    
    data.frame(
      cell = cl,
      bp_imbalanced = bp_imbalanced,
      frac_genome_imbalanced = bp_imbalanced / genome_size_bp,
      pct_genome_imbalanced = 100 * bp_imbalanced / genome_size_bp
    )
  })
  
  bind_rows(out)
}


############################3
## % genome imbalance

library(dplyr)
library(GenomicRanges)
library(IRanges)
library(purrr)

compute_percent_genome_imbalanced <- function(cnv_df, genome_size_bp = 3.1e9) {
  
  required_cols <- c("cell_name", "chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(cnv_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  segs_df <- cnv_df %>%
    transmute(
      cell = as.character(cell_name),
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    filter(!is.na(cell), !is.na(chr), !is.na(start), !is.na(end), end >= start)
  
  split_segs <- split(segs_df, segs_df$cell)
  
  out <- lapply(names(split_segs), function(cl) {
    x <- split_segs[[cl]]
    
    gr <- GenomicRanges::GRanges(
      seqnames = x$chr,
      ranges = IRanges(start = x$start, end = x$end)
    )
    
    gr_red <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
    bp_imbalanced <- sum(width(gr_red))
    
    data.frame(
      cell = cl,
      bp_imbalanced = bp_imbalanced,
      frac_genome_imbalanced = bp_imbalanced / genome_size_bp,
      pct_genome_imbalanced = 100 * bp_imbalanced / genome_size_bp
    )
  })
  
  bind_rows(out)
}



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
  left_join(cnv_filtered %>%select(cell_name, dataset), by = c("cell" = "cell_name"))

genome_segments <- genome_burden %>%
  left_join(cnv_filtered %>%select(cell_name, dataset), by = c("cell" = "cell_name"))


head(genome_burden_expr)



library(ggplot2)
cols <- c(
  "#1b9e77",
  "#d95f02",
  "#7570b3",
  "#e7298a"
)

g <- ggplot(expression_segments, aes(x = dataset, y = pct_genome_imbalanced, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = dataset), width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "",
    y = "% whole expressed genes imbalanced"
  )

p <- ggplot(genome_segments, aes(x = dataset, y = pct_genome_imbalanced, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = dataset), width = 0.15, alpha = 0.5, size = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "",
    y = "% genome imbalanced"
  )      

ggsave("Edo_whole_genome_imbalance.pdf",
       p,
       width = 10,
       height = 8,
       dpi = 400)

ggsave("Edo_genome_expressed_imbalance.pdf",
       g ,
       width = 10,
       height = 8,
       dpi = 400)
