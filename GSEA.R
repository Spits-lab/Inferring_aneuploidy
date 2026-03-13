#' @title GSEA
#'
#' @description
#' 
#' Functions that perform a GSEA pseudo-bulk RNA-seq analysis
#' from scRNA-seq
#' 
#' @author Pedro Granjo
#' @date 13-03-2026
#' 
#' 


#Set Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading functions developed for the workflow
source("clusters_cell_type_processing_functions.R")

# Package groups
cran_packages <- c(
  "dplyr"
)

bioc_packages <- c(
  "edgeR"
)


# Install missing CRAN and Bioconductor packages
install_if_missing(cran_packages, install.packages)

# Load BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install_if_missing(bioc_packages, BiocManager::install)

# Combine all for loading
all_packages <- c(cran_packages, bioc_packages)

# Load all quietly
invisible(lapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}))


############################################################################## -
############### Karyotype Labelling and Pseudo-Bulk Compression ##############
############################################################################## -

label_karyotype_from_aneu_table <- function(seu, cnv_filtered, cell_col = "cell_name_new") {
  cnv <- as.data.table(cnv_filtered)
  if (!(cell_col %in% names(cnv))) stop("cnv_filtered must contain column: ", cell_col)
  
  aneu_cells <- unique(cnv[[cell_col]])
  all_cells  <- colnames(seu)
  
  # Align to Seurat cells only (drop CNV cells not in object)
  aneu_cells <- intersect(aneu_cells, all_cells)
  
  seu$karyotype <- ifelse(all_cells %in% aneu_cells, "aneuploid", "euploid")
  seu$karyotype <- factor(seu$karyotype, levels = c("euploid", "aneuploid"))
  seu
}


make_pseudobulk <- function(seu,
                            group_vars = c("lineage", "karyotype"),
                            assay = "RNA",
                            slot = "counts",
                            K = 5,
                            min_cells_per_group = 25) {
  counts <- GetAssayData(seu, assay = assay, slot = slot)
  stopifnot(inherits(counts, "dgCMatrix"))
  
  meta <- as.data.table(seu@meta.data, keep.rownames = "cell_id")
  if (!all(group_vars %in% names(meta))) {
    stop("Missing group vars in Seurat metadata: ",
         paste(setdiff(group_vars, names(meta)), collapse = ", "))
  }
  
  meta <- meta[complete.cases(meta[, ..group_vars])]
  meta[, group_id := do.call(paste, c(.SD, sep = " | ")), .SDcols = group_vars]
  
  # filter small groups
  gs <- meta[, .N, by = group_id]
  keep <- gs[N >= min_cells_per_group, group_id]
  meta <- meta[group_id %in% keep]
  if (nrow(meta) == 0) stop("No groups left after filtering by min_cells_per_group.")
  
  # deterministic split within each group
  setorder(meta, group_id, cell_id)
  meta[, rep_id := ((seq_len(.N) - 1L) %% K) + 1L, by = group_id]
  meta[, pb_id := paste0(group_id, " || rep", rep_id)]
  
  pb_ids <- unique(meta$pb_id)
  
  sum_cols_sparse <- function(cell_ids) {
    idx <- match(cell_ids, colnames(counts))
    idx <- idx[!is.na(idx)]
    Matrix::rowSums(counts[, idx, drop = FALSE])
  }
  
  pb_list <- lapply(pb_ids, function(pid) {
    cells <- meta[pb_id == pid, cell_id]
    sum_cols_sparse(cells)
  })
  
  pb_mat <- do.call(cbind, pb_list)
  rownames(pb_mat) <- rownames(counts)
  colnames(pb_mat) <- pb_ids
  pb_counts <- Matrix::Matrix(pb_mat, sparse = TRUE)
  
  pb_meta <- unique(meta[, c(group_vars, "rep_id", "pb_id"), with = FALSE])
  setkey(pb_meta, pb_id)
  
  list(pb_counts = pb_counts, pb_meta = pb_meta)
}


################################################################################# -
#################### GSEA Set up - Cell type based approach #####################
################################################################################# -


run_edger_by_lineage <- function(pb_counts, pb_meta,
                                 min_reps = 2,
                                 cell_type_col = "cell_type",
                                 filter_min_cpm = 1,
                                 filter_min_n = 2) {
  pb_meta <- as.data.table(pb_meta)
  
  cell_type_cols <- colnames(pb_meta) %in% c(cell_type_col)
  lineages <- sort(unique(pb_meta[, ..cell_type_col][[1]]))
  
  res_list <- setNames(vector("list", length(lineages)), lineages)
  
  for (lin in lineages) {
    smeta <- pb_meta[get(cell_type_col) == lin]
    tab <- table(smeta$karyotype)
    
    if (!all(c("euploid", "aneuploid") %in% names(tab))) next
    if (any(tab[c("euploid", "aneuploid")] < min_reps)) next
    
    cols <- smeta$pb_id
    Y <- DGEList(counts = pb_counts[, cols, drop = FALSE])
    
    group <- factor(smeta$karyotype, levels = c("euploid", "aneuploid"))
    design <- model.matrix(~ group)
    
    keep <- rowSums(cpm(Y) >= filter_min_cpm) >= filter_min_n
    Y <- Y[keep, , keep.lib.sizes = FALSE]
    
    Y <- calcNormFactors(Y, method = "TMM")
    Y <- estimateDisp(Y, design, robust = TRUE)
    fit <- glmQLFit(Y, design, robust = TRUE)
    
    qlf <- glmQLFTest(fit, coef = "groupaneuploid")
    tt <- topTags(qlf, n = Inf)$table
    tt <- as.data.table(tt, keep.rownames = "gene")
    
    tt[, `:=`(
      lineage = lin,
      n_pb_euploid = as.integer(tab["euploid"]),
      n_pb_aneuploid = as.integer(tab["aneuploid"]),
      contrast = "aneuploid_vs_euploid"
    )]
    
    res_list[[lin]] <- tt
  }
  
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  de_all <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  list(de_by_lineage = res_list, de_all = de_all)
}



get_hallmark_sets <- function(species = "Homo sapiens") {
  msigdbr(species = species, category = "H") %>%
    as.data.table() %>%
    split(by = "gs_name", keep.by = FALSE) %>%
    lapply(function(dt) unique(dt$gene_symbol))
}



run_fgsea_for_lineage <- function(de_dt, pathways,
                                  minSize = 10, maxSize = 500,
                                  score_col = c("logFC", "FDR")) {
  de_dt <- as.data.table(de_dt)
  
  # Ensure columns exist
  if (!all(score_col %in% names(de_dt))) stop("DE table missing: ", paste(setdiff(score_col, names(de_dt)), collapse = ", "))
  if (!("gene" %in% names(de_dt))) stop("DE table must include `gene` (SYMBOL).")
  
  # Build ranked vector: signed -log10(P)
  p <- -log10(pmax(de_dt$PValue, .Machine$double.xmin))
  ranks <- sign(de_dt$logFC) * p
  names(ranks) <- de_dt$gene
  
  # Remove duplicates (keep max absolute rank)
  # (fgsea expects unique names; duplicates can happen with weird gene symbols)
  ranks_dt <- data.table(gene = names(ranks), rank = as.numeric(ranks))
  ranks_dt <- ranks_dt[!is.na(gene) & gene != ""]
  ranks_dt <- ranks_dt[ , .SD[which.max(abs(rank))], by = gene]
  ranks <- ranks_dt$rank
  names(ranks) <- ranks_dt$gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  fg <- fgsea::fgsea(pathways = pathways, stats = ranks, minSize = minSize, maxSize = maxSize)
  as.data.table(fg)[order(padj, -abs(NES))]
}



plot_fgsea_bar <- function(fgsea_res,
                              padj_cutoff = 0.05,
                              n_max = 20,
                              title = "TE — Hallmark enrichment (fgsea)",
                              show_labels = TRUE) {
  library(dplyr)
  library(ggplot2)
  
  df <- fgsea_res %>%
    filter(!is.na(padj), padj <= padj_cutoff) %>%
    mutate(min_p = .Machine$double.xmin,
           neg_log10_padj = -log10(pmax(padj, min_p))) %>%
    arrange(desc(abs(NES))) %>%       
    slice_head(n = n_max) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        annotate("text", x = 0, y = 0,
                 label = paste0("No significant pathways (padj ≤ ", padj_cutoff, ").")) +
        labs(title = title)
    )
  }
  
  p <- ggplot(df, aes(x = NES, y = pathway, fill = neg_log10_padj)) +
    geom_col(width = 0.8) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "none") +
    labs(
      title = title,
      x = "NES",
      y = NULL,
      fill = expression(-log[10]("padj"))
    )
  
  if (show_labels) {
    p <- p +
      geom_text(
        aes(label = paste0("padj=", signif(padj, 2))),
        hjust = ifelse(df$NES >= 0, -0.05, 1.05),
        size = 3
      ) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5.5, 40, 5.5, 5.5))  # right margin for labels
  }
  
  p
}


