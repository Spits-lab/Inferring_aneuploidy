#' @title GSEA
#'
#' @description
#' 
#' Functions that perform a GSEA pseudo-bulk RNA-seq analysis
#' from scRNA-seq
#' 
#' @author Pedro Granjo
#' @date 02-09-2026
#' 
#' 



# Package groups
cran_packages <- c(
  "dplyr","Seurat", "msigdbr"
)

bioc_packages <- c(
  "edgeR","GSVA"
)

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


############################################################################## -
############### Karyotype Labelling and Pseudo-Bulk Compression ##############
############################################################################## -

#' Preprocess a Seurat object for single-cell GSVA
#'
#' Applies SCTransform normalization, principal component analysis (PCA),
#' and UMAP dimensionality reduction to a Seurat object.
#'
#' @param seu A Seurat object that has already been Qced.
#'
#' @return A Seurat object with SCTransform normalization, PCA, and UMAP added.
#'
#' @details
#' This function assumes that QC has already been performed before
#' calling it.
preprocess_seurat <- function(seu) {
  
  #this function assumes that Seurat has been QCed before
  
  seu <- SCTransform(seu, verbose = FALSE)
  
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  
  return(seu)
}


#' Prepare an expression matrix from a Seurat object
#'
#' Extracts expression data from a specified assay in a Seurat object and
#' filters out genes expressed in too few cells.
#'
#' @param seu A Seurat object.
#' @param assay Character string specifying the assay to extract data from.
#' @param min_pct_cells Minimum fraction of cells in which a gene must be
#'   expressed to be retained. 
#'
#' @return A numeric matrix of gene expression values with genes in rows and
#'   cells in columns.
prepare_expression_matrix <- function(seu,
                                      assay = "SCT",
                                      min_pct_cells = 0.05) {
  
  expr_mat <- GetAssayData(seu, assay = assay, slot = "data")
  expr_mat <- as.matrix(expr_mat)
  
  # Remove very sparse genes
  keep_genes <- rowSums(expr_mat > 0) > (min_pct_cells * ncol(expr_mat))
  expr_mat <- expr_mat[keep_genes, ]
  
  return(expr_mat)
}


#' Retrieve MSigDB Hallmark gene sets
#'
#' Downloads Hallmark gene sets from the Molecular Signatures Database (MSigDB)
#' using the \pkg{msigdbr} package and returns them as a list suitable for
#' enrichment analysis tools such as GSVA, ssGSEA, or AUCell.
#'
#' @param species Character string specifying the species for which gene sets
#'   should be retrieved. Default is \code{"Homo sapiens"}.
#'
#' @return A named list of gene sets. Each element corresponds to a Hallmark
#'   pathway and contains a character vector of gene symbols.
get_hallmark_sets <- function(species = "Homo sapiens") {
  msigdbr(species = species, category = "H") %>%
    as.data.table() %>%
    split(by = "gs_name", keep.by = FALSE) %>%
    lapply(function(dt) unique(dt$gene_symbol))
}



#' Run single-cell GSVA on an expression matrix
#'
#' Computes gene set enrichment scores for each cell using the GSVA package.
#'
#' @param expr_mat A numeric expression matrix with genes in rows and cells in columns.
#' @param gene_sets A list of gene sets, where each element is a character vector
#'   of gene symbols.
#' @param method Character string specifying the enrichment method.
#' @param parallel.sz Integer specifying the number of parallel workers.
#' 
#' @return A matrix of gene set enrichment scores with gene sets in rows and
#'   cells in columns.
#'
#'
#' @export
run_scgsva <- function(expr_mat,
                       gene_sets,
                       method = "ssgsea",
                       parallel.sz = 4) {
  
  gsva_res <- gsva(expr = expr_mat,
                   gset.idx.list = gene_sets,
                   method = method,
                   kcdf = "Gaussian",
                   parallel.sz = parallel.sz)
  
  return(gsva_res)
}


#' Add GSVA results to Seurat metadata
#'
#' Transposes a GSVA result matrix and adds the enrichment scores as metadata
#' columns to a Seurat object.
#'
#' @param seu A Seurat object.
#' @param gsva_res A matrix of GSVA enrichment scores with gene sets in rows
#'   and cells in columns.
#'
#' @return A Seurat object with GSVA enrichment scores added to the metadata.
#'
#' @details
#' The GSVA result matrix is transposed before being added so that rows match
#' cells in the Seurat object.
add_gsva_to_seurat <- function(seu, gsva_res) {
  
  gsva_res <- t(gsva_res)
  seu <- AddMetaData(seu, metadata = gsva_res)
  
  return(seu)
}


#' Run a full single-cell GSVA pipeline on a Seurat object
#'
#' Preprocesses a Seurat object, extracts an expression matrix, computes GSVA
#' enrichment scores, and adds them back to the Seurat object's metadata.
#'
#' @param seu A Seurat object.
#' @param gene_sets A list of gene sets, where each element is a character vector
#'   of gene symbols per each Pathway
#' @param assay Character string specifying the assay to use for expression
#'   extraction
#'
#' @return A Seurat object with GSVA enrichment scores added as metadata per cell
scgsva_pipeline <- function(seu,
                            gene_sets,
                            assay = "SCT") {
  
  seu <- preprocess_seurat(seu)
  
  expr_mat <- prepare_expression_matrix(seu, assay = assay)
  
  gsva_res <- run_scgsva(expr_mat, gene_sets)
  
  seu <- add_gsva_to_seurat(seu, gsva_res)
  
  return(seu)
}



#' Run linear models between Hallmark pathway scores and BCR signal per cell type
#'
#' Fits linear models testing the association between Hallmark pathway GSVA scores
#' and BCR signal (\code{btotal}) within each cell type. The function evaluates
#' statistical significance and model fit and returns both per-pathway and
#' aggregated results.
#'
#' @param gsva_meta A data frame or matrix containing GSVA pathway scores, where rows correspond to cells and columns
#'   include Hallmark pathway scores.
#' @param hallmark_cols Character vector specifying the column names corresponding to Hallmark pathway scores.
#' @param btotal Named numeric vector containing the BCR signal values per cell.
#' @param cell_type Character vector specifying the cell type for each cell.
#' @param min_cells Minimum number of cells required within a cell type to fit a model.
#' @param r2_threshold Minimum R-squared value used to flag models with
#'   sufficient explanatory power.
#' @param padj_method Multiple testing correction method used by
#'
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








