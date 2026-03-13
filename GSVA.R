#' @title GSVA
#'
#' @description
#' 
#' Functions that perform a ssGSEA or GSVA on scRNA-seq
#' 
#' @author Pedro Granjo
#' @date 13-03-2026
#' 

#Set Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Loading functions developed for the workflow
source("clusters_cell_type_processing_functions.R")

# Package groups
cran_packages <- c(
  "Seurat", "msigdbr"
)

bioc_packages <- c(
  "GSVA"
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




