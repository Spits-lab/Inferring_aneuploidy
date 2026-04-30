# =============================================================================
# infercnv_make_objects.R
# Functions to create inferCNV objects for within and across comparisons
# =============================================================================

BASE <- "C:/Users/pmgra/Documents/GitHub/Inferring_aneuploidy/R/"

source(file.path(BASE, "Infercnv_utils.R"))


# =============================================================================
# INTERNAL: single object builders
# =============================================================================

#' Build one within-celltype inferCNV object
#'
#' Reference group = one of A/B/C
#' Query          = the other two groups
#'
#' @param counts_mx       genes x cells raw count matrix
#' @param split_metadata  metadata subset for one cell type, with split_group col
#' @param ref_group       "A", "B", or "C"
#' @param gene_order_file path to gene order file (hg38/mm10 etc.)
#' @param chr_exclude     chromosomes to exclude (default c("MT","Y"))
#' @param min_max_counts  c(min, max) counts per cell filter
#'
#' @return inferCNV object
.build_within_object <- function(counts_mx,
                                 split_metadata,
                                 ref_group,
                                 gene_order_file,
                                 chr_exclude    = c("MT", "Y"),
                                 min_max_counts = c(100, 1e6)) {
  
  cells     <- split_metadata$cell_name
  sub_counts <- counts_mx[, cells, drop = FALSE]
  
  annot <- build_annotations_df(
    cell_names   = cells,
    group_labels = split_metadata$split_group
  )
  
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix         = sub_counts,
    annotations_file          = annot,
    gene_order_file           = gene_order_file,
    chr_exclude               = chr_exclude,
    ref_group_names           = ref_group,
    min_max_counts_per_cell   = min_max_counts
  )
  
  return(obj)
}


#' Build one across-celltype inferCNV object
#'
#' Query  = all cells of query_type
#' Reference = all cells of ref_type (full cell type, no splitting)
#'
#' @param counts_mx       genes x cells raw count matrix
#' @param metadata        full metadata data.frame
#' @param cell_type_col   column name for cell type labels
#' @param query_type      cell type string for the query
#' @param ref_type        cell type string for the reference
#' @param gene_order_file path to gene order file
#' @param chr_exclude     chromosomes to exclude
#' @param min_max_counts  c(min, max) counts per cell filter
#'
#' @return inferCNV object
.build_across_object <- function(counts_mx,
                                 metadata,
                                 cell_type_col,
                                 query_type,
                                 ref_type,
                                 gene_order_file,
                                 chr_exclude    = c("MT", "Y"),
                                 min_max_counts = c(100, 1e6)) {
  
  # Subset metadata to only query + reference cell type
  sub_meta <- metadata[metadata[[cell_type_col]] %in% c(as.character(query_type), as.character(ref_type)), ,
                       drop = FALSE]
  
  cells      <- sub_meta$cell_name
  sub_counts <- counts_mx[, cells, drop = FALSE]
  
  # Labels are just the cell type names directly
  annot <- build_annotations_df(
    cell_names   = cells,
    group_labels = sub_meta[[cell_type_col]]
  )
  
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix         = sub_counts,
    annotations_file          = annot,
    gene_order_file           = gene_order_file,
    chr_exclude               = chr_exclude,
    ref_group_names           = ref_type,
    min_max_counts_per_cell   = min_max_counts
  )
  
  return(obj)
}


# =============================================================================
# INTERNAL: mode-level builders
# =============================================================================

#' Build all within-celltype objects across all cell types
#'
#' Returns nested list: list[[cell_type]][[ref_group]] = inferCNV object
#'
.build_all_within <- function(counts_mx,
                               metadata,
                               cell_type_col,
                               gene_order_file,
                               chr_exclude,
                               min_max_counts,
                               n_splits_within) {
  
  cell_types         <- unique(metadata[[cell_type_col]])
  all_split_metadata <- list()
  
  message("\n── Building WITHIN objects ──────────────────────────────────────")
  message(sprintf("   n_splits_within = %d  →  refs: %s",
                  n_splits_within,
                  paste(LETTERS[seq_len(n_splits_within)],
                        collapse = ", ")))
  
  objects <- setNames(
    lapply(cell_types, function(ct) {
      
      number_of_cells <- nrow(
        metadata[metadata[[cell_type_col]] == ct, ])
      
      if (number_of_cells >= 100) {
        
        message("\nCell type: ", ct)
        
        # Split cells for this cell type
        split_meta <- make_splits(
          metadata      = metadata,
          cell_type_col = cell_type_col,
          cell_type_val = ct,
          n_splits      = n_splits_within
        )
        
        # Store split registry
        all_split_metadata[[ct]] <<- split_meta
        
        # One inferCNV object per reference group
        refs <- LETTERS[seq_len(n_splits_within)]
        
        objs <- setNames(
          lapply(refs, function(ref) {
            
            others <- setdiff(refs, ref)
            message(sprintf("  Building ref=%s vs (%s)",
                            ref,
                            paste(others,
                                  collapse = "+")))
            
            tryCatch(
              .build_within_object(
                counts_mx       = counts_mx,
                split_metadata  = split_meta,
                ref_group       = ref,
                gene_order_file = gene_order_file,
                chr_exclude     = chr_exclude,
                min_max_counts  = min_max_counts
              ),
              error = function(e) {
                warning(sprintf(
                  "Failed for %s ref=%s: %s",
                  ct, ref, conditionMessage(e)))
                NULL
              }
            )
          }), refs)
        
        Filter(Negate(is.null), objs)  # ← return value for if block
        
      } else {
        
        message(sprintf(
          "\nCell type: %s skipped — low cells (%d)",
          ct, number_of_cells))
        
        cell_types <<- cell_types[
          !(cell_types == ct)]
      }
      
    }), cell_types)
  
  # Combine split metadata
  split_metadata_combined <- do.call(
    rbind, all_split_metadata)
  rownames(split_metadata_combined) <- NULL
  
  return(list(
    objects        = objects,
    split_metadata = split_metadata_combined
  ))
}


#' Build all across-celltype objects across all cell types
#'
#' Returns nested list: list[[query_type]][[ref_type]] = inferCNV object
#'
.build_all_across <- function(counts_mx,
                              metadata,
                              cell_type_col,
                              gene_order_file,
                              chr_exclude,
                              min_max_counts) {
  
  cell_types <- unique(metadata[[cell_type_col]])
  message("\n── Building ACROSS objects ──────────────────────────────────────")
  
  result <- setNames(lapply(cell_types, function(query) {
    
    # Reference = every other cell type
    ref_types <- setdiff(cell_types, query)
    message("\nQuery: ", query, " | References: ",
            paste(ref_types, collapse = ", "))
    
    objs <- setNames(lapply(ref_types, function(ref) {
      
      message("  Building query=", query, " vs ref=", ref)
      
      tryCatch(
        .build_across_object(
          counts_mx       = counts_mx,
          metadata        = metadata,
          cell_type_col   = cell_type_col,
          query_type      = query,
          ref_type        = ref,
          gene_order_file = gene_order_file,
          chr_exclude     = chr_exclude,
          min_max_counts  = min_max_counts
        ),
        error = function(e) {
          warning("Failed to build across object for query=", query,
                  " ref=", ref, ": ", conditionMessage(e))
          NULL
        }
      )
      
    }), ref_types)
    
    Filter(Negate(is.null), objs)
    
  }), cell_types)
  
  return(result)
}


# =============================================================================
# MAIN: exported function
# =============================================================================

#' Create all inferCNV objects for within and/or across comparisons
#'
#' @param counts_mx       genes x cells raw count matrix
#'                        (e.g. from GetAssayData(seurat_obj, layer = "counts"))
#' @param metadata        data.frame with required columns:
#'                          - 'cell_name': must match colnames(counts_mx)
#'                          - cell_type_col: cell type labels
#' @param cell_type_col   string, name of column in metadata with cell types
#'                        (default "cell_type")
#' @param gene_order_file path to inferCNV gene order file
#'                        (e.g. hg38_gencode_v27.txt)
#' @param mode            one of "within", "across", or "both" (default "both")
#' @param chr_exclude     chromosomes to exclude (default c("MT","Y"))
#' @param min_max_counts  c(min, max) counts per cell (default c(100, 1e6))
#'
#' @return named list with elements 'within_cell_type' and/or 'across_cell_type'
#'         each containing nested lists of inferCNV objects:
#'           within_cell_type[[cell_type]][[ref_group]]   (ref_group: A/B/C)
#'           across_cell_type[[query_type]][[ref_type]]
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(infercnv)
#'
#' counts_mx <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
#'
#' metadata <- data.frame(
#'   cell_name = colnames(seurat_obj),
#'   cell_type = seurat_obj$cell_type
#' )
#'
#' obj_list <- make_infercnv_objects(
#'   counts_mx       = counts_mx,
#'   metadata        = metadata,
#'   cell_type_col   = "cell_type",
#'   gene_order_file = "/path/to/hg38_gencode_v27.txt",
#'   mode            = "both"
#' )
#'
#' saveRDS(obj_list, "infercnv_objcomp.rds")
#' }
make_infercnv_objects <- function(counts_mx,
                                  metadata,
                                  cell_type_col   = "cell_type",
                                  gene_order_file,
                                  mode            = "both",
                                  chr_exclude     = c("MT", "Y"),
                                  min_max_counts  = c(100, 1e6),
                                  n_splits_within) {
  
  # ── Input checks ────────────────────────────────────────────────────────
  if (!mode %in% c("within", "across", "both")) {
    stop("mode must be one of: 'within', 'across', 'both'. Got: '", mode, "'")
  }
  
  if (!file.exists(gene_order_file)) {
    stop("gene_order_file not found: ", gene_order_file)
  }
  
  if (!is.matrix(counts_mx) && !inherits(counts_mx, "dgCMatrix")) {
    stop("counts_mx must be a matrix or dgCMatrix (sparse matrix).")
  }
  
  # ── Validate metadata ────────────────────────────────────────────────────
  message("Validating metadata...")
  metadata   <- validate_metadata(metadata, counts_mx, cell_type_col)
  
  # Align counts to metadata (drop cells not in metadata)
  keep_cells <- intersect(colnames(counts_mx), metadata$cell_name)
  counts_mx  <- counts_mx[, keep_cells, drop = FALSE]
  metadata   <- metadata[metadata$cell_name %in% keep_cells, , drop = FALSE]
  
  # ── Build objects ────────────────────────────────────────────────────────
  result <- list()
  
  if (mode %in% c("within", "both")) {
    result$within_cell_type <- .build_all_within(
      counts_mx       = counts_mx,
      metadata        = metadata,
      cell_type_col   = cell_type_col,
      gene_order_file = gene_order_file,
      chr_exclude     = chr_exclude,
      min_max_counts  = min_max_counts,
      n_splits_within
    )
  }
  
  if (mode %in% c("across", "both")) {
    result$across_cell_type <- .build_all_across(
      counts_mx       = counts_mx,
      metadata        = metadata,
      cell_type_col   = cell_type_col,
      gene_order_file = gene_order_file,
      chr_exclude     = chr_exclude,
      min_max_counts  = min_max_counts
    )
  }
  
  # ── Summary ──────────────────────────────────────────────────────────────
  message("\n── Summary ───────────────────────────────────────────────────────")
  
  if (!is.null(result$within_cell_type)) {
    total_within <- sum(sapply(result$within_cell_type$objects, length))
    message("Within objects built: ", total_within,
            " (", length(result$within_cell_type$objects), " cell types x ",length(result$within_cell_type$objects[[1]]),  " refs)")
  }
  
  if (!is.null(result$across_cell_type)) {
    total_across <- sum(sapply(result$across_cell_type, length))
    message("Across objects built: ", total_across)
  }
  
  message("\nDone. Save with: saveRDS(obj_list, 'infercnv_objcomp.rds')")
  
  return(result)
}
