#' Validate metadata for inferCNV pipeline
#'
#' @param metadata     data.frame with at least cell_name and cell_type_col
#' @param counts_mx    raw counts matrix (genes x cells)
#' @param cell_type_col string, name of the column containing cell type labels
#'
#' @return invisibly returns TRUE if all checks pass; stops with message if not
validate_metadata <- function(metadata,
                               counts_mx,
                               cell_type_col,
                               min_cells = 90) {
  
  # ── 1. Required columns present ───────────────────────────────────────────
  required_cols <- c("cell_name", cell_type_col)
  missing_cols  <- setdiff(required_cols,
                            colnames(metadata))
  
  if (length(missing_cols) > 0) {
    stop(
      "metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nRequired: 'cell_name' and '",
      cell_type_col, "'"
    )
  }
  
  # ── 2. cell_name must be unique ───────────────────────────────────────────
  if (any(duplicated(metadata$cell_name))) {
    n_dup <- sum(duplicated(metadata$cell_name))
    stop(
      n_dup, " duplicated cell_name(s) found. ",
      "Each cell must appear exactly once."
    )
  }
  
  # ── 3. cell_name must match colnames of counts_mx ─────────────────────────
  cells_in_meta   <- metadata$cell_name
  cells_in_counts <- colnames(counts_mx)
  
  only_in_meta   <- setdiff(cells_in_meta,
                              cells_in_counts)
  only_in_counts <- setdiff(cells_in_counts,
                              cells_in_meta)
  
  if (length(only_in_meta) > 0) {
    stop(
      length(only_in_meta),
      " cell(s) in metadata not found in counts_mx. ",
      "First few: ",
      paste(head(only_in_meta, 5), collapse = ", ")
    )
  }
  
  if (length(only_in_counts) > 0) {
    warning(
      length(only_in_counts),
      " cell(s) in counts_mx not in metadata ",
      "— will be dropped. First few: ",
      paste(head(only_in_counts, 5), collapse = ", ")
    )
  }
  
  # ── 4. cell_type_col must not have NA ─────────────────────────────────────
  n_na <- sum(is.na(metadata[[cell_type_col]]))
  if (n_na > 0) {
    stop(
      n_na, " NA value(s) in '", cell_type_col, "'. ",
      "All cells must have a cell type label."
    )
  }
  
  # ── 5. Report cell type sizes and filter ──────────────────────────────────
  type_counts <- table(metadata[[cell_type_col]])
  
  # Identify cell types below threshold
  below_min <- names(type_counts)[
    type_counts < min_cells]
  above_min <- names(type_counts)[
    type_counts >= min_cells]
  
  message("Cell type composition:")
  for (ct in names(type_counts)) {
    n    <- type_counts[[ct]]
    flag <- if (n < min_cells)
              sprintf(" [REMOVED: %d < %d minimum]",
                      n, min_cells)
            else ""
    message("  ", ct, ": ", n, " cells", flag)
  }
  
  # Remove cell types below threshold
  if (length(below_min) > 0) {
    message(sprintf(
      "\nRemoving %d cell type(s) with < %d cells: %s",
      length(below_min),
      min_cells,
      paste(below_min, collapse = ", ")
    ))
    
    metadata <- metadata %>%
      dplyr::filter(
        !(!!sym(cell_type_col) %in% below_min)
      )
    
    message(sprintf(
      "Remaining: %d cells across %d cell type(s)",
      nrow(metadata),
      length(above_min)
    ))
  }
  
  # ── 6. At least 2 cell types remaining ───────────────────────────────────
  n_types <- length(unique(
    metadata[[cell_type_col]]))
  
  if (n_types < 2) {
    warning(
      "Only ", n_types,
      " cell type(s) remaining after filtering. ",
      "Across-cell-type comparisons not possible."
    )
  }
  
  if (n_types == 0) {
    stop("No cell types remaining after filtering. ",
         "Lower min_cells threshold.")
  }
  
  # Return filtered metadata
  return(metadata)
}


#' Add A/B/C random split column to metadata for a single cell type
#'
#' Splits cells of one cell type into three roughly equal groups.
#'
#' @param metadata      data.frame (full metadata, not pre-subsetted)
#' @param cell_type_col string, column name for cell type
#' @param cell_type_val string, which cell type to split
#'
#' @return data.frame subset for that cell type with added 'split_group' column
make_splits <- function(metadata, cell_type_col, cell_type_val, n_splits = 3) {
  sub <- metadata[metadata[[cell_type_col]] == cell_type_val, , drop = FALSE]

  if (!is.numeric(n_splits) || n_splits < 2 || n_splits != round(n_splits)) {
    stop("n_splits must be an integer >= 2. Got: ", n_splits)
  }
  n_splits <- as.integer(n_splits)
  
  if (n_splits > 26) {
    stop("n_splits cannot exceed 26 (only 26 capital letters available). ",
         "Got: ", n_splits)
  }
  
  sub <- metadata[metadata[[cell_type_col]] == cell_type_val, , drop = FALSE]
  n   <- nrow(sub)
  
  if (n < n_splits) {
    stop(
      "Cell type '", cell_type_val, "' has only ", n, " cell(s) but ",
      "n_splits = ", n_splits, ". Need at least ", n_splits, " cells."
    )
  }
  
  # ── Generate letter labels ─────────────────────────────────────────────────
  labels <- LETTERS[seq_len(n_splits)]  
  
  # ── Compute group sizes ───────────────────────────────────────────────────
  # Each of the first (n_splits - 1) groups gets floor(n / n_splits) cells.
  # The last group gets whatever remains so the total always equals n.
  base_size  <- floor(n / n_splits)
  group_sizes <- rep(base_size, n_splits)
  group_sizes[n_splits] <- n - base_size * (n_splits - 1L)
  
  # ── Build label vector and shuffle ────────────────────────────────────────
  split_labels <- c()
  for (i in seq_len(n_splits)) {
    split_labels <- c(split_labels, rep(labels[i], group_sizes[i]))
  }

  sub$split_group <- split_labels
  
  size_summary <- paste(
    mapply(function(lbl, sz) sprintf("%s=%d", lbl, sz), labels, group_sizes),
    collapse = " | "
  )
  message(sprintf("  Split '%s' (n=%d, %d groups): %s",
                  cell_type_val, n, n_splits, size_summary))
  return(sub)
}


#' Build annotations data.frame for inferCNV
#'
#' inferCNV requires a single-column data.frame where:
#'   - rownames = cell barcodes
#'   - single column = group label
#'
#' @param cell_names character vector of cell barcodes
#' @param group_labels character vector of group labels (same length)
#'
#' @return data.frame suitable for inferCNV annotations_file argument
build_annotations_df <- function(cell_names, group_labels) {
  
  if (length(cell_names) != length(group_labels)) {
    stop("cell_names and group_labels must have the same length.")
  }
  
  annot <- data.frame(group = group_labels, row.names = cell_names,
                                   stringsAsFactors = FALSE)
  return(annot)
}