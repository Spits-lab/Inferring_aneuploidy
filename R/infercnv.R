#' @title infercnv
#'
#' @description
#' 
#' Functions that focus on the preparation or the processing of data to be runned through CNV
#' 
#' @author Pedro Granjo
#' @date 02-09-2026
#'



#' Validate metadata for inferCNV pipeline
#'
#' @param metadata     data.frame with at least cell_name and group_col
#' @param counts_mx    raw counts matrix (genes x cells)
#' @param group_col string, name of the column containing cell group labels
#'
#' @return invisibly returns TRUE if all checks pass; stops with message if not
validate_metadata <- function(metadata,
                              counts_mx,
                              group_col,
                              min_cells = 90) {
  
  # ── 1. Required columns present ───────────────────────────────────────────
  required_cols <- c("cell_name", group_col)
  missing_cols  <- setdiff(required_cols,
                           colnames(metadata))
  
  if (length(missing_cols) > 0) {
    stop(
      "metadata is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nRequired: 'cell_name' and '",
      group_col, "'"
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
  
  # ── 4. group_col must not have NA ─────────────────────────────────────
  n_na <- sum(is.na(metadata[[group_col]]))
  if (n_na > 0) {
    stop(
      n_na, " NA value(s) in '", group_col, "'. ",
      "All cells must have a group label."
    )
  }
  
  # ── 5. Report group sizes and filter ──────────────────────────────────
  type_counts <- table(metadata[[group_col]])
  
  # Identify groups below threshold
  below_min <- names(type_counts)[
    type_counts < min_cells]
  above_min <- names(type_counts)[
    type_counts >= min_cells]
  
  message("Groups composition:")
  for (ct in names(type_counts)) {
    n    <- type_counts[[ct]]
    flag <- if (n < min_cells)
      sprintf(" [REMOVED: %d < %d minimum]",
              n, min_cells)
    else ""
    message("  ", ct, ": ", n, " cells", flag)
  }
  
  # Remove groups below threshold
  if (length(below_min) > 0) {
    message(sprintf(
      "\nRemoving %d group(s) with < %d cells: %s",
      length(below_min),
      min_cells,
      paste(below_min, collapse = ", ")
    ))
    
    metadata <- metadata %>%
      dplyr::filter(
        !(!!sym(group_col) %in% below_min)
      )
    
    message(sprintf(
      "Remaining: %d cells across %d group(s)",
      nrow(metadata),
      length(above_min)
    ))
  }
  
  # ── 6. At least 2 groups remaining ───────────────────────────────────
  n_types <- length(unique(
    metadata[[group_col]]))
  
  if (n_types < 2) {
    warning(
      "Only ", n_types,
      " group(s) remaining after filtering. ",
      "Across-cell-type comparisons not possible."
    )
  }
  
  if (n_types == 0) {
    stop("No groups remaining after filtering. ",
         "Lower min_cells threshold.")
  }
  
  # Return filtered metadata
  return(metadata)
}


#' Add A/B/C random split column to metadata for a single groups
#'
#' Splits cells of one cluster into three roughly equal groups.
#'
#' @param metadata      data.frame (full metadata, not pre-subsetted)
#' @param group_col string, column name for cluster/cell type
#' @param subset_group_val string, which clusters to split
#'
#' @return data.frame subset for that clusters with added 'split_group' column
make_splits <- function(
    metadata,
    group_col,
    subset_group_val,
    n_splits  = 3,
    counts_mx,
    seed      = 42
) {
  
  # ── Validate ───────────────────────────────────────────────────────────────

  if (!is.numeric(n_splits) ||
      n_splits < 2 ||
      n_splits != round(n_splits)) {
    stop("n_splits must be an integer >= 2.")
  }
  
  n_splits <- as.integer(n_splits)
  if (n_splits > 26) stop("n_splits cannot exceed 26.")
  
  labels <- LETTERS[seq_len(n_splits)]
  
  # ── Subset metadata ────────────────────────────────────────────────────────
  sub <- metadata[
    metadata[[group_col]] == subset_group_val, ,
    drop = FALSE]
  n          <- nrow(sub)
  cell_names <- sub$cell_name
  
  if (n < n_splits) {
    stop(sprintf(
      "Group '%s' has %d cells but n_splits=%d.",
      subset_group_val, n, n_splits
    ))
  }
  
  target_size <- floor(n / n_splits)
  
  message(sprintf(
    "  Split '%s' (n=%d): target=%d cells/group via PC1",
    subset_group_val, n, target_size
  ))
  
  # ── Subset counts ──────────────────────────────────────────────────────────
  missing <- setdiff(cell_names, colnames(counts_mx))
  if (length(missing) > 0) {
    stop(sprintf("%d cells not in counts_mx.", length(missing)))
  }
  counts_sub <- counts_mx[, cell_names, drop = FALSE]
  
  # ── Seurat pipeline — PCA only ────────────────────────────────────────────
  suppressMessages({
    sobj <- Seurat::CreateSeuratObject(
      counts       = counts_sub,
      min.cells    = 0,
      min.features = 0
    )
    sobj <- Seurat::NormalizeData(sobj,        verbose = FALSE)
    sobj <- Seurat::FindVariableFeatures(sobj, verbose = FALSE)
    sobj <- Seurat::ScaleData(sobj,            verbose = FALSE)
    
    sobj <- Seurat::RunPCA(
      sobj,
      npcs     = 20,
      verbose  = FALSE,
      seed.use = seed
    )
  })

  
  # ── Sort by PC1 and split into equal groups ────────────────────────────────
  pca_coords <- Seurat::Embeddings(sobj, "pca")
  
  # Order cells by PC1
  pc1_order <- order(pca_coords[cell_names, 1])
  
  # Equal group sizes
  group_sizes <- rep(target_size, n_splits)
  group_sizes[n_splits] <- n -
    target_size * (n_splits - 1L)
  
  # Assign labels positionally along PC1
  split_labels <- rep(NA_character_, n)
  
  for (g in seq_len(n_splits)) {
    start_idx <- if (g == 1) 1L else
      sum(group_sizes[seq_len(g - 1)]) + 1L
    end_idx <- sum(group_sizes[seq_len(g)])
    split_labels[pc1_order[start_idx:end_idx]] <- labels[g]
  }
  
  sub$split_group <- split_labels
  
  # ── Validate no NAs ────────────────────────────────────────────────────────
  na_count <- sum(is.na(sub$split_group))
  if (na_count > 0) {
    stop(sprintf(
      "%d cells have NA split_group — check PC1 ordering.",
      na_count
    ))
  }
  
  # ── Report ─────────────────────────────────────────────────────────────────
  final_sizes <- table(sub$split_group)
  final_imbal <- max(as.integer(final_sizes)) -
                 min(as.integer(final_sizes))
  
  size_summary <- paste(
    mapply(
      function(lbl, sz) sprintf("%s=%d", lbl, sz),
      names(final_sizes),
      as.integer(final_sizes)
    ),
    collapse = " | "
  )
  
  message(sprintf(
    "  Done '%s': %s (imbalance=%d cells)",
    subset_group_val,
    size_summary,
    final_imbal
  ))
  
  return(sub)
}


# ── Subfunction 1: PC1-based splitting (original) ─────────────────────────────
.split_by_pc1 <- function(
    sub,
    cell_names,
    pca_coords,
    n_splits,
    labels
) {
  
  n           <- nrow(sub)
  target_size <- floor(n / n_splits)
  pc1_order   <- order(pca_coords[cell_names, 1])
  group_sizes <- rep(target_size, n_splits)
  group_sizes[n_splits] <- n - target_size * (n_splits - 1L)
  split_labels <- rep(NA_character_, n)
  
  for (g in seq_len(n_splits)) {
    start_idx <- if (g == 1) 1L else
      sum(group_sizes[seq_len(g - 1)]) + 1L
    end_idx <- sum(group_sizes[seq_len(g)])
    split_labels[pc1_order[start_idx:end_idx]] <- labels[g]
  }
  
  sub$split_group <- split_labels
  return(sub)
}

# ── Subfunction 2: Donor-aware splitting ──────────────────────────────────────
.split_by_clonal <- function(
    sub,
    cell_names,
    pca_coords,
    n_splits,
    labels,
    clonal_col,
    donor_col = NULL
) {
  
  n           <- nrow(sub)
  target_size <- floor(n / n_splits)
  
  # ── Median PC1 per clonal group ───────────────────────────────────────────
  clonal_summary <- sub %>%
    dplyr::mutate(
      pc1 = pca_coords[cell_name, 1]
    ) %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(
        c(clonal_col,
          if (!is.null(donor_col)) donor_col)
      ))
    ) %>%
    dplyr::summarise(
      n_cells    = dplyr::n(),
      median_pc1 = median(pc1, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::arrange(median_pc1)
  
  cat("  Clonal groups sorted by PC1:\n")
  print(clonal_summary)
  
  n_clonal    <- nrow(clonal_summary)
  split_sizes <- rep(0L, n_splits)
  
  # ── STRICT assignment: one pass, consecutive PC1 ──────────────────────────
  # Each clonal group assigned to exactly ONE split
  # Move to next split only when current is full
  # Last split absorbs all remainder
  # Donor mixing is SECONDARY — only logged, never overrides
  
  clonal_to_split <- rep(NA_character_, n_clonal)
  current_split   <- 1L
  
  # Track donors per split for logging only
  split_donors <- vector("list", n_splits)
  for (s in seq_len(n_splits)) {
    split_donors[[s]] <- character(0)
  }
  
  for (i in seq_len(n_clonal)) {
    
    grp_size  <- clonal_summary$n_cells[i]
    grp_donor <- if (!is.null(donor_col))
      clonal_summary[[donor_col]][i] else NULL
    
    # ── Advance split if current is full ──────────────────────────────────
    # Only advance if not last split
    if (current_split < n_splits &&
        split_sizes[current_split] >= target_size) {
      current_split <- current_split + 1L
      cat(sprintf(
        "  → Advancing to split %s (size reached %d)\n",
        labels[current_split],
        split_sizes[current_split - 1L]
      ))
    }
    
    # ── Always assign to current split ────────────────────────────────────
    # No donor logic can override this
    clonal_to_split[i]          <- labels[current_split]
    split_sizes[current_split]  <- split_sizes[current_split] +
                                    grp_size
    
    # Log donor placement (informational only)
    if (!is.null(grp_donor)) {
      donor_already_in_split <- grp_donor %in%
        split_donors[[current_split]]
      if (donor_already_in_split) {
        cat(sprintf(
          "  ℹ Donor %s already in split %s — ",
          "accepted (clonal integrity takes priority)\n",
          grp_donor, labels[current_split]
        ))
      }
      split_donors[[current_split]] <- c(
        split_donors[[current_split]], grp_donor)
    }
    
    cat(sprintf(
      "  Assigned %s → split %s (split now %d cells)\n",
      clonal_summary[[clonal_col]][i],
      labels[current_split],
      split_sizes[current_split]
    ))
  }
  
  clonal_summary$split_group <- clonal_to_split
  
  # ── Report ────────────────────────────────────────────────────────────────
  cat("\n  Final assignment:\n")
  print(clonal_summary %>%
    dplyr::select(
      dplyr::all_of(
        c(clonal_col,
          if (!is.null(donor_col)) donor_col,
          "n_cells", "median_pc1", "split_group")
      )
    )
  )
  
  cat("\n  Cells per split:\n")
  for (s in seq_along(labels)) {
    rows      <- clonal_summary$split_group == labels[s]
    n_in      <- sum(clonal_summary$n_cells[rows],
                     na.rm = TRUE)
    donors_in <- if (!is.null(donor_col))
      unique(clonal_summary[[donor_col]][rows]) else NULL
    cat(sprintf(
      "  Split %s: %d cells%s\n",
      labels[s], n_in,
      if (!is.null(donors_in))
        paste0(" | donors: ",
               paste(donors_in, collapse = ", "))
      else ""
    ))
  }
  
  # ── Map back to cells ──────────────────────────────────────────────────────
  sub <- sub %>%
    dplyr::left_join(
      clonal_summary %>%
        dplyr::select(
          dplyr::all_of(clonal_col),
          split_group
        ),
      by = clonal_col
    )
  
  # ── Clonal integrity check — FATAL if broken ──────────────────────────────
  broken <- sub %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(clonal_col))
    ) %>%
    dplyr::summarise(
      n_splits_used = dplyr::n_distinct(split_group),
      .groups       = "drop"
    ) %>%
    dplyr::filter(n_splits_used > 1)
  
  if (nrow(broken) > 0) {
    stop(
      "CRITICAL: Clonal integrity VIOLATED for ",
      nrow(broken), " groups: ",
      paste(broken[[clonal_col]], collapse = ", ")
    )
  }
  
  message("  Clonal integrity: ALL ",
          clonal_col, " groups intact ✅")
  
  return(sub)
}


# ── Mother function ───────────────────────────────────────────────────────────
make_splits <- function(
    metadata,
    group_col,
    subset_group_val,
    n_splits   = 3,
    counts_mx,
    seed       = 42,
    clonal_col = NULL,  # if NULL → PC1 only
    donor_col  = NULL   # only used if clonal_col set
) {
  
  # ── Validate ───────────────────────────────────────────────────────────────
  if (!is.numeric(n_splits) ||
      n_splits < 2 ||
      n_splits != round(n_splits)) {
    stop("n_splits must be an integer >= 2.")
  }
  n_splits <- as.integer(n_splits)
  if (n_splits > 26) stop("n_splits cannot exceed 26.")
  labels <- LETTERS[seq_len(n_splits)]
  
  # ── Subset metadata ────────────────────────────────────────────────────────
  sub <- metadata[
    metadata[[group_col]] == subset_group_val, ,
    drop = FALSE
  ]
  n          <- nrow(sub)
  cell_names <- sub$cell_name
  
  if (n < n_splits) {
    stop(sprintf(
      "Group '%s' has %d cells but n_splits=%d.",
      subset_group_val, n, n_splits
    ))
  }
  
  # ── Validate clonal/donor cols ────────────────────────────────────────────
  if (!is.null(clonal_col) &&
      !clonal_col %in% colnames(sub)) {
    stop("clonal_col '", clonal_col,
         "' not found in metadata.")
  }
  if (!is.null(donor_col) &&
      !donor_col %in% colnames(sub)) {
    stop("donor_col '", donor_col,
         "' not found in metadata.")
  }
  
  # ── Subset counts ──────────────────────────────────────────────────────────
  missing <- setdiff(cell_names, colnames(counts_mx))
  if (length(missing) > 0) {
    stop(sprintf("%d cells not in counts_mx.", length(missing)))
  }
  counts_sub <- counts_mx[, cell_names, drop = FALSE]
  
  # ── PCA ───────────────────────────────────────────────────────────────────
  message(sprintf(
    "  Split '%s' (n=%d, mode=%s)",
    subset_group_val, n,
    if (!is.null(clonal_col)) "clonal-aware" else "PC1"
  ))
  
  suppressMessages({
    sobj <- Seurat::CreateSeuratObject(
      counts       = counts_sub,
      min.cells    = 0,
      min.features = 0
    )
    sobj <- Seurat::NormalizeData(sobj,        verbose = FALSE)
    sobj <- Seurat::FindVariableFeatures(sobj, verbose = FALSE)
    sobj <- Seurat::ScaleData(sobj,            verbose = FALSE)
    sobj <- Seurat::RunPCA(
      sobj,
      npcs     = 20,
      verbose  = FALSE,
      seed.use = seed
    )
  })
  
  pca_coords <- Seurat::Embeddings(sobj, "pca")
  
  # ── Dispatch to subfunction ───────────────────────────────────────────────
  sub <- if (!is.null(clonal_col)) {
    
    .split_by_clonal(
      sub        = sub,
      cell_names = cell_names,
      pca_coords = pca_coords,
      n_splits   = n_splits,
      labels     = labels,
      clonal_col = clonal_col,
      donor_col  = donor_col
    )
    
  } else {
    
    .split_by_pc1(
      sub        = sub,
      cell_names = cell_names,
      pca_coords = pca_coords,
      n_splits   = n_splits,
      labels     = labels
    )
  }
  
  # ── Validate no NAs ────────────────────────────────────────────────────────
  na_count <- sum(is.na(sub$split_group))
  if (na_count > 0) {
    stop(sprintf(
      "%d cells have NA split_group.",
      na_count
    ))
  }
  
  # ── Final report ──────────────────────────────────────────────────────────
  final_sizes <- table(sub$split_group)
  final_imbal <- max(as.integer(final_sizes)) -
                 min(as.integer(final_sizes))
  
  size_summary <- paste(
    mapply(
      function(lbl, sz) sprintf("%s=%d", lbl, sz),
      names(final_sizes),
      as.integer(final_sizes)
    ),
    collapse = " | "
  )
  
  message(sprintf(
    "  Done '%s': %s (imbalance=%d cells)",
    subset_group_val,
    size_summary,
    final_imbal
  ))
  
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



# =============================================================================
# INTERNAL: single object builders
# =============================================================================

#' Build one within-celltype inferCNV object
#'
#' Reference group = one of A/B/C
#' Query          = the other two groups
#'
#' @param counts_mx       genes x cells raw count matrix
#' @param split_metadata  metadata subset for one group, with split_group col
#' @param ref_group       "A", "B", or "C"
#' @param gene_order_file path to gene order file (hg38/mm10 etc.)
#' @param chromosomes_to_exclude     chromosomes to exclude (default c("MT","Y"))
#' @param min_max_counts  c(min, max) counts per cell filter
#'
#' @return inferCNV object
.build_within_object <- function(counts_mx,
                                 split_metadata,
                                 ref_group,
                                 gene_order_file,
                                 chromosomes_to_exclude    = c("MT", "Y"),
                                 min_max_counts = c(100, 1e6)) {
  
  
  cells     <- split_metadata$cell_name
  sub_counts <- counts_mx[, cells, drop = FALSE]
  
  annot <- build_annotations_df(
    cell_names   = cells,
    group_labels = split_metadata$split_group
  )
  
  obj <- tryCatch(
    infercnv::CreateInfercnvObject(
      raw_counts_matrix         = sub_counts,
      annotations_file          = annot,
      gene_order_file           = gene_order_file,
      chr_exclude               = chromosomes_to_exclude,
      ref_group_names           = ref_group,
      min_max_counts_per_cell   = min_max_counts
    ),
    error = function(e) {
      message(sprintf(
        "  ERROR in CreateInfercnvObject — ref=%s: %s",
        ref_group,
        conditionMessage(e)
      ))
      NULL
    },
    warning = function(w) {
      message(sprintf(
        "  WARNING in CreateInfercnvObject — ref=%s: %s",
        ref_group,
        conditionMessage(w)
      ))
      # Warnings don't stop execution — invoke restart to continue
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(obj)) {
    message(sprintf(
      "  CreateInfercnvObject returned NULL for ref=%s — skipping",
      ref_group
    ))
    return(NULL)
  }
  
  return(obj)
}


#' Build one across-celltype inferCNV object
#'
#' Query  = all cells of query_type
#' Reference = all cells of ref_type (full group, no splitting)
#'
#' @param counts_mx       genes x cells raw count matrix
#' @param metadata        full metadata data.frame
#' @param group_col   column name for group labels
#' @param query_type      group string for the query
#' @param ref_type        group string for the reference
#' @param gene_order_file path to gene order file
#' @param chromosomes_to_exclude     chromosomes to exclude
#' @param min_max_counts  c(min, max) counts per cell filter
#'
#' @return inferCNV object
.build_across_object <- function(counts_mx,
                                 metadata,
                                 group_col,
                                 query_type,
                                 ref_type,
                                 gene_order_file,
                                 chromosomes_to_exclude    = c("MT", "Y"),
                                 min_max_counts = c(100, 1e6)) {
  
  # Subset metadata to only query + reference group
  sub_meta <- metadata[metadata[[group_col]] %in% c(as.character(query_type), as.character(ref_type)), ,
                       drop = FALSE]
  
  cells      <- sub_meta$cell_name
  sub_counts <- counts_mx[, cells, drop = FALSE]
  
  # Labels are just the group names directly
  annot <- build_annotations_df(
    cell_names   = cells,
    group_labels = sub_meta[[group_col]]
  )
  
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix         = sub_counts,
    annotations_file          = annot,
    gene_order_file           = gene_order_file,
    chr_exclude               = chromosomes_to_exclude,
    ref_group_names           = ref_type,
    min_max_counts_per_cell   = min_max_counts
  )
  
  return(obj)
}


# =============================================================================
# INTERNAL: mode-level builders
# =============================================================================

#' Build all within-celltype objects across all groups
#'
#' Returns nested list: list[[cell_group]][[ref_group]] = inferCNV object
#'
.build_all_within <- function(counts_mx,
                              metadata,
                              group_col,
                              gene_order_file,
                              chromosomes_to_exclude,
                              min_max_counts,
                              n_splits_within,
                              clonal_col,
                              donor_col) {
  
  group_clusters         <- unique(metadata[[group_col]])
  all_split_metadata <- list()
  
  message("\n── Building WITHIN objects ──────────────────────────────────────")
  message(sprintf("   n_splits_within = %d  →  refs: %s",
                  n_splits_within,
                  paste(LETTERS[seq_len(n_splits_within)],
                        collapse = ", ")))
  
  objects <- setNames(
    lapply(group_clusters, function(ct) {
      
      number_of_cells <- nrow(
        metadata[metadata[[group_col]] == ct, ])
      
      if (number_of_cells >= 100) {
        
        message("\nGroup cell: ", ct)
        
        # Split cells for this group
        split_meta <- make_splits(
          metadata      = metadata,
          group_col = group_col,
          subset_group_val = ct,
          n_splits      = n_splits_within,
          counts_mx     = counts_mx,
          clonal_col = clonal_col,
          donor_col  = donor_col
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
                chromosomes_to_exclude     = chromosomes_to_exclude,
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
        
        Filter(Negate(is.null), objs) 
        
      } else {
        
        message(sprintf(
          "\nCell Group: %s skipped — low cells (%d)",
          ct, number_of_cells))
        
        group_clusters <<- group_clusters[
          !(group_clusters == ct)]
      }
      
    }), group_clusters)
  
  # Combine split metadata
  split_metadata_combined <- do.call(
    rbind, all_split_metadata)
  rownames(split_metadata_combined) <- NULL
  
  return(list(
    objects        = objects,
    split_metadata = split_metadata_combined
  ))
}


#' Build all across-celltype objects across all cell groups
#'
#' Returns nested list: list[[query_type]][[ref_type]] = inferCNV object
#'
.build_all_across <- function(counts_mx,
                              metadata,
                              group_col,
                              gene_order_file,
                              chromosomes_to_exclude,
                              min_max_counts) {
  
  group_clusters <- unique(metadata[[group_col]])
  message("\n── Building ACROSS objects ──────────────────────────────────────")
  
  result <- setNames(lapply(group_clusters, function(query) {
    
    # Reference = every other cell groups
    ref_types <- setdiff(group_clusters, query)
    message("\nQuery: ", query, " | References: ",
            paste(ref_types, collapse = ", "))
    
    objs <- setNames(lapply(ref_types, function(ref) {
      
      message("  Building query=", query, " vs ref=", ref)
      
      tryCatch(
        .build_across_object(
          counts_mx       = counts_mx,
          metadata        = metadata,
          group_col   = group_col,
          query_type      = query,
          ref_type        = ref,
          gene_order_file = gene_order_file,
          chromosomes_to_exclude     = chromosomes_to_exclude,
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
    
  }), group_clusters)
  
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
#'                          - group_col: cell group labels
#' @param group_col   string, name of column in metadata with cell groups
#'                        (default "cell_group")
#' @param gene_order_file path to inferCNV gene order file
#'                        (e.g. hg38_gencode_v27.txt)
#' @param mode            one of "within", "across", or "both" (default "both")
#' @param chromosomes_to_exclude     chromosomes to exclude (default c("MT","Y"))
#' @param min_max_counts  c(min, max) counts per cell (default c(100, 1e6))
#'
#' @return named list with elements 'within_cell_group' and/or 'across_cell_group'
#'         each containing nested lists of inferCNV objects:
#'           within_cell_group[[cell_group]][[ref_group]]   (ref_group: A/B/C)
#'           across_cell_group[[query_type]][[ref_type]]
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
#'   cell_group = seurat_obj$cell_group
#' )
#'
#' obj_list <- make_infercnv_objects(
#'   counts_mx       = counts_mx,
#'   metadata        = metadata,
#'   group_col   = "cell_group",
#'   gene_order_file = "/path/to/hg38_gencode_v27.txt",
#'   mode            = "both"
#' )
#'
#' saveRDS(obj_list, "infercnv_objcomp.rds")
#' }
make_infercnv_objects <- function(counts_mx,
                                  metadata,
                                  group_col   = "cell_type",
                                  gene_order_file,
                                  mode            = "both",
                                  chromosomes_to_exclude     = c("MT", "Y"),
                                  min_max_counts  = c(100, 1e6),
                                  n_splits_within,
                                  clonal_col,
                                 donor_col) {
  
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
  metadata   <- validate_metadata(metadata, counts_mx, group_col)
  
  # Align counts to metadata (all cells should in the metadata after validate metadata check for this status)
  keep_cells <- intersect(colnames(counts_mx), metadata$cell_name)
  counts_mx  <- counts_mx[, keep_cells, drop = FALSE]
  metadata   <- metadata[metadata$cell_name %in% keep_cells, , drop = FALSE]
  
  # ── Build objects ────────────────────────────────────────────────────────
  result <- list()
  
  if (mode %in% c("within", "both")) {
    result$within_cell_group <- .build_all_within(
      counts_mx       = counts_mx,
      metadata        = metadata,
      group_col   = group_col,
      gene_order_file = gene_order_file,
      chromosomes_to_exclude     = chromosomes_to_exclude,
      min_max_counts  = min_max_counts,
      n_splits_within,
       clonal_col,
       donor_col
    )
  }
  
  if (mode %in% c("across", "both")) {
    result$across_cell_group <- .build_all_across(
      counts_mx       = counts_mx,
      metadata        = metadata,
      group_col   = group_col,
      gene_order_file = gene_order_file,
      chromosomes_to_exclude     = chromosomes_to_exclude,
      min_max_counts  = min_max_counts
    )
  }
  
  # ── Summary ──────────────────────────────────────────────────────────────
  message("\n── Summary ───────────────────────────────────────────────────────")
  
  if (!is.null(result$within_cell_group)) {
    total_within <- sum(sapply(result$within_cell_group$objects, length))
    message("Within objects built: ", total_within,
            " (", length(result$within_cell_group$objects), " cell groups x ",length(result$within_cell_group$objects[[1]]),  " refs)")
  }
  
  if (!is.null(result$across_cell_group)) {
    total_across <- sum(sapply(result$across_cell_group, length))
    message("Across objects built: ", total_across)
  }
  
  message("\nDone. Save with: saveRDS(obj_list, 'infercnv_objcomp.rds')")
  
  return(result)
}




# =============================================================================
# infercnv_run.R
# Function to run inferCNV on all objects produced by make_infercnv_objects()
# =============================================================================


#' Run inferCNV on all objects in a nested list
#'
#' Expects the nested list structure produced by make_infercnv_objects():
#'
#' Output directories are created automatically:

#'   {base_outdir}/across/{query_type}/ref_{ref_type}/
#'
#' @param infercnv_obj_list nested list from make_infercnv_objects()
#' @param base_outdir       root output directory (will be created if needed)
#' @param cutoff            minimum average counts per gene for reference cells
#'                          (default 0.1 — correct for RNA-seq, use 1 for 10x)
#' @param cluster_by_groups logical, cluster cells within groups (default TRUE)
#' @param analysis_mode     one of "subclusters", "samples", "cells"
#'                          (default "subclusters")
#' @param window_length     smoothing window length (default 140)
#' @param plot_steps        logical, save intermediate plots (default FALSE)
#' @param no_plot           logical, skip final heatmap (default TRUE —
#'                          set FALSE if you want plots)
#' @param no_prelim_plot    logical, skip preliminary plot (default TRUE)
#' @param plot_probabilities logical (default FALSE)
#' @param diagnostics       logical (default FALSE)
#' @param inspect_subclusters logical (default FALSE)
#' @param resume_if_exists  logical, skip runs where out_dir already has
#'                          run.final.infercnv_obj (default TRUE)
#'
#' @return invisibly returns a data.frame log of all runs with status
#'         (success / failed / skipped)
#'
#' @examples
#' \dontrun{
#' obj_list <- readRDS("infercnv_objcomp.rds")
#'
#' run_log <- run_infercnv_objects(
#'   infercnv_obj_list = obj_list,
#'   base_outdir       = "/path/to/output/",
#'   no_plot           = FALSE   # set TRUE on HPC to skip heavy plotting
#' )
#'
#' # Check what failed
#' run_log[run_log$status == "failed", ]
#' }
run_infercnv_objects <- function(infercnv_obj_list,
                                 base_outdir,
                                 cutoff               = 0.1,
                                 cluster_by_groups    = TRUE,
                                 analysis_mode        = "subclusters",
                                 window_length        = 140,
                                 plot_steps           = FALSE,
                                 no_plot              = TRUE,
                                 no_prelim_plot       = TRUE,
                                 plot_probabilities   = FALSE,
                                 diagnostics          = FALSE,
                                 inspect_subclusters  = FALSE,
                                 resume_if_exists     = TRUE) {
  
  # ── Input checks ────────────────────────────────────────────────────────
  valid_modes <- c("within_cell_group", "across_cell_group")
  found_modes <- intersect(names(infercnv_obj_list), valid_modes)
  
  if (length(found_modes) == 0) {
    stop(
      "infercnv_obj_list does not contain expected modes. ",
      "Expected names: 'within_cell_group' and/or 'across_cell_group'. ",
      "Got: ", paste(names(infercnv_obj_list), collapse = ", ")
    )
  }
  
  dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Run log ──────────────────────────────────────────────────────────────
  run_log <- data.frame(
    mode      = character(),
    cell_group = character(),
    comp      = character(),
    out_dir   = character(),
    status    = character(),
    message   = character(),
    stringsAsFactors = FALSE
  )
  
  # ── Map mode names to output folder names ────────────────────────────────
  mode_folder_map <- c(
    within_cell_group = "within",
    across_cell_group = "across"
  )
  processed_objects <- list()
  # ── Loop ─────────────────────────────────────────────────────────────────
  for (mode in found_modes) {
    
    mode_folder  <- mode_folder_map[[mode]]
    mode_objects <- infercnv_obj_list[[mode]]
    processed_objects[[mode]] <- list()
    
    for (cell_group in names(mode_objects)) {
      
      type_objects <- mode_objects[[cell_group]]
      
      for (comp in names(type_objects)) {
        
        infer_obj <- type_objects[[comp]]
        
        # NULL guard — object may have failed during creation
        if (is.null(infer_obj)) {
          message("Skipping NULL object: ", mode, " / ", cell_group, " / ", comp)
          run_log <- rbind(run_log, data.frame(
            mode = mode, cell_group = cell_group, comp = comp,
            out_dir = NA, status = "skipped_null", message = "Object is NULL",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        # Build output directory
        outdir <- file.path(base_outdir, mode_folder, cell_group, comp)
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        
        # Resume check — skip if already completed
        final_obj_path <- file.path(outdir, "run.final.infercnv_obj")
        if (resume_if_exists && file.exists(final_obj_path)) {
          message("Skipping (already complete): ",
                  mode_folder, "/", cell_group, "/", comp)
          processed_objects[[mode]][[cell_group]][[comp]] <- readRDS(final_obj_path)
           
          run_log <- rbind(run_log, data.frame(
            mode = mode, cell_group = cell_group, comp = comp,
            out_dir = outdir, status = "skipped_exists",
            message = "run.final.infercnv_obj already present",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        message("\n── Running: ", mode_folder, " / ", cell_group, " / ", comp)
        message("   Output: ", outdir)
        
        status  <- "success"
        err_msg <- ""
    
        tryCatch({
          options(scipen = 100)
          processed_obj = infercnv::run(
            infercnv_obj          = infer_obj,
            out_dir               = outdir,
            cutoff                = cutoff,
            cluster_by_groups     = cluster_by_groups,
            HMM                   = F,
            denoise               = F,
            analysis_mode         = analysis_mode,
            output_format         = NA,
            no_plot               = no_plot,
            #remove_genes_at_chr_ends = T,
            no_prelim_plot        = no_prelim_plot,
            window_length         = window_length,
            plot_probabilities    = plot_probabilities,
            plot_steps            = plot_steps,
            diagnostics           = diagnostics,
            inspect_subclusters   = inspect_subclusters
          )
          
        processed_objects[[mode]][[cell_group]][[comp]] <- processed_obj
        
        message("   Done: ", mode_folder, " / ", cell_group, " / ", comp)
          
        }, error = function(e) {
          status  <<- "failed"
          err_msg <<- conditionMessage(e)
          warning("FAILED: ", mode_folder, "/", cell_group, "/", comp,
                  "\n  Error: ", err_msg)
        })
        
        
        run_log <- rbind(run_log, data.frame(
          mode      = mode,
          cell_group = cell_group,
          comp      = comp,
          out_dir   = outdir,
          status    = status,
          message   = err_msg,
          stringsAsFactors = FALSE
        ))
        
        # Clean up memory between runs
        rm(infer_obj)
        gc()
        
      } # comp loop
    } # cell_group loop
  } # mode loop
  
  # ── Final summary ────────────────────────────────────────────────────────
  message("\n── Run Summary ───────────────────────────────────────────────────")
  print(table(run_log$status))
  
  failed <- run_log[run_log$status == "failed", ]
  if (nrow(failed) > 0) {
    message("\nFailed runs:")
    print(failed[, c("mode", "comp", "message")])
  }
  
  # Save log to base_outdir
  log_path <- file.path(base_outdir, "run_log.tsv")
  write.table(run_log, log_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("\nRun log saved to: ", log_path)
  
  list(
  run_log           = run_log,
  processed_objects = processed_objects)
}



run_infercnv_pipeline <- function(
    
  # ---- make_infercnv_objects parameters ----------------------------------
  counts_mx,
  metadata,
  group_col   = "cell_type",
  gene_order_file,
  mode            = c("within", "across"),
  chromosomes_to_exclude     = c("MT", "Y"),
  min_max_counts  = c(100, 1e6),
  n_splits_within = 3,
  clonal_col = NULL,
  donor_col = NULL,
  # ---- run_infercnv_objects parameters -----------------------------------
  base_outdir,
  cutoff          = 0.1,
  cluster_by_groups = TRUE,
  analysis_mode   = "subclusters",
  window_length   = 140,
  no_plot         = TRUE,
  resume_if_exists = TRUE
  
) {
  
  mode <- match.arg(mode)
  
  # ---- Input validation ---------------------------------------------------
  if (!is.matrix(counts_mx) && !inherits(counts_mx, "dgCMatrix")) {
    stop("counts_mx must be a matrix or sparse matrix (dgCMatrix).")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame.")
  }
  if (!group_col %in% colnames(metadata)) {
    stop("group_col '", group_col, "' not found in metadata.")
  }
  if (!file.exists(gene_order_file)) {
    stop("gene_order_file not found: ", gene_order_file)
  }
  if (!dir.exists(base_outdir)) {
    message("base_outdir does not exist — creating: ", base_outdir)
    dir.create(base_outdir, recursive = TRUE)
  }
  if (length(min_max_counts) != 2L) {
    stop("min_max_counts must be a numeric vector of length 2: c(min, max).")
  }
  
  message(sprintf(paste0(
    "inferCNV pipeline starting:\n",
    "  Mode:           %s\n",
    "  Cells:          %d\n",
    "  Genes:          %d\n",
    "  Cell Groups:     %s\n",
    "  Output dir:     %s"
  ),
  mode,
  ncol(counts_mx),
  nrow(counts_mx),
  paste(unique(metadata[[group_col]]), collapse = ", "),
  base_outdir
  ))
  
  t_start <- proc.time()
  
  # ---- Step 1: make inferCNV objects --------------------------------------
  message("\n[1/2] Creating inferCNV objects...")
  t_make_start <- proc.time()
  
  obj_list <- make_infercnv_objects(
    counts_mx       = counts_mx,
    metadata        = metadata,
    group_col   = group_col,
    gene_order_file = gene_order_file,
    mode            = mode,
    chromosomes_to_exclude     = chromosomes_to_exclude,
    min_max_counts  = min_max_counts,
    n_splits_within = n_splits_within,
    clonal_col,
  donor_col
  )
  
  t_make_end <- proc.time()
  
  # ---- Lightweight sanity check between steps -----------------------------
  if (is.null(obj_list)) {
    stop("make_infercnv_objects() returned NULL — check inputs.")
  }
  if (!"objects" %in% names(obj_list[["within_cell_group"]]) &&
      mode %in% c("within")) {
    stop(
      "Expected 'objects' element in obj_list$within_cell_group. ",
      "make_infercnv_objects() may have failed silently."
    )
  }
  
  message(sprintf(
    "Objects created in %.1f seconds.",
    (t_make_end - t_make_start)[["elapsed"]]
  ))
  
  split_metadata <- NULL  # default — only populated for within/both
  
  if (mode %in% c("within")) {
    
    if (!"split_metadata" %in% names(obj_list[["within_cell_group"]])) {
      stop(
        "Expected 'split_metadata' in obj_list$within_cell_group. ",
        "Check make_infercnv_objects() output structure."
      )
    }
    
    split_metadata <- obj_list[["within_cell_group"]][["split_metadata"]]
    obj_list[["within_cell_group"]] <- obj_list[["within_cell_group"]][["objects"]]
    
    message(sprintf(
      "Within mode: split_metadata extracted (%d cell group splits)",
      length(split_metadata)
    ))
  }
  
  
  # ---- Step 2: run inferCNV -----------------------------------------------
  t_run_start <- proc.time()
  
  run_result  <- run_infercnv_objects(
    infercnv_obj_list = obj_list,
    base_outdir       = base_outdir,
    cutoff            = cutoff,
    cluster_by_groups = cluster_by_groups,
    analysis_mode     = analysis_mode,
    window_length     = window_length,
    no_plot           = no_plot,
    resume_if_exists  = resume_if_exists
  )
  
 t_run_end <- proc.time()
  
  run_log           <- run_result$run_log
  processed_objects <- run_result$processed_objects
  
  # ---- Lightweight sanity check on run_log --------------------------------
  if (is.null(run_log)) {
    warning("run_infercnv_objects() returned NULL run_log — check output directory.")
  }
  
  t_end     <- proc.time()
  runtime   <- (t_end - t_start)[["elapsed"]]
  
  message(sprintf(paste0(
    "\nPipeline complete:\n",
    "  Make time:   %.1f seconds\n",
    "  Run time:    %.1f seconds\n",
    "  Total time:  %.1f seconds (%.1f minutes)"
  ),
  (t_make_end - t_make_start)[["elapsed"]],
  (t_run_end  - t_run_start)[["elapsed"]],
  runtime,
  runtime / 60
  ))
  
  
  list(
    obj_list       = processed_objects,
    run_log        = run_log,
    metadata       = split_metadata,
    runtime        = runtime
  )
}



#' Convert a wide expression matrix to long format
#'
#' Reshapes a gene-by-cell expression table into a long-format data frame with
#' one row per gene-cell combination.
#'
#' @param expr_df A data frame or matrix with genes as rows and cells as columns.
melt_expr_to_long <- function(expr_df, cell_prefix = "cell_") {
  
  if (is.null(rownames(expr_df)) || 
      all(rownames(expr_df) == as.character(seq_len(nrow(expr_df))))) {
    stop("expr_df has no meaningful rownames — gene names expected as rownames.")
  }
  
  expr_df |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(
      cols      = -gene,
      names_to  = "cell_name",
      values_to = "state_raw"
    ) |>
    # Strip prefix from cell_name after melting if present — silently skips
    # cells that don't carry the prefix, so mixed naming is handled safely
    dplyr::mutate(
      cell_name = stringr::str_remove(cell_name, paste0("^", cell_prefix))
    )
}




#' Attach genomic coordinates to long-format gene data
#'
#' Merges long-format gene-level data with a gene order table containing genomic
#' coordinates.
#'
#' @param long_df A long-format data frame containing a \code{gene} column.
#' @param gene_order A data frame containing gene coordinates. Row names are
#'   assumed to store gene names.
#'
#' @return A merged data frame containing gene-level values and genomic
#'   coordinate columns.
attach_gene_order <- function(long_df, gene_order) {
  
  gene_order["gene"] <- rownames(gene_order)
  
  required <- c("gene", "chr", "start", "stop")
  if (!all(required %in% colnames(gene_order))) {
    stop("gene_order missing required columns")
  }
  
  merged <- long_df |>
    inner_join(by = "gene",
               y = gene_order)
  
  if (nrow(merged) == 0) {
    stop("Merge produced empty table — gene names may not match between expr and gene_order.")
  }
  if (nrow(merged) < nrow(long_df)) {
    warning(sprintf(
      "Inner join dropped %d rows — some genes have no coordinate annotation.",
      nrow(long_df) - nrow(merged)
    ))
  }
  
  merged
}




#' Discretize raw CNV scores into gain, loss, or neutral states
#' 
#' Converts a continuous inferCNV score into a categorical CNV state using
#' mean +/- k * SD thresholds computed globally across the input data frame.
#' 
#' 
#' Global thresholding is appropriate here because inferCNV output values are
#' already centred by the tool's internal smoothing relative to the reference
#' signal. The neutral state dominates the distribution and anchors the global
#' mean near the neutral baseline, so gains and losses appear as symmetric
#' tail deviations. Per-cell thresholding would artificially force every cell
#' to have calls even in largely euploid cells.
#'
#' @param df A data frame containing a numeric \code{state_raw} column.
#' @param k Numeric multiplier for the standard deviation cutoff. Default is 1.5.
discretize_cnv_state_infer_cnv <- function(df, k =1.5) {
  if (!"state_raw" %in% colnames(df)) {
    stop("Missing required column: state_raw")
  }
  if (!is.numeric(df$state_raw)) {
    stop("state_raw must be numeric")
  }
  if (all(is.na(df$state_raw))) {
    stop("state_raw is entirely NA — check upstream expression extraction")
  }
  
  mu <- mean(df$state_raw, na.rm = TRUE)
  sigma <- sd(df$state_raw, na.rm = TRUE)
  
  upper <- mu + k * sigma
  lower <- mu - k * sigma
  df |>
    mutate(
      state = case_when(
        state_raw < lower ~ "loss",
        state_raw > upper ~ "gain",
        TRUE ~ "neutral"
      )
    )
}




#' Load and prepare inferCNV reference outputs
#'
#' Converts inferCNV outputs into a unified long-format table with genomic
#' coordinates, discrete CNV states, and reference labels.
#'
#' @param infercnv_list A named list of inferCNV objects. Each object is loaded and the 1 and 2 are taken
#'
#' @return A tibble containing gene-level CNV calls across all references.
load_and_prepare_infercnv_reference <- function(infercnv_list, k = 1.5) {
  refs <- names(infercnv_list)
  
  res <- lapply(seq_along(infercnv_list), function(i) {
    
    infercnv_obj <- infercnv_list[[i]]
    reference <- refs[i]
    melt_expr_to_long(as.data.frame(infercnv_obj[[1]])) |>
      attach_gene_order(infercnv_obj[[2]]) |>
      discretize_cnv_state_infer_cnv(k = k) |>
      dplyr::mutate(reference = reference)
     
    
  })
  
  dplyr::bind_rows(res)
}



#' Discover and load inferCNV run outputs
#'
#' Finds inferCNV result files matching a pattern inside one or more reference
#' directories and loads them as R objects with \code{readRDS()}.
#'
#' @param base_dir Optional base directory containing the reference-specific
#'   subdirectories.
#' @param ref_dirs Character vector of reference directory names or paths.
#' @param pattern Regular expression used to match inferCNV result files.
#' 
#' @return A named list of loaded inferCNV objects, one per reference directory.
#'
#' @details
#' For each entry in \code{ref_dirs}, the function searches for files matching
#' \code{pattern}. Matching files are loaded with \code{readRDS()} and stored in
#' a list using the reference name as the list element name.
discover_infercnv_runs <- function(base_dir = NULL,
                                   ref_dirs,
                                   pattern = "^run\\.final") {
  ref_paths <- if (!is.null(base_dir)) file.path(base_dir, ref_dirs) else ref_dirs
  names(ref_paths) <- ref_dirs
  
  runs <- lapply(ref_paths, function(ref_path){
    files <- list.files(ref_path, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) stop(sprintf("No files matching '%s' found in %s", pattern, ref_path))
    if (length(files)  > 1) stop(sprintf("Multiple files matching '%s' found in %s", pattern, ref_path))
    inferobj <- readRDS(files[[1]])
  })
  
}


load_infercnv_data <- function(base_dir, ref_dirs, pattern = "^run\\.final", k = 1.5) {
  
  # Discover and load runs — already written as discover_infercnv_runs()
  infer_objs <- discover_infercnv_runs(
    base_dir = base_dir,
    ref_dirs = ref_dirs,
    pattern  = pattern
  )
  
  # Convert to standard schema — already written
  infer_objs_1 <- lapply(infer_objs, function(x) {
    list(x@expr.data, x@gene_order)
  })
  
  return(list(gene_level_df = load_and_prepare_infercnv_reference(infer_objs_1, k = k), gene_order_df = infer_objs_1[[1]][[2]]))
  
}





