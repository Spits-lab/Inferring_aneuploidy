#' @title cnv_processing
#'
#' @description
#' 
#' Functions that perform further processing from the CNV performing a final filtering and Confidence Score in a Single Tool Approach
#' 
#' @author Pedro Granjo
#' @date 02-09-2026
#'


# Package groups
cran_packages <- c(
  "dplyr", "tidyr", "data.table", "cowplot", "igraph", "BiocManager", "purrr"
)


bioc_packages <- c(
  "GenomicRanges", "IRanges"
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


#' Add mode and cell type metadata to pipeline output tables
#'
#' Iterates over a nested results list (mode → cell_type → tables)
#' and adds mode and cell_type columns to each output dataframe.
#' Uses imap to access names at each level — two levels deep matching
#' the results list structure.
#'
#' @param results Nested list output of run_single_tool_pipeline.
#' @param mode_col Name of the mode column to add. Default "mode".
#'
#' @return Same nested list structure with metadata columns added to
#'   all dataframe elements.
add_metadata <- function(
    df,
    mode_name      = "within",
    metadata
) {
  
  difference <- setdiff(df$cell_name, metadata$cell_name)
  
  if(length(difference) > 0){
    stop("you should not have empty cell type values. It will disrupt the clustering analysis later on")
  }
  df |> dplyr::mutate(
            mode = mode_name
            ) |> left_join(metadata, by = "cell_name")
      }

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


#' Collapse consecutive gene-level CNV calls into segments
#'
#' Merges adjacent non-neutral gene-level CNV calls into larger genomic segments
#' within each reference, cell, and chromosome.
#'
#' @param gene_cnv_df A data frame of gene-level CNV calls.
#'
#' @return A tibble of collapsed CNV segments.
collapse_genes_to_cnv_segments <- function(gene_cnv_df) {

  required <- c("reference","cell_name","chr","start","stop","state")
  if (!all(required %in% colnames(gene_cnv_df))) {
    stop("Missing required columns")
  }
  
  dt <- as.data.table(gene_cnv_df)
  setorder(dt, reference, cell_name, chr, start)

  
  #number of rows before
  n_before <- nrow(dt)
  # Identify state breaks
  dt[, prev_state := data.table::shift(state), by = .(reference, cell_name, chr)]
  dt[, new_block := 
       state == "neutral" |
       prev_state == "neutral" |
       state != prev_state |
       is.na(prev_state)
  ]
  
  dt[, block_id := cumsum(new_block), by = .(reference, cell_name, chr)]

  segments <- dt[state != "neutral",
                 .(
                   start = min(start),
                   stop  = max(stop),
                   state = (state),
                   raw_state = mean(state_raw),
                   n_genes   = .N,
                   genes     = paste(gene,   # ← gene names
                   collapse = ";"),
                   genes_total_length = sum(
                     stop - start)
                   ),
                 by = .(reference, cell_name, chr, block_id)
  ]
  
  segments[, block_id := NULL]
  
  segments <- unique(segments)
  
  n_after <- nrow(segments)
  
  if (n_after >= n_before) {
    warning(sprintf(
      "Number of rows did not decrease after collapsing genes: before = %d, after = %d",
      n_before, n_after
    ))
  } else {
    message(sprintf(
      "Collapsed genes successfully: before = %d rows, after = %d rows",
      n_before, n_after
    ))
  }
  
  return(as_tibble(segments))
}


filt_remove_refs_cells <- function(df, metadata, filter_seq_mb, mode,
                                   remove_ref = TRUE) {
  
  # ---- Pre-filter ---------------------------------------------------------
  n_input <- nrow(df)
  
  df <- df |>
    dplyr::arrange(reference, cell_name, chr, cnv_state, start) |>
    dplyr::mutate(
      cnv_length    = as.numeric(end) - as.numeric(start) + 1,
      cnv_length_mb = cnv_length / 1e6
    ) |>
    dplyr::filter(cnv_length_mb > filter_seq_mb)
  
  n_postfilter <- nrow(df)
  
  message(sprintf(paste0(
    "Length filter (> %.0f Mb):\n",
    "  Input segments:    %d\n",
    "  Retained:          %d\n",
    "  Removed:           %d (%.1f%%)"
  ),
  filter_seq_mb,
  n_input,
  n_postfilter,
  n_input - n_postfilter,
  100 * (n_input - n_postfilter) / n_input
  ))
  
  # ---- Join metadata ------------------------------------------------------
  df_joined <- df |>
    dplyr::left_join(metadata, by = "cell_name")
  
  unmatched <- sum(is.na(df_joined$reference))
  if (unmatched > 0L) {
    warning(sprintf(
      "%d cell_name(s) did not match metadata — these cells have no group assignment.\n  Unmatched: %s",
      unmatched,
      paste(df$cell_name[is.na(df_joined$reference)], collapse = ", ")
    ))
  }
  
  # ---- Remove reference cells ---------------------------------------------
  # In within mode: cells from the same group as the reference are removed
  #   — comparing a group against itself inflates similarity
  # In across mode: cells whose cell type IS the reference are removed
  #   — a cell type cannot serve as its own reference
  if (remove_ref) {
    
    n_before_ref_removal <- nrow(df_joined)
  
    if (mode == "within") {
    
    df_joined <- df_joined |>
      dplyr::filter(!(reference == split_group)) |>
      dplyr::select(-dplyr::any_of("split_group"))
    
    message(sprintf(paste0(
      "Reference removal (mode = within):\n",
      "  Removing cells belonging to the same group as their reference\n",
      "  Before: %d segments\n",
      "  After:  %d segments\n",
      "  Removed: %d segments (self-reference cells)"
    ),
    n_before_ref_removal,
    nrow(df_joined),
    n_before_ref_removal - nrow(df_joined)
    ))
    
    } else if (mode == "across") {
    
    df_joined <- df_joined |>
      dplyr::filter(!(reference == cell_type)) |>
      dplyr::select(-dplyr::any_of("split_group"))
    
    message(sprintf(paste0(
      "Reference removal (mode = across):\n",
      "  Removing cells where their cell type is the reference\n",
      "  — a cell type cannot be compared against itself\n",
      "  Before: %d segments\n",
      "  After:  %d segments\n",
      "  Removed: %d segments (same cell type as reference)"
    ),
    n_before_ref_removal,
    nrow(df_joined),
    n_before_ref_removal - nrow(df_joined)
    ))
    }
  
  # ---- Final check --------------------------------------------------------
  if (nrow(df_joined) == 0L) {
    stop(sprintf(
      "No segments remaining after reference removal in mode = '%s'.\n",
      "  Check that reference labels in df match group labels in metadata.",
      mode
    ))
  }
  
  message(sprintf(paste0(
    "Pipeline ready:\n",
    "Total of %d segments retained\n",
    "In %d cells\n",
    "mode = '%s'"),
    nrow(df_joined),
    dplyr::n_distinct(df_joined$cell_name),
    mode
  ))}
  
  return(df_joined)
}

#' Merge nearby CNV segments
#'
#' Merges CNV segments that are close together within the same reference, cell,
#' chromosome, and CNV state.
#'
#' @param df A data frame of CNV segments.
#' @param max_gap Maximum genomic gap allowed between two segments for merging.
#'   Default is 100000.
#'
#' @return A data frame of merged CNV regions.
merge_nearby_regions <- function(df, max_gap = 100000L) {
  # ---- Input validation ---------------------------------------------------
  required <- c("reference", "cell_name", "chr", "state", "start", "stop")
  if (!all(required %in% colnames(df))) {
    stop("Missing required columns")
  }
  
  if (any(df$start >= df$stop)) {
    stop("Invalid CNV intervals: start >= stop")
  }
  
  if (any(is.na(df[, required]))) {
    stop("Null values detected in CNV table")
  }
  

  n_postfilter <- nrow(df)
  # ---- Core merging logic -------------------------------------------------
  merged_df <- df %>%
    group_by(reference, cell_name, chr, state) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(
      gap = start - lag(stop),
      new_block = is.na(gap) | gap > max_gap,
      merge_id = cumsum(new_block)
    ) %>%
    group_by(reference, cell_name, chr, state, merge_id) %>%
    summarise(
      start      = min(start),
      end        = max(stop),
      n_segments = n(),
      raw_state  = mean(raw_state),
      .groups    = "drop"
    ) %>%
    dplyr::rename(cnv_state = state) %>%
    dplyr::select(-merge_id) %>%
    dplyr::arrange(reference, cell_name, chr, cnv_state, start)
  

  # ---- Sanity checks ------------------------------------------------------
  n_merged <- nrow(merged_df)
  # Check 1: merge summary reporting
  message(sprintf(paste0(
    "Merge summary:\n",
    "  Input:                    %d rows\n",
    "  After merging:            %d rows\n",
    "  Reduction:                %d rows (%.1f%%)"
  ),
  n_postfilter,
  n_merged,
  n_postfilter - n_merged,
  100 * (n_postfilter - n_merged) / n_postfilter
  ))
  
  # Check 2: warn if merging produced no reduction
  if (n_merged >= n_postfilter) {
    warning(sprintf(
      "No reduction after merging: before = %d, after = %d. using the max_gap = %d.",
      as.integer(n_postfilter),
      as.integer(n_merged),
      as.integer(max_gap)
    ))
  }
  
  # Check 3: empty output guard
  if (n_merged == 0L) {
    stop("No segments remain after merging. Check max_gap and filter_seq_mb.")
  }
  
  # Check 4: coordinate integrity after merging
  if (any(merged_df$start >= merged_df$end, na.rm = TRUE)) {
    stop("Merged segments have start >= end. Check summarise logic.")
  }
  
  # Check 5: no NA in key output columns
  key_cols  <- c("reference", "cell_name", "chr", "cnv_state", "start", "end")
  na_counts <- colSums(is.na(merged_df[, key_cols]))
  if (any(na_counts > 0L)) {
    stop(
      "NA values in output columns after merging: ",
      paste(names(na_counts[na_counts > 0L]), collapse = ", ")
    )
  }
  
  # Check 6: n_segments should always be >= 1
  if (any(merged_df$n_segments < 1L, na.rm = TRUE)) {
    stop("Merged segments with n_segments < 1 detected. Check summarise logic.")
  }
  
  
  merged_df
}



#' Compute pairwise overlap scores using a named strategy
#'
#' Acts as the single entry point for all overlap methods. Individual strategies
#' are defined as internal nested functions — adding a new strategy means adding
#' a nested function here and registering it in the registry, nothing else changes
#' in the pipeline.
#'
#' @param q_start,q_end Integer vectors. Query segment coordinates.
#' @param s_start,s_end Integer vectors. Subject segment coordinates.
#' @param method Character string naming the overlap strategy to use.
#'   One of "reciprocal", "jaccard", "symmetric_reciprocal".
#'
#' @return Numeric vector of overlap scores in [0, 1], same length as inputs.
compute_overlap <- function(q_start, q_end, s_start, s_end, method = "reciprocal") {
  
  # --- Internal strategy definitions ---------------------------------------
  # Each function shares the same signature and returns a numeric vector in [0,1]
  # Add new strategies here as additional nested functions + registry entry
  
  .reciprocal <- function(q_start, q_end, s_start, s_end) {
    intersection_len <- pmax(0L, pmin(q_end, s_end) - pmax(q_start, s_start) + 1L)
    q_len            <- q_end - q_start + 1L
    s_len            <- s_end - s_start + 1L
    intersection_len / pmax(q_len, s_len)
  }
  
  .jaccard <- function(q_start, q_end, s_start, s_end) {
    intersection_len <- pmax(0L, pmin(q_end, s_end) - pmax(q_start, s_start) + 1L)
    q_len            <- q_end - q_start + 1L
    s_len            <- s_end - s_start + 1L
    union_len        <- q_len + s_len - intersection_len
    intersection_len / union_len
  }
  
  .symmetric_reciprocal <- function(q_start, q_end, s_start, s_end) {
    intersection_len <- pmax(0L, pmin(q_end, s_end) - pmax(q_start, s_start) + 1L)
    q_len            <- q_end - q_start + 1L
    s_len            <- s_end - s_start + 1L
    (intersection_len / q_len + intersection_len / s_len) / 2
  }
  
  # --- Internal registry ---------------------------------------------------
  # Maps method name strings to their corresponding internal functions
  
  .registry <- list(
    reciprocal           = .reciprocal,
    jaccard              = .jaccard,
    symmetric_reciprocal = .symmetric_reciprocal
  )
  
  # --- Validate and dispatch -----------------------------------------------
  
  valid_methods <- names(.registry)
  if (!method %in% valid_methods) {
    stop(
      "Unknown overlap method: '", method, "'. ",
      "Valid options are: ", paste(valid_methods, collapse = ", ")
    )
  }
  
  .registry[[method]](q_start, q_end, s_start, s_end)
}





find_maximal_cliches <- function(q_pass, s_pass, grp){
  n   <- nrow(grp)
  # ---- Maximal clique assignment ------------------------------------------
  # Maximal cliques guarantee every pair within a group directly passes
  # the overlap threshold — no transitive chaining
  total_duplicated <- 0L
  connected_idx <- unique(c(q_pass, s_pass))
  isolated_idx  <- setdiff(seq_len(n), connected_idx)
  
  g <- igraph::graph_from_edgelist(
    cbind(as.character(q_pass), as.character(s_pass)),
    directed = FALSE
  )
  
  cliques <- igraph::max_cliques(g)
  
  # Build a mapping: local segment index → vector of clique IDs it belongs to
  # A segment in one clique gets one ID
  # A segment in multiple cliques gets one row per clique — duplicated
  segment_clique_map <- vector("list", n)
  
  for (k in seq_along(cliques)) {
    members <- as.integer(names(cliques[[k]]))
    for (m in members) {
      segment_clique_map[[m]] <- c(segment_clique_map[[m]], k)
    }
  }
  
  # Count duplicated segments — those appearing in more than one clique
  n_duplicated <- sum(vapply(
    segment_clique_map[connected_idx],
    function(x) max(0L, length(x) - 1L),
    integer(1)
  ))
  
  
  total_duplicated <- total_duplicated + n_duplicated
  
  # Build output rows for this group
  # Each segment is emitted once per clique it belongs to
  # Isolated segments (no passing overlap partner) are dropped — NA equiv_id
  group_output <- vector("list", length(connected_idx))
  unresolved_rows <- vector("list", length(connected_idx)) 
  for (j in seq_along(connected_idx)) {
    seg_idx    <- connected_idx[j]
    clique_ids <- segment_clique_map[[seg_idx]]
    
    if (is.null(clique_ids) || length(clique_ids) == 0L) {
      # Connected but assigned to no clique — should not happen, guard only
      unresolved_row <- grp[seg_idx, ]
      unresolved_row$removal_reason <- "unresolved_clique"
      unresolved_rows[[j]] <- unresolved_row
      next
    }
    
    # One row per clique membership, with globally unique equiv ID
    seg_rows <- grp[rep(seg_idx, length(clique_ids)), ]
    seg_rows$local_clique_id <- clique_ids
    group_output[[j]] <- seg_rows
  }
  
  isolated_rows <- grp[isolated_idx, ]
  isolated_rows$removal_reason <- "non_overlap (isolated)"
  
  if (!is.null(unresolved_rows)){
    unresolved_rows <- do.call(rbind, unresolved_rows)
    removed_log <- bind_rows(unresolved_rows,isolated_rows)
  }else{
    removed_log <- bind_rows(isolated_rows)
  }
  
  if (length(connected_idx) + n_duplicated != nrow(dplyr::bind_rows(group_output))) {
    rows_groups<-  nrow(dplyr::bind_rows(group_output))
    stop(sprintf(
      "Somethings is WRONG. Check this numbers they should add up!  
       connected=%d, isolated=%d, n_duplicated=%d, group_output_rows=%d",
      length(connected_idx),
      length(isolated_idx),
      n_duplicated,
      rows_groups
    ))
  }
  
  return(list(rows =  dplyr::bind_rows(group_output), 
              n_duplicated = total_duplicated, removed = removed_log))
  }




process_cnv_cluster <- function(grp,overlap_method,min_overlap){
  
  n           <- nrow(grp)

  # Single-segment group — trivially its own equivalence class
  if (n == 1L) {
    removed_log <- grp
    removed_log$removal_reason <- "single_sequence"
    grp$local_clique_id <- seq_len(n)
    return(list(rows = grp[0, ], n_duplicated = 0L, removed = removed_log))
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = grp$chr,
    ranges   = IRanges::IRanges(start = grp$start, end = grp$end),
    strand   = "*"
  )
  names(gr) <- as.character(seq_len(n))
  
  hits <- GenomicRanges::findOverlaps(gr, gr, type = "any", select = "all")
  hits <- hits[S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits)]
  
  if (length(hits) == 0L) {
    removed_log <- grp
    removed_log$removal_reason <- "0 hits"
    grp$local_clique_id <- seq_len(n)
    return(list(rows = grp[0,], n_duplicated = 0L, removed = removed_log))
  }
  
  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)

  scores <- compute_overlap(
    q_start = grp$start[q_idx],
    q_end   = grp$end[q_idx],
    s_start = grp$start[s_idx],
    s_end   = grp$end[s_idx],
    method  = overlap_method
  )
  
  passing <- scores >= min_overlap
  q_pass  <- q_idx[passing]
  s_pass  <- s_idx[passing]
  
  
  if (length(q_pass) == 0L) {
    # Overlaps exist but none pass threshold — all separate classes
    removed_log <- grp
    removed_log$removal_reason <- "non overall (threshold based)"
    grp$local_clique_id <- seq_len(n)
    return(list(rows = grp[0,], n_duplicated = 0L, removed = removed_log))
  }

  res <- find_maximal_cliches(q_pass, s_pass, grp)
  
  
  return(res)
  
  
  return(res)
  
}




assign_cnv_equivalence <- function(
    df,
    min_overlap = 0.5,
    overlap_method         = "reciprocal",
    filter_seq_mb          = 7,
    parallel               = FALSE,
    by_columns = c("cell_name", "chr", "cnv_state"),
    n_cores = 1L
) {
  
  # ---- Input validation ---------------------------------------------------
  missing_cols  <- setdiff(c(by_columns, "start", "end"), colnames(df))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (any(df$start > df$end, na.rm = TRUE)) {
    stop("Detected rows where start > end. Check upstream segmentation.")
  }
  

  cnv_missing_collumns <- setdiff(c("cnv_length_mb","cnv_length"), colnames(df))
  if(length(cnv_missing_collumns) > 0L){
    df <- df |>
      dplyr::mutate(
        cnv_length    = end - start + 1L,
        cnv_length_mb = cnv_length / 1e6
      )
  }
  
  # ---- Pre-filter ---------------------------------------------------------
  if(filter_seq_mb > 0L) {
    n_rows_input <- nrow(df)
    df <- df |>
      dplyr::filter(cnv_length_mb > filter_seq_mb)
    
    n_rows_postfilter <- nrow(df)
    n_rows_removed    <- n_rows_input - n_rows_postfilter
    
    message(sprintf(paste0(
      "Pre-filter:\n", 
      "%d input rows\n",
      "%d retained\n", 
      "%.1f%% removed\n",
      "Segments length > %.0f Mb"),
      n_rows_input,
      n_rows_postfilter,
      100 * n_rows_removed / n_rows_input,
      filter_seq_mb
    ))
    
    if (n_rows_postfilter == 0L) {
      stop(sprintf(
        "No segments remain after length filter (> %.0f Mb). Consider lowering filter_seq_mb (currently %.0f).",
        filter_seq_mb,
        filter_seq_mb
      ))
    }
  }
  n_rows_postfilter <- nrow(df)
  
  # ---- Split into groups --------------------------------------------------
  group_indices <- df |>
    dplyr::mutate(.row_idx = dplyr::row_number()) |>
    dplyr::group_by(dplyr::across(all_of(by_columns))) |>
    dplyr::group_split()
  
  
  message(sprintf("Processing %d groups (%s combinations)",
                  length(group_indices), paste(by_columns, collapse = " x ")))
  
  # ---- Step 1: process groups ---------------------------------------------
  
  if (parallel) {
    results <- BiocParallel::bplapply(
      group_indices,
      process_cnv_cluster,
      overlap_method         = overlap_method,
      min_overlap = min_overlap,
      BPPARAM = BiocParallel::MulticoreParam(workers = n_cores)
    )
  } else {
    results <- lapply(
      group_indices,
      process_cnv_cluster,
      overlap_method         = overlap_method,
      min_overlap = min_overlap
    )
  }
  
  # ---- Step 2: assign composite equiv IDs ---------------------------------
  # Composite key = cell_name|chr|cnv_state|local_clique_id
  # Unique by construction — no global counter needed
  # Human readable — carries its own context for debugging
  n_duplicated_vec <- vapply(results, `[[`, integer(1), "n_duplicated")
  total_duplicated <- sum(n_duplicated_vec)

  
  result <- purrr::map(results, ~ {
    .x$rows$cnv_equiv_id <- do.call(
      paste,
      c(
        as.list(.x$rows[, c(by_columns, "local_clique_id")]),
        sep = "|"
      )
    )
    .x$rows$local_clique_id <- NULL
    .x$rows
  }) |>
    dplyr::bind_rows() |>
    dplyr::select(-.row_idx)
  
  removed_log <- purrr::map(results, ~ {
    .x$removed}) %>%
    dplyr::bind_rows()
  
  
  
  # ---- Sanity checks on output --------------------------------------------
  
  # Check 1: row count
  # Expected: post-filter rows + duplicated rows from multi-clique segments
  # Isolated segments are dropped so: result rows = retained + duplicated
  
  n_rows_result   <- nrow(result)
  n_rows_expected <- n_rows_postfilter + total_duplicated
  
  # Count isolated segments — those that were filtered out during clique
  # assignment (no passing overlap partner)
  
  n_removed <- n_rows_postfilter - (n_rows_result - total_duplicated)
  
  if (n_removed != nrow(removed_log)) {
    stop(sprintf(
      "Removed items: output has %d rows supposedly removed while removed log has %d. Check find_maximal_cliques.",
      n_removed,
      nrow(removed_log)
    ))
  }
  
  
  message(sprintf(paste0(
    "Row accounting:\n",
    "  Post-filter input:                                   %d\n",
    "  Isolated (dropped):                                  %d\n",
    "  Single Sequence (dropped):                           %d\n",
    "  Groups with Overlap bellow the threshold (dropped):  %d\n",
    "  0 hits sequences                                     %d\n",
    "  Total of segments (dropped):                         %d\n",
    "  Duplicated (added):                                  %d\n",
    "  Final output:                                        %d"
  ),
  n_rows_postfilter,
  sum(removed_log$removal_reason == "non_overlap (isolated)"),
  sum(removed_log$removal_reason == "single_sequence"),
  sum(removed_log$removal_reason == "non overall (threshold based)"),
  sum(removed_log$removal_reason == "0 hits"),
  nrow(removed_log),
  total_duplicated,
  n_rows_result
  ))
  
  if (n_removed < 0L) {
    stop(sprintf(
      "Row count inconsistency: output has %d more rows than expected. ",
      "Check find_maximal_cliques for duplicate row generation errors.",
      abs(n_removed)
    ))
  }
  
  # Check 2: no NA equiv IDs
  n_na_equiv <- sum(is.na(result$cnv_equiv_id))
  if (n_na_equiv > 0L) {
    stop(sprintf(
      "%d rows have NA cnv_equiv_id after assignment. ",
      "Check process_cnv_cluster for unhandled cases.",
      n_na_equiv
    ))
  }
  
  # Check 3: every cnv_equiv_id within a group
  # should correspond to at least one row — no phantom IDs
  equiv_counts <- result |>
    dplyr::group_by( dplyr::across(all_of(c(by_columns,"cnv_equiv_id")))) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  
  if (any(equiv_counts$n == 0L)) {
    warning("Some cnv_equiv_id values have zero rows — this should not happen.")
  }
  
  message(sprintf(
    "Equivalence complete: %d output rows, %d duplicated across cliques",
    n_rows_result,
    total_duplicated
  ))
  
  list(results_id = result, removed_log = removed_log)
}


#' Summarize reference support for CNV equivalence groups
#'
#' Aggregates CNV equivalence groups and counts how many references support each
#' event within each cell.
#'
#' @param df A data frame containing CNV equivalence assignments.
#'
#' @return A data frame with one row per cell and equivalence group, including
#'   genomic span and reference support counts.
summarize_cnv_support <- function(df) {
  required <- c(
    "cnv_equiv_id", "reference", "cell_name",
    "chr", "cnv_state", "start", "end"
  )
  stopifnot(all(required %in% colnames(df)))
  
  df |>
    group_by(cell_name,chr,cnv_state,cnv_equiv_id) |>
    summarise(
      start        = min(start),
      end          = max(end),
      cnv_length = end - start + 1,
      cnv_length_mb = cnv_length / 1e6,
      n_references = n_distinct(reference),
      references   = paste(sort(unique(reference)), collapse = ","),
      raw_state    = mean(raw_state),
      genes = paste(
        unique(trimws(unlist(
          strsplit(genes, ";")
        ))),
        collapse = ";"
      ),
      n_genes = dplyr::n_distinct(
        trimws(unlist(strsplit(genes, ";")))
      ),
      .groups = "drop"
    ) %>%
    filter(!is.na(cnv_equiv_id))
}



#' Filter CNV events by reference support
#'
#' Keeps only CNV events supported by at least a minimum number of references.
#'
#' @param cnv_events A data frame containing an \code{n_references} column.
#' @param min_references Minimum number of references required to retain an event.
#'
#' @return A filtered data frame of CNV events.
filter_cnv_events <- function(cnv_events, min_references = 2) {
  
  if (!"n_references" %in% colnames(cnv_events)) {
    stop("Missing required column: n_references")
  }
  
  cnv_events |>
    filter(n_references >= min_references)
}


identify_duplicated_segments <- function(
    all_segments,
    consistent,
    coord_cols = c("cell_name", "chr", "cnv_state", "start", "end",
                   "cnv_length", "cnv_length_mb", "reference")
) {
  
  # ---- Validation ---------------------------------------------------------
  missing_seg  <- setdiff(coord_cols, colnames(all_segments))
  missing_cons <- setdiff(c(coord_cols, "cnv_equiv_id"), c(colnames(consistent), "reference"))
  
  if (length(missing_seg) > 0L) {
    stop("all_segments missing columns: ", paste(missing_seg, collapse = ", "))
  }
  if (length(missing_cons) > 0L) {
    stop("consistent missing columns: ", paste(missing_cons, collapse = ", "))
  }
  
  # ---- Find duplicated segments -------------------------------------------
  # Segments appearing under more than one cnv_equiv_id for same coordinates
  # Filtered to only those present in supported events
  duplicated_segments <- all_segments |>
    dplyr::group_by(dplyr::across(dplyr::all_of(coord_cols))) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::ungroup() |>
    dplyr::filter(cnv_equiv_id %in% consistent$cnv_equiv_id) #to ensure 
  
  not_duplicated <- consistent |>
    dplyr::filter(!(cnv_equiv_id %in% duplicated_segments$cnv_equiv_id)) |>
    dplyr::mutate(merge_group_id = cnv_equiv_id)
  
  # ---- Split for processing -----------------------------------------------
  duplicate_list <- duplicated_segments |>
    dplyr::group_by(dplyr::across(dplyr::all_of(coord_cols))) |>
    dplyr::group_split()
  
  # ---- Sanity checks ------------------------------------------------------
  if (!all(purrr::map_lgl(duplicate_list, ~ nrow(.x) >= 2L))) {
    stop("Some duplicate groups have fewer than 2 members — check coord_cols definition.")
  }
  
  message(sprintf(paste0(
    "Duplicated segment summary:                                 \n",
    "  Duplicate segments (add into multiple calls):  %d\n",
    "  Duplicate groups to resolve:                   %d"
  ),
  length(unique(duplicated_segments$cnv_equiv_id)),
  length(duplicate_list)
  ))
  
  list(
    duplicated      = duplicated_segments,
    not_duplicated  = not_duplicated,
    duplicate_list  = duplicate_list
  )
}


resolve_duplicate_overlaps <- function(
    df,
    consistent,
    min_overlap    = 0.6,
    overlap_method = "reciprocal",
    clique_mode    = c("connected", "complete")
) {
  
  clique_mode <- match.arg(clique_mode)
  
  grp <- consistent[consistent$cnv_equiv_id %in% df$cnv_equiv_id, ]
  
  if (nrow(grp) == 0L) {
    warning("No matching events found in consistent for this duplicate group.")
    return(grp)
  }
  
  n <- nrow(grp)
  
  gr <- GenomicRanges::GRanges(
    seqnames = grp$chr,
    ranges   = IRanges::IRanges(start = grp$start, end = grp$end),
    strand   = "*"
  )
  
  hits <- GenomicRanges::findOverlaps(gr, gr, type = "any", select = "all")
  hits <- hits[S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits)]
  
  # No overlaps at all — every event keeps its own ID
  if (length(hits) == 0L) {
    grp$merge_group_id <- grp$cnv_equiv_id
    return(grp)
  }
  
  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)
  
  scores <- compute_overlap(
    q_start = grp$start[q_idx],
    q_end   = grp$end[q_idx],
    s_start = grp$start[s_idx],
    s_end   = grp$end[s_idx],
    method  = overlap_method
  )
  
  passing <- scores >= min_overlap
  q_pass  <- q_idx[passing]
  s_pass  <- s_idx[passing]
  
  # Preassign  every row starts with its own ID
  grp$merge_group_id <- grp$cnv_equiv_id
  
  # ---- Graph construction -------------------------------------------------
  g <- igraph::graph_from_edgelist(
    cbind(as.character(q_pass), as.character(s_pass)),
    directed = FALSE
  )
  
  # Add isolated vertices explicitly
  n_connected <- length(unique(c(q_pass, s_pass)))
  g <- igraph::add_vertices(g, n_connected - igraph::vcount(g))
  
  # ---- Component assignment -----------------------------------------------
  if (clique_mode == "connected") {
    comps       <- igraph::components(g)$membership
    node_groups <- split(
      as.integer(names(comps)),
      as.integer(comps)
    )
    
  } else {
    cliques     <- igraph::max_cliques(g)
    node_groups <- lapply(cliques, function(cl) as.integer(names(cl)))
  }
  
  # Initialise collector before the loop
  ambiguous_nodes <- list()
  # ---- Build merge_group_id -----------------------------------------------
  # ---- Assign merge_group_id — same logic for both modes ------------------
  for (node_idx in node_groups) {
    if (length(node_idx) > 1L) {
      clique_nums <- stringr::str_extract(
        grp$cnv_equiv_id[node_idx], "[^|]+$"
      )
      prefix    <- stringr::str_remove(
        grp$cnv_equiv_id[node_idx[1]], "\\|[^|]+$"
      )
      merged_id <- paste0(
        prefix, "|",
        paste(sort(unique(clique_nums)), collapse = "_")
      )
      # First clique wins for ambiguous segments in complete mode
      already_assigned <- node_idx[
        grp$merge_group_id[node_idx] != grp$cnv_equiv_id[node_idx]
      ]
      
      if (length(already_assigned) > 0L) {
        ambiguous_nodes <- c(ambiguous_nodes, list(data.frame(
          node_idx     = already_assigned,
          cnv_equiv_id = grp$cnv_equiv_id[already_assigned],
          kept_id      = grp$merge_group_id[already_assigned],
          skipped_id   = merged_id
        )))
      }
      
      unassigned <- node_idx[
        grp$merge_group_id[node_idx] == grp$cnv_equiv_id[node_idx]
      ]
      
      grp$merge_group_id[unassigned] <- merged_id
    }
  }
  
  # Report once after loop — only if ambiguous nodes were found
  if (length(ambiguous_nodes) > 0L) {
    ambiguous_df <- dplyr::bind_rows(ambiguous_nodes)
    message(sprintf(
      "%d segments appeared in multiple overlap and they were assigned to the first overlap cluster.",
      nrow(ambiguous_df)
    ))
  }
  grp
}


resolve_shared_cliques <- function(
    all_segments,
    consistent,
    coord_cols     = c("cell_name", "chr", "cnv_state", "start", "end",
                       "cnv_length", "cnv_length_mb", "reference"),
    min_overlap    = 0.6,
    overlap_method = "reciprocal",
    clique_mode    = c("connected", "complete"),
    parallel       = FALSE,
    n_cores        = 1L
) {
  
  clique_mode <- match.arg(clique_mode)
  
  # ---- Step 1: identify duplicates ----------------------------------------
  identified <- identify_duplicated_segments(
    all_segments = all_segments,
    consistent   = consistent,
    coord_cols   = coord_cols
  )
  
  not_duplicated <- identified$not_duplicated
  duplicate_list <- identified$duplicate_list
  
  if (length(duplicate_list) == 0L) {
    message("No duplicated segments found — returning consistent table unchanged.")
    consistent$merge_group_id <- consistent$cnv_equiv_id
    return(consistent)
  }
  
  # ---- Step 2: resolve overlaps -------------------------------------------
  if (parallel) {
    resolved_list <- furrr::future_map(
      duplicate_list,
      resolve_duplicate_overlaps,
      consistent     = consistent,
      min_overlap    = min_overlap,
      overlap_method = overlap_method,
      clique_mode    = clique_mode,
      .options       = furrr::furrr_options(seed = TRUE)
    )
  } else {
    resolved_list <- purrr::map(
      duplicate_list,
      resolve_duplicate_overlaps,
      consistent     = consistent,
      min_overlap    = min_overlap,
      overlap_method = overlap_method,
      clique_mode    = clique_mode
    )
  }
  
  
  # ---- Step 3: summarise merged groups ------------------------------------
  new_merge_id_duplicates <- dplyr::bind_rows(resolved_list) |>
    dplyr::group_by(merge_group_id, chr, cnv_state) |>
    dplyr::summarise(
      cell_name    = dplyr::first(cell_name),
      start        = min(start),
      end          = max(end),
      references   = paste(
        sort(unique(unlist(strsplit(references, ",")))),
        collapse = ","
      ),
      n_references = length(unique(unlist(strsplit(references, ",")))),
      # Recompute length from updated coordinates
      cnv_length    = end - start + 1L,
      cnv_length_mb = (end - start + 1L) / 1e6,
      raw_state     = mean(raw_state),
      genes = paste(
        unique(trimws(unlist(
          strsplit(genes, ";")
        ))),
        collapse = ";"
      ),
      n_genes = dplyr::n_distinct(
        trimws(unlist(strsplit(genes, ";")))
      ),
      .groups       = "drop"
    )
  
  # ---- Step 4: combine with non-duplicated --------------------------------
  combined <- dplyr::bind_rows(not_duplicated, new_merge_id_duplicates) |>
    dplyr::select(-cnv_equiv_id)
  
  # ---- Sanity checks ------------------------------------------------------
  n_merged    <- nrow(new_merge_id_duplicates)
  n_not_dup   <- nrow(not_duplicated)
  n_combined  <- nrow(combined)
  n_original  <- nrow(consistent)
  
  # Count how many groups were actually merged vs kept separate
  n_truly_merged <- sum(stringr::str_detect(
    new_merge_id_duplicates$merge_group_id, "_"
  ))
  n_kept_separate <- n_merged - n_truly_merged
  
  message(sprintf(paste0(
    "Resolution summary:\n",
    "  Original supported events:                      %d\n",
    "  Non-duplicated (unchanged):                     %d\n",
    "  Duplicate groups resolved:                      %d\n",
    "    → merged into one event (within each group):  %d\n",
    "    → kept separate:                              %d\n",
    "  Final combined events:                          %d"
  ),
  n_original,
  n_not_dup,
  length(duplicate_list),
  n_truly_merged,
  n_kept_separate,
  n_combined
  ))
  
  if (n_combined > n_original) {
    warning(sprintf(
      "Combined table (%d rows) exceeds original (%d rows). ",
      "Check for unintended row duplication in resolve_duplicate_overlaps.",
      n_combined, n_original
    ))
  }
  
  combined
}


# ── Core density computation ───────────────────────────────────────────────────
# Used by both filter and merge steps
compute_segment_density <- function(
    df,
    coding_gr,
    expressed_gr,         
    coding_expressed_set,  
    start_col = "start",
    stop_col  = "stop"
) {
  
  if (!"genes" %in% colnames(df)) {
    stop(
      "df must contain 'genes' column from ",
      "collapse_genes_to_cnv_segments()."
    )
  }
  
  seg_gr <- GenomicRanges::GRanges(
    seqnames = df$chr,
    ranges   = IRanges::IRanges(
      start = df[[start_col]],
      end   = df[[stop_col]]
    )
  )
  seg_gr$seg_idx <- seq_len(nrow(df))
  
  # ── n_total_protein_coding ────────────────────────────────────────────────
  # ALL coding + lncRNA genes in segment from GTF
  # coordinate-based, not expression-based
  hits_coding <- GenomicRanges::findOverlaps(
    seg_gr, coding_gr, type = "any")
  
  coding_counts <- if (length(hits_coding) > 0) {
    data.frame(
      seg_idx = S4Vectors::queryHits(hits_coding),
      gene    = coding_gr$gene[
        S4Vectors::subjectHits(hits_coding)]
    ) %>%
      dplyr::group_by(seg_idx) %>%
      dplyr::summarise(
        n_total_protein_coding = dplyr::n(),
        .groups = "drop"
      )
  } else {
    data.frame(
      seg_idx                = integer(0),
      n_total_protein_coding = integer(0)
    )
  }
  
  # ── n_expressing_genes ────────────────────────────────────────────────────
  # Genes in @gene_order overlapping segment
  # = genes inferCNV used (passed cutoff)
  # same for all cells — sample level ✅
  hits_expressed <- GenomicRanges::findOverlaps(
    seg_gr, expressed_gr, type = "any")
  
  expressed_counts <- if (length(hits_expressed) > 0) {
    data.frame(
      seg_idx = S4Vectors::queryHits(hits_expressed)
    ) %>%
      dplyr::count(seg_idx,
                   name = "n_expressing_genes")
  } else {
    data.frame(
      seg_idx            = integer(0),
      n_expressing_genes = integer(0)
    )
  }
  
  # ── n_coding_expressed ────────────────────────────────────────────────────
  # Per-cell: genes in THIS cell's gene_list
  # that are also in coding_expressed_set
  # coding_expressed_set = intersect(
  #   coding genes from GTF,
  #   expressed genes from @gene_order
  # )
  # → protein coding genes that inferCNV
  #   actually detected in this sample ✅
  # → then filtered to which ones THIS
  #   cell has in its segment ✅
  
  result <- df %>%
    dplyr::mutate(
      seg_idx   = seq_len(dplyr::n()),
      gene_list = strsplit(genes, ";")
    ) %>%
    dplyr::left_join(coding_counts,
                     by = "seg_idx") %>%
    dplyr::left_join(expressed_counts,
                     by = "seg_idx") %>%
    dplyr::mutate(
      
      n_total_protein_coding = tidyr::replace_na(
        n_total_protein_coding, 0L),
      n_expressing_genes     = tidyr::replace_na(
        n_expressing_genes, 0L),
      
      cnv_length_mb = (
        .data[[stop_col]] -
        .data[[start_col]] + 1) / 1e6,
      
      # Per-cell: how many genes in THIS
      # cell's segment are coding AND expressed
      # in this sample ✅
      n_coding_expressed = purrr::map_int(
        gene_list,
        function(g) {
          if (length(g) == 0) return(0L)
          sum(g %in% coding_expressed_set)
        }
      ),
      
      coding_density = round(
        n_total_protein_coding / cnv_length_mb, 1),
      
      expressed_density = round(
        n_expressing_genes / cnv_length_mb, 1),
      
      # pct of coding genes in this segment
      # that THIS cell expressed ✅
      pct_coding_expressed = round(
        100 * n_coding_expressed /
          pmax(n_total_protein_coding, 1), 1)
      
    ) %>%
    dplyr::select(-seg_idx, -gene_list)
  
  return(result)
}




# ── Apply density filter ───────────────────────────────────────────────────────
apply_density_filter <- function(
    df,
    pct_max          = 45,
    pct_floor        = 30,
    min_expr_density = 1.5,
    min_coding_density = 1.0
) {
  df %>%
    dplyr::mutate(
      pct_threshold = pmax(
        pct_floor,
        pct_max / sqrt(pmax(expressed_density, 0.01))
      ),
      passes_coding_density    = coding_density    >= min_coding_density,
      passes_expressing_density = expressed_density >= min_expr_density,
      passes_pct               = pct_coding_expressed >= pct_threshold
    ) %>%
    dplyr::filter(
      passes_coding_density,
      passes_expressing_density,
      passes_pct
    ) %>%
    dplyr::select(
      -pct_threshold,
      -passes_coding_density,
      -passes_expressing_density,
      -passes_pct
    )
}



# ── Pre-merge filter ───────────────────────────────────────────────────────────
filter_segments_by_gene_density <- function(
    collapse_df,
    gene_order,
    coding_gr,
    coding_expressed_set,
    pct_max            = 45,
    pct_floor          = 30,
    min_expr_density   = 1.5,
    min_coding_density = 1.0
) {
  
  # Build expressed_gr from gene_order
  expressed_table <- as.data.frame(gene_order)
  expressed_table$gene <- rownames(expressed_table)
  colnames(expressed_table) <- c("chr", "start",
                                  "stop", "gene")
  expressed_gr <- GenomicRanges::GRanges(
    seqnames = expressed_table$chr,
    ranges   = IRanges::IRanges(
      start = expressed_table$start,
      end   = expressed_table$stop
    )
  )
  expressed_gr$gene <- expressed_table$gene
  
  cat("  expressed_gr:",
      length(expressed_gr), "genes\n")
  
  n_before <- nrow(collapse_df)
  
  if (nrow(coding_counts) == 0) {
    message("  No coding genes found — skipping")
    return(collapse_df)
  }
  
  result <- collapse_df %>%
    compute_segment_density(
      coding_gr            = coding_gr,
      expressed_gr         = expressed_gr,
      coding_expressed_set = coding_expressed_set
    ) %>%
    apply_density_filter(
      pct_max            = pct_max,
      pct_floor          = pct_floor,
      min_expr_density   = min_expr_density,
      min_coding_density = min_coding_density
    )
  
  message(sprintf(
    "  Pre-merge filter: %d → %d (removed %d)",
    n_before, nrow(result),
    n_before - nrow(result)
  ))
  
  return(result)
}


filter_segments_by_gene_density <- function(
    collapse_df,
    gene_order,
    coding_gr,
    coding_expressed_set,
    pct_max            = 45,
    pct_floor          = 30,
    min_expr_density   = 1.5,
    min_coding_density = 1.0
) {
  
  # Build expressed_gr from gene_order
  expressed_table <- as.data.frame(gene_order)
  expressed_table$gene <- rownames(expressed_table)
  colnames(expressed_table) <- c("chr", "start",
                                  "stop", "gene")
  expressed_gr <- GenomicRanges::GRanges(
    seqnames = expressed_table$chr,
    ranges   = IRanges::IRanges(
      start = expressed_table$start,
      end   = expressed_table$stop
    )
  )
  expressed_gr$gene <- expressed_table$gene
  
  cat("  expressed_gr:", length(expressed_gr), "genes\n")
  
  n_before <- nrow(collapse_df)
  
  # ── Compute density ───────────────────────────────────────────────────────
  result <- compute_segment_density(
    df                   = collapse_df,
    coding_gr            = coding_gr,
    expressed_gr         = expressed_gr,
    coding_expressed_set = coding_expressed_set
  )
  
  # ── Checking if there are any results whatsoever ────────────────────────────
  if (all(result$n_total_protein_coding == 0)) {
    message(
      "  No coding genes overlap any segment — ",
      "skipping filter, returning unchanged.\n",
      "  Likely chromosome naming mismatch: ",
      "check seg_gr vs coding_gr seqnames"
    )
    return(collapse_df)
  }
  
  # ── Apply filter ──────────────────────────────────────────────────────────
  result <- result %>%
    apply_density_filter(
      pct_max            = pct_max,
      pct_floor          = pct_floor,
      min_expr_density   = min_expr_density,
      min_coding_density = min_coding_density
    )
  
  message(sprintf(
    "  Pre-merge filter: %d → %d (removed %d)",
    n_before, nrow(result),
    n_before - nrow(result)
  ))
  
  return(result)
}
# ── Merge + post-merge density check ──────────────────────────────────────────
merge_and_density_check <- function(
    df,
    gene_order,
    coding_gr,
    coding_expressed_set,
    max_gap_mb         = 10,
    pct_max            = 45,
    pct_floor          = 30,
    min_expr_density   = 1.5,
    min_coding_density = 1.0
) {
  
  # Build expressed_gr
  expressed_table <- as.data.frame(gene_order)
  expressed_table$gene <- rownames(expressed_table)
  colnames(expressed_table) <- c("chr", "start",
                                  "stop", "gene")
  expressed_gr <- GenomicRanges::GRanges(
    seqnames = expressed_table$chr,
    ranges   = IRanges::IRanges(
      start = expressed_table$start,
      end   = expressed_table$stop
    )
  )
  
  n_before   <- nrow(df)
  max_gap_bp <- max_gap_mb * 1e6
  
  cat("=== Merge + density check ===\n")
  cat("  max_gap:     ", max_gap_mb, "Mb\n")
  cat("  min_density: ", min_expr_density, "\n")
  cat("  Input:       ", n_before, "\n")
  
  # ── Step 1: Merge ──────────────────────────────────────────────────────────
  cat("\n[1] Merging within", max_gap_mb, "Mb...\n")
  
  merged_df <- df %>%
    dplyr::group_by(
      reference, cell_name, chr, state
    ) %>%
    dplyr::arrange(start, .by_group = TRUE) %>%
    dplyr::mutate(
      gap_bp    = start - dplyr::lag(stop),
      new_block = is.na(gap_bp) |
                  gap_bp > max_gap_bp,
      merge_id  = cumsum(new_block)
    ) %>%
    dplyr::group_by(
      reference, cell_name, chr,
      state, merge_id
    ) %>%
    dplyr::summarise(
      start      = min(start),
      end       = max(stop),
      n_segments = dplyr::n(),
      raw_state  = mean(raw_state, na.rm = TRUE),
      genes      = paste(
      unique(trimws(unlist(
        strsplit(genes, ";")
      ))),
      collapse = ";"
    ),
    n_genes    = dplyr::n_distinct(trimws(unlist(
      strsplit(genes, ";")
    ))),
    genes_total_length = sum(
      genes_total_length, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::select(-merge_id) %>%
    dplyr::rename(cnv_state = state)
  
  cat("  After merge:", nrow(merged_df), "\n")
  
  # ── Step 2: Post-merge density check ──────────────────────────────────────
  cat("\n[2] Post-merge density check...\n")
  
  result <- merged_df %>%
    compute_segment_density(
      coding_gr            = coding_gr,
      expressed_gr         = expressed_gr,
      coding_expressed_set = coding_expressed_set,
      stop_col  = "end"
    ) %>%
    apply_density_filter(
      pct_max            = pct_max,
      pct_floor          = pct_floor,
      min_expr_density   = min_expr_density,
      min_coding_density = min_coding_density
    )
  
  message(sprintf(paste0(
    "Merge + density check:\n",
    "  Input:         %d\n",
    "  After merge:   %d\n",
    "  After density: %d\n",
    "  Total removed: %d (%.1f%%)"
  ),
  n_before, nrow(merged_df), nrow(result),
  n_before - nrow(result),
  100 * (n_before - nrow(result)) / n_before
  ))
  
  return(result)
}



add_gene_length_stats <- function(
    df,
    expressed_gr,
    start_col = "start",
    end_col   = "end"
) {
  
  cat("=== Adding gene length stats ===\n")
  cat("  Input rows:", nrow(df), "\n")
  
  seg_gr <- GenomicRanges::GRanges(
    seqnames = df$chr,
    ranges   = IRanges::IRanges(
      start = df[[start_col]],
      end   = df[[end_col]]
    )
  )
  seg_gr$seg_idx <- seq_len(nrow(df))
  
  hits <- GenomicRanges::findOverlaps(
    seg_gr, expressed_gr, type = "any")
  
  if (length(hits) == 0) {
    message("  No overlaps — returning zeros")
    return(df %>%
      dplyr::mutate(
        genes_total_length = 0L,
        gene_coverage_pct  = 0
      ))
  }
  
  s_idx <- S4Vectors::subjectHits(hits)
  
  
  gene_stats <- data.frame(
    seg_idx     = S4Vectors::queryHits(hits),
    gene_length = GenomicRanges::width(
      expressed_gr)[s_idx]
  ) %>%
    dplyr::group_by(seg_idx) %>%
    dplyr::summarise(
      genes_total_length = sum(gene_length,
                               na.rm = TRUE),
      .groups = "drop"
    )
  
  df %>%
    dplyr::mutate(seg_idx = seq_len(dplyr::n())) %>%
    dplyr::left_join(gene_stats, by = "seg_idx") %>%
    dplyr::mutate(
      genes_total_length = tidyr::replace_na(
        genes_total_length, 0L),
     
      gene_coverage_pct  = round(
        100 * genes_total_length /
          pmax((.data[[end_col]] -
                .data[[start_col]] + 1L), 1), 1)
    ) %>%
    dplyr::select(-seg_idx)
}


compute_all_segment_stats <- function(
    df,
    coding_gr,
    gene_order,
    coding_expressed_set,
    start_col = "start",
    stop_col  = "stop"
) {
  
  expressed_table <- as.data.frame(gene_order)
  expressed_table$gene <- rownames(expressed_table)
  colnames(expressed_table) <- c("chr", "start",
                                "stop", "gene") 
  expressed_gr <- GenomicRanges::GRanges(
  seqnames = expressed_table$chr,
  ranges   = IRanges::IRanges(
    start = expressed_table$start,
    end   = expressed_table$stop))
    
 expressed_gr$gene <- expressed_table$gene
 
 cat("=== Computing all segment stats ===\n")
  
  # ── Step 1: Gene density (coding + expressed) ─────────────────────────────
  cat("\n[1] Gene density stats...\n")
  result <- compute_segment_density(
    df                   = df,
    coding_gr            = coding_gr,
    expressed_gr         = expressed_gr,
    coding_expressed_set = coding_expressed_set,
    start_col            = start_col,
    stop_col             = stop_col
  )
  
  # ── Step 2: Gene length stats ─────────────────────────────────────────────
  cat("\n[2] Gene length stats...\n")
  print(colnames(result))
  result <- add_gene_length_stats(
    df           = result,
    expressed_gr = expressed_gr,
    start_col    = start_col,
    end_col      = stop_col
  )
  
  cat("\n=== Done ===\n")
  cat("  Rows:                    ", nrow(result), "\n")
  cat("  Median expressed_density:",
      median(result$expressed_density,
             na.rm = TRUE), "\n")
  cat("  Median coding_density:   ",
      median(result$coding_density,
             na.rm = TRUE), "\n")
  cat("  Median gene_coverage_pct:",
      median(result$gene_coverage_pct,
             na.rm = TRUE), "\n")
  
  return(result)
}


#' Run a fast CNV event consolidation pipeline
#'
#' Runs a CNV processing workflow from gene-level calls to merged CNV events
#' with reference support summaries.
#'
#' @param gene_level_df A gene-level CNV data frame.
#' @param max_gap Maximum genomic gap allowed when merging nearby segments.
#' @param min_overlap Minimum reciprocal overlap for equivalence
#'   assignment.
#' @param min_references Minimum number of references required to keep a CNV.
#' @param overlap_method Select the overlap method
run_fast_cnv_pipeline <- function(
    gene_level_df,
    max_gap = 100000,
    min_overlap_consistent_calls = 0.5,
    min_overlap_multiple_nodes = 0.6,
    filter_seq_mb_init = 5,
    filter_seq_mb_equiv = 7,
    min_references = 2,
    overlap_method_equiv_cnv_call_merge = "reciprocal",
    overlap_method_equiv_cnv_after_filter = "reciprocal",
    parallel = F,
    cores = 1L,
    clique_mode_consistent = "connected",
    removed_log_return = F,
    mode = "within",
    metadata,
    remove_ref = T,
    gene_order           = NULL,
    coding_gr            = NULL,
    coding_expressed_set = NULL,
    pct_max              = 45,
    pct_floor            = 30,
    min_expr_density     = 1.5,
    min_coding_density   = 1.0,
    max_gap_mb           = 10
) {

  has_density_params <- !is.null(gene_order) &&
                        !is.null(coding_gr)  &&
                        !is.null(coding_expressed_set)
  
  if (!has_density_params) {
    stop(
      "gene_order, coding_gr and coding_expressed_set ",
      "are required for run_fast_cnv_pipeline.\n",
      "These are needed for gene density filter ",
      "and post-merge density check."
    )
  }
  
  
  message("→ Collapsing genes to segments")
  segments <- collapse_genes_to_cnv_segments(gene_cnv_df = gene_level_df)
  
  
    segments <- filter_segments_by_gene_density(
    collapse_df          = segments,
    gene_order           = gene_order,
    coding_gr            = coding_gr,
    coding_expressed_set = coding_expressed_set,
    pct_max              = pct_max,
    pct_floor            = pct_floor,
    min_expr_density     = min_expr_density,
    min_coding_density   = min_coding_density
  )
  
    merged <- merge_and_density_check(
    df                   = segments,
    gene_order           = gene_order,
    coding_gr            = coding_gr,
    coding_expressed_set = coding_expressed_set,
    max_gap_mb           = max_gap_mb,
    pct_max              = pct_max,
    pct_floor            = pct_floor,
    min_expr_density     = min_expr_density,
    min_coding_density   = min_coding_density
  )
  
  message("→ Removing reference cells")
  filt_segments <- filt_remove_refs_cells(merged, metadata, filter_seq_mb = filter_seq_mb_init, mode,
                                          remove_ref = remove_ref)
  
  message("→ Assigning CNV equivalence")
  equiv <- assign_cnv_equivalence(
    df = merged,
    min_overlap = min_overlap_consistent_calls,
    overlap_method         = overlap_method_equiv_cnv_call_merge,
    filter_seq_mb          = filter_seq_mb_equiv,
    parallel               = parallel,
    n_cores = cores
  )
  
  table_with_equiv_id <- equiv$results_id
  
 
  cnv_events <- summarize_cnv_support(table_with_equiv_id)
  
  message("→ Filtering CNVs by reference support")
  
  supported_events <- filter_cnv_events(
    cnv_events,
    min_references = min_references
  )

  final_consistent_events <- resolve_shared_cliques(
    all_segments = table_with_equiv_id, 
    consistent = supported_events, 
    coord_cols = c("cell_name", "chr", "cnv_state", "start", "end","cnv_length", "cnv_length_mb", "reference"),
    min_overlap    = min_overlap_multiple_nodes,
    overlap_method = overlap_method_equiv_cnv_after_filter,
    clique_mode    = clique_mode_consistent,
    parallel       = parallel,
    n_cores        = cores
  )
  
  final_consistent_events <- compute_all_segment_stats(
    final_consistent_events,
    coding_gr,
    gene_order,
    coding_expressed_set,
    start_col = "start",
    stop_col  = "end")
    
  if(removed_log_return){
    list(
      cnvs_per_segment   =  table_with_equiv_id,  # one row per segment, with equiv IDs - IDs which tells which CNVs overlap each other
      cnvs_summarized    =  cnv_events,      # one row per equiv group, with support counts
      cnvs_supported     =  supported_events,    # filtered to min_references threshold
      removed_log        = equiv$removed_log,    #segments that were removed for numerous reason
      cnvs_supported_overlaped = add_metadata(df = final_consistent_events, mode_name = mode,  metadata =metadata)  #segments that are in multiple calls are assess for overlap and if they pass the threshold the segment increases
    )
  } else{
    list(
      cnvs_per_segment   = table_with_equiv_id,  # one row per segment, with equiv IDs - IDs which tells which CNVs overlap each other
      cnvs_summarized    =  cnv_events,          # one row per equiv group, with support counts
      cnvs_supported     =  supported_events,    # filtered to min_references threshold
      cnvs_supported_overlaped = add_metadata(df = final_consistent_events, mode_name = mode,  metadata = metadata)  #segments that are in multiple calls are assess for overlap and if they pass the threshold the segment increases
    )
  }

}















############################################################################################################-
############################# Across and Within Cell type Integration Approach##############################
############################################################################################################-

std_events <- function(tbl, dataset, mode){
  tbl %>%
    mutate(cell_type = dataset, mode = mode,
           ds_cell = paste(cell_type, cell_name, sep="|"))
}














