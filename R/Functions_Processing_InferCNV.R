#' @title Score_system
#'
#' @description
#' 
#' Functions that perform further processing from the CNV performing a final filtering and Confidence Score in a Single Tool Approach
#' 
#' @author Pedro Granjo
#' @date 13-03-2026
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


load_infercnv_data <- function(base_dir, ref_dirs, pattern = "^run\\.final") {
  
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
  
  load_and_prepare_infercnv_reference(infer_objs_1)
}


load_tool_data <- function(
    tool    = c("infercnv", "scevan", "copykat"),
    base_dir,
    ref_dirs,
    pattern
) {
  tool <- match.arg(tool)
  
  switch(tool,
         "infercnv" = load_infercnv_data(base_dir, ref_dirs, pattern)#,
         #"scevan"   = load_scevan_data(base_dir, ref_dirs, pattern),
         #"copykat"  = load_copykat_data(base_dir, ref_dirs, pattern)
  )
  
}


process_tool_cnv_runs <- function(
    base_dir,
    mode                                = c("within", "across"),
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
    metadata
) {
  
  # ---- Validate base_dir --------------------------------------------------
  if (!dir.exists(base_dir)) {
    stop("base_dir does not exist: ", base_dir)
  }
  
  # ---- Validate modes directories exist -----------------------------------
  missing_mode_dirs <- mode[!dir.exists(file.path(base_dir, mode))]
  if (length(missing_mode_dirs) > 0L) {
    stop(
      "Mode directories not found in base_dir: ",
      paste(missing_mode_dirs, collapse = ", ")
    )
  }
  
  # ---- Run pipeline per mode and cell type --------------------------------
  cell_dir   <- file.path(base_dir, mode)
  cell_types <- list.files(cell_dir)
    
  if (length(cell_types) == 0L) {
      warning(sprintf("No cell types found in %s — skipping mode '%s'.", 
                      cell_dir, mode))
      return(NULL)
  }
    
  ct_results <- purrr::map(cell_types, \(ct) {
    
      ct_dir   <- file.path(cell_dir, ct)
      refs_dir <- list.files(ct_dir, full.names = T)
      refs_dir <- list.files(ct_dir)[dir.exists(file.path(refs_dir))]

      if (length(refs_dir) == 0L) {
        warning(sprintf("No reference directories found in %s — skipping.", ct_dir))
        return(NULL)
      }
      
      message(sprintf("Processing mode='%s', cell_type='%s', refs=%s",
                      mode, ct, paste(refs_dir, collapse = ", ")))
      
      # Load tool data
      tools_data <- load_tool_data(
        tool     = tool,
        base_dir = ct_dir,
        ref_dirs = refs_dir,
        pattern  = pattern
      )
      
      # Run CNV pipeline
      run_fast_cnv_pipeline(
        gene_level_df                         = tools_data,
        max_gap                               = max_gap,
        min_overlap_consistent_calls          = min_overlap_consistent_calls,
        min_overlap_multiple_nodes            = min_overlap_multiple_nodes,
        filter_seq_mb_init                    =  filter_seq_mb_init,
        filter_seq_mb_equiv                   = filter_seq_mb_equiv,
        min_references                        = min_references,
        overlap_method_equiv_cnv_call_merge   = overlap_method_equiv_cnv_call_merge,
        overlap_method_equiv_cnv_after_filter = overlap_method_equiv_cnv_after_filter,
        parallel                              = parallel,
        cores                                 = cores,
        clique_mode_consistent                = clique_mode_consistent,
        removed_log_return                    = removed_log_return,
        metadata = metadata,
        mode = mode
      )
    })
    
    names(ct_results) <- cell_types
    ct_results
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
load_and_prepare_infercnv_reference <- function(infercnv_list) {
  
  refs <- names(infercnv_list)
  
  res <- lapply(seq_along(infercnv_list), function(i) {
    
    infercnv_obj <- infercnv_list[[i]]
    reference <- refs[i]
    melt_expr_to_long(as.data.frame(infercnv_obj[[1]])) |>
      attach_gene_order(infercnv_obj[[2]]) |>
      discretize_cnv_state_infer_cnv() |>
      dplyr::mutate(reference = reference)
  })
  
  dplyr::bind_rows(res)
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
                   state = (state)
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


filt_remove_refs_cells <- function(df, metadata, filter_seq_mb, mode) {
  
  # ---- Pre-filter ---------------------------------------------------------
  n_input <- nrow(df)
  
  df <- df |>
    dplyr::arrange(reference, cell_name, chr, state, start) |>
    dplyr::mutate(
      cnv_length    = as.numeric(stop) - as.numeric(start) + 1,
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
  
  n_before_ref_removal <- nrow(df_joined)
  
  if (mode == "within") {
    
   # df_joined <- df_joined |>
    #  dplyr::filter(!(reference == split_group)) |>
     # dplyr::select(-dplyr::any_of("split_group"))
    
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
  ))
  
  df_joined
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
merge_nearby_regions <- function(df, max_gap = 100000L, debug = FALSE) {
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
  warning(sprintf(
    "No reduction after merging: before = %d, after = %d. using the max_gap = %d.",
    as.integer(n_postfilter),
    as.integer(n_merged),
    as.integer(max_gap)
  ))
  
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
  
  # Check 7: overlap detection 
  if (debug) {
    overlap <- merged_df |>
      dplyr::group_by(reference, cell_name, chr, cnv_state) |>
      dplyr::mutate(overlaps_prev = start <= dplyr::lag(end)) |>
      dplyr::filter(!is.na(overlaps_prev) & overlaps_prev)
    
    if (nrow(overlap) > 0L) {
      stop(sprintf(
        "%d overlapping segment pair(s) detected after merging. Check merging logic.",
        nrow(overlap)
      ))
    }
    message("Debug: no overlapping segments detected after merging.")
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
  n   <- nrow(grp)

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

  res <- find_maximal_cliches(q_pass, s_pass,grp)
  
  
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
        "No segments remain after length filter (> %.0f Mb). ",
        "Consider lowering filter_seq_mb (currently %.0f).",
        filter_seq_mb, filter_seq_mb
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
    metadata 
) {
  
  message("→ Collapsing genes to segments")
  segments <- collapse_genes_to_cnv_segments(gene_cnv_df = gene_level_df)
  
  message("→ Removing reference cells")
  filt_segments <- filt_remove_refs_cells(segments, metadata, filter_seq_mb = filter_seq_mb_init, mode)
  
  message("→ Merging nearby CNVs")
  merged <- merge_nearby_regions(df = filt_segments, max_gap = max_gap)

 
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

  final_consistent_events <-resolve_shared_cliques(
    all_segments = table_with_equiv_id, 
    consistent = supported_events, 
    coord_cols = c("cell_name", "chr", "cnv_state", "start", "end","cnv_length", "cnv_length_mb", "reference"),
    min_overlap    = min_overlap_multiple_nodes,
    overlap_method = overlap_method_equiv_cnv_after_filter,
    clique_mode    = clique_mode_consistent,
    parallel       = parallel,
    n_cores        = cores
  )
  
    
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






#' Compute overlap in base pairs between two intervals
#'
#' @param a_start,a_end Start and end coordinates of the first interval.
#' @param b_start,b_end Start and end coordinates of the second interval.
#'
#' @return Integer overlap length in base pairs.
overlap_bp <- function(a_start, a_end, b_start, b_end) {
  max(0, min(a_end, b_end) - max(a_start, b_start))
}



########################################################


#' Classify a CNV by chromosome arm overlap
#'
#' Assigns a CNV event to an arm-level category based on overlap with
#' chromosome arm annotations.
#'
#' @param cnv_row A single-row data frame representing one CNV.
#' @param chromosome_arms A data frame describing chromosome arm intervals.
classify_single_cnv <- function(cnv_row, chromosome_arms) {
  
  required_cnv  <- c("chr", "start", "end")
  required_arms <- c("chr", "arm_start", "arm_end", "arm")
  
  missing_cnv  <- setdiff(required_cnv,  colnames(cnv_row))
  missing_arms <- setdiff(required_arms, colnames(chromosome_arms))
  
  if (length(missing_cnv)  > 0) stop("cnv_row missing columns: ",        paste(missing_cnv,  collapse = ", "))
  if (length(missing_arms) > 0) stop("chromosome_arms missing columns: ", paste(missing_arms, collapse = ", "))
  
  
  arms <- chromosome_arms[chromosome_arms$chr == cnv_row$chr, ]
  
  hit_p <- FALSE
  hit_q <- FALSE
  hit_c <- FALSE
  
  if (nrow(arms) == 0) {
    return(NA_character_)
  }
  
  for (j in seq_len(nrow(arms))) {
    
    ov <- overlap_bp(
      cnv_row$start, cnv_row$end,
      arms$arm_start[j], arms$arm_end[j]
    )
    
    if (ov > 0) {
      if (arms$arm[j] == "p")   hit_p <- TRUE
      if (arms$arm[j] == "q")   hit_q <- TRUE
      if (arms$arm[j] == "cen") hit_c <- TRUE
    }
  }
  
  if (hit_p & hit_c & hit_q) {
    "p_centromere_q"
  } else if (hit_p & hit_c) {
    "p_centromere"
  } else if (hit_c & hit_q) {
    "centromere_q"
  } else if (hit_p) {
    "p_arm"
  } else if (hit_q) {
    "q_arm"
  } else {
    NA_character_
  }
}

#' Classify CNV events by chromosome arm
#'
#' Applies arm classification to each CNV in a data frame.
#'
#' @param cnv_df A CNV data frame.
#' @param chromosome_arms A chromosome arm annotation data frame.
#'
#' @return The input CNV data frame with an added arm_class column.
classify_cnv_arms <- function(cnv_df, chromosome_arms) {
  
  required_cnv  <- c("chr", "start", "end")
  required_arms <- c("chr", "arm_start", "arm_end", "arm")
  
  missing_cnv  <- setdiff(required_cnv,  colnames(cnv_df))
  missing_arms <- setdiff(required_arms, colnames(chromosome_arms))
  
  if (length(missing_cnv)  > 0) stop("cnv_row missing columns: ",        paste(missing_cnv,  collapse = ", "))
  if (length(missing_arms) > 0) stop("chromosome_arms missing columns: ", paste(missing_arms, collapse = ", "))
  
  # Validate arm labels once here — not inside the per-row function
  valid_arms    <- c("p", "q", "cen")
  unexpected    <- setdiff(unique(chromosome_arms$arm), valid_arms)
  if (length(unexpected) > 0) {
    warning(
      "Unexpected arm labels in chromosome_arms: ",
      paste(unexpected, collapse = ", "),
      ". Expected: ", paste(valid_arms, collapse = ", ")
    )
  }
  
  if (nrow(cnv_df) == 0L) {
    warning("cnv_df is empty — returning with arm_class column set to NA.")
    cnv_df$arm_class <- NA_character_
    return(cnv_df)
  }
  
  arms_by_chr <- split(chromosome_arms, chromosome_arms$chr)
  
  cnv_df$arm_class <- vapply(
    seq_len(nrow(cnv_df)),
    function(i) {

      chr_arms <- arms_by_chr[[as.character(cnv_df[i,]$chr)]]
      if (is.null(chr_arms)) return(NA_character_)
      classify_single_cnv(cnv_df[i, ], chr_arms)
    },
    character(1)
  )

  na_rate <- mean(is.na(cnv_df$arm_class))
  if (na_rate > 0.1) {
    warning(sprintf(
      "%.1f%% of CNVs could not be arm-classified — check chromosome naming convention.",
      na_rate * 100
    ))}
    
  return(cnv_df)
}


#' Calculate chromosome and arm-level CNV coverage percentages
#'
#' Computes the percentage of chromosome-wide, p-arm, and q-arm span covered by
#' each CNV event.
#'
#' @param cnv_df A CNV data frame containing CNV coordinates and arm classes.
#' @param chromosome_arms A chromosome arm annotation table.
#' @return A data frame with additional percentage columns including the percentages for whole chrmossome, as well as p and q arm
calculate_cnv_arm_percentages <- function(cnv_df, chromosome_arms) {
  
  # --- Find chromosomes with CNVs ---
  chromosomes_with_cnv <- unique(cnv_df$chr)
  
  # Optional: check which chromosomes are missing
  all_chromosomes <- unique(chromosome_arms$chr)
  missing_chr <- setdiff(all_chromosomes, chromosomes_with_cnv)
  if(length(missing_chr) > 0){
    message("Skipping chromosomes with no CNV data: ", paste(missing_chr, collapse = ", "))
  }
  
  cnv_arm_percentage <- lapply(unique(cnv_df$chr), function(x){
    
    chr_subset <- cnv_df %>% filter(chr == x)
    
    chromosome_arms_subset <- chromosome_arms %>% filter(chr == x)
    
    whole_length <- chromosome_arms_subset %>%
      summarise(arm_length = max(arm_end) - min(arm_start)) %>%
      pull(arm_length)
    
    p_length <- chromosome_arms_subset %>%
      filter(arm == "p") %>%
      pull(arm_length)
    q_length <- chromosome_arms_subset %>%
      filter(arm == "q") %>%
      pull(arm_length)
    
    chr_subset %>%
      mutate(
        whole_chromosome_gain = ifelse(cnv_state == "gain",
                                       round(cnv_length / whole_length * 100,2),
                                       NA_real_),
        whole_chromosome_loss = ifelse(cnv_state == "loss",
                                       round(cnv_length / whole_length * 100,2),
                                       NA_real_),
        p_arm_gain = ifelse(cnv_state == "gain" & arm_class == "p_arm",
                            round(cnv_length / p_length*100,2),
                            NA_real_),
        p_arm_loss = ifelse(cnv_state == "loss" & arm_class == "p_arm",
                            round(cnv_length / p_length*100,2),
                            NA_real_),
        q_arm_gain = ifelse(cnv_state == "gain" & arm_class == "q_arm",
                            round(cnv_length / q_length*100,2),
                            NA_real_),
        q_arm_loss = ifelse(cnv_state == "loss" & arm_class == "q_arm",
                            round(cnv_length / q_length*100,2),
                            NA_real_)
      )
  })
  
  # Remove NULLs from skipped chromosomes
  cnv_total <- do.call(rbind, cnv_arm_percentage)
  
  return(cnv_total)
}



#' Add chromosome arm classification and percentages to a CNV table
#'
#' Validates coordinate columns and annotates CNVs with chromosome arm classes.
#'
#' @param main_df A CNV data frame.
#' @param chromosome_arms A chromosome arm annotation table.
#' @param chr_col Name of the chromosome column.
#' @param start_col Name of the interval start column.
#' @param end_col Name of the interval end column.
#'
#' @return The input data frame with arm classification added.
add_chromosome_info <- function(main_df,
                                chromosome_arms,
                                chr_col = "chr",
                                start_col = "start",
                                end_col = "end") {
  
  required_cols <- c(chr_col, start_col, end_col)
  
  if (!all(required_cols %in% colnames(main_df))) {
    stop("Missing required CNV coordinate columns.")
  }
  
  cnv_df <- classify_cnv_arms(main_df, chromosome_arms)
  
  missing_chr <- setdiff(cnv_df$chr,chromosome_arms$chr)
  
  if (any(missing_chr)) {
    message(message("Skipping chromosomes which are not present in cnv_df: ", paste(missing_chr, collapse = ", ")))
  }
  
  cnv_df <- cnv_df %>%
    filter(chr %in% chromosome_arms$chr)
  
  whole_chr_info <- calculate_cnv_arm_percentages(cnv_df, chromosome_arms)
  return(whole_chr_info) 
}




############################################################################################################-
############################# Across and Within Cell type Integration Approach##############################
############################################################################################################-

std_events <- function(tbl, dataset, mode){
  tbl %>%
    mutate(cell_type = dataset, mode = mode,
           ds_cell = paste(cell_type, cell_name, sep="|"))
}



compute_cell_sizes <- function(
    metadata,
    group_cols,
    cell_col = "cell_name"
) {
  
  missing_cols <- setdiff(c(group_cols, cell_col), colnames(metadata))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }
  
  metadata |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_total_cells = dplyr::n_distinct(.data[[cell_col]]),
      .groups       = "drop"
    )
}


run_full_cnv_pipeline <- function(
    
  # ---- Entry point -------------------------------------------------------
  start_from        = c("block1", "block2", "block3", "block4"),
  precomputed       = list(),
  save_intermediate = FALSE,
  outdir            = NULL,
  
  # ---- Block 1 -----------------------------------------------------------
  counts_mx         = NULL,
  metadata          = NULL,
  cell_type_col     = "cell_type",
  gene_order_file   = NULL,
  modes              = c("within", "across"),
  chr_exclude       = c("MT", "Y"),
  min_max_counts    = c(100, 1e6),
  n_splits_within   = 3,
  base_outdir       = NULL,
  cutoff            = 0.1,
  cluster_by_groups = TRUE,
  HMM               = FALSE,
  denoise           = TRUE,
  analysis_mode     = "subclusters",
  window_length     = 140,
  no_plot           = TRUE,
  resume_if_exists  = TRUE,
  
  # ---- Block 2 -----------------------------------------------------------
  base_dir                              = NULL,
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
  
  # ---- Block 3 -----------------------------------------------------------
  chromosome_arms   = NULL,
  group_cols        = NULL,
  cell_col          = "cell_name",
  chr_col           = "chr",
  start_col         = "start",
  end_col           = "end",
  
  # ---- Block 4 -----------------------------------------------------------
  by                          = NULL,
  sample_col                  = NULL,
  overlap_method              = "reciprocal",
  min_overlap                 = 0.8,
  boundaries_mb               = c(25, 10),
  base_fraction               = 0.05,
  step                        = 0.05,
  min_cap_threshold           = 5L,
  max_cap_threshold           = 25L,
  total_chromosome_permission = 65
) {
  
  start_from <- match.arg(start_from)
  modes       <- match.arg(modes)
  
  # ---- Validate save_intermediate -----------------------------------------
  if (save_intermediate && is.null(outdir)) {
    stop("outdir must be provided when save_intermediate = TRUE.")
  }
  if (save_intermediate && !dir.exists(outdir)) {
    message("Creating outdir: ", outdir)
    dir.create(outdir, recursive = TRUE)
  }
  
  # ---- Determine which blocks to run --------------------------------------
  valid_blocks  <- c("block1", "block2", "block3", "block4")
  blocks_to_run <- valid_blocks[
    which(valid_blocks == start_from):length(valid_blocks)
  ]
  
  if (start_from != "block1"  && length(precomputed) == 0L && start_from != "block2") {
    stop(
      "start_from = '", start_from, "' but precomputed is empty.\n",
      "Provide required inputs in precomputed list."
    )
  }
  
  # ---- Initialise ---------------------------------------------------------
  results  <- list(block1 = NULL, block2 = NULL,
                   block3 = NULL, block4 = NULL)
  summaries <- list()
  t_start   <- proc.time()
  
  message(sprintf(
    "=== CNV Pipeline | start_from = %s | tool = %s ===",
    start_from, tool
  ))
  
  # =========================================================================
  # BLOCK 1 — inferCNV creation and execution
  # =========================================================================
  if ("block1" %in% blocks_to_run) {
    
    if (is.null(counts_mx))       stop("counts_mx required for block1.")
    if (is.null(metadata))        stop("metadata required for block1.")
    if (is.null(gene_order_file)) stop("gene_order_file required for block1.")
    if (is.null(base_outdir))     stop("base_outdir required for block1.")
    
    message("\n[1/4] Running inferCNV pipeline...")
    t1 <- proc.time()
    
    results$block1 <- run_infercnv_pipeline(
      counts_mx         = counts_mx,
      metadata          = metadata,
      cell_type_col     = cell_type_col,
      gene_order_file   = gene_order_file,
      mode              = modes,
      chr_exclude       = chr_exclude,
      min_max_counts    = min_max_counts,
      n_splits_within   = n_splits_within,
      base_outdir       = base_outdir,
      cutoff            = cutoff,
      cluster_by_groups = cluster_by_groups,
      HMM               = HMM,
      denoise           = denoise,
      analysis_mode     = analysis_mode,
      window_length     = window_length,
      no_plot           = no_plot,
      resume_if_exists  = resume_if_exists
    )
    
    
    if (!is.null(results$block1$metadata)) {
      metadata <- results$block1$metadata
      message("metadata enriched with split_group column (within mode).")
    }
    
    # Pass base_outdir as base_dir for block2
    base_dir <- base_outdir
    
    summaries$block1 <- list(runtime_s = (proc.time() - t1)[["elapsed"]])
    
    if (save_intermediate) {
      path <- file.path(outdir, "block1_infercnv.rds")
      saveRDS(results$block1, path)
      message("  Saved: ", path)
    }
  }
  
  # =========================================================================
  # BLOCK 2 — Load CNV calls and extract supported events
  # =========================================================================
  if ("block2" %in% blocks_to_run) {
    
    # Re-entry: resolve base_dir from precomputed
    if (start_from == "block2") {
      if (is.null(base_dir)) stop("base_dir required for block2.")
      if (is.null(metadata)) stop("metadata required for block2.")
    }
    
    modes <- if(modes == "both") c("within", "across") else modes
    
    message("\n[2/4] Loading and processing CNV calls...")
    t2 <- proc.time()
    
    
    full_results <- process_tool_cnv_runs(
      base_dir                              = base_dir,
      mode                                  = modes,
      tool                                  = tool,
      pattern                               = pattern,
      max_gap                               = max_gap,
      min_overlap_consistent_calls          = min_overlap_consistent_calls,
      min_overlap_multiple_nodes            = min_overlap_multiple_nodes,
      filter_seq_mb_init                    = filter_seq_mb_init,
      filter_seq_mb_equiv                   = filter_seq_mb_equiv,
      min_references                        = min_references,
      overlap_method_equiv_cnv_call_merge   = overlap_method_equiv_cnv_call_merge,
      overlap_method_equiv_cnv_after_filter = overlap_method_equiv_cnv_after_filter,
      parallel                              = parallel,
      cores                                 = cores,
      clique_mode_consistent                = clique_mode_consistent,
      removed_log_return                    = removed_log_return,
      metadata                              = metadata
    )
    
    # Extract supported events — the key output for downstream blocks
    supported_events <- purrr::map(full_results, \(x) {
      x[["cnvs_supported_overlaped"]]
    }) |>
      purrr::compact() |>
      dplyr::bind_rows()
    
    if (nrow(supported_events) == 0L) {
      stop("No supported events after block2 — check pipeline parameters.")
    }
    
    results$block2 <- list(
      supported_events = supported_events,
      full_results     = full_results
    )
    
    summaries$block2 <- list(
      n_events  = nrow(supported_events),
      n_cells   = dplyr::n_distinct(supported_events$cell_name),
      runtime_s = (proc.time() - t2)[["elapsed"]]
    )
    
    message(sprintf(
      "  %d supported events across %d cells",
      summaries$block2$n_events,
      summaries$block2$n_cells
    ))
    
    if (save_intermediate) {
      path <- file.path(outdir, "block2_supported_events.rds")
      saveRDS(supported_events, path)
      message("  Saved: ", path)
    }
  }
  
  # =========================================================================
  # BLOCK 3 — Annotate CNV events + compute cell sizes
  # =========================================================================
  if ("block3" %in% blocks_to_run) {
    
    if (is.null(chromosome_arms)) stop("chromosome_arms required for block3.")
    if (is.null(group_cols))      stop("group_cols required for block3.")
    
    # Resolve supported_events
    supported_events <- if (start_from == "block3") {
      if (!"supported_events" %in% names(precomputed)) {
        stop(
          "precomputed$supported_events required when start_from = 'block3'.\n",
          "Can be a data frame or path to RDS file."
        )
      }
      se <- precomputed$supported_events
      if (is.character(se)) {
        if (!file.exists(se)) stop("RDS path not found: ", se)
        message("  Loading supported_events from: ", se)
        readRDS(se)
      } else {
        se
      }
    } else {
      results$block2$supported_events
    }
    
    if (is.null(metadata)) stop("metadata required for block3.")
    
    message("\n[3/4] Annotating CNV events...")
    t3 <- proc.time()
    
    cnv_annotated <- add_chromosome_info(
      supported_events,
      chromosome_arms,
      chr_col   = chr_col,
      start_col = start_col,
      end_col   = end_col
    )
    
    cell_sizes <- compute_cell_sizes(
      metadata   = metadata,
      group_cols = group_cols,
      cell_col   = cell_col
    )
    
    results$block3 <- list(
      cnv_annotated = cnv_annotated,
      cell_sizes    = cell_sizes
    )
    
    summaries$block3 <- list(
      n_annotated      = nrow(cnv_annotated),
      n_groups         = nrow(cell_sizes),
      cell_size_range  = range(cell_sizes$n_total_cells),
      runtime_s        = (proc.time() - t3)[["elapsed"]]
    )
    
    message(sprintf(
      "  %d annotated events | %d groups | cell range: %d-%d",
      summaries$block3$n_annotated,
      summaries$block3$n_groups,
      summaries$block3$cell_size_range[1],
      summaries$block3$cell_size_range[2]
    ))
    
    if (save_intermediate) {
      path <- file.path(outdir, "block3_annotated.rds")
      saveRDS(results$block3, path)
      message("  Saved: ", path)
    }
  }
  
  # =========================================================================
  # BLOCK 4 — Cluster and score CNV loci
  # =========================================================================
  if ("block4" %in% blocks_to_run) {
    
    # Resolve inputs
    if (start_from == "block4") {
      if (!all(c("cnv_annotated", "cell_sizes") %in% names(precomputed))) {
        stop(
          "precomputed must contain 'cnv_annotated' and 'cell_sizes' ",
          "when start_from = 'block4'."
        )
      }
      cnv_annotated <- precomputed$cnv_annotated
      cell_sizes    <- precomputed$cell_sizes
    } else {
      cnv_annotated <- results$block3$cnv_annotated
      cell_sizes    <- results$block3$cell_sizes
    }
    
    # Default by and sample_col to group_cols if not specified
    if (is.null(by)) {
      by <- group_cols
      message("  by defaulting to group_cols: ",
              paste(group_cols, collapse = ", "))
    }
    if (is.null(sample_col)) {
      sample_col <- group_cols[1]
      message("  sample_col defaulting to: ", sample_col)
    }
    
    message("\n[4/4] Clustering and scoring CNV loci...")
    t4 <- proc.time()
    
    clustered_events <- run_cnv_locus_analysis(
      cnv_annotated,
      by             = by,
      overlap_method = overlap_method,
      min_ovelap     = min_overlap,
      sample_col     = sample_col,
      cell_col       = cell_col
    )
    
    scored_events <- score_cnv_clusters(
      summary_df               = clustered_events$cnv_locus_summary,
      clustered_events         = clustered_events$clustered_events,
      cell_sizes               = cell_sizes,
      by_union                 = sample_col,
      boundaries_mb            = boundaries_mb,
      base_fraction            = base_fraction,
      step                     = step,
      threshold_method         = "auto",
      threshold_mode           = "fractions",
      min_cap_threshold        = min_cap_threshold,
      max_cap_threshold        = max_cap_threshold,
      total_chromosome_permission = total_chromosome_permission,
      round_fun                = ceiling
    )
    
    results$block4 <- list(
      clustered_events = clustered_events,
      scored_events    = scored_events
    )
    
    summaries$block4 <- list(
      n_clustered = nrow(clustered_events$clustered_events),
      n_scored    = nrow(scored_events),
      runtime_s   = (proc.time() - t4)[["elapsed"]]
    )
    
    message(sprintf(
      "  %d clustered events | %d scored loci",
      summaries$block4$n_clustered,
      summaries$block4$n_scored
    ))
    
    if (save_intermediate) {
      path <- file.path(outdir, "block4_clustered_scored.rds")
      saveRDS(results$block4, path)
      message("  Saved: ", path)
    }
  }
  
  # ---- Pipeline summary ---------------------------------------------------
  total_runtime <- (proc.time() - t_start)[["elapsed"]]
  
  message(sprintf(paste0(
    "\n=== Pipeline complete ===\n",
    "  Blocks run:    %s\n",
    "  Total runtime: %.1f seconds (%.1f minutes)"
  ),
  paste(blocks_to_run, collapse = " → "),
  total_runtime,
  total_runtime / 60
  ))
  
  c(
    results,
    list(summary = list(
      blocks_run    = blocks_to_run,
      per_block     = summaries,
      total_runtime = total_runtime
    ))
  )
}









