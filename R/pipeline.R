#' @title pipeline
#'
#' @description
#' 
#' Functions that centralise all steps of the reference-independent
#' single-cell CNV karyotyping pipeline. The pipeline is organised
#' into four sequential parts that represent the 4 major processing steps which the current pipeline performs
#'
#' @author Pedro Granjo
#' @date 02-09-2026
#' 
#' 

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



load_tool_data <- function(
    tool    = c("infercnv", "scevan", "copykat"),
    base_dir,
    ref_dirs,
    pattern,
    k
) {
  tool <- match.arg(tool)
  
  switch(tool,
         "infercnv" = load_infercnv_data(base_dir,  ref_dirs = ref_dirs, pattern = pattern, k = k),
         #"scevan"   = load_scevan_data(base_dir, ref_dirs, pattern),
         #"copykat"  = load_copykat_data(base_dir, ref_dirs, pattern)
  )
  
}



#' Execute CNV inference tool
#'
#' Dispatches to the correct tool-specific pipeline based on the tool parameter.
#' Currently supports inferCNV. SCEVAN and CopyKat are planned for future integration.
#'
#' @param tool Character. CNV tool to run. Currently only "infercnv" supported.
#' @inheritParams run_infercnv_pipeline
#'
#' @return Named list with tool output. Structure depends on tool:
#'   inferCNV: list(obj_list, run_log, metadata, runtime)
run_cnv_tool <- function(
    tool              = c("infercnv", "scevan", "copykat"),
    counts_mx,
    metadata,
    group_col     = "cell_type",
    gene_order_file   = NULL,
    chromosomes_to_exclude       = c("MT", "Y"),
    min_max_counts    = c(100, 1e6),
    n_splits_within   = 3,
    base_outdir,
    cutoff            = 0.1,
    cluster_by_groups = TRUE,
    denoise           = TRUE,
    analysis_mode     = "subclusters",
    window_length     = 140L,
    no_plot           = TRUE,
    resume_if_exists  = TRUE,
    clonal_col = NULL,
    donor_col = NULL
) {
  
  tool <- match.arg(tool)
  
  switch(tool,
         
         "infercnv" = run_infercnv_pipeline(
           counts_mx         = counts_mx,
           metadata          = metadata,
           group_col     = group_col,
           gene_order_file   = gene_order_file,
           mode              = "within",
           chromosomes_to_exclude       = chromosomes_to_exclude,
           min_max_counts    = min_max_counts,
           n_splits_within   = n_splits_within,
           base_outdir       = base_outdir,
           cutoff            = cutoff,
           cluster_by_groups = cluster_by_groups,
           analysis_mode     = analysis_mode,
           window_length     = window_length,
           no_plot           = no_plot,
           resume_if_exists  = resume_if_exists,
           clonal_col = clonal_col,
           donor_col = donor_col
         ),
         
         "scevan" = stop(
           "SCEVAN integration not yet implemented. ",
           "Use tool = 'infercnv' for now."
         ),
         
         "copykat" = stop(
           "CopyKat integration not yet implemented. ",
           "Use tool = 'infercnv' for now."
         )
  )
}




process_tool_cnv_runs <- function(
    base_dir,
    mode                                  = "within",
    tool                                  = "infercnv",
    pattern                               = NULL,
    max_gap                               = 100000,
    min_overlap_consistent_calls          = 0.5,
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
    metadata,
    remove_ref                            = TRUE,
    k                                     = 1.5,
    # ── Gene density params ───────────────────────────────────────────────
    coding_gr            = NULL,
    coding_expressed_set = NULL,
    pct_max              = 45,
    pct_floor            = 30,
    min_expr_density     = 1.5,
    min_coding_density   = 1.0,
    max_gap_mb           = 10
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
  group_clusters <- list.files(cell_dir)
  
  if (length(group_clusters) == 0L) {
    stop(sprintf(
      paste0(
        "No cell types found in '%s' for mode '%s'.\n",
        "Check that inferCNV has been run and results exist:\n",
        "  %s/{cell_group}/{A,B,C}/run.final.infercnv_obj"
      ),
      cell_dir, mode, cell_dir
    ))
  }
  
  ct_results <- purrr::map(group_clusters, \(ct) {
    
    ct_dir   <- file.path(cell_dir, ct)
    refs_dir <- list.files(ct_dir, full.names = T)
    refs_dir <- list.files(ct_dir)[dir.exists(file.path(refs_dir))]
    
    if (length(refs_dir) == 0L) {
      warning(sprintf("No reference directories found in %s — skipping.", ct_dir))
      return(NULL)
    }
    
    message(sprintf("Processing mode='%s', ='%s', refs=%s",
                    mode, ct, paste(refs_dir, collapse = ", ")))
    
    # Load tool data
    tools_data <- load_tool_data(
      tool     = tool,
      base_dir = ct_dir,
      ref_dirs = refs_dir,
      pattern  = pattern,
      k = k
    )
    
    # Run CNV pipeline
    run_fast_cnv_pipeline(
      gene_level_df                         = tools_data[[1]],
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
      mode = mode,
      remove_ref = remove_ref,
      # ── Gene density params ──────────────────────────────────────────
      coding_gr                 = coding_gr,
      min_expr_density          = min_expr_density,
      gene_order           = tools_data[[2]],
      coding_expressed_set = coding_expressed_set,
      pct_max              = pct_max,
      pct_floor            = pct_floor,
      min_coding_density   = min_coding_density,
      max_gap_mb           = max_gap_mb
      )
     
    
  })
  
  names(ct_results) <- group_clusters
  ct_results
}


#' Resolve metadata for Block2 in single execution mode
#'
#' In single mode, one split_metadata.rds exists at the within/ level.
#' If metadata already in memory — use directly.
#' Otherwise search base_dir/within/ for split_metadata.rds.
#'
#' @param metadata  data.frame or NULL
#' @param base_dir  path to inferCNV output base directory
#' @param outdir    path to pipeline output directory
#'
#' @return validated metadata data.frame
.resolve_metadata_single <- function(metadata,
                                     base_dir,
                                     outdir,
                                     group_col = "cell_type") {
  
  if (!is.null(metadata) && "split_group" %in% colnames(metadata)) {
    # Already in memory with split_group — use directly
    message("  Using provided metadata (single mode).")
    return(metadata)
  }
  
  if (!is.null(metadata) && !"split_group" %in% colnames(metadata)) {
    message(
      "  Provided metadata missing 'split_group'.\n",
      "  Searching for split_metadata.rds saved by Block1..."
    )
    metadata <- NULL
  }
  
  # Search outdir first, fall back to base_dir/within/
  search_paths <- c(
    file.path(outdir,    "split_metadata.rds"),
    file.path(base_dir,  "within", "split_metadata.rds")
  )
  
  search_paths <- search_paths[!sapply(search_paths, is.null)]
  
  meta_path <- NULL
  for (path in search_paths) {
    if (file.exists(path)) {
      meta_path <- path
      break
    }
  }
  
  if (is.null(meta_path)) {
    stop(
      "No split_metadata.rds found in single mode.\n",
      "Searched:\n",
      paste(" ", search_paths, collapse = "\n"), "\n",
      "Either:\n",
      "  1. Provide --metadata_path with split_group column\n",
      "  2. Re-run Block1 (single mode) to generate split_metadata.rds"
    )
  }
  
  message("  Loading split_metadata from: ", meta_path)
  metadata <- readRDS(meta_path)
  
  # Validate
  if (!"split_group" %in% colnames(metadata)) {
    stop(
      "Loaded split_metadata missing 'split_group' column.\n",
      "Re-run Block1 to regenerate it."
    )
  }
  
  n_na <- sum(is.na(metadata$split_group))
  if (n_na > 0L) {
    stop(sprintf(
      "%d cells have NA split_group — re-run Block1.",
      n_na
    ))
  }
  
  message(sprintf(
    "  Resolved metadata: %d cells, %d cell types, groups: %s",
    nrow(metadata),
    dplyr::n_distinct(metadata[[group_col]]),
    paste(sort(unique(metadata$split_group)), collapse = ", ")
  ))
  
  return(metadata)
}


#' Resolve metadata for Block2 in array execution mode
#'
#' In array mode, one split_metadata.rds exists per cell type directory.
#' Auto-discovers all files recursively and combines via bind_rows.
#' Validates no cell overlap between files.
#'
#' @param metadata  data.frame or NULL
#' @param base_dir  path to inferCNV output base directory
#' @param outdir    path to pipeline output directory
#'
#' @return validated combined metadata data.frame
.resolve_metadata_array <- function(metadata,
                                    base_dir,
                                    outdir,
                                    group_col = "cell_type") {
  
  if (!is.null(metadata) && "split_group" %in% colnames(metadata)) {
    message("  Using provided metadata (array mode).")
    return(metadata)
  }
  
  if (!is.null(metadata) && !"split_group" %in% colnames(metadata)) {
    message(
      "  Provided metadata missing 'split_group'.\n",
      "  Searching for split_metadata.rds files saved by Block1 array..."
    )
    metadata <- NULL
  }
  
  # Search outdir first, fall back to base_dir
  # Array mode: multiple files expected, one per cell type
  search_dirs <- c(outdir, base_dir)
  search_dirs <- search_dirs[!sapply(search_dirs, is.null)]
  meta_files  <- character(0)
  
  for (search_dir in search_dirs) {
    meta_files <- list.files(
      search_dir,
      pattern    = "split_metadata.rds",
      recursive  = TRUE,
      full.names = TRUE
    )
    if (length(meta_files) > 0L) {
      message(sprintf(
        "  Found %d split_metadata.rds file(s) in: %s",
        length(meta_files), search_dir
      ))
      break
    }
  }
  
  if (length(meta_files) == 0L) {
    stop(
      "No split_metadata.rds files found in array mode.\n",
      "Searched in:\n",
      paste(" ", search_dirs, collapse = "\n"), "\n",
      "Either:\n",
      "  1. Provide --metadata_path with split_group column\n",
      "  2. Re-run Block1 array job to generate per-cell-type files"
    )
  }
  
  purrr::walk(meta_files, ~ message("  Loading: ", .x))
  meta_list <- purrr::map(meta_files, readRDS)
  
  # ---- Validate no cell overlap between files -----------------------------
  all_cells <- unlist(purrr::map(meta_list, ~ .x$cell_name))
  dup_cells <- all_cells[duplicated(all_cells)]
  
  if (length(dup_cells) > 0L) {
    stop(sprintf(
      "%d cell(s) appear in multiple split_metadata files:\n%s\n",
      "Check that Block1 array tasks processed distinct cell types.",
      length(dup_cells),
      paste(head(dup_cells, 10), collapse = ", ")
    ))
  }
  
  metadata <- dplyr::bind_rows(meta_list)
  
  # ---- Validate after combining -------------------------------------------
  if (!"split_group" %in% colnames(metadata)) {
    stop(
      "Combined split_metadata missing 'split_group' column.\n",
      "Re-run Block1 array to regenerate files."
    )
  }
  
  n_na <- sum(is.na(metadata$split_group))
  if (n_na > 0L) {
    stop(sprintf(
      "%d cells have NA split_group — check Block1 output for: %s",
      n_na,
      paste(unique(mmetadata[[group_col]][is.na(metadata$split_group)]),
            collapse = ", ")
    ))
  }
  
  message(sprintf(
    "  Combined metadata: %d cells, %d cell types, groups: %s",
    nrow(metadata),
    dplyr::n_distinct(metadata[[group_col]]),
    paste(sort(unique(metadata$split_group)), collapse = ", ")
  ))
  
  return(metadata)
}

#' Resolve and validate metadata for Block2
#'
#' Dispatches to the correct resolution strategy based on execution_mode.
#' After resolution validates metadata against disk structure.
#'
#' @param metadata        data.frame or NULL
#' @param base_dir        path to inferCNV output base directory
#' @param outdir          path to pipeline output directory
#' @param execution_mode  "single" or "array"
#' @param pattern         regex to match inferCNV result files
#'
#' @return validated metadata data.frame with split_group column
resolve_block2_metadata <- function(metadata,
                                    base_dir,
                                    outdir,
                                    execution_mode = c("single", "array"),
                                    pattern        = "^run\\.final",
                                    group_col = "cell_type") {
  
  execution_mode <- match.arg(execution_mode)
  
  message(sprintf(
    "\n── Resolving Block2 metadata (execution_mode = '%s') ───────────",
    execution_mode
  ))
  
  # ---- Dispatch to correct resolver ---------------------------------------
  metadata <- switch(
    execution_mode,
    "single" = .resolve_metadata_single(
      metadata = metadata,
      base_dir = base_dir,
      outdir   = outdir,
      group_col = group_col
    ),
    "array"  = .resolve_metadata_array(
      metadata = metadata,
      base_dir = base_dir,
      outdir   = outdir,
      group_col = group_col
    )
  )
  
  # ---- Validate against disk — always, regardless of mode ----------------
  validate_metadata_against_disk(
    metadata = metadata,
    base_dir = base_dir,
    pattern  = pattern,
    group_col = group_col
  )
  
  return(metadata)
}


#' Validate metadata against existing inferCNV output directories
#'
#' Checks that:
#' 1. split_group column exists
#' 2. No NA in split_group
#' 3. All cell types in metadata have a corresponding directory
#' 4. All reference directories (A/B/C) exist per cell type
#' 5. All cells in metadata are accounted for across cell types
#'
#' @param metadata Data frame with split_group and cell_type columns
#' @param base_dir Path to inferCNV output base directory
#' @param mode     "within" or "across"
#' @param pattern  Pattern to match inferCNV result files
#'
#' @return invisibly returns TRUE if all checks pass
validate_metadata_against_disk <- function(
    metadata,
    base_dir,
    mode          = "within",
    pattern       = "^run\\.final",
    group_col     = "cell_type", 
    cell_name_col = "cell_name" 
) {
  
  mode_dir <- file.path(base_dir, mode)
  issues   <- list()
  
  # ---- Check 1 — split_group column exists --------------------------------
  if (!"split_group" %in% colnames(metadata)) {
    stop(
      "metadata is missing 'split_group' column.\n",
      "Provide split_metadata.rds saved by Block 1."
    )
  }

  # ---- Check 2 — no NA in split_group ------------------------------------
  na_cells <- metadata[[cell_name_col]][is.na(metadata$split_group)]
  if (length(na_cells) > 0L) {
    issues <- c(issues, sprintf(
      "%d cells have NA split_group: %s",
      length(na_cells),
      paste(head(na_cells, 5), collapse = ", ")
    ))
  }
  
  # ---- Check 3 — group_col column exists ----------------------------------
  if (!group_col %in% colnames(metadata)) {
    stop(sprintf(
      "metadata is missing '%s' column.\n",
      "Set group_col to the correct column name.",
      group_col
    ))
  }
  
  # ---- Check 4 — mode directory exists ------------------------------------
  if (!dir.exists(mode_dir)) {
    stop(sprintf(
      "Mode directory not found: %s\nRun Block 1 first.",
      mode_dir
    ))
  }
  
  # ---- Check 5 — group directories exist on disk -------------------------
  group_values <- unique(metadata[[group_col]])
  group_dirs   <- file.path(mode_dir, group_values)
  missing_dirs <- group_values[!dir.exists(group_dirs)]
  
  if (length(missing_dirs) > 0L) {
    issues <- c(issues, sprintf(
      "Group directories missing in %s: %s",
      mode_dir,
      paste(missing_dirs, collapse = ", ")
    ))
  }
  
  if (length(missing_dirs) == length(group_values)) {
    stop(sprintf(
      "No group directories found in %s.\nExpected: %s\nRun Block 1 first.",
      mode_dir,
      paste(group_values, collapse = ", ")
    ))
  }
  
  # ---- Check 6 — reference directories exist per group -------------------
  present_groups <- group_values[dir.exists(group_dirs)]
  
  for (grp in present_groups) {
    
    grp_dir  <- file.path(mode_dir, grp)
    ref_dirs <- list.files(grp_dir, full.names = FALSE)
    ref_dirs <- ref_dirs[dir.exists(file.path(grp_dir, ref_dirs))]
    
    if (length(ref_dirs) == 0L) {
      issues <- c(issues, sprintf(
        "No reference directories found for group '%s' in %s",
        grp, grp_dir
      ))
      next
    }
    
    # ---- Check 7 — inferCNV result files exist per reference ------------
    for (ref in ref_dirs) {
      ref_path     <- file.path(grp_dir, ref)
      result_files <- list.files(ref_path,
                                 pattern    = pattern,
                                 full.names = FALSE)
      if (length(result_files) == 0L) {
        issues <- c(issues, sprintf(
          "No result files matching '%s' found for %s/%s",
          pattern, grp, ref
        ))
      }
    }
    
    # ---- Check 8 — split_group values match reference directories -------
    groups_meta <- sort(unique(
      metadata$split_group[metadata[[group_col]] == grp]
    ))
    groups_disk  <- sort(ref_dirs)
    missing_refs <- setdiff(groups_meta, groups_disk)
    extra_refs   <- setdiff(groups_disk, groups_meta)
    
    if (length(missing_refs) > 0L) {
      issues <- c(issues, sprintf(
        "Group '%s': split_group values %s in metadata but no directory on disk",
        grp, paste(missing_refs, collapse = ", ")
      ))
    }
    if (length(extra_refs) > 0L) {
      issues <- c(issues, sprintf(
        "Group '%s': directories %s on disk but not in metadata split_group",
        grp, paste(extra_refs, collapse = ", ")
      ))
    }
  }
  
  # ---- Check 9 — all cells accounted for ---------------------------------
  if (length(present_groups) > 0L) {
    n_cells_meta <- nrow(metadata)
    n_cells_disk <- sum(sapply(present_groups, function(grp) {
      nrow(metadata[metadata[[group_col]] == grp, ])
    }))
    
    if (n_cells_meta != n_cells_disk) {
      issues <- c(issues, sprintf(
        "Cell count mismatch: %d cells in metadata but %d accounted for across groups",
        n_cells_meta, n_cells_disk
      ))
    }
  }
  
  # ---- Report all issues at once -----------------------------------------
  if (length(issues) > 0L) {
    stop(sprintf(
      "%d validation issue(s) found:\n%s",
      length(issues),
      paste(seq_along(issues), issues,
            sep = ". ", collapse = "\n")
    ))
  }
  
  # ---- Summary on success -------------------------------------------------
  message(sprintf(paste0(
    "Metadata validation passed:\n",
    "  Cells:        %d\n",
    "  Groups (%s):  %d (%s)\n",
    "  Split groups: %s"
  ),
  nrow(metadata),
  group_col,
  dplyr::n_distinct(metadata[[group_col]]),
  paste(sort(unique(metadata[[group_col]])), collapse = ", "),
  paste(sort(unique(metadata$split_group)), collapse = ", ")
  ))
  
  invisible(TRUE)
}



#' Run Block 1 for all cell types in a single job
#'
#' @param counts_mx       genes x cells count matrix
#' @param metadata        full metadata data.frame
#' @param group_col   column name for cell types
#' @param gene_order_file path to gene order file
#' @param base_outdir     base output directory
#' @param tool            CNV tool to use
#' @param chromosomes_to_exclude     chromosomes to exclude
#' @param min_max_counts  c(min, max) counts per cell
#' @param n_splits_within number of reference splits
#' @param cutoff          inferCNV cutoff
#' @param window_length   smoothing window
#' @param no_plot         skip plots
#' @param resume_if_exists skip completed runs
#'
#' @return list with obj_list and split_metadata
.run_block1_single <- function(
    counts_mx,
    metadata,
    group_col     = "cell_type",
    gene_order_file,
    base_outdir,
    tool              = "infercnv",
    chromosomes_to_exclude       = c("MT", "Y"),
    min_max_counts    = c(100, 1e6),
    n_splits_within   = 3L,
    cutoff            = 0.1,
    cluster_by_groups = TRUE,
    analysis_mode     = "subclusters",
    window_length     = 140L,
    no_plot           = TRUE,
    resume_if_exists  = TRUE,
    clonal_col = NULL,
    donor_col = NULL
) {
  
  message("\n── Block1 single mode — all cell types ──────────────────────────")
  
  # ---- Run all cell types together ----------------------------------------
  result <- run_cnv_tool(
    tool              = tool,
    counts_mx         = counts_mx,
    metadata          = metadata,
    group_col     = group_col,
    gene_order_file   = gene_order_file,
    chromosomes_to_exclude       = chromosomes_to_exclude,
    min_max_counts    = min_max_counts,
    n_splits_within   = n_splits_within,
    base_outdir       = base_outdir,
    cutoff            = cutoff,
    cluster_by_groups = cluster_by_groups,
    analysis_mode     = analysis_mode,
    window_length     = window_length,
    no_plot           = no_plot,
    resume_if_exists  = resume_if_exists,
    clonal_col = clonal_col,
    donor_col = donor_col
    
  )
  
  # ---- Validate and extract split_metadata --------------------------------
  split_metadata <- result$metadata
  
  if (is.null(split_metadata)) {
    stop(
      "split_metadata not returned by run_cnv_tool().\n",
      "Block2 requires this — check run_cnv_tool() output."
    )
  }
  
  if (!"split_group" %in% colnames(split_metadata)) {
    stop(
      "split_metadata missing 'split_group' column.\n",
      "Check run_cnv_tool() output."
    )
  }
  
  message(sprintf(
    "Block1 complete: %d cells, %d cell types, groups: %s",
    nrow(split_metadata),
    dplyr::n_distinct(split_metadata[[group_col]]),
    paste(sort(unique(split_metadata$split_group)), collapse = ", ")
  ))
  
  # ---- Save split_metadata at within/ level -------------------------------
  # Single mode saves one file at base level
  # Block2 single mode loads from here
  within_dir <- file.path(base_outdir, "within")
  dir.create(within_dir, recursive = TRUE, showWarnings = FALSE)
  
  meta_path <- file.path(within_dir, "split_metadata.rds")
  saveRDS(split_metadata, meta_path)
  message("split_metadata saved to: ", meta_path)
  
  return(result)
}




#' Run Block 1 for one cell type (array job task)
#'
#' Designed to be called once per array task.
#' Subsets metadata to the specified cell type before running.
#'
#' @param cell_group       cell type to process for this task
#' @param counts_mx       genes x cells count matrix
#' @param metadata        full metadata data.frame (all cell types)
#' @param group_col   column name for cell types
#' @param gene_order_file path to gene order file
#' @param base_outdir     base output directory
#' @param tool            CNV tool to use
#' @param chromosomes_to_exclude     chromosomes to exclude
#' @param min_max_counts  c(min, max) counts per cell
#' @param n_splits_within number of reference splits
#' @param cutoff          inferCNV cutoff
#' @param window_length   smoothing window
#' @param no_plot         skip plots
#' @param resume_if_exists skip completed runs
#'
#' @return list with obj_list and split_metadata for this cell type
.run_block1_array <- function(
    cell_group,
    counts_mx,
    metadata,
    group_col     = "cell_type",
    gene_order_file,
    base_outdir,
    tool              = "infercnv",
    chromosomes_to_exclude       = c("MT", "Y"),
    min_max_counts    = c(100, 1e6),
    n_splits_within   = 3L,
    cutoff            = 0.1,
    cluster_by_groups = TRUE,
    analysis_mode     = "subclusters",
    window_length     = 140L,
    no_plot           = TRUE,
    resume_if_exists  = TRUE
) {
  
  message(sprintf(
    "\n── Block1 array mode — cell type: %s ───────────────────────────",
    cell_group
  ))
  
  # ---- Validate cell_group -------------------------------------------------
  if (length(cell_group) != 1L) {
    stop(
      "cell_group must be a single value in array mode.\n",
      "For multiple cell types use execution_mode = 'single'."
    )
  }
  
  available_types <- unique(metadata[[group_col]])
  
  if (!cell_group %in% available_types) {
    stop(sprintf(
      "Cell Cluster '%s' not found in metadata column '%s'.\nAvailable: %s",
      cell_group,
      group_col,
      paste(available_types, collapse = ", ")
    ))
  }
  
  # ---- Subset metadata to this cell type ----------------------------------
  metadata_ct <- metadata[metadata[[group_col]] %in% cell_group, ,
                          drop = FALSE]
  
  if (nrow(metadata_ct) == 0L) {
    stop(sprintf(
      "No cells found for cell type '%s' after subsetting.", cell_group
    ))
  }
  
  # ---- Validate all cells present in counts matrix ------------------------
  cells_missing <- setdiff(metadata_ct$cell_name, colnames(counts_mx))
  
  if (length(cells_missing) > 0L) {
    stop(sprintf(
      "%d cell(s) in metadata not found in counts matrix:\n%s",
      length(cells_missing),
      paste(head(cells_missing, 10), collapse = ", ")
    ))
  }
  
  message(sprintf(
    "  Subsetting to %d cells for cell group '%s'",
    nrow(metadata_ct), cell_group
  ))
  
  # ---- Output directory for this cell type --------------------------------
  outdir_ct <- file.path(base_outdir, "within", cell_group)
  dir.create(outdir_ct, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Run tool for this cell type ----------------------------------------
  result <- run_cnv_tool(
    tool              = tool,
    counts_mx         = counts_mx,
    metadata          = metadata_ct, 
    group_col     = group_col,
    gene_order_file   = gene_order_file,
    chromosomes_to_exclude       = chromosomes_to_exclude,
    min_max_counts    = min_max_counts,
    n_splits_within   = n_splits_within,
    base_outdir       = outdir_ct,
    cutoff            = cutoff,
    cluster_by_groups = cluster_by_groups,
    analysis_mode     = analysis_mode,
    window_length     = window_length,
    no_plot           = no_plot,
    resume_if_exists  = resume_if_exists,
    clonal_col = clonal_col,
    donor_col = donor_col
  )
  
  # ---- Validate and save split_metadata -----------------------------------
  split_metadata <- result$metadata
  
  if (is.null(split_metadata)) {
    stop(
      "split_metadata not returned by run_cnv_tool().\n",
      "Block2 requires this — check run_cnv_tool() output."
    )
  }
  
  # Validate all cells for this cell type are in split_metadata
  cells_missing_meta <- setdiff(
    metadata_ct$cell_name,
    split_metadata$cell_name
  )
  
  if (length(cells_missing_meta) > 0L) {
    warning(sprintf(
      "%d cell(s) missing from split_metadata (filtered by inferCNV):\n%s",
      length(cells_missing_meta),
      paste(head(cells_missing_meta, 10), collapse = ", ")
    ))
  }
  
  # Array mode saves per cell type directory
  # Block2 array mode discovers these via recursive list.files()
  meta_path <- file.path(outdir_ct, "split_metadata.rds")
  saveRDS(split_metadata, meta_path)
  
  message(sprintf(
    "  split_metadata saved: %d cells, groups: %s → %s",
    nrow(split_metadata),
    paste(sort(unique(split_metadata$split_group)), collapse = ", "),
    meta_path
  ))
  
  return(result)
}



run_full_cnv_pipeline <- function(
    # ---- Entry point -------------------------------------------------------
    start_from        = c("block1", "block2", "block3", "block4"),
    execution_mode    = c("single", "array"),
    save_intermediate = FALSE,
    workdir            = NULL,
    
    # ---- Block 1 -----------------------------------------------------------
    counts_mx         = NULL,
    metadata          = NULL,
    group_clusters_col     = "cell_type",
    gene_order_file   = NULL,
    chromosomes_to_exclude       = c("MT", "Y"),
    min_max_counts    = c(100, 1e6),
    n_splits_within   = 3,
    cutoff            = 0.1,
    resume_if_exists  = TRUE,
    tool_outdir = NULL,
    k_discrete = 1.5,
    remove_ref = T,
    clonal_col = NULL,
    donor_col = NULL,
    # ---- Block 2 -----------------------------------------------------------
    tool                                  = "infercnv",
    pattern                               = "^run\\.final",
    max_gap                               = 100000,
    min_overlap_consistent_calls          = 0.75,
    min_overlap_multiple_nodes            = 0.6,
    filter_seq_mb_init                    = 5,
    filter_seq_mb_equiv                   = 7,
    min_references                        = 2,
    parallel                              = FALSE,
    cores                                 = 1L,
    clique_mode_consistent                = "connected",
    removed_log_return                    = FALSE,
    min_coding_density   = 1.0,
    max_gap_mb           = 10,
    
    # ---- Block 3 -----------------------------------------------------------
    supported_events  = NULL,
    chromosome_arms   = NULL,
    group_cols        = NULL,
    cell_col          = "cell_name",
    chr_col           = "chr",
    start_col         = "start",
    end_col           = "end",
    
    
    # ---- Block 4 -----------------------------------------------------------
    cnv_annotated = NULL,
    cell_sizes = NULL,
    by                          = NULL,
    sample_col                  = NULL,
    overlap_method              = "reciprocal",
    min_overlap                 = 0.8,
    min_required_cells           = 5L,
    k_threshold_growth = 1.5,
    sensitivity_floor_mb = 20,
    p_arm_permission     = 60,
    q_arm_permission     = 60,  
    whole_chr_permission = 65,
    coding_gr                 = NULL,
    coding_expressed_set      = NULL,
    pct_max                   = 45,
    pct_floor                 = 30,
    min_expr_density          = 1.7
) {
  
  start_from     <- match.arg(start_from)
  execution_mode <- match.arg(execution_mode)
  
  # ---- Validate save_intermediate -----------------------------------------
  if (save_intermediate && is.null(workdir)) {
    stop("workdir must be provided when save_intermediate = TRUE.")
  }
  if (save_intermediate && !dir.exists(workdir)) {
    message("Creating Working Directory: ", workdir)
    dir.create(workdir, recursive = TRUE)
  }
  
  # ---- Determine which blocks to run --------------------------------------
  valid_blocks  <- c("block1", "block2", "block3", "block4")
  blocks_to_run <- valid_blocks[
    which(valid_blocks == start_from):length(valid_blocks)
  ]
  
  # ---- Initialise ---------------------------------------------------------
  results   <- list(block1 = NULL, block2 = NULL,
                    block3 = NULL, block4 = NULL)
  summaries <- list()
  t_start   <- proc.time()
  
  message(sprintf(
    "=== CNV Pipeline | start_from = %s | execution_mode = %s | tool = %s ===",
    start_from, execution_mode, tool
  ))
  
  # =========================================================================
  # BLOCK 1 — CNV tool creation and execution (single mode only)
  # =========================================================================
  if ("block1" %in% blocks_to_run) {
    
    if (is.null(counts_mx))       stop("counts_mx required for block1.")
    if (is.null(metadata))        stop("metadata required for block1.")
    if (is.null(gene_order_file)) stop("gene_order_file required for block1.")
    if (is.null(workdir))     stop("workdir required for block1.")
    if(is.null(tool_outdir)) stop("tool_outdir required for block1.")
    
    message("\n[1/4] Running Block1 (single mode — all cell types)...")
    t1 <- proc.time()
    
    results$block1 <- .run_block1_single(
      counts_mx         = counts_mx,
      metadata          = metadata,
      group_col     = group_clusters_col,
      gene_order_file   = gene_order_file,
      base_outdir       = tool_outdir,
      tool              = tool,
      chromosomes_to_exclude       = chromosomes_to_exclude,
      min_max_counts    = min_max_counts,
      n_splits_within   = n_splits_within,
      cutoff            = cutoff,
      cluster_by_groups = TRUE,
      analysis_mode     = "subclusters",
      no_plot           = TRUE,
      resume_if_exists  = resume_if_exists,
      clonal_col = clonal_col,
      donor_col = donor_col
    )
    
    
    summaries$block1 <- list(
      runtime_s = (proc.time() - t1)[["elapsed"]]
    )
    
    message(sprintf(
      "  Block1 complete in %.1f minutes",
      summaries$block1$runtime_s / 60
    ))
    
    if (save_intermediate) {
      path <- file.path(workdir, "block1_results.rds")
      saveRDS(results$block1, path)
      message("  Saved: ", path)
    }
  }
  # =========================================================================
  # BLOCK 2 — Load CNV calls and extract supported events
  # =========================================================================
  if ("block2" %in% blocks_to_run) {
    
    message("\n[2/4] Loading and processing CNV calls...")
    t2 <- proc.time()
    
    # ---- Determine Block2 execution mode ----------------------------------
    # complete  → Block1 just ran, metadata in memory
    # single    → Block1 ran as single job previously, one metadata file
    # array     → Block1 ran as array job, N metadata files per cell type
    
    execution_mode_b2 <- if (!is.null(results$block1)) {
      "complete"
    } else {
      execution_mode
    }
    
    message(sprintf("  Block2 mode: %s", execution_mode_b2))
    
    # ---- Metadata resolution ----------------------------------------------
    if (execution_mode_b2 == "complete") {
      
      message("  Using metadata from Block1 (in memory).")
      metadata <- results$block1$metadata
      
    } else {
      
      # Fragmented — Block1 ran separately, load metadata from disk
      if (is.null(workdir)) {
        stop("workdir required when start_from = 'block2'.")
      }
      
      metadata <- resolve_block2_metadata(
        metadata       = metadata,
        base_dir       = tool_outdir,
        outdir         = workdir,
        execution_mode = execution_mode_b2,
        pattern        = pattern,
        group_col = group_clusters_col
      )
    }
    
    # ---- Execute Block2 ---------------------------------------------------
    # Two paths depending on whether Block1 objects are in memory
    
    if (!is.null(results$block1) &&
        !is.null(results$block1$obj_list)) {
      
      within_cell_group <- results$block1$obj_list[["within_cell_group"]]
      within_objs <- within_cell_group[["objects"]] %||% within_cell_group
      within_objs <- Filter(function(x) !is.null(x) && length(x) > 0L, within_objs)
      
    if (length(within_objs) == 0L) {
    stop(
      "No valid cell types found in Block1 output.\n",
      "All cell types were skipped — check min_cells threshold."
    )
  }
  
    message(sprintf(
    "  Processing %d cell type(s) from memory: %s",
    length(within_objs),
    paste(names(within_objs), collapse = ", ")
    ))
  
      
     full_results <- purrr::imap(within_objs, \(ct_obj_list, ct_name) {
    message(sprintf("  Processing cell_group = '%s'", ct_name))
    message(sprintf("  ct_obj_list class: %s", class(ct_obj_list)))
    message(sprintf("  ct_obj_list names: %s", paste(names(ct_obj_list), collapse = ", ")))
    message(sprintf("  first element class: %s", class(ct_obj_list[[1]])))
  
    # Convert inferCNV objects to gene-level data frame
    ct_prepared <- lapply(ct_obj_list, function(obj) {
        list(obj@expr.data, obj@gene_order)
    })
  
  
  gene_level_df <- load_and_prepare_infercnv_reference(ct_prepared, k = k_discrete)
  
  run_fast_cnv_pipeline(
    gene_level_df                         = gene_level_df,
    max_gap                               = max_gap,
    min_overlap_consistent_calls          = min_overlap_consistent_calls,
    min_overlap_multiple_nodes            = min_overlap_multiple_nodes,
    filter_seq_mb_init                    = filter_seq_mb_init,
    filter_seq_mb_equiv                   = filter_seq_mb_equiv,
    min_references                        = min_references,
    overlap_method_equiv_cnv_call_merge   = overlap_method,
    overlap_method_equiv_cnv_after_filter = overlap_method,
    parallel                              = parallel,
    cores                                 = cores,
    clique_mode_consistent                = clique_mode_consistent,
    removed_log_return                    = removed_log_return,
    metadata                              = metadata,
    mode                                  = "within",
    remove_ref = remove_ref,

    coding_gr                 = coding_gr,
    gene_order                = ct_prepared[[1]][[2]],
    coding_expressed_set      = coding_expressed_set,
    pct_max                   = pct_max,
    pct_floor                 = pct_floor,
    min_expr_density          = min_expr_density,
    min_coding_density        = min_coding_density,
    max_gap_mb                = max_gap_mb)
  })
      
    } else {
      
      full_results <- process_tool_cnv_runs(
        base_dir                              = tool_outdir,
        mode                                  = "within",
        tool                                  = tool,
        pattern                               = pattern,
        max_gap                               = max_gap,
        min_overlap_consistent_calls          = min_overlap_consistent_calls,
        min_overlap_multiple_nodes            = min_overlap_multiple_nodes,
        filter_seq_mb_init                    = sensitivity_floor_mb - 7.5,
        filter_seq_mb_equiv                   = 0,
        min_references                        = min_references,
        overlap_method_equiv_cnv_call_merge   = overlap_method,
        overlap_method_equiv_cnv_after_filter = overlap_method,
        parallel                              = parallel,
        cores                                 = cores,
        clique_mode_consistent                = clique_mode_consistent,
        removed_log_return                    = removed_log_return,
        metadata                              = metadata,
         k = k_discrete,
        coding_gr                 = coding_gr,
        coding_expressed_set      = coding_expressed_set,
        pct_max                   = pct_max,
        pct_floor                 = pct_floor,
        min_expr_density     = min_expr_density,
        min_coding_density   = min_coding_density,
        max_gap_mb           = max_gap_mb
      )
      
      if (is.null(full_results)) {
        stop(
          "No results from process_tool_cnv_runs.\n",
          "Check that inferCNV output exists in: ", workdir
        )
      }
    }
    
    # ---- Extract supported events — same regardless of path ---------------
    supported_events <- purrr::map(full_results, \(ct_results) {
      ct_results[["cnvs_supported_overlaped"]]
    }) |>
      purrr::compact() |>
      dplyr::bind_rows()
    
    if (nrow(supported_events) == 0L) {
      stop(
        "No supported events after Block2.\n",
        "Check pipeline parameters — consider lowering:\n",
        "  --min_references (current: ",               min_references, ")\n",
        "  --filter_seq_mb_equiv (current: ",          filter_seq_mb_equiv, ")\n",
        "  --min_overlap_consistent_calls (current: ", min_overlap_consistent_calls, ")"
      )
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
      "  %d supported events across %d cells (%.1f seconds)",
      summaries$block2$n_events,
      summaries$block2$n_cells,
      summaries$block2$runtime_s
    ))
    
    if (save_intermediate) {
      path <- file.path(workdir, "block2_supported_events.rds")
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
    
    # ---- Resolve supported_events -----------------------------------------
    supported_events <- if (start_from == "block3") {
      
      if (is.null(supported_events)) {
        stop(
          "supported_events required when start_from = 'block3'.\n",
          "Provide a data frame or path to a saved RDS file."
        )
      }
      if (is.character(supported_events)) {
        if (!file.exists(supported_events)) {
          stop("RDS path not found: ", supported_events)
        }
        message("  Loading supported_events from: ", supported_events)
        readRDS(supported_events)
      } else {
        supported_events
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
    
    if (!sample_col %in% group_cols) {
    message(
    "  sample_col '", sample_col,
    "' not in group_cols — adding it for cell_sizes"
  )
  group_cols_for_sizes <- unique(c(group_cols, sample_col))
} else {
  group_cols_for_sizes <- group_cols
}

    cell_sizes <- compute_cell_sizes(
      metadata   = metadata,
      group_cols = group_cols_for_sizes,
      cell_col   = cell_col
    )
    
    results$block3 <- list(
      cnv_annotated = cnv_annotated,
      cell_sizes    = cell_sizes
    )
    
    summaries$block3 <- list(
      n_annotated     = nrow(cnv_annotated),
      n_groups        = nrow(cell_sizes),
      cell_size_range = range(cell_sizes$n_total_cells),
      runtime_s       = (proc.time() - t3)[["elapsed"]]
    )
    
    message(sprintf(
      "  %d annotated events | %d groups | cell range: %d-%d",
      summaries$block3$n_annotated,
      summaries$block3$n_groups,
      summaries$block3$cell_size_range[1],
      summaries$block3$cell_size_range[2]
    ))
    
    if (save_intermediate) {
      path <- file.path(workdir, "block3_annotated.rds")
      saveRDS(results$block3, path)
      message("  Saved: ", path)
    }
  }
  
  # =========================================================================
  # BLOCK 4 — Cluster and score CNV loci
  # =========================================================================
  if ("block4" %in% blocks_to_run) {
    
    # ---- Resolve inputs ---------------------------------------------------
    cnv_annotated <- if (start_from == "block4") {
      if (is.null(cnv_annotated)) {
        stop("cnv_annotated required when start_from = 'block4'.")
      }
      cnv_annotated
    } else {
      results$block3$cnv_annotated
    }
    
    cell_sizes <- if (start_from == "block4") {
      if (is.null(cell_sizes)) {
        stop("cell_sizes required when start_from = 'block4'.")
      }
      cell_sizes
    } else {
      results$block3$cell_sizes
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
      summary_df                  = clustered_events$cnv_locus_summary,
      clustered_events            = clustered_events$clustered_events,
      cell_sizes                  = cell_sizes,
      by_union                    = sample_col,
      chromosome_arms             = chromosome_arms,
      min_required_cells           = min_required_cells,
      round_fun                   = ceiling,
      k = k_threshold_growth,
      sensitivity_floor_mb = sensitivity_floor_mb,
      p_arm_permission     = p_arm_permission,
      q_arm_permission     = q_arm_permission,  
      whole_chr_permission = whole_chr_permission  
    )
    
    scored_events <-  deduplicate_cnv_cells(scored_events,
                               sample_col   = sample_col,
                               cell_col     = cell_col)
                               
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
      path <- file.path(workdir, "block4_clustered_scored.rds")
      saveRDS(results$block4, path)
      message("  Saved: ", path)
    }
  }
  
  # ---- Pipeline summary ---------------------------------------------------
  total_runtime <- (proc.time() - t_start)[["elapsed"]]
  
  message(sprintf(paste0(
    "\n=== Pipeline complete ===\n",
    "  Blocks run:      %s\n",
    "  Execution mode:  %s\n",
    "  Total runtime:   %.1f seconds (%.1f minutes)"
  ),
  paste(blocks_to_run, collapse = " → "),
  execution_mode,
  total_runtime,
  total_runtime / 60
  ))
  
  f<- list(
    all_results = results,
    list(summary = list(
      blocks_run     = blocks_to_run,
      execution_mode = execution_mode,
      per_block      = summaries,
      total_runtime  = total_runtime
    ))
  )
}