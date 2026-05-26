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
    pattern
) {
  tool <- match.arg(tool)
  
  switch(tool,
         "infercnv" = load_infercnv_data(base_dir, ref_dirs, pattern)#,
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
    cell_type_col     = "cell_type",
    gene_order_file   = NULL,
    chr_exclude       = c("MT", "Y"),
    min_max_counts    = c(100, 1e6),
    n_splits_within   = 3,
    base_outdir,
    cutoff            = 0.1,
    cluster_by_groups = TRUE,
    HMM               = FALSE,
    denoise           = TRUE,
    analysis_mode     = "subclusters",
    window_length     = 140,
    no_plot           = TRUE,
    resume_if_exists  = TRUE
) {
  
  tool <- match.arg(tool)
  
  switch(tool,
         
         "infercnv" = run_infercnv_pipeline(
           counts_mx         = counts_mx,
           metadata          = metadata,
           cell_type_col     = cell_type_col,
           gene_order_file   = gene_order_file,
           mode              = "within",
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
    stop(sprintf(
      paste0(
        "No cell types found in '%s' for mode '%s'.\n",
        "Check that inferCNV has been run and results exist:\n",
        "  %s/{cell_type}/{A,B,C}/run.final.infercnv_obj"
      ),
      cell_dir, mode, cell_dir
    ))
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





run_full_cnv_pipeline <- function(
    
  # ---- Entry point -------------------------------------------------------
  start_from        = c("block1", "block2", "block3", "block4"),
  save_intermediate = FALSE,
  outdir            = NULL,
  # ---- Block 1 -----------------------------------------------------------
  counts_mx         = NULL,
  metadata          = NULL,
  cell_type_col     = "cell_type",
  gene_order_file   = NULL,
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
  parallel                              = FALSE,
  cores                                 = 1L,
  clique_mode_consistent                = "connected",
  removed_log_return                    = FALSE,
  supported_events = NULL,
  
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
    
    results$block1 <- run_cnv_tool(
      tool              = tool,
      counts_mx         = counts_mx,
      metadata          = metadata,
      cell_type_col     = cell_type_col,
      gene_order_file   = gene_order_file,
      mode              = c("within"),
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
    
    
    if (!is.null(results$block1$obj_list$within_cell_type$split_metadata)) {
      metadata <- results$block1$obj_list$within_cell_type$split_metadata
      message("metadata enriched with split_group column (within mode).")
    }
    
    # Pass base_outdir as base_dir for block2
    base_dir <- base_outdir
    
    summaries$block1 <- list(runtime_s = (proc.time() - t1)[["elapsed"]])
    
    if (save_intermediate) {
      meta_path <- file.path(base_dir, "split_metadata.rds")
      saveRDS(meta_path)
      message("Enriched metadata saved")
      path <- file.path(outdir, "block1_infercnv.rds")
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
    
    # ---- Fast path — Block 1 objects already in memory --------------------
    if (!is.null(results$block1) && !is.null(results$block1$obj_list)) {
      
      within_cell_type <- results$block1$obj_list[["within_cell_type"]]
      
      if (is.null(within_cell_type)) {
        stop("within_cell_type not found in obj_list.")
      }
      
      # Extract split_metadata — enriched metadata with split_group
      within_objs <- within_cell_type[["objects"]]
      
      if (is.null(metadata)) {
        stop("split_metadata not found in within_cell_type.")
      }
      if (is.null(within_objs)) {
        stop("objects not found in within_cell_type.")
      }
      
      message(sprintf(
        "  split_metadata loaded: %d cells across %d groups",
        nrow(metadata),
        dplyr::n_distinct(metadata$split_group)
      ))
      
      full_results <- purrr::imap(within_objs, \(ct_df, ct_name) {
        
        message(sprintf(
          "  Processing mode='within', cell_type='%s'", ct_name
        ))
        
        # no load_and_prepare_infercnv_reference needed
        run_fast_cnv_pipeline(
          gene_level_df                         = ct_df,
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
          mode                                  = "within"
        )
      })
    } else {
      
      # Re-entry: resolve base_dir and metadata
      if (start_from == "block2") {
        
        if (is.null(base_dir)) stop("base_dir required when start_from = 'block2'.")
        
        # Try loading split_metadata from saved RDS if metadata not provided
        if (is.null(metadata)) {
          meta_path <- file.path(base_dir, "split_metadata.rds")
          if (file.exists(meta_path)) {
            message("  Loading split_metadata from: ", meta_path)
            metadata <- readRDS(meta_path)
          } else {
            stop(
              "metadata required for block2. ",
              "No split_metadata.rds found in base_dir either.\n",
              "Provide metadata directly or re-run Block 1 to generate it."
            )
          }
        }
      }
      
      full_results <- process_tool_cnv_runs(
        base_dir                              = base_dir,
        mode                                  = "within",
        tool                                  = tool,
        pattern                               = pattern,
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
        metadata                              = metadata
      )
      
      if (is.null(full_results)) {
        message("No results returned from process_tool_cnv_runs — block2 returning NULL.")
        results$block2 <- NULL
        summaries$block2 <- list(n_events = 0L, n_cells = 0L,
                                 runtime_s = (proc.time() - t2)[["elapsed"]])
        return(c(results, list(summary = list(
          blocks_run    = blocks_to_run,
          per_block     = summaries,
          total_runtime = (proc.time() - t_start)[["elapsed"]]
        ))))
      }
    
    }
    
    # ---- Extract supported events — same regardless of path ---------------
    supported_events <- purrr::map(full_results, \(ct_results) {
      ct_results[["cnvs_supported_overlaped"]]
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
      "  %d supported events across %d cells (%.1f seconds)",
      summaries$block2$n_events,
      summaries$block2$n_cells,
      summaries$block2$runtime_s
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
