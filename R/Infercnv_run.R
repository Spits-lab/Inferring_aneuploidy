# =============================================================================
# infercnv_run.R
# Function to run inferCNV on all objects produced by make_infercnv_objects()
# =============================================================================


#' Run inferCNV on all objects in a nested list
#'
#' Expects the nested list structure produced by make_infercnv_objects():
#'   list$within_cell_type[[cell_type]][[ref_group]]
#'   list$across_cell_type[[query_type]][[ref_type]]
#'
#' Output directories are created automatically:
#'   {base_outdir}/within/{cell_type}/ref_{ref_group}/
#'   {base_outdir}/across/{query_type}/ref_{ref_type}/
#'
#' @param infercnv_obj_list nested list from make_infercnv_objects()
#' @param base_outdir       root output directory (will be created if needed)
#' @param cutoff            minimum average counts per gene for reference cells
#'                          (default 0.1 — correct for RNA-seq, use 1 for 10x)
#' @param cluster_by_groups logical, cluster cells within groups (default TRUE)
#' @param HMM               logical, run HMM CNV prediction (default FALSE)
#' @param denoise           logical, apply denoising (default TRUE)
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
                                 HMM                  = FALSE,
                                 denoise              = TRUE,
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
  valid_modes <- c("within_cell_type", "across_cell_type")
  found_modes <- intersect(names(infercnv_obj_list), valid_modes)
  
  if (length(found_modes) == 0) {
    stop(
      "infercnv_obj_list does not contain expected modes. ",
      "Expected names: 'within_cell_type' and/or 'across_cell_type'. ",
      "Got: ", paste(names(infercnv_obj_list), collapse = ", ")
    )
  }
  
  dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Run log ──────────────────────────────────────────────────────────────
  run_log <- data.frame(
    mode      = character(),
    cell_type = character(),
    comp      = character(),
    out_dir   = character(),
    status    = character(),
    message   = character(),
    stringsAsFactors = FALSE
  )
  
  # ── Map mode names to output folder names ────────────────────────────────
  mode_folder_map <- c(
    within_cell_type = "within",
    across_cell_type = "across"
  )
  
  # ── Loop ─────────────────────────────────────────────────────────────────
  for (mode in found_modes) {
    
    mode_folder  <- mode_folder_map[[mode]]
    mode_objects <- infercnv_obj_list[[mode]]
    
    for (cell_type in names(mode_objects)) {
      
      type_objects <- mode_objects[[cell_type]]
      
      for (comp in names(type_objects)) {
        
        infer_obj <- type_objects[[comp]]
        
        # NULL guard — object may have failed during creation
        if (is.null(infer_obj)) {
          message("Skipping NULL object: ", mode, " / ", cell_type, " / ", comp)
          run_log <- rbind(run_log, data.frame(
            mode = mode, cell_type = cell_type, comp = comp,
            out_dir = NA, status = "skipped_null", message = "Object is NULL",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        # Build output directory
        outdir <- file.path(base_outdir, mode_folder, cell_type, comp)
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        
        # Resume check — skip if already completed
        final_obj_path <- file.path(outdir, "run.final.infercnv_obj")
        if (resume_if_exists && file.exists(final_obj_path)) {
          message("Skipping (already complete): ",
                  mode_folder, "/", cell_type, "/", comp)
          run_log <- rbind(run_log, data.frame(
            mode = mode, cell_type = cell_type, comp = comp,
            out_dir = outdir, status = "skipped_exists",
            message = "run.final.infercnv_obj already present",
            stringsAsFactors = FALSE
          ))
          next
        }
        
        message("\n── Running: ", mode_folder, " / ", cell_type, " / ", comp)
        message("   Output: ", outdir)
        
        status  <- "success"
        err_msg <- ""
        
        tryCatch({
          options(scipen = 100)
          infercnv::run(
            infercnv_obj          = infer_obj,
            out_dir               = outdir,
            cutoff                = cutoff,
            cluster_by_groups     = cluster_by_groups,
            HMM                   = HMM,
            denoise               = denoise,
            analysis_mode         = analysis_mode,
            output_format         = NA,
            no_plot               = no_plot,
            no_prelim_plot        = no_prelim_plot,
            window_length         = window_length,
            plot_probabilities    = plot_probabilities,
            plot_steps            = plot_steps,
            diagnostics           = diagnostics,
            inspect_subclusters   = inspect_subclusters
          )
          
          message("   Done: ", mode_folder, " / ", cell_type, " / ", comp)
          
        }, error = function(e) {
          status  <<- "failed"
          err_msg <<- conditionMessage(e)
          warning("FAILED: ", mode_folder, "/", cell_type, "/", comp,
                  "\n  Error: ", err_msg)
        })
        
        run_log <- rbind(run_log, data.frame(
          mode      = mode,
          cell_type = cell_type,
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
    } # cell_type loop
  } # mode loop
  
  # ── Final summary ────────────────────────────────────────────────────────
  message("\n── Run Summary ───────────────────────────────────────────────────")
  print(table(run_log$status))
  
  failed <- run_log[run_log$status == "failed", ]
  if (nrow(failed) > 0) {
    message("\nFailed runs:")
    print(failed[, c("mode", "cell_type", "comp", "message")])
  }
  
  # Save log to base_outdir
  log_path <- file.path(base_outdir, "run_log.tsv")
  write.table(run_log, log_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("\nRun log saved to: ", log_path)
  
  invisible(run_log)
}



run_infercnv_pipeline <- function(
    
  # ---- make_infercnv_objects parameters ----------------------------------
  counts_mx,
  metadata,
  cell_type_col   = "cell_type",
  gene_order_file,
  mode            = c("within", "across"),
  chr_exclude     = c("MT", "Y"),
  min_max_counts  = c(100, 1e6),
  n_splits_within = 3,
  
  # ---- run_infercnv_objects parameters -----------------------------------
  base_outdir,
  cutoff          = 0.1,
  cluster_by_groups = TRUE,
  HMM             = FALSE,
  denoise         = TRUE,
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
  if (!cell_type_col %in% colnames(metadata)) {
    stop("cell_type_col '", cell_type_col, "' not found in metadata.")
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
    "  Cell types:     %s\n",
    "  Output dir:     %s"
  ),
  mode,
  ncol(counts_mx),
  nrow(counts_mx),
  paste(unique(metadata[[cell_type_col]]), collapse = ", "),
  base_outdir
  ))
  
  t_start <- proc.time()
  
  # ---- Step 1: make inferCNV objects --------------------------------------
  message("\n[1/2] Creating inferCNV objects...")
  t_make_start <- proc.time()
  
  obj_list <- make_infercnv_objects(
    counts_mx       = counts_mx,
    metadata        = metadata,
    cell_type_col   = cell_type_col,
    gene_order_file = gene_order_file,
    mode            = mode,
    chr_exclude     = chr_exclude,
    min_max_counts  = min_max_counts,
    n_splits_within = n_splits_within
  )
  
  t_make_end <- proc.time()
  
  # ---- Lightweight sanity check between steps -----------------------------
  if (is.null(obj_list)) {
    stop("make_infercnv_objects() returned NULL — check inputs.")
  }
  if (!"objects" %in% names(obj_list[["within_cell_type"]]) &&
      mode %in% c("within")) {
    stop(
      "Expected 'objects' element in obj_list$within_cell_type. ",
      "make_infercnv_objects() may have failed silently."
    )
  }
  
  message(sprintf(
    "Objects created in %.1f seconds.",
    (t_make_end - t_make_start)[["elapsed"]]
  ))
  
  split_metadata <- NULL  # default — only populated for within/both
  
  if (mode %in% c("within")) {
    
    if (!"split_metadata" %in% names(obj_list[["within_cell_type"]])) {
      stop(
        "Expected 'split_metadata' in obj_list$within_cell_type. ",
        "Check make_infercnv_objects() output structure."
      )
    }
    
    split_metadata <- obj_list[["within_cell_type"]][["split_metadata"]]
    obj_list[["within_cell_type"]] <- obj_list[["within_cell_type"]][["objects"]]
    
    message(sprintf(
      "Within mode: split_metadata extracted (%d cell type splits)",
      length(split_metadata)
    ))
  }
  
  
  # ---- Step 2: run inferCNV -----------------------------------------------
  t_run_start <- proc.time()
  
  run_log <- run_infercnv_objects(
    infercnv_obj_list = obj_list,
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
  
  t_run_end <- proc.time()
  
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
    obj_list       = obj_list,
    run_log        = run_log,
    metadata       = split_metadata,
    runtime        = runtime
  )
}


