#' @title Score_system
#'
#' @description
#' 
#' Functions that perform further processing from the CNV performing a final filtering and Confidence Score in a Single Tool Approach
#' 
#' @author Pedro Granjo
#' @date 13-03-2026
#'



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

# Package groups
cran_packages <- c(
  "dplyr", "igraph", "BiocManager", "purrr", "tidyr"
)

bioc_packages <- c(
  "GenomicRanges", "IRanges"
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


#' Cluster CNV events within groups
#'
#' defined by one or more metadata columns.
#'
#' @param df Data frame containing CNV events.
#' @param by Character vector specifying grouping columns (e.g. dataset,
#' sample, or patient).
#' @param cluster_mode Clustering strategy passed to `cluster_cnv_events()`.
#' @param min_ovelap Minimum reciprocal overlap threshold.
#'
#' @details
#' The function splits the input data frame according to the grouping
#' variables and clusters CNV events independently within each subset.
#'
#' @return
#' Data frame with clustered CNV events including the `cnv_equiv_id` column.
cluster_cnv_events_by <- function(df, by = NULL, overlap_method = "reciprocal", min_ovelap = 0.75, parallel = F, n_cores = 1L){
  
  cols_to_remove <- intersect(c("merge_group_id", "cnv_equiv_id"), colnames(df))
  
  if (length(cols_to_remove) > 0L) {
    df <- df |> dplyr::select(-dplyr::all_of(cols_to_remove))
  }
  
  
  required_cols <- c("chr", "start", "end", "cnv_state")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  mandatory_group <- c("chr", "cnv_state")
  
  if (is.null(by)) {
    by_columns <- mandatory_group 
    return(out)
  } else{
    by_columns <- c(by,mandatory_group)
  }
  
  
  missing_by <- setdiff(by, colnames(df))
  if (length(missing_by) > 0) {
    stop("Grouping columns not found: ", paste(missing_by, collapse = ", "))
  }
  
  res <- assign_cnv_equivalence(
    df,
    min_overlap = min_ovelap,
    overlap_method         = overlap_method,
    filter_seq_mb          = 0,
    parallel               = parallel,
    by_columns = by_columns,
    n_cores = n_cores
  )
  
  return(res)
}
  
  


#' Summarise clustered CNV loci
#'
#' Computes summary statistics for clustered CNV loci including genomic
#' span, number of events, number of cells, and number of samples.
#'
#' @param df Data frame containing clustered CNV events.
#' @param by Optional grouping columns.
#' @param sample_col Column identifying samples or datasets.
#' @param cell_col Column identifying cells.
#' @param mode_col Optional column indicating event origin classification
#' (e.g. "within", "across", "both").
#'
#' @details
#' The function aggregates CNV clusters and computes:
#'
#' * genomic locus boundaries
#' * number of contributing CNV events
#' * number of unique cells
#' * number of samples
#'
#' If `mode_col` is present additional counts per mode are reported.
#'
#' @return
#' Data frame summarising CNV loci.
summarise_cnv_loci <- function(df, by = NULL,
                               sample_col = "dataset",
                               cell_col = "ds_cell",
                               mode_col = "mode") {
  
  required_cols <- c("cnv_equiv_id", "chr", "cnv_state", "start", "end", by)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!any(by %in% colnames(df))) {
    stop("Missing sample column: ", sample_col)
  }
  
  if (!cell_col %in% colnames(df)) {
    stop("Missing cell column: ", cell_col)
  }
  
  grouping_vars <- c(by, "cnv_equiv_id", "chr", "cnv_state")
  
  df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars)))
  
  # Base summary always computed
  out <- df %>%
    dplyr::summarise(
      locus_start = min(start),
      locus_end = max(end),
      locus_width = locus_end - locus_start + 1,
      locus_width_mb = locus_width / 1e6,
      n_events = dplyr::n(),
      n_cells = dplyr::n_distinct(.data[[cell_col]]),
      n_samples = dplyr::n_distinct(.data[[sample_col]]),
      .groups = "drop"
    )
  
  out %>%
    dplyr::arrange(dplyr::desc(n_cells), dplyr::desc(n_events))
}

#' Run CNV locus analysis workflow
#'
#' Performs CNV event clustering followed by locus-level summarisation.
#'
#' @param df Data frame containing CNV events.
#' @param by Optional grouping variables.
#' @param min_ovelap Minimum reciprocal overlap threshold.
#' @param cluster_mode Clustering strategy ("connected" or "complete").
#' @param sample_col Column identifying samples.
#'
#' @return
#' A list containing:
#'
#' * `clustered_events` — CNV events with cluster assignments
#' * `cnv_locus_summary` — summarised CNV loci
#'
run_cnv_locus_analysis <- function(df, by = NULL, min_ovelap = 0.75, sample_col, cell_col, overlap_method = "reciprocal",
                                   parallel = F, n_cores = 1L, removed_log_retur = F){
  
  clustered <- cluster_cnv_events_by(
    df = df,
    by = by,
    overlap_method = overlap_method,
    min_ovelap = min_ovelap,
    parallel = parallel, n_cores = n_cores
  )
  
  

  clustered_table_with_equiv_id <- clustered$results_id
  
  
  summary_tbl <- summarise_cnv_loci(
    df = clustered_table_with_equiv_id,
    by = by,
    sample_col = sample_col,
    cell_col = cell_col
  )
  
  
  if(removed_log_retur){
    list(
      clustered_events = clustered_table_with_equiv_id,
      remove_log = clustered$removed_log,
      cnv_locus_summary = summary_tbl
    )
  } else{
    list(
      clustered_events = clustered_table_with_equiv_id,
      cnv_locus_summary = summary_tbl
    )
  }

}




#' Filter CNV loci based on minimum cell thresholds
#'
#' Removes CNV loci with insufficient supporting cells and filters
#' likely artefactual events near centromeres.
#'
#' @param df CNV locus summary data frame.
#' @param min_cells_keep Global minimum cell threshold.
#' @param threshold_df Optional table containing group-specific thresholds.
#' @param group_col Column used for threshold grouping.
#'
#' @return
#' Filtered CNV locus data frame.
filter_cnv_loci <- function(
    total_chromosome_permission  = 70,
    clustered_events = NULL
) {
  
  # ---- Input validation ---------------------------------------------------
  
  required_cols <- c(
    "n_cells",
    "arm_class",
    "whole_chromosome_gain",
    "whole_chromosome_loss"
  )
  
  missing_cols <- setdiff(required_cols, colnames(clustered_events))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # ---- Threshold resolution -----------------------------------------------
  
  n_before <- nrow(clustered_events)
    
    out <- clustered_events |>
      dplyr::mutate(
        max_chr_event = pmax(
          whole_chromosome_gain,
          whole_chromosome_loss,
          na.rm = TRUE
        ),
        pass_cells      = case_when(is.na(effective_threshold) ~ F,
                                    n_cells >= effective_threshold ~ T),
        pass_centromere = !(arm_class == "p_centromere_q" &
                              max_chr_event < total_chromosome_permission)
      ) |>
      dplyr::filter(pass_cells, pass_centromere) |>
      dplyr::select(-all_of(c("max_chr_event",
                       "pass_cells", "pass_centromere", "effective_threshold")))
    
  
  # ---- Reporting ----------------------------------------------------------
  n_after         <- nrow(out)
  n_removed_cells <- length(unique(clustered_events["cell_name"])) - length(unique(out["cell_name"]))
  
  message(sprintf(paste0(
    "Cell filter summary:\n",
    "  Input:                    %d rows\n",
    "  Retained:                 %d rows\n",
    "  Removed (total):          %d rows (%.1f%%)\n",
    "  Centromere threshold:     %.0f%%"
  ),
  n_before,
  n_after,
  n_before - n_after,
  100 * (n_before - n_after) / n_before,
  total_chromosome_permission
  ))
  
  return(out)
}



#' Assign confidence levels to CNV clusters
#'
#' Classifies CNV loci as high or low confidence based on
#' minimum cell counts and locus length.
#'
#' @param summary_df Data frame containing CNV locus summaries.
#' @param high_min_cells Minimum number of cells required for high confidence.
#' @param high_threshold_df Optional table defining group-specific thresholds.
#' @param group_col Column used for joining group-specific thresholds.
#' @param min_length_mb Minimum CNV length (in Mb) required for high confidence.
#' @param length_col Column representing CNV length.
#' @param keep_low Logical indicating whether low-confidence clusters
#' should be retained.
#'
#' @return
#' Data frame with additional columns:
#'
#' * `confidence`
#' * `pass_cells`
#' * `pass_length`
#' * `low_reason`
classify_cnv_loci <- function(
    df
){
  df <- df %>% dplyr::mutate( tier = as.numeric(gsub("^tier_","",tier )),
    confidence = case_when(tier == 1 ~ "High", tier == 2 ~ "Low"))
  return(df)
}


make_tier_definitions <- function(
    base_fraction = 0.05,
    step          = 0.03,
    boundaries_mb   = c(25,50),
    fractions     = NULL,
    n_cells       = NULL,
    mode          = c("fractions", "number")
) {
  
  mode <- match.arg(mode)
  boundaries_mb <- sort(boundaries_mb,decreasing = T)
  # ---- Validate boundaries_mb -----------------------------------------------
  if (!all(diff(boundaries_mb) < 0)) {
    stop(
      "boundaries_mb must be strictly decreasing — ",
      "tier 1 should have the largest minimum size.\n",
      "  Got: ", paste(boundaries_mb, collapse = ", ")
    )
  }
  if (any(boundaries_mb < 0)) {
    stop("boundaries_mb values must be non-negative.")
  }
  
  
  # ---- Fraction resolution ------------------------------------------------
  if(mode == "fractions"){
    
    if (!is.null(n_cells)) {
      warning("n_cells ignored when selected mode is 'fractions'.")
    }
    
    if (!is.null(fractions)) {
      
      # Manual fractions supplied — validate
      if (any(fractions <= 0) || any(fractions >= 1)) {
        stop(
          "All fractions must be in (0, 1).\n",
          "  Got: ", paste(round(fractions, 4), collapse = ", ")
        )
      }
      
      if (!all(diff(fractions) > 0)) {
        stop(
          "fractions must be strictly increasing — ",
          "tier 1 (largest CNVs) should have the smallest fraction.\n",
          "  Got: ", paste(round(fractions, 4), collapse = ", ")
        )
      }
      n_tiers <- length(fractions)
      
    } else {
      n_tiers <- length(boundaries_mb)
      # Auto-generate fractions from base + step
      fractions <- base_fraction + (0:(n_tiers - 1)) * step
      
      if (any(fractions <= 0) || any(fractions >= 1)) {
        stop(sprintf(
          "Generated fractions outside (0, 1) — adjust base_fraction or step.\n  Got: %s",
          paste(round(fractions, 4), collapse = ", ")
        ))
      }
      
      
    } 
    out <- data.frame(
      tier        = paste0("tier_", seq_len(n_tiers)),
      boundaries_mb = boundaries_mb,
      fraction    = fractions,
      stringsAsFactors = FALSE
    ) 
  } else if(mode == "number"){
    
    if (!is.null(fractions)) {
      warning("fractions ignored when mode = 'number'.")
    }
    
    if (is.null(n_cells)) {
      stop("n_cells must be provided when mode = 'number'.")
    }
    
    n_tier <- length(n_cells)
    
    if (any(n_cells < 1L)) {
      stop("n_cells values must be >= 1.")
    }
    
    if (n_tiers > 1L && !all(diff(n_cells) > 0)) {
      stop(
        "n_cells must be strictly increasing — ",
        "tier 1 (largest CNVs) should require the fewest cells.\n",
        "  Got: ", paste(n_cells, collapse = ", ")
      )
    }
    
    out <- data.frame(
      tier        = paste0("tier_", seq_len(n_tiers)),
      boundaries_mb = boundaries_mb,
      n_cells     = as.integer(n_cells),
      stringsAsFactors = FALSE
    )
    
  }
  
  # ---- Report -------------------------------------------------------------
  message(sprintf(
    "Tier definitions created (%d tiers, mode = '%s'):", n_tiers, mode
  ))
  
  if (mode == "fractions") {
    message(paste(
      sprintf("  tier_%d:\n min_size = %.0f Mb\n, fraction = %.3f",
              seq_len(n_tiers), boundaries_mb, fractions),
      collapse = "\n"
    ))
  } else {
    message(paste(
      sprintf("  tier_%d:\n  min_size = %.0f Mb\n, n_cells = %d",
              seq_len(n_tiers), boundaries_mb, n_cells),
      collapse = "\n"
    ))
  }
  
  out
}


resolve_tier_thresholds <- function(method = c("auto", "single", "manual"),
                                    base_fraction = 0.05,
                                    step          = 0.03,
                                    boundaries_mb   = c(25,50),
                                    fractions     = NULL,
                                    n_cells       = NULL,
                                    mode          = c("fractions", "number")
                                    
) {
  
  method <- match.arg(method)
  mode   <- match.arg(mode)
  
  # ---- Method dispatch ----------------------------------------------------
  if (method == "automated") {
    warning(
      "Automated method uses fraction-based approach only. ",
      "For number-based thresholds use method = 'vector' and supply n_cells."
    )
    mode     <- "fractions"
    fractions <- NULL
    n_cells   <- NULL
  }
  
  if (method == "fixed") {
    message("Fixed threshold — single tier applied uniformly.")
    
    if (!is.null(fractions) && !is.null(n_cells)) {
      stop("Provide either fractions or n_cells, not both.")
    }
    if (is.null(fractions) && is.null(n_cells)) {
      stop("Provide either a fraction or n_cells value for the fixed threshold.")
    }
    
    mode    <- if (!is.null(fractions)) "fractions" else "number"
  }
  
  out <- make_tier_definitions(
    base_fraction = base_fraction,
    step          = step,
    boundaries_mb   = boundaries_mb,
    fractions     = fractions,
    n_cells       = n_cells,
    mode          = mode
  )
  
  out
}


prepare_cnv_thresholds <- function(
    summary_df,
    clustered_events,
    by_union = NULL,
    method = c("auto", "single", "manual"),
    max_tiers = 2,
    base_fraction = 0.05,
    step          = 0.03,
    cell_sizes = NULL,
    boundaries_mb   = c(100, 50, 20),
    fractions     = NULL,
    n_cells       = NULL,
    mode          = c("fractions", "number"),
    min_cap_threshold = 2,
    max_cap_threshold = 25,
    round_fun = ceiling
){
  if(is.null(cell_sizes)){
    stop("Your dataframe with total cell numbers is not coming through")
  }
  
  mode <- match.arg(mode)
  
  required_prepared_cols <- c(
    "cnv_equiv_id", "cnv_length_mb", "n_cells", "n_total_cells"
  )
  
  already_merged <- all(required_prepared_cols %in% colnames(clustered_events))
  
  if (already_merged) {
    message("Input already contains required columns — skipping merge step.")
    merged_df <- clustered_events
    
  } else {
    # Validate summary_df is provided
    if (is.null(summary_df)) {
      stop(
        "summary_df must be provided when clustered_events is missing columns: ",
        paste(setdiff(required_prepared_cols, colnames(clustered_events)), collapse = ", ")
      )
    }
    
    
    join_cols <- c(by_union, "cnv_equiv_id")
    merged_df <- clustered_events %>%
      dplyr::left_join(
        summary_df %>% dplyr::select(dplyr::all_of(join_cols), n_cells),
        by = join_cols
      ) %>%
      dplyr::left_join(cell_sizes, by = by_union
      )
    
    if (anyNA(merged_df$n_total_cells)) {
      missing <- length(unique(merged_df[["cnv_equiv_id"]][is.na(merged_df$n_total_cells)]))
      stop(sprintf("No cell count found in %d groups: ", missing))
    }
  }
  
  if(length(boundaries_mb) > max_tiers){
    stop(sprintf("You are trying to make more than %d groups.",
                 max_tiers))
  }
  
  tier_table <- resolve_tier_thresholds(
    method = method,
    base_fraction = base_fraction,
    step          = step,
    boundaries_mb   = boundaries_mb,
    fractions     = fractions,
    n_cells       = n_cells,
    mode          = mode)
  
  threshold_col  <- if ("fraction" %in% colnames(tier_table)) "fraction" else "n_cells"
  tier_sorted    <- tier_table[order(tier_table$boundaries_mb, decreasing = FALSE), ]
  boundaries_asc <- tier_sorted$boundaries_mb
  thresholds_asc <- tier_sorted[[threshold_col]]
  
  thresholded_df <- merged_df |>
    dplyr::mutate(
      interval_idx  = findInterval(cnv_length_mb, boundaries_asc),
      safe_idx      = dplyr::if_else(interval_idx == 0L, NA_integer_, interval_idx),
      raw_threshold = thresholds_asc[safe_idx],
      tier          = dplyr::if_else(
        is.na(safe_idx),
        "below_threshold",
        tier_sorted$tier[safe_idx]
      )
    ) |>
    dplyr::select(-interval_idx, -safe_idx)
  
  # Apply effective_threshold computation outside mutate
  # threshold_col is a scalar — evaluate once, apply vectorised
  if (threshold_col == "fraction") {
    thresholded_df <- thresholded_df |>
      dplyr::mutate(
        effective_threshold = round_fun(pmin(max_cap_threshold, pmax(raw_threshold * n_total_cells, min_cap_threshold)))
      )
  } else if (threshold_col == "number"){
    thresholded_df <- thresholded_df |>
      dplyr::mutate(
        effective_threshold = round_fun(pmin(max_cap_threshold, pmax(raw_threshold, min_cap_threshold)))
      )
  }
  thresholded_df <- thresholded_df |>
    dplyr::select(-raw_threshold)
  return(thresholded_df)
}



score_cnv_clusters <- function(
    summary_df,
    clustered_events,
    cell_sizes,
    by_union,
    boundaries_mb        = c(50, 19),
    base_fraction        = 0.05,
    step                 = 0.03,
    fractions            = NULL,
    n_cells              = NULL,
    threshold_method     = c("auto", "single", "manual"),
    threshold_mode       = c("fractions", "number"),
    min_cap_threshold    = 2L,
    max_cap_threshold    = 25L,
    max_tiers            = 2L,
    total_chromosome_permission = 70,
    round_fun = ceiling
) {
  
  threshold_method <- match.arg(threshold_method)
  threshold_mode <- match.arg(threshold_mode)
  
  # ---- Step 1: prepare thresholds — merge only if needed ------------------
  thresholded_df <- prepare_cnv_thresholds(
    summary_df        = summary_df,
    clustered_events  = clustered_events,
    by_union          = by_union,
    cell_sizes        = cell_sizes,
    method            = threshold_method,
    max_tiers         = max_tiers,
    base_fraction     = base_fraction,
    step              = step,
    boundaries_mb     = boundaries_mb,
    fractions         = fractions,
    n_cells           = n_cells,
    mode              =  threshold_mode,
    min_cap_threshold = min_cap_threshold,
    max_cap_threshold = max_cap_threshold,
    round_fun         = ceiling
  )
  
  # ---- Step 2: filter -----------------------------------------------------
  filtered_df <- filter_cnv_loci(
    clustered_events                          = thresholded_df,
    total_chromosome_permission = total_chromosome_permission
  )
  
  # ---- Step 3: classify confidence ----------------------------------------
  # thresholded_df already has effective_threshold — pass directly
  # No second prepare needed — filter only removes rows, columns intact
  classify_cnv_loci(
    df         = filtered_df)
}