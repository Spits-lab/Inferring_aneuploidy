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
  
  
  
  

