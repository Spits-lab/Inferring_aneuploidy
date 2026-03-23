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
  "dplyr", "tidyr", "data.table", "ggplot2", "patchwork", "cowplot"
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
merge_nearby_regions <- function(df, max_gap = 100000) {
  # ---- schema validation (minimal, mirrors Python intent) ----
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
  
  # ---- ensure genomic ordering ----
  df <- df |>
    arrange(reference, cell_name, chr, state, start)
  
  # ---- core merging logic ----
  merged_df <- df |>
    group_by(reference, cell_name, chr, state) |>
    arrange(start, .by_group = TRUE) |>
    mutate(
      gap = start - lag(stop),
      new_block = is.na(gap) | gap > max_gap,
      merge_id = cumsum(new_block)
    ) |>
    group_by(reference, cell_name, chr, state, merge_id) |>
    summarise(
      cnv_state  = dplyr::first(state),
      start      = min(start),
      end        = max(stop),
      n_segments = n(),
      .groups    = "drop"
    )
  
  # ---- post-merge sanity check ----
  merged_df <- merged_df |>
    arrange(reference, cell_name, chr, cnv_state, start)
  
  overlap <- merged_df |>
    group_by(reference, cell_name, chr, cnv_state) |>
    mutate(overlap = start <= lag(end)) |>
    filter(!is.na(overlap) & overlap)
  
  if (nrow(overlap) > 0) {
    stop("Overlapping CNVs detected within the same CNV state")
  }
  
  # ---- sanity check: row reduction ----
  n_before <- nrow(df)
  n_after  <- nrow(merged_df)
  
  if (n_after < n_before) {
    message(sprintf(
      "CNV segments merged successfully: before = %d rows, after = %d rows",
      n_before, n_after
    ))
  } else {
    warning(sprintf(
      "No reduction in CNV rows after merging: before = %d, after = %d",
      n_before, n_after
    ))
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
    intersection_len / pmin(q_len, s_len)
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



#' Assign CNV equivalence IDs using overlap-graph connected components
#'
#' For each (cell_name, chr, cnv_state) group, segments that mutually satisfy
#' the reciprocal overlap threshold are assigned the same cnv_equiv_id via
#' connected-component clustering. Overlap is computed on ORIGINAL segment
#' coordinates, not an expanding window, to avoid transitive chaining bias.
#'
#' @param df A data frame with columns: cell_name, chr, cnv_state, start, end.
#' @param min_reciprocal_overlap Minimum reciprocal overlap to consider two
#'   segments equivalent. Default 0.5.
#'
#' @return The input data frame with an added integer column `cnv_equiv_id`.
#'   IDs are globally unique across all groups.
assign_cnv_equivalence <- function(df, min_reciprocal_overlap = 0.5, overlap_method  = "reciprocal") {
  
  required_cols <- c("cell_name", "chr", "cnv_state", "start", "end")
  missing_cols  <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate coordinate integrity before doing any interval arithmetic
  if (any(df$start > df$end, na.rm = TRUE)) {
    stop("Detected rows where start > end. Check upstream segmentation.")
  }
  
  # Sort once globally — findOverlaps requires sorted GRanges
  df <- df |>
    dplyr::arrange(cell_name, chr, cnv_state, start)
  
  # Global equiv ID counter — incremented per connected component across
  # all groups so IDs are unique in the output dataframe
  global_equiv_counter <- 0L
  df$cnv_equiv_id      <- NA_integer_
  
  # Split into (cell_name x chr x cnv_state) groups
  # We retain row indices into the original df to avoid the which() lookup
  group_indices <- df |>
    dplyr::mutate(.row_idx = dplyr::row_number()) |>
    dplyr::group_by(cell_name, chr, cnv_state) |>
    dplyr::group_split()
  
  for (grp in group_indices) {
    
    n <- nrow(grp)
    
    # Single-segment group: trivially its own equivalence class
    if (n == 1L) {
      global_equiv_counter <- global_equiv_counter + 1L
      df$cnv_equiv_id[grp$.row_idx] <- global_equiv_counter
      next
    }
    
    # Build a GRanges for this group using original (unextended) coordinates
    # We use seq_along as names so we can recover indices after findOverlaps
    gr <- GenomicRanges::GRanges(
      seqnames = grp$chr,
      ranges   = IRanges::IRanges(start = grp$start, end = grp$end),
      strand   = "*"
    )
    names(gr) <- as.character(seq_len(n))
    
    # Find all pairwise overlaps within this group
    # type = "any" catches partial overlaps; we filter by RO threshold below
    hits <- GenomicRanges::findOverlaps(gr, gr, type = "any", select = "all")
    
    # Remove self-hits
    hits <- hits[S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits)]
    
    if (length(hits) == 0L) {
      # No overlaps at all — every segment is its own equivalence class
      new_ids <- seq(
        from = global_equiv_counter + 1L,
        by   = 1L,
        length.out = n
      )
      df$cnv_equiv_id[grp$.row_idx] <- new_ids
      global_equiv_counter <- global_equiv_counter + n
      next
    }
    
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    
    # Replace the previous overlap_fn() call with:
    scores <- compute_overlap(
      q_start = grp$start[q_idx],
      q_end   = grp$end[q_idx],
      s_start = grp$start[s_idx],
      s_end   = grp$end[s_idx],
      method  = overlap_method
    )
    
    # Keep only pairs that meet the reciprocal overlap threshold
    passing <- scores >= min_reciprocal_overlap
    q_pass  <- q_idx[passing]
    s_pass  <- s_idx[passing]
    
    if (length(q_pass) == 0L) {
      # Overlaps exist but none pass the threshold — all separate classes
      new_ids <- seq(
        from = global_equiv_counter + 1L,
        by   = 1L,
        length.out = n
      )
      df$cnv_equiv_id[grp$.row_idx] <- new_ids
      global_equiv_counter <- global_equiv_counter + n
      next
    }
    
    # Connected components on the overlap graph
    # Each node is a segment (1..n), edges connect segments that pass RO
    # Components define equivalence classes — no sequential walk, no chaining
    g <- igraph::graph_from_edgelist(
      cbind(q_pass, s_pass),
      directed = FALSE
    )
    
    # Add isolated vertices (segments with no passing overlaps) so every
    # segment gets a component assignment
    g <- igraph::add_vertices(g, n - igraph::vcount(g))
    
    component_ids <- igraph::components(g)$membership
    
    # Offset component IDs by the global counter to ensure global uniqueness
    df$cnv_equiv_id[grp$.row_idx] <- component_ids + global_equiv_counter
    global_equiv_counter <- global_equiv_counter + max(component_ids)
  }
  
  df
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


#' Run a fast CNV event consolidation pipeline
#'
#' Runs a CNV processing workflow from gene-level calls to merged CNV events
#' with reference support summaries.
#'
#' @param gene_level_df A gene-level CNV data frame.
#' @param max_gap Maximum genomic gap allowed when merging nearby segments.
#' @param min_reciprocal_overlap Minimum reciprocal overlap for equivalence
#'   assignment.
#' @param min_references Minimum number of references required to keep a CNV.
#' @param overlap_method Select the overlap method
run_fast_cnv_pipeline <- function(
    gene_level_df,
    max_gap = 100000,
    min_reciprocal_overlap = 0.5,
    min_references = 2,
    overlap_method = "reciprocal"
) {
  
  message("→ Collapsing genes to segments")
  segments <- collapse_genes_to_cnv_segments(gene_cnv_df = gene_level_df)
 
  message("→ Merging nearby CNVs")
  merged <- merge_nearby_regions(df = segments, max_gap = max_gap)
 
  message("→ Assigning CNV equivalence")
  equiv <- assign_cnv_equivalence(
    df = merged,
    min_reciprocal_overlap = min_reciprocal_overlap,
    overlap_method = "reciprocal"
  )
  
  cnv_events <- summarize_cnv_support(equiv)
  
  message("→ Filtering CNVs by reference support")
  
  supported_events <- filter_cnv_events(
    cnv_events,
    min_references = min_references
  )
  
  
  list(
    cnvs_per_segment   = equiv,             # one row per segment, with equiv IDs - IDs which tells which CNVs overlap each other
    cnvs_summarized    = cnv_events,        # one row per equiv group, with support counts
    cnvs_supported     = supported_events  # filtered to min_references threshold
  )
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
      chr_arms <- arms_by_chr[[cnv_df$chr[i]]]
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





#' Plot the density distribution of CNV lengths
#'
#' Creates a density plot of CNV lengths in megabases, optionally restricted
#' to a single CNV state such as gain or loss.
#'
#' @param dt A data frame containing at least the columns cnv_length_mb
#'   and optionally cnv_state.
#' @param state Optional character string specifying a CNV state to subset,
#'   such as "gain" or "loss". If NULL, all rows are used.
#' @param thresholds Numeric vector of CNV length thresholds to display as
#'   vertical dashed lines.
#' @param fill_color Fill color for the density polygon. Default is
#' @param title Plot title. Default is "CNV length distribution".
#'
#' @return A ggplot2 object.
plot_cnv_density <- function(
    dt,
    state = NULL,
    thresholds = c(5, 25, 50),
    fill_color = "grey40",
    title = "CNV length distribution"
) {
  
  # Subset efficiently (no copy)
  if (!is.null(state)) {
    dt_sub <- dt[dt$cnv_state == state,]
  } else {
    dt_sub <- dt
  }
  
  ggplot(dt_sub, aes(x = cnv_length_mb)) +
    geom_density(fill = fill_color, alpha = 0.4) +
    geom_vline(
      xintercept = thresholds,
      linetype = "dashed",
      alpha = 0.6
    ) +
    labs(
      title = title,
      x = "CNV length (Mb)",
      y = "Density"
    ) +
    theme_minimal()
}



#' Plot overall and state-specific CNV length distributions
#'
#' Builds three density plots showing CNV length distributions for all CNVs,
#' gain events, and loss events, then stacks them vertically.
#'
#' @param dt A data frame containing CNV length information, including
#'   cnv_length_mb and cnv_state.
#' @param thresholds Numeric vector of thresholds to display as dashed vertical
#'   lines in each panel. Default is c(5, 25, 50).
#'
#' @return A combined patchwork plot object.
plot_all_cnv_distributions <- function(
    dt,
    thresholds = c(5, 25, 50)
) {
  
  p_overall <- plot_cnv_density(
    dt = dt,
    state = NULL,
    thresholds = thresholds,
    fill_color = "grey40",
    title = "Overall CNV length distribution"
  )
  
  p_gain <- plot_cnv_density(
    dt = dt,
    state = "gain",
    thresholds = thresholds,
    fill_color = "steelblue",
    title = "Gain CNV length distribution"
  )
  
  p_loss <- plot_cnv_density(
    dt = dt,
    state = "loss",
    thresholds = thresholds,
    fill_color = "firebrick",
    title = "Loss CNV length distribution"
  )
  
  p_overall / p_gain / p_loss
}




#' Plot density distributions by level and type
#'
#' Creates a density plot for percentage values within a selected level,
#' grouped by type.
#'
#' @param level_name Character string specifying which level in
#'   plot_long$level to plot.
#' @param plot_long A long-format data frame containing at least the columns
#'   level, type, and percentage.
#' @param threshold Numeric threshold shown as a dashed vertical line.
#'
#' @return A ggplot2 object.
#'
#' @details
#' The function also computes the mean percentage per type, although the
#' current vertical reference line uses the supplied threshold value.
make_density_plot <- function(level_name,plot_long, threshold) {
  
  df_sub <- filter(plot_long, level == level_name)
  
  # compute means per type
  mean_df <- df_sub %>%
    group_by(type) %>%
    summarise(mean_value = mean(percentage, na.rm = TRUE),
              .groups = "drop")
  
  ggplot(df_sub,
         aes(x = percentage, fill = type, color = type)) +
    geom_density(alpha = 0.3, linewidth = 1) +
    
    # vertical mean lines (like abline)
    geom_vline(data = mean_df,
               aes(xintercept = threshold),
               linetype = "dashed",
               linewidth = 1.2,
               show.legend = FALSE) +
    
    labs(
      title = level_name,
      x = "Percentage",
      y = "Density"
    ) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top")
}





################################################################
## Functions to Process Information for conjoined heatmap#######
#################################################################
#' Prepare genome-wide chromosome and arm coordinates
#'
#' Builds chromosome-level cumulative coordinates and lookup tables for
#' chromosome arms, enabling mapping of chromosome-local intervals into a
#' genome-wide coordinate system.
#'
#' @param chromosome_arms A data frame containing chromosome arm annotation.
#'   It should include at least \code{chr}, \code{arm}, \code{arm_start},
#'   \code{arm_end}, and \code{arm_length}.
#'
#' @return A named list containing information the different fraction from the chromossome
#'
prepare_genome_structure <- function(chromosome_arms) {
  
  chromosome_lengths <- chromosome_arms %>%
    group_by(chr) %>%
    summarise(chr_length = sum(arm_length), .groups = "drop") %>%
    mutate(
      chr_num = suppressWarnings(as.numeric(gsub("chr", "", chr)))
    ) %>%
    arrange(chr_num, chr) %>%
    select(-chr_num) %>%
    mutate(
      chr_start = lag(cumsum(chr_length), default = 0),
      chr_end   = chr_start + chr_length
    )
  
  # Precompute arm lookup tables (vectorized, no match needed later)
  arm_lookup <- chromosome_arms %>%
    left_join(chromosome_lengths, by = "chr") %>%
    mutate(
      genome_arm_start = chr_start + arm_start,
      genome_arm_end   = chr_start + arm_end
    )
  
  p_lookup <- arm_lookup %>%
    filter(arm == "p") %>%
    select(chr,
           p_genome_start = genome_arm_start,
           p_genome_end   = genome_arm_end)
  
  q_lookup <- arm_lookup %>%
    filter(arm == "q") %>%
    select(chr,
           q_genome_start = genome_arm_start,
           q_genome_end   = genome_arm_end)
  
  list(
    chromosome_lengths = chromosome_lengths,
    p_lookup = p_lookup,
    q_lookup = q_lookup,
    arm_lookup = arm_lookup
  )
}

#' Map CNV intervals to genome-wide coordinates
#'
#' Converts chromosome-level CNV coordinates into cumulative genome-wide
#' coordinates for visualization. The function joins chromosome and chromosome-arm
#' lookup tables, computes genome-wide start and end positions, optionally
#' expands large events to full chromosome or arm boundaries for plotting,
#' and encodes CNV state numerically.
#'
#' @param cnv_filtered A data frame of CNV events
#' @param genome_structure
#' @param threshold Optional numeric threshold used to expand plotted CNV
#'   intervals to whole-chromosome or chromosome-arm boundaries.
#' @return A data frame with additional columns
map_cnv_to_genome <- function(cnv_filtered,
                              genome_structure,
                              threshold = NULL) {
  
  chromosome_lengths <- genome_structure$chromosome_lengths
  p_lookup <- genome_structure$p_lookup
  q_lookup <- genome_structure$q_lookup
  
  
  cnv_mapped <- cnv_filtered %>%
    left_join(chromosome_lengths, by = "chr") %>%
    left_join(p_lookup, by = "chr") %>%
    left_join(q_lookup, by = "chr") %>%
    mutate(
      genome_start = chr_start + start,
      genome_end   = chr_start + end
    )
  
  # If threshold is NULL → no override logic
  if (is.null(threshold)) {
    
    cnv_mapped <- cnv_mapped %>%
      mutate(
        genome_start_plot = genome_start,
        genome_end_plot   = genome_end
      )
    
  } else {
    
    cnv_mapped <- cnv_mapped %>%
      mutate(
        genome_start_plot = case_when(
          whole_chromosome_gain > threshold |
            whole_chromosome_loss > threshold ~ chr_start,
          
          p_arm_gain > threshold |
            p_arm_loss > threshold ~ p_genome_start,
          
          q_arm_gain > threshold |
            q_arm_loss > threshold ~ q_genome_start,
          
          TRUE ~ genome_start
        ),
        genome_end_plot = case_when(
          whole_chromosome_gain > threshold |
            whole_chromosome_loss > threshold ~ chr_end,
          
          p_arm_gain > threshold |
            p_arm_loss > threshold ~ p_genome_end,
          
          q_arm_gain > threshold |
            q_arm_loss > threshold ~ q_genome_end,
          
          TRUE ~ genome_end
        )
      )
  }
  
  cell_order <- cnv_filtered %>%
    distinct(cell_name, embryo, stage) %>%
    arrange(embryo, stage)
  
  cnv_mapped <- cnv_mapped %>%
    mutate(
      cell_name = factor(cell_name, levels = cell_order$cell_name),
      cell_id = as.numeric(cell_name),
      cnv_state_numeric = case_when(
        cnv_state == "gain" ~ 1,
        cnv_state == "loss" ~ -1,
        TRUE ~ 0
      )
    ) %>%
    arrange(embryo, cell_id)
}

#' Plot a chromosome ideogram
#'
#' Creates a simple ideogram-style plot of chromosome arms aligned to genome-wide
#' coordinates.
#'
#' @param genome_structure Output from prepare_genome_structure().
#' @param max_cell_id Maximum cell index used to place the ideogram above the
#'   heatmap.
#' @param arm_colors Named vector of colors for chromosome arms.
#'
#' @return A ggplot2 object.
plot_ideogram <- function(genome_structure,
                          max_cell_id,
                          arm_colors = c("p"="#4DBBD5",
                                         "cen"="black",
                                         "q"="#E64B35")) {
  
  chromosome_lengths <- genome_structure$chromosome_lengths
  arm_lookup <- genome_structure$arm_lookup
  
  arm_plot <- arm_lookup %>%
    mutate(
      ymin = max_cell_id + 1,
      ymax = ymin + 1
    )
  
  ggplot(arm_plot) +
    geom_rect(aes(
      xmin = genome_arm_start,
      xmax = genome_arm_end,
      ymin = ymin,
      ymax = ymax,
      fill = arm
    ), color = NA) +
    scale_fill_manual(values = arm_colors, name = "Arm") +
    theme_void() +
    scale_x_continuous(
      breaks = chromosome_lengths$chr_start,
      labels = chromosome_lengths$chr,
      expand = c(0,0),
      limits = c(0, max(chromosome_lengths$chr_end) +1e6)
    )
}



#' Plot a CNV heatmap across the genome
#'
#' Draws a heatmap-like representation of CNV intervals across genome-wide
#' coordinates and cells.
#'
#' @param cnv_mapped A mapped CNV data frame produced by
#'   map_cnv_to_genome().
#' @param chromosome_lengths Chromosome coordinate table from
#'   prepare_genome_structure().
#' @param boundary_lines Optional numeric vector of horizontal boundary positions.
#' @param show_x_labels Logical; whether chromosome labels should be displayed on
#'   the x-axis.
#'
#' @return A ggplot2 object.
plot_cnv_heatmap <- function(cnv_mapped,
                             chromosome_lengths,
                             boundary_lines = NULL,
                             show_x_labels = FALSE) {
  
  max_plot_idx <- max(cnv_mapped$plot_idx)
  
  p <- ggplot(cnv_mapped) +
    geom_rect(aes(
      xmin = genome_start_plot,
      xmax = genome_end_plot,
      ymin = plot_idx - 0.5,
      ymax = plot_idx + 0.5,
      fill = cnv_state_numeric
    )) +
    geom_segment(
      data = chromosome_lengths,
      aes(x = chr_end,
          xend = chr_end,
          y = 0,
          yend = max_plot_idx),
      color = "black",
      linewidth = 0.3
    ) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0
    ) +
    scale_x_continuous(
      breaks = chromosome_lengths$chr_start,
      labels = chromosome_lengths$chr,
      expand = c(0,0),
      limits = c(0, max(chromosome_lengths$chr_end))
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = if (show_x_labels)
        element_text(size = 10, angle = 45, hjust = 1)
      else element_blank()
    ) +
    labs(x = "Chromosome", y = "Cells", fill = "CNV State")
 
  if (!is.null(boundary_lines) & length(boundary_lines)  < length(cnv_mapped$plot_idx)/2 +50) { #Change this later
    p <- p + geom_hline(
      yintercept = boundary_lines,
      color = "black",
      linewidth = 0.1,
      alpha = 0.7
    )
  }
  
  p
}



#' Combine a CNV heatmap with an ideogram
#'
#' Adds an ideogram plot above a CNV heatmap using a grob annotation.
#'
#' @param heatmap_plot A heatmap plot created by plot_cnv_heatmap().
#' @param ideogram_plot An ideogram plot created by plot_ideogram().
#' @param cnv_mapped A mapped CNV data frame used to determine y-axis placement.
#' @param remove_legends Logical; whether to remove legends from both plots.
#'
#' @return A combined ggplot2 object.
assemble_heatmap_with_ideogram <- function(heatmap_plot,
                                           ideogram_plot,
                                           cnv_mapped,
                                           remove_legends = T) {
  
  if (remove_legends) {
    heatmap_plot  <- heatmap_plot  + theme(legend.position = "none")
    ideogram_plot <- ideogram_plot + theme(legend.position = "none")
  }
  
  # Convert ideogram to grob
  ideogram_grob <- ggplotGrob(ideogram_plot)
  
  # Add ideogram as annotation to heatmap
  p_combined <- heatmap_plot +
    annotation_custom(
      grob = ideogram_grob,
      xmin = -Inf, xmax = Inf,
      ymin = max(cnv_mapped$plot_idx) + 1,  
      ymax = max(cnv_mapped$plot_idx) + 3
    )
  return(p_combined)
}



