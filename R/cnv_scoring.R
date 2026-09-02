#' @title cnv_scoring
#'
#' @description
#' 
#' Functions that perform further processing from the CNV performing a final filtering and Confidence Score in a Single Tool Approach
#' 
#' @author Pedro Granjo
#' @date 02-09-2026
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
  if (!is.null(by)) {
    missing_by <- setdiff(by, colnames(df))
    if (length(missing_by) > 0) {
      stop("Grouping columns not found: ",
           paste(missing_by, collapse = ", "))
    }
  }
  
  # ── Validate cell and sample columns ──────────────────────────────────────
  if (!cell_col %in% colnames(df)) {
    stop("Missing cell column: ", cell_col)
  }
  if (!sample_col %in% colnames(df)) {
    stop("Missing sample column: ", sample_col)
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
run_cnv_locus_analysis <- function(
    df,
    by               = NULL,
    min_ovelap       = 0.75,
    sample_col,
    cell_col,
    overlap_method   = "reciprocal",
    parallel         = FALSE,
    n_cores          = 1L,
    removed_log_retur = FALSE) {
  
  # ── Step 1: Cluster CNV events ──────────────────────────────────────────────
  clustered <- cluster_cnv_events_by(
    df             = df,
    by             = by,
    overlap_method = overlap_method,
    min_ovelap     = min_ovelap,
    parallel       = parallel,
    n_cores        = n_cores
  )
  clustered_table_with_equiv_id <- clustered$results_id
  
  # ── Step 2: First summary ───────────────────────────────────────────────────
  summary_tbl <- summarise_cnv_loci(
    df         = clustered_table_with_equiv_id,
    by         = by,
    sample_col = sample_col,
    cell_col   = cell_col
  )
  
  summary_tbl <- summary_tbl %>%
    dplyr::rename(
      start         = locus_start,
      end           = locus_end,
      cnv_length    = locus_width,
      cnv_length_mb = locus_width_mb
    )
  
    if (removed_log_retur) {
      return(list(
        clustered_events  = clustered_table_with_equiv_id,
        remove_log        = clustered$removed_log,
        cnv_locus_summary = summary_tbl
      ))
    } else {
      return(list(
        clustered_events  = clustered_table_with_equiv_id,
        cnv_locus_summary = summary_tbl
      ))
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
    p_arm_permission     = 70,
    q_arm_permission     = 70,
    clustered_events     = NULL,
    whole_chr_permission = 60
) {
  
  required_cols <- c(
    "n_cells",
    "arm_class",
    "whole_chromosome_gain",
    "whole_chromosome_loss",
    "p_arm_pct",
    "q_arm_pct",
    "p_centromere_pct",
    "centromere_q_pct"
  )
  
  missing_cols <- setdiff(required_cols,
                          colnames(clustered_events))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ",
         paste(missing_cols, collapse = ", "),
         "\nDid you run add_arm_percentages() first?")
  }
  
  n_before <- nrow(clustered_events)
  
  out <- clustered_events %>%
    dplyr::mutate(
      max_chr_event = pmax(
        whole_chromosome_gain,
        whole_chromosome_loss,
        na.rm = TRUE
      ),
      pass_cells = dplyr::case_when(
        is.na(effective_threshold)     ~ FALSE,
        n_cells >= effective_threshold ~ TRUE,
        TRUE                           ~ FALSE
      ),
      
      pass_centromere = dplyr::case_when(
        
        # p_arm: always pass
        # no centromere concern
        arm_class == "p_arm" ~ TRUE,
        
        # q_arm: always pass
        # no centromere concern
        arm_class == "q_arm" ~ TRUE,
        
        # p_centromere_q: both arms must pass
        # compensates for arm length differences
        arm_class == "p_centromere_q" ~
          dplyr::coalesce(
            p_arm_pct >= p_arm_permission &
              q_arm_pct >= q_arm_permission & max_chr_event >= whole_chr_permission,
            FALSE
          ),
        
        # p_centromere: p arm only
        # centromere length removed
        arm_class == "p_centromere" ~
          dplyr::coalesce(
            p_centromere_pct >= p_arm_permission,
            FALSE
          ),
        
        # centromere_q: q arm only
        # centromere length removed
        arm_class == "centromere_q" ~
          dplyr::coalesce(
            centromere_q_pct >= q_arm_permission,
            FALSE
          ),
        
        # Any other arm class: pass
        TRUE ~ TRUE
      )
      
    ) %>%
    dplyr::filter(pass_cells,
                  pass_centromere) %>%
    dplyr::select(
      -dplyr::all_of(c("pass_cells",
                       "pass_centromere",
                       "effective_threshold", "max_chr_event"))
    )
  
  n_after <- nrow(out)
  
  message(sprintf(paste0(
    "Cell filter summary:\n",
    "  Input:                    %d rows\n",
    "  Retained:                 %d rows\n",
    "  Removed (total):          %d rows (%.1f%%)\n",
    "  P arm threshold:          %.0f%%\n",
    "  Q arm threshold:          %.0f%%"
  ),
  n_before,
  n_after,
  n_before - n_after,
  100 * (n_before - n_after) / n_before,
  p_arm_permission,
  q_arm_permission
  ))
  
  return(out)
}



standardise_cluster_boundaries <- function(
    clustered_events,
    chromosome_arms,
    sample_col = "sample",
    cell_col   = "cell_name") {
  
  cat("Standardising cluster boundaries...\n")
  
  n_before <- nrow(clustered_events)
  
  # Compute union boundaries per cluster
  cluster_boundaries <- clustered_events %>%
    dplyr::group_by(cnv_equiv_id) %>%
    dplyr::summarise(
      start_std         = min(.data$start, na.rm = TRUE),
      end_std           = max(.data$end,   na.rm = TRUE),
      genes_cluster     = paste(
        unique(trimws(unlist(
          strsplit(genes, ";")
        ))),
        collapse = ";"
      ),
      n_genes_cluster   = dplyr::n_distinct(
        trimws(unlist(strsplit(genes, ";")))
      ),
      genes_total_length_cluster    = mean(
        genes_total_length, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cnv_length_std    = end_std - start_std,
      cnv_length_mb_std = cnv_length_std / 1e6,
      
      genes_total_length_mb_cluster = 
        genes_total_length_cluster / 1e6,
        
      genes_coverage_pct_cluster    = round(
        100 * genes_total_length_cluster /
          pmax(cnv_length_std, 1), 1)
    )
  
  cat("Clusters standardised:",
      nrow(cluster_boundaries), "\n")
  
  # Check how many clusters had boundary changes
  changed <- clustered_events %>%
    dplyr::left_join(cluster_boundaries,
                     by = "cnv_equiv_id") %>%
    dplyr::summarise(
      n_start_changed = sum(start != start_std,
                            na.rm = TRUE),
      n_end_changed   = sum(end   != end_std,
                            na.rm = TRUE)
    )
  
  cat("Rows with start changed:",
      changed$n_start_changed, "\n")
  cat("Rows with end changed:  ",
      changed$n_end_changed, "\n")
  
  # Apply standardised boundaries
  result <- clustered_events %>%
    dplyr::left_join(cluster_boundaries,
                     by = "cnv_equiv_id") %>%
    dplyr::mutate(
      start          = start_std,
      end            = end_std,
      cnv_length     = cnv_length_std,
      cnv_length_mb  = cnv_length_mb_std,
      genes              = genes_cluster,
      n_genes            = n_genes_cluster,
      genes_total_length = genes_total_length_cluster,
      genes_total_length_mb = genes_total_length_mb_cluster,
      genes_coverage_pct = genes_coverage_pct_cluster
    ) %>%
    dplyr::select(-start_std, -end_std,
      -cnv_length_std, -cnv_length_mb_std,
      -genes_cluster, -n_genes_cluster,
      -genes_total_length_cluster,
      -genes_total_length_mb_cluster,
      -genes_coverage_pct_cluster)
  
  cat("Rows before:", n_before, "\n")
  cat("Rows after: ", nrow(result), "\n")
  
  # Recompute arm percentages
  # since boundaries changed
  cat("Recomputing arm percentages...\n")
  result <- add_arm_percentages(
    cnv_df          = result,
    chromosome_arms = chromosome_arms
  )
  
  return(result)
}




deduplicate_cnv_cells <- function(
    clustered_events,
    sample_col   = "sample",
    cell_col     = "cell_name",
    min_overlap  = 0.75) {
  
  cat("Deduplicating cells across clusters...\n")
  cat("Min overlap threshold:", min_overlap, "\n")
  
  n_before <- nrow(clustered_events)
  
  # Find cells appearing in multiple clusters
  # on same chr + cnv_state
  cell_cluster_map <- clustered_events %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(
        c(cell_col, sample_col,
          "chr", "cnv_state")
      ))
    ) %>%
    dplyr::summarise(
      n_clusters   = dplyr::n(),
      equiv_ids    = list(cnv_equiv_id),
      starts       = list(start),
      ends         = list(end),
      cnv_lengths  = list(cnv_length),
      .groups      = "drop"
    ) %>%
    dplyr::filter(n_clusters > 1)
  
  cat("Cells in multiple clusters:",
      nrow(cell_cluster_map), "\n")
  
  if (nrow(cell_cluster_map) == 0) {
    cat("No duplicates found\n")
    return(clustered_events)
  }
  
  # For each duplicated cell
  # find which cluster to keep
  cells_to_remove <- lapply(
    seq_len(nrow(cell_cluster_map)),
    function(i) {
      
      row      <- cell_cluster_map[i, ]
      ids      <- unlist(row$equiv_ids)
      starts   <- unlist(row$starts)
      ends     <- unlist(row$ends)
      lengths  <- unlist(row$cnv_lengths)
      n        <- length(ids)
      
      # Compute pairwise reciprocal overlap
      # between all cluster pairs for this cell
      keep_id <- ids[1]  # default: keep first
      
      if (n == 2) {
        
        # Overlap between two clusters
        overlap_start <- max(starts)
        overlap_end   <- min(ends)
        overlap_len   <- max(0,
                             overlap_end -
                               overlap_start)
        
        # Reciprocal overlap
        recip_overlap <- overlap_len /
          max(lengths)
        
        if (recip_overlap >= min_overlap) {
          # Same event — keep cluster
          # with largest span
          keep_id <- ids[which.max(lengths)]
          remove_ids <- ids[ids != keep_id]
          
          return(data.frame(
            cell      = row[[cell_col]],
            sample    = row[[sample_col]],
            chr       = row$chr,
            cnv_state = row$cnv_state,
            remove_id = remove_ids,
            keep_id   = keep_id,
            overlap   = recip_overlap,
            action    = "merged"
          ))
        } else {
          # Different events — keep both
          return(data.frame(
            cell      = row[[cell_col]],
            sample    = row[[sample_col]],
            chr       = row$chr,
            cnv_state = row$cnv_state,
            remove_id = NA_character_,
            keep_id   = NA_character_,
            overlap   = recip_overlap,
            action    = "kept_both"
          ))
        }
        
      } else {
        
        # More than 2 clusters
        # compute all pairwise overlaps
        # keep cluster with most total overlap
        overlap_scores <- sapply(
          seq_len(n), function(j) {
            other_starts <- starts[-j]
            other_ends   <- ends[-j]
            overlaps <- pmax(
              0,
              pmin(ends[j], other_ends) -
                pmax(starts[j], other_starts)
            )
            sum(overlaps) / lengths[j]
          }
        )
        
        keep_id    <- ids[which.max(overlap_scores)]
        remove_ids <- ids[ids != keep_id]
        
        return(data.frame(
          cell      = row[[cell_col]],
          sample    = row[[sample_col]],
          chr       = row$chr,
          cnv_state = row$cnv_state,
          remove_id = remove_ids,
          keep_id   = keep_id,
          overlap   = max(overlap_scores),
          action    = "merged_multiple"
        ))
      }
    }
  )
  
  dedup_log <- dplyr::bind_rows(cells_to_remove)
  
  # IDs to remove per cell
  remove_pairs <- dedup_log %>%
    dplyr::filter(!is.na(remove_id)) %>%
    dplyr::select(
      !!cell_col    := cell,
      !!sample_col  := sample,
      chr,
      cnv_state,
      remove_id
    )
  
  cat("Cluster assignments to remove:",
      nrow(remove_pairs), "\n")
  cat("Action summary:\n")
  print(table(dedup_log$action,
              useNA = "always"))
  
  if (nrow(remove_pairs) == 0) {
    cat("No cells removed\n")
    return(clustered_events)
  }
  
  # Remove duplicate assignments
  result <- clustered_events %>%
    dplyr::anti_join(
      remove_pairs %>%
        dplyr::rename(cnv_equiv_id = remove_id),
      by = c(cell_col, sample_col,
             "chr", "cnv_state",
             "cnv_equiv_id")
    )
  
  n_after <- nrow(result)
  
  cat("Rows before dedup:", n_before, "\n")
  cat("Rows after dedup: ", n_after,  "\n")
  cat("Rows removed:     ",
      n_before - n_after, "\n")
  
  attr(result, "dedup_log") <- dedup_log
  
  return(result)
}


add_arm_percentages <- function(cnv_df,
                                chromosome_arms) {
  
  chromosomes_with_cnv <- unique(cnv_df$chr)
  
  result_list <- lapply(
    chromosomes_with_cnv,
    function(x) {
      
      chr_subset <- cnv_df %>%
        dplyr::filter(chr == x)
      
      chr_arms <- chromosome_arms %>%
        dplyr::filter(chr == x)
      
      # ── Extract arm boundaries ─────────────────────────────────────────────
      p_arm <- chr_arms %>% dplyr::filter(arm == "p")
      q_arm <- chr_arms %>% dplyr::filter(arm == "q")
      cen   <- chr_arms %>% dplyr::filter(arm == "cen")
      
      p_start  <- if (nrow(p_arm) == 0) NA_real_ else
        as.numeric(p_arm$arm_start[1])
      p_end    <- if (nrow(p_arm) == 0) NA_real_ else
        as.numeric(p_arm$arm_end[1])
      p_length <- if (nrow(p_arm) == 0) NA_real_ else
        as.numeric(p_arm$arm_length[1])
      
      q_start  <- if (nrow(q_arm) == 0) NA_real_ else
        as.numeric(q_arm$arm_start[1])
      q_end    <- if (nrow(q_arm) == 0) NA_real_ else
        as.numeric(q_arm$arm_end[1])
      q_length <- if (nrow(q_arm) == 0) NA_real_ else
        as.numeric(q_arm$arm_length[1])
      
      cen_length <- if (nrow(cen) == 0) 0 else
        as.numeric(cen$arm_length[1])
      
      chr_subset %>%
        dplyr::mutate(
          
          # ── Overlap of CNV with p arm ──────────────────────────────────────
          overlap_p = pmax(0,
            pmin(end,   p_end)   -
            pmax(start, p_start)
          ),
          
          # ── Overlap of CNV with q arm ──────────────────────────────────────
          overlap_q = pmax(0,
            pmin(end,   q_end)   -
            pmax(start, q_start)
          ),
          
          # ── p_arm_pct ─────────────────────────────────────────────────────
          p_arm_pct = dplyr::case_when(
            arm_class == "p_centromere_q" &
              !is.na(p_length) & p_length > 0 ~
              round(overlap_p / p_length * 100, 2),
            TRUE ~ NA_real_
          ),
          
          # ── q_arm_pct ─────────────────────────────────────────────────────
          q_arm_pct = dplyr::case_when(
            arm_class == "p_centromere_q" &
              !is.na(q_length) & q_length > 0 ~
              round(overlap_q / q_length * 100, 2),
            TRUE ~ NA_real_
          ),
          
          # ── p_centromere_pct ──────────────────────────────────────────────
          p_centromere_pct = dplyr::case_when(
            arm_class == "p_centromere" &
              !is.na(p_length) & p_length > 0 ~
              round(
                pmax(0, overlap_p - cen_length) /
                  p_length * 100, 2),
            TRUE ~ NA_real_
          ),
          
          # ── centromere_q_pct ──────────────────────────────────────────────
          centromere_q_pct = dplyr::case_when(
            arm_class == "centromere_q" &
              !is.na(q_length) & q_length > 0 ~
              round(
                pmax(0, overlap_q - cen_length) /
                  q_length * 100, 2),
            TRUE ~ NA_real_
          )
        ) %>%
        dplyr::select(-overlap_p, -overlap_q)
    }
  )
  
  dplyr::bind_rows(result_list)
}





prepare_cnv_thresholds <- function(
    summary_df,
    clustered_events,
    by_union = NULL,
    cell_sizes = NULL,
    k, 
    sensitivity_floor_mb = 20,
    min_required_cells   = 3,
    round_fun             = ceiling
){
  
  if (is.null(cell_sizes)) {
    stop("Your dataframe with total cell numbers is not coming through")
  }
  
  if (missing(k) || is.null(k) || !is.numeric(k) || length(k) != 1L) {
    stop(
      "k must be supplied as a single numeric value.\n",
      "Derive it from one trusted calibration point:\n",
      "  k = known_good_threshold / sqrt(known_good_n_cells)"
    )
  }
  
  if (sensitivity_floor_mb <= 0) {
    stop("sensitivity_floor_mb must be positive.")
  }
  
  
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
    
    cell_sizes_clean <- cell_sizes %>%
      dplyr::select(
        dplyr::all_of(by_union),
        n_total_cells
      )
      
    merged_df <- clustered_events %>%
      dplyr::left_join(
        summary_df %>% dplyr::select(dplyr::all_of(join_cols), n_cells),
        by = join_cols
      ) %>%
      dplyr::left_join(cell_sizes_clean, by = by_union
      )
    
    if (anyNA(merged_df$n_total_cells)) {
      missing <- length(unique(merged_df[["cnv_equiv_id"]][is.na(merged_df$n_total_cells)]))
      stop(sprintf("No cell count found in %d groups: ", missing))
    }
  }
  
  merged_df <- merged_df %>%
    dplyr::filter(cnv_length_mb > sensitivity_floor_mb)
  
  if (nrow(merged_df) == 0L) {
    stop(sprintf(
      "No events remain above the sensitivity floor (%.1f Mb).\n",
      sensitivity_floor_mb
    ))
  }
  
  thresholded_df <- merged_df %>%
    dplyr::mutate(
      size_factor          = cnv_length_mb / sensitivity_floor_mb,
      group_threshold       = k * sqrt(n_total_cells),
      effective_threshold = dplyr::case_when(
      cnv_length_mb >= 110 ~ 1L,
      TRUE ~ round_fun(
        pmax(min_required_cells, group_threshold / size_factor)
      )
    )
    ) %>%
    dplyr::select(-size_factor, -group_threshold)
  
  
  message(sprintf(paste0(
    "Thresholds computed:\n",
    "  k:                     %.4f\n",
    "  sensitivity_floor_mb:  %.1f\n",
    "  min_required_cells:    %d\n",
    "  Loci scored:            %d\n",
    "  effective_threshold range: %d - %d"
  ),
  k,
  sensitivity_floor_mb,
  min_required_cells,
  nrow(thresholded_df),
  min(thresholded_df$effective_threshold, na.rm = TRUE),
  max(thresholded_df$effective_threshold, na.rm = TRUE)
  ))
  
  return(thresholded_df)
}



score_cnv_clusters <- function(
    summary_df,
    clustered_events,
    cell_sizes,
    by_union,
    chromosome_arms,
    k,
    n_cells              = NULL,
    min_required_cells    = 2L,
    p_arm_permission     = 60,
    q_arm_permission     = 60,  
    whole_chr_permission = 65,  
    sensitivity_floor_mb = 20,
    round_fun = ceiling,
    sample_col = "sample",
    cell_coll = "cell_name"
) {
  
  # ---- Step 1: prepare thresholds — merge only if needed ------------------
  thresholded_df <- prepare_cnv_thresholds(
    summary_df         = summary_df,
    clustered_events   = clustered_events,
    by_union           = by_union,
    cell_sizes         = cell_sizes,
    k                  = k,
    sensitivity_floor_mb = sensitivity_floor_mb,
    min_required_cells = min_required_cells,
    round_fun          = ceiling
  )
  
  # Add percentage of p and q for filtering based on total length
  thresholded_df <- add_arm_percentages(
    cnv_df          = thresholded_df,
    chromosome_arms = chromosome_arms
  )
  
  
  # ---- Step 2: filter -----------------------------------------------------
  filtered_df <- filter_cnv_loci(
    clustered_events     = thresholded_df,
    p_arm_permission     = p_arm_permission,
    q_arm_permission     = q_arm_permission,
    whole_chr_permission = whole_chr_permission
  )
  
  scored_std <- standardise_cluster_boundaries(
    clustered_events = filtered_df,
    chromosome_arms  = chromosome_arms,
    sample_col       = sample_col,
    cell_col         = cell_col
  )
}
