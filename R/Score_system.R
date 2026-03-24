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


#' Cluster CNV events by reciprocal overlap
#'
#' Groups CNV events into clusters based on genomic
#' reciprocal overlap and CNV state. Events are first grouped by chromosome
#' and CNV state, and then clustered using graph connectivity derived from
#' reciprocal overlaps. A sequence similarity can be achieve by assessing sequences within subclusters
#'
#' @param df Data frame containing CNV events.
#' @param min_reciprocal_overlap Minimum reciprocal overlap required between two events
#' to be considered connected. Default is 0.75.
refine_complete_subclusters <- function(x, min_reciprocal_overlap = 0.75) {
  
  if (nrow(x) <= 1) {
    x$subcluster_id <- 1L
    return(x)
  }
  
  gr <- GRanges(
    seqnames = x$chr,
    ranges = IRanges(start = x$start, end = x$end)
  )
  
  hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  keep <- q != s
  q <- q[keep]
  s <- s[keep]
  
  keep <- q < s
  q <- q[keep]
  s <- s[keep]
  
  if (length(q) == 0) {
    x$subcluster_id <- seq_len(nrow(x))
    return(x)
  }
  
  inter <- width(pintersect(gr[q], gr[s]))
  recip_ok <- (inter / width(gr[q]) >= min_reciprocal_overlap) &
    (inter / width(gr[s]) >= min_reciprocal_overlap)
  
  q <- q[recip_ok]
  s <- s[recip_ok]
  
  
  # adjacency matrix for strict pairwise overlap
  adj <- matrix(FALSE, nrow = nrow(x), ncol = nrow(x))
  diag(adj) <- TRUE
  
  if (length(q) > 0) {
    adj[cbind(q, s)] <- TRUE
    adj[cbind(s, q)] <- TRUE
  }
  
  ord <- order(x$start, x$end)
  subcluster <- integer(nrow(x))
  groups <- list()
  k <- 0L
  
  for (i in ord) {
    assigned <- FALSE
    
    if (length(groups) > 0) {
      for (g in seq_along(groups)) {
        members <- groups[[g]]
        if (all(adj[i, members])) {
          groups[[g]] <- c(members, i)
          subcluster[i] <- g
          assigned <- TRUE
        }
      }
    }
    
    if (!assigned) {
      k <- k + 1L
      groups[[k]] <- i
      subcluster[i] <- k
    }
  }
  
  x$subcluster_id <- subcluster
  x
}

#' Cluster CNV events within groups
#'
#' This function is applied separately to groups of CNV events
#' defined by one or more metadata columns.
#'
#' @param df Data frame containing CNV events.
#' @param by Character vector specifying grouping columns (e.g. dataset, embryo, cell type).
#' @param cluster_mode Clustering strategy passed to connected (graph connectivity derived from
#' reciprocal overlaps or complete if we want subgroups within the clusters connectivity)
#' @param min_reciprocal_overlap Minimum reciprocal overlap threshold.
#'
#' @details
#' The function splits the input data frame according to the grouping
#' variables and clusters CNV events independently within each subset.
#'
#' @return
#' Data frame with clustered CNV events including the `cluster_id` column.
cluster_cnv_events <- function(df, min_reciprocal_overlap = 0.75,
                               cluster_mode = c("connected", "complete")) {
  
  cluster_mode <- match.arg(cluster_mode)
  
  stopifnot(all(c("event_id","chr","start","end","cnv_state") %in% colnames(df)))
  
  df <- df %>% mutate(cluster_key = paste(chr, cnv_state, sep = "|"))
  split_df <- split(df, df$cluster_key)
  
  out <- vector("list", length(split_df))
  names(out) <- names(split_df)
  
  for (nm in names(split_df)) {
    x <- split_df[[nm]]
    
    if (nrow(x) == 1) {
      x$cluster_id <- paste0(nm, "|1")
      out[[nm]] <- x
      next
    }
    
    gr <- GRanges(
      seqnames = x$chr,
      ranges   = IRanges(start = x$start, end = x$end),
      event_id = x$event_id
    )
    
    hits <- findOverlaps(gr, gr, ignore.strand = TRUE)
    q <- queryHits(hits)
    s <- subjectHits(hits)
    
    keep <- q != s
    q <- q[keep]
    s <- s[keep]
    
    if (length(q) == 0) {
      x$cluster_id <- paste0(nm, "|", seq_len(nrow(x)))
      out[[nm]] <- x
      next
    }
    
    keep <- q < s
    q <- q[keep]
    s <- s[keep]
    
    inter <- width(pintersect(gr[q], gr[s]))
    recip_ok <- (inter / width(gr[q]) >= min_reciprocal_overlap) &
      (inter / width(gr[s]) >= min_reciprocal_overlap)
    
    q <- q[recip_ok]
    s <- s[recip_ok]
    
    if (length(q) == 0) {
      x$cluster_id <- paste0(nm, "|", seq_len(nrow(x)))
      out[[nm]] <- x
      next
    }
    
    # first pass: connected components
    g <- igraph::graph_from_edgelist(cbind(as.character(q), as.character(s)), directed = FALSE)
    comps <- igraph::components(g)$membership
    
    cluster_num <- seq_len(nrow(x))
    cluster_num[as.integer(names(comps))] <- comps
    x$cluster_id <- paste0(nm, "|", cluster_num)
    
    if (cluster_mode == "complete") {
      # refine inside each connected component
      x <- x %>%
        group_split(cluster_id, .keep = TRUE) %>%
        lapply(refine_complete_subclusters, min_reciprocal_overlap = min_reciprocal_overlap) %>%
        bind_rows()
      
      x <- x %>%
        mutate(cluster_id = paste0(cluster_id, "_", subcluster_id)) %>%
        select(-subcluster_id)
    }
    
    out[[nm]] <- x
  }
  
  bind_rows(out)
}




#' Cluster CNV events within groups
#'
#' Applies `cluster_cnv_events()` separately to groups of CNV events
#' defined by one or more metadata columns.
#'
#' @param df Data frame containing CNV events.
#' @param by Character vector specifying grouping columns (e.g. dataset,
#' sample, or patient).
#' @param cluster_mode Clustering strategy passed to `cluster_cnv_events()`.
#' @param min_reciprocal_overlap Minimum reciprocal overlap threshold.
#'
#' @details
#' The function splits the input data frame according to the grouping
#' variables and clusters CNV events independently within each subset.
#'
#' @return
#' Data frame with clustered CNV events including the `cluster_id` column.
cluster_cnv_events_by <- function(df, by = NULL, cluster_mode = c("connected", "complete"), min_reciprocal_overlap = 0.75) {
  cluster_mode <- match.arg(cluster_mode)
  required_cols <- c("event_id", "chr", "start", "end", "cnv_state")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.null(by) || length(by) == 0) {
    out <- cluster_cnv_events(df, min_reciprocal_overlap = min_reciprocal_overlap, cluster_mode = cluster_mode)
    return(out)
  }
  
  missing_by <- setdiff(by, colnames(df))
  if (length(missing_by) > 0) {
    stop("Grouping columns not found: ", paste(missing_by, collapse = ", "))
  }
  
  split_keys <- interaction(df[, by], drop = TRUE, sep = "___")
  split_list <- split(df, split_keys)
  
  out <- purrr::imap(split_list, function(subdf, nm) {
    res <- cluster_cnv_events(subdf, min_reciprocal_overlap = min_reciprocal_overlap, cluster_mode = cluster_mode)
    res
  }) %>%
    dplyr::bind_rows()
  out
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
  
  required_cols <- c("cluster_id", "chr", "cnv_state", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!sample_col %in% colnames(df)) {
    stop("Missing sample column: ", sample_col)
  }
  
  if (!cell_col %in% colnames(df)) {
    stop("Missing cell column: ", cell_col)
  }
  
  grouping_vars <- c(by, "cluster_id", "chr", "cnv_state")
  
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
  
  # Add mode-specific summaries only if mode exists
  if (mode_col %in% colnames(df)) {
    mode_summary <- df %>%
      dplyr::summarise(
        n_within = sum(.data[[mode_col]] == "within", na.rm = TRUE),
        n_across = sum(.data[[mode_col]] == "across", na.rm = TRUE),
        n_both = sum(.data[[mode_col]] == "both", na.rm = TRUE),
        cells_within = dplyr::n_distinct(.data[[cell_col]][.data[[mode_col]] == "within"]),
        cells_across = dplyr::n_distinct(.data[[cell_col]][.data[[mode_col]] == "across"]),
        cells_both = dplyr::n_distinct(.data[[cell_col]][.data[[mode_col]] == "both"]),
        samples_within = dplyr::n_distinct(.data[[sample_col]][.data[[mode_col]] == "within"]),
        samples_across = dplyr::n_distinct(.data[[sample_col]][.data[[mode_col]] == "across"]),
        samples_both = dplyr::n_distinct(.data[[sample_col]][.data[[mode_col]] == "both"]),
        .groups = "drop"
      )
    
    out <- out %>%
      dplyr::left_join(mode_summary, by = grouping_vars)
  }
  
  out %>%
    dplyr::arrange(dplyr::desc(n_cells), dplyr::desc(n_events))
}

#' Run CNV locus analysis workflow
#'
#' Performs CNV event clustering followed by locus-level summarisation.
#'
#' @param df Data frame containing CNV events.
#' @param by Optional grouping variables.
#' @param min_reciprocal_overlap Minimum reciprocal overlap threshold.
#' @param cluster_mode Clustering strategy ("connected" or "complete").
#' @param sample_col Column identifying samples.
#'
#' @return
#' A list containing:
#'
#' * `clustered_events` — CNV events with cluster assignments
#' * `cnv_locus_summary` — summarised CNV loci
#'
run_cnv_locus_analysis <- function(df, by = NULL, min_reciprocal_overlap = 0.75, cluster_mode = c("connected", "complete"),sample_col) {
  cluster_mode <- match.arg(cluster_mode)
  
  clustered <- cluster_cnv_events_by(
    df = df,
    by = by,
    min_reciprocal_overlap = min_reciprocal_overlap,
   cluster_mode = cluster_mode
  )

  summary_tbl <- summarise_cnv_loci(
    df = clustered,
    by = by,
    sample_col = sample_col
  )
  
  
  list(
    clustered_events = clustered,
    cnv_locus_summary = summary_tbl
  )
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
filter_cnv_loci_min_cells <- function(df,
                                      min_cells_keep = 5,
                                      threshold_df = NULL,
                                      group_col = NULL) {
  
  required_cols <- c(
    "n_cells",
    "arm_class",
    "whole_chromosome_gain",
    "whole_chromosome_loss"
  )
  
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- df
  
  # Optional warning if group_col is passed but no threshold table is used
  if (!is.null(group_col) && is.null(threshold_df)) {
    warning("group_col was provided but threshold_df is NULL; using global min_cells_keep.")
  }
  
  # If per-group thresholds are provided, join them
  if (!is.null(threshold_df)) {
    
    if (is.null(group_col)) {
      stop("group_col must be provided when threshold_df is used.")
    }
    
    if (!group_col %in% colnames(out)) {
      stop("group_col not found in df: ", group_col)
    }
    
    if (!group_col %in% colnames(threshold_df)) {
      stop("group_col not found in threshold_df: ", group_col)
    }
    
    if (!"min_cells_keep" %in% colnames(threshold_df)) {
      stop("threshold_df must contain column `min_cells_keep`.")
    }
    
    # Avoid duplicated min_cells_keep column if function is reused
    if ("min_cells_keep" %in% colnames(out)) {
      out <- out %>% dplyr::select(-min_cells_keep)
    }
    
    out <- out %>%
      dplyr::left_join(
        threshold_df %>%
          dplyr::select(dplyr::all_of(group_col), min_cells_keep),
        by = group_col
      )
    
    if (any(is.na(out$min_cells_keep))) {
      stop("Some rows have no min_cells_keep after joining threshold_df.")
    }
    
  } else {
    # Global threshold
    out <- out %>%
      dplyr::mutate(min_cells_keep = min_cells_keep)
  }
  
  out %>%
    dplyr::filter(
      n_cells >= min_cells_keep,
      !(arm_class == "p_centromere_q" &
          pmax(whole_chromosome_gain, whole_chromosome_loss, na.rm = TRUE) < 70)
    )
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
assign_cluster_confidence_binary <- function(
    summary_df,
    high_min_cells = 5,
    high_threshold_df = NULL,
    group_col = NULL,
    min_length_mb = 10,
    length_col = c("cnv_length_mb", "locus_width_mb"),
    keep_low = TRUE
){
  
  length_col <- match.arg(length_col)
  
  req_cols <- c("n_cells", length_col)
  missing_cols <- setdiff(req_cols, colnames(summary_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in summary_df: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- summary_df
  
  
  
  if (!is.null(high_threshold_df)) {
    if (is.null(group_col)) {
      stop("group_col must be provided when high_threshold_df is used")
    }
    
    if (!group_col %in% colnames(out)) {
      stop("group_col not found in summary_df: ", group_col)
    }
    
    if (!group_col %in% colnames(high_threshold_df)) {
      stop("group_col not found in high_threshold_df: ", group_col)
    }
    
    if (!"high_min_cells" %in% colnames(high_threshold_df)) {
      stop("high_threshold_df must contain column `high_min_cells`")
    }
    
  
    out <- out %>%
      dplyr::left_join(
        high_threshold_df %>% dplyr::select(dplyr::all_of(group_col), high_min_cells),
        by = group_col
      )
    
    if (any(is.na(out$high_min_cells))) {
      stop("Some rows have no high_min_cells after joining high_threshold_df")
    }
    
  } else {
    out <- out %>%
      dplyr::mutate(high_min_cells = high_min_cells)
  }
  
  out <- out %>%
    dplyr::mutate(
      pass_cells = n_cells >= high_min_cells,
      pass_length = .data[[length_col]] >= min_length_mb,
      confidence = dplyr::case_when(
        pass_cells & pass_length ~ "high",
        TRUE ~ "low"
      ),
      low_reason = dplyr::case_when(
        confidence == "high" ~ NA_character_,
        !pass_cells ~ "failed_cells",
        !pass_length ~ "failed_length",
        TRUE ~ "other"
      )
    )
  
  if (!keep_low) {
    out <- out %>% dplyr::filter(confidence == "high")
  }
  
  out %>%
    dplyr::mutate(confidence = factor(confidence, levels = c("high", "low")))
  
  
}

#' Determine CNV detection thresholds per cell type
#'
#' Calculates minimum cell thresholds required to retain CNV loci for
#' different cell types based on dataset size or user-defined rules.
#'
#' @param celltype_sizes Table containing cell types and their total
#' cell counts.
#' @param method Threshold calculation strategy.
#' @param fixed_value Fixed threshold used when `method = "fixed"`.
#' @param fraction Fraction of total cells used when `method = "fraction"`.
#' @param manual_thresholds Named vector of manually specified thresholds.
#' @param min_threshold Minimum allowed threshold. (Minimum cap)
#' @param max_threshold Maximum allowed threshold. (Maximum cap)
#' @param round_fun Rounding function used when computing thresholds.
#'
#' @return
#' Data frame with calculated `min_cells_keep` per cell type. We can use this function for filtering and also produced degree of confidence scores
resolve_celltype_thresholds <- function(celltype_sizes,
                                        method = c("fixed", "fraction", "median", "mean", "manual"),
                                        fixed_value = 5,
                                        fraction = 0.05,
                                        manual_thresholds = NULL,
                                        min_threshold = 3,
                                        max_threshold = 20,
                                        round_fun = ceiling) {
  
  method <- match.arg(method)
  
  if (is.vector(celltype_sizes) && !is.null(names(celltype_sizes))) {
    celltype_sizes <- data.frame(
      cell_type = names(celltype_sizes),
      n_total_cells = as.numeric(celltype_sizes),
      stringsAsFactors = FALSE
    )
  }
  
  required_cols <- c("cell_type", "n_total_cells")
  missing_cols <- setdiff(required_cols, colnames(celltype_sizes))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in celltype_sizes: ",
         paste(missing_cols, collapse = ", "))
  }
  
  out <- celltype_sizes
  
  if (method == "fixed") {
    out$min_cells_keep <- fixed_value
    
  } else if (method == "fraction") {
    raw_thr <- round_fun(out$n_total_cells * fraction)
    out$min_cells_keep <- pmin(max_threshold, pmax(min_threshold, raw_thr))
    
  } else if (method == "median") {
    thr <- round_fun(median(out$n_total_cells, na.rm = TRUE))
    out$min_cells_keep <- pmin(max_threshold, pmax(min_threshold, thr))
    
  } else if (method == "mean") {
    thr <- round_fun(mean(out$n_total_cells, na.rm = TRUE))
    out$min_cells_keep <- pmin(max_threshold, pmax(min_threshold, thr))
    
  } else if (method == "manual") {
    if (is.null(manual_thresholds)) {
      stop("manual_thresholds must be provided when method = 'manual'")
    }
    
    if (is.null(names(manual_thresholds))) {
      stop("manual_thresholds must be a named vector")
    }
    
    out$min_cells_keep <- manual_thresholds[out$cell_type]
    
    if (any(is.na(out$min_cells_keep))) {
      missing_types <- out$cell_type[is.na(out$min_cells_keep)]
      stop("Missing manual thresholds for: ", paste(missing_types, collapse = ", "))
    }
    
    out$min_cells_keep <- pmin(max_threshold, pmax(min_threshold, out$min_cells_keep))
  }
  
  out
}


#' Score and filter CNV clusters
#'
#' Integrates CNV cluster summaries with event-level data to filter loci
#' and assign confidence scores.
#'
#' @param summary_df CNV locus summary table.
#' @param clustered_events Event-level CNV table containing cluster IDs 
#' (with their respective subcluster depending on settings)
#' @param by_union Optional grouping variables used for joining. (eg., cell type, embryo)
#' @param min_cells_keep Minimum number of cells required to keep a locus.
#' @param min_threshold_df Optional per-group thresholds.
#' @param high_min_cells Minimum cells required for high confidence.
#' @param high_threshold_df Optional group-specific thresholds for
#' high confidence classification.
#' @param threshold_group_col Column defining threshold groups.
#' @param min_length_mb Minimum locus length required for high confidence.
#' @param length_col Column representing CNV length.
#' @param keep_low Whether to keep low-confidence clusters.
#'
#' @return
#' Data frame of scored CNV clusters.
score_cnv_clusters <- function(
    summary_df,
    clustered_events,
    by_union = NULL,
    min_cells_keep = 3,
    min_threshold_df = NULL,
    high_min_cells = 3,
    high_threshold_df = NULL,
    threshold_group_col = NULL,
    min_length_mb = 5,
    length_col = c("cnv_length_mb", "locus_width_mb"),
    keep_low = TRUE
){
  
  length_col <- match.arg(length_col)
  

  join_cols <- c(by_union, "cluster_id")
  
  missing_summary <- setdiff(join_cols, colnames(summary_df))
  if (length(missing_summary) > 0) {
    stop("Missing join columns in summary_df: ",
         paste(missing_summary, collapse = ", "))
  }
  
  missing_events <- setdiff(join_cols, colnames(clustered_events))
  if (length(missing_events) > 0) {
    stop("Missing join columns in clustered_events: ",
         paste(missing_events, collapse = ", "))
  }
  
  merged_df <- clustered_events %>%
    dplyr::left_join(
      summary_df %>% dplyr::select(dplyr::all_of(join_cols), n_cells),
      by = join_cols
    )

  filtered_df <- filter_cnv_loci_min_cells(
    df = merged_df,
    min_cells_keep = min_cells_keep,
    threshold_df = min_threshold_df,
    group_col = threshold_group_col
  )
  
  assign_cluster_confidence_binary(
    summary_df = filtered_df,
    high_min_cells = high_min_cells,
    high_threshold_df = high_threshold_df,
    group_col = threshold_group_col,
    min_length_mb = min_length_mb,
    length_col = length_col,
    keep_low = keep_low)
}


