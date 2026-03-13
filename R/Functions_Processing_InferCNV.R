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
  
  runs <- list()
  
  for (ref in ref_dirs) {
    if( !is.null(base_dir)){
      ref_path <- file.path(base_dir, ref)
    } else{
      ref_path <- file.path(ref)
    }
    
    
    files <- list.files(
      ref_path,
      pattern = pattern,
      full.names = TRUE
    )
    
    
    
    if (length(files) == 0) {
      stop(sprintf("No run.final files found in %s", ref_path))
    }
    
    
    inferobj <- readRDS(files)
    
    
    runs[[ref]] <- inferobj
  }
  
  runs
}


#' Convert a wide expression matrix to long format
#'
#' Reshapes a gene-by-cell expression table into a long-format data frame with
#' one row per gene-cell combination.
#'
#' @param expr_df A data frame or matrix with genes as rows and cells as columns.
melt_expr_to_long <- function(expr_df) {
  expr_df |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(
      cols = -gene,
      names_to = "cell_name",
      values_to = "state_raw"
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
  
  if (nrow(merged) == 0 & nrow(merged) == nrow(long_df)) {
    stop("Merge produced empty table — gene names mismatch?")
  }

  merged
}


#' Discretize raw CNV scores into gain, loss, or neutral states
#'
#' Converts a continuous CNV score into a categorical CNV state using
#' mean-plus-minus-\code{k}-standard-deviation thresholds.
#'
#' @param df A data frame containing a numeric \code{state_raw} column.
#' @param k Numeric multiplier for the standard deviation cutoff. Default is 1.5.
discretize_cnv_state <- function(df, k =1.5) {
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
      discretize_cnv_state() |>
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

#' Compute reciprocal overlap between two intervals
#'
#' Calculates the reciprocal overlap between two genomic intervals as the
#' minimum overlap fraction relative to each interval length.
#'
#' @param start1,end1 Start and end coordinates for the first interval.
#' @param start2,end2 Start and end coordinates for the second interval.
#'
#' @return A numeric value between 0 and 1.
reciprocal_overlap <- function(start1, end1, start2, end2) {
  
  overlap <- max(0, min(end1, end2) - max(start1, start2))
  
  if (overlap == 0) {
    return(0)
  }
  
  len1 <- end1 - start1
  len2 <- end2 - start2
  
  min(overlap / len1, overlap / len2)
}


#' Assign equivalence groups to CNV events
#'
#' Groups CNV segments into equivalence classes within each cell, chromosome,
#' and CNV state based on reciprocal overlap.
#'
#' @param df A data frame containing CNV segments.
#' @param min_reciprocal_overlap Minimum reciprocal overlap required for two
#'   consecutive segments to belong to the same equivalence group.
assign_cnv_equivalence <- function(df, min_reciprocal_overlap = 0.5) {
  
  required <- c("cell_name", "chr", "cnv_state", "start", "end")
  if (!all(required %in% colnames(df))) {
    stop("Missing required columns")
  }
  
  df <- df |>
    arrange(cell_name, chr, cnv_state, start)
  
  df$cnv_equiv_id <- NA_integer_
  
  equiv_counter <- 0
  
  #This split in least for each cell, each chromossome and each state
  grouped <- df |>
    group_split(cell_name, chr, cnv_state)
  
  for (group in grouped) {
    
    current_start <- NULL
    current_end   <- NULL
    
    for (i in seq_len(nrow(group))) {
      
      row_idx <- which(df$cell_name == group$cell_name[i] &
                         df$chr       == group$chr[i] &
                         df$cnv_state == group$cnv_state[i] &
                         df$start     == group$start[i] &
                         df$end       == group$end[i])[1]
      
      if (is.null(current_start)) {
        # first CNV in this group
        equiv_counter <- equiv_counter + 1
        current_start <- group$start[i]
        current_end   <- group$end[i]
      } else {
        
        ro <- reciprocal_overlap(
          current_start, current_end,
          group$start[i], group$end[i]
        )
        
        if (ro < min_reciprocal_overlap) {
          # start new equivalence group
          equiv_counter <- equiv_counter + 1
          current_start <- group$start[i]
          current_end   <- group$end[i]
        } else {
          # extend current equivalence interval
          current_end <- max(current_end, group$end[i])
        }
      }
      
      df$cnv_equiv_id[row_idx] <- equiv_counter
    }
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
    group_by(cell_name, cnv_equiv_id) |>
    summarise(
      chr          = dplyr::first(chr),
      cnv_state    = dplyr::first(cnv_state),
      start        = min(start),
      end          = max(end),
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
  
  stopifnot("n_references" %in% colnames(cnv_events))
  
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
run_fast_cnv_pipeline <- function(
    gene_level_df,
    max_gap = 100000,
    min_reciprocal_overlap = 0.5,
    min_references = 2
) {
  
  message("→ Collapsing genes to segments")
  segments <- collapse_genes_to_cnv_segments(gene_level_df)
  dim(segments)
 
  
  message("→ Merging nearby CNVs")
  merged <- merge_nearby_regions(segments, max_gap)
  
  message("→ Assigning CNV equivalence")
  equiv <- assign_cnv_equivalence(
    merged,
    min_reciprocal_overlap
  )
  
  cnv_events <- summarize_cnv_support(equiv)
  
  n_refs <- length(unique(equiv$reference))
  
  message("→ Filtering CNVs by reference support")
  
  supported_events <- filter_cnv_events(
    cnv_events,
    min_references = 2
  )
  
  
  list(
    cnvs_all        = equiv,
    cnvs_filtered   = cnv_events,
    cnv_support_tbl = supported_events
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

#' Remove a prefix from cell names
#'
#' Removes a regular expression prefix from cell identifiers.
#'
#' @param cell_names Character vector of cell names.
#' @param prefix Regular expression prefix to remove. Default is
#'   \code{"^cell_"}.
#'
#' @return Character vector of cleaned cell names.
clean_cell_name <- function(cell_names, prefix = "^cell_") {
  sub(prefix, "", cell_names)
}



#' Add parsed metadata and CNV length columns
#'
#' Adds cleaned cell names, embryo labels, stage labels, and CNV length columns
#' to a CNV data frame.
#'
#' @param df A data frame containing at least \code{cell_name}, \code{start},
#'   and \code{end}.
#'
#' @return The input data frame with added columns:
#' \code{cell_name_new}, \code{embryo}, \code{stage}, \code{cnv_length}, and
#' \code{cnv_length_mb}.
add_additional_columns <- function(df){
  # Remove the "cell_" prefix first
  df$cell_name_new <- sub("^cell_", "", df$cell_name)
  
  # Extract embryo (everything except the last dot and cell number)
  df$embryo <- sub("^(.+)\\.[0-9]+$", "\\1", df$cell_name_new)
  
  # Extract stage (first part before the first dot)
  df$stage <- sub("^([^\\.]+)\\..*$", "\\1", df$cell_name_new)
  

  df <- df %>%
    mutate(
      cnv_length = end - start + 1,
      cnv_length_mb = cnv_length / 1e6
    )
  
  
  return(df)
}



########################################################
#Stats Function

#' Add parsed metadata and CNV length columns
#'
#' Adds cleaned cell names, embryo labels, stage labels, and CNV length columns
#' to a CNV data frame.
#'
#' @param df A data frame containing at least "cell_name", "start", "end"
#'
count_stage_embryo <- function(df, count_name = "count") {
  
  df %>%
    dplyr::distinct(stage, embryo, cell_name) %>%
    dplyr::group_by(stage, embryo) %>%
    dplyr::summarise(!!count_name := dplyr::n(), .groups = "drop")
}

#' Count CNV-positive cells by stage and embryo
#'
#' Extracts stage, embryo, and cell information from a CNV table and counts
#' unique cells per stage-embryo combination.
#'
#' @param df A CNV data frame containing "stage","embryo", "cell_name"
#'
#' @return A summary data frame with a cnv_count column.
#'
get_cnv_counts <- function(df) {
  
  cnv_df <- df %>%
    dplyr::select(stage, embryo, cell_name)
  
  count_stage_embryo(cnv_df, count_name = "cnv_count")
}

#' Count Seurat cells by stage and embryo
#'
#' Parses Seurat cell names into embryo and stage labels and counts cells per
#' stage-embryo combination.
#'
#' @param seurat_obj A Seurat object.
#'
#' @return A summary data frame with a seurat_count column.
get_seurat_counts <- function(seurat_obj) {
  
  seurat_cells <- data.frame(
    cell_name = colnames(seurat_obj),
    stringsAsFactors = FALSE
  )
  
  seurat_cells$cell_name_new <- sub("^cell_", "", seurat_cells$cell_name)
  seurat_cells$embryo <- sub("^(.+)\\.[0-9]+$", "\\1", seurat_cells$cell_name_new)
  seurat_cells$stage <- sub("^([^\\.]+)\\..*$", "\\1", seurat_cells$cell_name_new)
  
  count_stage_embryo(seurat_cells, count_name = "seurat_count")
}

#' Merge CNV and Seurat cell count summaries
#'
#' Combines embryo-stage count summaries from CNV data and Seurat data.
#'
#' @param cnv_counts A data frame of CNV-based counts.
#' @param seurat_counts A data frame of Seurat-based counts.
#' @param drop_missing_cnv Logical; if TRUE, rows with missing cnv_count are removed
#'
#' @return A merged count table.
merge_counts <- function(cnv_counts, seurat_counts, drop_missing_cnv = TRUE) {
  
  merged <- dplyr::full_join(
    cnv_counts,
    seurat_counts,
    by = c("stage", "embryo")
  )
  
  if (drop_missing_cnv) {
    merged <- merged[!is.na(merged$cnv_count), ]
  }
  
  return(merged)
}


#' Summarise counts by stage
#'
#' Aggregates embryo-level CNV and Seurat counts into stage-level totals.
#'
#' @param merged_counts A merged count data frame containing stage,cnv_count, seurat_count columns.
#'
#' @return A stage-level summary data frame.
summarise_by_stage <- function(merged_counts) {
  
  merged_counts %>%
    dplyr::group_by(stage) %>%
    dplyr::summarise(
      cnv_total = sum(cnv_count, na.rm = TRUE),
      seurat_total = sum(seurat_count, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Build embryo-level and stage-level CNV summaries
#'
#' Generates embryo-level and stage-level summaries comparing CNV-positive cells
#' with all Seurat cells.
#'
#' @param dfs A CNV data frame.
#' @param seurat_obj A Seurat object.
make_embryo_stage_summary <- function(dfs, seurat_obj) {
  
  cnv_counts <- get_cnv_counts(dfs)

  seurat_counts <- get_seurat_counts(seurat_obj)

  merged_counts <- merge_counts(cnv_counts, seurat_counts)

  stage_summary <- summarise_by_stage(merged_counts)
  
  
  return(list(
    embryo_level = merged_counts,
    stage_level  = stage_summary
  ))
}



#' Classify a CNV by chromosome arm overlap
#'
#' Assigns a CNV event to an arm-level category based on overlap with
#' chromosome arm annotations.
#'
#' @param cnv_row A single-row data frame representing one CNV.
#' @param chromosome_arms A data frame describing chromosome arm intervals.
classify_single_cnv <- function(cnv_row, chromosome_arms) {
  
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
  
  cnv_df$arm_class <- vapply(
    seq_len(nrow(cnv_df)),
    function(i) classify_single_cnv(cnv_df[i, ], chromosome_arms),
    character(1)
  )
  
  return(cnv_df)
}

#' Add chromosome arm classification to a CNV table
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
  
  classify_cnv_arms(main_df, chromosome_arms)
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



