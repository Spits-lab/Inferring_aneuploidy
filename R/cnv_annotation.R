#' @title cnv_annotation
#'
#' @description
#' 
#' Functions that helps in the annotation different CNVs with their chromosome information (arm, and coverage)
#' 
#' @author Pedro Granjo
#' @date 02-09-2026
#'



#' Compute overlap in base pairs between two intervals
#'
#' @param a_start,a_end Start and end coordinates of the first interval.
#' @param b_start,b_end Start and end coordinates of the second interval.
#'
#' @return Integer overlap length in base pairs.
overlap_bp <- function(a_start, a_end, b_start, b_end) {
  max(0, min(a_end, b_end) - max(a_start, b_start))
}



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
    cnv_df$arm_class <- character(0)
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
  
  if (chr_col != "chr" || start_col != "start" || end_col != "end") {
    main_df <- main_df |>
      dplyr::rename(
        chr   = dplyr::all_of(chr_col),
        start = dplyr::all_of(start_col),
        end   = dplyr::all_of(end_col)
      )
  }
  
  
  cnv_df <- classify_cnv_arms(main_df, chromosome_arms)
  
  missing_chr <- setdiff(cnv_df$chr, chromosome_arms$chr)
  
  if (length(missing_chr) > 0L) {
    message(
      "Skipping chromosomes not present in chromosome_arms: ",
      paste(missing_chr, collapse = ", ")
    )
  }
  
  cnv_df <- cnv_df |>
    dplyr::filter(chr %in% chromosome_arms$chr)
  
  whole_chr_info <- calculate_cnv_arm_percentages(cnv_df, chromosome_arms)
  return(whole_chr_info) 
}


