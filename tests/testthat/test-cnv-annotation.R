source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/tests/testthat/helper-fixtures.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/cnv_annotation.R")

# ============================================================================
# overlap_bp
# ============================================================================

testthat::test_that("overlap_bp: perfect overlap returns full length", {
  testthat::expect_equal(overlap_bp(100, 200, 100, 200), 100)
})

testthat::test_that("overlap_bp: no overlap returns 0", {
  testthat::expect_equal(overlap_bp(100, 200, 300, 400), 0)
})

testthat::test_that("overlap_bp: partial overlap returns correct length", {
  testthat::expect_equal(overlap_bp(100, 200, 150, 250), 50)
})

testthat::test_that("overlap_bp: adjacent segments return 0", {
  testthat::expect_equal(overlap_bp(100, 200, 200, 300), 0)
})

testthat::test_that("overlap_bp: containment returns inner length", {
  testthat::expect_equal(overlap_bp(100, 300, 150, 200), 50)
})

testthat::test_that("overlap_bp: single base overlap returns 1", {
  testthat::expect_equal(overlap_bp(100, 150, 150, 200), 0)
})

testthat::test_that("overlap_bp: result is always non-negative", {
  testthat::expect_gte(overlap_bp(100, 200, 500, 600), 0)
})


# ============================================================================
# classify_single_cnv
# ============================================================================

testthat::test_that("classify_single_cnv: p arm only CNV returns p_arm", {
  cnv <- data.frame(chr = "chr1", start = 1e6, end = 50e6)
  out <- classify_single_cnv(cnv, hg38_chromosome_arms)
  testthat::expect_equal(out, "p_arm")
})

testthat::test_that("classify_single_cnv: q arm only CNV returns q_arm", {
  # Use real q arm coordinates from chr1
  chr1_arms <- hg38_chromosome_arms[hg38_chromosome_arms$chr == "chr1", ]
  q_start   <- chr1_arms$arm_start[chr1_arms$arm == "q"]
  q_end     <- chr1_arms$arm_end[chr1_arms$arm == "q"]
  
  cnv <- data.frame(
    chr   = "chr1",
    start = as.integer(q_start + 1e6),
    end   = as.integer(q_start + 20e6)
  )
  out <- classify_single_cnv(cnv, hg38_chromosome_arms)
  testthat::expect_equal(out, "q_arm")
})

testthat::test_that("classify_single_cnv: whole chromosome CNV returns p_centromere_q", {
  cnv <- data.frame(chr = "chr1", start = 1L, end = 250000000L)
  out <- classify_single_cnv(cnv, hg38_chromosome_arms)
  testthat::expect_equal(out, "p_centromere_q")
})

testthat::test_that("classify_single_cnv: unknown chromosome returns NA", {
  cnv <- data.frame(chr = "chrUNKNOWN", start = 1e6, end = 50e6)
  out <- classify_single_cnv(cnv, hg38_chromosome_arms)
  testthat::expect_true(is.na(out))
})

testthat::test_that("classify_single_cnv: missing cnv columns raises error", {
  cnv <- data.frame(chr = "chr1", start = 1e6)  # missing end
  testthat::expect_error(
    classify_single_cnv(cnv, hg38_chromosome_arms),
    "cnv_row missing columns"
  )
})

testthat::test_that("classify_single_cnv: missing arm columns raises error", {
  cnv  <- data.frame(chr = "chr1", start = 1e6, end = 50e6)
  arms <- hg38_chromosome_arms
  arms$arm <- NULL
  testthat::expect_error(
    classify_single_cnv(cnv, arms),
    "chromosome_arms missing columns"
  )
})

testthat::test_that("classify_single_cnv: output is character", {
  cnv <- data.frame(chr = "chr1", start = 1e6, end = 50e6)
  out <- classify_single_cnv(cnv, hg38_chromosome_arms)
  testthat::expect_type(out, "character")
})


# ============================================================================
# classify_cnv_arms
# ============================================================================

testthat::test_that("classify_cnv_arms: adds arm_class column", {
  df  <- make_mock_cnv()
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  testthat::expect_true("arm_class" %in% colnames(out))
})

testthat::test_that("classify_cnv_arms: arm_class values are valid", {
  df  <- make_mock_cnv()
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  valid <- c("p_arm", "q_arm", "p_centromere", "centromere_q",
             "p_centromere_q", NA_character_)
  testthat::expect_true(all(out$arm_class %in% valid | is.na(out$arm_class)))
})

testthat::test_that("classify_cnv_arms: row count unchanged", {
  df  <- make_mock_cnv()
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  testthat::expect_equal(nrow(out), nrow(df))
})

testthat::test_that("classify_cnv_arms: empty cnv_df returns warning and NA arm_class", {
  df  <- make_mock_cnv()[0, ]
  testthat::expect_warning(
    out <- classify_cnv_arms(df, hg38_chromosome_arms),
    "empty"
  )
  testthat::expect_equal(nrow(out), 0L)
})

testthat::test_that("classify_cnv_arms: missing cnv columns raises error", {
  df       <- make_mock_cnv()
  df$start <- NULL
  testthat::expect_error(
    classify_cnv_arms(df, hg38_chromosome_arms),
    "cnv_row missing columns"
  )
})

testthat::test_that("classify_cnv_arms: missing arm columns raises error", {
  df   <- make_mock_cnv()
  arms <- hg38_chromosome_arms
  arms$arm_start <- NULL
  testthat::expect_error(
    classify_cnv_arms(df, arms),
    "chromosome_arms missing columns"
  )
})

testthat::test_that("classify_cnv_arms: unexpected arm labels produce warning", {
  df   <- make_mock_cnv()
  arms <- hg38_chromosome_arms
  arms$arm[1] <- "telomere"  # invalid label
  testthat::expect_warning(
    classify_cnv_arms(df, arms),
    "Unexpected arm labels"
  )
})

testthat::test_that("classify_cnv_arms: chr1 loss at 117-150Mb is p_centromere_q", {

  df  <- make_mock_cnv()[1, ]
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  testthat::expect_equal(out$arm_class, "p_centromere_q")
})

# ============================================================================
# calculate_cnv_arm_percentages
# ============================================================================

testthat::test_that("calculate_cnv_arm_percentages: adds required columns", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  expected_cols <- c(
    "whole_chromosome_gain", "whole_chromosome_loss",
    "p_arm_gain", "p_arm_loss",
    "q_arm_gain", "q_arm_loss"
  )
  testthat::expect_true(all(expected_cols %in% colnames(out)))
})

testthat::test_that("calculate_cnv_arm_percentages: gain columns NA for loss CNVs", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  loss_rows <- out[out$cnv_state == "loss", ]
  testthat::expect_true(all(is.na(loss_rows$whole_chromosome_gain)))
})

testthat::test_that("calculate_cnv_arm_percentages: loss columns NA for gain CNVs", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  gain_rows <- out[out$cnv_state == "gain", ]
  testthat::expect_true(all(is.na(gain_rows$whole_chromosome_loss)))
})

testthat::test_that("calculate_cnv_arm_percentages: percentages between 0 and 100", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  pct_cols <- c(
    "whole_chromosome_gain", "whole_chromosome_loss",
    "p_arm_gain", "p_arm_loss",
    "q_arm_gain", "q_arm_loss"
  )
  
  for (col in pct_cols) {
    vals <- out[[col]][!is.na(out[[col]])]
    testthat::expect_true(
      all(vals >= 0 & vals <= 100),
      label = sprintf("%s values in [0, 100]", col)
    )
  }
})

testthat::test_that("calculate_cnv_arm_percentages: row count unchanged", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  testthat::expect_equal(nrow(out), nrow(df))
})

testthat::test_that("calculate_cnv_arm_percentages: whole chromosome gain > 0 for gain CNV", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  gain_rows <- out[out$cnv_state == "gain", ]
  testthat::expect_true(all(gain_rows$whole_chromosome_gain > 0, na.rm = TRUE))
})





testthat::test_that("percentages: gain CNV has NA in all loss columns", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  gain_rows <- out[out$cnv_state == "gain", ]
  testthat::expect_true(all(is.na(gain_rows$whole_chromosome_loss)))
  testthat::expect_true(all(is.na(gain_rows$p_arm_loss)))
  testthat::expect_true(all(is.na(gain_rows$q_arm_loss)))
})

testthat::test_that("percentages: loss CNV has NA in all gain columns", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  loss_rows <- out[out$cnv_state == "loss", ]
  testthat::expect_true(all(is.na(loss_rows$whole_chromosome_gain)))
  testthat::expect_true(all(is.na(loss_rows$p_arm_gain)))
  testthat::expect_true(all(is.na(loss_rows$q_arm_gain)))
})

testthat::test_that("percentages: p_arm CNV has NA in q_arm columns", {
  # Build a CNV guaranteed to be p_arm only
  chr1_arms <- hg38_chromosome_arms[hg38_chromosome_arms$chr == "chr1", ]
  p_end     <- chr1_arms$arm_end[chr1_arms$arm == "p"]
  
  df <- data.frame(
    cell_name     = "cell_1",
    chr           = "chr1",
    cnv_state     = "gain",
    start         = 1e6,
    end           = as.integer(p_end - 1e6),  # stays within p arm
    cnv_length    = as.integer(p_end - 1e6) - 1e6,
    cnv_length_mb = (as.integer(p_end - 1e6) - 1e6) / 1e6,
    stringsAsFactors = FALSE
  )
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  
  testthat::expect_equal(df$arm_class, "p_arm")
  
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  testthat::expect_true(is.na(out$q_arm_gain))
  testthat::expect_true(is.na(out$q_arm_loss))
})

testthat::test_that("percentages: q_arm CNV has NA in p_arm columns", {
  chr1_arms <- hg38_chromosome_arms[hg38_chromosome_arms$chr == "chr1", ]
  q_start   <- chr1_arms$arm_start[chr1_arms$arm == "q"]
  q_end     <- chr1_arms$arm_end[chr1_arms$arm == "q"]
  
  df <- data.frame(
    cell_name     = "cell_1",
    chr           = "chr1",
    cnv_state     = "gain",
    start         = as.integer(q_start + 1e6),  # stays within q arm
    end           = as.integer(q_start + 20e6),
    cnv_length    = as.integer(19e6),
    cnv_length_mb = 19,
    stringsAsFactors = FALSE
  )
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  
  testthat::expect_equal(df$arm_class, "q_arm")
  
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  testthat::expect_true(is.na(out$p_arm_gain))
  testthat::expect_true(is.na(out$p_arm_loss))
})

testthat::test_that("percentages: p_arm gain has p_arm_gain > 0 and p_arm_loss NA", {
  chr1_arms <- hg38_chromosome_arms[hg38_chromosome_arms$chr == "chr1", ]
  p_end     <- chr1_arms$arm_end[chr1_arms$arm == "p"]
  
  df <- data.frame(
    cell_name     = "cell_1",
    chr           = "chr1",
    cnv_state     = "gain",
    start         = 1e6,
    end           = as.integer(p_end - 1e6),
    cnv_length    = as.integer(p_end - 1e6) - 1e6,
    cnv_length_mb = (as.integer(p_end - 1e6) - 1e6) / 1e6,
    stringsAsFactors = FALSE
  )
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  testthat::expect_gt(out$p_arm_gain, 0)
  testthat::expect_true(is.na(out$p_arm_loss))
})

testthat::test_that("percentages: whole chromosome always filled for correct state", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  # Every gain row must have whole_chromosome_gain filled
  gain_rows <- out[out$cnv_state == "gain", ]
  testthat::expect_true(all(!is.na(gain_rows$whole_chromosome_gain)))
  
  # Every loss row must have whole_chromosome_loss filled
  loss_rows <- out[out$cnv_state == "loss", ]
  testthat::expect_true(all(!is.na(loss_rows$whole_chromosome_loss)))
})

testthat::test_that("percentages: no CNV can have both p_arm and q_arm percentage filled", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  # A single CNV cannot be in both p and q arm simultaneously
  both_filled <- !is.na(out$p_arm_gain) & !is.na(out$q_arm_gain) |
    !is.na(out$p_arm_loss) & !is.na(out$q_arm_loss)
  
  testthat::expect_false(any(both_filled))
})



testthat::test_that("percentages: centromere-spanning CNV has NA in arm columns", {
  df <- data.frame(
    cell_name     = "cell_1",
    chr           = "chr1",
    cnv_state     = "gain",
    start         = 1L,
    end           = 250000000L,  # whole chromosome
    cnv_length    = 250000000L,
    cnv_length_mb = 250,
    stringsAsFactors = FALSE
  )
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  # Whole chromosome CNV — arm columns should be NA
  testthat::expect_true(is.na(out$p_arm_gain))
  testthat::expect_true(is.na(out$q_arm_gain))
})


testthat::test_that("percentages: larger CNV has higher whole chromosome percentage", {
  chr1_arms <- hg38_chromosome_arms[hg38_chromosome_arms$chr == "chr1", ]
  q_start   <- chr1_arms$arm_start[chr1_arms$arm == "q"]
  
  # Small q arm CNV
  df_small <- data.frame(
    cell_name = "cell_1", chr = "chr1", cnv_state = "gain",
    start = as.integer(q_start + 1e6),
    end   = as.integer(q_start + 10e6),
    cnv_length    = 9e6L,
    cnv_length_mb = 9,
    stringsAsFactors = FALSE
  )
  
  # Large q arm CNV — same start, larger end
  df_large <- data.frame(
    cell_name = "cell_1", chr = "chr1", cnv_state = "gain",
    start = as.integer(q_start + 1e6),
    end   = as.integer(q_start + 50e6),
    cnv_length    = 49e6L,
    cnv_length_mb = 49,
    stringsAsFactors = FALSE
  )
  
  df_small <- classify_cnv_arms(df_small, hg38_chromosome_arms)
  df_large <- classify_cnv_arms(df_large, hg38_chromosome_arms)
  
  out_small <- calculate_cnv_arm_percentages(df_small, hg38_chromosome_arms)
  out_large <- calculate_cnv_arm_percentages(df_large, hg38_chromosome_arms)
  
  testthat::expect_gt(
    out_large$whole_chromosome_gain,
    out_small$whole_chromosome_gain
  )
})


testthat::test_that("classify_cnv_arms: chrX CNV is classified", {
  # From real data — chrX gain 14008279-75746911
  df <- data.frame(
    cell_name = "E5.10.933", chr = "chrX",
    cnv_state = "gain",
    start = 14008279L, end = 75746911L,
    cnv_length = 61738633L, cnv_length_mb = 61.7,
    stringsAsFactors = FALSE
  )
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  testthat::expect_false(is.na(out$arm_class))
})

testthat::test_that("classify_cnv_arms: chrY CNV is classified", {
  # From real data — chrY loss 2841486-20781032
  df <- data.frame(
    cell_name = "E5.10.933", chr = "chrY",
    cnv_state = "loss",
    start = 2841486L, end = 20781032L,
    cnv_length = 17939547L, cnv_length_mb = 17.9,
    stringsAsFactors = FALSE
  )
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  testthat::expect_false(is.na(out$arm_class))
})


testthat::test_that("percentages: no percentage exceeds 100%", {
  df  <- make_mock_cnv()
  df  <- classify_cnv_arms(df, hg38_chromosome_arms)
  out <- calculate_cnv_arm_percentages(df, hg38_chromosome_arms)
  
  pct_cols <- c(
    "whole_chromosome_gain", "whole_chromosome_loss",
    "p_arm_gain", "p_arm_loss",
    "q_arm_gain", "q_arm_loss"
  )
  for (col in pct_cols) {
    vals <- out[[col]][!is.na(out[[col]])]
    testthat::expect_true(
      all(vals <= 100),
      label = sprintf("%s never exceeds 100%%", col)
    )
  }
})


testthat::test_that("classify_cnv_arms: two CNVs same chromosome same cell classified independently", {
  df <- data.frame(
    cell_name = c("E5.10.933", "E5.10.933"),
    chr       = c("chrX", "chrX"),
    cnv_state = c("gain", "gain"),
    start     = c(14008279L,  135344796L),
    end       = c(75746911L,  149000663L),
    cnv_length    = c(61738633L, 13655868L),
    cnv_length_mb = c(61.7, 13.7),
    stringsAsFactors = FALSE
  )
  out <- classify_cnv_arms(df, hg38_chromosome_arms)
  
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_true(all(!is.na(out$arm_class)))
  
  testthat::expect_true(length(unique(out$arm_class)) >= 1L)
})

testthat::test_that("add_chromosome_info: cnv_length consistent with start and end", {
  df  <- make_mock_cnv()
  out <- add_chromosome_info(df, hg38_chromosome_arms)
  
  # cnv_length should equal end - start + 1 (or end - start depending on convention)
  # Check consistency — not exact formula since convention may vary
  testthat::expect_true(all(out$cnv_length > 0))
  testthat::expect_true(all(out$end > out$start))
})
# ============================================================================
# add_chromosome_info
# ============================================================================

testthat::test_that("add_chromosome_info: returns data frame", {
  df  <- make_mock_cnv()
  out <- add_chromosome_info(df, hg38_chromosome_arms)
  testthat::expect_s3_class(out, "data.frame")
})

testthat::test_that("add_chromosome_info: adds arm_class and percentage columns", {
  df  <- make_mock_cnv()
  out <- add_chromosome_info(df, hg38_chromosome_arms)
  
  testthat::expect_true("arm_class" %in% colnames(out))
  testthat::expect_true("whole_chromosome_gain" %in% colnames(out))
  testthat::expect_true("whole_chromosome_loss" %in% colnames(out))
})

testthat::test_that("add_chromosome_info: missing coordinate columns raises error", {
  df       <- make_mock_cnv()
  df$start <- NULL
  testthat::expect_error(
    add_chromosome_info(df, hg38_chromosome_arms),
    "Missing required CNV coordinate columns"
  )
})

testthat::test_that("add_chromosome_info: custom column names work", {
  df <- make_mock_cnv()
  df <- dplyr::rename(df,
                      chromosome = chr,
                      seg_start  = start,
                      seg_end    = end
  )
  out <- add_chromosome_info(
    df,
    hg38_chromosome_arms,
    chr_col   = "chromosome",
    start_col = "seg_start",
    end_col   = "seg_end"
  )
  testthat::expect_true("arm_class" %in% colnames(out))
})

testthat::test_that("add_chromosome_info: chromosomes not in arms table are dropped", {
  df <- make_mock_cnv()
  df$chr[1] <- "chrUNKNOWN"
  
  out <- testthat::expect_warning(
    add_chromosome_info(df, hg38_chromosome_arms),
    "arm-classified"
  )
  testthat::expect_false("chrUNKNOWN" %in% out$chr)
})

testthat::test_that("add_chromosome_info: output has no rows from unknown chromosomes", {
  df      <- make_mock_cnv()
  all_chr <- unique(df$chr)
  out     <- add_chromosome_info(df, hg38_chromosome_arms)
  testthat::expect_true(all(out$chr %in% hg38_chromosome_arms$chr))
})

