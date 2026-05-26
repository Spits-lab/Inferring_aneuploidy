source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/tests/testthat/helper-fixtures.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/cnv_scoring.R")


testthat::test_that("make_tier_definitions: fractions mode returns correct structure", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    base_fraction = 0.05,
    step          = 0.03,
    mode          = "fractions"
  )
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_true("tier"          %in% colnames(out))
  testthat::expect_true("boundaries_mb" %in% colnames(out))
  testthat::expect_true("fraction"      %in% colnames(out))
  testthat::expect_false("n_cells"      %in% colnames(out))
})

testthat::test_that("make_tier_definitions: number mode returns correct structure", {
  out <- make_tier_definitions(
    boundaries_mb = c(25,50),
    n_cells       = c(3L, 8L),
    mode          = "number"
  )
  testthat::expect_true("n_cells"  %in% colnames(out))
  testthat::expect_false("fraction" %in% colnames(out))
  testthat::expect_equal(nrow(out), 2L)
})

testthat::test_that("make_tier_definitions: boundaries sorted descending automatically", {
  out <- make_tier_definitions(
    boundaries_mb = c(25, 50),
    mode          = "fractions"
  )
  testthat::expect_true(all(diff(out$boundaries_mb) < 0))
})

testthat::test_that("make_tier_definitions: tier names follow tier_N convention", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    mode          = "fractions"
  )
  testthat::expect_equal(out$tier, c("tier_1", "tier_2"))
})

testthat::test_that("make_tier_definitions: auto fractions increase with step", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    base_fraction = 0.05,
    step          = 0.03,
    mode          = "fractions"
  )
  testthat::expect_equal(out$fraction[1], 0.05)
  testthat::expect_equal(out$fraction[2], 0.08)
})

testthat::test_that("make_tier_definitions: auto fractions increase with step (2)", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    base_fraction = 0.05,
    step          = 0.05,
    mode          = "fractions"
  )
  testthat::expect_equal(out$fraction[1], 0.05)
  testthat::expect_equal(out$fraction[2], 0.1)
})



testthat::test_that("make_tier_definitions: manual fractions used when provided", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    fractions     = c(0.03, 0.10),
    mode          = "fractions"
  )
  testthat::expect_equal(out$fraction, c(0.03, 0.10))
})

testthat::test_that("make_tier_definitions: non-increasing fractions auto-sorted with warning", {
  testthat::expect_warning(
    out <- make_tier_definitions(
      boundaries_mb = c(50, 25),
      fractions     = c(0.10, 0.03),
      mode          = "fractions"
    ),
    "sorting automatically"  
  )
  testthat::expect_true(all(diff(out$fraction) > 0))
  testthat::expect_equal(out$fraction, c(0.03, 0.10))
})



testthat::test_that("make_tier_definitions: fractions must be in (0, 1)", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      fractions     = c(0, 0.5),
      mode          = "fractions"
    ),
    "must be in"
  )
})

testthat::test_that("make_tier_definitions: negative boundaries raise error", {
  testthat::expect_error(
    make_tier_definitions(boundaries_mb = c(50, -5)),
    "non-negative"
  )
})

testthat::test_that("make_tier_definitions: number mode requires n_cells", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      mode          = "number"
    ),
    "n_cells must be provided"
  )
})

testthat::test_that("make_tier_definitions: number mode n_cells must be >= 1", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      n_cells       = c(0L, 5L),
      mode          = "number"
    ),
    ">= 1"
  )
})


testthat::test_that("make_tier_definitions: number mode warns when fractions provided", {
  testthat::expect_warning(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      n_cells       = c(3L, 8L),
      fractions     = c(0.05, 0.10),
      mode          = "number"
    ),
    "ignored"
  )
})

testthat::test_that("make_tier_definitions: character n_cells raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      n_cells       = c("three", "eight"),
      mode          = "number"
    ),
    "numeric"
  )
})

testthat::test_that("make_tier_definitions: logical n_cells raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      n_cells       = c(TRUE, FALSE),
      mode          = "number"
    ),
    "numeric"
  )
})

testthat::test_that("make_tier_definitions: numeric n_cells accepted and coerced to integer", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    n_cells       = c(3.0, 8.0),  # numeric not integer — should work
    mode          = "number"
  )
  testthat::expect_type(out$n_cells, "integer")
  testthat::expect_equal(out$n_cells, c(3L, 8L))
})

testthat::test_that("make_tier_definitions: n_cells stored as integer in number mode", {
  out <- make_tier_definitions(
    boundaries_mb = c(50, 25),
    n_cells       = c(3L, 8L),
    mode          = "number"
  )
  testthat::expect_type(out$n_cells, "integer")
})


testthat::test_that("make_tier_definitions: n_cells ignored warns in fractions mode", {
  testthat::expect_warning(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      n_cells       = c(3L, 8L),
      mode          = "fractions"
    ),
    "ignored"
  )
})


testthat::test_that("make_tier_definitions: character fractions raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      fractions     = c("0.05", "0.08"),
      mode          = "fractions"
    ),
    "numeric"
  )
})


testthat::test_that("make_tier_definitions: logical fractions raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      fractions     = c(TRUE, FALSE),
      mode          = "fractions"
    ),
    "numeric"
  )
})


testthat::test_that("make_tier_definitions: negative step raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      base_fraction = 0.05,
      step          = -0.03,
      mode          = "fractions"
    ),
    "positive"
  )
})

testthat::test_that("make_tier_definitions: step >= 1 raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      base_fraction = 0.05,
      step          = 1.5,
      mode          = "fractions"
    ),
    "less than 1"
  )
})


testthat::test_that("make_tier_definitions: step that pushes fraction above 1 raises error", {
  # base_fraction = 0.8, step = 0.3 → second fraction = 1.1 → invalid
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      base_fraction = 0.80,
      step          = 0.30,
      mode          = "fractions"
    ),
    "outside \\(0, 1\\)"
  )
})

testthat::test_that("make_tier_definitions: non-numeric step raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      step          = "big",
      mode          = "fractions"
    ),
    "numeric"
  )
})


testthat::test_that("make_tier_definitions: fraction >= 0.25 produces warning", {
  testthat::expect_warning(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      fractions     = c(0.25, 0.50),
      mode          = "fractions"
    ),
    "25%"
  )
})

testthat::test_that("make_tier_definitions: auto-generated fraction >= 0.25 produces warning", {
  # base_fraction = 0.25 → first fraction already at warning threshold
  testthat::expect_warning(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      base_fraction = 0.25,
      step          = 0.03,
      mode          = "fractions"
    ),
    "25%"
  )
})


testthat::test_that("make_tier_definitions: base_fraction <= 0 raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      base_fraction = 0,
      mode          = "fractions"
    ),
    "\\(0, 1\\)"
  )
})

testthat::test_that("make_tier_definitions: base_fraction >= 1 raises error", {
  testthat::expect_error(
    make_tier_definitions(
      boundaries_mb = c(50, 25),
      base_fraction = 1.5,
      mode          = "fractions"
    ),
    "\\(0, 1\\)"
  )
})




testthat::test_that("make_tier_definitions: single boundary produces one tier", {
  out <- make_tier_definitions(
    boundaries_mb = c(25),
    mode          = "fractions"
  )
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$tier, "tier_1")
})



# ============================================================================
# resolve_tier_thresholds
# ============================================================================

testthat::test_that("resolve_tier_thresholds: auto method returns tier table", {
  out <- resolve_tier_thresholds(
    method        = "auto",
    boundaries_mb = c(50, 25),
    base_fraction = 0.05,
    step          = 0.03
  )
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true("tier"     %in% colnames(out))
  testthat::expect_true("fraction" %in% colnames(out))
})

testthat::test_that("resolve_tier_thresholds: manual method uses provided fractions", {
  out <- resolve_tier_thresholds(
    method        = "manual",
    boundaries_mb = c(50, 25),
    fractions     = c(0.03, 0.10),
    mode          = "fractions"
  )
  testthat::expect_equal(out$fraction, c(0.03, 0.10))
})

testthat::test_that("resolve_tier_thresholds: single method with fraction", {
  out <- resolve_tier_thresholds(
    method        = "single",
    boundaries_mb = c(25),
    fractions     = c(0.05),
    mode          = "fractions"
  )
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$fraction, 0.05)
})

testthat::test_that("resolve_tier_thresholds: invalid method raises error", {
  testthat::expect_error(
    resolve_tier_thresholds(method = "invalid"),
    "should be one of"
  )
})


testthat::test_that("resolve_tier_thresholds: more than 2 boundaries raises error", {
  testthat::expect_error(
    resolve_tier_thresholds(
      method        = "auto",
      boundaries_mb = c(100, 50, 25), 
      base_fraction = 0.05,
      step          = 0.03
    ),
    "maximum of 2"
  )
})

testthat::test_that("resolve_tier_thresholds: exactly 2 boundaries accepted", {
  testthat::expect_no_error(
    resolve_tier_thresholds(
      method        = "auto",
      boundaries_mb = c(50, 25),
      base_fraction = 0.05,
      step          = 0.03
    )
  )
})

testthat::test_that("resolve_tier_thresholds: single boundary accepted", {
  testthat::expect_no_error(
    resolve_tier_thresholds(
      method        = "single",
      boundaries_mb = c(25),
      fractions     = c(0.05),
      mode          = "fractions"
    )
  )
})




# ============================================================================
# prepare_cnv_thresholds
# ============================================================================

# Minimal input for prepare_cnv_thresholds


testthat::test_that("prepare_cnv_thresholds: adds tier column", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df       = NULL,
    clustered_events = df,
    by_union         = "embryo",
    cell_sizes       = make_mock_cell_sizes(),
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    mode             = "fractions"
  )
  testthat::expect_true("tier" %in% colnames(out))
})

testthat::test_that("prepare_cnv_thresholds: adds effective_threshold column", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df       = NULL,
    clustered_events = df,
    by_union         = "embryo",
    cell_sizes       = make_mock_cell_sizes(),
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    mode             = "fractions"
  )
  testthat::expect_true("effective_threshold" %in% colnames(out))
})

testthat::test_that("prepare_cnv_thresholds: tier_1 assigned to large CNVs", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df       = NULL,
    clustered_events = df,
    by_union         = "embryo",
    cell_sizes       = make_mock_cell_sizes(),
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    mode             = "fractions"
  )
  # 86.7 Mb → above 25 Mb → tier_1
  large_cnv <- out[out$cnv_length_mb > 25, ]
  testthat::expect_true(all(large_cnv$tier == "tier_1"))
})



testthat::test_that("prepare_cnv_thresholds: below_threshold for small CNVs", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df       = NULL,
    clustered_events = df,
    by_union         = "embryo",
    cell_sizes       = make_mock_cell_sizes(),
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    mode             = "fractions"
  )
  # 9.8 Mb → below 10 Mb → below_threshold
  small_cnv <- out[out$cnv_length_mb < 10, ]
  testthat::expect_true(all(small_cnv$tier == "below_threshold"))
})

testthat::test_that("prepare_cnv_thresholds: effective_threshold respects min_cap", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df        = NULL,
    clustered_events  = df,
    by_union          = "embryo",
    cell_sizes        = make_mock_cell_sizes(),
    boundaries_mb     = c(25, 10),
    base_fraction     = 0.05,
    step              = 0.03,
    mode              = "fractions",
    min_cap_threshold = 5L
  )
  non_na <- out$effective_threshold[!is.na(out$effective_threshold)]
  testthat::expect_true(all(non_na >= 5L))
})

testthat::test_that("prepare_cnv_thresholds: effective_threshold respects max_cap", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df        = NULL,
    clustered_events  = df,
    by_union          = "embryo",
    cell_sizes        = make_mock_cell_sizes(),
    boundaries_mb     = c(25, 10),
    base_fraction     = 0.05,
    step              = 0.03,
    mode              = "fractions",
    max_cap_threshold = 10L
  )
  non_na <- out$effective_threshold[!is.na(out$effective_threshold)]
  testthat::expect_true(all(non_na <= 10L))
})

testthat::test_that("prepare_cnv_thresholds: missing cell_sizes raises error", {
  df <- make_mock_prepare_input()
  testthat::expect_error(
    prepare_cnv_thresholds(
      summary_df       = NULL,
      clustered_events = df,
      by_union         = "embryo",
      cell_sizes       = NULL,
      boundaries_mb    = c(25, 10)
    ),
    "cell numbers"
  )
})

testthat::test_that("prepare_cnv_thresholds: works with single boundary", {
  df  <- make_mock_prepare_input()
  out <- prepare_cnv_thresholds(
    summary_df       = NULL,
    clustered_events = df,
    by_union         = "embryo",
    cell_sizes       = make_mock_cell_sizes(),
    boundaries_mb    = c(25),      
    base_fraction    = 0.05,
    step             = 0.03,
    mode             = "fractions",
    max_tiers        = 1L     
  )
  testthat::expect_true("tier" %in% colnames(out))
  testthat::expect_true(
    all(out$tier %in% c("tier_1", "below_threshold"))
  )
  # CNVs above 25Mb → tier_1
  testthat::expect_true(
    all(out$tier[out$cnv_length_mb >= 25] == "tier_1")
  )
  # CNVs below 25Mb → below_threshold
  testthat::expect_true(
    all(out$tier[out$cnv_length_mb < 25] == "below_threshold")
  )
})

testthat::test_that("prepare_cnv_thresholds: exceeding max_tiers raises error", {
  df <- make_mock_prepare_input()
  testthat::expect_error(
    prepare_cnv_thresholds(
      summary_df       = NULL,
      clustered_events = df,
      by_union         = "embryo",
      cell_sizes       = make_mock_cell_sizes(),
      boundaries_mb    = c(50, 25, 10),
      max_tiers        = 2L
    ),
    "more than"
  )
})

testthat::test_that("prepare_cnv_thresholds: skips merge when columns already present", {
  df  <- make_mock_prepare_input()  
  testthat::expect_message(
    out <- prepare_cnv_thresholds(
      summary_df       = NULL,
      clustered_events = df,
      by_union         = "embryo",
      cell_sizes       = make_mock_cell_sizes(),
      boundaries_mb    = c(25, 10),
      base_fraction    = 0.05,
      step             = 0.03,
      mode             = "fractions"
    ),
    "skipping merge"
  )
})


# ============================================================================
# filter_cnv_loci
# ============================================================================

make_mock_filter_input <- function() {
  data.frame(
    cell_name             = c("E5.10.945", "E5.10.949", "E5.10.934", "E5.10.936"),
    chr                   = c("chr1",  "chr19", "chrX",  "chrX"),
    cnv_state             = c("loss",  "gain",  "gain",  "gain"),
    n_cells               = c(2L,      2L,      1L,      1L),
    arm_class             = c("q_arm", "p_arm", "q_arm", "p_centromere_q"),
    whole_chromosome_gain = c(NA,      14.12,   19.93,   55.55),
    whole_chromosome_loss = c(7.11,    NA,      NA,      NA),
    effective_threshold   = c(2L,      3L,      5L,      2L),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("filter_cnv_loci: removes rows below effective_threshold", {
  df  <- make_mock_filter_input()
  # n_cells = 1 for chrX gain → below effective_threshold = 5
  out <- filter_cnv_loci(
    clustered_events            = df,
    total_chromosome_permission = 70
  )
  testthat::expect_true(all(out$n_cells >= out$effective_threshold,
                            na.rm = TRUE))
})

testthat::test_that("filter_cnv_loci: removes p_centromere_q below permission threshold", {
  df  <- make_mock_filter_input()
  out <- filter_cnv_loci(
    clustered_events            = df,
    total_chromosome_permission = 70
  )
  testthat::expect_false("E5.10.936" %in% out$cell_name)
})

testthat::test_that("filter_cnv_loci: keeps p_centromere_q above permission threshold", {
  df       <- make_mock_filter_input()
  df$whole_chromosome_gain[df$arm_class == "p_centromere_q"] <- 80
  df$n_cells[df$arm_class == "p_centromere_q"] <- 5L
  out <- filter_cnv_loci(
    clustered_events            = df,
    total_chromosome_permission = 70
  )
  testthat::expect_true("E5.10.936" %in% out$cell_name)
})

testthat::test_that("filter_cnv_loci: removes effective_threshold from output", {
  df  <- make_mock_filter_input()
  out <- filter_cnv_loci(
    clustered_events            = df,
    total_chromosome_permission = 70
  )
  testthat::expect_false("effective_threshold" %in% colnames(out))
})

testthat::test_that("filter_cnv_loci: missing required columns raises error", {
  df      <- make_mock_filter_input()
  df$n_cells <- NULL
  testthat::expect_error(
    filter_cnv_loci(clustered_events = df),
    "Missing required columns"
  )
})

testthat::test_that("filter_cnv_loci: NA effective_threshold rows removed", {
  df <- make_mock_filter_input()
  df$effective_threshold[1] <- NA
  out <- filter_cnv_loci(
    clustered_events            = df,
    total_chromosome_permission = 70
  )
  testthat::expect_false("E5.10.945" %in% out$cell_name)
})

# ============================================================================
# classify_cnv_loci
# ============================================================================

make_mock_classify_input <- function() {
  data.frame(
    cell_name = c("E5.10.945", "E5.10.949", "E5.10.934"),
    chr       = c("chr1",      "chr19",     "chrX"),
    cnv_state = c("loss",      "gain",      "gain"),
    tier      = c("tier_1",    "tier_2",    "tier_1"),
    n_cells   = c(5L,          3L,          2L),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("classify_cnv_loci: adds confidence column", {
  df  <- make_mock_classify_input()
  out <- classify_cnv_loci(df)
  testthat::expect_true("confidence" %in% colnames(out))
})

testthat::test_that("classify_cnv_loci: tier_1 classified as High", {
  df  <- make_mock_classify_input()
  out <- classify_cnv_loci(df)
  tier1_rows <- out[out$tier == 1, ]
  testthat::expect_true(all(tier1_rows$confidence == "High"))
})

testthat::test_that("classify_cnv_loci: tier_2 classified as Low", {
  df  <- make_mock_classify_input()
  out <- classify_cnv_loci(df)
  tier2_rows <- out[out$tier == 2, ]
  testthat::expect_true(all(tier2_rows$confidence == "Low"))
})

testthat::test_that("classify_cnv_loci: confidence only contains High or Low", {
  df  <- make_mock_classify_input()
  out <- classify_cnv_loci(df)
  testthat::expect_true(all(out$confidence %in% c("High", "Low", NA)))
})

testthat::test_that("classify_cnv_loci: row count unchanged", {
  df  <- make_mock_classify_input()
  out <- classify_cnv_loci(df)
  testthat::expect_equal(nrow(out), nrow(df))
})

testthat::test_that("classify_cnv_loci: tier converted to numeric", {
  df  <- make_mock_classify_input()
  out <- classify_cnv_loci(df)
  testthat::expect_type(out$tier, "double")
})


testthat::test_that("filter_cnv_loci: all rows fail threshold returns empty data frame", {
  df <- make_mock_filter_input()
  # Set effective_threshold above n_cells for all rows
  df$effective_threshold <- 100L
  df$arm_class           <- "q_arm"  # no centromere complication
  
  # Expect error since n_before > 0 but 100 * (n_before - n_after) / n_before
  # division should still work — just returns empty df
  out <- filter_cnv_loci(
    clustered_events            = df,
    total_chromosome_permission = 70
  )
  testthat::expect_equal(nrow(out), 0L)
  testthat::expect_true("cnv_state" %in% colnames(out))  # structure preserved
})

testthat::test_that("score_cnv_clusters: high fraction may return empty result", {

  out <- suppressWarnings(score_cnv_clusters(
    summary_df       = make_mock_locus_summary(),
    clustered_events = make_mock_clustered_events(),
    cell_sizes       = make_mock_cell_sizes(),
    by_union         = "embryo",
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.99,  
    step             = 0.003,
    threshold_method = "auto",
    threshold_mode   = "fractions",
    min_cap_threshold = 2L,
    max_cap_threshold = 25L,
    total_chromosome_permission = 70
  ))
  # Should return empty data frame not an error
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(nrow(out), 0L)
})

# ============================================================================
# summarise_cnv_loci
# ============================================================================

make_mock_summarise_input <- function() {
  data.frame(
    embryo       = c("E5.10", "E5.10", "E5.10", "E5.10"),
    cnv_equiv_id = c("E5.10|chr1|loss|1",  "E5.10|chr1|loss|1",
                     "E5.10|chr19|gain|1", "E5.10|chr19|gain|1"),
    chr          = c("chr1",  "chr1",  "chr19", "chr19"),
    cnv_state    = c("loss",  "loss",  "gain",  "gain"),
    start        = c(231223763L, 231332753L, 301444L, 301444L),
    end          = c(248919946L, 248919946L, 8577577L, 10115372L),
    cell_name    = c("E5.10.945", "E5.10.949", "E5.10.945", "E5.10.949"),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("summarise_cnv_loci: one row per equiv group", {
  df  <- make_mock_summarise_input()
  out <- summarise_cnv_loci(
    df         = df,
    by         = "embryo",
    sample_col = "embryo",
    cell_col   = "cell_name"
  )
  testthat::expect_equal(nrow(out), 2L)
})

testthat::test_that("summarise_cnv_loci: locus_start is min start", {
  df  <- make_mock_summarise_input()
  out <- summarise_cnv_loci(
    df         = df,
    by         = "embryo",
    sample_col = "embryo",
    cell_col   = "cell_name"
  )
  chr1_row <- out[out$chr == "chr1", ]
  testthat::expect_equal(chr1_row$locus_start, 231223763L)
})

testthat::test_that("summarise_cnv_loci: locus_end is max end", {
  df  <- make_mock_summarise_input()
  out <- summarise_cnv_loci(
    df         = df,
    by         = "embryo",
    sample_col = "embryo",
    cell_col   = "cell_name"
  )
  chr19_row <- out[out$chr == "chr19", ]
  testthat::expect_equal(chr19_row$locus_end, 10115372L)
})

testthat::test_that("summarise_cnv_loci: n_cells counts distinct cells", {
  df  <- make_mock_summarise_input()
  out <- summarise_cnv_loci(
    df         = df,
    by         = "embryo",
    sample_col = "embryo",
    cell_col   = "cell_name"
  )
  testthat::expect_equal(out$n_cells[out$chr == "chr1"], 2L)
})

testthat::test_that("summarise_cnv_loci: missing required columns raises error", {
  df      <- make_mock_summarise_input()
  df$chr  <- NULL
  testthat::expect_error(
    summarise_cnv_loci(df, by = "embryo",
                       sample_col = "embryo", cell_col = "cell_name"),
    "Missing required columns"
  )
})

testthat::test_that("summarise_cnv_loci: missing cell_col raises error", {
  df <- make_mock_summarise_input()
  testthat::expect_error(
    summarise_cnv_loci(df, by = "embryo",
                       sample_col = "embryo",
                       cell_col   = "nonexistent_col"),
    "Missing cell column"
  )
})

testthat::test_that("summarise_cnv_loci: locus_width_mb computed correctly", {
  df  <- make_mock_summarise_input()
  out <- summarise_cnv_loci(
    df         = df,
    by         = "embryo",
    sample_col = "embryo",
    cell_col   = "cell_name"
  )
  chr1_row <- out[out$chr == "chr1", ]
  expected_mb <- (248919946L - 231223763L + 1L) / 1e6
  testthat::expect_equal(chr1_row$locus_width_mb, expected_mb,
                         tolerance = 1e-3)
})

# ============================================================================
# score_cnv_clusters — integration test
# ============================================================================

testthat::test_that("score_cnv_clusters: returns data frame", {
  out <- score_cnv_clusters(
    summary_df       = make_mock_locus_summary(),
    clustered_events = make_mock_clustered_events(),
    cell_sizes       = make_mock_cell_sizes(),
    by_union         = "embryo",
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    threshold_method = "auto",
    threshold_mode   = "fractions",
    min_cap_threshold = 2L,
    max_cap_threshold = 25L,
    total_chromosome_permission = 70
  )
  testthat::expect_s3_class(out, "data.frame")
})

testthat::test_that("score_cnv_clusters: output has confidence column", {
  out <- score_cnv_clusters(
    summary_df       = make_mock_locus_summary(),
    clustered_events = make_mock_clustered_events(),
    cell_sizes       = make_mock_cell_sizes(),
    by_union         = "embryo",
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    threshold_method = "auto",
    threshold_mode   = "fractions",
    min_cap_threshold = 2L,
    max_cap_threshold = 25L,
    total_chromosome_permission = 70
  )
  testthat::expect_true("confidence" %in% colnames(out))
})

testthat::test_that("score_cnv_clusters: confidence values are High or Low only", {
  out <- score_cnv_clusters(
    summary_df       = make_mock_locus_summary(),
    clustered_events = make_mock_clustered_events(),
    cell_sizes       = make_mock_cell_sizes(),
    by_union         = "embryo",
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    threshold_method = "auto",
    threshold_mode   = "fractions",
    min_cap_threshold = 2L,
    max_cap_threshold = 25L,
    total_chromosome_permission = 70
  )
  testthat::expect_true(all(out$confidence %in% c("High", "Low", NA)))
})

testthat::test_that("score_cnv_clusters: output has fewer rows than input after filtering", {
  input <- make_mock_clustered_events()
  out <- score_cnv_clusters(
    summary_df       = make_mock_locus_summary(),
    clustered_events = input,
    cell_sizes       = make_mock_cell_sizes(),
    by_union         = "embryo",
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    threshold_method = "auto",
    threshold_mode   = "fractions",
    min_cap_threshold = 2L,
    max_cap_threshold = 25L,
    total_chromosome_permission = 70
  )
  testthat::expect_lte(nrow(out), nrow(input))
})

testthat::test_that("score_cnv_clusters: no effective_threshold in output", {
  out <- score_cnv_clusters(
    summary_df       = make_mock_locus_summary(),
    clustered_events = make_mock_clustered_events(),
    cell_sizes       = make_mock_cell_sizes(),
    by_union         = "embryo",
    boundaries_mb    = c(25, 10),
    base_fraction    = 0.05,
    step             = 0.03,
    threshold_method = "auto",
    threshold_mode   = "fractions",
    min_cap_threshold = 2L,
    max_cap_threshold = 25L,
    total_chromosome_permission = 70
  )
  testthat::expect_false("effective_threshold" %in% colnames(out))
})