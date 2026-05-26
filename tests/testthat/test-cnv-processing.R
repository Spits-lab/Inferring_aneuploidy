
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/tests/testthat/helper-fixtures.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/cnv_processing.R")

# ============================================================================
# compute_overlap
# ============================================================================

testthat::test_that("compute_overlap: perfect overlap returns 1", {
  score <- compute_overlap(
    q_start = 100L, q_end = 200L,
    s_start = 100L, s_end = 200L,
    method  = "reciprocal"
  )
  testthat::expect_equal(score, 1)
})

testthat::test_that("compute_overlap: no overlap returns 0", {
  score <- compute_overlap(
    q_start = 100L, q_end = 200L,
    s_start = 300L, s_end = 400L,
    method  = "reciprocal"
  )
  testthat::expect_equal(score, 0)
})

testthat::test_that("compute_overlap: reciprocal uses larger segment as denominator", {
  # q = 26bp, s = 101bp, intersection = 26bp
  # reciprocal = 26 / max(26, 101) = 26/101
  score <- compute_overlap(
    q_start = 150L, q_end = 175L,
    s_start = 100L, s_end = 200L,
    method  = "reciprocal"
  )
  testthat::expect_equal(score, 26 / 101)
})

testthat::test_that("compute_overlap: jaccard is symmetric", {
  score_ab <- compute_overlap(100L, 200L, 150L, 250L, method = "jaccard")
  score_ba <- compute_overlap(150L, 250L, 100L, 200L, method = "jaccard")
  testthat::expect_equal(score_ab, score_ba)
})

testthat::test_that("compute_overlap: jaccard perfect overlap returns 1", {
  score <- compute_overlap(100L, 200L, 100L, 200L, method = "jaccard")
  testthat::expect_equal(score, 1)
})

testthat::test_that("compute_overlap: symmetric_reciprocal is symmetric", {
  score_ab <- compute_overlap(100L, 200L, 150L, 250L,
                              method = "symmetric_reciprocal")
  score_ba <- compute_overlap(150L, 250L, 100L, 200L,
                              method = "symmetric_reciprocal")
  testthat::expect_equal(score_ab, score_ba)
})

testthat::test_that("compute_overlap: all methods return values in [0, 1]", {
  methods <- c("reciprocal", "jaccard", "symmetric_reciprocal")
  for (m in methods) {
    score <- compute_overlap(100L, 200L, 150L, 250L, method = m)
    testthat::expect_gte(score, 0, label = paste(m, ">= 0"))
    testthat::expect_lte(score, 1, label = paste(m, "<= 1"))
  }
})

testthat::test_that("compute_overlap: unknown method raises error", {
  testthat::expect_error(
    compute_overlap(100L, 200L, 100L, 200L, method = "invalid"),
    "Unknown overlap method"
  )
})

testthat::test_that("compute_overlap: vectorised input returns correct length", {
  scores <- compute_overlap(
    q_start = c(100L, 200L, 300L),
    q_end   = c(150L, 250L, 350L),
    s_start = c(125L, 225L, 400L),
    s_end   = c(175L, 275L, 450L),
    method  = "reciprocal"
  )
  testthat::expect_length(scores, 3L)
})

testthat::test_that("compute_overlap: real coordinates from embryo data", {
  # chr1 loss E5.10.933: 117863485-150476566
  # same segment from ref B — perfect overlap expected
  score <- compute_overlap(
    q_start = 117863485L, q_end = 150476566L,
    s_start = 117863485L, s_end = 150476566L,
    method  = "reciprocal"
  )
  testthat::expect_equal(score, 1)
})

# ============================================================================
# collapse_genes_to_cnv_segments
# ============================================================================

testthat::test_that("collapse_genes_to_cnv_segments: returns fewer rows than input", {
  df  <- make_mock_gene_level()
  # Keep only gain state to ensure collapsing happens
  df  <- df[df$state != "neutral", ]
  out <- collapse_genes_to_cnv_segments(df)
  testthat::expect_lt(nrow(out), nrow(df))
})

testthat::test_that("collapse_genes_to_cnv_segments: neutral rows are removed", {
  df  <- make_mock_gene_level()
  out <- collapse_genes_to_cnv_segments(df)
  testthat::expect_false("neutral" %in% out$state)
})

testthat::test_that("collapse_genes_to_cnv_segments: output has required columns", {
  df  <- make_mock_gene_level()
  out <- collapse_genes_to_cnv_segments(df)
  required <- c("reference", "cell_name", "chr", "start", "stop", "state")
  testthat::expect_true(all(required %in% colnames(out)))
})

testthat::test_that("collapse_genes_to_cnv_segments: start <= stop in all rows", {
  df  <- make_mock_gene_level()
  out <- collapse_genes_to_cnv_segments(df)
  testthat::expect_true(all(out$start <= out$stop))
})

testthat::test_that("collapse_genes_to_cnv_segments: missing required columns raises error", {
  df      <- make_mock_gene_level()
  df$chr  <- NULL
  testthat::expect_error(
    collapse_genes_to_cnv_segments(df),
    "Missing required columns"
  )
})

testthat::test_that("collapse_genes_to_cnv_segments: consecutive same-state genes collapsed", {
  # Three consecutive gain genes on chr1 for same cell and reference
  # Should collapse to one segment
  df <- data.frame(
    gene      = c("G1",   "G2",   "G3"),
    cell_name = c("cell_1", "cell_1", "cell_1"),
    state_raw = c(1.5,    1.5,    1.5),
    chr       = c("chr1", "chr1", "chr1"),
    start     = c(1000L,  2000L,  3000L),
    stop      = c(1500L,  2500L,  3500L),
    state     = c("gain", "gain", "gain"),
    reference = c("A",    "A",    "A"),
    stringsAsFactors = FALSE
  )
  out <- collapse_genes_to_cnv_segments(df)
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$start, 1000L)
  testthat::expect_equal(out$stop,  3500L)
})

testthat::test_that("collapse_genes_to_cnv_segments: different states not collapsed", {
  df <- data.frame(
    gene      = c("G1",     "G2"),
    cell_name = c("cell_1", "cell_1"),
    state_raw = c(1.5,      0.5),
    chr       = c("chr15",  "chr15"),
    start     = c(1000000L, 2000000L),
    stop      = c(1500000L, 2500000L),
    state     = c("gain",   "loss"),
    reference = c("A",      "A"),
    stringsAsFactors = FALSE
  )
  # No reduction expected — different states cannot collapse
  # Warning is correct behaviour here
  testthat::expect_warning(
    out <- collapse_genes_to_cnv_segments(df),
    "did not decrease"
  )
  testthat::expect_equal(nrow(out), 2L)
})

testthat::test_that("collapse_genes_to_cnv_segments: different cells not collapsed", {
  df <- data.frame(
    gene      = c("G1",     "G1"),
    cell_name = c("cell_1", "cell_2"),
    state_raw = c(1.5,      1.5),
    chr       = c("chr15",  "chr15"),
    start     = c(1000000L, 1000000L),
    stop      = c(1500000L, 1500000L),
    state     = c("gain",   "gain"),
    reference = c("A",      "A"),
    stringsAsFactors = FALSE
  )
  testthat::expect_warning(
    out <- collapse_genes_to_cnv_segments(df),
    "did not decrease"
  )
  testthat::expect_equal(nrow(out), 2L)
})

# ============================================================================
# merge_nearby_regions
# ============================================================================

testthat::test_that("merge_nearby_regions: nearby segments are merged", {
  df <- make_mock_segments()
  # Add two nearby segments for same cell/ref/chr/state
  df_extra <- data.frame(
    reference = c("A",         "A"),
    cell_name = c("E5.10.933", "E5.10.933"),
    chr       = c("chr5",      "chr5"),
    state     = c("gain",      "gain"),
    start     = c(92151L,      1345100L),  # gap = ~1.25Mb < max_gap
    stop      = c(1345099L,    6757048L),
    stringsAsFactors = FALSE
  )
  df  <- rbind(df, df_extra)
  out <- merge_nearby_regions(df, max_gap = 2000000L)
  testthat::expect_lt(nrow(out), nrow(df))
})

testthat::test_that("merge_nearby_regions: distant segments not merged", {
  df  <- make_mock_segments()
  testthat::expect_warning(
    out <- merge_nearby_regions(df, max_gap = 100L),
    "No reduction"
  )
  testthat::expect_equal(nrow(out), nrow(df))
})

testthat::test_that("merge_nearby_regions: output column end exists", {
  df <- make_mock_segments()
  # Add two nearby segments that will merge
  df_extra <- data.frame(
    reference  = c("A",         "A"),
    cell_name  = c("E5.10.933", "E5.10.933"),
    chr        = c("chr15",     "chr15"),
    state      = c("gain",      "gain"),
    start      = c(92151L,      500000L),
    stop       = c(499999L,     6757048L),
    stringsAsFactors = FALSE
  )
  df  <- rbind(df, df_extra)
  out <- merge_nearby_regions(df, max_gap = 100000L)
  testthat::expect_true("end" %in% colnames(out))
})

testthat::test_that("merge_nearby_regions: start <= end in all rows", {
  df <- make_mock_segments()
  df_extra <- data.frame(
    reference  = c("A",         "A"),
    cell_name  = c("E5.10.933", "E5.10.933"),
    chr        = c("chr15",     "chr15"),
    state      = c("gain",      "gain"),
    start      = c(92151L,      500000L),
    stop       = c(499999L,     6757048L),
    stringsAsFactors = FALSE
  )
  df  <- rbind(df, df_extra)
  out <- merge_nearby_regions(df, max_gap = 100000L)
  testthat::expect_true(all(out$start <= out$end))
})

testthat::test_that("merge_nearby_regions: missing required columns raises error", {
  df      <- make_mock_segments()
  df$chr  <- NULL
  testthat::expect_error(
    merge_nearby_regions(df),
    "Missing required columns"
  )
})

testthat::test_that("merge_nearby_regions: start >= stop raises error", {
  df         <- make_mock_segments()
  df$stop[1] <- df$start[1] - 1L
  testthat::expect_error(
    merge_nearby_regions(df),
    "start >= stop"
  )
})

testthat::test_that("merge_nearby_regions: merged segment spans full range", {
  df <- data.frame(
    reference = c("A",    "A"),
    cell_name = c("cell_1", "cell_1"),
    chr       = c("chr1", "chr1"),
    state     = c("gain", "gain"),
    start     = c(1000L,  1500L),
    stop      = c(1200L,  2000L),
    stringsAsFactors = FALSE
  )
  out <- merge_nearby_regions(df, max_gap = 1000L)
  testthat::expect_equal(out$start, 1000L)
  testthat::expect_equal(out$end,   2000L)
})

testthat::test_that("merge_nearby_regions: gain and loss segments not merged even when adjacent", {
  df <- data.frame(
    reference = c("A",      "A"),
    cell_name = c("cell_1", "cell_1"),
    chr       = c("chr15",  "chr15"),
    state     = c("gain",   "loss"),
    start     = c(1000000L, 1500000L),
    stop      = c(1200000L, 2000000L),
    stringsAsFactors = FALSE
  )
  # No reduction expected — different states never merge
  # Warning is correct behaviour
  testthat::expect_warning(
    out <- merge_nearby_regions(df, max_gap = 1000000L),
    "No reduction"
  )
  testthat::expect_equal(nrow(out), 2L)
  # Both states preserved
  testthat::expect_true("gain" %in% out$cnv_state)
  testthat::expect_true("loss" %in% out$cnv_state)
})



testthat::test_that("merge_nearby_regions: same state adjacent segments merged without warning", {
  df <- data.frame(
    reference = c("A",      "A"),
    cell_name = c("cell_1", "cell_1"),
    chr       = c("chr15",  "chr15"),
    state     = c("gain",   "gain"),
    start     = c(1000000L, 1500000L),
    stop      = c(1200000L, 2000000L),
    stringsAsFactors = FALSE
  )
  # Reduction expected — same state merges → no warning
  testthat::expect_no_warning(
    out <- merge_nearby_regions(df, max_gap = 1000000L)
  )
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$start, 1000000L)
  testthat::expect_equal(out$end,   2000000L)
})


testthat::test_that("merge_nearby_regions: same state different chromosomes not merged", {
  df <- data.frame(
    reference = c("A",      "A"),
    cell_name = c("cell_1", "cell_1"),
    chr       = c("chr15",  "chr22"),  # different chromosomes
    state     = c("gain",   "gain"),
    start     = c(1000000L, 1000000L),
    stop      = c(1200000L, 1200000L),
    stringsAsFactors = FALSE
  )
  testthat::expect_warning(
    out <- merge_nearby_regions(df, max_gap = 1000000L),
    "No reduction"
  )
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_true(all(c("chr15", "chr22") %in% out$chr))
})

testthat::test_that("merge_nearby_regions: same state different references not merged", {
  df <- data.frame(
    reference = c("A",      "B"),
    cell_name = c("cell_1", "cell_1"),
    chr       = c("chr15",  "chr15"),
    state     = c("gain",   "gain"),
    start     = c(1000000L, 1500000L),
    stop      = c(1200000L, 2000000L),
    stringsAsFactors = FALSE
  )
  testthat::expect_warning(
    out <- merge_nearby_regions(df, max_gap = 1000000L),
    "No reduction"
  )
  testthat::expect_equal(nrow(out), 2L)
})


testthat::test_that("merge_nearby_regions: same state different cells not merged", {
  df <- data.frame(
    reference = c("A",      "A"),
    cell_name = c("cell_1", "cell_2"),
    chr       = c("chr15",  "chr15"),
    state     = c("gain",   "gain"),
    start     = c(1000000L, 1500000L),
    stop      = c(1200000L, 2000000L),
    stringsAsFactors = FALSE
  )
  testthat::expect_warning(
    out <- merge_nearby_regions(df, max_gap = 1000000L),
    "No reduction"
  )
  testthat::expect_equal(nrow(out), 2L)
})


testthat::test_that("merge_nearby_regions: three adjacent segments merge into one", {
  df <- data.frame(
    reference = c("A",      "A",      "A"),
    cell_name = c("cell_1", "cell_1", "cell_1"),
    chr       = c("chr15",  "chr15",  "chr15"),
    state     = c("gain",   "gain",   "gain"),
    start     = c(1000000L, 1500000L, 2000000L),
    stop      = c(1200000L, 1800000L, 2500000L),
    stringsAsFactors = FALSE
  )
  out <- merge_nearby_regions(df, max_gap = 500000L)
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$start, 1000000L)
  testthat::expect_equal(out$end,   2500000L)
})

testthat::test_that("merge_nearby_regions: three adjacent segments merge into one", {
  df <- data.frame(
    reference = c("A",      "A",      "A"),
    cell_name = c("cell_1", "cell_1", "cell_1"),
    chr       = c("chr15",  "chr15",  "chr15"),
    state     = c("gain",   "gain",   "gain"),
    start     = c(1000000L, 1500000L, 2000000L),
    stop      = c(1200000L, 1800000L, 2500000L),
    stringsAsFactors = FALSE
  )
  out <- merge_nearby_regions(df, max_gap = 500000L)
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$start, 1000000L)
  testthat::expect_equal(out$end,   2500000L)
})


testthat::test_that("merge_nearby_regions: n_segments reflects number of merged segments", {
  df <- data.frame(
    reference = c("A",      "A",      "A"),
    cell_name = c("cell_1", "cell_1", "cell_1"),
    chr       = c("chr15",  "chr15",  "chr15"),
    state     = c("gain",   "gain",   "gain"),
    start     = c(1000000L, 1500000L, 2000000L),
    stop      = c(1200000L, 1800000L, 2500000L),
    stringsAsFactors = FALSE
  )
  out <- merge_nearby_regions(df, max_gap = 500000L)
  testthat::expect_equal(out$n_segments, 3L)
})


testthat::test_that("compute_overlap: single base overlap returns small positive value", {
  score <- compute_overlap(
    q_start = 1000000L, q_end = 2000000L,
    s_start = 2000000L, s_end = 3000000L,
    method  = "reciprocal"
  )
  testthat::expect_gt(score, 0)
  testthat::expect_lt(score, 0.01)
})

testthat::test_that("compute_overlap: reciprocal penalises containment correctly", {
  # Small segment fully inside large — reciprocal uses max length
  score_contained <- compute_overlap(
    q_start = 1000000L, q_end = 1100000L,   # 100kb
    s_start = 1000000L, s_end = 10000000L,  # 9Mb
    method  = "reciprocal"
  )
  # Should be low — small segment contained in large
  testthat::expect_lt(score_contained, 0.05)
})


testthat::test_that("compute_overlap: jaccard and reciprocal differ for asymmetric segments", {
  # q = 100-160Mb (60Mb), s = 130-180Mb (50Mb)
  # intersection = 30Mb
  # reciprocal = 30 / max(60, 50) = 30/60 = 0.500
  # jaccard    = 30 / (60 + 50 - 30) = 30/80 = 0.375
  # These differ — good test case
  
  score_rec <- compute_overlap(
    100000000L, 160000000L,
    130000000L, 180000000L,
    method = "reciprocal"
  )
  score_jac <- compute_overlap(
    100000000L, 160000000L,
    130000000L, 180000000L,
    method = "jaccard"
  )
  
  # Verify expected values explicitly
  testthat::expect_equal(score_rec, 30000001 / 60000001,
                         tolerance = 1e-6)
  testthat::expect_equal(score_jac, 30000001 / 80000001,
                         tolerance = 1e-6)
  
  # Confirm they differ
  testthat::expect_false(isTRUE(all.equal(score_rec, score_jac)))
})


testthat::test_that("filter_cnv_events: references column preserved after filtering", {
  df  <- make_mock_supported_events()
  out <- filter_cnv_events(df, min_references = 2L)
  testthat::expect_true("references" %in% colnames(out))
  testthat::expect_true(all(grepl("A|B|C", out$references)))
})



testthat::test_that("filter_cnv_events: all columns preserved after filtering", {
  df  <- make_mock_supported_events()
  out <- filter_cnv_events(df, min_references = 2L)
  testthat::expect_equal(colnames(out), colnames(df))
})



testthat::test_that("summarize_cnv_support: references column lists all references", {
  df <- data.frame(
    cnv_equiv_id = c("cell1|chr15|loss|1", "cell1|chr15|loss|1",
                     "cell1|chr15|loss|1"),
    reference    = c("A", "B", "C"),
    cell_name    = c("cell_1", "cell_1", "cell_1"),
    chr          = c("chr15",  "chr15",  "chr15"),
    cnv_state    = c("loss",   "loss",   "loss"),
    start        = c(100L,     100L,     100L),
    end          = c(500L,     500L,     500L),
    stringsAsFactors = FALSE
  )
  out <- summarize_cnv_support(df)
  testthat::expect_true(grepl("A", out$references))
  testthat::expect_true(grepl("B", out$references))
  testthat::expect_true(grepl("C", out$references))
})



testthat::test_that("summarize_cnv_support: cnv_length_mb matches end - start + 1", {
  df <- data.frame(
    cnv_equiv_id = c("cell1|chr15|loss|1", "cell1|chr15|loss|1"),
    reference    = c("A", "B"),
    cell_name    = c("cell_1", "cell_1"),
    chr          = c("chr15",  "chr15"),
    cnv_state    = c("loss",   "loss"),
    start        = c(10000000L, 10000000L),
    end          = c(50000000L, 50000000L),
    stringsAsFactors = FALSE
  )
  out <- summarize_cnv_support(df)
  
  # cnv_length uses inclusive coordinates: end - start + 1
  expected_mb <- (50000000L - 10000000L + 1L) / 1e6
  testthat::expect_equal(out$cnv_length_mb, expected_mb)
})



testthat::test_that("collapse_genes_to_cnv_segments: gain and loss on same chr kept separate", {
  df <- data.frame(
    gene      = c("G1",     "G2"),
    cell_name = c("cell_1", "cell_1"),
    state_raw = c(1.5,      0.5),
    chr       = c("chr15",  "chr15"),
    start     = c(1000000L, 5000000L),
    stop      = c(2000000L, 6000000L),
    state     = c("gain",   "loss"),
    reference = c("A",      "A"),
    stringsAsFactors = FALSE
  )
  testthat::expect_warning(
    out <- collapse_genes_to_cnv_segments(df),
    "did not decrease"
  )
  testthat::expect_equal(nrow(out), 2L)
  testthat::expect_true("gain" %in% out$state)
  testthat::expect_true("loss" %in% out$state)
})



testthat::test_that("collapse_genes_to_cnv_segments: output start always <= stop", {
  df  <- make_mock_gene_level()
  df  <- df[df$state != "neutral", ]
  out <- collapse_genes_to_cnv_segments(df)
  testthat::expect_true(all(out$start <= out$stop))
})



testthat::test_that("collapse_genes_to_cnv_segments: collapsed segment spans full gene range", {
  df <- data.frame(
    gene      = c("G1",     "G2",     "G3"),
    cell_name = c("cell_1", "cell_1", "cell_1"),
    state_raw = c(1.5,      1.5,      1.5),
    chr       = c("chr15",  "chr15",  "chr15"),
    start     = c(1000000L, 2000000L, 3000000L),
    stop      = c(1500000L, 2500000L, 3500000L),
    state     = c("gain",   "gain",   "gain"),
    reference = c("A",      "A",      "A"),
    stringsAsFactors = FALSE
  )
  out <- collapse_genes_to_cnv_segments(df)
  testthat::expect_equal(out$start, 1000000L)
  testthat::expect_equal(out$stop,  3500000L)
})


# ============================================================================
# filter_cnv_events
# ============================================================================

testthat::test_that("filter_cnv_events: keeps events with sufficient references", {
  df  <- make_mock_supported_events()
  out <- filter_cnv_events(df, min_references = 2L)
  testthat::expect_true(all(out$n_references >= 2L))
})

testthat::test_that("filter_cnv_events: removes events below threshold", {
  df <- make_mock_supported_events()
  df$n_references[1] <- 1L  # force one to fail
  out <- filter_cnv_events(df, min_references = 2L)
  testthat::expect_equal(nrow(out), nrow(df) - 1L)
})

testthat::test_that("filter_cnv_events: min_references = 1 keeps all rows", {
  df  <- make_mock_supported_events()
  out <- filter_cnv_events(df, min_references = 1L)
  testthat::expect_equal(nrow(out), nrow(df))
})

testthat::test_that("filter_cnv_events: missing n_references raises error", {
  df               <- make_mock_supported_events()
  df$n_references  <- NULL
  testthat::expect_error(
    filter_cnv_events(df),
    "Missing required column: n_references"
  )
})

testthat::test_that("filter_cnv_events: all below threshold returns 0 rows", {
  df              <- make_mock_supported_events()
  df$n_references <- 1L
  out             <- filter_cnv_events(df, min_references = 3L)
  testthat::expect_equal(nrow(out), 0L)
})

# ============================================================================
# summarize_cnv_support
# ============================================================================

testthat::test_that("summarize_cnv_support: returns one row per equiv group", {
  # Build input with multiple rows per equiv group (from different references)
  df <- data.frame(
    cnv_equiv_id = c("cell1|chr1|loss|1", "cell1|chr1|loss|1",
                     "cell1|chrX|gain|1"),
    reference    = c("A", "B", "A"),
    cell_name    = c("cell_1", "cell_1", "cell_1"),
    chr          = c("chr1",   "chr1",   "chrX"),
    cnv_state    = c("loss",   "loss",   "gain"),
    start        = c(100L,     105L,     200L),
    end          = c(500L,     490L,     800L),
    stringsAsFactors = FALSE
  )
  out <- summarize_cnv_support(df)
  testthat::expect_equal(nrow(out), 2L)
})

testthat::test_that("summarize_cnv_support: n_references counts distinct references", {
  df <- data.frame(
    cnv_equiv_id = c("cell1|chr1|loss|1", "cell1|chr1|loss|1"),
    reference    = c("A", "B"),
    cell_name    = c("cell_1", "cell_1"),
    chr          = c("chr1",   "chr1"),
    cnv_state    = c("loss",   "loss"),
    start        = c(100L,     100L),
    end          = c(500L,     500L),
    stringsAsFactors = FALSE
  )
  out <- summarize_cnv_support(df)
  testthat::expect_equal(out$n_references, 2L)
})

testthat::test_that("summarize_cnv_support: start is min, end is max across references", {
  df <- data.frame(
    cnv_equiv_id = c("cell1|chr1|loss|1", "cell1|chr1|loss|1"),
    reference    = c("A",  "B"),
    cell_name    = c("cell_1", "cell_1"),
    chr          = c("chr1",   "chr1"),
    cnv_state    = c("loss",   "loss"),
    start        = c(100L,     90L),   # B starts earlier
    end          = c(500L,     510L),  # B ends later
    stringsAsFactors = FALSE
  )
  out <- summarize_cnv_support(df)
  testthat::expect_equal(out$start, 90L)
  testthat::expect_equal(out$end,   510L)
})

testthat::test_that("summarize_cnv_support: missing required columns raises error", {
  df <- data.frame(
    cnv_equiv_id = "cell1|chr1|loss|1",
    reference    = "A",
    cell_name    = "cell_1",
    chr          = "chr1"
    # missing cnv_state, start, end
  )
  testthat::expect_error(summarize_cnv_support(df))
})

# ============================================================================
# assign_cnv_equivalence — scenario-based tests
# ============================================================================

testthat::test_that("equivalence scenario 1: perfect A+B+C overlap → one clique", {
  df  <- make_scenario_perfect_overlap()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  # All three references → one equiv group
  testthat::expect_equal(
    dplyr::n_distinct(out$results_id$cnv_equiv_id), 1L
  )
  # Summarised → n_references = 3
  summary <- summarize_cnv_support(out$results_id)
  testthat::expect_equal(summary$n_references, 3L)
})

testthat::test_that("equivalence scenario 2: A+B overlap, C isolated → C in removed_log", {
  df  <- make_scenario_ab_only()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  # C is isolated — should appear in removed_log
  testthat::expect_true(nrow(out$removed_log) > 0L)
  # One equiv group for A+B
  testthat::expect_equal(
    dplyr::n_distinct(out$results_id$cnv_equiv_id), 1L
  )
})

testthat::test_that("equivalence scenario 3: A+B, A+C overlap but not B+C → two cliques, A duplicated", {
  df  <- make_scenario_ac_ab_no_bc()
  
  # Verify overlap scores match expectation
  score_ab <- compute_overlap(10000000L, 100000000L,
                              10000000L, 65000000L,  "reciprocal")
  score_ac <- compute_overlap(10000000L, 100000000L,
                              55000000L, 100000000L, "reciprocal")
  score_bc <- compute_overlap(10000000L, 65000000L,
                              55000000L, 100000000L, "reciprocal")
  
  # Confirm geometry before testing equivalence
  testthat::expect_gte(score_ab, 0.5)   # A+B pass
  testthat::expect_gte(score_ac, 0.5)   # A+C pass
  testthat::expect_lt(score_bc,  0.5)   # B+C fail
  
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.5,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  
  # Two cliques → two distinct equiv IDs
  testthat::expect_equal(
    dplyr::n_distinct(out$results_id$cnv_equiv_id), 2L
  )
  # A appears in both cliques → more rows than input
  testthat::expect_gt(nrow(out$results_id), nrow(df))
})

testthat::test_that("equivalence scenario 4: no overlap → all isolated → all removed", {
  df  <- make_scenario_no_overlap()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  # All isolated — no results
  testthat::expect_equal(nrow(out$results_id), 0L)
  # All in removed_log
  testthat::expect_equal(nrow(out$removed_log), nrow(df))
})

testthat::test_that("equivalence scenario 5: single reference → removed as single_sequence", {
  df  <- make_scenario_single_ref()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  testthat::expect_equal(nrow(out$results_id), 0L)
  testthat::expect_true(
    any(out$removed_log$removal_reason == "single_sequence")
  )
})

testthat::test_that("equivalence scenario 6: A+B pass threshold, A+C fail → one clique", {
  df  <- make_scenario_ab_pass_ac_fail()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  # Only A+B form a clique
  testthat::expect_equal(
    dplyr::n_distinct(out$results_id$cnv_equiv_id), 1L
  )
  summary <- summarize_cnv_support(out$results_id)
  testthat::expect_equal(summary$n_references, 2L)
})

# ---- General assign_cnv_equivalence properties ----------------------------

testthat::test_that("assign_cnv_equivalence: no NA equiv IDs in results", {
  df  <- make_scenario_perfect_overlap()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  testthat::expect_false(any(is.na(out$results_id$cnv_equiv_id)))
})

testthat::test_that("assign_cnv_equivalence: missing required columns raises error", {
  df      <- make_scenario_perfect_overlap()
  df$chr  <- NULL
  testthat::expect_error(
    assign_cnv_equivalence(df, filter_seq_mb = 0),
    "Missing required columns"
  )
})

testthat::test_that("assign_cnv_equivalence: start > end raises error", {
  df        <- make_scenario_perfect_overlap()
  df$end[1] <- df$start[1] - 1L
  testthat::expect_error(
    assign_cnv_equivalence(df, filter_seq_mb = 0),
    "start > end"
  )
})

testthat::test_that("assign_cnv_equivalence: returns list with results_id and removed_log", {
  df  <- make_scenario_perfect_overlap()
  out <- assign_cnv_equivalence(
    df,
    min_overlap    = 0.75,
    overlap_method = "reciprocal",
    filter_seq_mb  = 0
  )
  testthat::expect_true("results_id"  %in% names(out))
  testthat::expect_true("removed_log" %in% names(out))
})

testthat::test_that("assign_cnv_equivalence: filter_seq_mb removes all segments raises error", {
  df <- make_scenario_perfect_overlap()
  testthat::expect_error(
    assign_cnv_equivalence(
      df,
      min_overlap    = 0.75,
      overlap_method = "reciprocal",
      filter_seq_mb  = 50
    ),
    "No segments remain"
  )
})

# ============================================================================
# add_metadata
# ============================================================================

testthat::test_that("add_metadata: adds mode column", {
  df  <- make_mock_supported_events()
  df  <- df[df$cell_name %in% make_mock_metadata()$cell_name, ]
  # Use subset that matches metadata
  df2 <- make_mock_supported_events()[1, ]
  df2$cell_name <- "E5.5.101"
  out <- add_metadata(df2, mode_name = "within", metadata = make_mock_metadata())
  testthat::expect_true("mode" %in% colnames(out))
  testthat::expect_equal(out$mode, "within")
})

testthat::test_that("add_metadata: joins cell_type from metadata", {
  df <- data.frame(
    cell_name = "E5.5.101",
    chr       = "chr1",
    cnv_state = "gain",
    stringsAsFactors = FALSE
  )
  out <- add_metadata(df, mode_name = "within", metadata = make_mock_metadata())
  testthat::expect_true("cell_type" %in% colnames(out))
  testthat::expect_equal(out$cell_type, "TE")
})

testthat::test_that("add_metadata: unmatched cell_names raise error", {
  df <- data.frame(
    cell_name = "UNKNOWN_CELL",
    chr       = "chr1",
    cnv_state = "gain",
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    add_metadata(df, mode_name = "within", metadata = make_mock_metadata()),
    "empty cell type values"
  )
})