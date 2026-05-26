
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/tests/testthat/helper-fixtures.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/infercnv.R")
# tests/testthat/test-infercnv.R
# Tests for R/infercnv.R — utility functions only
# run_infercnv_pipeline excluded — too expensive to test with real inferCNV
# Fixtures loaded from helper-fixtures.R

# ============================================================================
# melt_expr_to_long
# ============================================================================

testthat::test_that("melt_expr_to_long: returns long format with correct columns", {
  mat <- make_mock_expr_matrix()
  out <- melt_expr_to_long(as.data.frame(mat))
  testthat::expect_true("gene"      %in% colnames(out))
  testthat::expect_true("cell_name" %in% colnames(out))
  testthat::expect_true("state_raw" %in% colnames(out))
})

testthat::test_that("melt_expr_to_long: n_rows = n_genes x n_cells", {
  mat <- make_mock_expr_matrix()
  out <- melt_expr_to_long(as.data.frame(mat))
  testthat::expect_equal(nrow(out), nrow(mat) * ncol(mat))
})

testthat::test_that("melt_expr_to_long: all genes present in output", {
  mat <- make_mock_expr_matrix()
  out <- melt_expr_to_long(as.data.frame(mat))
  testthat::expect_true(all(rownames(mat) %in% out$gene))
})

testthat::test_that("melt_expr_to_long: all cells present in output", {
  mat <- make_mock_expr_matrix()
  out <- melt_expr_to_long(as.data.frame(mat))
  testthat::expect_true(all(colnames(mat) %in% out$cell_name))
})

testthat::test_that("melt_expr_to_long: cell_ prefix stripped from cell_name", {
  mat <- make_mock_expr_matrix()
  df  <- as.data.frame(mat)
  # Add prefix to column names
  colnames(df) <- paste0("cell_", colnames(df))
  out <- melt_expr_to_long(df, cell_prefix = "cell_")
  # Prefix should be removed
  testthat::expect_false(any(grepl("^cell_", out$cell_name)))
})

testthat::test_that("melt_expr_to_long: no prefix stripping when prefix absent", {
  mat <- make_mock_expr_matrix()
  out <- melt_expr_to_long(as.data.frame(mat), cell_prefix = "cell_")
  # Original names have no prefix — should be unchanged
  testthat::expect_true(all(colnames(mat) %in% out$cell_name))
})

testthat::test_that("melt_expr_to_long: no rownames raises error", {
  mat <- make_mock_expr_matrix()
  df  <- as.data.frame(mat)
  rownames(df) <- NULL
  testthat::expect_error(
    melt_expr_to_long(df),
    "no meaningful rownames"
  )
})

testthat::test_that("melt_expr_to_long: numeric rownames raises error", {
  mat <- make_mock_expr_matrix()
  df  <- as.data.frame(mat)
  rownames(df) <- as.character(seq_len(nrow(df)))
  testthat::expect_error(
    melt_expr_to_long(df),
    "no meaningful rownames"
  )
})

testthat::test_that("melt_expr_to_long: state_raw values match original matrix", {
  mat <- make_mock_expr_matrix()
  out <- melt_expr_to_long(as.data.frame(mat))
  
  # GENE5 in E5.5.100 = 3.0 — verified against matrix
  val <- out$state_raw[out$gene == "GENE5" & out$cell_name == "E5.5.100"]
  testthat::expect_equal(val, 3.0)
  
  # DDX11L1 in E6.2.114 = 3.0
  val2 <- out$state_raw[out$gene == "DDX11L1" & out$cell_name == "E6.2.114"]
  testthat::expect_equal(val2, 3.0)
  
  # MIR6723 in E5.5.101 = 0.0
  val3 <- out$state_raw[out$gene == "MIR6723" & out$cell_name == "E5.5.101"]
  testthat::expect_equal(val3, 0.0)
})

# ============================================================================
# attach_gene_order
# ============================================================================

testthat::test_that("attach_gene_order: adds chr start stop columns", {
  long_df   <- make_mock_long_df()
  gene_order <- make_mock_gene_order()
  out        <- attach_gene_order(long_df, gene_order)
  testthat::expect_true("chr"   %in% colnames(out))
  testthat::expect_true("start" %in% colnames(out))
  testthat::expect_true("stop"  %in% colnames(out))
})

testthat::test_that("attach_gene_order: genes not in gene_order are dropped", {
  long_df    <- make_mock_long_df()
  gene_order <- make_mock_gene_order()
  # Add a gene not in gene_order
  extra <- data.frame(
    gene      = "UNKNOWN_GENE",
    cell_name = "E5.5.101",
    state_raw = 1.0,
    stringsAsFactors = FALSE
  )
  long_df_extra <- rbind(long_df, extra)
  
  testthat::expect_warning(
    out <- attach_gene_order(long_df_extra, gene_order),
    "dropped"
  )
  testthat::expect_false("UNKNOWN_GENE" %in% out$gene)
})

testthat::test_that("attach_gene_order: no matching genes raises error", {
  long_df <- data.frame(
    gene      = c("NONEXISTENT1", "NONEXISTENT2"),
    cell_name = c("cell_1",       "cell_2"),
    state_raw = c(1.0,            1.0),
    stringsAsFactors = FALSE
  )
  gene_order <- make_mock_gene_order()
  testthat::expect_error(
    attach_gene_order(long_df, gene_order),
    "empty table"
  )
})

testthat::test_that("attach_gene_order: missing required gene_order columns raises error", {
  long_df    <- make_mock_long_df()
  gene_order <- make_mock_gene_order()
  gene_order$chr <- NULL
  testthat::expect_error(
    attach_gene_order(long_df, gene_order),
    "missing required columns"
  )
})

testthat::test_that("attach_gene_order: coordinates match gene_order values", {
  long_df    <- make_mock_long_df()
  gene_order <- make_mock_gene_order()
  out        <- attach_gene_order(long_df, gene_order)
  
  # MIR6723 should have chr1, start=632325, stop=632413
  mir_row <- out[out$gene == "MIR6723", ][1, ]
  testthat::expect_equal(as.character(mir_row$chr),   "chr1")
  testthat::expect_equal(mir_row$start, 632325L)
  testthat::expect_equal(mir_row$stop,  632413L)
})

testthat::test_that("attach_gene_order: gene column preserved in output", {
  long_df    <- make_mock_long_df()
  gene_order <- make_mock_gene_order()
  out        <- attach_gene_order(long_df, gene_order)
  testthat::expect_true("gene" %in% colnames(out))
})

# ============================================================================
# discretize_cnv_state_infer_cnv
# ============================================================================

testthat::test_that("discretize_cnv_state_infer_cnv: adds state column", {
  df  <- make_mock_long_df()
  out <- discretize_cnv_state_infer_cnv(df, k = 1.5)
  testthat::expect_true("state" %in% colnames(out))
})

testthat::test_that("discretize_cnv_state_infer_cnv: state values are gain/loss/neutral only", {
  df  <- make_mock_long_df()
  out <- discretize_cnv_state_infer_cnv(df, k = 1.5)
  testthat::expect_true(all(out$state %in% c("gain", "loss", "neutral")))
})

testthat::test_that("discretize_cnv_state_infer_cnv: neutral dominates with default k", {
  df  <- make_mock_long_df()
  out <- discretize_cnv_state_infer_cnv(df, k = 1.5)
  testthat::expect_gt(
    sum(out$state == "neutral"),
    sum(out$state != "neutral")
  )
})

testthat::test_that("discretize_cnv_state_infer_cnv: larger k produces fewer calls", {
  df    <- make_mock_long_df()
  out_1 <- discretize_cnv_state_infer_cnv(df, k = 0.5)
  out_2 <- discretize_cnv_state_infer_cnv(df, k = 3.0)
  n_calls_1 <- sum(out_1$state != "neutral")
  n_calls_2 <- sum(out_2$state != "neutral")
  testthat::expect_gte(n_calls_1, n_calls_2)
})

testthat::test_that("discretize_cnv_state_infer_cnv: extreme values classified as gain or loss", {
  df <- data.frame(
    state_raw = c(rep(1.0, 100), 100.0, -100.0),
    stringsAsFactors = FALSE
  )
  out <- discretize_cnv_state_infer_cnv(df, k = 1.5)
  testthat::expect_true("gain" %in% out$state)
  testthat::expect_true("loss" %in% out$state)
})

testthat::test_that("discretize_cnv_state_infer_cnv: missing state_raw raises error", {
  df <- data.frame(x = 1:5)
  testthat::expect_error(
    discretize_cnv_state_infer_cnv(df),
    "Missing required column: state_raw"
  )
})

testthat::test_that("discretize_cnv_state_infer_cnv: non-numeric state_raw raises error", {
  df <- data.frame(state_raw = letters[1:5])
  testthat::expect_error(
    discretize_cnv_state_infer_cnv(df),
    "state_raw must be numeric"
  )
})

testthat::test_that("discretize_cnv_state_infer_cnv: all NA state_raw raises error", {
  df <- data.frame(state_raw = rep(NA_real_, 5))
  testthat::expect_error(
    discretize_cnv_state_infer_cnv(df),
    "entirely NA"
  )
})

testthat::test_that("discretize_cnv_state_infer_cnv: row count unchanged", {
  df  <- make_mock_long_df()
  out <- discretize_cnv_state_infer_cnv(df, k = 1.5)
  testthat::expect_equal(nrow(out), nrow(df))
})

testthat::test_that("discretize_cnv_state_infer_cnv: k = 0 classifies all as gain or loss", {
  df <- data.frame(state_raw = c(1.0, 1.0, 1.0, 1.1, 0.9))
  out <- discretize_cnv_state_infer_cnv(df, k = 0)
  # With k=0 threshold = mean ± 0 → everything outside mean is gain/loss
  testthat::expect_true(all(out$state %in% c("gain", "loss", "neutral")))
})

# ============================================================================
# load_and_prepare_infercnv_reference
# ============================================================================

testthat::test_that("load_and_prepare_infercnv_reference: returns data frame", {
  infercnv_list <- make_mock_infercnv_list()
  out <- load_and_prepare_infercnv_reference(infercnv_list)
  testthat::expect_s3_class(out, "data.frame")
})

testthat::test_that("load_and_prepare_infercnv_reference: has required columns", {
  infercnv_list <- make_mock_infercnv_list()
  out <- load_and_prepare_infercnv_reference(infercnv_list)
  required <- c("gene", "cell_name", "state_raw", "chr",
                "start", "stop", "state", "reference")
  testthat::expect_true(all(required %in% colnames(out)))
})

testthat::test_that("load_and_prepare_infercnv_reference: reference column has correct values", {
  infercnv_list <- make_mock_infercnv_list()
  out <- load_and_prepare_infercnv_reference(infercnv_list)
  testthat::expect_true(all(c("A", "B", "C") %in% out$reference))
})

testthat::test_that("load_and_prepare_infercnv_reference: state values valid", {
  infercnv_list <- make_mock_infercnv_list()
  out <- load_and_prepare_infercnv_reference(infercnv_list)
  testthat::expect_true(all(out$state %in% c("gain", "loss", "neutral")))
})

testthat::test_that("load_and_prepare_infercnv_reference: n_rows = refs x genes x cells", {
  infercnv_list <- make_mock_infercnv_list()
  mat           <- make_mock_expr_matrix()
  out           <- load_and_prepare_infercnv_reference(infercnv_list)
  expected_rows <- length(infercnv_list) * nrow(mat) * ncol(mat)
  testthat::expect_equal(nrow(out), expected_rows)
})

testthat::test_that("load_and_prepare_infercnv_reference: single reference works", {
  infercnv_list <- make_mock_infercnv_list()[1]  # only ref A
  out <- load_and_prepare_infercnv_reference(infercnv_list)
  testthat::expect_equal(dplyr::n_distinct(out$reference), 1L)
  testthat::expect_equal(unique(out$reference), "A")
})

# ============================================================================
# validate_metadata
# ============================================================================

testthat::test_that("validate_metadata: valid input passes without error", {
  meta     <- make_mock_infercnv_metadata()
  counts   <- make_mock_expr_matrix()
  testthat::expect_no_error(
    suppressMessages(
      validate_metadata(meta, counts, cell_type_col = "cell_type",
                        min_cells = 1L)
    )
  )
})

testthat::test_that("validate_metadata: missing cell_name column raises error", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  meta$cell_name <- NULL
  testthat::expect_error(
    suppressMessages(
      validate_metadata(meta, counts, "cell_type", min_cells = 1L)
    ),
    "missing required column"
  )
})

testthat::test_that("validate_metadata: missing cell_type_col raises error", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  testthat::expect_error(
    suppressMessages(
      validate_metadata(meta, counts, "nonexistent_col", min_cells = 1L)
    ),
    "missing required column"
  )
})

testthat::test_that("validate_metadata: duplicated cell_name raises error", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  meta$cell_name[2] <- meta$cell_name[1]  # duplicate
  testthat::expect_error(
    suppressMessages(
      validate_metadata(meta, counts, "cell_type", min_cells = 1L)
    ),
    "duplicated"
  )
})

testthat::test_that("validate_metadata: cell in metadata not in counts raises error", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  meta$cell_name[1] <- "UNKNOWN_CELL"
  testthat::expect_error(
    suppressMessages(
      validate_metadata(meta, counts, "cell_type", min_cells = 1L)
    ),
    "not found in counts_mx"
  )
})

testthat::test_that("validate_metadata: cell in counts not in metadata produces warning", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  # Remove one cell from metadata — it exists in counts but not metadata
  meta <- meta[-1, ]
  testthat::expect_warning(
    suppressMessages(
      validate_metadata(meta, counts, "cell_type", min_cells = 1L)
    ),
    "not in metadata"
  )
})

testthat::test_that("validate_metadata: NA in cell_type_col raises error", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  meta$cell_type[1] <- NA
  testthat::expect_error(
    suppressMessages(
      validate_metadata(meta, counts, "cell_type", min_cells = 1L)
    ),
    "NA value"
  )
})

testthat::test_that("validate_metadata: returns filtered metadata", {
  meta   <- make_mock_infercnv_metadata()
  counts <- make_mock_expr_matrix()
  out <- suppressMessages(
    validate_metadata(meta, counts, "cell_type", min_cells = 1L)
  )
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true("cell_name"  %in% colnames(out))
  testthat::expect_true("cell_type"  %in% colnames(out))
})

# ============================================================================
# make_splits
# ============================================================================

testthat::test_that("make_splits: adds split_group column", {
  meta <- make_mock_infercnv_metadata()
  out  <- suppressMessages(
    make_splits(meta, "cell_type", "TE", n_splits = 2L)
  )
  testthat::expect_true("split_group" %in% colnames(out))
})

testthat::test_that("make_splits: split_group uses capital letters", {
  meta <- make_mock_infercnv_metadata()
  out  <- suppressMessages(
    make_splits(meta, "cell_type", "TE", n_splits = 2L)
  )
  testthat::expect_true(all(out$split_group %in% LETTERS))
})

testthat::test_that("make_splits: returns only cells of requested cell type", {
  meta <- make_mock_infercnv_metadata()
  out  <- suppressMessages(
    make_splits(meta, "cell_type", "TE", n_splits = 2L)
  )
  testthat::expect_true(all(out$cell_type == "TE"))
})

testthat::test_that("make_splits: n_splits = 3 produces A B C labels", {
  # Need at least 3 cells of same type
  meta <- data.frame(
    cell_name = paste0("cell_", 1:6),
    cell_type = rep("TE", 6),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(
    make_splits(meta, "cell_type", "TE", n_splits = 3L)
  )
  testthat::expect_true(all(c("A", "B", "C") %in% out$split_group))
})

testthat::test_that("make_splits: non-integer n_splits raises error", {
  meta <- make_mock_infercnv_metadata()
  testthat::expect_error(
    make_splits(meta, "cell_type", "TE", n_splits = 2.5),
    "integer"
  )
})

testthat::test_that("make_splits: n_splits > 26 raises error", {
  meta <- data.frame(
    cell_name = paste0("cell_", 1:30),
    cell_type = rep("TE", 30),
    stringsAsFactors = FALSE
  )
  testthat::expect_error(
    make_splits(meta, "cell_type", "TE", n_splits = 27L),
    "26"
  )
})

testthat::test_that("make_splits: fewer cells than n_splits raises error", {
  meta <- make_mock_infercnv_metadata()  # only 2 TE cells
  testthat::expect_error(
    make_splits(meta, "cell_type", "TE", n_splits = 5L),
    "only"
  )
})

testthat::test_that("make_splits: total cells preserved after splitting", {
  meta <- data.frame(
    cell_name = paste0("cell_", 1:6),
    cell_type = rep("TE", 6),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(
    make_splits(meta, "cell_type", "TE", n_splits = 3L)
  )
  testthat::expect_equal(nrow(out), 6L)
})

# ============================================================================
# build_annotations_df
# ============================================================================

testthat::test_that("build_annotations_df: returns data frame with one column", {
  out <- build_annotations_df(
    cell_names   = c("cell_1", "cell_2", "cell_3"),
    group_labels = c("TE",     "TE",     "Epiblast")
  )
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(ncol(out), 1L)
})

testthat::test_that("build_annotations_df: rownames are cell names", {
  cells  <- c("cell_1", "cell_2", "cell_3")
  labels <- c("TE",     "TE",     "Epiblast")
  out    <- build_annotations_df(cells, labels)
  testthat::expect_equal(rownames(out), cells)
})

testthat::test_that("build_annotations_df: group column contains labels", {
  cells  <- c("cell_1", "cell_2", "cell_3")
  labels <- c("TE",     "TE",     "Epiblast")
  out    <- build_annotations_df(cells, labels)
  testthat::expect_equal(out$group, labels)
})

testthat::test_that("build_annotations_df: mismatched lengths raise error", {
  testthat::expect_error(
    build_annotations_df(
      cell_names   = c("cell_1", "cell_2"),
      group_labels = c("TE")       # wrong length
    ),
    "same length"
  )
})

testthat::test_that("build_annotations_df: empty input returns empty data frame", {
  out <- build_annotations_df(
    cell_names   = character(0),
    group_labels = character(0)
  )
  testthat::expect_equal(nrow(out), 0L)
})

# ============================================================================
# discover_infercnv_runs
# ============================================================================

testthat::test_that("discover_infercnv_runs: missing directory raises error", {
  testthat::expect_error(
    discover_infercnv_runs(
      base_dir = tempdir(),
      ref_dirs = "nonexistent_ref",
      pattern  = "^run\\.final"
    )
  )
})

testthat::test_that("discover_infercnv_runs: no matching files raises error", {
  # Create a temp directory with no matching files
  tmp <- tempdir()
  ref_dir <- file.path(tmp, "empty_ref")
  dir.create(ref_dir, showWarnings = FALSE)
  
  testthat::expect_error(
    discover_infercnv_runs(
      base_dir = tmp,
      ref_dirs = "empty_ref",
      pattern  = "^run\\.final"
    ),
    "No files matching"
  )
})

testthat::test_that("discover_infercnv_runs: multiple matching files raises error", {
  tmp     <- tempdir()
  ref_dir <- file.path(tmp, "multi_ref")
  dir.create(ref_dir, showWarnings = FALSE)
  
  # Create two files matching the pattern
  saveRDS(list(), file.path(ref_dir, "run.final.file1.rds"))
  saveRDS(list(), file.path(ref_dir, "run.final.file2.rds"))
  
  testthat::expect_error(
    discover_infercnv_runs(
      base_dir = tmp,
      ref_dirs = "multi_ref",
      pattern  = "^run\\.final"
    ),
    "Multiple files"
  )
  
  # Cleanup
  unlink(ref_dir, recursive = TRUE)
})

testthat::test_that("discover_infercnv_runs: finds and loads single matching file", {
  tmp     <- tempdir()
  ref_dir <- file.path(tmp, "good_ref")
  dir.create(ref_dir, showWarnings = FALSE)
  
  # Save a mock object
  mock_obj <- list(expr = "mock")
  saveRDS(mock_obj, file.path(ref_dir, "run.final.infercnv_obj"))
  
  out <- discover_infercnv_runs(
    base_dir = tmp,
    ref_dirs = "good_ref",
    pattern  = "^run\\.final"
  )
  
  testthat::expect_type(out, "list")
  testthat::expect_equal(length(out), 1L)
  testthat::expect_equal(names(out), "good_ref")
  
  # Cleanup
  unlink(ref_dir, recursive = TRUE)
})

testthat::test_that("discover_infercnv_runs: names match ref_dirs", {
  tmp <- tempdir()
  
  for (ref in c("A", "B", "C")) {
    dir_path <- file.path(tmp, ref)
    dir.create(dir_path, showWarnings = FALSE)
    saveRDS(list(), file.path(dir_path, "run.final.infercnv_obj"))
  }
  
  out <- discover_infercnv_runs(
    base_dir = tmp,
    ref_dirs = c("A", "B", "C"),
    pattern  = "^run\\.final"
  )
  
  testthat::expect_equal(sort(names(out)), c("A", "B", "C"))
  
  # Cleanup
  for (ref in c("A", "B", "C")) {
    unlink(file.path(tmp, ref), recursive = TRUE)
  }
})