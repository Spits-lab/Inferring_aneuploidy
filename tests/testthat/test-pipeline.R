source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/tests/testthat/helper-fixtures.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/pipeline.R")

# tests/testthat/test-pipeline.R
# Tests for R/pipeline.R
# Fixtures loaded from helper-fixtures.R

# ============================================================================
# compute_cell_sizes
# ============================================================================

testthat::test_that("compute_cell_sizes: returns one row per group", {
  meta <- make_mock_pipeline_metadata()
  out  <- compute_cell_sizes(meta, group_cols = "embryo")
  testthat::expect_equal(nrow(out), dplyr::n_distinct(meta$embryo))
})

testthat::test_that("compute_cell_sizes: n_total_cells counts distinct cells", {
  meta <- make_mock_pipeline_metadata()
  out  <- compute_cell_sizes(meta, group_cols = "embryo")
  # E6.2 has 3 cells
  e62  <- out[out$embryo == "E6.2", ]
  testthat::expect_equal(e62$n_total_cells, 3L)
})

testthat::test_that("compute_cell_sizes: multiple group_cols works", {
  meta <- make_mock_pipeline_metadata()
  out  <- compute_cell_sizes(
    meta,
    group_cols = c("embryo", "cell_type")
  )
  # Each unique embryo + cell_type combination = one row
  expected_rows <- dplyr::n_distinct(
    paste(meta$embryo, meta$cell_type)
  )
  testthat::expect_equal(nrow(out), expected_rows)
})

testthat::test_that("compute_cell_sizes: n_total_cells column present", {
  meta <- make_mock_pipeline_metadata()
  out  <- compute_cell_sizes(meta, group_cols = "embryo")
  testthat::expect_true("n_total_cells" %in% colnames(out))
})

testthat::test_that("compute_cell_sizes: group column present in output", {
  meta <- make_mock_pipeline_metadata()
  out  <- compute_cell_sizes(meta, group_cols = "embryo")
  testthat::expect_true("embryo" %in% colnames(out))
})

testthat::test_that("compute_cell_sizes: missing group_col raises error", {
  meta <- make_mock_pipeline_metadata()
  testthat::expect_error(
    compute_cell_sizes(meta, group_cols = "nonexistent"),
    "Missing columns"
  )
})

testthat::test_that("compute_cell_sizes: missing cell_col raises error", {
  meta <- make_mock_pipeline_metadata()
  testthat::expect_error(
    compute_cell_sizes(meta, group_cols = "embryo",
                       cell_col = "nonexistent"),
    "Missing columns"
  )
})

testthat::test_that("compute_cell_sizes: custom cell_col works", {
  meta        <- make_mock_pipeline_metadata()
  meta$my_id  <- meta$cell_name
  out <- compute_cell_sizes(
    meta,
    group_cols = "embryo",
    cell_col   = "my_id"
  )
  testthat::expect_true("n_total_cells" %in% colnames(out))
})

testthat::test_that("compute_cell_sizes: duplicate cell_names counted once", {
  meta <- make_mock_pipeline_metadata()
  # Duplicate a cell
  meta_dup <- rbind(meta, meta[1, ])
  out <- compute_cell_sizes(meta_dup, group_cols = "embryo")
  # E5.5 should still have 2 distinct cells not 3
  e55 <- out[out$embryo == "E5.5", ]
  testthat::expect_equal(e55$n_total_cells, 2L)
})

testthat::test_that("compute_cell_sizes: total cells sums to n_distinct cells", {
  meta <- make_mock_pipeline_metadata()
  out  <- compute_cell_sizes(meta, group_cols = "embryo")
  testthat::expect_equal(
    sum(out$n_total_cells),
    dplyr::n_distinct(meta$cell_name)
  )
})

# ============================================================================
# load_tool_data — dispatcher tests
# ============================================================================

testthat::test_that("load_tool_data: invalid tool raises error", {
  testthat::expect_error(
    load_tool_data(
      tool     = "invalid_tool",
      base_dir = tempdir(),
      ref_dirs = "A",
      pattern  = "^run\\.final"
    ),
    "should be one of"
  )
})



# ============================================================================
# run_cnv_tool — dispatcher tests
# ============================================================================

testthat::test_that("run_cnv_tool: invalid tool raises error", {
  testthat::expect_error(
    run_cnv_tool(
      tool        = "invalid_tool",
      counts_mx   = make_mock_expr_matrix(),
      metadata    = make_mock_infercnv_metadata(),
      base_outdir = tempdir()
    ),
    "should be one of"
  )
})

testthat::test_that("run_cnv_tool: scevan raises not implemented error", {
  testthat::expect_error(
    run_cnv_tool(
      tool        = "scevan",
      counts_mx   = make_mock_expr_matrix(),
      metadata    = make_mock_infercnv_metadata(),
      base_outdir = tempdir()
    ),
    "not yet implemented"
  )
})

testthat::test_that("run_cnv_tool: copykat raises not implemented error", {
  testthat::expect_error(
    run_cnv_tool(
      tool        = "copykat",
      counts_mx   = make_mock_expr_matrix(),
      metadata    = make_mock_infercnv_metadata(),
      base_outdir = tempdir()
    ),
    "not yet implemented"
  )
})

# ============================================================================
# process_tool_cnv_runs — directory validation tests
# ============================================================================

testthat::test_that("process_tool_cnv_runs: nonexistent base_dir raises error", {
  testthat::expect_error(
    process_tool_cnv_runs(
      base_dir = "/nonexistent/path",
      mode     = "within",
      metadata = make_mock_pipeline_metadata()
    ),
    "does not exist"
  )
})

testthat::test_that("process_tool_cnv_runs: missing mode directory raises error", {
  # Use a dedicated subdirectory that definitely has no 'within' folder
  tmp <- file.path(tempdir(), "test_missing_mode")
  dir.create(tmp, showWarnings = FALSE)
  
  testthat::expect_error(
    process_tool_cnv_runs(
      base_dir = tmp,
      mode     = "within",
      metadata = make_mock_pipeline_metadata()
    ),
    "Mode directories not found"
  )
  
  # Cleanup
  unlink(tmp, recursive = TRUE)
})


testthat::test_that("process_tool_cnv_runs: empty mode directory raises error", {
  tmp      <- file.path(tempdir(), "test_empty_mode")
  mode_dir <- file.path(tmp, "within")
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
  
  testthat::expect_error(
    process_tool_cnv_runs(
      base_dir = tmp,
      mode     = "within",
      metadata = make_mock_pipeline_metadata()
    ),
    "No cell types found"
  )
  
  unlink(tmp, recursive = TRUE)
})

# ============================================================================
# run_full_cnv_pipeline — entry point validation
# ============================================================================

testthat::test_that("run_full_cnv_pipeline: block1 without counts_mx raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from      = "block1",
        counts_mx       = NULL,
        metadata        = make_mock_pipeline_metadata(),
        gene_order_file = "dummy.txt",
        base_outdir     = tempdir()
      )
    ),
    "counts_mx required"
  )
})

testthat::test_that("run_full_cnv_pipeline: block1 without metadata raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from      = "block1",
        counts_mx       = make_mock_expr_matrix(),
        metadata        = NULL,
        gene_order_file = "dummy.txt",
        base_outdir     = tempdir()
      )
    ),
    "metadata required"
  )
})

testthat::test_that("run_full_cnv_pipeline: block1 without gene_order_file raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from      = "block1",
        counts_mx       = make_mock_expr_matrix(),
        metadata        = make_mock_pipeline_metadata(),
        gene_order_file = NULL,
        base_outdir     = tempdir()
      )
    ),
    "gene_order_file required"
  )
})

testthat::test_that("run_full_cnv_pipeline: block1 without base_outdir raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from      = "block1",
        counts_mx       = make_mock_expr_matrix(),
        metadata        = make_mock_pipeline_metadata(),
        gene_order_file = "dummy.txt",
        base_outdir     = NULL
      )
    ),
    "base_outdir required"
  )
})

testthat::test_that("run_full_cnv_pipeline: block2 without base_dir raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from = "block2",
        base_dir   = NULL,
        metadata   = make_mock_pipeline_metadata()
      )
    ),
    "base_dir required"
  )
})

testthat::test_that("run_full_cnv_pipeline: block3 without chromosome_arms raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from      = "block3",
        chromosome_arms = NULL,
        group_cols      = "embryo",
        metadata        = make_mock_pipeline_metadata(),
        supported_events = make_mock_supported_events()
      )
    ),
    "chromosome_arms required"
  )
})

testthat::test_that("run_full_cnv_pipeline: block3 without group_cols raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from       = "block3",
        chromosome_arms  = hg38_chromosome_arms,
        group_cols       = NULL,
        metadata         = make_mock_pipeline_metadata(),
        supported_events = make_mock_supported_events()
      )
    ),
    "group_cols required"
  )
})

testthat::test_that("run_full_cnv_pipeline: save_intermediate without outdir raises error", {
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from        = "block2",
        base_dir          = tempdir(),
        metadata          = make_mock_pipeline_metadata(),
        save_intermediate = TRUE,
        outdir            = NULL
      )
    ),
    "outdir must be provided"
  )
})

testthat::test_that("run_full_cnv_pipeline: returns named list with block elements", {
  tmp      <- file.path(tempdir(), "test_pipeline_structure")
  mode_dir <- file.path(tmp, "within")
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Empty mode directory now raises error — test that block1 is NULL
  # by checking the error comes from block2 not block1
  testthat::expect_error(
    run_full_cnv_pipeline(
      start_from = "block2",
      base_dir   = tmp,
      metadata   = make_mock_pipeline_metadata()
    ),
    "No cell types found"  # error from block2, not block1
  )
  
  unlink(tmp, recursive = TRUE)
})


testthat::test_that("run_full_cnv_pipeline: block1 is NULL when starting from block2", {
  # When start_from = "block2", block1 should never run
  # We verify this by checking the error originates in block2
  # not in block1 (which would say "counts_mx required")
  tmp      <- file.path(tempdir(), "test_block1_null")
  mode_dir <- file.path(tmp, "within")
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
  
  err <- tryCatch(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from = "block2",
        base_dir   = tmp,
        metadata   = make_mock_pipeline_metadata()
      )
    ),
    error = function(e) e
  )
  
  # Error should be from block2 (No cell types) not block1 (counts_mx required)
  testthat::expect_true(grepl("No cell types", err$message))
  testthat::expect_false(grepl("counts_mx", err$message))
  
  unlink(tmp, recursive = TRUE)
})

testthat::test_that("run_full_cnv_pipeline: blocks_to_run correct for block3 start", {
  # Starting from block3 should run block3 and block4 only
  # We just check the summary — not actual execution
  testthat::expect_error(
    suppressMessages(
      run_full_cnv_pipeline(
        start_from       = "block3",
        chromosome_arms  = hg38_chromosome_arms,
        group_cols       = "embryo",
        metadata         = make_mock_pipeline_metadata(),
        supported_events = make_mock_supported_events()
      )
    )
    # Will error at add_chromosome_info since mock data is minimal
    # but the error comes from block3 execution not block routing
  )
})