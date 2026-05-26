#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R
# Reference-independent CNV karyotyping pipeline
# Compatible with SLURM and Snakemake (Pattern A — optparse)
#
# Usage:
#   Rscript run_pipeline.R --base_dir /path/to/infercnv --metadata_path /path/to/metadata.rds
#   Rscript run_pipeline.R --help
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# =============================================================================
# Define options
# =============================================================================

option_list <- list(
  
  # ---- General ------------------------------------------------------------
  optparse::make_option(
    "--start_from",
    type    = "character",
    default = "block2",
    help    = "Block to start from: block1, block2, block3, block4 [default: %default]"
  ),
  optparse::make_option(
    "--outdir",
    type    = "character",
    default = NULL,
    help    = "Output directory for pipeline results [required]"
  ),
  optparse::make_option(
    "--save_intermediate",
    action  = "store_true",
    default = FALSE,
    help    = "Save intermediate block outputs to outdir [default: %default]"
  ),
  
  # ---- Metadata -----------------------------------------------------------
  optparse::make_option(
    "--metadata_path",
    type    = "character",
    default = NULL,
    help    = "Path to metadata RDS file [required]"
  ),
  
  # ---- Block 2 ------------------------------------------------------------
  optparse::make_option(
    "--base_dir",
    type    = "character",
    default = NULL,
    help    = "Path to inferCNV output folder [required for block2]"
  ),
  optparse::make_option(
    "--tool",
    type    = "character",
    default = "infercnv",
    help    = "CNV tool: infercnv, scevan, copykat [default: %default]"
  ),
  optparse::make_option(
    "--pattern",
    type    = "character",
    default = "^run\\.final",
    help    = "Regex pattern to match inferCNV result files [default: %default]"
  ),
  optparse::make_option(
    "--max_gap",
    type    = "integer",
    default = 100000L,
    help    = "Maximum gap for merging nearby segments [default: %default]"
  ),
  optparse::make_option(
    "--min_overlap_consistent_calls",
    type    = "double",
    default = 0.75,
    help    = "Minimum overlap for equivalence assignment [default: %default]"
  ),
  optparse::make_option(
    "--min_overlap_multiple_nodes",
    type    = "double",
    default = 0.6,
    help    = "Minimum overlap for multi-node merging [default: %default]"
  ),
  optparse::make_option(
    "--filter_seq_mb_init",
    type    = "double",
    default = 5.0,
    help    = "Minimum segment length before merging in Mb [default: %default]"
  ),
  optparse::make_option(
    "--filter_seq_mb_equiv",
    type    = "double",
    default = 7.0,
    help    = "Minimum segment length before equivalence in Mb [default: %default]"
  ),
  optparse::make_option(
    "--min_references",
    type    = "integer",
    default = 2L,
    help    = "Minimum references required to retain a CNV [default: %default]"
  ),
  optparse::make_option(
    "--parallel",
    action  = "store_true",
    default = FALSE,
    help    = "Use parallel processing [default: %default]"
  ),
  optparse::make_option(
    "--cores",
    type    = "integer",
    default = 1L,
    help    = "Number of cores when parallel = TRUE [default: %default]"
  ),
  
  # ---- Block 3 ------------------------------------------------------------
  optparse::make_option(
    "--group_cols",
    type    = "character",
    default = "embryo",
    help    = "Comma-separated grouping columns e.g. 'embryo' or 'embryo,cell_type' [default: %default]"
  ),
  optparse::make_option(
    "--chromosome_arms_path",
    type    = "character",
    default = NULL,
    help    = "Path to chromosome arms RDS file. If NULL uses hg38 built-in [default: %default]"
  ),
  
  # ---- Block 4 ------------------------------------------------------------
  optparse::make_option(
    "--boundaries_mb",
    type    = "character",
    default = "25,10",
    help    = "Comma-separated size boundaries in Mb e.g. '25,10' [default: %default]"
  ),
  optparse::make_option(
    "--base_fraction",
    type    = "double",
    default = 0.05,
    help    = "Base fraction for frequency threshold [default: %default]"
  ),
  optparse::make_option(
    "--step",
    type    = "double",
    default = 0.05,
    help    = "Fraction step between tiers [default: %default]"
  ),
  optparse::make_option(
    "--min_cap_threshold",
    type    = "integer",
    default = 5L,
    help    = "Minimum cell count floor [default: %default]"
  ),
  optparse::make_option(
    "--max_cap_threshold",
    type    = "integer",
    default = 25L,
    help    = "Maximum cell count cap [default: %default]"
  ),
  optparse::make_option(
    "--total_chromosome_permission",
    type    = "double",
    default = 65.0,
    help    = "Centromere filter threshold percent [default: %default]"
  ),
  optparse::make_option(
    "--min_overlap",
    type    = "double",
    default = 0.8,
    help    = "Minimum overlap for locus clustering [default: %default]"
  )
)

# =============================================================================
# Parse arguments
# =============================================================================

opt <- optparse::parse_args(
  optparse::OptionParser(
    option_list = option_list,
    description = paste(
      "Reference-independent single-cell CNV karyotyping pipeline.",
      "Runs blocks 1-4 starting from --start_from.",
      "See run_full_cnv_pipeline() documentation for details."
    )
  )
)

# =============================================================================
# Validate required arguments
# =============================================================================

if (is.null(opt$outdir)) {
  stop("--outdir is required")
}

if (is.null(opt$metadata_path)) {
  stop("--metadata_path is required")
}

if (!file.exists(opt$metadata_path)) {
  stop("metadata_path not found: ", opt$metadata_path)
}

if (opt$start_from %in% c("block2") && is.null(opt$base_dir)) {
  stop("--base_dir is required when --start_from = block2")
}

# =============================================================================
# Load dependencies
# =============================================================================

message("Loading pipeline functions...")

# Source all pipeline files — adjust paths as needed
source("R/cnv_annotation.R")
source("R/cnv_processing.R")
source("R/cnv_scoring.R")
source("R/infercnv.R")
source("R/pipeline.R")

# =============================================================================
# Parse vector arguments
# =============================================================================

# group_cols: "embryo,cell_type" → c("embryo", "cell_type")
group_cols <- strsplit(opt$group_cols, ",")[[1]]
group_cols <- trimws(group_cols)

# boundaries_mb: "25,10" → c(25, 10)
boundaries_mb <- as.numeric(strsplit(opt$boundaries_mb, ",")[[1]])

if (any(is.na(boundaries_mb))) {
  stop(
    "Could not parse --boundaries_mb: '", opt$boundaries_mb, "'. ",
    "Expected comma-separated numbers e.g. '25,10'"
  )
}

# =============================================================================
# Load inputs
# =============================================================================

message("Loading metadata from: ", opt$metadata_path)
metadata <- readRDS(opt$metadata_path)

# Load chromosome arms
if (!is.null(opt$chromosome_arms_path)) {
  message("Loading chromosome arms from: ", opt$chromosome_arms_path)
  chromosome_arms <- readRDS(opt$chromosome_arms_path)
} else {
  message("Using built-in hg38 chromosome arms")
  chromosome_arms <- hg38_chromosome_arms
}

# =============================================================================
# Create output directory
# =============================================================================

if (!dir.exists(opt$outdir)) {
  message("Creating output directory: ", opt$outdir)
  dir.create(opt$outdir, recursive = TRUE)
}

# =============================================================================
# Run pipeline
# =============================================================================

message(sprintf(paste0(
  "\n=== Pipeline configuration ===\n",
  "  start_from:     %s\n",
  "  base_dir:       %s\n",
  "  group_cols:     %s\n",
  "  boundaries_mb:  %s\n",
  "  base_fraction:  %.3f\n",
  "  outdir:         %s"
),
opt$start_from,
opt$base_dir %||% "N/A",
paste(group_cols, collapse = ", "),
paste(boundaries_mb, collapse = ", "),
opt$base_fraction,
opt$outdir
))

results <- run_full_cnv_pipeline(
  
  start_from        = opt$start_from,
  save_intermediate = opt$save_intermediate,
  outdir            = opt$outdir,
  
  # Block 2
  base_dir                     = opt$base_dir,
  tool                         = opt$tool,
  pattern                      = opt$pattern,
  max_gap                      = opt$max_gap,
  min_overlap_consistent_calls = opt$min_overlap_consistent_calls,
  min_overlap_multiple_nodes   = opt$min_overlap_multiple_nodes,
  filter_seq_mb_init           = opt$filter_seq_mb_init,
  filter_seq_mb_equiv          = opt$filter_seq_mb_equiv,
  min_references               = opt$min_references,
  parallel                     = opt$parallel,
  cores                        = opt$cores,
  metadata                     = metadata,
  
  # Block 3
  chromosome_arms = chromosome_arms,
  group_cols      = group_cols,
  
  # Block 4
  boundaries_mb               = boundaries_mb,
  base_fraction               = opt$base_fraction,
  step                        = opt$step,
  min_cap_threshold           = opt$min_cap_threshold,
  max_cap_threshold           = opt$max_cap_threshold,
  total_chromosome_permission = opt$total_chromosome_permission,
  min_overlap                 = opt$min_overlap
)

# =============================================================================
# Save outputs
# =============================================================================

out_path <- file.path(opt$outdir, "pipeline_results.rds")
message("Saving full pipeline results to: ", out_path)
saveRDS(results, out_path)

message(sprintf(
  "\nDone. Results saved to: %s",
  out_path
))