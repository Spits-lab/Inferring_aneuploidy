#!/usr/bin/env Rscript
# =============================================================================
# run_pipeline.R
# Reference-independent CNV karyotyping pipeline
# Compatible with SLURM and Snakemake (Pattern A — optparse)
#
# Usage:
#   Single job full pipeline:
#     Rscript run_pipeline.R --execution-mode single --start-from block1 ...
#   Resume block2-4 after single job block1:
#     Rscript run_pipeline.R --execution-mode single --start-from block2 ...
#   Resume block2-4 after array job block1:
#     Rscript run_pipeline.R --execution-mode array --start-from block2 ...
#   Rscript run_pipeline.R --help
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# =============================================================================
# Resolve project root from script location
# =============================================================================

args        <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))

project_dir <- dirname(dirname(script_path)) 

# =============================================================================
# Define options
# =============================================================================

option_list <- list(
  # ---- General ------------------------------------------------------------
  optparse::make_option(
    "--execution-mode",
    type    = "character",
    default = "single",
    help    = paste(
      "Execution mode:",
      "'single' — all cell types in one job,",
      "'array'  — block1 ran as array job.",
      "[default: %default]"
    )
  ),
  optparse::make_option(
  "--remove-reference",
  type    = "logical",
  default = TRUE,
  help    = "Remove reference cells from analysis [default: %default]"
  ),
  
optparse::make_option(
  "--min-coding-density",
  type    = "double",
  default = 1.0,
  help    = "Minimal Number of Coding Genes Density [default: %default]"
  ),
  optparse::make_option(
  "--max-gap-mb",
  type    = "double",
  default = 10,
  help    = "Maximal Gap between Segments [default: %default]"
  ),
  optparse::make_option(
  "--clonal-col",
  type    = "character",
  default = NULL,
  help    = paste(
    "Metadata column defining clonal origin.",
    "Cells from the same group will never be split",
    "across reference groups (e.g. embryo, clone_id).",
    "[default: NULL — PC1 splitting only]"
  )
),
optparse::make_option(
  "--donor-col",
  type    = "character",
  default = NULL,
  help    = paste(
    "Metadata column defining donor.",
    "Used to mix donors across splits when",
    "--clonal-col is set (e.g. donor, patient).",
    "[default: NULL]"
  )
),
  
optparse::make_option(
  "--coding-genes-path",
  type    = "character",
  default = NULL,
  help    = "Path to coding genes RDS file"
),
optparse::make_option(
  "--pct-max",
  type    = "double",
  default = 45,
  help    = "Maximum pct threshold for gene density filter [default: %default]"
),
optparse::make_option(
  "--pct-floor",
  type    = "double",
  default = 30,
  help    = "Minimum pct threshold for gene density filter [default: %default]"
),
optparse::make_option(
  "--min-expr-density",
  type    = "double",
  default = 1.5,
  help    = "Minimum expressed gene density per Mb [default: %default]"
),
  optparse::make_option(
    "--start-from",
    type    = "character",
    default = "block2",
    help    = "Block to start from: block1, block2, block3, block4 [default: %default]"
  ),
  optparse::make_option(
    "--workdir",
    type    = "character",
    default = NULL,
    help    = "Output directory for pipeline results [required]"
  ),
  optparse::make_option(
    "--save-intermediate",
    action  = "store_true",
    default = FALSE,
    help    = "Save intermediate block outputs to workdir [default: %default]"
  ),
  
  # ---- Metadata -----------------------------------------------------------
  optparse::make_option(
    "--metadata-path",
    type    = "character",
    default = NULL,
    help    = paste(
      "Path to metadata RDS file.",
      "block1: full metadata with cell_type column.",
      "block2 single: split_metadata with split_group column.",
      "block2 array: optional — auto-discovered if not provided."
    )
  ),
  
  # ---- Block 1 (single mode only) -----------------------------------------
  optparse::make_option(
    "--counts-path",
    type    = "character",
    default = NULL,
    help    = "Path to counts matrix RDS [required for block1]"
  ),
  optparse::make_option(
    "--gene-order-file",
    type    = "character",
    default = NULL,
    help    = "Path to gene order file [required for block1]"
  ),
  optparse::make_option(
    "--tool-outdir",
    type    = "character",
    default = NULL,
    help    = "Base output directory for inferCNV results [required for block1]"
  ),
  optparse::make_option(
    "--group-clusters-col",
    type    = "character",
    default = "cell_type",
    help    = "Metadata column for cell types [default: %default]"
  ),
  optparse::make_option(
    "--n-splits-within",
    type    = "integer",
    default = 3L,
    help    = "Number of reference splits [default: %default]"
  ),
  optparse::make_option(
    "--cutoff",
    type    = "double",
    default = 0.1,
    help    = "inferCNV cutoff [default: %default]"
  ),
  optparse::make_option(
    "--chr-exclude",
    type    = "character",
    default = "MT,Y",
    help    = "Comma-separated chromosomes to exclude [default: %default]"
  ),
  optparse::make_option(
    "--resume-if-exists",
    action  = "store_true",
    default = TRUE,
    help    = "Skip completed runs [default: %default]"
  ),
  
  # ---- Block 2 ------------------------------------------------------------
  optparse::make_option(
    "--base-dir",
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
    "--max-gap",
    type    = "integer",
    default = 100000L,
    help    = "Maximum gap for merging nearby segments [default: %default]"
  ),
  optparse::make_option(
    "--min-overlap-consistent-calls",
    type    = "double",
    default = 0.75,
    help    = "Minimum overlap for equivalence assignment [default: %default]"
  ),
  optparse::make_option(
    "--min-overlap-multiple-nodes",
    type    = "double",
    default = 0.6,
    help    = "Minimum overlap for multi-node merging [default: %default]"
  ),
  optparse::make_option(
    "--filter-seq-mb-init",
    type    = "double",
    default = 5.0,
    help    = "Minimum segment length before merging in Mb [default: %default]"
  ),
    optparse::make_option(
    "--overlap-method",
    type    = "character",
    default = "jaccard",
    help    = "Type of overlap method implemented"
  ),
  optparse::make_option(
    "--filter-seq-mb-equiv",
    type    = "double",
    default = 7.0,
    help    = "Minimum segment length before equivalence in Mb [default: %default]"
  ),
  optparse::make_option(
    "--min-references",
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
    "--group-cols",
    type    = "character",
    default = "cell_type",
    help    = "Comma-separated grouping columns [default: %default]"
  ),
  optparse::make_option(
    "--chromosome-arms-path",
    type    = "character",
    default = NULL,
    help    = "Path to chromosome arms RDS. If NULL uses hg38 built-in [default: %default]"
  ),
  
  # ---- Block 4 ------------------------------------------------------------
    optparse::make_option(
    "--by-col",
    type    = "character",
    default = "cell_type",
    help    = "How do you want to cluster your segments"
  ),
  optparse::make_option(
  "--cnv-annotated",
  type    = "character",
  default = NULL,
  help    = "Mandatory Variable if you want to start the pipeline in later stages [default: %default]"
),
  optparse::make_option(
  "--cell-sizes",
  type    = "character",
  default = NULL,
  help    = "Number of cells per cell type within a data frame[default: %default]"
),
  
  optparse::make_option(
  "--k-value",
  type    = "double",
  default = NULL,
  help    = paste(
    "Scaling constant for group-size threshold:",
    "group_threshold = k * sqrt(n_total_cells).",
    "Must be derived from a trusted calibration point",
    "(k = known_good_threshold / sqrt(known_good_n_cells)).",
    "No default — required."
  )
    ),
  optparse::make_option(
  "--sensitivity-floor-mb",
  type    = "double",
  default = 20,
  help    = paste(
    "Event size (Mb) at or below which events are filtered",
    "out entirely (below empirical detection limit)",
    "[default: %default]"
    )
  ),
  optparse::make_option(
    "--cell-id-col",
    type    = "character",
    default = "cell_name",
    help    = "Column which has the Unique ID for each cell"
  ),
  optparse::make_option("--sample-col",
    type = "character",
    default = "cell_type",
    help = "Column that we want to use to separate the cells, e.g (embryo or cell type)"
    ),
  optparse::make_option(
    "--min-required-cells",
    type    = "integer",
    default = 5L,
    help    = "Minimum cell count floor [default: %default]"
  ),
  optparse::make_option(
  "--p-arm-permission",
  type    = "double",
  default = 70,
  help    = "Minimum percentage of p arm covered for p_centromere events [default: %default]"
),
optparse::make_option(
  "--q-arm-permission",
  type    = "double",
  default = 70,
  help    = "Minimum percentage of q arm covered for centromere_q events [default: %default]"
),
optparse::make_option(
  "--whole-chr-permission",
  type    = "double",
  default = 60,
  help    = "Minimum percentage of both arms covered for p_centromere_q events [default: %default]"
),
optparse::make_option(
  "--k-interval",
  type    = "double",
  default = 60,
  help    = " Turning continous values of gain and loss into discrete ones"
),
  optparse::make_option(
    "--min-overlap",
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
      "Runs blocks starting from --start-from.",
      "Use --execution-mode to specify single or array job context."
    )
  )
)

# =============================================================================
# Validate execution_mode
# =============================================================================


if (!opt$`execution-mode` %in% c("single", "array")) {
  stop(
    "--execution-mode must be 'single' or 'array'.\n",
    "Got: '", opt$`execution-mode`, "'"
  )
}

# =============================================================================
# Validate required arguments per start_from + execution_mode
# =============================================================================

if (is.null(opt$workdir)) {
  stop("--workdir is required")
}

if (opt$`start-from` == "block1") {
  if (opt$`execution-mode` == "array") {
    stop(
      "--start-from block1 not supported with --execution-mode array.\n",
      "For array mode Block1 use: Rscript scripts/run_block1_array.R"
    )
  }
  if (is.null(opt$`counts-path`))      stop("--counts-path required for block1")
  if (is.null(opt$`metadata-path`))    stop("--metadata-path required for block1")
  if (is.null(opt$`gene-order-file`))  stop("--gene-order-file required for block1")
  if (is.null(opt$`tool-outdir`))      stop("--tool-outdir required for block1")
}


if (is.null(opt$`k-value`)) {
  stop(
    "--k-value is required.\n",
    "Derive it from a trusted calibration point:\n",
    "  k = known_good_threshold / sqrt(known_good_n_cells)"
  )
}


if (opt$`start-from` == "block2" && is.null(opt$`tool-outdir`)) {
  stop("--base-dir required when --start-from = block2")
}


# Parse by_col as vector
# same pattern as group_cols
by_col <- if (!is.null(opt$`by-col`) &&
               nzchar(opt$`by-col`)) {
  trimws(strsplit(opt$`by-col`, ",")[[1]])
} else {
  NULL
}


# =============================================================================
# Load dependencies
# =============================================================================

message("Loading pipeline functions...")
source(file.path(project_dir, "R", "cnv_annotation.R"))
source(file.path(project_dir, "R", "cnv_processing.R"))
source(file.path(project_dir, "R", "cnv_scoring.R"))
source(file.path(project_dir, "R", "infercnv.R"))
source(file.path(project_dir, "R", "pipeline.R"))

# =============================================================================
# Parse vector arguments
# Note: optparse converts hyphens to underscores in opt$ names
# e.g. --group-cols → opt$`group_cols`
# =============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

chr_exclude   <- as.character(trimws(strsplit(opt$`chr-exclude` %||% "MT,Y", ",")[[1]]))
group_cols    <- trimws(strsplit(opt$`group-cols`, ",")[[1]])


# =============================================================================
# Type conversion and validation
# SLURM passes everything as strings
# convert to correct types before pipeline
# =============================================================================

# ── Integers ─────────────────────────────────────────────────────────────────
opt$`n-splits-within`  <- as.integer(opt$`n-splits-within`)
opt$`max-gap`          <- as.integer(opt$`max-gap`)
opt$`min-references`   <- as.integer(opt$`min-references`)
opt$`min-required-cells` <- as.integer(opt$`min-required-cells`)
opt$cores              <- as.integer(opt$cores)

# ── Doubles ───────────────────────────────────────────────────────────────────
opt$cutoff                          <- as.double(opt$cutoff)
opt$`min-overlap`                   <- as.double(opt$`min-overlap`)
opt$`min-overlap-consistent-calls`  <- as.double(opt$`min-overlap-consistent-calls`)
opt$`min-overlap-multiple-nodes`    <- as.double(opt$`min-overlap-multiple-nodes`)
opt$`filter-seq-mb-init`            <- as.double(opt$`filter-seq-mb-init`)
opt$`filter-seq-mb-equiv`           <- as.double(opt$`filter-seq-mb-equiv`)
opt$`p-arm-permission`              <- as.double(opt$`p-arm-permission`)
opt$`q-arm-permission`              <- as.double(opt$`q-arm-permission`)
opt$`whole-chr-permission`          <- as.double(opt$`whole-chr-permission`)
k_interval                          <- as.double(opt$`k-interval`)

opt$`min-coding-density` <- as.double(opt$`min-coding-density`)
opt$`max-gap-mb` <- as.double(opt$`max-gap-mb`)
  
opt$`k-value`               <- as.double(opt$`k-value`)
opt$`sensitivity-floor-mb`  <- as.double(opt$`sensitivity-floor-mb`)

# ── Logical ───────────────────────────────────────────────────────────────────
# isTRUE handles "TRUE"/"FALSE" strings from SLURM
opt$`remove-reference`  <- isTRUE(as.logical(opt$`remove-reference`))
opt$`save-intermediate` <- isTRUE(as.logical(opt$`save-intermediate`))
opt$`resume-if-exists`  <- isTRUE(as.logical(opt$`resume-if-exists`))
opt$parallel            <- isTRUE(as.logical(opt$parallel))

# ── Validate no NAs after conversion ─────────────────────────────────────────
check_no_na <- function(val, name) {
  if (any(is.na(val))) {
    stop("Failed to parse '", name,
         "' — got NA after type conversion.\n",
         "Check the value passed from SLURM.")
  }
}

check_no_na(opt$`n-splits-within`,  "--n-splits-within")
check_no_na(opt$`max-gap`,          "--max-gap")
check_no_na(opt$`min-references`,   "--min-references")
check_no_na(opt$`min-required-cells`,"--min-required-cells")
check_no_na(opt$cutoff,             "--cutoff")
check_no_na(opt$`min-overlap`,      "--min-overlap")
check_no_na(k_interval,             "--k-interval")
check_no_na(opt$`k-value`,              "--k-value")
check_no_na(opt$`sensitivity-floor-mb`, "--sensitivity-floor-mb")
check_no_na(opt$`min-coding-density`, "--min-coding-density")
check_no_na(opt$`max-gap-mb`, "--max-gap-mb")


cat("Type conversion complete\n")



# =============================================================================
# Load inputs
# =============================================================================

# Metadata — required for block1, optional for block2 array mode
metadata <- NULL
if (!is.null(opt$`metadata-path`)) {
  if (!file.exists(opt$`metadata-path`)) {
    stop("metadata-path not found: ", opt$`metadata-path`)
  }
  message("Loading metadata from: ", opt$`metadata-path`)
  metadata <- readRDS(opt$`metadata-path`)
}

# Counts matrix — block1 only
counts_mx <- NULL
if (!is.null(opt$`counts-path`)) {
  if (!file.exists(opt$`counts-path`)) {
    stop("counts-path not found: ", opt$`counts-path`)
  }
  message("Loading counts matrix from: ", opt$`counts-path`)
  counts_mx <- readRDS(opt$`counts-path`)
}

if (opt$`start-from` %in% c("block3", "block4")) {
  
  cnv_annotated <- NULL
  if (!is.null(opt$`cnv-annotated`)) {
    if (!file.exists(opt$`cnv-annotated`)) {
      stop("cnv-annotated not found: ",
           opt$`cnv-annotated`)
    }
    cnv_annotated <- readRDS(opt$`cnv-annotated`)
    message("Loaded cnv_annotated: ",
            nrow(cnv_annotated), " rows")
  }
  
  cell_sizes <- NULL
  if (!is.null(opt$`cell-sizes`)) {
    if (!file.exists(opt$`cell-sizes`)) {
      stop("cell-sizes not found: ",
           opt$`cell-sizes`)
    }
    cell_sizes <- readRDS(opt$`cell-sizes`)
    message("Loaded cell_sizes: ",
            nrow(cell_sizes), " rows")
  }
  
} else {
  # block1 or block2 — not needed
  cnv_annotated <- NULL
  cell_sizes    <- NULL
}




# Chromosome arms
chromosome_arms <- NULL

if (!is.null(opt$`chromosome-arms-path`)) {
  chromosome_arms <- readRDS(opt$`chromosome-arms-path`)
} else {
  chrom_path <- file.path(project_dir, "data", "hg38_chromosome_arms.rds")
  if (!file.exists(chrom_path)) {
    stop(
      "No --chromosome-arms-path provided and built-in not found at:\n",
      chrom_path, "\n",
      "Either provide --chromosome-arms-path or save hg38_chromosome_arms.rds to data/"
    )
  }
  message("Loading built-in hg38 chromosome arms from: ", chrom_path)
  chromosome_arms <- readRDS(chrom_path)
}


# ── Build coding_gr and coding_expressed_set if needed ────────────────────────
coding_gr            <- NULL
coding_expressed_set <- NULL


  
  if (is.null(opt$`coding-genes-path`)) {
    stop("--coding-genes-path required")
  }
  if (!file.exists(opt$`coding-genes-path`)) {
    stop("coding-genes-path not found: ", opt$`coding-genes-path`)
  }
  
  cat("Loading coding genes from:",
      opt$`coding-genes-path`, "\n")
  coding_genes <- readRDS(opt$`coding-genes-path`)
  
  coding_gr <- GenomicRanges::GRanges(
    seqnames = coding_genes$chr,
    ranges   = IRanges::IRanges(
      start = coding_genes$start,
      end   = coding_genes$end
    )
  )
  coding_gr$gene <- coding_genes$gene_name
  
  coding_expressed_set <- unique(coding_genes$gene_name)
  
  print(coding_gr)
  cat("Coding genes loaded:", length(coding_expressed_set), "\n")



# =============================================================================
# Create output directory
# =============================================================================

if (!dir.exists(opt$workdir)) {
  message("Creating output directory: ", opt$workdir)
  dir.create(opt$workdir, recursive = TRUE)
}


# =============================================================================
# Print input summary
# =============================================================================

message(paste0(
  "\n=== Submission parameters ===\n",
  "  EXECUTION_MODE:          ", opt$`execution-mode`,            "\n",
  "  START_FROM:              ", opt$`start-from`,                "\n",
  "  WORKDIR:                 ", opt$workdir,                     "\n",
  "\n--- Paths ---\n",
  "  METADATA_PATH:           ", opt$`metadata-path`        %||% "NULL", "\n",
  "  COUNTS_PATH:             ", opt$`counts-path`          %||% "NULL", "\n",
  "  TOOL_OUTDIR:             ", opt$`tool-outdir`          %||% "NULL", "\n",
  "  CHROMOSOME_ARMS:         ", opt$`chromosome-arms-path` %||% "built-in hg38", "\n",
  "  GENE_ORDER:              ", opt$`gene-order-file`      %||% "NULL", "\n",
  "  TOOL:                    ", opt$tool,                        "\n",
  "  N_SPLITS:                ", opt$`n-splits-within`,           "\n",
  "  CUTOFF:                  ", opt$cutoff,                      "\n",
  "  REMOVE_REFERENCE:        ", opt$`remove-reference`,          "\n",
  "  PATTERN:                 ", opt$pattern,                     "\n",
  "  MAX_GAP:                 ", opt$`max-gap`,                   "\n",
  "  MIN_OVERLAP_CONSISTENT:  ", opt$`min-overlap-consistent-calls`, "\n",
  "  MIN_OVERLAP_NODES:       ", opt$`min-overlap-multiple-nodes`,"\n",
  "  FILTER_SEQ_MB_INIT:      ", opt$`filter-seq-mb-init`,        "\n",
  "  FILTER_SEQ_MB_EQUIV:     ", opt$`filter-seq-mb-equiv`,       "\n",
  "  MIN_REFERENCES:          ", opt$`min-references`,            "\n",
  "  GROUP_COLS:              ", paste(group_cols, collapse = ", "), "\n",
  "  K_VALUE_THRESHOLD:                 ", opt$`k-value`, "\n",
  "  SENSITIVITY_FLOOR_MB:    ", opt$`sensitivity-floor-mb`, "\n",
  "  MIN_REQUIRED_CELLS:      ", opt$`min-required-cells`, "\n",
  "  BY_COL:                  ", paste(by_col %||% "NULL", collapse = ", "), "\n",
  "  SAMPLE_COL:              ", opt$`sample-col`  %||% "NULL",  "\n",
  "  CELL_COL:                ", opt$`cell-id-col`,               "\n",
  "  P_ARM_PERMISSION:        ", opt$`p-arm-permission`,   "%\n",
  "  Q_ARM_PERMISSION:        ", opt$`q-arm-permission`,   "%\n",
  "  WHOLE_CHR_PERMISSION:    ", opt$`whole-chr-permission`, "%\n",
  "  MIN_OVERLAP:             ", opt$`min-overlap`,               "\n",
   " OVERLAP METHOD:          ", opt$`overlap-method`,           "\n",
  "============================="
))

# =============================================================================
# Run pipeline
# =============================================================================



print(opt$`clonal-col`)
print(opt$`donor-col`)

results <- run_full_cnv_pipeline(
  
  # ---- General ─────────────────────────────────────────────────────────────
  start_from        = opt$`start-from`,
  execution_mode    = opt$`execution-mode`,
  save_intermediate = isTRUE(opt$`save-intermediate`),
  workdir           = opt$workdir,
  
  # ---- Block 1 ─────────────────────────────────────────────────────────────
  counts_mx              = counts_mx,
  metadata               = metadata,
  group_clusters_col     = opt$`group-clusters-col`,
  gene_order_file        = opt$`gene-order-file`,
  chromosomes_to_exclude = chr_exclude,
  n_splits_within        = opt$`n-splits-within`,
  tool_outdir            = opt$`tool-outdir`,
  tool                   = opt$tool,
  cutoff                 = opt$cutoff,
  remove_ref             = isTRUE(opt$`remove-reference`),  
  resume_if_exists       = isTRUE(opt$`resume-if-exists`),
  clonal_col = opt$`clonal-col`,
  donor_col = opt$`donor-col`,
  
  # ---- Block 2 ─────────────────────────────────────────────────────────────
  pattern                      = opt$pattern,
  max_gap                      = opt$`max-gap`,
  k_discrete                   = k_interval,
  min_overlap_consistent_calls = opt$`min-overlap-consistent-calls`,
  min_overlap_multiple_nodes   = opt$`min-overlap-multiple-nodes`,
  filter_seq_mb_init           = opt$`filter-seq-mb-init`,
  filter_seq_mb_equiv          = opt$`filter-seq-mb-equiv`,
  min_references               = opt$`min-references`,
  parallel                     = isTRUE(opt$parallel),
  cores                        = opt$cores,
  pct_max                   = opt$`pct-max`,
  pct_floor                 = opt$`pct-floor`,
  min_expr_density          = opt$`min-expr-density`,
  min_coding_density        = opt$`min-coding-density`,
  max_gap_mb                = opt$`max-gap-mb`,
  coding_gr                 = coding_gr,
  coding_expressed_set      = coding_expressed_set,
  
  # ---- Block 3 ─────────────────────────────────────────────────────────────
  chromosome_arms = chromosome_arms,
  group_cols      = group_cols,
  
  # ---- Block 4 ─────────────────────────────────────────────────────────────
  cell_sizes           = cell_sizes,
  cnv_annotated        = cnv_annotated,
  by                   = by_col,
  sample_col         =    opt$`sample-col`,
  cell_col             = opt$`cell-id-col`,
  min_required_cells   = opt$`min-required-cells`,
  p_arm_permission     = opt$`p-arm-permission`,
  q_arm_permission     = opt$`q-arm-permission`,
  whole_chr_permission = opt$`whole-chr-permission`,
  min_overlap          = opt$`min-overlap`,
  sensitivity_floor_mb = opt$`sensitivity-floor-mb`,
  k_threshold_growth =  opt$`k-value`,
  overlap_method       = opt$`overlap-method`)




# =============================================================================
# Save outputs
# =============================================================================

out_path <- file.path(opt$workdir, "pipeline_results.rds")
message("Saving full pipeline results to: ", out_path)
saveRDS(results, out_path)

message(sprintf("\nDone. Results saved to: %s", out_path))