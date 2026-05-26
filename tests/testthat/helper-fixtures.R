

source("~/GitHub/cnvscfinder/R/cnv_annotation.R")
source("~/GitHub/cnvscfinder/R/cnv_processing.R")
source("~/GitHub/cnvscfinder/R/cnv_scoring.R")
source("~/GitHub/cnvscfinder/R/gsva_gsea.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/infercnv.R")
source("C:/Users/pmgra/Documents/GitHub/cnvscfinder/R/pipeline.R")


hg38_chromosome_arms <- readRDS("C:/Users/pmgra/Documents/GitHub/cnvscfinder/tests/testthat/fixtures/hg38_chromosome_arms.rds")



make_mock_cnv <- function() {
  data.frame(
    cell_name     = c("E5.10.933", "E5.10.933", "E5.10.933", "E5.10.934"),
    chr           = c("chr1",      "chrX",      "chrX",       "chr3"),
    cnv_state     = c("loss", "gain",      "gain",       "gain"),
    start         = c(117863485L,  135344796L,  14008279L,    169772831L),
    end           = c(150476566L,  149000663L,  75746911L,    184368596L),
    cnv_length    = c(32613082L, 13655868L,   61738633L,    14595766L),
    cnv_length_mb = c(32.6, 13.7,        61.7,         14.6),
    n_references  = c(2L, 2L,         2L,            2L),
    references    = c("A,B", "A,B",       "A,B",         "A,B"),
    mode          = c("within", "within",   "within",      "within"),
    cell_type     = c("TE",  "TE",      "TE",          "TE"),
    embryo        = c("E5.10",  "E5.11",   "E5.10",       "E5.10"),
    stringsAsFactors = FALSE
  )
}




make_mock_segments <- function() {
  data.frame(
    reference = c("A", "A", "A", "A", "A"),
    cell_name = c("E5.10.933", "E5.10.933", "E5.10.933",
                  "E5.10.934", "E5.10.934"),
    chr       = c("chr1",  "chr1",      "chr13",
                  "chr3",  "chr3"),
    start     = c(117863485L, 201359008L, 18816718L,
                  50350695L,  169772831L),
    stop      = c(150476566L, 202809470L, 75549439L,
                  51983409L,  183684477L),
    state     = c("loss", "gain", "loss",
                  "loss", "gain"),
    stringsAsFactors = FALSE
  )
}



make_mock_metadata <- function() {
  data.frame(
    cell_name   = c("E5.5.101", "E5.5.100", "E6.2.114",
                    "E6.2.104", "E6.2.107", "E7.2.166"),
    cell_type   = c("TE",       "TE",       "Epiblast",
                    "Epiblast", "Epiblast", "TE"),
    split_group = c("A",       "B",        "A",
                    "B",       "C",        "A"),
    embryo      = c("E5.5",     "E5.5",     "E6.2",
                    "E6.2",     "E6.2",     "E7.2"),
    stringsAsFactors = FALSE
  )
}



make_mock_pipeline_metadata <- function() {
  data.frame(
    cell_name  = c("E5.5.101", "E5.5.100", "E6.2.114",
                   "E6.2.104", "E6.2.107", "E7.2.166"),
    cell_type  = c("TE",       "TE",       "Epiblast",
                   "Epiblast", "Epiblast", "TE"),
    embryo     = c("E5.5",     "E5.5",     "E6.2",
                   "E6.2",     "E6.2",     "E7.2"),
    split_group = c("A",       "B",        "A",
                    "B",       "C",        "A"),
    stringsAsFactors = FALSE
  )
}



make_mock_gene_level <- function() {
  data.frame(
    gene = c(
      "GENE1", "GENE1", "GENE1",
      "GENE2", "GENE2", "GENE2",
      "GENE3", "GENE3", "GENE3",
      "GENE1", "GENE1", "GENE1",
      "GENE2", "GENE2", "GENE2",
      "GENE3", "GENE3", "GENE3",
      "GENE1", "GENE1", "GENE1", "GENE1",
      "GENE4", "GENE5", "GENE6", "GENE7",
      "GENE8", "GENE9",
      "GENE10", "GENE11"
    ),
    
    cell_name = c(
      rep("E5.5.101", 9),
      rep("E5.5.100", 9),
      rep("E6.2.114", 4),
      rep("E6.2.107", 8)
    ),
    
    state_raw = c(
      rep(1.50, 9),   
      rep(0.60, 9), 
      rep(1.01, 4),  
      c(1.50, 1.01, 1.50, 1.50,   
        1.50, 1.50, 0.60, 0.60)  
    ),
    
    chr = c(
      rep("chr15", 9), 
      rep("chr15", 9),   
      rep("chr15", 4),  
      c("chr15", "chr15", "chr15", "chr15",  
        "chr7",  "chr7",                     
        "chr22", "chr22")                   
    ),
    
    start = c(
      1000000L, 1000000L, 1000000L,  
      2000000L, 2000000L, 2000000L,  
      3000000L, 3000000L, 3000000L,
      1000000L, 1000000L, 1000000L,
      2000000L, 2000000L, 2000000L,
      3000000L, 3000000L, 3000000L,
      1000000L, 2000000L, 3000000L, 3450000L,
      1000000L, 2000000L, 3000000L, 3500000L,
      5000000L, 6000000L,
      8000000L, 9000000L
    ),
    
    stop = c(
      1500000L, 1500000L, 1500000L,
      2500000L, 2500000L, 2500000L,
      3500000L, 3500000L, 3500000L,
      1500000L, 1500000L, 1500000L,
      2500000L, 2500000L, 2500000L,
      3500000L, 3500000L, 3500000L,
      1500000L, 2500000L, 3500000L, 4000000L,
      1500000L, 2500000L, 3500000L, 4000000L,
      5500000L, 6500000L,
      8500000L, 9500000L
    ),
    
    state = c(
      rep("gain", 9),    
      rep("loss", 9),     
      rep("neutral", 4),  
      c("gain", "neutral", "gain", "gain",   
        "gain", "gain",                     
        "loss", "loss")                      
    ),
    
    reference = c(
      "A", "B", "C",
      "A", "B", "C",
      "A", "B", "C",
      "A", "B", "C",
      "A", "B", "C",
      "A", "B", "C",

      "A", "B", "C", "C",
      rep("A", 8)
    ),
    
    stringsAsFactors = FALSE
  )
}


make_mock_segments <- function() {
  data.frame(
    reference = c("A",         "A",         "B",         "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.934", "E5.10.933", "E5.10.934", "E5.10.933"),
    chr       = c("chr1",      "chr3",      "chr1",      "chr3",      "chr1"),
    start     = c(117863485L,  169772831L,  117863485L,  169772831L,  117863485L),
    stop      = c(150476566L,  183684477L,  150476566L,  183684477L,  150476566L),
    state     = c("loss",      "gain",      "loss",      "gain",      "loss"),
    stringsAsFactors = FALSE
  )
}


make_mock_merged <- function() {
  data.frame(
    reference = c("A",         "A",         "B",         "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.934", "E5.10.933", "E5.10.934", "E5.10.933"),
    chr       = c("chr1",      "chr3",      "chr1",      "chr3",      "chr1"),
    cnv_state = c("loss",      "gain",      "loss",      "gain",      "loss"),
    start     = c(117863485L,  169772831L,  117863485L,  169772831L,  117863485L),
    end       = c(150476566L,  183684477L,  150476566L,  183684477L,  150476566L),
    n_segments = c(1L, 1L, 1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}


make_scenario_perfect_overlap <- function() {
  data.frame(
    reference = c("A",         "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.933", "E5.10.933"),
    chr       = c("chr1",      "chr1",      "chr1"),
    cnv_state = c("loss",      "loss",      "loss"),
    start     = c(117863485L,  117863485L,  117863485L),
    end       = c(150476566L,  150476566L,  150476566L),
    stringsAsFactors = FALSE
  )
}

make_scenario_ab_only <- function() {
  data.frame(
    reference = c("A",         "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.933", "E5.10.933"),
    chr       = c("chr1",      "chr1",      "chr1"),
    cnv_state = c("loss",      "loss",      "loss"),
    start     = c(117863485L,  117863485L,  200000000L),  # C is far away
    end       = c(150476566L,  150476566L,  210000000L),
    stringsAsFactors = FALSE
  )
}

make_scenario_ac_ab_no_bc <- function() {
  data.frame(
    reference = c("A",         "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.933", "E5.10.933"),
    chr       = c("chr15",     "chr15",     "chr15"),
    cnv_state = c("loss",      "loss",      "loss"),
    start     = c(10000000L,   10000000L,   55000000L),
    end       = c(100000000L,  65000000L,   100000000L),
    stringsAsFactors = FALSE
  )
}



make_scenario_no_overlap <- function() {
  data.frame(
    reference = c("A",        "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.933", "E5.10.933"),
    chr       = c("chr1",      "chr1",      "chr1"),
    cnv_state = c("loss",      "loss",      "loss"),
    start     = c(10000000L,   100000000L,  200000000L),  # all far apart
    end       = c(20000000L,   110000000L,  210000000L),
    stringsAsFactors = FALSE
  )
}


make_scenario_ab_pass_ac_fail <- function() {
  data.frame(
    reference = c("A",         "B",         "C"),
    cell_name = c("E5.10.933", "E5.10.933", "E5.10.933"),
    chr       = c("chr1",      "chr1",      "chr1"),
    cnv_state = c("loss",      "loss",      "loss"),
    # A and B: near perfect overlap → passes 0.75
    # A and C: very small overlap → fails 0.75
    start     = c(100000000L,  100000000L,  149000000L),
    end       = c(150000000L,  148000000L,  160000000L),
    stringsAsFactors = FALSE
  )
}


make_mock_supported_events <- function() {
  data.frame(
    cell_name     = c("E5.10.933", "E5.10.933", "E5.10.933", "E5.10.934"),
    chr           = c("chr1",      "chrX",      "chrX",       "chr3"),
    cnv_state     = c("loss",      "gain",      "loss",       "gain"),
    cnv_equiv_id  = c("E5.10.933|chr1|loss|1",
                      "E5.10.933|chrX|gain|1",
                      "E5.10.933|chrY|loss|1",
                      "E5.10.934|chr3|gain|1"),
    start         = c(117863485L,  14008279L,   2841486L,    169772831L),
    end           = c(150476566L,  75746911L,   20781032L,   184368596L),
    cnv_length    = c(32613082L,   61738633L,   17939547L,   14595766L),
    cnv_length_mb = c(32.6,        61.7,        17.9,        14.6),
    n_references  = c(2L,          2L,           2L,           2L),
    references    = c("A,B",       "A,B",        "A,B",        "A,B"),
    merge_group_id = c("E5.10.933|chr1|loss|1",
                       "E5.10.933|chrX|gain|1",
                       "E5.10.933|chrX|loss|1",
                       "E5.10.934|chr3|gain|1"),
    stringsAsFactors = FALSE
  )
}

make_scenario_single_ref <- function() {
  data.frame(
    reference = c("A"),
    cell_name = c("E5.10.933"),
    chr       = c("chr15"),
    cnv_state = c("loss"),
    start     = c(117863485L),
    end       = c(150476566L),
    stringsAsFactors = FALSE
  )
}


make_mock_clustered_events <- function() {
  data.frame(
    cell_name             = c("E5.10.945", "E5.10.949", "E5.10.945",
                              "E5.10.949", "E5.10.934", "E5.10.936"),
    chr                   = c("chr1",  "chr1",  "chr19",
                              "chr19", "chrX",  "chrX"),
    cnv_state             = c("loss",  "loss",  "gain",
                              "gain",  "gain",  "gain"),
    start                 = c(231223763L, 231332753L, 301444L,
                              301444L,    117897813L, 13035618L),
    end                   = c(248919946L, 248919946L, 8577577L,
                              10115372L,  149000663L, 99721158L),
    cnv_length            = c(17696184L,  17587194L,  8276134L,
                              9813929L,   31102851L,  86685541L),
    cnv_length_mb         = c(17.7,       17.6,       8.3,
                              9.8,        31.1,       86.7),
    n_references          = c(2L, 2L, 2L, 2L, 2L, 2L),
    references            = c("A,B", "A,B", "A,B", "A,B", "A,B", "A,B"),
    mode                  = c("within", "within", "within",
                              "within", "within", "within"),
    cell_type             = c("TE", "TE", "TE", "TE", "TE", "TE"),
    embryo                = c("E5.10", "E5.10", "E5.10",
                              "E5.10", "E5.10", "E5.10"),
    split_group           = c("C", "C", "C", "C", "C", "C"),
    arm_class             = c("q_arm",       "q_arm",       "p_arm",
                              "p_arm",       "q_arm",       "p_centromere_q"),
    whole_chromosome_gain = c(NA,    NA,    14.12, 16.74, 19.93, 55.55),
    whole_chromosome_loss = c(7.11,  7.06,  NA,    NA,    NA,    NA),
    p_arm_gain            = c(NA,    NA,    33.78, 40.06, NA,    NA),
    p_arm_loss            = c(NA,    NA,    NA,    NA,    NA,    NA),
    q_arm_gain            = c(NA,    NA,    NA,    NA,    33.22, NA),
    q_arm_loss            = c(14.27, 14.18, NA,    NA,    NA,    NA),
    cnv_equiv_id          = c("E5.10|chr1|loss|1",  "E5.10|chr1|loss|1",
                              "E5.10|chr19|gain|1", "E5.10|chr19|gain|1",
                              "E5.10|chrX|gain|1",  "E5.10|chrX|gain|2"),
    stringsAsFactors = FALSE
  )
}

make_mock_locus_summary <- function() {
  data.frame(
    embryo        = c("E5.10",             "E5.10",              "E5.10",             "E5.10"),
    cnv_equiv_id  = c("E5.10|chr1|loss|1", "E5.10|chr19|gain|1", "E5.10|chrX|gain|1", "E5.10|chrX|gain|2"),
    chr           = c("chr1",  "chr19", "chrX",  "chrX"),
    cnv_state     = c("loss",  "gain",  "gain",  "gain"),
    locus_start   = c(231223763L, 301444L,    117897813L, 13035618L),
    locus_end     = c(248919946L, 10115372L,  149000663L, 99721158L),
    locus_width   = c(17696184L,  9813929L,   31102851L,  86685541L),
    locus_width_mb = c(17.7,      9.8,        31.1,       86.7),
    n_events      = c(2L, 2L, 1L, 1L),
    n_cells       = c(2L, 2L, 1L, 1L),
    n_samples     = c(1L, 1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}


make_mock_cell_sizes <- function() {
  data.frame(
    embryo        = c("E5.10", "E5.11", "E6.2"),
    n_total_cells = c(20L,     15L,     25L),
    stringsAsFactors = FALSE
  )
}


make_mock_thresholded <- function() {
  df <- make_mock_clustered_events()
  df$n_total_cells       <- 20L
  df$tier                <- c("tier_2", "tier_2", "below_threshold",
                              "below_threshold", "tier_1", "tier_1")
  df$effective_threshold <- c(5L, 5L, NA, NA, 3L, 3L)
  df
}



make_mock_prepare_input <- function() {
  data.frame(
    cnv_equiv_id  = c("E5.10|chr1|loss|1",  "E5.10|chr19|gain|1",
                      "E5.10|chrX|gain|1",  "E5.10|chrX|gain|2"),
    embryo        = c("E5.10", "E5.10", "E5.10", "E5.10"),
    cnv_length_mb = c(17.7,    9.8,     31.1,    86.7),
    n_cells       = c(2L,      2L,      1L,      1L),
    n_total_cells = c(20L,     20L,     20L,     20L),
    stringsAsFactors = FALSE
  )
}


make_mock_cell_sizes <- function() {
  data.frame(
    embryo        = c("E5.10", "E5.11", "E6.2"),
    n_total_cells = c(20L,     15L,     25L),
    stringsAsFactors = FALSE
  )
}




# ============================================================================
# inferCNV fixtures
# ============================================================================

# Mock expression matrix — genes x cells (sparse-like, use regular matrix)
make_mock_expr_matrix <- function() {
  genes <- c("DDX11L1", "WASH7P", "MIR6723", "GENE4", "GENE5")
  cells <- c("E5.5.101", "E5.5.100", "E6.2.114", "E6.2.104")
  
  mat <- matrix(
    data     = c(0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0, 0.0,
                 0.0, 3.0, 3.0, 2.0,
                 1.5, 1.5, 1.5, 1.5,
                 0.8, 0.8, 1.2, 1.2),
    nrow     = 5,
    ncol     = 4,
    dimnames = list(genes, cells)
  )
  mat
}

# Mock gene order — matches real structure: rownames=gene, cols=chr/start/stop
make_mock_gene_order <- function() {
  df <- data.frame(
    chr   = c("chr1",   "chr1",   "chr1",   "chr15",  "chr15"),
    start = c(11869L,   14404L,   632325L,  1000000L, 2000000L),
    stop  = c(14409L,   29570L,   632413L,  1500000L, 2500000L),
    row.names = c("DDX11L1", "WASH7P", "MIR6723", "GENE4", "GENE5"),
    stringsAsFactors = FALSE
  )
  df
}

# Mock long format — output of melt_expr_to_long
make_mock_long_df <- function() {
  data.frame(
    gene      = c("MIR6723", "MIR6723", "MIR6723", "MIR6723",
                  "GENE4",   "GENE4",   "GENE4",   "GENE4"),
    cell_name = c("E5.5.101", "E5.5.100", "E6.2.114", "E6.2.104",
                  "E5.5.101", "E5.5.100", "E6.2.114", "E6.2.104"),
    state_raw = c(0.0,  3.0,  3.0,  2.0,
                  1.5,  1.5,  1.5,  1.5),
    stringsAsFactors = FALSE
  )
}

# Mock infercnv_list — input to load_and_prepare_infercnv_reference
# Each element is list(expr_matrix, gene_order)
make_mock_infercnv_list <- function() {
  list(
    A = list(make_mock_expr_matrix(), make_mock_gene_order()),
    B = list(make_mock_expr_matrix(), make_mock_gene_order()),
    C = list(make_mock_expr_matrix(), make_mock_gene_order())
  )
}

# Mock metadata for validate_metadata and make_splits
make_mock_infercnv_metadata <- function() {
  data.frame(
    cell_name = c("E5.5.101", "E5.5.100", "E6.2.114", "E6.2.104"),
    cell_type = c("TE",       "TE",       "Epiblast", "Epiblast"),
    stringsAsFactors = FALSE
  )
}











