#!/bin/bash
# =============================================================================
# submit_pipeline.sh
# Submit CNV pipeline to SLURM
#
# Usage:
#   Single job:  sbatch submit_pipeline.sh
#   Array job:   sbatch --array=1-N submit_pipeline.sh
# =============================================================================

#SBATCH --job-name=cnv_pipeline
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/cnv_pipeline_%j_%a.out
#SBATCH --error=logs/cnv_pipeline_%j_%a.err

# ---- Load modules -----------------------------------------------------------
module load R/4.3.0

# ---- Paths ------------------------------------------------------------------
BASE_DIR="/scratch/brussel/vo/000/bvo00016/vsc11567/infercnv_results/d7"
METADATA="/scratch/brussel/vo/000/bvo00016/vsc11567/data/metadata.rds"
OUTDIR="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_results/d7"
SCRIPT="scripts/run_pipeline.R"

# ---- Array job support ------------------------------------------------------
# If running as array job, read dataset from sample sheet
# Sample sheet format: one dataset path per line
# e.g. samples.txt:
#   /scratch/.../d7
#   /scratch/.../d10
#   /scratch/.../d14

if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
    SAMPLE_SHEET="config/samples.txt"
    BASE_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_SHEET}")
    DATASET=$(basename "${BASE_DIR}")
    OUTDIR="${OUTDIR}/${DATASET}"
    echo "Array task ${SLURM_ARRAY_TASK_ID}: processing ${DATASET}"
fi

# ---- Create directories -----------------------------------------------------
mkdir -p "${OUTDIR}"
mkdir -p logs

# ---- Run pipeline -----------------------------------------------------------
echo "Starting CNV pipeline at $(date)"
echo "  base_dir: ${BASE_DIR}"
echo "  outdir:   ${OUTDIR}"

Rscript "${SCRIPT}" \
    --start_from    "block2" \
    --base_dir      "${BASE_DIR}" \
    --metadata_path "${METADATA}" \
    --group_cols    "embryo" \
    --boundaries_mb "25,10" \
    --base_fraction 0.05 \
    --step          0.05 \
    --outdir        "${OUTDIR}"

EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Pipeline completed successfully at $(date)"
else
    echo "Pipeline FAILED with exit code ${EXIT_CODE} at $(date)"
    exit ${EXIT_CODE}
fi