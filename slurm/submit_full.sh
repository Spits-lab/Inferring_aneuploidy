#!/bin/bash
#SBATCH --job-name=cnv_pipeline
#SBATCH --time=00:15:00
#SBATCH --mem=64G
#SBATCH --output=/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/logs/pipeline_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework/logs/pipeline_%j.err

module purge
module load inferCNV/1.18.1-foss-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load R-bundle-CRAN/2023.12-foss-2023a


# ---- Required — must be passed at submission time -----------------------
# sbatch --export=COUNTS_PATH=...,METADATA_PATH=...,WORKDIR=... ssubmit_full.sh

if [[ -z "${SCRIPT}" ]];           then echo "ERROR: SCRIPT not set";           exit 1; fi
if [[ -z "${COUNTS_PATH}" ]];      then echo "ERROR: COUNTS_PATH not set";      exit 1; fi
if [[ -z "${METADATA_PATH}" ]];    then echo "ERROR: METADATA_PATH not set";    exit 1; fi
if [[ -z "${GENE_ORDER}" ]];       then echo "ERROR: GENE_ORDER not set";       exit 1; fi
if [[ -z "${TOOL_OUTDIR}" ]];      then echo "ERROR: TOOL_OUTDIR not set";      exit 1; fi
if [[ -z "${WORKDIR}" ]];          then echo "ERROR: WORKDIR not set";          exit 1; fi
if [[ -z "${CHROMOSSOME_PATH}" ]]; then echo "ERROR: CHROMOSSOME_PATH not set"; exit 1; fi
CLONAL_COL="${CLONAL_COL:-NULL}"
DONOR_COL="${DONOR_COL:-NULL}"


mkdir -p logs

echo "=== Submission parameters ==="
echo "  EXECUTION_MODE:              ${EXECUTION_MODE}"
echo "  START_FROM:                  ${START_FROM}"
echo "  COUNTS_PATH:                 ${COUNTS_PATH}"
echo "  METADATA_PATH:               ${METADATA_PATH}"
echo "  WORKDIR:                     ${WORKDIR}"
echo "  GROUP_COLS:                  ${GROUP_COLS}"
echo "  N_SPLITS:                    ${N_SPLITS}"
echo "  REMOVE_REFERENCE:            ${REMOVE_REFERENCE}"
echo "  SENSITIVITY FLOOR MB:        ${SENSITIVITY_FLOOR_MB}"
echo "  P_ARM_PERMISSION:            ${P_ARM_PERMISSION}"
echo "  Q_ARM_PERMISSION:            ${Q_ARM_PERMISSION}"
echo "  WHOLE_CHR_PERMISSION:        ${WHOLE_CHR_PERMISSION}"
echo "  MIN_OVERLAP:                 ${MIN_OVERLAP}"
echo "  CUTOFF:                      ${CUTOFF}"
echo "  SAMPLE_COL:                  ${SAMPLE_COL}"
echo "  BY:                          ${BY_COL}"
echo "  CELL_COL:                    ${CELL_COL}"
echo "  Discrete K Value:            ${K_DIS_VALUE}"
echo "  Frequency Threshold K Value: ${K_FRE_VALUE}"
echo "  Cell Group Approach:         ${CELL_GROUP_CLUSTER}"
echo "  OVERLAP Method:              ${OVERLAP_METHOD}"
echo "  APPLY_GENE_DENSITY_FILTER:   ${APPLY_GENE_DENSITY_FILTER:-FALSE}"
echo "  CODING_GENES_PATH:           ${CODING_GENES_PATH:-NA}"
echo "  PCT_MAX:                     ${PCT_MAX:-45}"
echo "  PCT_FLOOR:                   ${PCT_FLOOR:-30}"
echo "  MIN_EXPR_DENSITY:            ${MIN_EXPR_DENSITY:-1.5}"
echo "  MIN_CODING_DENSITY:          ${MIN_CODING_DENSITY}"
echo "  MAX_GAP_MB:                  ${MAX_GAP_MB}"
echo "============================="

Rscript ${SCRIPT} \
  --execution-mode               "${EXECUTION_MODE}" \
  --start-from                   "${START_FROM}" \
  --counts-path                  "${COUNTS_PATH}" \
  --metadata-path                "${METADATA_PATH}" \
  --gene-order-file              "${GENE_ORDER}" \
  --tool-outdir                  "${TOOL_OUTDIR}" \
  --workdir                      "${WORKDIR}" \
  --chromosome-arms-path         "${CHROMOSSOME_PATH}" \
  --k-interval                   "${K_DIS_VALUE}" \
  --k-value                      "${K_FRE_VALUE}" \
  --sensitivity-floor-mb         "${SENSITIVITY_FLOOR_MB}" \
  --min-required-cells           "${MIN_REQUIRED_CELLS}" \
  --group-cols                   "${GROUP_COLS}" \
  --n-splits-within               "${N_SPLITS}" \
  --remove-reference             "${REMOVE_REFERENCE}" \
  --p-arm-permission             "${P_ARM_PERMISSION}" \
  --q-arm-permission             "${Q_ARM_PERMISSION}" \
  --whole-chr-permission         "${WHOLE_CHR_PERMISSION}" \
  --min-overlap-consistent-calls "${MIN_OVERLAP_CONSISTENT}" \
  --min-overlap-multiple-nodes   "${MIN_OVERLAP_NODES}" \
  --min-overlap                  "${MIN_OVERLAP}" \
  --min-references               "${MIN_REFERENCES}" \
  --max-gap                      "${MAX_GAP}" \
  --cutoff                       "${CUTOFF}" \
  --tool                         "${TOOL}" \
  --sample-col                   "${SAMPLE_COL}" \
  --by-col                       "${BY_COL}" \
  --cell-id-col                  "${CELL_COL}" \
  --group-clusters-col           "${CELL_GROUP_CLUSTER}" \
  --overlap-method               "${OVERLAP_METHOD}" \
  --cnv-annotated                "${CNV_ANNOTATED_PATH}" \
  --cell-sizes                   "${CELL_SIZES_PATH}" \
  --coding-genes-path         "${CODING_GENES_PATH}" \
  --pct-max                   "${PCT_MAX}" \
  --pct-floor                 "${PCT_FLOOR}" \
  --min-expr-density          "${MIN_EXPR_DENSITY}" \
  --min-coding-density        "${MIN_CODING_DENSITY}" \
  --max-gap-mb                "${MAX_GAP_MB}" \
  --donor-col                 "${DONOR_COL}" \
  --clonal-col                "${CLONAL_COL}" \
  
echo "Done: $(date)"

EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Pipeline completed at $(date)"
else
    echo "Pipeline FAILED with exit code ${EXIT_CODE}"
    exit ${EXIT_CODE}
fi




