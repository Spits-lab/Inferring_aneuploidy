#!/bin/bash
# slurm/run_round3_benchmark.sh
BASE="/scratch/brussel/vo/000/bvo00016/vsc11567/cnv_framework"

module purge
module load inferCNV/1.18.1-foss-2023a
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load R-bundle-CRAN/2023.12-foss-2023a

sbatch \
  --time=02:25:00 \
  --mem=64G \
  --job-name="r3_benchmark" \
  --output="${BASE}/logs/r3_norefbenchmark_%j.out" \
  --error="${BASE}/logs/r3_norefinfercnv_%j.err" \
  --wrap="
    module purge
    module load inferCNV/1.18.1-foss-2023a
    module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
    module load R-bundle-CRAN/2023.12-foss-2023a
    Rscript ${BASE}/scripts/run_block_test.R
  "

echo "Submitted ✅"
echo "Monitor: tail -f ${BASE}/logs/r3_benchmark_*.out"