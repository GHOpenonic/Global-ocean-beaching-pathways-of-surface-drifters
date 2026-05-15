#!/bin/bash
#SBATCH -p pi_abodner
#SBATCH -w node4005
#SBATCH -N 1
#SBATCH --job-name=watershed_stability_3regime_200_lower_diag
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
#SBATCH --cpus-per-task=125
#SBATCH --mem=350G
#SBATCH --time=4:00:00

set -euo pipefail

mkdir -p logs

BEACH_INPUT="/home/codycruz/drifters_watersheds/undrogued_beach.parquet"
OUTDIR="/home/codycruz/drifters_watersheds/watershed_stability_outputs/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"

source /home/codycruz/drifters_watersheds/.venv/bin/activate

python watershed_clustering_stability_regimes.py \
  --beach-parquet "${BEACH_INPUT}" \
  --outdir "${OUTDIR}" \
  --n-workers 125 \
  --threads-per-worker 1 \
  --trajectory-history-weights "0.1,0.15" \
  --trajectory-sample-counts "50" \
  --regime1-min-cluster-sizes "20,21,22,23,24,25,26,27,28,29" \
  --regime1-cluster-selection-epsilons "6,6.5,7,7.5,8,8.5,9,9.5,10,10.5" \
  --regime2-min-cluster-sizes "20,21,22,23,24,25,26,27,28,29" \
  --regime2-cluster-selection-epsilons "11,11.5,12,12.5,13,13.5,14,14.5,15.5,16" \
  --regime3-min-cluster-sizes "30,31,32,33,34,35,36,37,38,39" \
  --regime3-cluster-selection-epsilons "11,11.5,12,12.5,13,13.5,14,14.5,15.5,16" \
  --coast-detour-weight 0.1 \
  --coast-grid-resolution 0.5 \
  --watershed-existence-jaccard-threshold 0.3 \
  --noise-grouped-threshold 0.3 \
  --min-consensus-cluster-size 25 \
  --min-consensus-existence-probability 0.3

  # --regime1-min-cluster-sizes "20,21,22,23,24,25,26,27,28,29" \
  # --regime1-cluster-selection-epsilons "1,1.5,2,2.5,3,3.5,4,4.5,5,5.5" \
  # --regime2-min-cluster-sizes "30,31,32,33,34,35,36,37,38,39" \
  # --regime2-cluster-selection-epsilons "6,6.5,7,7.5,8,8.5,9,9.5,10,10.5" \
  # --regime3-min-cluster-sizes "40,41,42,43,44,45,46,47,48,49" \
  # --regime3-cluster-selection-epsilons "11,11.5,12,12.5,13,13.5,14,14.5,15.5,16" \

  #   --regime1-min-cluster-sizes "30,31,32,33,34,35,36,37,38,39" \
  # --regime1-cluster-selection-epsilons "1,1.5,2,2.5,3,3.5,4,4.5,5,5.5" \
  # --regime2-min-cluster-sizes "40,41,42,43,44,45,46,47,48,49" \
  # --regime2-cluster-selection-epsilons "1,1.5,2,2.5,3,3.5,4,4.5,5,5.5" \
  # --regime3-min-cluster-sizes "40,41,42,43,44,45,46,47,48,49" \
  # --regime3-cluster-selection-epsilons "6,6.5,7,7.5,8,8.5,9,9.5,10,10.5" \