#!/bin/bash
#SBATCH -p moon
#SBATCH --job-name=qq_ns_5261
#SBATCH --array=1-100
#SBATCH --mail-user=chiara.signorilesign@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH -o /users/ttp/signorile/Desktop/DY-offshell/mixed_offshell/logs/OLD_Z/qq_ns_%A_%a.stdout
#SBATCH -e /users/ttp/signorile/Desktop/DY-offshell/mixed_offshell/logs/OLD_Z/qq_ns_%A_%a.stderr

# Go to working directory
cd /users/ttp/signorile/Desktop/DY-offshell/mixed_offshell/f90

# Use SLURM array index as seed
SEED=${SLURM_ARRAY_TASK_ID}

echo "Running with seed ${SEED}"

./mixed_ll_Z_old \
  -corr nnlo \
  -ch ns_ga \
  -sec rr_5261 \
  -vegasNc0 5000000 \
  -vegasNc1 50000000 \
  -vegasIt0 15 \
  -vegasIt1 5 \
  -seed ${SEED}
