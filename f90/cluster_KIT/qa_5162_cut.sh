#!/bin/bash
#SBATCH -p moon
#SBATCH --mail-user=chiara.signorilesign@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH -o /users/ttp/signorile/Desktop/off_Z/mixed_ll/f90/cluster_results/qa_5162_cut.stdout
#SBATCH -e /users/ttp/signorile/Desktop/off_Z/mixed_ll/f90/cluster_results/qa_5162_cut.stderr

cd /users/ttp/signorile/Desktop/off_Z/mixed_ll/f90
./mixed_ll -corr nnlo -ch qa -sec rr5162 -ptlep_cut 20 -ylep_cut 2.5 -vegasNc0 1000000 -vegasNc1 50000000 -vegasIt0 10 -vegasIt1 5
