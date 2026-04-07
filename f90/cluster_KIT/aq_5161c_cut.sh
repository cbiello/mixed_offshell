#!/bin/bash
#SBATCH -p moon
#SBATCH --mail-user=chiara.signorilesign@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH -o /users/ttp/signorile/Desktop/off_Z/mixed_ll/f90/cluster_results/aq_5161c_cut.stdout
#SBATCH -e /users/ttp/signorile/Desktop/off_Z/mixed_ll/f90/cluster_results/aq_5161c_cut.stderr

cd /users/ttp/signorile/Desktop/off_Z/mixed_ll/f90
./mixed_ll -corr nnlo -ch aq -sec rr5161c -vegasNc0 2000000 -vegasNc1 50000000 -vegasIt0 10 -vegasIt1 5 -ptlep_cut 20 -ylep_cut 2.5
