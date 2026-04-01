#!/bin/bash
#SBATCH -p moon
#SBATCH --mail-user=chiara.signorilesign@gmail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH -o /users/ttp/signorile/Desktop/DY-offshell/mixed_offshell/qq_ns_5161a_Z_old.stdout
#SBATCH -e /users/ttp/signorile/Desktop/DY-offshell/mixed_offshell/qq_ns_5161a_Z_old.stderr

cd /users/ttp/signorile/Desktop/DY-offshell/mixed_offshell/f90
./mixed_ll -corr nnlo -ch ns_ga -sec rr_5262a -vegasNc0 5000000 -vegasNc1 500000000 -vegasIt0 20 -vegasIt1 10 
