#!/bin/bash
# One sector on the KIT cluster (SLURM): a warm-up job writes the Vegas grid, then a job array runs the production
# seeds 1..n on that grid (It1 > 1: Vegas keeps refining it between iterations); the array starts when the warm-up has finished.
#   cluster_KIT/submit.sh <boson> <corr> <ch> <sec> [nseeds=100] [Nc0=5000000] [It0=15] [Nc1=50000000] [It1=1]
# Environment: PART partition (default moon), BIN binary (default: a copy of mixed_ll in bin/ named by its checksum, so a rebuild
# does not change queued jobs), LOGS log directory (default ../logs; logs and grids are named after the analysis in makefile.inc). Histograms: histo/*_s<seed>.histo.
set -e
b=$1 corr=$2 ch=$3 sec=$4 n=${5:-100} nc0=${6:-5000000} it0=${7:-15} nc1=${8:-50000000} it1=${9:-1}
f90=$(cd "$(dirname "$0")/.." && pwd); part=${PART:-moon}; ana=$(sed -n 's/^analysis=//p' "$f90/makefile.inc")
logs=$(realpath -m "${LOGS:-$f90/../logs}/$ana/$b")
mkdir -p "$logs" "$f90/grid" "$f90/bin"
bin=$(cd "$f90" && realpath -m "${BIN:-bin/mixed_ll_$(md5sum < mixed_ll | cut -c1-8)}"); [ -e "$bin" ] || cp "$f90/mixed_ll" "$bin"
g=${ana}_${b}_${corr}_${ch}_${sec}; run="$bin -boson $b -corr $corr -ch $ch -sec $sec -gridfile $g"
w=$(sbatch --parsable -p "$part" -c 1 -D "$f90" -J "w_$g" -o "$logs/${g}_warmup.out" \
    --wrap "$run -writegrid -vegasNc0 $nc0 -vegasIt0 $it0 -vegasNc1 1000 -vegasIt1 1 </dev/null")
p=$(sbatch --parsable -p "$part" -c 1 -D "$f90" -J "p_$g" -a "1-$n" -d "afterok:$w" -o "$logs/${g}_s%a.out" \
    --wrap "$run -readgrid -vegasNc1 $nc1 -vegasIt1 $it1 -seed \$SLURM_ARRAY_TASK_ID </dev/null")
echo "$g: warm-up $w, seeds 1-$n $p ($bin, $part)"
