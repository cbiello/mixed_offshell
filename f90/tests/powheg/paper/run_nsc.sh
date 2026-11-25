#!/bin/bash
# usage: run_nsc.sh <W_incl|W_ATLASfid> <wp|wm> <corr> <ch> <sec> <seed> [flags]
# one mixed_ll sector, one seed; mu = mT/2 with the three scale columns (mu, mu/2, 2 mu)
R=/scratch2/cbiello/mixed_runs
a=$1; b=$2; c=$3; ch=$4; s=$5; sd=$6; shift 6
n0=1000000; n1=5000000
[ "$c $s" = "nloewk v" ] && { n0=300000; n1=1000000; }
d=prod/nsc/$a/$b; mkdir -p $R/$d/log          # -histodir holds 30 characters: relative path from the build directory
cd $R/build_$a && ./mixed_ll -boson $b -corr $c -ch $ch -sec $s -seed $sd -dynscale .true. -histodir ../$d \
   -vegasIt0 5 -vegasNc0 $n0 -vegasIt1 10 -vegasNc1 $n1 "$@" </dev/null > $R/$d/log/${c}_${ch}_${s}_s$sd.log 2>&1
echo "nsc $a $b $c $ch $s $sd exit $? $(date +%H:%M)"
