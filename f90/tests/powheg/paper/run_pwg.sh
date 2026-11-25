#!/bin/bash
# usage: run_pwg.sh <lo|nloqcd|nloew> <wp|wm> <scale factor> <seed> [key=val ...]
# one POWHEG seed, setup A and B in one run (nscsetup 3), mu = xi mT/2
R=/scratch2/cbiello/mixed_runs; P=/scratch/cbiello/POWHEG-BOX-RES/MiNNLOPS_res
o=$1; b=$2; xi=$3; sd=$4; shift 4
case $o in
  lo)     code=NSC-lnuQCD; keys="bornonly=1 ncall1=300000 itmx1=3 ncall2=2000000 itmx2=5";;
  nloqcd) code=NSC-lnuQCD; keys="ncall1=300000 itmx1=3 ncall2=1000000 itmx2=5";;
  nloew)  code=NSC-lnuEW;  keys="ncall1=100000 itmx1=3 ncall2=600000 itmx2=6";;
esac
dy=1; [ $b = wm ] && dy=-1
d=$R/prod/pwg/$b/${o}_mu$xi/s$sd; rm -rf $d; mkdir -p $d; cp $R/pwg/base.input $d/powheg.input
declare -A K; for kv in chooseDY=$dy renscfact=$xi facscfact=$xi $keys "$@" iseed=$sd; do K[${kv%%=*}]=${kv#*=}; done   # later keys win
for k in "${!K[@]}"; do echo "$k ${K[$k]}" >> $d/powheg.input; done
[ -n "$FIXMU" ] && sed -i "s/^fixedscale 0/fixedscale 1/" $d/powheg.input
cd $d && $P/$code/pwhg_main_paper > run.log 2>&1
echo "pwg $o $b $xi $sd exit $? $(date +%H:%M)"
