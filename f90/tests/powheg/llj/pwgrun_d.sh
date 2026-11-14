#!/bin/bash
# usage: pwgrun_d.sh <tag> <code dir> <nseeds | first-last> <grid file> key=val ...
# as pwgrun.sh, but every seed starts from an existing importance-sampling grid (pwgxgrid-btl.dat):
# needed for runs with btildeborn=0 (NLO correction only), whose own grid stage sees a zero integrand
tag=$1; code=$2; n=$3; grid=$4; shift 4
case $n in *-*) a=${n%-*}; b=${n#*-};; *) a=1; b=$n;; esac
for i in $(seq $a $b); do d=$tag/s$i; rm -rf $d; mkdir -p $d; cp base.input $d/powheg.input; cp $grid $d/pwgxgrid-btl.dat
  for kv in "$@" iseed=$i; do echo "${kv/=/ }" >> $d/powheg.input; done
  ( cd $d && $code/pwhg_main > run.log 2>&1 ) & done; wait
