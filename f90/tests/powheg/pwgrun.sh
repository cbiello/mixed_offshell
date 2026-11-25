#!/bin/bash
# usage: [PWG=<binary>] pwgrun.sh <tag> <code dir> <nseeds | first-last> key=val ...   (parallel independent seeds)
tag=$1; code=$2; n=$3; shift 3
case $n in *-*) a=${n%-*}; b=${n#*-};; *) a=1; b=$n;; esac
for i in $(seq $a $b); do d=$tag/s$i; rm -rf $d; mkdir -p $d; cp base.input $d/powheg.input
  for kv in "$@" iseed=$i; do echo "${kv/=/ }" >> $d/powheg.input; done
  ( cd $d && $code/${PWG:-pwhg_main} > run.log 2>&1 ) & done; wait
