#!/bin/bash
# usage: pwgrun.sh <tag> <code dir> <nseeds> key=val ...   (parallel independent seeds)
tag=$1; code=$2; n=$3; shift 3
for i in $(seq 1 $n); do d=$tag/s$i; rm -rf $d; mkdir -p $d; cp base.input $d/powheg.input
  for kv in "$@" iseed=$i; do echo "${kv/=/ }" >> $d/powheg.input; done
  ( cd $d && $code/pwhg_main > run.log 2>&1 ) & done; wait
