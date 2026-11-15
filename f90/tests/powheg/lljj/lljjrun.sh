#!/bin/bash
# usage: lljjrun.sh <boson> <outfile> <sector list> <nc0> <nc1> [extra flags]  -- four-quark sectors with two jets (13 TeV, fixed scale MZ)
b=$1; out=$2; list=$3; nc0=$4; nc1=$5; shift 5; EXTRA="$*"
run(){ r=$(./mixed_ll -boson $b -corr $1 -ch $2 -sec $3 -pref jj${b}_ -histodir jjhisto_$b -ecoll 13000 -dynscale .false. -mu 91.1876 \
   -vegasIt0 5 -vegasNc0 $nc0 -vegasIt1 5 -vegasNc1 $nc1 $EXTRA </dev/null 2>&1 | grep -a "# Result\|# Error\|STOP" | awk '{print $3}' | tr '\n' ' ')
  echo "$1 $2 $3 $r"; }
export -f run; export b EXTRA nc0 nc1
mkdir -p jjhisto_$b
cat $list | xargs -P 10 -L 1 bash -c 'run $0 $1 $2' | sort > $out
