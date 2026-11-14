#!/bin/bash
# usage: lljrun.sh <boson> <outfile> <nc0> <nc1> [extra flags]  -- all sectors of NLO EW l l' + jet (13 TeV, fixed scale MZ)
# default technical mass rescaling: the photon recombination needs events generated below the analysis cuts
b=$1; out=$2; nc0=$3; nc1=$4; shift 4; EXTRA="$*"
run(){ r=$(./mixed_ll -boson $b -corr $1 -ch $2 -sec $3 -pref llj${b}_ -histodir lljhisto_$b -ecoll 13000 -dynscale .false. -mu 91.1876 \
   -vegasIt0 5 -vegasNc0 $nc0 -vegasIt1 5 -vegasNc1 $nc1 $EXTRA </dev/null 2>&1 | grep -a "# Result\|# Error\|STOP" | awk '{print $3}' | tr '\n' ' ')
  echo "$1 $2 $3 $r"; }
export -f run; export b EXTRA nc0 nc1
mkdir -p lljhisto_$b
cat tests/powheg/llj/sectors.txt | xargs -P 12 -L 1 bash -c 'run $0 $1 $2' | sort > $out
