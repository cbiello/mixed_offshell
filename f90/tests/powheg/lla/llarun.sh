#!/bin/bash
# usage: llarun.sh <boson> <outfile> <nc0> <nc1> [extra flags]  -- all sectors of NLO QCD l l' gamma (13 TeV, fixed scale MZ)
# use nc0 >= 1500000, nc1 >= 5000000 and no technical mass rescaling: with less the sectors come out low (see README)
b=$1; out=$2; nc0=$3; nc1=$4; shift 4; EXTRA="$*"
run(){ r=$(./mixed_ll -boson $b -corr $1 -ch $2 -sec $3 -pref lla${b}_ -histodir llahisto_$b -ecoll 13000 -dynscale .false. -mu 91.1876 -qmin_rsc 1 \
   -vegasIt0 5 -vegasNc0 $nc0 -vegasIt1 5 -vegasNc1 $nc1 $EXTRA </dev/null 2>&1 | grep -a "# Result\|# Error\|STOP" | awk '{print $3}' | tr '\n' ' ')
  echo "$1 $2 $3 $r"; }
export -f run; export b EXTRA nc0 nc1
mkdir -p llahisto_$b
cat tests/powheg/lla/sectors.txt | xargs -P 12 -L 1 bash -c 'run $0 $1 $2' | sort > $out
