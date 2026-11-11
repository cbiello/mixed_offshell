#!/bin/bash
# usage: nscrun.sh <outfile> <boson> "<corr ch sec>" ...   (13 TeV, fixed scale 91.1876)
out=$1; b=$2; shift 2
run(){ r=$(./mixed_ll -boson $b -corr $1 -ch $2 -sec $3 -pref nsc${b}_ -histodir nschisto -ecoll 13000 -dynscale .false. -mu 91.1876 \
   -vegasIt0 5 -vegasNc0 500000 -vegasIt1 10 -vegasNc1 2000000 $EXTRA </dev/null 2>&1 | grep -a "# Result\|# Error" | awk '{print $3}' | tr '\n' ' ')
  echo "$1 $2 $3 $r"; }
export -f run; export b EXTRA
mkdir -p nschisto
printf "%s\n" "$@" | xargs -P 16 -L 1 bash -c 'run $0 $1 $2' | sort > $out
python3 - $out <<'P'
import sys,math
s=e=0
for l in open(sys.argv[1]):
    v=l.split(); print(f"  {v[0]:7s}{v[1]:4s}{v[2]:9s}{float(v[3]):14.5f} +- {float(v[4]):9.5f} fb"); s+=float(v[3]); e+=float(v[4])**2
print(f"  sum                 {s:14.5f} +- {math.sqrt(e):9.5f} fb")
P
