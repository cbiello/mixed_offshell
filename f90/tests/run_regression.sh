#!/bin/bash
# Same-seed, low-statistics run of every sector; compares with the stored reference.
# usage (from f90/):  tests/run_regression.sh z|wm|wp [update]
# Z needs analysis=dflt, W needs analysis=W_ATLAS (see the header of the reference).
boson=${1:?usage: tests/run_regression.sh z|wm|wp [update]}
T=$(dirname "$0"); ref=$T/ref_$boson.txt; new=$T/new_$boson.txt
list=sectors_list.dat; [ $boson != z ] && list=sectors_list_$boson.dat
tmp=tests/tmp_$boson; mkdir -p $tmp
run() { out=$(timeout 900 ./mixed_ll $* -boson $boson -histodir $tmp -pref reg_ \
          -vegasIt0 2 -vegasNc0 20000 -vegasIt1 1 -vegasNc1 20000 < /dev/null 2>&1)
  r=$(echo "$out" | grep -aE "# Result:|STOP|runtime error" | head -1 | sed 's/# Result: *//' | tr -s ' ')
  [ -z "$r" ] && r=$(echo "$out" | grep -av '^ *$' | tail -1 | tr -s ' ' | cut -c1-80)
  echo "$2 $4 $6 : $r"; }
export -f run; export boson tmp
{ echo "# analysis: $(grep -o 'ANALYSIS=[^ ]*' config.cfg | tr -d '\\"')"
  sed 's/ *$//' $list | xargs -P $(nproc) -L 1 bash -c 'run "$@"' _ | sort; } > $new
rm -rf $tmp grid/reg_*
[ "$2" = update ] && { mv $new $ref; echo "reference updated: $ref"; exit 0; }
# numbers agree if the relative difference is below 1e-7
paste -d'|' $ref $new | awk -F'|' '
  { n=split($1,a," "); split($2,b," ")
    if ($1==$2) next
    x=a[5]; y=b[5]
    if (a[1]" "a[2]" "a[3]==b[1]" "b[2]" "b[3] && x+0==x && y+0==y && x!=0) { d=(x-y)/x; if (d<0) d=-d; if (d<1e-7) next }
    bad++; print "DIFF  ref: " $1; print "      new: " $2 }
  END { if (bad) { print bad " sectors differ"; exit 1 } else print "all sectors agree with the reference" }'
