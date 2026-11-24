#!/bin/bash
# usage: f90/tests/checklim_run.sh <boson> <npts> <csvdir>   (from any directory; PLAN C1.5, see tests/CHECKLIM.md)
# -checklim_csv for every NLO real, RV, RR and ONLO sector of sectors_list_<boson>.dat (not four-quark, nf, ns s_oewk_is)
b=$1; n=$2; out=$(realpath -m "$3"); mkdir -p "$out"
cd "$(dirname "$0")/.." || exit 1
grep -h -- '-sec \(r_\|rv\|rr_\|s_o\)' sectors_list_$b.dat | awk '{print $2, $4, $6}' | grep -v 'ns_qq\|nf_is\|^nnlo ns s_oewk_is$' | sort -u |
  xargs -P 4 -L 1 bash -c "./mixed_ll -boson $b -corr \$0 -ch \$1 -sec \$2 -checklim_csv -checklim_npts $n </dev/null 2>&1 | grep -a '^[PRU],' > '$out'/$b.\$0.\$1.\$2.csv"
