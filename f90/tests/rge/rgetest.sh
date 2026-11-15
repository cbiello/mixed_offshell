#!/bin/bash
# RGE test of the integrated (S) sectors: with PDFs and alpha_s frozen at mu0, the slope in log(mu^2) of the mixed sectors
# of a channel must equal minus the lower-order results with one PDF replaced by a splitting function convolved with it.
# ns: the two-loop sector is not run; its slope is minus the row V (delta(1-x) part of the two-loop splitting function).
# usage (from f90/): tests/rge/rgetest.sh <z|wp|wm> <gq|aq|ag|ns> [nc0 nc1]     then: tests/rge/rge.py tests/rge/out_<boson>_<channel>
b=$1; c=$2; nc0=${3:-200000}; nc1=${4:-500000}; mu=91.1876; out=tests/rge/out_${b}_$c; mkdir -p $out
run() { p=rg${1}_${2}_${3}_; shift 3; ./mixed_ll -boson $b "$@" -dynscale .false. -mu $mu -pdfscale $mu -pref $p -histodir $out \
          -vegasIt0 3 -vegasNc0 $nc0 -vegasIt1 3 -vegasNc1 $nc1 < /dev/null > /dev/null 2>&1; rm -f grid/$p*; }
lhs() { for ch in $1; do for s in $2; do echo "L $ch $s -corr nnlo -ch $ch -sec $s"; done; done; }
rhs() { for ch in $3; do for s in $4; do echo "$1 $ch $s -corr $2 -ch $ch -sec $s -rge $5"; done; done; }
export -f run; export b mu out nc0 nc1
case $c in
gq) lhs "gq qg" "s_oqcd s_oewk_is s_oewk_fs_53 s_oewk_fs_54 s_12 s_vewk"
    rhs A nloewk ns "r_is s v r_fs_53 r_fs_54" qg; rhs B nloqcd "gq qg" "r_is s" qq; rhs D lo ns na qg11;;
aq) lhs "aq qa" "s_vqcd s_12 s_oqcd s_oewk"
    rhs A nloqcd ns "r_is s v" qa; rhs B nloewk "aq qa" "r_is s" qqs; rhs D lo ns na qa11;;
ag) lhs "ag ga" "s_12 s_oqcd s_oewk"
    rhs A nloqcd "gq qg" "r_is s" qa; rhs B nloewk "aq qa" "r_is s" qg;;
ns) lhs ns "s_12 s_vqcd s_vewk s_oqcd s_oewk_is_ra s_oewk_fs_53 s_oewk_fs_54 s_qqb s_qq"
    rhs A nloewk ns "r_is s v r_fs_53 r_fs_54" qqs; rhs B nloqcd ns "r_is s v" qq; rhs D lo ns na qq11; rhs E lo ns na qb11; rhs V lo ns na qq11d;;
esac | xargs -P $(nproc) -L 1 bash -c 'run "$@"' _
