#!/bin/bash
# Evaluate W_NNLO_polecheck_cleaned.nb and compare its non-q qb' channels with the Fortran. Logs in checks/logs/.
#   verify_W_polecheck_{gq,aq_ag,4q}: notebook S-sector structures = Fortran (transcribed in W_polecheck_fortran_defs.wl)
#   verify_W_polecheck_rv:           ONLO[FLVQEDfin] (gq), ONLO[FLVQCDfin] (photon-quark) = rvewk_is / rvqcd_is
#   verify_W_polecheck_compiled:     W_polecheck_fortran_defs.wl = compiled Fortran (needs a built f90/, see run_verify_fortran_conventions.sh)
# usage: math/subtraction-notebooks/checks/run_verify_W_polecheck.sh [wolframscript]
set -e
WS=${1:-/opt/software/mathematica/14.3/Executables/wolframscript}
C=$(cd "$(dirname "$0")" && pwd); NB=$(dirname "$C"); TOP=$(cd "$NB/../.." && pwd)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
H=$TOP/f90/src/IntSub/mod_hoppet_nnlo.f90
python3 "$C/extract_kernels_fortran.py" "$H" > "$T/k.m"
python3 "$C/extract_4q_fortran.py" "$H" > "$T/k4q.m"
# compiled values: private copies of the HOPPET kernels + the objects of mixed_ll
( cd "$TOP/f90"
  { echo "module kern_copy"; sed -n 2,10p $H; echo "  implicit none"; echo "contains"
    for f in mySub_gq_func mySub_gq_Lmu_func mySub_gq_b_func mySub_gq_b_Lmu_func mySub_TCgq_func mySub_TCgq_Lmu_func mySub_TCgq_Lmu2_func \
             mySub_ns_qq_func mySub_ns_qq_Lmu_func mySub_ns_qqb_func mySub_ns_qqb_Lmu_func; do
      awk -v f="$f" '$0 ~ "^ *function "f"\\(" {p=1} p {print} p && $0 ~ "end function "f" *$" {p=0}' $H; done
    echo "end module kern_copy"; } > "$T/kern_copy.f90"
  link=$(make -n -W src/main.f90 | grep gfortran | tail -1 | tr ' ' '\n' | grep -E '^(objects/|-L|-l|-I|"-Wl)' | tr -d '"' | tr '\n' ' ')
  gfortran -cpp @config.cfg -O2 -ffree-line-length-1000 -I modules -J "$T" "$T/kern_copy.f90" "$C/fortran_polecheck_driver.f90" $link -o "$T/driver" )
( cd "$T" && ./driver && sed -i 's/E\([-+]\)/*^\1/g' fortran_polecheck_values.m )
NBDIR=$NB NBOUT=$T "$WS" -file "$C/run_W_NNLO_polecheck.wls" 2>&1 | grep -v configuaration > "$C/logs/run_W_NNLO_polecheck.log"
for n in gq aq_ag 4q rv compiled; do
  NBOUT=$T FKERN=$T/k.m F4Q=$T/k4q.m FVALS=$T/fortran_polecheck_values.m "$WS" -file "$C/verify_W_polecheck_${n}_fortran.wls" 2>&1 \
    | grep -v -e configuaration -e meprec -e '^ ' -e '^$' -e indet > "$C/logs/verify_W_polecheck_${n}_fortran.log"
done
cat "$C"/logs/verify_W_polecheck_*_fortran.log
