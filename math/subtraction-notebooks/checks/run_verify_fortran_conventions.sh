#!/bin/bash
# usage: checks/run_verify_fortran_conventions.sh <outdir>   (from math/subtraction-notebooks; needs a built f90/mixed_ll and wolframscript)
# fortran_conventions.wl against the compiled Fortran: the hoppet kernels are copied verbatim (private in mod_hoppet_nnlo),
# everything else is called from the objects of mixed_ll. Prints relative differences, all ~1e-15.
set -e
out=$(realpath -m $1); mkdir -p $out; nb=$PWD; cd ../../f90
K=src/IntSub/mod_hoppet_nnlo.f90
{ echo "module kern_copy"; sed -n 2,10p $K; echo "  implicit none"; echo "contains"
  for f in mySub_ns_func mySub_ns_Lmu_func PqqNLO_func PqqNLO_Lmu_func PqqNLO_b_func PqqNLO_b_Lmu_func mySub_TCqq_func mySub_TCqq_Lmu_func mySub_TCqq_Lmu2_func; do
    awk -v f="$f" '$0 ~ "^ *function "f"\\(" {p=1} p {print} p && $0 ~ "end function "f"$" {p=0}' $K; done
  echo "end module kern_copy"; } > $out/kern_copy.f90
link=$(make -n -W src/main.f90 | grep gfortran | tail -1 | tr ' ' '\n' | grep -E '^(objects/|-L|-l|-I|"-Wl)' | tr -d '"' | tr '\n' ' ')
gfortran -cpp @config.cfg -O3 -ffree-line-length-1000 -I modules -J $out $out/kern_copy.f90 $nb/checks/fortran_conventions_driver.f90 $link -o $out/driver
cd $out && ./driver && sed -i 's/E\([-+]\)/*^\1/g' fortran_values.m
NBDIR=$nb FVALS=$out/fortran_values.m wolframscript -file $nb/checks/verify_fortran_conventions.wls
