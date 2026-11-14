# NLO EW corrections to l l' + jet: mixed_ll against POWHEG `NSC-lnujQCD`

Setup: 13 TeV, PDF 325100, fixed scale 91.1876 GeV, G_mu scheme, inputs of `base.input`.
Channels: q qb', g q, q g (no photon-induced channels, no four-quark interference).
Cuts (analysis `V_llj_nloew` here, `pwhg_analysis-NSCj.f` in POWHEG):
a photon is recombined with the nearest charged lepton or final quark if Delta R < 0.1 (never with a gluon);
jet = final quark (dressed) or gluon with pT > 30, no rapidity cut;
W: `W_ATLAS` cuts on the dressed lepton; Z: mll > 200, pT(l) > 20, |y(l)| < 2.5, pT(l+) pT(l-) > 35^2.

## Run

    # mixed_ll, built with analysis=V_llj_nloew, from f90/
    tests/powheg/llj/lljrun.sh wp out_wp.txt 1500000 5000000      # also wm, z
    # POWHEG (build and keys: README in NSC-lnujQCD); LO and correction separately
    ../pwgrun.sh lo_1 <path>/NSC-lnujQCD 10 chooseDY=1 bornonly=1 ncall1=200000 itmx1=3 ncall2=500000 itmx2=6
    ../pwgrun.sh grid_1 <path>/NSC-lnujQCD 1 chooseDY=1 ncall1=40000 itmx1=3 ncall2=2000 itmx2=1
    ./pwgrun_d.sh d_1 <path>/NSC-lnujQCD 20 grid_1/s1/pwgxgrid-btl.dat chooseDY=1 btildeborn=0 ncall2=100000 itmx2=4

`sectors.txt`: the hard events of `nloqcd r_is` are the LO; `rvewk_is` the EW loop; `ns_ga`, `gq`, `qg` `rr_*`
the real photon emission; `s_oqcd` the integrated counterterms and collinear renormalisation.
The default technical mass rescaling is kept (the recombination needs events generated below the cuts).
The Z `rvewk_is` sectors cost 7 ms per point (16 times the W ones): run them with fewer points over several seeds.

## Results (2026-11-14)

| | channel | POWHEG [fb] | mixed_ll [fb] | difference [fb] |
|---|---|---|---|---|
| LO W+ | all | 555.51 +- 0.20 | 555.24 +- 0.09 | -0.27 +- 0.22 |
| LO W+ | q qb' | 307.54 +- 0.12 | 307.36 +- 0.08 | -0.18 +- 0.14 |
| LO W+ | gluon | 247.97 +- 0.12 | 247.88 +- 0.04 | -0.08 +- 0.13 |
| LO W- | all | 365.99 +- 0.10 | 365.75 +- 0.05 | -0.24 +- 0.11 |
| LO W- | q qb' | 187.10 +- 0.06 | 186.95 +- 0.05 | -0.15 +- 0.07 |
| LO W- | gluon | 178.89 +- 0.07 | 178.80 +- 0.03 | -0.09 +- 0.07 |
| LO Z | all | 585.61 +- 0.18 | 585.25 +- 0.09 | -0.36 +- 0.21 |
| LO Z | q qb' | 291.37 +- 0.07 | 291.13 +- 0.08 | -0.24 +- 0.11 |
| LO Z | gluon | 294.24 +- 0.14 | 294.11 +- 0.05 | -0.12 +- 0.15 |
| correction W+ | all | -38.86 +- 0.08 | -38.82 +- 0.02 | 0.04 +- 0.09 |
| correction W+ | q qb' | -22.36 +- 0.04 | -22.36 +- 0.02 | 0.00 +- 0.05 |
| correction W+ | gluon | -16.50 +- 0.05 | -16.46 +- 0.01 | 0.04 +- 0.05 |
| correction W- | all | -25.01 +- 0.05 | -24.93 +- 0.01 | 0.08 +- 0.05 |
| correction W- | q qb' | -13.15 +- 0.03 | -13.13 +- 0.01 | 0.02 +- 0.03 |
| correction W- | gluon | -11.86 +- 0.03 | -11.79 +- 0.01 | 0.06 +- 0.03 |
| correction Z | all | -18.63 +- 0.07 | -18.64 +- 0.01 | -0.00 +- 0.07 |
| correction Z | q qb' | -8.76 +- 0.03 | -8.81 +- 0.01 | -0.05 +- 0.03 |
| correction Z | gluon | -9.86 +- 0.05 | -9.83 +- 0.01 | 0.03 +- 0.05 |

q qb' and gluon refer to the incoming partons (`rate_qq`, `rate_g` in POWHEG).
POWHEG: 10 seeds (LO), 20 seeds (correction; Z: 28 seeds, error-weighted, 16 further seeds did not finish in 12 h).
mixed_ll: 1.5M / 5M points (trees; LO of W: two seeds), 200k / 500k (`rvewk_is`; Z: plus six seeds of 20k / 50k).

The corrections agree within 2 sigma, i.e. within 0.02% of the LO. The LO of mixed_ll is 0.05% below POWHEG for the
three bosons (about 2 sigma each): not understood, to be checked with `-qmin_rsc 1` and tighter POWHEG generation cuts.

Two errors in the integrated sectors were found with this comparison (W only, both fixed):

- `ns s_oqcd` used a soft function with the wrong sign of the lepton-quark charge terms (`calG_ONLOQCD_ns`, now
  replaced by `genGnloQCD_cc`): the q qb' correction was W+ -26.86, W- -15.71 fb.
- `calG_ONLOQCD_gq` filled only the (g,q) element, so its term was missing in `qg s_oqcd`
  (W+ -0.081 instead of 0.359 fb; gq 0.358).

mixed_ll sectors:

| corr ch sec | W+ [fb] | W- [fb] | Z [fb] |
|---|---|---|---|
| nloqcd gq r_is | 123.9421 +- 0.0265 | 89.3846 +- 0.0184 | 147.0620 +- 0.0350 |
| nloqcd ns r_is | 307.3563 +- 0.0776 | 186.9505 +- 0.0458 | 291.1322 +- 0.0809 |
| nloqcd qg r_is | 123.9427 +- 0.0264 | 89.4115 +- 0.0183 | 147.0511 +- 0.0348 |
| nnlo gq rr_5161a | -0.0878 +- 0.0001 | -0.0757 +- 0.0001 | -0.2239 +- 0.0001 |
| nnlo gq rr_5161b | -0.1395 +- 0.0003 | -0.0917 +- 0.0003 | -0.1813 +- 0.0004 |
| nnlo gq rr_5161c | -0.2159 +- 0.0002 | -0.1916 +- 0.0002 | -0.5653 +- 0.0004 |
| nnlo gq rr_5161d | -0.1373 +- 0.0005 | -0.1120 +- 0.0006 | -0.4336 +- 0.0007 |
| nnlo gq rr_5162 | -0.0916 +- 0.0012 | -0.1689 +- 0.0005 | -1.0066 +- 0.0010 |
| nnlo gq rr_5163 | -0.3463 +- 0.0006 | -0.2842 +- 0.0016 | -1.1793 +- 0.0020 |
| nnlo gq rr_5164 | -0.3144 +- 0.0022 | -0.2483 +- 0.0003 | -1.0747 +- 0.0019 |
| nnlo gq rvewk_is | -7.2608 +- 0.0073 | -5.1930 +- 0.0049 | -1.8108 +- 0.0021 |
| nnlo gq s_oqcd | 0.3584 +- 0.0009 | 0.4591 +- 0.0003 | 1.5623 +- 0.0007 |
| nnlo ns rvewk_is | -19.2261 +- 0.0199 | -11.1469 +- 0.0110 | -5.9718 +- 0.0054 |
| nnlo ns s_oqcd | -0.9493 +- 0.0049 | -0.6184 +- 0.0023 | 4.9524 +- 0.0031 |
| nnlo ns_ga rr_5161a | -0.0265 +- 0.0014 | -0.0126 +- 0.0009 | -0.3668 +- 0.0009 |
| nnlo ns_ga rr_5161c | -0.2962 +- 0.0013 | -0.1783 +- 0.0003 | -0.6637 +- 0.0007 |
| nnlo ns_ga rr_5162 | -0.1737 +- 0.0011 | -0.0665 +- 0.0008 | -0.9575 +- 0.0012 |
| nnlo ns_ga rr_5163 | -0.3260 +- 0.0004 | -0.2163 +- 0.0019 | -0.9395 +- 0.0015 |
| nnlo ns_ga rr_5164 | -0.2686 +- 0.0026 | -0.2100 +- 0.0003 | -0.9726 +- 0.0017 |
| nnlo ns_ga rr_5261 | -0.1738 +- 0.0008 | -0.0666 +- 0.0005 | -0.9536 +- 0.0012 |
| nnlo ns_ga rr_5262a | -0.0281 +- 0.0006 | -0.0131 +- 0.0005 | -0.3656 +- 0.0013 |
| nnlo ns_ga rr_5262c | -0.2961 +- 0.0005 | -0.1784 +- 0.0003 | -0.6628 +- 0.0009 |
| nnlo ns_ga rr_5263 | -0.3257 +- 0.0004 | -0.2156 +- 0.0014 | -0.9350 +- 0.0016 |
| nnlo ns_ga rr_5264 | -0.2706 +- 0.0022 | -0.2099 +- 0.0003 | -0.9751 +- 0.0019 |
| nnlo qg rr_5261 | -0.0935 +- 0.0010 | -0.1689 +- 0.0004 | -1.0061 +- 0.0008 |
| nnlo qg rr_5262a | -0.0879 +- 0.0001 | -0.0754 +- 0.0001 | -0.2237 +- 0.0001 |
| nnlo qg rr_5262b | -0.1398 +- 0.0003 | -0.0915 +- 0.0002 | -0.1809 +- 0.0004 |
| nnlo qg rr_5262c | -0.2156 +- 0.0002 | -0.1918 +- 0.0002 | -0.5655 +- 0.0004 |
| nnlo qg rr_5262d | -0.1387 +- 0.0006 | -0.1112 +- 0.0006 | -0.4323 +- 0.0007 |
| nnlo qg rr_5263 | -0.3464 +- 0.0005 | -0.2823 +- 0.0019 | -1.1830 +- 0.0019 |
| nnlo qg rr_5264 | -0.3128 +- 0.0030 | -0.2487 +- 0.0003 | -1.0751 +- 0.0022 |
| nnlo qg rvewk_is | -7.2450 +- 0.0073 | -5.1784 +- 0.0049 | -1.8114 +- 0.0020 |
| nnlo qg s_oqcd | 0.3592 +- 0.0008 | 0.4602 +- 0.0004 | 1.5657 +- 0.0007 |
