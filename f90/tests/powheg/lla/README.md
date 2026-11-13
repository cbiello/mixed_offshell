# NLO QCD corrections to l l' gamma: mixed_ll against POWHEG `NSC-lnuaEW`

Setup: 13 TeV, PDF 325100, fixed scale 91.1876 GeV, G_mu scheme, inputs of `base.input`.
Cuts (analysis `V_lla_nloqcd` here, `pwhg_analysis-NSCa.f` in POWHEG):
photon pT > 50, |y| < 2.5, Delta R(photon, charged lepton) > 0.1, Delta R(photon, final quark) > 0.4;
pT(l l') > 30; W: `W_ATLAS` cuts on the bare lepton; Z: mll > 200, pT(l) > 20, |y(l)| < 2.5, pT(l+) pT(l-) > 35^2.

## Run

    # mixed_ll, built with analysis=V_lla_nloqcd, from f90/
    tests/powheg/lla/llarun.sh wp out_wp.txt 1500000 5000000      # also wm, z; then sum the fourth column
    # POWHEG (see the README in NSC-lnuaEW for the build and the generation cuts); base.input has withdamp 0
    ../pwgrun.sh nlo_1 <path>/NSC-lnuaEW 20 chooseDY=1 ncall1=200000 itmx1=3 ncall2=500000 itmx2=6

The sectors are listed in `sectors.txt`: the hard events of `nloewk ns r_*` are the LO; `rvqcd_*` the virtual;
`ns_ga`, `gq`, `qg` `rr_*` the real; `s_oewk_*` the integrated counterterms and collinear renormalisation.

## Results (2026-11-13)

| | POWHEG [fb] | mixed_ll [fb] | difference [fb] |
|---|---|---|---|
| NLO QCD W+ | 6.1937 +- 0.0075 | 6.1928 +- 0.0021 | -0.0009 +- 0.0078 |
| NLO QCD W- | running | 3.7968 +- 0.0012 | |
| NLO QCD Z | running | 13.1031 +- 0.0031 | |
| LO W+ | 5.0308 +- 0.0015 | 5.0291 +- 0.0009 | -0.0017 +- 0.0018 |
| LO W- | 3.0387 +- 0.0011 | 3.0376 +- 0.0005 | -0.0011 +- 0.0012 |
| LO Z | 10.5800 +- 0.0026 | 10.5766 +- 0.0014 | -0.0034 +- 0.0029 |

W- and Z with the damping off are being rerun; with the gluon PDF switched off W- already agrees (3.4048 +- 0.0011 against 3.4066 +- 0.0039).

Three things were needed to get here (all in the recipe now):

- `-qmin_rsc 1`: the default technical mass rescaling (20) makes Vegas unreliable in the initial-state photon sectors
  for this hard-photon setup (`nloewk ns r_is` W+: 1.495 to 1.517 fb over seeds, quoted error 0.0015).
- At least 1.5M / 5M points per iteration (warm-up / main): with 300k / 1M many sectors come out 1-3% low with
  errors that do not show it.
- POWHEG without damping (`withdamp 0`, now the default of `NSC-lnuaEW`): with damping the hard real emission with the
  photon close to the final quark goes to the remnant, whose integral scatters between 0 and 57 fb over seeds; the
  W+ rate was 6.1622 +- 0.0052, i.e. 0.5% low.

How the last point was found (W+, fb; useful as partial checks):

| Quantity | mixed_ll | POWHEG |
|---|---|---|
| q qb' channel only (gluon PDF off) | 5.6531 +- 0.0018 | 5.6586 +- 0.0052 (with damping) |
| real only, 30 < pT(l l') < 50 | 0.9843 +- 0.0013 | 0.9870 +- 0.0043 (with damping) |
| real, final gluon pT > 20 | 2.4284 +- 0.0009 | 2.4286 +- 0.0046 |
| real, final quark pT > 20 | 1.1865 +- 0.0005 | 1.1870 +- 0.0017 (with damping: 1.160 +- 0.008) |
| real, final quark pT > 50 | 0.7363 +- 0.0003 | 0.7365 +- 0.0012 (with damping: 0.710 +- 0.007) |

mixed_ll sectors (two seeds, 1.5M / 5M points, `-qmin_rsc 1`):

| corr ch sec | W+ [fb] | W- [fb] | Z [fb] |
|---|---|---|---|
| nloewk ns r_fs_53 | 0.36285 +- 0.00010 | 1.89574 +- 0.00041 | 3.41338 +- 0.00074 |
| nloewk ns r_fs_54 | 3.14873 +- 0.00071 | 0.21854 +- 0.00006 | 3.48379 +- 0.00075 |
| nloewk ns r_is | 1.51754 +- 0.00063 | 0.92337 +- 0.00030 | 3.67943 +- 0.00090 |
| nnlo gq rr_5161a | 0.01367 +- 0.00001 | 0.01093 +- 0.00001 | 0.02333 +- 0.00001 |
| nnlo gq rr_5161b | 0.02722 +- 0.00002 | 0.02205 +- 0.00002 | 0.04259 +- 0.00003 |
| nnlo gq rr_5161c | 0.01131 +- 0.00004 | 0.01283 +- 0.00003 | 0.03772 +- 0.00008 |
| nnlo gq rr_5161d | 0.02279 +- 0.00002 | 0.02060 +- 0.00002 | 0.06217 +- 0.00004 |
| nnlo gq rr_5162 | 0.02673 +- 0.00009 | 0.01994 +- 0.00005 | 0.02846 +- 0.00012 |
| nnlo gq rr_5163 | 0.02639 +- 0.00007 | 0.06109 +- 0.00020 | 0.12270 +- 0.00034 |
| nnlo gq rr_5164 | 0.08093 +- 0.00024 | 0.02087 +- 0.00005 | 0.11074 +- 0.00031 |
| nnlo gq s_oewk_fs_53 | 0.00494 +- 0.00002 | 0.02685 +- 0.00006 | 0.04346 +- 0.00010 |
| nnlo gq s_oewk_fs_54 | 0.04435 +- 0.00009 | 0.00115 +- 0.00001 | 0.03414 +- 0.00010 |
| nnlo gq s_oewk_is | 0.01150 +- 0.00007 | -0.00033 +- 0.00004 | 0.01490 +- 0.00012 |
| nnlo ns rvqcd_fs_53 | -0.06553 +- 0.00002 | -0.37460 +- 0.00008 | -0.66801 +- 0.00015 |
| nnlo ns rvqcd_fs_54 | -0.62248 +- 0.00014 | -0.03927 +- 0.00001 | -0.68185 +- 0.00015 |
| nnlo ns rvqcd_is | -0.25691 +- 0.00011 | -0.15327 +- 0.00005 | -0.66747 +- 0.00016 |
| nnlo ns s_oewk_fs_53 | 0.13953 +- 0.00018 | 0.66938 +- 0.00067 | 1.25715 +- 0.00113 |
| nnlo ns s_oewk_fs_54 | 1.12780 +- 0.00109 | 0.08398 +- 0.00010 | 1.29372 +- 0.00118 |
| nnlo ns s_oewk_is_ra | 0.61581 +- 0.00077 | 0.37790 +- 0.00043 | 1.50014 +- 0.00118 |
| nnlo ns_ga rr_5161a | -0.08210 +- 0.00015 | -0.05444 +- 0.00007 | -0.16717 +- 0.00018 |
| nnlo ns_ga rr_5161c | -0.00241 +- 0.00058 | -0.00153 +- 0.00022 | -0.00423 +- 0.00073 |
| nnlo ns_ga rr_5162 | -0.02635 +- 0.00023 | -0.01630 +- 0.00015 | -0.06057 +- 0.00028 |
| nnlo ns_ga rr_5163 | -0.00147 +- 0.00013 | -0.02440 +- 0.00027 | -0.01687 +- 0.00065 |
| nnlo ns_ga rr_5164 | -0.04427 +- 0.00058 | -0.00146 +- 0.00006 | -0.02527 +- 0.00074 |
| nnlo ns_ga rr_5261 | -0.02662 +- 0.00023 | -0.01637 +- 0.00010 | -0.06102 +- 0.00030 |
| nnlo ns_ga rr_5262a | -0.08237 +- 0.00013 | -0.05446 +- 0.00007 | -0.16759 +- 0.00017 |
| nnlo ns_ga rr_5262c | -0.00277 +- 0.00053 | -0.00178 +- 0.00023 | -0.00336 +- 0.00054 |
| nnlo ns_ga rr_5263 | -0.00156 +- 0.00011 | -0.02474 +- 0.00035 | -0.01676 +- 0.00053 |
| nnlo ns_ga rr_5264 | -0.04433 +- 0.00064 | -0.00156 +- 0.00006 | -0.02531 +- 0.00080 |
| nnlo qg rr_5261 | 0.02681 +- 0.00009 | 0.01987 +- 0.00005 | 0.02873 +- 0.00019 |
| nnlo qg rr_5262a | 0.01366 +- 0.00001 | 0.01092 +- 0.00001 | 0.02334 +- 0.00001 |
| nnlo qg rr_5262b | 0.02723 +- 0.00002 | 0.02210 +- 0.00002 | 0.04259 +- 0.00003 |
| nnlo qg rr_5262c | 0.01134 +- 0.00004 | 0.01283 +- 0.00003 | 0.03760 +- 0.00007 |
| nnlo qg rr_5262d | 0.02278 +- 0.00002 | 0.02060 +- 0.00002 | 0.06217 +- 0.00004 |
| nnlo qg rr_5263 | 0.02645 +- 0.00008 | 0.06121 +- 0.00021 | 0.12320 +- 0.00042 |
| nnlo qg rr_5264 | 0.08077 +- 0.00024 | 0.02094 +- 0.00004 | 0.11080 +- 0.00031 |
| nnlo qg s_oewk_fs_53 | 0.00489 +- 0.00001 | 0.02680 +- 0.00006 | 0.04339 +- 0.00010 |
| nnlo qg s_oewk_fs_54 | 0.04433 +- 0.00009 | 0.00114 +- 0.00001 | 0.03412 +- 0.00010 |
| nnlo qg s_oewk_is | 0.01159 +- 0.00007 | -0.00032 +- 0.00004 | 0.01487 +- 0.00012 |
