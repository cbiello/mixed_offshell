# CHECKS — validated comparisons

Record of the checks that passed, with the numbers. Open tasks are in [PLAN.md](PLAN.md).

## POWHEG at LO and NLO, W+ and W- (2026-11-12)

Codes: `NSC-lnuQCD`, `NSC-lnuEW` (POWHEG-BOX-RES, fixed order) against `mixed_ll`.
Recipe, scripts, input and sector numbers: [f90/tests/powheg/README.md](f90/tests/powheg/README.md).

Setup: pp at 13 TeV, PDF 325100 (NNPDF31 luxqed), fixed scale mu_R = mu_F = 91.1876 GeV, G_mu scheme,
complex masses, MW = 80.398, GammaW = 2.1054, MZ = 91.1876, GammaZ = 2.4952 GeV, unit CKM.
Cuts (`W_ATLAS`): mT in (200,5000) GeV, pT(l) > 65, pT(nu) > 85, pT(l) pT(nu) > 35^2, |y(l)|, |y(nu)| < 2.4;
charged lepton dressed with a photon within Delta R < 0.1.

| | POWHEG [fb] | mixed_ll [fb] | POWHEG - mixed_ll [fb] |
|---|---|---|---|
| LO W+ | 1417.2 +- 1.5 | 1416.41 +- 0.18 | 0.8 +- 1.5 |
| LO W- | 990.73 +- 0.47 | 990.88 +- 0.12 | -0.15 +- 0.49 |
| NLO QCD W+ (ns + gq + qg) | 1737.51 +- 0.61 | 1736.68 +- 0.31 | 0.83 +- 0.69 |
| NLO QCD W- (ns + gq + qg) | 1203.14 +- 0.44 | 1203.02 +- 0.18 | 0.12 +- 0.48 |
| NLO EW W+ (ns + aq + qa) | 1327.29 +- 0.25 | 1326.99 +- 0.21 | 0.30 +- 0.33 |
| NLO EW W- (ns + aq + qa) | 929.24 +- 0.25 | 929.23 +- 0.13 | 0.01 +- 0.28 |

All six agree within 1.4 standard deviations; NLO EW to 0.02%.

Notes on the POWHEG side:
- `NSC-lnuEW/openloops2.f` needs the real and the virtual divided by `st_alpha/(2pi)` (the public POWHEG-BOX-RES
  multiplies by alpha_s/2pi also in QED-only mode) and the scale passed to OpenLoops (`set_parameter("mu",..)`).
  With the pushed version (`eb7c4b3`) the NLO total is wrong by a large factor.
- Run with `detailedNLO 1`, otherwise the photon is not identified and the lepton is not dressed
  (W+ undressed: 1305.28 +- 0.76 fb).
- Use `pwhg_analysis-NSC.f` (in `f90/tests/powheg/`); the committed `pwhg_analysis-W_ATLAS.f` is disabled.

## Limits of the RV and RR sectors of the l l' + photon QCD corrections (2026-11-12)

Tool: `./mixed_ll ... -checkint -vegasIt0 3 -vegasNc0 30000` (integrability scan, see STRUCTURE.md).
A limit passes if the median of |f| d^k over 100 points falls with the depth d (k = number of variables at the edge).

Sectors (63): `ns rvqcd_is`, `rvqcd_fs_53/54`; `aq`/`qa rvqcd_is`; RR `ns_ga` (10), `gq` (7), `qg` (7), `aq` (10), `qa` (10), `ag` (7), `ga` (7).

| Boson | Sectors | Limits per sector | Result |
|---|---|---|---|
| Z (`dflt`) | 63 | 16 (RV), 49 (RR) | all pass |
| W+ (`W_ATLAS`) | 63 | same | all pass |
| W- (`W_ATLAS`) | 63 | same | all pass |

Before the fixes the scan failed for W in `ns rvqcd_is` (photon collinear to leg 2) and in the neutrino sector of
`ns rvqcd_fs` (PLAN RV6a-RV7b). Sensitivity, from planted bugs in `aq rr_5161a`: a wrong charge in any counterterm
and a 1% error in the soft or in a nested counterterm are all flagged.

Z, old code vs general code after the conversion of `aq`, `qa`, `ag`, `ga` (RR) and `aq`, `qa` (`rvqcd_is`):
`tests/run_regression.sh z` reproduces the old numbers (1e-7) in every sector except `ag rr_5262a` (2.5e-6, rounding amplified
by the grid adaptation; 6e-11 with one warm-up iteration).

Pairwise check of `aq`/`qa rvqcd_is` for W (`-checklim`): median |(H+CT)/H| = 1e-2, 1e-3, 1e-4 at depth 1e-4, 1e-6, 1e-8.

## W+ q qbar' integrated subtraction: Fortran against `W.nb` (2026-11-12)

Scripts in `math/subtraction-notebooks/checks/` (Mathematica 14.3, `/opt/software/mathematica/14.3/Executables/wolframscript`).

1. `run_W_nb.wls` re-evaluates the derivation of `W.nb` (RR + RV + VV + PDF counterterms, q qb' channel, general charges)
   with the `aux/` files: the 1/eps^4, 1/eps^3, 1/eps^2 and 1/eps coefficients are all 0; the finite remainder is saved.
2. `verify_Wp_qq_fortran.wls` contains a transcription of the Fortran: the `Vcharge == 1` branches of `s_12`, `s_vqcd`,
   `s_vewk`, the functions `fin_ns_vqcd_cc`, `genGnloQCD_cc`, `legfinalcalG_cc`, `fin_ns_v7_cc` and the kernels
   `mySub_ns`, `PqqNLO(_b)(_Lmu)`, `mySub_TCqq(_Lmu,_Lmu2)`, `Pqq_NLO`, `PqqAP_0_R`. It is compared with the notebook
   remainder, with charge conservation Qq + Qqp = Qlp, structure by structure:

| Structure | Compared | Difference |
|---|---|---|
| FLVqcd elastic, leg 1, leg 2 | coefficients of D0..D3 and the regular part at three z | 0 |
| FLVew elastic, leg 1, leg 2 | same | 0 |
| FLM double boost | all D(z1) D(z2) coefficients on a 2x2 grid | 0 |
| FLM single boost, leg 1 and leg 2 | D0..D3 and regular part | 0 |
| FLM elastic | number at a random point (Born kinematics, eta12 -> 1) | 1e-30 (the limit parameter) |

Sensitivity: changing one coefficient in `legfinalcalG_cc` or in `fin_ns_v7_cc` gives differences of order 1.
Not covered: the ONLO terms (`s_oqcd`, `s_oewk_*`), which involve the partition functions; W- (not in the notebook, to be
obtained with 3 <-> 4); the identification of the kinematic variables in the code, which is shared with the Z implementation.

## NLO QCD corrections to l l' gamma against POWHEG (2026-11-13)

All sectors with a resolved photon (39 per boson: LO from `nloewk ns r_*`, `ns rvqcd_*`, RR `ns_ga`/`gq`/`qg`,
`s_oewk_*` of `ns`/`gq`/`qg`) against the new POWHEG generator `NSC-lnuaEW`. Setup, cuts, recipe and the full sector
table: [f90/tests/powheg/lla/README.md](f90/tests/powheg/lla/README.md).

| | POWHEG [fb] | mixed_ll [fb] | difference [fb] |
|---|---|---|---|
| NLO QCD W+ | 6.1937 +- 0.0075 | 6.1928 +- 0.0021 | -0.0009 +- 0.0078 |
| NLO QCD W- | running | 3.7968 +- 0.0012 | |
| NLO QCD Z | running | 13.1031 +- 0.0031 | |
| LO W+ | 5.0308 +- 0.0015 | 5.0291 +- 0.0009 | -0.0017 +- 0.0018 |
| LO W- | 3.0387 +- 0.0011 | 3.0376 +- 0.0005 | -0.0011 +- 0.0012 |
| LO Z | 10.5800 +- 0.0026 | 10.5766 +- 0.0014 | -0.0034 +- 0.0029 |

W- and Z with the damping off are being rerun; with the gluon PDF switched off W- already agrees (3.4048 +- 0.0011 against 3.4066 +- 0.0039).
The check covers the W fixes of `rvqcd` (RV6a-RV7b), the RR sectors and the resolved-photon S sectors for q qb, gq
and qg. The `_gen` trees with g + gamma agree with the stored MadGraph points for W+ and W- up to one constant factor
per coupling order (1.00069 at alpha^2, 1.00115 at alpha^3: inputs of the stored points; not followed up).

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

## W- q qbar' integrated subtraction (2026-11-13)

`s_12`, `s_vqcd`, `s_vewk` use one W branch for both signs with `ilept_chg`, `Qlept_chg` (W+: 4, +1; W-: 3, -1): the
derivation of `W.nb` does not depend on which leg carries the charge. CP test (`tests/cptest.sh`): W- in p p against W+
with antiproton PDFs (`-cptest`), 14 TeV, `W_ATLAS`:

| Sector | W- [fb] | W+, pbar pbar [fb] | relative difference |
|---|---|---|---|
| lo ns | 1096.73 +- 0.21 | 1096.98 +- 0.21 | 2e-4 |
| nloqcd ns v | -208.199 +- 0.041 | -208.251 +- 0.040 | 3e-4 |
| nnlo ns s_vqcd | 0.71284 +- 0.00019 | 0.71282 +- 0.00019 | 4e-5 |
| nnlo ns s_12 | -1.7546 +- 0.0010 | -1.7553 +- 0.0010 | 4e-4 |
| nnlo ns s_vewk | two runs, combined difference -0.009 +- 0.019 on -22.5 | | 4e-4 |

Resolved-photon S sectors (`ns s_oewk_is_ra`, `s_oewk_fs_53/54`; `gq`, `qg` `s_oewk_is`, `s_oewk_fs_53/54`): the QCD kernels
in the code are the ONLO-EWK terms of the notebooks (q qb: `W.nb`, CF [Ptilde + (2 D0 - (1+z)) log(eta) omega] and
CF (2 pi^2/3 - 3 lmu); gq: `Z.nb`, TR [Ptilde_gq + P_qg log(eta/2) omega]); all nine pass the `-checkint` scan for W+ and W-.
`gq s_oewk_fs_53/54` was zero at every point (two debug lines overwrote the integration variables): fixed, Z gq = qg mirror
(0.0925 fb each in the l l gamma setup).
