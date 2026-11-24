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
| NLO QCD W+ | 6.1896 +- 0.0049 | 6.1928 +- 0.0021 | 0.0032 +- 0.0053 |
| NLO QCD W- | 3.8011 +- 0.0038 | 3.7968 +- 0.0012 | -0.0043 +- 0.0040 |
| NLO QCD Z | 13.1107 +- 0.0112 | 13.1031 +- 0.0031 | -0.0076 +- 0.0116 |
| LO W+ | 5.0308 +- 0.0015 | 5.0291 +- 0.0009 | -0.0017 +- 0.0018 |
| LO W- | 3.0387 +- 0.0011 | 3.0376 +- 0.0005 | -0.0011 +- 0.0012 |
| LO Z | 10.5800 +- 0.0026 | 10.5766 +- 0.0014 | -0.0034 +- 0.0029 |

POWHEG NLO: mean over 50 (W+) and 30 (W-, Z) seeds without damping; seeds more than 0.2 fb from the median are dropped (1 for W+, 2 for W-, none for Z: single large-weight points).

With the gluon PDF switched off W- also agrees (3.4048 +- 0.0011 against 3.4066 +- 0.0039).
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

## NLO EW corrections to l l' + jet against POWHEG (2026-11-14)

Sum of the 33 sectors with one resolved jet (`nloqcd r_is` as LO, `rvewk_is`, RR `ns_ga`/`gq`/`qg`, `s_oqcd`)
against the new POWHEG generator `NSC-lnujQCD`, channels q qb', gq, qg; photon recombined with leptons and quarks
(Delta R 0.1), jet pT > 30. Setup, recipe and sector table: [f90/tests/powheg/llj/README.md](f90/tests/powheg/llj/README.md).

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

The corrections agree within 2 sigma (0.02% of the LO); the LO of mixed_ll is 0.05% low (2 sigma), still to be understood.
Before the fixes the q qb' correction was W+ -26.86, W- -15.71 fb. Found with this check:

- `ns s_oqcd` (W): wrong sign of the lepton-quark charge terms in `calG_ONLOQCD_ns`; now `genGnloQCD_cc`.
- `qg s_oqcd` (W): `calG_ONLOQCD_gq` filled only the (g,q) element. Mirror after the fix: gq 0.3584(9), qg 0.3592(8) fb (W+).
- `gq`/`qg rvewk_is` for W: new, collinear limit with the extra -2 pi^2 Q_i (Q3+Q4); integrability scan passes.
- `qg rvewk_is` (Z, old code): quark and antiquark rows of `qg_lumi_splitb` swapped; mirror after the fix gq -1.8129(26), qg -1.8085(25) fb.

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

## Four-quark channels (2026-11-14)

W+/W- `ns_qqb`, `ns_qq` (`rr_5161a/c`, `rr_5262a/c`): hard matrix element from OpenLoops `pplnjj_ew`, triple-collinear
counterterm from the Born `_gen` with the squared charge of the emitting leg. `-checklim` (pointwise TC limit): the
difference closes as sqrt(depth) in all 16 sectors (slowest points, W- `ns_qqb rr_5262a`: 0.92, 0.52, 0.10 at 1e-6, 1e-8, 1e-10).
The `-checkint` scan is not sensitive here (logarithmic singularity only: a mutated counterterm passes it).

Integrated `s_qqb`, `s_qq`: the generic code run for Z reproduces the legacy `s_qqb` to 14 digits. Legacy Z `s_qq` had quark
and antiquark exchanged in the Born of the leg-1 term; per-leg mirror test with a charge-asymmetric cut: -2.707e-3 vs
-0.703e-3 fb before, -2.716e-3 vs -2.711e-3 fb after the fix.

Hard events with two jets against POWHEG `NSC-lnujjMIX` (13 TeV, both partons pT > 30, Delta R > 0.4; `f90/tests/powheg/lljj/`):

| | initial state | POWHEG [fb] | mixed_ll [fb] | difference [fb] |
|---|---|---|---|---|
| W+ | opposite sign | -0.5683 +- 0.0008 | -0.5690 +- 0.0010 | -0.0008 +- 0.0013 |
| W+ | same sign | 0.4628 +- 0.0042 | 0.4598 +- 0.0013 | -0.0031 +- 0.0044 |
| W- | opposite sign | -0.3864 +- 0.0008 | -0.3862 +- 0.0009 | 0.0002 +- 0.0012 |
| W- | same sign | 0.2644 +- 0.0035 | 0.2614 +- 0.0006 | -0.0030 +- 0.0036 |
| Z | opposite sign | -0.4496 +- 0.0006 | -0.4488 +- 0.0007 | 0.0007 +- 0.0009 |
| Z | same sign | 0.3080 +- 0.0028 | 0.3086 +- 0.0008 | 0.0006 +- 0.0030 |

Found with this check: Z `ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w` registered the OpenLoops processes with u and d exchanged
(shifts of 3 to 6 sigma per channel in this setup; in the regression run `ns_qqpb_w` goes from -0.199 to -0.051 fb). Fixed.
The integrated `s_qqb`, `s_qq` and the counterterms are not covered by this benchmark.

Further checks of the W four-quark sectors (2026-11-14, `W_ATLAS`):

- Regression: Z, W+, W- all sectors agree with the references after the fixes.
- `-checkint` scan: the 8 RR sectors and `s_qqb`, `s_qq` pass for W+ and W- (weak test for RR, see above).
- CP test (`tests/cptest.sh`, W- in p p against W+ in pbar pbar):

| Sector | W- [fb] | W+, pbar pbar [fb] | difference / error |
|---|---|---|---|
| ns s_qqb | 0.140551 +- 0.000035 | 0.140535 +- 0.000035 | 0.3 |
| ns s_qq | -0.004520 +- 0.000002 | -0.004518 +- 0.000002 | -0.7 |
| ns_qqb rr_5161a + c | -0.20815 +- 0.00068 | -0.20907 +- 0.00065 | 1.0 |
| ns_qqb rr_5262a + c | -0.20854 +- 0.00073 | -0.20975 +- 0.00068 | 1.2 |
| ns_qq rr_5161a, c, 5262a, c | 0.0273, 0.0662, 0.0271, 0.0672 | 0.0276, 0.0677, 0.0275, 0.0675 | -1.1, -2.2, -1.2, -0.5 |
| ns_qqb_w rr | -0.0468 +- 0.0069 | -0.0542 +- 0.0056 | 0.8 |
| ns_qqp_w rr | 0.1261 +- 0.0016 | 0.1263 +- 0.0010 | -0.1 |

  For `ns_qqb` only the sum a + c is CP symmetric: the sectors order the quark (5) against the antiquark (6), and CP exchanges
  the two (a alone differs by 6%, 8 sigma).
- Leg mirror of the integrated sectors (convolution on leg 1 against leg 2, symmetric cuts): `s_qqb` W+ 0.105941(64) vs
  0.105994(65), W- 0.070265(42) vs 0.070241(41); `s_qq` W+ -3.8955(5)e-4 vs -3.8895(5)e-4, W- -2.2623(22)e-3 vs -2.2575(21)e-3 fb.

Consistency of the W counterterms and integrated sectors with the Z ones (2026-11-15):

- The W code run for Z reproduces the old Z code at the regression points: TC counterterm of the 8 RR sectors to 12 digits,
  `s_qqb` and `s_qq` to 14 digits. The old Z branches are removed, Z now runs the same code as W there.
- Kernels against the notebook (`math/subtraction-notebooks/checks/verify_4q_fortran.wls`): `mySub_ns_qqb`, `mySub_ns_qq` and
  their `_Lmu` parts equal `FINITEqqbqqb`, `FINITEqqqq` of `Z.nb` (ratio 1 to 15 digits at six values of z).
- The scale term of `s_qq` is -2 P_A(z), P_A = 4(1-z) + 2(1+z) log z + 2 p_qq(-z) S2(z): the known two-loop q -> qbar
  splitting function with CF (CF - CA/2) -> 2 CF Q^2, as required by the collinear renormalisation.

With the pointwise TC limit and the POWHEG comparison of the hard events, what is left untested for W is only what is common
to Z: that the notebook kernels are the integral of the TC counterterm (`Z.nb` holds the results, not the derivation).

## Integrated sectors of gq/qg and of the photon-induced channels for W (2026-11-15)

New W code: gq, qg `s_12`, `s_vewk`; aq, qa `s_vqcd`, `s_12`, `s_oewk` (qa also `s_oqcd`); ag, ga `s_12`, `s_oqcd`, `s_oewk`.
Z regression unchanged. No benchmark and no notebook for these sectors: the checks below test the wiring, not the formulas.
The gq/qg `s_12` formula for W is not derived (PLAN S10).

CP test (W- in p p against W+ in pbar pbar, 4 seeds each, `W_ATLAS`, 14 TeV), difference / error:

| Sector | gq | qg | aq | qa | ag | ga |
|---|---|---|---|---|---|---|
| s_12 | -0.2 | -0.0 | -1.1 | 0.6 | 0.3 | -0.2 |
| s_vewk | 0.8 | -1.0 | | | | |
| s_vqcd | | | -0.4 | -1.1 | | |
| s_oqcd | | | 1.3 | -1.5 | 0.5 | -1.4 |
| s_oewk | | | 1.8 | -0.8 | 1.0 | 0.2 |

Mirror (W-, fb): gq `s_12` 0.13523(30), qg 0.13485(30); gq `s_vewk` 1.6511(4), qg 1.6509(4); aq `s_12` -5.237(4)e-3, qa -5.237(4)e-3;
ag `s_12` 3.278(1)e-3, ga 3.277(1)e-3.

Seeds: before I18 (PLAN) runs with different `-seed` shared the generator of the lepton rapidity sign. With 8 seeds gq `s_12`
came out 0.13595(30) and qg 0.13404(30) (4.5 sigma apart); each charge structure of the boosted term taken alone was symmetric.
After the fix the two agree (numbers above). The four-quark CP result for `ns_qqb` (only a + c symmetric) is unchanged with independent seeds.

## Scale dependence of the integrated sectors (RGE test, 2026-11-15)

With PDFs and alpha_s frozen at mu0 (`-pdfscale`), the slope c1 in log mu^2 of the mixed sectors of a channel must equal -r0,
r0 = lower-order results with one PDF replaced by a splitting function convolved with it (`-rge`): P_qg, P_qq, P_q gamma at
one loop times the NLO QCD / NLO EW sectors, and the mixed kernels P^(1,1) times LO. The finite remainders (`v`, `rvewk`,
`rvqcd`, the hard events) do not depend on mu, so the test covers every S sector of the channel at once.
Both sides come from the scale-variation weights of single runs, at the same points. Recipe: `f90/tests/rge/rgetest.sh`, `rge.py`.
The two-loop sector cannot run here (library not linked): its slope is taken as minus the delta(1-x) part of P_qq^(1,1) (row V).

Calibration at NLO: gq -104.89 vs -104.884; Z ns QCD 65.330 vs 65.324, QED 0.9643 vs 0.9642.

Mixed sectors, fb per unit of log mu^2 (`W_ATLAS`, 14 TeV, mu0 = MZ, one seed):

| Channel | W+ c1 | W+ -r0 | W- c1 | W- -r0 | Z c1 | Z -r0 |
|---|---|---|---|---|---|---|
| gq + qg | 6.6211 | 6.6122 | 5.2229 | 5.2117 | 4.163 | 4.182 |
| aq + qa | 0.0653 | 0.0668 | 0.03676 | 0.03682 | | |
| ag + ga | -0.04335 | -0.04323 | -0.04318 | -0.04330 | | |
| q qb (ns) | -3.778 | -3.838 | -3.188 | -3.192 | -1.927 | -1.923 |

- W+ ns: the gap is noise of `s_oewk_fs_54` (slope -0.412 here; -0.465, -0.499, -0.496, -0.506 with four other seeds).
  With the mean of the five seeds c1 = -3.842.
- Z gq + qg without `s_vewknf`; the 0.45% gap is within the noise of `s_oewk_fs`.
- Single sectors that can be isolated: W+ gq + qg `s_vewk` 5.231 vs 5.230 (P_qg times the EW virtual); Z `s_qq` -0.00851 vs -0.00850
  (P_q qb^(1,1)); the log^2 mu coefficient of W+ gq + qg also matches.

What it tests: all coefficients of log mu in the S sectors, including the products of one-loop kernels with the finite soft
and virtual functions at boosted kinematics. For gq/qg `s_12` of W this includes the W-specific term of `calG_gq_cc`
(lepton-quark charge correlation), which enters r0 through the NLO EW sectors convolved with P_qg.
What it does not test: terms without log mu (constants such as the zeta2 terms, the finite parts of the triple-collinear kernels).

## Integrated subtraction of gq and of the photon channels: Fortran against `Z.nb` (2026-11-15)

Scripts in `math/subtraction-notebooks/checks/`, logs in `checks/logs/`.

1. `verify_Z_gq_aq_ag_fortran.wls` (Z). The kernels `mySub_gq*`, `mySub_TCgq*`, `mySub_TCaa*`, `mySub_TCag*` are read from
   `mod_hoppet_nnlo.f90` by `extract_kernels_fortran.py`; the assembly of `s_12`, `s_vewk`, `s_vqcd` and `calG_QqQl` are
   transcribed. Compared with `FINITEgq`, `FINITEγq`, `FINITEγg` of `Z.nb`, all terms without a resolved parton,
   at six values of z (40 digits). Every difference is 0:

   | Channel | Terms |
   |---|---|
   | gq | `s_vewk` kernel; `s_12` double boost, triple-collinear kernel with its log and log^2 parts, soft-photon term |
   | photon-quark | `s_vqcd` kernel; `s_12` double boost, single boost q qb, single boost gamma gamma |
   | photon-gluon | `s_12` double boost, single boost gamma gamma |

   No term of the notebook is left unassigned. The terms with a resolved parton (`s_oqcd`, `s_oewk*`) are in item 5.

2. What this means for W. Photon channels: the W code uses the same kernel tables, with the squared charge of the quark
   from the photon splitting and without the gamma gamma terms; no other charge enters these terms, so `s_12` and `s_vqcd`
   of aq, qa, ag, ga are covered. gq, qg: `s_vewk`, the double boost and the triple-collinear kernel are covered in the
   same way (charge of one quark line squared). The soft-photon term is the only W-specific formula.

3. `verify_W_gq_fortran.wls`: the W soft-photon term (`calG_gq_cc`) written for two charged leptons and reduced to Z
   charges reproduces the `Z.nb` single boost: the Qq^2 and Qq Qe terms agree identically, log(z) term included, and what
   is left is a function of E3, E4, eta34, L only (the lepton-lepton correlation, absent for W). With the opposite sign of
   the lepton charge the Qq Qe terms do not close.

4. `derive_gq_general_charges.wls`: the gq derivation of `notes/fabrizio/subtraction_f.wl` (RR, RV and PDF counterterms of the
   nested subtraction, written with the charges q1..q4 of the four legs) is evaluated without specifying the charges,
   only q1 + q2 + q3 + q4 = 0.
   - All poles of the F(z.p1,p2) coefficient cancel for general charges.
   - Z charges: the finite part equals the Fortran single boost (`calG_QqQl` and the `mySub_TCgq` kernels), difference 0.
   - W charges (charged lepton on leg 4 or on leg 3, any q1, q2): the finite part equals the Fortran W formula
     (`calG_gq_cc` and the same kernels), difference 0 in each of the structures q1^2, q1 q2, q2^2 at five values of z.
   - The double boost carries q2^2 only and the virtual term no charge, as in the code.
   - Sensitivity: without the -4 zeta2 Q_i^2 term of `calG_gq_cc` the difference is of order 20.

   This closes PLAN S10: the W formula of gq/qg `s_12` is derived, constants included. The derivation uses the charges of
   the external legs only (soft and collinear photon emission, Catani operator), which is all that enters for W as well.

5. `derive_aq_ag_general_charges.wls`: the photon-quark and photon-gluon derivations of the same file, charges left general.
   - All poles cancel, in the terms with a resolved parton and in the others.
   - Charges: the terms with a resolved quark from the photon (`s_oewk`) carry no charge besides the amplitude; the terms
     with the photon splitting integrated (`s_oqcd`, single and double boost, QCD virtual) carry only the squared charge of
     the quark from the photon splitting; the gamma gamma term is the only one with another charge (absent for W).
   - Finite parts of the resolved-parton terms against the Fortran (transcribed), each difference 0:

   | Sector | Term | Fortran |
   |---|---|---|
   | aq `s_oqcd` | P~_qg(z) with the damped log(eta15/2), times Nc Q^2 | `intsub` of `xsect_nnlo_s_aq_oqcd` |
   | aq `s_oewk`, boosted | CF [P~_qq(z) with the damped log(eta25/2)], plus distributions included | `Pqq_NLO`, `PqqAP_0_R` |
   | aq `s_oewk`, elastic | CF [G - 3/2 log(mu^2/E^2)] | `calG_CF_hard` |
   | ag `s_oqcd` | Nc Q^2 P~_qg(z) | `tildePqgNLO` |
   | ag `s_oewk` | TR [P~_qg(z) with the damped log(eta25/2)] | `intsub` of `xsect_nnlo_s_ag_oewk` |

   Sensitivity: 13/2 -> 6 in the elastic function gives a difference of 2/3. qa, ga are the mirrors (aq = qa, ag = ga tests).
   Not in this comparison: the local counterterms of the resolved parton inside these sectors (covered by the scan below).

## Integrability scan of `rvewk` and of the S sectors with a resolved parton (2026-11-16)

`-checkint -vegasIt0 3 -vegasNc0 30000` on `rvewk_is` (ns, gq, qg), ns `s_oqcd`, `s_oewk_is_ra`, `s_oewk_fs_53/54`,
gq/qg `s_oqcd`, `s_oewk_is`, `s_oewk_fs_53/54`, aq/qa/ag/ga `s_oqcd`, `s_oewk`: 23 sectors for Z, W+, W-.
All 69 runs pass (10 to 14 limits each closing, no FAIL; the remaining lines have too few non-zero points to judge).
This completes PLAN C1.6.

## Tree amplitudes against MadGraph (2026-11-16)

`./mixed_ll -boson <z|wp|wm> -corr <lo|nloqcd|nloew|qcdew> -ch ns -sec na -MGcheck -cm_scheme false`, 20 phase-space
points per channel, tolerance 1e-6.

| Boson | lo | nloqcd | nloew | qcdew (g + photon, all crossings incl. ag, ga) |
|---|---|---|---|---|
| W+ | 40/40 | 120/120 | 120/120 | 240/240 |
| W- | 40/40 | 120/120 | 120/120 | 240/240 |
| Z | 80/80 | 240/240 | 240/240 | 440/440 |

The stored MadGraph points were made with real couplings. With the default complex-mass scheme every W entry is off by
the same factor (0.99931470 = 1/(1 + Gamma_W^2/M_W^2) at lo and nloqcd, 0.99885166 with a photon, spread 1.5e-9), which is
the coupling normalisation only. This closes PLAN RR9.

## One convention for notebooks and Fortran: q qb' against `W.nb` (2026-11-24)

`math/subtraction-notebooks/fortran_conventions.wl` fixes the Fortran convention for the notebooks (legs `p1..p4`,
`Qq = Q_IS(al)`, `Qqp = Q_IS(be)`, `Q3`, `Q4`, charged lepton `Qlept` on leg `ilept`), maps `W.nb` (`fromWnb[ilept]`)
and the all-outgoing charges of `notes/fabrizio` (`fromAllOut`), and writes the q qb' finite remainder per sector with
the Fortran routines and `ilept` as a parameter. `checks/bench_fortran_conventions_W.wls` (Mathematica 13.2, `W.nb`
evaluated with `run_W_nb.wls`): all ten structures agree for W+ (`ilept = 4`) and W- (`ilept = 3`), general charges.
Sensitivity: `W.nb` charges read as incoming, or W- without the leg swap, give O(1) differences. Log in `checks/logs/`.

Second link, `fortran_conventions.wl` against the compiled Fortran (`checks/run_verify_fortran_conventions.sh`): a driver
linked to the objects of `mixed_ll` evaluates the hoppet kernels of `s_12`, `s_vqcd`, `s_vewk` (pieces REAL, VIRT, DELTA;
copied verbatim, they are private), `Pqq_NLO`, `PqqAP_0_R` and `fin_ns_vqcd_cc`, `genGnloQCD_cc`, `legfinalcalG_cc`,
`fin_ns_v7_cc` for W+ (u db, db u, c sb), W- (d ub, ub d) and a charge-violating pair. All agree to 4e-15; a 1e-3 change
of one kernel or the lepton on the wrong leg gives 2e-4 and 0.2. The assembly in `mod_xsects_nnlo_s_ns.f90` (tables
A-D, per-leg `Q_IS`, `ilept`, `mu2ref = mll^2` for the hoppet logs, unboosted `EC` in `kinematics_lo_z`) was read against
`finW`. Not covered here: the hoppet convolution itself and the z integration of the plus distributions (same code as Z).

## Limits of all RR and ONLO sectors with `-checklim` (2026-11-24)

`check_limits_rr` and `check_limits_onlo` (`f90/src/Checks/mod_check_int.f90`, PLAN C1.3) are called from `vegas_integrate`
with the sector itself. For a random point that passes the cuts, every limit of the sector (RR: S5, S6, X = TC or C5,
Y = C, C5, C6 or C56 and combinations; ONLO at fixed z: S, C1, C2, SC1, SC2) is set to 1e-1..1e-10 and H and all
counterterms are printed side by side, with the sign of H. Below: W+, one point per sector, largest relative difference
between the terms that must cancel, at depth 1e-3 and 1e-6 (soft: falls as x, collinear: as sqrt x).
Command, from `f90/`: `yes '' | ./mixed_ll -boson wp -corr nnlo -ch gq -sec rr_5161c -checklim | grep -av 'overriding input\|limit depth'`

```
ag_rr_5161           X: 4e-02 -> 2e-03  Y: 3e-03 -> 7e-05
ag_rr_5162           X: 9e-02 -> 3e-03  Y: 3e-01 -> 1e-03  X+Y: 9e-02 -> 3e-03
ag_rr_5261           X: 3e-02 -> 9e-04  Y: 4e-02 -> 1e-03  X+Y: 6e-04 -> 5e-05
ag_rr_5262a          Y: 7e-03 -> 5e-04
ag_rr_5262b          only the gamma gamma TC counterterm, zero for W
ag_rr_5262c          Y: 1e-02 -> 4e-04
ag_rr_5262d          only the gamma gamma TC counterterm, zero for W
aq_rr_5161a          S5: 6e-03 -> 6e-06  X: 8e-03 -> 3e-04  S5+X: 1e-02 -> 3e-04  Y: 1e-02 -> 5e-04  S5+Y: 1e-02 -> 5e-04  X+Y: 3e-02 -> 1e-03  S5+X+Y: 2e-02 -> 8e-04
aq_rr_5161b          S5: 9e-03 -> 9e-06  X: 3e-02 -> 9e-04  S5+X: 4e-02 -> 9e-04  Y: 2e-02 -> 2e-05  S5+Y: 2e-02 -> 2e-05  X+Y: 3e-02 -> 9e-04  S5+X+Y: 4e-02 -> 4e-03
aq_rr_5161c          S5: 4e-03 -> 4e-06  X: 1e-02 -> 5e-04  S5+X: 2e-02 -> 5e-04
aq_rr_5161d          S5: 9e-03 -> 9e-06  X: 3e-02 -> 9e-04  S5+X: 4e-02 -> 9e-04  Y: 2e-02 -> 2e-05  S5+Y: 2e-02 -> 2e-05  X+Y: 3e-02 -> 2e-03  S5+X+Y: 4e-02 -> 2e-03
aq_rr_5162           S5: 4e-03 -> 4e-06
aq_rr_5261           S5: 6e-03 -> 6e-06  X: 9e-03 -> 3e-04  S5+X: 8e-03 -> 3e-04  Y: 2e-02 -> 7e-04  S5+Y: 1e-02 -> 2e-04  X+Y: 4e-03 -> 2e-04  S5+X+Y: 5e-03 -> 7e-05
aq_rr_5262a          S5: 5e-03 -> 5e-06
aq_rr_5262b          S5: 7e-03 -> 7e-06  Y: 3e-03 -> 3e-06  S5+Y: 1e-02 -> 1e-05
aq_rr_5262c          S5: 1e-03 -> 1e-06  Y: 7e-03 -> 2e-04  S5+Y: 9e-03 -> 2e-04
aq_rr_5262d          S5: 7e-03 -> 7e-06  Y: 1e-03 -> 1e-06  S5+Y: 7e-03 -> 7e-06
ga_rr_5161a          Y: 2e-02 -> 6e-04
ga_rr_5161b          only the gamma gamma TC counterterm, zero for W
ga_rr_5161c          Y: 2e-03 -> 8e-05
ga_rr_5161d          only the gamma gamma TC counterterm, zero for W
ga_rr_5162           X: 2e-02 -> 7e-04  Y: 3e-02 -> 1e-03  X+Y: 1e-02 -> 5e-04
ga_rr_5261           X: 2e-01 -> 8e-03  Y: 2e-02 -> 5e-04  X+Y: 2e-01 -> 8e-03
ga_rr_5262           X: 1e-01 -> 4e-03  Y: 5e-03 -> 3e-05
gq_rr_5161a          S6: 5e-03 -> 5e-06  X: 1e-02 -> 3e-04  S6+X: 4e-03 -> 3e-04
gq_rr_5161b          S6: 7e-03 -> 8e-06  X: 2e-02 -> 5e-04  S6+X: 2e-02 -> 4e-04  Y: 2e-03 -> 2e-06  S6+Y: 7e-03 -> 8e-06  X+Y: 2e-02 -> 5e-04  S6+X+Y: 2e-02 -> 3e-04
gq_rr_5161c          S6: 8e-03 -> 8e-06  X: 5e-02 -> 2e-03  S6+X: 5e-02 -> 1e-03  Y: 1e-02 -> 5e-04  S6+Y: 2e-02 -> 6e-04  X+Y: 5e-02 -> 2e-03  S6+X+Y: 6e-02 -> 1e-03
gq_rr_5161d          S6: 1e-02 -> 1e-05  X: 8e-03 -> 3e-04  S6+X: 2e-02 -> 1e-04  Y: 5e-04 -> 5e-07  S6+Y: 1e-02 -> 1e-05  X+Y: 5e-03 -> 1e-04  S6+X+Y: 1e-02 -> 9e-05
gq_rr_5162           S6: 3e-03 -> 3e-06  X: 8e-02 -> 3e-03  S6+X: 5e-02 -> 2e-03  Y: 7e-02 -> 2e-03  S6+Y: 7e-02 -> 2e-03  X+Y: 7e-02 -> 2e-03  S6+X+Y: 7e-02 -> 2e-03
gq_rr_5163           S6: 2e-03 -> 2e-06  X: 4e-02 -> 1e-03  S6+X: 1e-02 -> 4e-04
gq_rr_5164           S6: 2e-03 -> 2e-06  X: 1e-01 -> 3e-03  S6+X: 6e-02 -> 2e-03  Y: 6e-03 -> 1e-04  S6+Y: 7e-03 -> 1e-04  X+Y: 1e-02 -> 2e-04  S6+X+Y: 1e-02 -> 2e-04
ns_ga_rr_5161a       S5: 3e-03 -> 3e-06  S6: 3e-03 -> 3e-06  S5+S6: 4e-03 -> 4e-06  X: 8e-02 -> 3e-03  S5+X: 8e-02 -> 3e-03  S6+X: 8e-02 -> 3e-03  S5+S6+X: 8e-02 -> 3e-03  Y: 6e-02 -> 2e-03  S5+Y: 6e-02 -> 2e-03  S6+Y: 6e-02 -> 2e-03  S5+S6+Y: 6e-02 -> 2e-03  X+Y: 3e-02 -> 8e-04  S5+X+Y: 2e-02 -> 5e-04  S6+X+Y: 3e-02 -> 7e-04  S5+S6+X+Y: 2e-02 -> 4e-04
ns_ga_rr_5161c       S5: 3e-03 -> 3e-06  S6: 3e-03 -> 3e-06  S5+S6: 4e-03 -> 4e-06  X: 9e-02 -> 3e-03  S5+X: 9e-02 -> 3e-03  S6+X: 9e-02 -> 3e-03  S5+S6+X: 9e-02 -> 3e-03  Y: 5e-02 -> 1e-03  S5+Y: 4e-02 -> 1e-03  S6+Y: 5e-02 -> 2e-03  S5+S6+Y: 4e-02 -> 1e-03  X+Y: 7e-02 -> 2e-03  S5+X+Y: 7e-02 -> 2e-03  S6+X+Y: 7e-02 -> 2e-03  S5+S6+X+Y: 7e-02 -> 2e-03
ns_ga_rr_5162        S5: 3e-03 -> 3e-06  S6: 3e-03 -> 3e-06  S5+S6: 5e-03 -> 5e-06  X: 2e-02 -> 6e-04  S5+X: 1e-03 -> 1e-06  S6+X: 1e-02 -> 4e-04  S5+S6+X: 3e-03 -> 3e-06  Y: 7e-02 -> 2e-03  S5+Y: 6e-02 -> 2e-03  S6+Y: 7e-02 -> 2e-03  S5+S6+Y: 6e-02 -> 2e-03  X+Y: 8e-02 -> 2e-03  S5+X+Y: 6e-02 -> 2e-03  S6+X+Y: 8e-02 -> 2e-03  S5+S6+X+Y: 6e-02 -> 2e-03
ns_ga_rr_5163        S5: 4e-03 -> 4e-06  S6: 1e-03 -> 1e-06  S5+S6: 4e-03 -> 4e-06  X: 2e-02 -> 6e-04  S5+X: 9e-04 -> 8e-07  S6+X: 2e-02 -> 5e-04  S5+S6+X: 1e-03 -> 1e-06
ns_ga_rr_5164        S5: 2e-03 -> 2e-06  S6: 2e-03 -> 2e-06  S5+S6: 3e-03 -> 3e-06  X: 6e-03 -> 2e-04  S5+X: 2e-03 -> 2e-06  S6+X: 6e-03 -> 2e-04  S5+S6+X: 3e-03 -> 3e-06  Y: 3e-03 -> 1e-04  S5+Y: 5e-03 -> 1e-04  S6+Y: 3e-03 -> 1e-04  S5+S6+Y: 4e-03 -> 1e-04  X+Y: 9e-03 -> 1e-04  S5+X+Y: 5e-03 -> 1e-04  S6+X+Y: 8e-03 -> 1e-04  S5+S6+X+Y: 5e-03 -> 1e-04
ns_ga_rr_5261        S5: 3e-03 -> 3e-06  S6: 3e-03 -> 3e-06  S5+S6: 4e-03 -> 4e-06  X: 1e-02 -> 5e-04  S5+X: 2e-03 -> 2e-06  S6+X: 1e-02 -> 4e-04  S5+S6+X: 3e-03 -> 3e-06  Y: 5e-02 -> 1e-03  S5+Y: 4e-02 -> 1e-03  S6+Y: 5e-02 -> 1e-03  S5+S6+Y: 4e-02 -> 1e-03  X+Y: 6e-02 -> 2e-03  S5+X+Y: 4e-02 -> 1e-03  S6+X+Y: 6e-02 -> 2e-03  S5+S6+X+Y: 4e-02 -> 1e-03
ns_ga_rr_5262a       S5: 3e-03 -> 3e-06  S6: 3e-03 -> 3e-06  S5+S6: 3e-03 -> 3e-06  X: 1e-01 -> 4e-03  S5+X: 1e-01 -> 4e-03  S6+X: 1e-01 -> 4e-03  S5+S6+X: 1e-01 -> 4e-03  Y: 8e-02 -> 2e-03  S5+Y: 8e-02 -> 2e-03  S6+Y: 8e-02 -> 2e-03  S5+S6+Y: 8e-02 -> 2e-03  X+Y: 3e-02 -> 7e-04  S5+X+Y: 2e-02 -> 5e-04  S6+X+Y: 3e-02 -> 6e-04  S5+S6+X+Y: 2e-02 -> 4e-04
ns_ga_rr_5262c       S5: 3e-03 -> 3e-06  S6: 3e-03 -> 3e-06  S5+S6: 4e-03 -> 4e-06  X: 1e-01 -> 4e-03  S5+X: 1e-01 -> 4e-03  S6+X: 1e-01 -> 4e-03  S5+S6+X: 1e-01 -> 4e-03  Y: 5e-02 -> 2e-03  S5+Y: 4e-02 -> 1e-03  S6+Y: 5e-02 -> 2e-03  S5+S6+Y: 4e-02 -> 1e-03  X+Y: 1e-01 -> 3e-03  S5+X+Y: 1e-01 -> 3e-03  S6+X+Y: 1e-01 -> 3e-03  S5+S6+X+Y: 1e-01 -> 3e-03
ns_ga_rr_5263        S5: 3e-03 -> 3e-06  S6: 1e-03 -> 1e-06  S5+S6: 3e-03 -> 3e-06  X: 2e-02 -> 5e-04  S5+X: 2e-03 -> 2e-06  S6+X: 1e-02 -> 4e-04  S5+S6+X: 3e-03 -> 3e-06
ns_ga_rr_5264        S5: 2e-03 -> 2e-06  S6: 2e-03 -> 2e-06  S5+S6: 3e-03 -> 3e-06  X: 6e-03 -> 2e-04  S5+X: 1e-03 -> 1e-06  S6+X: 6e-03 -> 2e-04  S5+S6+X: 2e-03 -> 2e-06  Y: 5e-03 -> 1e-04  S5+Y: 5e-03 -> 1e-04  S6+Y: 3e-03 -> 1e-04  S5+S6+Y: 5e-03 -> 1e-04  X+Y: 8e-03 -> 1e-04  S5+X+Y: 4e-03 -> 1e-04  S6+X+Y: 7e-03 -> 1e-04  S5+S6+X+Y: 4e-03 -> 1e-04
qa_rr_5161a          S5: 5e-03 -> 5e-06
qa_rr_5161b          S5: 6e-03 -> 6e-06  Y: 2e-03 -> 2e-06  S5+Y: 5e-03 -> 5e-06
qa_rr_5161c          S5: 1e-03 -> 1e-06  Y: 8e-03 -> 3e-04  S5+Y: 9e-03 -> 3e-04
qa_rr_5161d          S5: 5e-03 -> 5e-06  Y: 4e-03 -> 4e-06  S5+Y: 7e-03 -> 7e-06
qa_rr_5162           S5: 6e-03 -> 6e-06  X: 2e-02 -> 7e-04  S5+X: 9e-03 -> 3e-04  Y: 2e-02 -> 1e-03  S5+Y: 1e-02 -> 6e-04  X+Y: 7e-03 -> 5e-04  S5+X+Y: 3e-03 -> 5e-04
qa_rr_5261           S5: 2e-03 -> 2e-06
qa_rr_5262a          S5: 5e-03 -> 5e-06  X: 5e-02 -> 2e-03  S5+X: 6e-02 -> 2e-03  Y: 9e-02 -> 3e-03  S5+Y: 9e-02 -> 3e-03  X+Y: 4e-02 -> 1e-03  S5+X+Y: 2e-02 -> 8e-04
qa_rr_5262b          S5: 8e-03 -> 8e-06  X: 2e-01 -> 6e-03  S5+X: 2e-01 -> 7e-03  Y: 3e-03 -> 3e-06  S5+Y: 8e-03 -> 8e-06  X+Y: 2e-01 -> 6e-03  S5+X+Y: 2e-01 -> 2e-03
qa_rr_5262c          S5: 3e-03 -> 3e-06  X: 9e-02 -> 3e-03  S5+X: 9e-02 -> 3e-03
qa_rr_5262d          S5: 8e-03 -> 8e-06  X: 2e-01 -> 6e-03  S5+X: 2e-01 -> 7e-03  Y: 9e-04 -> 9e-07  S5+Y: 8e-03 -> 8e-06  X+Y: 2e-01 -> 7e-03  S5+X+Y: 2e-01 -> 3e-03
qg_rr_5261           S6: 4e-03 -> 4e-06  X: 1e-01 -> 4e-03  S6+X: 1e-01 -> 3e-03  Y: 4e-02 -> 1e-03  S6+Y: 4e-02 -> 1e-03  X+Y: 5e-02 -> 1e-03  S6+X+Y: 5e-02 -> 1e-03
qg_rr_5262a          S6: 6e-03 -> 6e-06  X: 3e-02 -> 9e-04  S6+X: 2e-02 -> 9e-04
qg_rr_5262b          S6: 8e-03 -> 8e-06  X: 2e-02 -> 5e-04  S6+X: 2e-02 -> 5e-04  Y: 2e-03 -> 2e-06  S6+Y: 8e-03 -> 8e-06  X+Y: 1e-02 -> 4e-04  S6+X+Y: 2e-02 -> 3e-04
qg_rr_5262c          S6: 9e-03 -> 9e-06  X: 7e-02 -> 2e-03  S6+X: 8e-02 -> 2e-03  Y: 2e-02 -> 7e-04  S6+Y: 2e-02 -> 7e-04  X+Y: 8e-02 -> 2e-03  S6+X+Y: 9e-02 -> 2e-03
qg_rr_5262d          S6: 1e-02 -> 1e-05  X: 2e-02 -> 5e-04  S6+X: 2e-02 -> 1e-04  Y: 4e-03 -> 4e-06  S6+Y: 1e-02 -> 1e-05  X+Y: 1e-02 -> 4e-04  S6+X+Y: 1e-02 -> 5e-05
qg_rr_5263           S6: 2e-03 -> 2e-06  X: 6e-03 -> 2e-04  S6+X: 4e-03 -> 1e-04
qg_rr_5264           S6: 3e-03 -> 3e-06  X: 5e-02 -> 2e-03  S6+X: 3e-02 -> 1e-03  Y: 6e-02 -> 2e-03  S6+Y: 6e-02 -> 2e-03  X+Y: 9e-03 -> 1e-04  S6+X+Y: 8e-03 -> 6e-05
ag_s_oewk            C1: 5e-02 -> 1e-03
ag_s_oqcd            C2: 6e-03 -> 3e-04
aq_s_oewk            C1: 4e-02 -> 2e-03
aq_s_oqcd            S: 1e+00 -> 1e-03  C1: 1e-02 -> 4e-04  C2: 6e-03 -> 2e-04  SC1: 2e-03 -> 2e-06  SC2: 1e-03 -> 1e-06
ga_s_oewk            C2: 1e-02 -> 2e-04
ga_s_oqcd            C1: 1e-02 -> 4e-04
gq_s_oewk_fs_53      S: 6e-03 -> 6e-06
gq_s_oewk_fs_54      S: 4e-03 -> 4e-06  C1: 8e-02 -> 1e-03  SC1: 5e-02 -> 9e-05
gq_s_oewk_is         S: 4e-03 -> 4e-06  C1: 4e-02 -> 1e-03  C2: 2e-02 -> 4e-04  SC1: 4e-02 -> 1e-03  SC2: 2e-02 -> 4e-04
gq_s_oqcd            C1: 1e-01 -> 3e-03
ns_s_oewk_fs_53      S: 6e-05 -> 6e-08
ns_s_oewk_fs_54      S: 1e-03 -> 1e-06  C1: 6e-03 -> 2e-04  SC1: 7e-03 -> 2e-04
ns_s_oewk_is_ra      S: 2e-03 -> 2e-06  C1: 5e-01 -> 5e-03  C2: 3e-01 -> 9e-03  SC1: 4e-02 -> 2e-03  SC2: 3e-01 -> 9e-03
ns_s_oqcd            C1: 2e-01 -> 2e-03  C2: 2e-01 -> 3e-03
qa_s_oewk            C2: 5e-01 -> 1e-02
qa_s_oqcd            S: 5e-03 -> 5e-06  C1: 1e-02 -> 3e-04  C2: 5e-02 -> 1e-04  SC1: 1e-03 -> 1e-06  SC2: 4e-03 -> 4e-06
qg_s_oewk_fs_53      S: 2e-03 -> 2e-06
qg_s_oewk_fs_54      S: 1e-03 -> 1e-06  C1: 4e-04 -> 2e-04  SC1: 2e-03 -> 2e-04
qg_s_oewk_is         S: 4e-03 -> 4e-06  C1: 4e-02 -> 1e-03  C2: 2e-02 -> 4e-04  SC1: 4e-02 -> 1e-03  SC2: 2e-02 -> 4e-04
qg_s_oqcd            C2: 2e-01 -> 3e-03
```

Four-quark (`xsect_checklim_nnlo4q`, W+ `ns_qqb rr_5161a`, depth 1e-6): |H+CT|/H = 5e-4 to 9e-4 on three points.
Not run: W- and Z (the binary has the W analysis compiled in), many points per sector (PLAN C1.8).

## Non-q qb' channels for W: `W_NNLO_polecheck_cleaned.nb` against the Fortran (2026-11-24)

`math/subtraction-notebooks/W_NNLO_polecheck_cleaned.nb` (Christian): W subtraction with general charges for gq, photon-quark,
photon-gluon and the four-quark channels, charged lepton on leg 4. `checks/run_verify_W_polecheck.sh` evaluates it headless
(`run_W_NNLO_polecheck.wls`, the `aux/` files replace the Dropbox paths, the q qb' part is skipped) and checks that the sum of the
Fortran sectors equals the notebook total: every structure of the finite remainder is matched to the Fortran piece that fills it
(a coverage line per channel confirms nothing is left over), general charges, five values of z, 40 digits.
Logs in `checks/logs/verify_W_polecheck_*_fortran.log`.

| Channel | Notebook structure -> Fortran sector, notebook - Fortran = 0 |
|---|---|
| gq | FLM[z p1] -> `s_12` (`calG_gq_cc` + `TCgq`, W+ and W-); FLM[z1 p1, z2 p2] -> `s_12` double boost; FLVew[z p1] -> `s_vewk`; ONLO FLM[pq] (+ log eta15 terms) -> `s_oqcd` elastic (W+, W-); ONLO FLM[p1, z p2, pq] -> `s_oqcd` boosted; ONLO FLM[z p1, p2, pa] -> `s_oewk_is` + `s_oewk_fs_53/54`; ONLO[FLVQEDfin] -> `rvewk_is` |
| photon-quark | `s_12` single and double boost; `s_vqcd`; `s_oqcd`; `s_oewk` boosted and elastic (`calG_CF_hard`); ONLO FLVQCDfin -> `rvqcd_is` |
| photon-gluon | both flavour orderings: `s_12` double boost, `s_oqcd`, `s_oewk` |
| four-quark | q qb' -> qb' qb' (two copies): `mySub_ns_qq`; q qb -> ...: `mySub_ns_qqb`; ratio 1 to 30 digits, log(mu) part included |

Real-virtual terms (`verify_W_polecheck_rv_fortran.wls`, from the notebook's own FVQEDgq, FLRVQCDγq, PqgRV, PqgNLO):
- Hard term: the Fortran finite remainder, OpenLoops eps^0 minus `get_fin_ewk_qqbllb_gen(5,2,3,4,...)` / `get_fin_qcd(2,5)`, equals
  the notebook definition (eps^0 of its pole term, alpha_B -> 1) for W+ and W-, both flavour orderings, and photon-quark; the
  four-point Born obeys the same relation. Legs, charge signs and the timelike pi^2 (pair 54 or 53) are what the notebook has.
- Counterterm: unintegrating the collinear terms of RVgq/RVγq gives the collinear limit of FLVQEDfin / FLVQCDfin implied by the
  notebook; its poles cancel and its finite part equals `col_rvewk_is_qg_gen` (W+: Q4 term, W-: Q3 term) and `col_rvqcd_is_qa`
  exactly. Fortran z = E5/E1 = 1 - z of the notebook; the notebook frame has E1 = E2 = EC. 5/3 pi^2 -> 3/2 pi^2 gives 0.73.

Compiled code (`verify_W_polecheck_compiled_fortran.wls`): the transcriptions used above (`W_polecheck_fortran_defs.wl`) agree with
`fortran_polecheck_driver.f90`, linked against the objects of `mixed_ll`: `calG_gq_cc`, `calG_ONLOQCD_gq` (W+, W-, both orderings,
two lmu), `calG_CF`, `col_rvewk_is_qg_gen`, `col_rvqcd_is_qa`, `Pqq_NLO`, `PqqAP_0_R`, `Pgq_spav` to 1e-16, HOPPET kernels to 2e-14.
Only the `get_fin` call arguments stay transcribed (private, next to the OpenLoops calls).

- Sector weights: the Fortran damps (`partition_nnlo_fact`, `partition_nlo_qed`, `calG_CF` damp1/damp2) are partitions of unity
  with the limits the notebook uses for omega tilde (kept symbolic there). Photon-quark elastic: the weight of
  log(eta52/(1-eta52)) is 1 - omega[{51,61},{5,6}] (the other choice differs by 0.2).
- Four-quark: the notebook charge is that of the incoming quark before the splitting (Qqpp); the Fortran squares the charge of
  the Born leg after it. They are equal in these channels (identical-flavour interference), as the chapter titles state.
- Sensitivity: dropping the -4 zeta2 Q1^2 term of `calG_gq_cc` gives 6.2; 13/2 -> 6 in `calG_CF_hard` gives 2/3; the wrong
  four-quark kernel gives O(1) to O(10^11).
- Not covered here: the mirrored channels qg, qa, ga (separate code; mirror tests above), the local counterterms inside ONLO
  (`-checkint` scans above), RR sectors, and the poles (the notebook has its own pole-check sections).
- `s_vewknf` and `rvewknf_is` stop for W ("not implemented") but are still listed in `sectors_list_wp.dat`/`_wm.dat`.

## W-aware limit checks (PLAN C1.4, 2026-11-24)

`check_limits_rr` and `check_limits_onlo` also test every limit in which a sector has no counterterm of its own (for W:
photon collinear to the neutrino, the gamma gamma -> l l triple-collinear limit; the soft gluon of `ns s_oqcd`, covered by
C1 + C2): x^k |H + counterterms| (k = variables at the edge) must fall between depths 1e-2 and 1e-8. One point per sector:

| Boson | Sectors | Limits with counterterm | Without: integrable | cut | SINGULAR |
|---|---|---|---|---|---|
| W+ | 78 (58 RR, 20 ONLO) | 440 | 620 | 4 | 1 |
| W- | 78 | 433 | 613 | 7 | 0 |

The W+ flag (`aq rr_5161b`, S6+X+Y) is double precision with three variables at the edge: x^3 |H + CT| falls to 1e-4 and
then rises; integrable on 7 further points (one cut). The first version tested H alone and flagged `ns s_oqcd` soft:
there the soft gluon is cancelled by C1 + C2 (1/(rho(1-rho)) = 1/rho + 1/(1-rho)), so H + counterterms is the right test.

## C1 limit of `qg s_oewk_fs_54` against the gq mirror (PLAN C1.9, 2026-11-24)

The C1.3 table (one point) showed 4e-4 -> 2e-4 between depths 1e-3 and 1e-6 for W+ `qg s_oewk_fs_54`. With 30 points
(28 series of a kinematic group with H and C1), median |H/C1 - 1| in the limit photon collinear to the charged lepton:

| Sector | 1e-3 | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 | 1e-10 |
|---|---|---|---|---|---|---|---|
| W+ `qg s_oewk_fs_54` | 1.0e-2 | 3.0e-3 | 8.9e-4 | 2.8e-4 | 9.0e-5 | 2.9e-5 | 2.9e-6 |
| W+ `gq s_oewk_fs_54` | 1.6e-2 | 6.5e-3 | 2.2e-3 | 7.1e-4 | 2.3e-4 | 7.2e-5 | 7.2e-6 |
| W- `qg s_oewk_fs_53` | 9.4e-3 | 2.7e-3 | 8.7e-4 | 2.7e-4 | 8.6e-5 | 2.7e-5 | 2.7e-6 |
| W- `gq s_oewk_fs_53` | 1.5e-2 | 4.8e-3 | 1.8e-3 | 5.9e-4 | 1.9e-4 | 6.0e-5 | 6.0e-6 |

All fall as sqrt x. In a few series H - C1 changes sign between 1e-3 and 1e-4 (qg W+: 3 of 28, gq W+: 2, qg W-: 0,
gq W-: 3); the point of the C1.3 table was one of them. Not a bug.

NLO real and RV sectors (16: `nloqcd`/`nloewk` `r_is`, `r_fs_53/54`, `rvewk_is`, `rvqcd_is`, `rvqcd_fs_53/54`) now fill
`limval_onlo` and `-checklim` sends them to `check_limits_onlo`. W+ and W-, one point per sector: 30 limits with counterterm,
all converge (soft as x, collinear as sqrt x); 50 limits without own counterterm, all integrable, among them the photon
collinear to the neutrino (`r_fs_53`, `rvqcd_fs_53` for W+; `_54` for W-) and the soft limit of `nloqcd r_is` (covered by C1 + C2).

## Intensive limit checks, W+ and W- (PLAN C1.5, C1.7, 2026-11-24)

`tests/checklim_run.sh wp|wm 200 csv` and `tests/checklim_analyse.py csv ref` (manual `f90/tests/CHECKLIM.md`):
94 sectors, 200 points each that pass the cuts, depths 1e-1..1e-10. W+: 1,866 limit/pair rows, all with at least 95% of
the points passing. W-: 1,881 rows, one below 95%: `qg rr_5262d` C56 integrability, one point, where the lepton dressing
is not collinear safe (PLAN team question 10). The remaining single-point flags are cut edges (the counterterm's Born
just below the mT cut, e.g. 199.9997 GeV, while H approaches it as sqrt x), double precision at triple and quadruple
limits, and sign crossings of A - B; the analysis does not count them (CHECKLIM.md).

C1.7 targets (median |A/B - 1| over pair rows and points; slope = median per-row Theil-Sen slope, expected 1/2 for
collinear limits; int = limits without own counterterm, graded integrable):

| Target | Sectors | Rows (pairs, int) | W+ min pass | W- min pass | W+ 1e-2 | 1e-4 | 1e-6 | 1e-8 | W- 1e-2 | 1e-4 | 1e-6 | 1e-8 | slopes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| RV6a C2 | ns `rvqcd_is` | 2, 0 | 1.000 | 0.989 | 3e-1 | 2e-2 | 2e-3 | 2e-4 | 3e-1 | 2e-2 | 2e-3 | 2e-4 | 0.50-0.53 |
| RV6b SC2 | ns `rvqcd_is` | 1, 0 | 1.000 | 1.000 | 2e-1 | 2e-2 | 2e-3 | 2e-4 | 2e-1 | 2e-2 | 2e-3 | 2e-4 | 0.50 |
| RV7a C1 | ns `rvqcd_fs_53/54` | 2, 1 | 1.000 | 1.000 | 5e-2 | 2e-3 | 2e-4 | 2e-5 | 5e-2 | 2e-3 | 2e-4 | 2e-5 | 0.52-0.58 |
| RV7b SC1 | ns `rvqcd_fs_53/54` | 1, 1 | 1.000 | 1.000 | 6e-2 | 2e-3 | 2e-4 | 2e-5 | 6e-2 | 2e-3 | 2e-4 | 2e-5 | 0.59-0.60 |
| RV3 C1, C2 | ns, gq, qg `rvewk_is` | 4, 2 | 0.990 | 0.995 | 5e-2 | 7e-3 | 8e-4 | 9e-5 | 4e-2 | 3e-3 | 5e-4 | 6e-5 | 0.46-0.54 |
| N2 | `nloewk` aq, qa `r_is` | 2, 8 | 1.000 | 1.000 | 1e-1 | 8e-3 | 9e-4 | 9e-5 | 2e-1 | 2e-2 | 2e-3 | 2e-4 | 0.49-0.50 |
| RR5 | `ns_ga`, gq, qg `rr_5163/5164/5263/5264` | 216, 56 | 1.000 | 1.000 | 3e-2 | 2e-3 | 2e-4 | 2e-5 | 3e-2 | 2e-3 | 2e-4 | 2e-5 | 0.5-1.0 |
| RR8 (TC, int only) | aq, qa, ag, ga RR | 0, 167-173 | 0.984 | 0.977 | | | | | | | | | |

In RV7a/b the neutrino sector (`fs_53` for W+, `fs_54` for W-) has no C1 counterterm and is integrable; the charged-lepton
sector converges as sqrt x. In RR5 the pairs S5 = S5C5 in the C5 limit agree to 2e-16 at every depth (slope 0, passing).
RR8: the lowest rows (97.7-99.5%) are TC with three or four variables at the edge, where double precision runs out.
