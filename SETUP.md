# SETUP — paper runs for W+ and W- (mixed QCD x EW, off-shell)

Common input for the production runs of `mixed_ll` and for the final POWHEG distribution
benchmarks (NLO QCD, NLO EW, NLO QCD l nu + photon, NLO EW l nu + jet). Every number below must be
identical in both codes. Status of the code: [PLAN.md](PLAN.md); checks: [CHECKS.md](CHECKS.md).

**Status: updated 2026-11-25** (energy, PDF, scale, cuts and benchmark objects decided). No open decision left (section 6).

Sources:
- [Z] Buccioni et al., *Mixed QCD-electroweak corrections to dilepton production at the LHC in the high invariant mass region*,
  [arXiv:2203.11237](https://arxiv.org/abs/2203.11237), Sec. 4. Numbers below are taken from the LaTeX source (`pheno.tex`).
- [ATLAS] ATLAS, *Measurement of double-differential charged-current Drell–Yan cross-sections at high transverse masses in pp
  collisions at sqrt(s) = 13 TeV*, [arXiv:2502.21088](https://arxiv.org/abs/2502.21088), JHEP 07 (2025) 026, Sec. 6.1, Table 1.

Two setups:
- **A — inclusive**: no cuts on the leptons apart from a lower cut on the l nu transverse mass.
- **B — ATLAS**: the fiducial region of [ATLAS] (high mT, W -> l nu, charges separately).

## 1. Common parameters (A and B)

The EW input follows [Z] and is the `mixed_ll` default. Energy and PDF differ from [Z] and must be passed explicitly.

| Parameter | Value | `mixed_ll` option (default) | POWHEG (`NSC-*`) |
|---|---|---|---|
| Collider | pp, sqrt(s) = **13 TeV** ([ATLAS]; [Z] used 13.6 TeV) | `-ecoll 13000` (default 14000!) | `ebeam1 = ebeam2 = 6500` |
| Boson | W+ -> nu e+, W- -> e- nubar; one lepton flavour, massless leptons | `-boson wp` / `-boson wm` | `chooseDY 1` / `-1` |
| PDF | `NNPDF40_nnlo_as_01180_qed`, member 0, LHAPDF **336100**, for all orders (LO included); [Z] used `NNPDF31_nnlo_as_0118_luxqed` | `-pdfset NNPDF40_nnlo_as_01180_qed -pdfmem 0` (not the default, see 5.7) | `lhans1 = lhans2 = 336100` |
| alpha_s | from the PDF set, alpha_s(mZ) = 0.118 | – | from the PDF |
| Flavours | 5 light flavours in the PDF, mb = 0; b-initiated W channels excluded (unit CKM) | `-mb 0` (default); `is_W_quark` | b excluded in the flavour lists |
| CKM | unit matrix | – | unit |
| EW scheme | G_mu, complex-mass scheme (with its O(alpha alpha_s) extension) | `-ew_scheme Gmu -cm_scheme true` (default) | `ewscheme 1`, complex masses |
| G_F | 1.16639e-5 GeV^-2 | `-Gf` (default) | same |
| mZ, Gamma_Z | 91.1876, 2.4952 GeV | `-mz`, `-gaz` (default) | same |
| mW, Gamma_W | 80.398, 2.1054 GeV | `-mw`, `-gaw` (default) | same |
| mH, Gamma_H | 125 GeV, 4.165 MeV | `-mh`, `-gah` (default) | same |
| mt, Gamma_t | 173.2 GeV, 0 | `-mt` (default) | `Gamma_t = 0` |
| alpha | derived: 1/132.277 (= 7.5598868e-3) | derived | fixed to 7.5598868e-3 in OpenLoops |
| Units | fb | `-units fb` | – |

### Leptons and photons

| | Value | Option |
|---|---|---|
| Charged lepton | **dressed**: photons within R(l, gamma) < 0.1 added to the lepton (E scheme, four-momentum sum); the neutrino is never dressed | `-R_phot 0.1` (default) |
| Distance | R = sqrt(dy^2 + dphi^2) (rapidity) | – |
| Beam photons | no special treatment | `-phot_beam false` (default) |
| All cuts and observables | in terms of the dressed lepton and the neutrino | – |
| Transverse mass | mT = sqrt(2 pT_l pT_nu (1 - cos dphi(l, nu))) | – |

[ATLAS] unfolds to Born level and provides `C_fsr` = dressed/Born (Delta R < 0.1) in HEPData: we compare dressed with dressed.

### Scales

| | Value | Option |
|---|---|---|
| Central | mu_R = mu_F = **mT / 2**, mT of the dressed lepton and the neutrino | `-dynscale true` (already mT/2 in the W analyses, see 5.2) |
| Variation | three points: mu_R = mu_F = mu x {1/2, 1, 2}, envelope (as [Z]) | `-xmuR`, `-xmuF` |
| Input-scheme variation | alpha(mZ) scheme with alpha(mZ) = 1/128, combined with the scale variation (as [Z]) | `-ew_scheme amz -alpha_mz 0.0078125` |

### Orders and channels

sigma = sigma^(0,0) + d sigma^(1,0) + d sigma^(0,1) + d sigma^(1,1), with (i,j) the order alpha_s^i alpha^j relative to LO.

| Order | Channels for W | Note |
|---|---|---|
| LO (0,0) | q qb' | no gamma gamma channel for W |
| NLO QCD (1,0) | q qb', g q, q g | |
| NLO EW (0,1) | q qb', gamma q, q gamma | photon-induced included |
| mixed (1,1) | q qb', g q, q g, gamma q, q gamma, gamma g, g gamma, four-quark (q qb', q q') | VV (two-loop) not yet in the code for W (PLAN.md section 3): the exact one will be added, no approximation in its place |
| not included | NNLO QCD (2,0), nf pieces | – |

### Technical settings (production)

| | Value | Why |
|---|---|---|
| Mass rescaling for resolved photons | `-qmin_rsc 1` | PLAN I17: the default 20 biases the IS photon sectors low |
| Statistics | >= 1.5M (grid) / 5M (production) points for sectors with a resolved photon; >= 4 independent seeds per sector, combined | PLAN I17: lower statistics are 1-3% low with underestimated errors |
| Seeds | `-seed` distinct per run | I18 |
| Cut-off parameters | `buff*`, damping exponents at the defaults; independence still to check (C4.4b) | – |
| Sectors | `sectors_list_wp.dat`, `sectors_list_wm.dat` | remove `s_oewk_is` draft first (S4) |

## 2. Setup A — inclusive

| Cut | Value |
|---|---|
| Transverse mass | mT > 200 GeV (dressed lepton) |
| Lepton pT, rapidity, missing pT, m_lnu | none |
| Generation cut | `-qmin 200 -qmax 13000` (m_lnu >= mT; to be confirmed, see 5.4) |
| Analysis | `analysis/mod_cut_histo_W_incl.f90` |

Transverse-mass windows for tables (as [Z] eq. 4.5, with mT in place of the mass):

| Window | mT |
|---|---|
| Phi(1) | 200 – 300 GeV |
| Phi(2) | 300 – 500 GeV |
| Phi(3) | 500 – 1500 GeV |
| Phi(4) | > 1500 GeV |

Distributions:

| Histogram | Range, bin |
|---|---|
| `rate` (total and per window) | – |
| `mlnu` | 200 – 1000 GeV, 20 GeV; 1000 – 3000 GeV, 100 GeV |
| `mT` | 200 – 1000 GeV, 20 GeV; 1000 – 3000 GeV, 100 GeV |
| `ptl` (charged lepton) | 0 – 1000 GeV, 20 GeV |
| `yl` (charged lepton) | -5 – 5, 0.2 |
| `ptlnu` | 0 – 500 GeV, 10 GeV |

## 3. Setup B — ATLAS fiducial ([ATLAS] Sec. 6.1)

| Cut | Value |
|---|---|
| Charged lepton | pT_l > 65 GeV, \|eta_l\| < 2.4 (no crack exclusion at particle level) |
| Neutrino | pT_nu > 85 GeV (no rapidity cut) |
| Transverse mass | 200 < mT < 5000 GeV (single-differential); 200 < mT < 2000 GeV (double-differential) |
| Jets, isolation | none |
| Lepton charge | part of the fiducial definition: W+ and W- separately |
| Generation cut | `-qmin 200` (m_lnu >= mT; to be confirmed, see 5.4) |
| Analysis | `analysis/mod_cut_histo_W_ATLASfid.f90` |

Distributions ([ATLAS] Sec. 6, Table 1, Sec. 8.4):

| Histogram | Binning |
|---|---|
| `rate` | – |
| d sigma / d mT | 12 bins in 200 – 5000 GeV; **edges from HEPData** (not in the paper text) — TODO |
| d^2 sigma / d mT d\|eta_l\| | mT in [200, 300], [300, 425]: \|eta\| edges 0, 0.2, ..., 2.4 (12 bins); mT in [425, 600], [600, 900]: 0, 0.4, ..., 2.4 (6 bins); mT in [900, 2000]: 0, 0.6, 1.2, 1.8, 2.4 (4 bins) |
| Charge asymmetry A_l = (d sigma+ - d sigma-) / (d sigma+ + d sigma-) | same bins as the double-differential (post-processing of the W+ and W- histograms) |
| Additional | `ptl`, `ptnu`, `etal`: 10 GeV / 0.1 bins (benchmark only) |

The ATLAS mT bins of the double-differential grouping are read off Table 1; confirm against HEPData.

## 4. POWHEG distribution benchmarks

Fixed-order POWHEG-BOX-RES (`testplots 1`), same input as section 1 and the same analysis (setup A or B).
Recipe and base input: `f90/tests/powheg/paper/`. `mixed_ll`: analyses `W_incl` (A) and `W_ATLASfid` (B) with
`-require none|photon|jet`. POWHEG: `pwhg_analysis-NSCpaper.f` with `nscsetup 1|2` (A, B) and `nscrequire 0|1|2`.

| Benchmark | POWHEG code | `mixed_ll` sectors | Analysis | Extra settings |
|---|---|---|---|---|
| NLO QCD | `NSC-lnuQCD` | `lo ns`, `nloqcd` ns/gq/qg | A / B | – |
| NLO EW | `NSC-lnuEW` | `lo ns`, `nloewk` ns/aq/qa | A / B | `detailedNLO 1` (else the photon has no id and the lepton is not dressed) |
| NLO QCD l nu + photon | `NSC-lnuaEW` | 39 sectors with a resolved photon (C3.6) | A / B + photon | resolved photon: pT_gamma > 50 GeV, \|y_gamma\| < 2.5, R(gamma, l) > 0.1, R(gamma, q) > 0.4; pT(l nu) > 30 GeV (as `V_lla_nloqcd`) |
| NLO EW l nu + jet | `NSC-lnujQCD` | 33 sectors (C3.7) | A / B + jet | jet = final quark (dressed) or gluon with pT_j > 30 GeV; photon never a jet (as `V_llj_nloew`) |
| four-quark (optional) | `NSC-lnujjMIX` | four-quark RR | A / B + 2 jets | as `V_lljj_4q` |

The photon and jet definitions of the validation runs (table above) are kept; no ATLAS-like objects.

## 5. Open points before production

1. **Scale in POWHEG.** mT/2 is added to `NSC-lnuQCD` and `NSC-lnuEW` (`f90/tests/powheg/paper/Born_phsp_mT_scale.patch`); `NSC-lnuaEW`, `NSC-lnujQCD` still have only mu = mZ.
2. **Generation cut.** With mT > 200 GeV the result must not depend on `-qmin`. First check at NLO: run `-qmin 100` and `-qmin 200` and compare (PLAN C4.8).
3. **mT bins of [ATLAS].** `W_ATLASfid` books mT in 10 GeV bins; the 12 ATLAS edges come from HEPData.

Done: `W_incl` and `W_ATLASfid` (no neutrino rapidity cut, pseudorapidity of the lepton, no product cut);
13 TeV and `NNPDF40_nnlo_as_01180_qed` are the defaults (install the set where it is missing).

## 6. Decisions (2026-11-25)

No open decision. 13 TeV for both setups; `NNPDF40_nnlo_as_01180_qed`; mu = mT/2 with the three-point variation;
setup A with mT > 200 GeV only, no resonance-region run; no approximate VV (the exact one will be added);
validation cuts kept for the POWHEG benchmarks.

## 7. Reference numbers

[Z] at 13.6 TeV, `NNPDF31_nnlo_as_0118_luxqed`, mu = m_ll/2, its fiducial cuts (m_ll > 200, pT_l > 30, sqrt(pT_l- pT_l+) > 35, \|y_l\| < 2.5), fb:

| | sigma^(0,0) | d sigma^(1,0) | d sigma^(0,1) | d sigma^(2,0) | d sigma^(1,1) |
|---|---|---|---|---|---|
| sum of channels | 1621.06 | 340.37 | -47.046 | 13.88 | -15.71 |

`mixed_ll` vs POWHEG at NLO with the current `W_ATLAS` analysis (13 TeV, fixed mu = mZ, neutrino \|y\| cut included):
`f90/tests/powheg/README.md` (C3.1, C3.2).
