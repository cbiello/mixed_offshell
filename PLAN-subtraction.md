# PLAN-subtraction — integrated subtraction ("S" sectors) for W+ and W-

Companion of [PLAN.md](PLAN.md) section 6. Checks that passed are recorded in [CHECKS.md](CHECKS.md).
Reference formulae: `math/subtraction-notebooks/Z.nb` (all channels, Z) and `W.nb` (q qb' channel, general
charges, poles cancelled). Scripts: `math/subtraction-notebooks/checks/`.

Legend: ✅ done and checked · 🟡 code exists, not validated · 🔧 broken · ❌ missing · ⏸ later

## 1. Conventions

Process q(1) qb'(2) -> l(3) l'(4) + X. Incoming charges `Q_IS(al)`, `Q_IS(be)` (antiquark: minus the quark charge).

| | Z | W+ | W- |
|---|---|---|---|
| final state (3,4) | e- e+ | nu e+ | e- nubar |
| charged lepton `ilept` | 3 and 4 | 4 | 3 |
| its outgoing charge `Qlept` | -1, +1 | +1 | -1 |
| charge conservation | Q1 + Q2 = 0 | Q1 + Q2 = +1 | Q1 + Q2 = -1 |

`W.nb` has one charged lepton with outgoing charge `Qe` (new notation `Qlp`) on leg 4. Nothing in the derivation
depends on the leg number, so **W- is the same set of functions with `Qlept = -1`, `ilept = 3`** (invariants and
angles of leg 3). The code functions already take `(Qlept, ilept)`.

Logs: `lmu = log(mu^2/4EC^2)` at fixed partonic energy; HOPPET tables use `log(mu^2/mll^2) = lmu - log z`.

## 2. Notebook pieces and where they live

| Notebook piece | ns | gq, qg | aq, qa | ag, ga |
|---|---|---|---|---|
| FLVqcd boosted + elastic | `s_vqcd` | – | `s_vqcd` | – |
| FLVew boosted + elastic | `s_vewk` | `s_vewk` | – | – |
| FLM double boost, single boost, elastic | `s_12` | `s_12` | `s_12` | `s_12` |
| ONLO-QCD: resolved gluon/quark x EW integrated CT | `s_oqcd` | `s_oqcd` | `s_oqcd` | `s_oqcd` |
| ONLO-EWK: resolved photon x QCD integrated CT | `s_oewk_is(_ra)`, `s_oewk_fs_53/54` | `s_oewk_is`, `s_oewk_fs_53/54` | `s_oewk` | `s_oewk` |
| four-quark | `s_qqb`, `s_qq` | – | – | – |

Kernels: `IntSub/mod_hoppet_nnlo.f90` (tables, with delta(1-z) endpoints that carry part of the elastic terms);
elastic functions: `IntSub/mod_int_sub_nnlo.f90`.

## 3. Validation methods

| ID | Method | Proves |
|---|---|---|
| M1 | Symbolic: Fortran transcribed to Mathematica against the finite remainder of the re-run notebook | formulae and bookkeeping of the elastic/boosted terms |
| M2 | `-checkint` scan | the integrand is integrable (local counterterms of the resolved emission) |
| M3 | CP test: W- in p p against W+ in pbar pbar (`-cptest`, PDFs q <-> qbar), same cuts | the W- wiring (lepton index, charge, ids), given W+ |
| M4 | NLO QCD l l' gamma: sum of the sectors with a resolved photon against POWHEG `NSC-lnuaEW` | everything in section 5 at once, independently of the notebooks |
| M5 | Z: old code against general code (`tests/run_regression.sh z`) | no regression |

## 4. Work packages

### 4a. q qb' (ns), W+

| ID | Task | Method | Status |
|---|---|---|---|
| A1 | `s_12`, `s_vqcd`, `s_vewk` against `W.nb` | M1 | ✅ 2026-11-12 |
| A2 | `s_oewk_is_ra`, `s_oewk_fs_53/54`: QCD kernels on the resolved-photon process against the ONLO-EWK terms of `W.nb` | M1 (kernels), M2, M4 | ✅ kernels, scan, POWHEG (2026-11-13) |
| A3 | `s_oqcd` against the ONLO-QCD terms of `W.nb` (`calG_ONLOQCD_ns` vs `genGnloQCD`, PLAN S3) | M1, M2 | 🟡 later, not needed for 5 |
| A4 | `s_oewk_is`: draft routine, use `_ra`; rewire and delete the draft (PLAN S4) | – | 🔧 |

### 4b. q qb' (ns), W-

| ID | Task | Method | Status |
|---|---|---|---|
| B1 | `s_12`, `s_vqcd`, `s_vewk`: one W branch for both signs, `Qlept`/`ilept` from the boson | M1 (same functions), M3, M5 | ✅ CP test (2026-11-13) |
| B2 | `s_oewk_is_ra`, `s_oewk_fs_53/54` for W- | M2, M3, M4 | ✅ scan, POWHEG |

### 4c. gq, qg

| ID | Task | Method | Status |
|---|---|---|---|
| C1 | `s_oewk_is`, `s_oewk_fs_53/54` for W+ and W-: kernels against the gq ONLO-EWK terms of `Z.nb` (no charges in them) | M1 (kernels), M2, M3, M4 | ✅ kernels, scan, POWHEG |
| C2 | `gq s_oewk_fs_53/54` is zero at every point, also for Z (PLAN S10b) | fix, then M5 against `qg` mirror | ✅ fixed: debug lines overwrote `xE`, `xRHO` |
| C3 | `s_12`, `s_vewk`, `s_oqcd` for W | – | ⏸ not needed for 5 |

### 4d. photon-induced (aq, qa, ag, ga)

Not needed for section 5 (no final-state photon). ⏸ Needs W formulae for these channels.

## 5. Target: NLO QCD corrections to l l' gamma

Sectors that contribute with a resolved photon (analysis `V_lla_nloqcd`):

| Piece | Sectors |
|---|---|
| LO | `nloewk ns` `r_is`, `r_fs_53`, `r_fs_54` (hard events only pass the photon cuts) |
| virtual | `nnlo ns` `rvqcd_is`, `rvqcd_fs_53`, `rvqcd_fs_54` |
| real | `nnlo ns_ga` `rr_*` (10), `gq` `rr_*` (7), `qg` `rr_*` (7) |
| integrated | `nnlo ns` `s_oewk_is_ra`, `s_oewk_fs_53/54`; `gq`, `qg` `s_oewk_is`, `s_oewk_fs_53/54` |

| ID | Task | Status |
|---|---|---|
| T1 | Matching analyses in both codes (photon cuts, pT(l l') cut) | ✅ `V_lla_nloqcd`, `pwhg_analysis-NSCa.f` |
| T2 | LO l l' gamma: mixed_ll against POWHEG, Z, W+, W- | ✅ agree at 0.03% (needs `-qmin_rsc 1`) |
| T3 | NLO QCD l l' gamma: mixed_ll against POWHEG, Z, W+, W- | ✅ W+ agrees at 0.1% (6.1928(21) vs 6.1937(75)); W-, Z: POWHEG rerun without damping in progress |

## 6. Order

A2 -> B1 -> M3 tool -> B2 -> C2 -> C1 -> T1..T3. A3, A4, C3, 4d afterwards.
