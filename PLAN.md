# PLAN — mixed QCD x EW NNLO for off-shell W+ and W-

Shared status board. Read [STRUCTURE.md](STRUCTURE.md) first for names and conventions.
Checks that passed, with numbers: [CHECKS.md](CHECKS.md).
Based on a read-through of branch `AI-era` (snapshot 2026-11-11). Nothing below was
established by running the NNLO sectors; items marked *suspected* come from reading the
code and must be confirmed before fixing.

**How to use**: pick a task, put your name and date in the Status column, tick it when the
listed acceptance checks pass for **both** W+ and W-. Keep the Z result unchanged
(old-code vs new-code build agreement, `OLDCODE=true/false`) at every step.

Legend: ✅ done and checked · 🟡 code exists, not validated · 🔧 partial / broken · ❌ missing · ⏸ out of scope for now

## 0. Summary

| Block | Z | W+ | W- |
|---|---|---|---|
| LO, NLO QCD, NLO EW (ns, gq, qg, aq, qa) | ✅ | ✅ POWHEG | ✅ POWHEG |
| VV (two-loop) | ✅ | ⏸ | ⏸ |
| RV, EW loop + gluon (`rvewk_is`: ns, gq, qg) | ✅ | ns 🟡 (limits ✅); gq, qg 🔧 integrand Z | same |
| RV, QCD loop + photon (`rvqcd_*` ns) | ✅ | ✅ limits | ✅ limits |
| RV, QCD loop + photon (aq, qa) | ✅ | ✅ limits | ✅ limits |
| RR `ns_ga` | ✅ | ✅ limits | ✅ limits |
| RR `gq`, `qg` | ✅ | ✅ limits | ✅ limits |
| RR `aq`, `qa`, `ag`, `ga` | ✅ | ✅ limits | ✅ limits |
| RR four-quark | ✅ | ❌ | ❌ |
| S ns: `s_12`, `s_vqcd`, `s_vewk` | ✅ | 🟡 | ❌ (`stop`) |
| S ns: `s_oqcd`, `s_oewk_*` | ✅ | 🟡 | 🟡 |
| S ns: `s_qqb`, `s_qq` | ✅ | ❌ | ❌ |
| S gq, qg: `s_oqcd`, `s_oewk_*` | ✅ | 🟡 | 🟡 |
| S gq, qg: `s_12`, `s_vewk` | ✅ | ❌ | ❌ |
| S aq, qa, ag, ga (all) | ✅ | ❌ | ❌ |
| Limit-scaling checks usable on all sectors | 🔧 | 🔧 | 🔧 |
| Fortran vs Mathematica integrated CT check | – | 🔧 (one term) | ❌ |

Important: sectors that are still Z-only have **no guard**. A W build runs them and
silently returns neutral-current numbers.

## 1. Infrastructure and hygiene (do first)

| ID | Task | Where | Status |
|---|---|---|---|
| I0 | Boson selected at run time: `-boson z/wm/wp` replaces the `-D_Vcharge` preprocessor flag (329 directives -> `if (Vcharge == ..)`), flag removed from makefile and `configure`, W output files tagged `_wm`/`_wp`. Verified: 51 sectors x (Z, W-, W+) reproduce the old builds digit for digit. Side fix: uninitialised `kin` in `xsect_nloewk_r_is_aq/_qa` | all of `src/`, `analysis/` | ✅ 2026-11-11 |
| I1 | Build fixes. Done in `f5196e2` (Chiara): `_Vcharge` passed by the makefile (now superseded by I0), stale `mod_amplitudes_loop_pplnu.o` removed, LO `allocate` guarded. Still local: `.o` dependencies in `get_dep.py` (needed for `make -j`), LHAPDF rpath in the makefile, `cfgcb.txt` | `f90/makefile`, `get_dep.py` | ✅ |
| I2 | README: `ppllj` added to the OpenLoops list, HOPPET and handyG listed, build (cfg keys) and run (`-boson`) sections | `README.md`, `f90/README` | ✅ 2026-11-11 |
| I3 | Guard on every Z-only integrand: `if (Vcharge /= 0) stop '<name>: not implemented for W'` right after the result is initialised (70 integrands: `aa`, `ag`, `ga`, `aq`, `qa` channels, four-quark, `rvewk`, `rvqcd_is_aq/qa`, `s_gq/s_qg` `s12`/`vewk`, `nf` pieces, VV). **Remove the line when the integrand is generalised.** Not guarded, because W work is in progress there: `rr_gq` (5), `rr_526k_qg`, `s_aq_oqcd` | `Sectors/` | ✅ 2026-11-11 |
| I4 | Debug code removed: `print`/`pause` in `rr_qg` (5262b/d, 526k), `s_gq` `oewk_fs`, `rvewk_is` (ns, gq, qg, with the unused `_gen` call); `print`/`stop` in `s_aq` `oqcd`; 'work in progress' print in `s_ns` `oewk_is`. Z runs again in these sectors. Seen afterwards (Z): `gq s_oewk_fs_53/54` crashes in Vegas `rebin` (I14); `ns s_oewk_is` (new branch unfinished) returns nonsense, use `s_oewk_is_ra` | `Sectors/` | ✅ 2026-11-11 |
| I5 | One global old/new switch: the ~40 local `oldcode` flags are replaced by the preprocessor flag `_oldcode` (`OLDCODE=true/false` in the cfg, default false; 206 blocks). Verified: new-code build reproduces 55 sectors x (Z, W-, W+) digit for digit, old-code build reproduces the all-`.true.` reference for Z. `aq s_oqcd`, formerly hard-wired to old code, now follows the flag (old and new agree to 9 digits) | `Sectors/`, `configure` | ✅ 2026-11-11 |
| I6 | Sector list synchronised with the dispatch in `Chooser/mod_run_proc_*.f90`: added `vvnf_full`, `rvewknf_is`, `s_vewknf` (ns, gq, qg) and `s_oewk_is_ra`; `sectors_list.dat` regenerated (`./mixed_ll -list sec`), so `ag rr_5161` and `ga rr_5262` are now right. Per-boson list: `./mixed_ll -list sec -boson wm` writes `sectors_list_wm.dat` (no `aa` channel for W); sectors not yet implemented for W stay in the list and stop with the I3 guard | `sectors_list.dat`, `Chooser/mod_sects_list.f90` | ✅ 2026-11-11 |
| I7 | Lepton ids: new `lep_ids(2)` (set by `-boson`); the 66 hard-wired `[id_el,-id_el]` in W-enabled integrands (all of `rr_gq`, `rr_526k_qg`, `rr_5i6ic_ns_ga`, `s_ns` `oqcd`/`oewk_is`, `s_qg`, `s_aq_oqcd`) now use it. Z unchanged; for W only `ns_ga rr_5262c`, `qg s_oewk_is`, `qg s_oewk_fs_53` move. Not done: collapsing the existing three-branch id blocks into one line | `Sectors/`, `Parms/mod_parms.f90` | ✅ 2026-11-11 |
| I7b | `%part` removed from `KinConfig` and from all integrands (419 lines); the W analyses (`W_ATLAS`, `W_ATLASNLOQCD`) identify the neutrino from `event%ids(3)` and stop if the ids are not set. `nloqcd` `v`/`s` integrands, which only set `%part`, now set `%ids`. Effect: Z identical; W- identical in all sectors that worked; W+ `nloewk aq/qa r_is` corrected (lepton and neutrino cuts were swapped): 1.465(8) / 1.462(6) fb; `rr_gq`, `rr_qg`, `aq s_oqcd` and W+ `ns s_12/s_vqcd/s_vewk` now return numbers instead of zero. `ns s_oewk_is` (unfinished, use `s_oewk_is_ra`) still has configurations without ids | `analysis/`, `Sectors/`, `GenProc/` | ✅ 2026-11-11 |
| I8 | Flavour selection in the tree amplitudes by integer tests: `is_W_pair(i,j)` and `is_W_quark(i)` in `mod_parms` replace the float comparison and the eight hand-written `cycle` blocks; `current_charge`, `Qgeneration` removed. b is now skipped for W in all single-quark amplitudes (was missing in `res_tree_ga_qg_gen`, `ga_aq`, `ga_qa`). 55 sectors x (W+, W-) identical. Open for the team: CKM beyond the unit matrix (question 1) | `Amplitudes/mod_amplitudes_tree_ppll.f90`, `Parms/` | ✅ 2026-11-11 |
| I9 | Dead W prototypes deleted: `xsect_lo_ns_wp/_wm`, `xsect_nloqcd_r_is_ns_wp/_wm`, `xsect_nloewk_r_is_ns_wp/_wm`, `qQpb_lumi_wp/_wm`, `get_coupl_w`, `mod_amplitudes_tree_ppnul.f90`, `mod_xsects_nnlo_rvqcd_original.f90`, `makefile_chiara`, two `*.bak`, commented chooser entries. Kept: `ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w` (Z four-quark channels with W exchange, now guarded for W) | – | ✅ 2026-11-11 |
| I10 | Run header prints the boson and the charges of legs 3,4 | `Parms/mod_proc_parms.f90` | ✅ 2026-11-11 |
| I11 | Regression suite: `tests/run_regression.sh z/wm/wp [update]` runs every sector of the list with the same seed at low statistics and compares with `tests/ref_*.txt` (relative 1e-7). References: Z with `dflt` (139 numbers), W+/W- with `W_ATLAS` (58 numbers, the rest stop as not implemented). About 15 s on 32 cores | `f90/tests/` | ✅ 2026-11-11 |
| I12 | Photons are not recombined with the neutrino (`recombine_la`, already there). Boson and analysis are checked at start: W stops unless the analysis is named `W_*`, `w_*` or `V_*`; `W_ATLAS` stops for Z. The analysis itself stays a compile-time choice | `Parms/mod_proc_parms.f90`, `analysis/` | ✅ 2026-11-11 |
| I13 | Build with `-finit-real=snan -ffpe-trap=invalid` and run all sectors (Z, W+, W-): all clean now. Fixed: 0/0 in the acceptance printout, PDF slot 7 divided before being set, unset slots of OpenLoops result arrays and of helicity arrays, a whole-array division in `rr_5261_qg`. None changed a number. With fewer than ~20k points Vegas itself can produce 0/0 in the grid refinement (statistics artefact) | `Sectors/`, `Amplitudes/`, `HoppetTools/` | ✅ 2026-11-11 |
| I14 | Vegas `rebin` stops with 'integrand is zero at every point' instead of running out of bounds | `Vegas/mod_my_vegas.f90` | ✅ 2026-11-11 |
| I15 | `V_minimal` fixed: the dynamic scale used a transverse mass that was never computed (scale = 0); now `mll/2`. For W the charged lepton and the neutrino are taken from `ids`. Debug print removed | `analysis/mod_cut_histo_V_minimal.f90` | ✅ 2026-11-11 |
| I16 | Z, old vs new code (same seed): all 55 tested sectors agree to at least 7 digits after I16a-c | – | ✅ 2026-11-11 |
| I16a | `nloewk aq s`, `qa s` (0.4% off): bug in the new code, Born multiplied by `Q2(fa)*Q2(fb)` (Q^4 for Z) instead of one Q^2 of the quark from the photon splitting. **Fixed**: loop replaced by `multiply_IS_charges_sq(res_lo_new,1)` (aq) and `(...,2)` (qa). Z new = old to all digits. W values change (W- aq: -3.7e-3 -> -2.3e-2 fb at test statistics) | `mod_xsects_nloewk_vs.f90` | ✅ 2026-11-11 |
| I16b | `nloewk qa r_is` (old = -19221 fb): **bug in the old branch**. In C1, `Pqq(z)/(1-z)` is applied inside the old branch and again after `#endif` (in `aq` the second one sits inside the new branch). With it once, old = new to 12 digits. New code is right | `mod_xsects_nloewk_r.f90` (`xsect_nloewk_r_is_qa`) | ✅ fixed (old = new, -9.06994348143 fb) |
| I16c | `ns_ga rr_5161a/5262a` (old = -58 fb): **old branch incomplete**, five defects: S6 never calls `res_tree_g_qqb`; S6 and S5S6 pass the unset `res_tmp_vect` instead of `res_tmp_vect_old` and leave the result in `respdf_vect` while the code reads `respdf_tmp`; TC, TCS6, TCS5 fill `respdf` instead of `respdf_tmp(:,1)`. With these repaired, old = new to 11 digits, so the **new code is validated for Z in sector a** | `mod_xsects_nnlo_rr_ns.f90` (`xsect_nnlo_rr_5i6ia_ns_ga`) | ✅ fixed (old = new to 11 digits) |
| I16d | `rvqcd_fs_53/54`, `rvqcd_is` (agree to 1e-7..1e-9): not a bug. Old code uses the analytic one-loop amplitude (`res_qcdloop_j`), new code uses OpenLoops; the difference is numerical precision | – | ✅ understood |
| I17 | Vegas reliability with a resolved photon: (a) the technical mass rescaling `-qmin_rsc 20` (default) biases the initial-state photon sectors low for hard-photon cuts (`nloewk ns r_is` W+: 1.495-1.517 fb over seeds); (b) with 300k/1M points many RR/S sectors are 1-3% low with underestimated errors. For now: `-qmin_rsc 1` and 1.5M/5M points. To do: check the production NNLO setup of the same sectors (`rr_5161/5262/5162/5261`, `rvqcd_is`, `s_oewk_is`) with several seeds, and decide the default | `Parms/mod_proc_parms.f90`, `Vegas/` | 🔧 found 2026-11-13 |

## 2. NLO (done; remaining loose ends)

| ID | Task | Where | Status |
|---|---|---|---|
| N1 | NLO QCD ns/gq/qg and NLO EW ns: implemented for W±, benchmarked against POWHEG (commit `346fdea`) | – | ✅ |
| N2 | NLO EW `aq`/`qa` `r_is` for W: limits checked with `-checklim` (2026-11-11). Z closes (ratio ~ sqrt(depth)); W does not, two bugs: (a) `qa` C2 uses the charge of leg 1, must be leg 2 (:1053); (b) `res_tree_a_aq_gen` / `res_tree_a_qa_gen` keep the b quark for W while the Born does not (missing `abs(.)==5` skip, see I8). Both fixed: all four W limits close, Z unchanged digit for digit. W+ aq/qa = 1.508(7)/1.494(6) fb, W- = 0.856(4)/0.847(3) fb (`W_ATLAS`) | `mod_xsects_nloewk_r.f90`, `mod_amplitudes_tree_ppll.f90` :484-550 | ✅ 2026-11-11 |
| N3 | NLO EW `s_aq`/`s_qa`: double charge factor, fixed, see I16a. Included in the POWHEG NLO EW comparison (C3.2) | `mod_xsects_nloewk_vs.f90` | ✅ |
| N4 | NLO QCD `gq`/`qg` integrated collinear: closed. The full NLO QCD (ns + gq + qg) agrees with POWHEG for W+ and W- at the 0.05% level (C3.1); gq + qg contribute 72 fb (W+), so they are tested to about 1% of their size | `mod_xsects_nloqcd_vs.f90` | ✅ 2026-11-11 |
| N5 | POWHEG comparison for W+ and W-: LO, NLO QCD (C3.1) and NLO EW (C3.2) agree | `f90/tests/powheg/` | ✅ 2026-11-12 |

## 3. Double virtual (out of scope for now)

| ID | Task | Status |
|---|---|---|
| V1 | Two-loop mixed amplitude for q qb' -> l nu, its interface, `vvfc/vvxf/vvnf` for W. `mod_xsects_nnlo_vv.f90` is Z only | ⏸ |
| V2 | `nf` pieces for W (`vvnf`, `rvewknf_is`, `s_vewknf`): decide whether they stay separate | ⏸ |

## 4. Real-virtual

### 4a. EW one-loop x real gluon

| ID | Task | Where | Status |
|---|---|---|---|
| RV1 | `rvewk_is` ns generalised for W± (2026-11-12): Hard, C1, C2 with `_gen` amplitudes, `get_respdf_gen`, `lep_ids`, per-leg charges (RV3); legacy under `_oldcode`; W guard removed. Z: new -0.71454867610 vs old -0.71454867562 fb (7e-10). W limits (`-checklim`, ~300 W_ATLAS points per leg): median \|H+CT\|/H 2e-3, 2e-4, 2e-5 at depth 1e-4, 1e-6, 1e-8 (W- and W+, C1 and C2), ~sqrt(depth); with the new charge logs switched off it stays at 9e-3. **Finite parts not validated**: matching with the integrated counterterm (W+ `s_12` term B has `-2 Qlept Q_IS log(s14/s12)`, same structure) still to check; W- integrated sectors missing (S2) | `mod_xsects_nnlo_rvewk.f90`, `mod_limits_rv.f90` | 🟡 limits ✅, finite ❌ |
| RV2 | `rvewk_is` gq, qg: same situation (`res_ewkloop_g_gq_gen`, `_g_qg_gen` ready and called as a test; integrand Z) | :172-380 | 🔧 |
| RV3 | One-loop collinear limits with general charges, ns: `col_rvewk_is_ns_gen`. The C1 (C2) counterterm is c_i^2 (splitting part) + 2 log(1-z) c_i [c3 log(s_i3/s12) + c4 log(s_i4/s12)], c = incoming charges (c3=-Q3, c4=-Q4), from the mismatch of the I operator between 5-point and reduced kinematics; reduces to the Z form for Q2=-Q1, Q4=-Q3. C2 uses the leg-2 charge (cf. RV6a). gq/qg (`col_rvewk_is_qg`) still to do | `SingLimits/mod_limits_rv.f90` | ✅ ns (2026-11-12) |
| RV4 | OpenLoops ids: W registers ids 3-4 only; the Z code reads `OL_id(j+3)`. Check the W+ ordering (`register_red_qcd` vs `register_red_ew` differ) | `OpenLoops/mod_ol_interface.f90` | 🟡 |
| RV5 | `rvewknf_is` for W (see V2) | – | ⏸ |

### 4b. QCD one-loop x real photon

| ID | Task | Where | Status |
|---|---|---|---|
| RV6 | `rvqcd_is` ns generalised and checked for W±: all limits close in the `-checkint` scan (Z, W+, W-) after RV6a, RV6b | `mod_xsects_nnlo_rvqcd.f90` | ✅ limits (2026-11-12) |
| RV7 | `rvqcd_fs_53/54` ns generalised and checked for W±: all limits close in the scan after RV7a, RV7b, both the charged-lepton and the neutrino sector | idem | ✅ limits (2026-11-12) |
| RV6a | C2 used the leg-1 charge: **confirmed by the scan** (5\|\|2 limit did not close for W+ and W-), fixed with `multiply_IS_charges_sq(res_tmp,2)` | idem | ✅ fixed |
| RV6b | SC2 was a copy of SC1 (leg-1 charge): fixed, computed with the leg-2 charge as at NLO EW | idem | ✅ fixed |
| RV7a | C used `Q_lep2` in both sectors: **confirmed** (the neutrino sector, `fs_53` for W+ and `fs_54` for W-, had a counterterm without a pole). Fixed with the charge of the emitting lepton, `Qfs2 = Q3^2` or `Q4^2`. W result of the neutrino sector: -17.6 -> 0.63 fb (W+), -14.0 -> 0.45 fb (W-) at test statistics | idem | ✅ fixed |
| RV7b | CS used `Q_lep2`: same fix as RV7a | idem | ✅ fixed |
| RV7c | Charged-lepton sector for W±: C, S, CS close in the scan | idem | ✅ |
| RV7d | **Verify for W**: sum of the three partitions (`rvqcd_is` + `fs_53` + `fs_54`) is independent of the damping/partition choice | `mod_partitions.f90` | ❌ to verify |
| RV8 | `rvqcd_is` aq, qa generalised for W±: `_gen` amplitudes, charge of the quark from the photon splitting, `transition`, `lep_ids`; the gamma gamma limit (C2 for aq, C1 for qa) only for Z; legacy under `_oldcode`. Z identical to the old code; W limits close (`-checklim`: 1e-2, 1e-3, 1e-4 at depth 1e-4, 1e-6, 1e-8; scan ok) | idem | ✅ limits (2026-11-12) |
| RV9 | Confirm that no additional RV sector is needed for W (e.g. photon-induced with EW loop, final-state EW-loop sectors): none exists for Z either | – | question |

## 5. Double real

General amplitudes exist for every channel with g + gamma or q + gamma in the final state
(`res_tree_ga_*_gen`). The work is in the integrands.

| ID | Task | Where | Status |
|---|---|---|---|
| RR1 | `ns_ga`, all 10 sectors: general code for W±. Fix ids (:995). Legacy branch of 5161a/5262a is broken (:433, :533 use an uninitialised array), so `oldcode=T` is not a valid reference there | `mod_xsects_nnlo_rr_ns.f90` :217-2637 | 🟡 |
| RR2 | `gq`, 7 sectors: general matrix elements and charges, but no `Vcharge` branch; ids are Z; `%part` unset. "agrees with oldcode=T", limits unchecked for W | `mod_xsects_nnlo_rr_gq.f90` | 🟡 |
| RR3 | `qg`, 7 sectors: as RR2; `pause` in 5262b/d and 5263/5264; ids Z in 526k | `mod_xsects_nnlo_rr_qg.f90` | 🟡 |
| RR4 | Scale-variation weights: `respdf = respdf_tmp(1,1)+respdf_tmp(1,2)` gave the central value to all weights in `gq`, `qg`, `aq`, `qa` (two sectors each); now `respdf_tmp(:,1)+respdf_tmp(:,2)`. Central values unchanged | idem | ✅ 2026-11-12 |
| RR5 | Neutrino partition (5164/5264 for W-, 5163/5263 for W+) in `ns_ga`, `gq`, `qg`: all limits close in the scan for W+ and W- | `ns_ga`, `gq`, `qg` | ✅ limits |
| RR6 | `aq`, `qa` (10 sectors each) generalised for W± by a scripted conversion: `_gen` amplitudes, `get_respdf_gen`, charge of the quark from the photon splitting, `lep_ids`, W guards removed; legacy under `_oldcode`. Z: all sectors reproduce the old numbers; W+ and W-: all limits close in the scan | `mod_xsects_nnlo_rr_aq.f90`, `_qa.f90` | ✅ limits (2026-11-12) |
| RR7 | `ag`, `ga` (7 sectors each) generalised for W±: flavour sums per leg with the new `collapse_leg` (quarks or antiquarks of one leg, optionally with their charge). Z reproduces the old numbers (`ag rr_5262a` to 6e-11 without grid adaptation, 2.5e-6 with it: rounding in a strongly cancelling sector); W limits close | `mod_xsects_nnlo_rr_ag.f90`, `_ga.f90`, `mod_aux_sectors.f90` | ✅ limits (2026-11-12) |
| RR8 | Limits built on gamma gamma -> l l: kept in the code, `res_treeAA_aa_gen` is zero for W, and the scan shows that the W integrand is integrable there without them | RR6, RR7 | ✅ |
| RR9 | `res_tree_ga_ag_gen`/`_ga_gen` for W: not compared with MadGraph yet, but all their collinear limits (onto validated lower-multiplicity amplitudes) close for W+ and W-. Fixed on the way: `res_tree_ga_aq_gen`/`_qa_gen` filled a spurious photon-photon entry | `Amplitudes/` | 🟡 limits ✅, MG check open |
| RR10 | Four-quark singular channels (`ns_qqb`, `ns_qq`, TC limit only): Z legacy amplitudes. Need q qb' -> l nu q'' qb'' amplitudes at O(alpha_s alpha) interference and the TC counterterm with W flavours | `mod_xsects_nnlo_rr_ns.f90` :2641-2883 | ❌ |
| RR11 | Four-quark finite channels (`ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w`): OpenLoops registration stops for W ("TO IMPLEMENT OL rr for W in the 4q case"). Enumerate W flavour channels, register, add luminosities (`qQpb_lumi_wp/_wm` exist, unused) | `mod_ol_interface.f90` :248-284, `rr_ns` :2885-3126 | ❌ |
| RR12 | List, on paper, all partonic channels at O(alpha_s alpha) for W+ and W- with their singular limits, and map each to a sector. Use it to confirm nothing is missing (notes exist for Z in `notes/chiara_federica`) | new note | ❌ |

## 6. Integrated subtraction ("S" sectors)

Mathematica input: `math/DY+_qq.wl` (W+, ns). Subtraction formulae for W+ and W- in all
channels are to be provided by Christian.

### 6a. ns channel

| ID | Task | Where | Status |
|---|---|---|---|
| S1 | `s_12`, `s_vqcd`, `s_vewk` for W+: **validated against `math/subtraction-notebooks/W.nb`** (2026-11-12). The notebook derivation was re-run (all poles vanish) and the Fortran W+ branches, the elastic functions (`fin_ns_vqcd_cc`, `genGnloQCD_cc`, `legfinalcalG_cc`, `fin_ns_v7_cc`) and the HOPPET kernels, transcribed to Mathematica, reproduce the full finite remainder in every structure (FLVqcd, FLVew, FLM double boost, two single boosts, elastic). Script: `math/subtraction-notebooks/checks/` | `mod_xsects_nnlo_s_ns.f90`, `IntSub/` | ✅ symbolic |
| S2 | `s_12`, `s_vqcd`, `s_vewk` for **W-**: done with one W branch for both signs (`ilept_chg`, `Qlept_chg`); CP test against W+ with antiproton PDFs passes (CHECKS.md) | idem | ✅ 2026-11-13 |
| S3 | `s_oqcd` W±: general, limits checked. Open: `calG_ONLOQCD_ns` and `genGnloQCD_cc` differ in the sign of the `Qlept*Q_IS` cross terms once the scale logs are aligned; the note in `math/tex` already flags they are not equal. Establish which is right in which context | :837-1528; `mod_int_sub_nnlo.f90` :263, :428 | 🟡 |
| S4 | `s_oewk_is`: the routine wired to `-sec s_oewk_is` is an unfinished draft (broken for Z too); the working one is `_raoul` under `s_oewk_is_ra`. Rewire, delete the draft | :1531-2491, :3181-4216; `Chooser/mod_run_proc_s_ns.f90` | 🔧 |
| S5 | `s_oewk_fs_53/54` (`_raoul`): general; "limits checked for W- only". Check W+; remove legacy `oewk_fs_5i`; `Q_lep2` left at :4861 | :4220-5019 | 🟡 |
| S6 | `s_qqb`, `s_qq`: Z only, no guard (goes with RR10) | :698-835 | ❌ |
| S7 | `s_vewknf` (see V2) | :3111 | ⏸ |

### 6b. gq, qg

| ID | Task | Where | Status |
|---|---|---|---|
| S8 | `s_oqcd`, `s_oewk_is`, `s_oewk_fs_53/54`: general; "agreement oldcode T/F for NC, limits work for + and - CC". The W elastic function `calG_ONLOQCD_gq` is never exercised by the NC comparison | `mod_xsects_nnlo_s_gq.f90`, `_qg.f90` | 🟡 |
| S9 | *Suspected*: `gq` passes `Lim_Ei(2)` (energy of the radiated parton) in the scale log of `calG_ONLOQCD_gq` (:499, :562), `qg` passes `Lim_Ei(1)` | idem | 🟡 |
| S10 | `s_12`, `s_vewk` (and `s_vewknf`): Z only | :66-255 in both | ❌ |
| S10b | `gq s_oewk_fs_53/54` zero at every point: two debug lines set `xx(xE)`, `xx(xRHO)` to fixed values. Removed; Z gq = qg mirror | `mod_xsects_nnlo_s_gq.f90` | ✅ 2026-11-13 |

### 6c. aq, qa, ag, ga

| ID | Task | Where | Status |
|---|---|---|---|
| S11 | `aq`: `s_12`, `s_vqcd`, `s_oewk` Z only; `s_oqcd` first pass (general amplitudes only, Z ids, debug `stop`) | `mod_xsects_nnlo_s_aq.f90` | ❌ |
| S12 | `qa`: everything Z only | `mod_xsects_nnlo_s_qa.f90` | ❌ |
| S13 | `ag`, `ga`: everything Z only. Terms built on gamma gamma -> l l vanish for W (cf. RR8) | `mod_xsects_nnlo_s_ag.f90`, `_ga.f90` | ❌ |
| S14 | HOPPET tables: `xPij_a` (aq/qa) carry `Qdn2/Qup2` weights; make them charge-free and apply charges per leg | `IntSub/mod_hoppet_nnlo.f90` :319, `mod_hoppet_nlo.f90` :87 | ❌ |

### 6d. Mathematica side

| ID | Task | Status |
|---|---|---|
| M1 | Add to `math/`: W- ns; W± for gq/qg, aq/qa, ag/ga, four-quark; ONLO functions | ❌ (Christian) |
| M2 | One notebook convention for charges/legs shared with the Fortran (`Qq, Qqp, Qlept, ilept`) | ❌ |

## 7. Validation (mandatory for every item above)

### C1. Scaling of local counterterms in singular limits
For each RR, RV and ONLO sector and each limit (S, C, SC, TC, double soft, ...):
(Hard + CT)/Hard must vanish with the expected power as the limit is approached.

| ID | Task | Status |
|---|---|---|
| C1.1 | New `-checkint`: integrability scan of any sector (takes points from a Vegas warm-up, sends subsets of the variables to the edge, prints the median of \|f\| d^k at three depths). Works for every integrand without wiring. Tested with planted bugs: detects a wrong charge and 1% errors, also in nested counterterms. `lim_proc` (pairwise scaling) still has no caller | ✅ scan (2026-11-12) |
| C1.2 | Fix wiring: `xsect_nloewk_r_is_ns` fills `limval_nlo_is` but the checker reads `FintNLO_ew`; `nloqcd r_is` routed to the wrong reader; `rr_qg` fills no `limval_nnlo` | 🔧 |
| C1.3 | Extend `xsect_selector` beyond RR `ns_ga`: RR gq/qg/aq/qa/ag/ga/4q, all RV, all ONLO | ❌ |
| C1.4 | Make it W-aware (which limits are singular per sector: no lepton-collinear pole on the neutrino; no gamma gamma TC) | 🔧 |
| C1.5 | Output one table per sector (limit, slope, grade); store reference tables for Z, W+, W- | ❌ |
| C1.6 | Scan run for all RV and RR sectors with a photon or from photon-induced channels (63 sectors: `ns` rvqcd, `aq`/`qa` rvqcd, RR `ns_ga`, `gq`, `qg`, `aq`, `qa`, `ag`, `ga`), Z, W+, W-: all pass. Still to run: `rvewk`, four-quark, S sectors | 🟡 |
| C1.7 | Priority targets for W (cannot be seen in Z): RV6a, RV6b, RV7a, RV7b, RV3, N2, RR5, RR8. Each needs a recorded scaling table for W+ and W- before it is ticked | ❌ |

### C2. Integrated counterterm left-over, Fortran vs Mathematica
| ID | Task | Status |
|---|---|---|
| C2.1 | Extend `math/validate_DYp7.wl` from term [7] to all of `FINITEDYp[1..7,60]` (symbolic, kernels included) | 🔧 |
| C2.2 | Numerical check: evaluate each finite-remainder term at fixed `(z, eta, mu)` in both codes (Fortran test driver + notebook export), per channel and boson | ❌ |
| C2.3 | Resolve S3 (`calG_ONLOQCD_ns` vs `genGnloQCD_cc`) and S9 with this check | ❌ |
| C2.4 | Same for W-, gq/qg, photon-induced and four-quark once M1 is in | ❌ |
| C2.5 | Store logs of the validation runs in the repo | ❌ |

### C3. POWHEG benchmarks

Christian's POWHEG-BOX-RES generators (`/scratch/cbiello/POWHEG-BOX-RES/MiNNLOPS_res/`, GitHub `cbiello/NSC-lnuEW`, `cbiello/NSC-lnuQCD`), used at fixed order.

| | `NSC-lnuQCD` | `NSC-lnuEW` |
|---|---|---|
| Process | pp -> l nu at NLO QCD | pp -> l nu at NLO EW, photon-induced channels included (`flg_with_em`, `flg_QEDonly`, photon PDF) |
| Boson | `chooseDY` +1 / -1 (W+ / W-); 0 = Z, flavour lists still commented out | same |
| Amplitudes | Born and real from OpenLoops `ppllj`; virtual analytic, in the nested-subtraction normalisation | all from OpenLoops `ppllj_ew` |
| Scheme | G_mu (`ewscheme 1`), complex-mass | same; alpha fixed to the OpenLoops value 7.5598868e-3, top width 0 |
| Cuts | generation cut `qmin`; analysis `pwhg_analysis-W_ATLAS.f` | same analysis file |
| Test input | `testrun-lhc`, 13 TeV, PDF 303600 | `testrun-lhc`, 13 TeV, PDF 325100 (NNPDF31 luxqed) |

The analysis is the twin of `analysis/mod_cut_histo_W_ATLAS.f90`: mT in [200,5000] GeV, pT(l) > 65, pT(nu) > 85, pT(l) pT(nu) > 35^2, |y| < 2.4; histograms `total`, `rate`, `massT` (10 GeV bins from 200 GeV).

| ID | Task | Status |
|---|---|---|
| C3.1 | NLO QCD vs `NSC-lnuQCD`, 13 TeV, fixed scale MZ, `W_ATLAS` cuts. W+: 1737.51(61) vs 1736.68(31) fb; W-: 1203.14(44) vs 1203.02(18) fb; LO agrees too. Setup, scripts and sector numbers in `f90/tests/powheg/README.md` | ✅ 2026-11-11 |
| C3.2 | NLO EW vs `NSC-lnuEW`, same setup, lepton dressed with ΔR 0.1. W+: 1327.29(25) vs 1326.99(21) fb; W-: 929.24(25) vs 929.23(13) fb. Needed three fixes in `NSC-lnuEW/openloops2.f` (real and virtual normalised by alpha_s/2pi as the public POWHEG-BOX-RES expects, scale passed to OpenLoops) and `detailedNLO 1`; see the README | ✅ 2026-11-12 |
| C3.3 | NLO QCD generator for l l' + photon: `NSC-lnuaEW` (POWHEG-BOX-RES, OpenLoops `pplla`, `ppllaj`), Z, W+, W- | ✅ 2026-11-13 |
| C3.4 | NLO EW generator for l l' + jet: `NSC-lnujQCD` (POWHEG-BOX-RES; l nu j for W+/W-, l l j for Z; photon-induced channels, resolved jet), built like `NSC-lnuaEW` | 🔧 next |
| C3.5 | Analyses in this code that match C3.3 and C3.4: resolved photon (cf. Z analysis `lla_nloqcd`) and resolved jet, for W+ and W- | ❌ |
| C3.6 | NLO QCD l l' gamma: sum of the 39 sectors with a resolved photon against `NSC-lnuaEW` (no damping): W+ 6.1928(21) vs 6.1937(75) fb; W-, Z being rerun. Needs `-qmin_rsc 1` and 1.5M/5M points (`f90/tests/powheg/lla/`) | ✅ 2026-11-13 |
| C3.7 | Compare C3.4 with the sum of the sectors here with one resolved jet (`rvewk` + RR + ONLO), channel by channel | ❌ |
| C3.8 | Script that runs and sums the right sectors for C3.1, C3.2, C3.6, C3.7 and prints the ratio to POWHEG | ❌ |

### C4. Further checks
| ID | Task | Status |
|---|---|---|
| C4.1 | `-MGcheck` for all `_gen` trees incl. g+gamma and four-quark, W±; add loop amplitudes vs OpenLoops/MadLoop points | 🔧 |
| C4.2 | Pole check (`-checkpoles`) for `_gen` one-loop wrappers | ❌ |
| C4.3 | Z regression: `OLDCODE=true` vs `false` builds and old published numbers after every change; script it (I11) | manual |
| C4.4 | Independence of technical parameters: `buff*`, damping exponents, partition choice; mu_R/mu_F logs | ❌ |
| C4.5 | W+ vs W- consistency: charge conjugation relates partonic W- to W+ with q <-> qb | ❌ |
| C4.6 | `qa` = mirror of `aq`, `ga` = mirror of `ag`, `qg` = mirror of `gq` (automated diff test with x1 <-> x2) | ❌ |
| C4.7 | Inclusive mixed corrections vs literature for on-shell/off-shell W where available | ❌ |

## 8. Suggested order

1. I1-I4, I7 (build, guards, debug removal, ids).
2. C1.1-C1.5 (scaling checker usable everywhere) and C2.1-C2.2 (math check harness).
3. ns channel to completion: RR1, RV6-RV7, RV1+RV3, S2-S5; then C1/C2 on ns.
4. gq/qg: RR2-RR5, RV2, S8-S10.
5. C3 POWHEG benchmarks: C3.1-C3.2 now (NLO is complete), C3.3-C3.5 in parallel, C3.6-C3.7 when the RV/RR sectors of a channel are done.
6. Photon-induced: N2-N3, RV8, RR6-RR9, S11-S14.
7. Four-quark: RR10-RR11, S6.
8. VV and nf pieces (section 3), then phenomenology.

## 9. Sector status matrix (fill as work proceeds)

Columns: code generalised / limits scale (C1) / finite part vs math (C2) — for W+ and W-.

| corr ch sec | W+ code | W+ C1 | W+ C2 | W- code | W- C1 | W- C2 |
|---|---|---|---|---|---|---|
| nnlo ns rvewk_is | 🟡 | ✅ | – | 🟡 | ✅ | – |
| nnlo ns rvqcd_is | 🟡 | | – | 🟡 | | – |
| nnlo ns rvqcd_fs_53/54 | 🟡 | | – | 🟡 | | – |
| nnlo gq/qg rvewk_is | 🔧 | | – | 🔧 | | – |
| nnlo aq/qa rvqcd_is | ❌ | | – | ❌ | | – |
| nnlo ns_ga rr (10 sectors) | 🟡 | | – | 🟡 | | – |
| nnlo gq rr (7) | 🟡 | | – | 🟡 | | – |
| nnlo qg rr (7) | 🟡 | | – | 🟡 | | – |
| nnlo aq/qa rr (10+10) | ❌ | | – | ❌ | | – |
| nnlo ag/ga rr (6+6) | ❌ | | – | ❌ | | – |
| nnlo ns 4q rr | ❌ | | – | ❌ | | – |
| nnlo ns s_12 / s_vqcd / s_vewk | 🟡 | – | 🔧 | ❌ | – | |
| nnlo ns s_oqcd | 🟡 | 🟡 | | 🟡 | 🟡 | |
| nnlo ns s_oewk_is(_ra) | 🟡 | | | 🟡 | 🟡 | |
| nnlo ns s_oewk_fs_53/54 | 🟡 | | | 🟡 | 🟡 | |
| nnlo ns s_qqb / s_qq | ❌ | – | | ❌ | – | |
| nnlo gq/qg s_oqcd, s_oewk_* | 🟡 | 🟡 | | 🟡 | 🟡 | |
| nnlo gq/qg s_12, s_vewk | ❌ | – | | ❌ | – | |
| nnlo aq/qa s_* | ❌ | | | ❌ | | |
| nnlo ag/ga s_* | ❌ | | | ❌ | | |

## 10. Open questions for the team

1. CKM: keep the unit matrix? b-initiated channels for W?
2. Do the `nf` sectors (`*nf`) and the `aa`-type terms have W analogues we must keep?
3. Which soft function is correct in the ONLO z=1 terms (S3), and which energy enters the gq log (S9)?
4. For sectors whose only TC limit is gamma gamma -> l l (RR8): drop the sector or keep it?
5. Four-quark W channels: amplitudes from OpenLoops (`pplnjj_ew` is public) or analytic?
6. POWHEG: exact common input for C3.1-C3.2 (PDF set, scales, `qmin`, widths); photon and jet definition for C3.3-C3.7.
7. Is any RV sector missing for W (RV9)?
