# PLAN — mixed QCD x EW NNLO for off-shell W+ and W-

Shared status board. Read [STRUCTURE.md](STRUCTURE.md) first for names and conventions.
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
| LO, NLO QCD, NLO EW (ns, gq, qg, aq, qa) | ✅ | ✅ (aq/qa 🟡) | ✅ (aq/qa 🟡) |
| VV (two-loop) | ✅ | ⏸ | ⏸ |
| RV, EW loop + gluon (`rvewk_is`: ns, gq, qg) | ✅ | 🔧 amplitudes ready, integrand Z | 🔧 same |
| RV, QCD loop + photon (`rvqcd_*` ns) | ✅ | 🟡 | 🟡 |
| RV, QCD loop + photon (aq, qa) | ✅ | ❌ | ❌ |
| RR `ns_ga` | ✅ | 🟡 | 🟡 |
| RR `gq`, `qg` | ✅ | 🟡 | 🟡 |
| RR `aq`, `qa`, `ag`, `ga` | ✅ | ❌ | ❌ |
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
| I7b | **`%part` is the field the W analyses use to tell lepton from neutrino** (`W_ATLAS`: `event%part(3)`), not `%ids`. It is unset in all of `rr_gq`, `rr_qg`, `s_aq_oqcd` (analysis prints 'where is the final state lepton?' and cuts on undefined pT, W integrand zero or garbage) and set to `[id_el,-id_el]` for both W in `nloewk` `aq`/`qa` `r_is` (:768-775, :953-960), which swaps the lepton and neutrino pT cuts for W+. Proposed: let the analyses read `event%ids(3)` (always set, see I7) and drop the lepton entries of `%part` | `analysis/`, `Sectors/` | ❌ to decide |
| I8 | Flavour selection in the tree amplitudes by integer tests: `is_W_pair(i,j)` and `is_W_quark(i)` in `mod_parms` replace the float comparison and the eight hand-written `cycle` blocks; `current_charge`, `Qgeneration` removed. b is now skipped for W in all single-quark amplitudes (was missing in `res_tree_ga_qg_gen`, `ga_aq`, `ga_qa`). 55 sectors x (W+, W-) identical. Open for the team: CKM beyond the unit matrix (question 1) | `Amplitudes/mod_amplitudes_tree_ppll.f90`, `Parms/` | ✅ 2026-11-11 |
| I9 | Dead W prototypes deleted: `xsect_lo_ns_wp/_wm`, `xsect_nloqcd_r_is_ns_wp/_wm`, `xsect_nloewk_r_is_ns_wp/_wm`, `qQpb_lumi_wp/_wm`, `get_coupl_w`, `mod_amplitudes_tree_ppnul.f90`, `mod_xsects_nnlo_rvqcd_original.f90`, `makefile_chiara`, two `*.bak`, commented chooser entries. Kept: `ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w` (Z four-quark channels with W exchange, now guarded for W) | – | ✅ 2026-11-11 |
| I10 | Run header prints `lepton charge = Qel` for every boson; print `Q3,Q4` | `Parms/mod_proc_parms.f90:426` | ❌ |
| I11 | Small regression suite: script that runs each sector for Z/W+/W- at low statistics with fixed seed and compares to stored numbers | new `f90/tests/` | ❌ |
| I12 | Analysis: do not recombine photons with the neutrino (`TODO_CB` item 4); make the analysis follow `-boson` (it is still chosen at compile time) | `analysis/mod_cut_histo_aux.f90` | ❌ |
| I13 | Uninitialised variables: several integrands sum `kin(:)` / `Fint(:)` entries that are only set in some branches (found in NLO EW aq/qa). Build once with `-finit-real=snan -ffpe-trap=invalid` and fix what trips | `Sectors/` | ❌ |
| I14 | Vegas `rebin` goes out of bounds (`mod_my_vegas.f90:388`) when the integrand is zero at every point or the warm-up is tiny; make it fail with a clear message | `Vegas/mod_my_vegas.f90` | ❌ |
| I15 | Analysis `V_minimal` returns zero in every sector for Z and W; fix or remove | `analysis/mod_cut_histo_V_minimal.f90` | ❌ |
| I16 | Z, old vs new code (same seed): agreement to ~10 digits in most sectors. The four exceptions were studied on 2026-11-11 (I16a-d). Verified fixes are ready as a patch, not applied | – | studied |
| I16a | `nloewk aq s`, `qa s` (0.4% off): bug in the new code, Born multiplied by `Q2(fa)*Q2(fb)` (Q^4 for Z) instead of one Q^2 of the quark from the photon splitting. **Fixed**: loop replaced by `multiply_IS_charges_sq(res_lo_new,1)` (aq) and `(...,2)` (qa). Z new = old to all digits. W values change (W- aq: -3.7e-3 -> -2.3e-2 fb at test statistics) | `mod_xsects_nloewk_vs.f90` | ✅ 2026-11-11 |
| I16b | `nloewk qa r_is` (old = -19221 fb): **bug in the old branch**. In C1, `Pqq(z)/(1-z)` is applied inside the old branch and again after `#endif` (in `aq` the second one sits inside the new branch). With it once, old = new to 12 digits. New code is right | `mod_xsects_nloewk_r.f90` (`xsect_nloewk_r_is_qa`) | ✅ fixed (old = new, -9.06994348143 fb) |
| I16c | `ns_ga rr_5161a/5262a` (old = -58 fb): **old branch incomplete**, five defects: S6 never calls `res_tree_g_qqb`; S6 and S5S6 pass the unset `res_tmp_vect` instead of `res_tmp_vect_old` and leave the result in `respdf_vect` while the code reads `respdf_tmp`; TC, TCS6, TCS5 fill `respdf` instead of `respdf_tmp(:,1)`. With these repaired, old = new to 11 digits, so the **new code is validated for Z in sector a** | `mod_xsects_nnlo_rr_ns.f90` (`xsect_nnlo_rr_5i6ia_ns_ga`) | ✅ fixed (old = new to 11 digits) |
| I16d | `rvqcd_fs_53/54`, `rvqcd_is` (agree to 1e-7..1e-9): not a bug. Old code uses the analytic one-loop amplitude (`res_qcdloop_j`), new code uses OpenLoops; the difference is numerical precision | – | ✅ understood |

## 2. NLO (done; remaining loose ends)

| ID | Task | Where | Status |
|---|---|---|---|
| N1 | NLO QCD ns/gq/qg and NLO EW ns: implemented for W±, benchmarked against POWHEG (commit `346fdea`) | – | ✅ |
| N2 | NLO EW `aq`/`qa` `r_is` for W: limits checked with `-checklim` (2026-11-11). Z closes (ratio ~ sqrt(depth)); W does not, two bugs: (a) `qa` C2 uses the charge of leg 1, must be leg 2 (:1053); (b) `res_tree_a_aq_gen` / `res_tree_a_qa_gen` keep the b quark for W while the Born does not (missing `abs(.)==5` skip, see I8). Both fixed: all four W limits close, Z unchanged digit for digit. W+ aq/qa = 1.508(7)/1.494(6) fb, W- = 0.856(4)/0.847(3) fb (`W_ATLAS`) | `mod_xsects_nloewk_r.f90`, `mod_amplitudes_tree_ppll.f90` :484-550 | ✅ 2026-11-11 |
| N3 | NLO EW `s_aq`/`s_qa`: double charge factor, fixed, see I16a. Redo the POWHEG comparison for the photon-induced channels (N5) | `mod_xsects_nloewk_vs.f90` | ✅ fixed, to re-benchmark |
| N4 | NLO QCD `gq`/`qg` integrated collinear: "still work in progress" in `f90/README`; confirm closed | `mod_xsects_nloqcd_vs.f90` | 🟡 |
| N5 | Redo the POWHEG comparison for W+ and W-, all channels including photon-induced: see C3.1, C3.2 | – | ❌ |

## 3. Double virtual (out of scope for now)

| ID | Task | Status |
|---|---|---|
| V1 | Two-loop mixed amplitude for q qb' -> l nu, its interface, `vvfc/vvxf/vvnf` for W. `mod_xsects_nnlo_vv.f90` is Z only | ⏸ |
| V2 | `nf` pieces for W (`vvnf`, `rvewknf_is`, `s_vewknf`): decide whether they stay separate | ⏸ |

## 4. Real-virtual

### 4a. EW one-loop x real gluon

| ID | Task | Where | Status |
|---|---|---|---|
| RV1 | `rvewk_is` ns: the **amplitudes are ready for W±** (`res_ewkloop_g_qqb_gen`: OpenLoops ids, W-charge finite part, flavour matrix; reduced `res_ewkloop_qqb_gen`; commit `b28b7a7`). The **integrand is still Z**: it calls the `_gen` routine only as a test (`print` + `pause`), then uses the legacy amplitude, `ns_lumi_splitb`, `[Qdn2,Qup2]`, Z ids. To do: hard, C1, C2 with `_gen` amplitudes, `get_respdf_gen`, per-leg charges, W ids | `mod_xsects_nnlo_rvewk.f90` :33-168 | 🔧 |
| RV2 | `rvewk_is` gq, qg: same situation (`res_ewkloop_g_gq_gen`, `_g_qg_gen` ready and called as a test; integrand Z) | :172-380 | 🔧 |
| RV3 | One-loop collinear limits with W charges. `col_rvewk_is_ns/_qg` use `-2 log(s_i3/s_i4)` times `Q_lep*[Qdn,Qup]`, which assumes `Q3=-Q4` and one quark charge. Re-derive for general `(Q_q, Q_q', Q_3, Q_4)` | `SingLimits/mod_limits_rv.f90` | ❌ |
| RV4 | OpenLoops ids: W registers ids 3-4 only; the Z code reads `OL_id(j+3)`. Check the W+ ordering (`register_red_qcd` vs `register_red_ew` differ) | `OpenLoops/mod_ol_interface.f90` | 🟡 |
| RV5 | `rvewknf_is` for W (see V2) | – | ⏸ |

### 4b. QCD one-loop x real photon

| ID | Task | Where | Status |
|---|---|---|---|
| RV6 | `rvqcd_is` ns: **integrand generalised** (ids/part for Z, W-, W+; `ol_res_qcdloop_a_qqb_gen`, `res_qcdloop_qqb_gen`, `get_qed_eik_gen`, `get_respdf_gen`). To verify for W: see RV6a, RV6b | `mod_xsects_nnlo_rvqcd.f90` :45-321 | 🟡 |
| RV7 | `rvqcd_fs_53/54` ns: **integrand generalised**, commit note "limits checked but NOT deeply". To verify for W: see RV7a, RV7b | idem :346-577 | 🟡 |
| RV6a | **Verify for W**: C2 uses the leg-1 charge, `multiply_IS_charges_sq(res_tmp,1)` (:223). Photon collinear to leg 2 should carry `Q_IS` of leg 2, as in NLO EW (`mod_xsects_nloewk_r.f90` :220). Invisible for Z (equal charges), so the old-vs-new comparison cannot detect it. Test: C1.6 scaling in the 5\|\|2 limit for W+ and W- | :223 | ❌ to verify |
| RV6b | **Verify for W**: SC2 is a copy of SC1, `respdf_tmp(:,3) = respdf_tmp(:,2)` with leg-1 charge (:277-279). NLO EW computes it separately with leg 2 (:276-278 there). Test: soft-collinear 5\|\|2 scaling for W± | :277-279 | ❌ to verify |
| RV7a | **Verify for W**: C uses `Q_lep2` for both sectors (:492). In the neutrino sector (`fs_54` for W-, `fs_53` for W+) the hard term has no collinear pole, so the CT should vanish; NLO EW uses `Qsq_FS(icoll-2)` (`mod_xsects_nloewk_r.f90` :479). Test: run the neutrino sector, check finiteness and collinear scaling | :492 | ❌ to verify |
| RV7b | **Verify for W**: CS uses `res1red*Q_lep2` (:530), same issue as RV7a. Test: soft-collinear scaling in the neutrino sector for W± | :530 | ❌ to verify |
| RV7c | **Verify for W**: charged-lepton sector (`fs_53` for W-, `fs_54` for W+): repeat the limit check with the scaling tool (C, S, CS), both signs | :346-577 | ❌ to verify |
| RV7d | **Verify for W**: sum of the three partitions (`rvqcd_is` + `fs_53` + `fs_54`) is independent of the damping/partition choice | `mod_partitions.f90` | ❌ to verify |
| RV8 | `rvqcd_is` aq, qa: **integrand is Z only** (no `oldcode`, no `Vcharge` branch, Z ids, `aq_lumi`). The general amplitudes `ol_res_qcdloop_a_aq_gen`/`_qa_gen` are ready and called, but the result is discarded (:633, :769). Generalise, then verify for W: collinear gamma -> q qb limit with the charge of the produced quark, both signs | idem :583-853 | ❌ |
| RV9 | Confirm that no additional RV sector is needed for W (e.g. photon-induced with EW loop, final-state EW-loop sectors): none exists for Z either | – | question |

## 5. Double real

General amplitudes exist for every channel with g + gamma or q + gamma in the final state
(`res_tree_ga_*_gen`). The work is in the integrands.

| ID | Task | Where | Status |
|---|---|---|---|
| RR1 | `ns_ga`, all 10 sectors: general code for W±. Fix ids (:995). Legacy branch of 5161a/5262a is broken (:433, :533 use an uninitialised array), so `oldcode=T` is not a valid reference there | `mod_xsects_nnlo_rr_ns.f90` :217-2637 | 🟡 |
| RR2 | `gq`, 7 sectors: general matrix elements and charges, but no `Vcharge` branch; ids are Z; `%part` unset. "agrees with oldcode=T", limits unchecked for W | `mod_xsects_nnlo_rr_gq.f90` | 🟡 |
| RR3 | `qg`, 7 sectors: as RR2; `pause` in 5262b/d and 5263/5264; ids Z in 526k | `mod_xsects_nnlo_rr_qg.f90` | 🟡 |
| RR4 | Scale-variation weights: `respdf = respdf_tmp(1,1)+respdf_tmp(1,2)` assigns the central value to all weights (`gq` :1200, :1330; `qg` :1236, :1353) | idem | 🔧 |
| RR5 | Neutrino partition (5164/5264 for W-, 5163/5263 for W+): confirm hard minus soft is integrable and collinear CTs vanish | `ns_ga`, `gq`, `qg` | 🟡 |
| RR6 | `aq`, `qa` (10 sectors each): pure Z code. Generalise: `_gen` amplitudes, `gen_lumi`, gamma -> q qb with the charge of the **produced** quark (isospin partner for W), ids | `mod_xsects_nnlo_rr_aq.f90`, `_qa.f90` | ❌ |
| RR7 | `ag`, `ga`: pure Z code. Generalise; the flavour sum `ndn*Qdn2+nup*Qup2` becomes per-leg (q and qb' carry different charges) | `mod_xsects_nnlo_rr_ag.f90`, `_ga.f90` | ❌ |
| RR8 | Limits built on gamma gamma -> l l (`res_treeAA_aa`: TC in `aq` 5262, `qa` 5161, `ag` 5262, `ga` 5161) do not exist for W. Decide per sector: drop, or keep the sector with soft/C5/C56 only | RR6, RR7 | ❌ |
| RR9 | Check `res_tree_ga_ag_gen`/`_ga_gen`: W+ and W- use the identical expression `res_tmp(1,3)*nup`, and the W branch of `res_tree_jj_ag_gen` uses `iconf1` only. Extend `-MGcheck` to cover it | `Amplitudes/` | 🟡 |
| RR10 | Four-quark singular channels (`ns_qqb`, `ns_qq`, TC limit only): Z legacy amplitudes. Need q qb' -> l nu q'' qb'' amplitudes at O(alpha_s alpha) interference and the TC counterterm with W flavours | `mod_xsects_nnlo_rr_ns.f90` :2641-2883 | ❌ |
| RR11 | Four-quark finite channels (`ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w`): OpenLoops registration stops for W ("TO IMPLEMENT OL rr for W in the 4q case"). Enumerate W flavour channels, register, add luminosities (`qQpb_lumi_wp/_wm` exist, unused) | `mod_ol_interface.f90` :248-284, `rr_ns` :2885-3126 | ❌ |
| RR12 | List, on paper, all partonic channels at O(alpha_s alpha) for W+ and W- with their singular limits, and map each to a sector. Use it to confirm nothing is missing (notes exist for Z in `notes/chiara_federica`) | new note | ❌ |

## 6. Integrated subtraction ("S" sectors)

Mathematica input: `math/DY+_qq.wl` (W+, ns). Subtraction formulae for W+ and W- in all
channels are to be provided by Christian.

### 6a. ns channel

| ID | Task | Where | Status |
|---|---|---|---|
| S1 | `s_12`, `s_vqcd`, `s_vewk` for W+: implemented (`FINITEDYp[1..7,60]`); only term [7] and the soft functions are validated symbolically | `mod_xsects_nnlo_s_ns.f90` :83-696, `IntSub/mod_int_sub_nnlo.f90` | 🟡 |
| S2 | Same three for **W-**: `stop` at :141, :520, :642. Helpers already take `(Qlept, ilept)`; add branches with `Qlept=-1`, `ilept=3`, `s13/s23`, ids `[id_el,-id_nue]`. Needs the W- Mathematica reference | idem | ❌ |
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
| C1.1 | Make the scaling checker reachable: `lim_proc` (`Chooser/mod_lim_proc.f90`) has no caller; add `-checkscaling` | ❌ |
| C1.2 | Fix wiring: `xsect_nloewk_r_is_ns` fills `limval_nlo_is` but the checker reads `FintNLO_ew`; `nloqcd r_is` routed to the wrong reader; `rr_qg` fills no `limval_nnlo` | 🔧 |
| C1.3 | Extend `xsect_selector` beyond RR `ns_ga`: RR gq/qg/aq/qa/ag/ga/4q, all RV, all ONLO | ❌ |
| C1.4 | Make it W-aware (which limits are singular per sector: no lepton-collinear pole on the neutrino; no gamma gamma TC) | 🔧 |
| C1.5 | Output one table per sector (limit, slope, grade); store reference tables for Z, W+, W- | ❌ |
| C1.6 | Run for every sector, Z/W+/W-; record in section 9 | ❌ |
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
| C3.1 | Redo the NLO QCD comparison with `NSC-lnuQCD`: W+ and W-, channels ns, gq, qg, rate and mT. Same PDF, scales, masses, widths and alpha in both codes (write the common input down here) | ❌ |
| C3.2 | Redo the NLO EW comparison with `NSC-lnuEW`: W+ and W-, channels ns, aq, qa. Needed again after the 2026-11-11 fixes of the photon-induced channels (I16a, N2) | ❌ |
| C3.3 | Build the NLO QCD generator for l l' + photon (l nu gamma for W+/W-, l l gamma for Z) in the same way as `NSC-lnuQCD` (Born with a resolved photon, OpenLoops reals, same analysis plus photon cuts) | ❌ |
| C3.4 | Build the NLO EW generator for l l' + jet (l nu j for W+/W-, l l j for Z) in the same way as `NSC-lnuEW` (OpenLoops `pplljj_ew`-type libraries, photon-induced channels, resolved jet) | ❌ |
| C3.5 | Analyses in this code that match C3.3 and C3.4: resolved photon (cf. Z analysis `lla_nloqcd`) and resolved jet, for W+ and W- | ❌ |
| C3.6 | Compare C3.3 with the sum of the O(alpha_s alpha) sectors here with one resolved photon (`rvqcd` + RR + ONLO), channel by channel | ❌ |
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
| nnlo ns rvewk_is | 🔧 | | – | 🔧 | | – |
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
