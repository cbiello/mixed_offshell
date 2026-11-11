# STRUCTURE — how the code is organised

Snapshot of branch `AI-era` (base commit `2edaf8e` + build fixes), written from a full
read of `f90/`. Companion to [PLAN.md](PLAN.md), which lists what is left for W+ and W-.
Line numbers refer to this snapshot and will drift.

## 1. What the code computes

`mixed_ll` computes fully differential off-shell Drell-Yan cross sections up to mixed
O(alpha_s alpha) with the nested soft-collinear subtraction scheme:

| `-boson` (`Vcharge`) | Process | Leptons (legs 3,4) | Charges `Q3,Q4` |
|---|---|---|---|
| `z` (0)  | pp -> gamma*/Z -> e- e+ | `[id_el,-id_el]`  | -1, +1 |
| `wm` (-1) | pp -> W- -> e- nubar    | `[id_el,-id_nue]` | -1, 0 |
| `wp` (+1) | pp -> W+ -> nu e+       | `[id_nue,-id_el]` | 0, +1 |

The boson is chosen **at run time** with `-boson z|wm|wp` (default `z`). It sets the integer
`Vcharge` in `Parms/mod_parms.f90` (and `lep_ids`, `Q3`, `Q4`); the code branches on
`if (Vcharge == ...)`. One binary serves all three bosons. Only the analysis is still fixed
at compile time. Output files of W runs carry a `_wm`/`_wp` tag after the analysis name.
Momentum convention everywhere: `q(1) q'(2) -> l(3) l'(4) [+ parton 5] [+ parton 6]`.

## 2. Repository layout

| Path | Content |
|---|---|
| `f90/` | The Fortran code, build system, run scripts |
| `f90/src/` | Sources (section 4) |
| `f90/analysis/` | Cuts and histograms; one `mod_cut_histo_<name>.f90` per analysis, chosen by `analysis=` in the cfg |
| `f90/aux/`, `f90/bin/` | Histogram combiners (`nsc-hsummer`, `nsc-combine_histos`, `rrsummer`, `cbsummer`) |
| `f90/cluster_run/`, `cluster_KIT/`, `cluster_results/` | Batch scripts and old logs |
| `f90/sectors_list.dat` | Static list of all `-corr/-ch/-sec` combinations (139 lines; not fully in sync with the chooser) |
| `math/` | `DY0_qq.wl`, `DY+_qq.wl`: analytic finite remainders (ns channel); `validate_DYp7.wl`, `run_validate.wls`: symbolic Fortran-vs-notebook check; `tex/DYp_implementation.tex`: term-by-term map |
| `notes/` | Derivation notes (Chiara/Federica: aq, ag, 4q channels; Fabrizio: Z subtraction notebook and pheno numbers; phase space) |
| `old_code/` | Reference Z-only code from Chiara, Fabrizio, Raoul |
| `2L/` | Two-loop amplitude material (Z) |
| `Analysis/`, `plots/`, `setup/` | Misc. older analyses, gnuplot, setup note |

## 3. Build and run

```
cd f90
./configure --cfgfile=cfgcb.txt     # writes config.cfg, makefile.inc, analysis/mod_cut_histo.f90
make remove && make -j16            # clean rebuild only when the analysis or flags change
./mixed_ll -boson <z|wm|wp> -corr <c> -ch <ch> -sec <s> [options]   # one sector per run
./mixed_ll -h                       # all options
```

- Dependencies: LHAPDF, HOPPET, handyG, OpenLoops with process libraries
  `pplla ppllj ppllj_ew pplljj_qcdew` (the W channels live in the same `ppll*` libraries;
  `ppllj` is required by `nloqcd v` but missing from `README.md`).
- Each run integrates one sector with Vegas and writes `histo/<analysis>_<corr>_<ch>_<sec>_...histo`.
  A physical result is the sum of all sectors, combined with the tools in `bin/`.
- Compile flags: `_oldcode` (section 8), `_withchecks` (extra limit bookkeeping), `_MPIbuild`, `_with2L` (two-loop Z amplitudes).

## 4. Source tree (`f90/src`)

| Dir | Modules | Role |
|---|---|---|
| `BasicDef` | `mod_types`, `mod_consts_dp`, `mod_mpi_common`, `mod_limvals` | kinds, constants, MPI state, storage for limit checks |
| `Parser` | `mod_parser` | command line / input file |
| `Parms` | `mod_parms`, `mod_proc_parms`, `mod_vegas_parms` | EW/QCD parameters, charges (`Q_IS(-5:7)`, `Q3`, `Q4`), flavour tests `is_W_pair`, `is_W_quark`, lepton ids `lep_ids`, CKM (unit), cuts, buffers |
| `GenProc` | `mod_process` | `type KinConfig` (below) |
| `AuxFunctions` | `mod_auxfunctions`, `mod_kinfunctions`, `mod_lumi`, `cpolylog_main` | spinors, polylogs, PDF calls |
| `HoppetTools` | `mod_hoppet_tools` | PDF tables (photon = flavour 7 internally) |
| `OpenLoops` | `mod_ol_interface` | registers the OpenLoops processes needed by the selected `(corr,ch,sec)` |
| `Amplitudes` | `mod_coupl`, `mod_amplitudes_tree_ppll`, `mod_amplitudes_loop_ppll`, `mod_amplitudes_twol_nf_ppll` | matrix elements (section 6) |
| `Kinematics` | `mod_kinematics_lo(_z)`, `_nlo(_z)`, `_lept`, `_gen`, `_nnlo_tc_is_eu_ac/abcd`, `_nnlo_dc_ii`, `_nnlo_dc_if`, `mod_aux_kinematics` | phase space + momenta in every singular limit. No `_Vcharge` dependence |
| `SingLimits` | `mod_splittings_bare`, `mod_eikonals`, `mod_limits_rv` | splitting functions, QED/QCD eikonals (`get_qed_eik_gen` is charge-general), one-loop collinear limits (Z-specific charges) |
| `IntSub` | `mod_hoppet_nlo`, `mod_hoppet_nnlo`, `mod_int_sub_nnlo` | integrated subtraction: convolution tables `xPij*` and soft/elastic functions `calG*`, `fin_*` |
| `Sectors` | `mod_xsects_*`, `mod_aux_sectors`, `mod_partitions` | **the integrands** (section 5), luminosity helpers, damping/partition functions |
| `Histograms` | `mod_histo`, `mod_proc_histo` | binning, 3 scale weights |
| `Vegas` | `mod_my_vegas`, `mod_vegas_mpi`, ... | integrator |
| `Chooser` | `mod_run_proc*`, `mod_sects_list`, `mod_lim_proc*` | maps `(corr,ch,sec)` to an integrand |
| `Checks` | `mod_checks`, `mod_check_amp`, `mod_check_lim`, `MG_output/` | amplitude and limit checks (section 9) |
| top | `main.f90`, `mod_chooser.f90` | init (parser, LHAPDF, HOPPET, OpenLoops, histos), then `run_proc` |

### `KinConfig` (`GenProc/mod_process.f90`)
One instance per kinematic configuration (hard process and every limit):
`AmpMom` (momenta for amplitudes), `LimMom`, `Lim_etaij`, `Lim_Ei`, `Lim_sij`, `Lim_z`,
`PartFrac`, `wgt`, `mu2ref`, `npart`, `flag`/`makecut`,
`ids(6)` (what the **analysis** sees: 0 for initial state, lepton ids from `lep_ids`, `id_g`, `id_a`, quark).
`ids` drives cuts, the lepton/neutrino assignment and photon recombination, so it must
be set in every configuration.

## 5. Sectors: naming and dispatch

A run is selected by `-corr`, `-ch`, `-sec`. Dispatch: `Chooser/mod_run_proc.f90` and its
`mod_run_proc_<x>.f90` helpers.

**Channels (`-ch`)**: initial state. `ns` = q qbar(') ; `gq`, `qg`; `aq`, `qa` (a = photon);
`ag`, `ga`; `aa`; RR sub-channels of ns: `ns_ga` (g + gamma emitted), `ns_qqb`, `ns_qq`
(four-quark, singular), `ns_qqb_w`, `ns_qqp_w`, `ns_qqpb_w` (four-quark, finite, OpenLoops).

**Sectors (`-sec`)**

| Order | `-sec` | Meaning | File |
|---|---|---|---|
| LO | `na` | Born | `mod_xsects_lo` |
| NLO | `r_is` | real, emission off initial state, minus local CTs | `mod_xsects_nloqcd_r`, `mod_xsects_nloewk_r` |
| NLO EW | `r_fs_53`, `r_fs_54` | real photon in the partition collinear to lepton 3 / 4 | `mod_xsects_nloewk_r` |
| NLO | `v`, `s` | virtual; integrated subtraction + PDF collinear | `mod_xsects_nlo*_vs` |
| NNLO VV | `vvnf`, `vvfc`, `vvxf`, `vvnf_full` | double virtual pieces (Z only, out of scope for now) | `mod_xsects_nnlo_vv` |
| NNLO RV | `rvewk_is` | EW one-loop x real gluon | `mod_xsects_nnlo_rvewk` |
| | `rvqcd_is`, `rvqcd_fs_53/54` | QCD one-loop x real photon | `mod_xsects_nnlo_rvqcd` |
| | `rvewknf_is` | nf part of EW loop x real gluon | `mod_xsects_nnlo_rvewk` |
| NNLO RR | `rr_5i6j[x]` | double real, see below | `mod_xsects_nnlo_rr_<ch>` |
| NNLO S | `s_12` | double-unresolved integrated terms on Born kinematics | `mod_xsects_nnlo_s_<ch>` |
| | `s_vqcd`, `s_vewk`, `s_vewknf` | integrated single-unresolved x one-loop | idem |
| | `s_oqcd` | "ONLO": resolved gluon/quark + integrated QED subtraction, with its own local CTs | idem |
| | `s_oewk_is`, `s_oewk_fs_53/54`, `s_oewk` | ONLO: resolved photon + integrated QCD subtraction | idem |
| | `s_qqb`, `s_qq` | integrated four-quark TC terms | `mod_xsects_nnlo_s_ns` |

**RR sector labels `5i6j`**: parton 5 is assigned to emitter `i`, parton 6 to emitter `j`
(1,2 = beams; 3,4 = leptons).
- `5161`, `5262`: both on the same beam leg, triple-collinear parametrisation;
  `a` = 6 more collinear, `c` = 5 more collinear, `b`,`d` = 5||6 (only when a 5||6 pole exists).
- `5162`, `5261`: double collinear, one parton on each beam.
- `5163`, `5164`, `5263`, `5264`: gluon/quark on a beam, photon collinear to lepton 3 or 4.
- `partition_nnlo_fact` splits the photon over 4 legs and the gluon over 2. It is
  charge-blind: for W the neutrino still gets a partition, whose collinear CTs vanish
  through `Q=0` but whose hard and soft pieces must still be run.

## 6. Matrix elements

Two generations coexist.

**Legacy (Z only)**: return `(2,2)` arrays `[q qbar / qbar q] x [dn, up]`, hard-wired
`Qdn/Qup/Q_lep`. E.g. `res_tree_qqb`, `res_tree_g_qqb`, `res_tree_ga_qqb`, `res_ewkloop_g_qqb`.

**General (`*_gen`, "New from Raoul", used for NC, CC+, CC-)**: return a full flavour matrix
`res(-5:7,-5:7)` (index 0 = gluon, 7 = photon). For W the entries with
`is_W_pair(i,j)` (charges sum to the W charge, same generation) are filled (unit CKM, no b).

| Process | Tree (`mod_amplitudes_tree_ppll`) | One-loop (`mod_amplitudes_loop_ppll`) |
|---|---|---|
| q qb' -> l l' | `res_tree_qqb_gen` | `res_qcdloop_qqb_gen`, `ol_res_qcdloop_qqb_gen`, `res_ewkloop_qqb_gen` |
| + g (ns, gq, qg) | `res_tree_g_qqb_gen`, `_g_gq_gen`, `_g_qg_gen` | `res_ewkloop_g_qqb_gen`, `_g_gq_gen`, `_g_qg_gen` (+ `_nf_gen`) |
| + gamma (ns, aq, qa) | `res_tree_a_qqb_gen`, `_a_aq_gen`, `_a_qa_gen` (W: MCFM `agamtree`, includes emission off the W) | `ol_res_qcdloop_a_qqb_gen`, `_a_aq_gen`, `_a_qa_gen` |
| + g gamma (ns, gq, qg, aq, qa, ag, ga) | `res_tree_ga_*_gen` (W: `ubdgmsq`) | – |
| four quarks | legacy `res_tree_qqb_qqb`, `res_tree_qq_qq`; OpenLoops for `ns_*_w` (Z only) | – |
| gamma gamma -> l l | `res_treeAA_aa(_gen)` (zero for W) | `res_ewkloopAA_aa` |
| two-loop | – | `mod_amplitudes_twol_nf_ppll`, `TwoLoop/` (Z only) |


Loop amplitudes come from OpenLoops. `mod_ol_interface.f90` registers per `(corr,ch,sec)`
and per `Vcharge`; process ids are stored in `OL_id(:)` with positions that differ between
Z and W.

## 7. Luminosities and charges (`Sectors/mod_aux_sectors.f90`)

| Legacy | General |
|---|---|
| `get_respdf(lumi,...)` with `ns_lumi`, `gq_lumi`, `aq_lumi`, `ag_lumi`, ... on `(2,2)` arrays | `get_respdf_gen`, `get_respdf_vect_gen`, `get_respdf_gen_mu`, `get_respdf_hoppet_gen` -> `gen_lumi`: sum of `f1(i) f2(j) res(i,j)` |
| `Qdn2`, `Qup2`, `Q_lep2`, module `charges` tables | `multiply_IS_charges_sq(res,leg)`, `multiply_FSQ_charges_sq`, `Qsq_FS=[Q3**2,Q4**2]`, `get_qed_eik_gen` with `Qvec=[-Q_IS(al),-Q_IS(be),Q3,Q4,...]` |
| – | `transition('g -> q' | 'ga -> q' | ...)` moves a flavour matrix to another initial state for collinear limits |

Soft and collinear counterterms need only the charges of the four external fermions; charge
conservation takes care of the W. Emission off the W enters only through the hard matrix
elements.

## 8. The `oldcode` switch

Generalised integrands contain both implementations, selected by the preprocessor:

```
#if (_oldcode == 1)
   ... legacy Z implementation ...
#else
   ... general implementation (Z, W-, W+) ...
#endif
```

`_oldcode` is one global compile-time flag: `OLDCODE=true|false` in the cfg file (default
`false`), written by `configure` as `-D_oldcode=1|0` into `config.cfg`. All channels are
therefore always on the same side. An old-code build refuses `-boson wm/wp`.
For `-boson z` the two builds must agree; this is the regression test for the general code.
Files with no `_oldcode` and no `Vcharge` branch are untouched Z code.

## 9. Check infrastructure

| Tool | How | State |
|---|---|---|
| `-MGcheck` | `Checks/mod_check_amp.f90`: `*_gen` tree amplitudes vs stored MadGraph points (`MG_output/`, Z, W-, W+) | works, tree level only, no four-quark |
| `-checklim -lim <x>` | `mod_xsects_checklim_*`: one phase-space point pushed to depth `x`, prints Hard, CT and ratio, then pauses. Needs `_withchecks=1` | NLO sectors, RV sectors, RR only `ns_qqb/ns_qq`; not W-aware |
| scaling checker | `Checks/mod_check_lim.f90`: scans depths 1e-4..1e-15, fits the slope, grades each limit | NLO (3 sectors) and RR `ns_ga`; **unreachable**: `lim_proc` has no caller |
| `-check` | `mod_checks.f90`, `icheck=6` hard-coded -> only the nf amplitude check | Z |
| `-checkpoles` | pole cancellation in legacy loop routines | Z only |
| old vs new code | two builds (`OLDCODE=true/false`), `-boson z`, same seed (section 8) | manual |
| `math/run_validate.wls` | symbolic check of W+ `fin_ns_v7_cc`, `genGnloQCD_cc`, `legfinalcalG_cc` | W+, ns `s_12` term [7] only |

## 10. Analyses

`analysis/mod_cut_histo_<name>.f90` provides `init_user_histo` and `cut_histo(event)`.
W analyses: `W_ATLAS` (high-mT ATLAS-like cuts), `W_ATLASNLOEW`, `W_ATLASNLOQCD`
(POWHEG benchmark setups), `w_dflt`. `mod_cut_histo_aux` has photon recombination
(`recombine_photons`), driven by `event%ids`.
