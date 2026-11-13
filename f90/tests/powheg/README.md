# POWHEG comparison at NLO (PLAN.md C3.1, C3.2)

Codes: `NSC-lnuQCD`, `NSC-lnuEW` in `POWHEG-BOX-RES/MiNNLOPS_res/`, used at fixed order (`testplots 1`, histograms in `pwg-NLO.top`).

## Common setup

13 TeV, PDF `NNPDF31_nnlo_as_0118_luxqed` (325100), alpha_s from the PDF, fixed scale mu = 91.1876 GeV,
G_mu scheme with complex masses: MW = 80.398, GammaW = 2.1054, MZ = 91.1876, GammaZ = 2.4952, G_F = 1.16639e-5,
unit CKM, no b. Cuts of `analysis/mod_cut_histo_W_ATLAS.f90`.

- POWHEG input: `base.input` plus `chooseDY 1|-1` and the statistics.
- POWHEG analysis: `pwhg_analysis-NSC.f` (copy it into the process directory). It mirrors `W_ATLAS`,
  including the dressing of the charged lepton with a photon within R = 0.1. To see the photon in a
  fixed-order run add `detailedNLO 1` (otherwise the radiated parton has id 0).
- mixed_ll: `-ecoll 13000 -dynscale .false. -mu 91.1876`, analysis `W_ATLAS`.

## Build of the POWHEG codes

    export PATH=<lhapdf>/bin:<fastjet>/bin:$PATH
    make ANALYSIS=W_ATLAS PWHGANAL="pwhg_bookhist-multi.o pwhg_analysis-NSC.o" \
         OLPATH=<openloops> LIBSOPENLOOPS="-Wl,-rpath,<openloops>/lib -L<openloops>/lib -lopenloops" \
         LHAPDFCXXFLAGS="-I<lhapdf>/include" pwhg_main

## Run

    ./pwgrun.sh qcd_nlo_wp <path>/NSC-lnuQCD 16 chooseDY=1 ncall1=300000 itmx1=3 ncall2=1000000 itmx2=4
    ./pwgcomb.py qcd_nlo_wp            # average of the 'rate' histogram over the seeds
    ./nscrun.sh out.txt wp "lo ns na" "nloqcd ns v" ...   # from f90/, sums the sectors

## Results (2026-11-11, NLO EW 2026-11-12)

| | POWHEG | mixed_ll | difference |
|---|---|---|---|
| LO W+ | 1417.2 +- 1.5 fb | 1416.41 +- 0.18 fb | -0.7 +- 1.5 |
| LO W- | 990.73 +- 0.47 fb | 990.88 +- 0.12 fb | +0.15 +- 0.49 |
| NLO QCD W+ (ns + gq + qg) | 1737.51 +- 0.61 fb | 1736.68 +- 0.31 fb | -0.83 +- 0.69 |
| NLO QCD W- (ns + gq + qg) | 1203.14 +- 0.44 fb | 1203.02 +- 0.18 fb | -0.12 +- 0.48 |
| NLO EW W+ (ns + aq + qa) | 1327.29 +- 0.25 fb | 1326.99 +- 0.21 fb | -0.30 +- 0.33 |
| NLO EW W- (ns + aq + qa) | 929.24 +- 0.25 fb | 929.23 +- 0.13 fb | -0.01 +- 0.28 |

mixed_ll sectors, W+ / W- in fb. NLO QCD: ns v -283.75 / -198.50, ns s 560.18 / 389.76, ns r_is -28.54 / -17.22,
gq s -3.78 / -5.92, gq r_is 39.94 / 24.96, qg s -3.78 / -5.92, qg r_is 39.98 / 24.98.
NLO EW: ns v -74.30 / -50.38, ns s -5.54 / -3.94, ns r_is -1.94 / -0.97, ns r_fs_53 -2.66 / -5.85,
ns r_fs_54 -7.56 / -1.98, aq s 0.0056 / 0.0088, aq r_is 1.281 / 0.725, qa s 0.0056 / 0.0088, qa r_is 1.282 / 0.724.

## NLO EW: fixes in `NSC-lnuEW` (2026-11-12)

`NSC-lnuEW` as pushed (`eb7c4b3`) gave about 38 pb total against 4.06 pb at LO. The public POWHEG-BOX-RES
(`365bb20`) multiplies the output of `setreal` and `setvirtual` by alpha_s/2pi also in QED-only mode
(`sigreal.f:1331`, `sigsoftvirt.f:394`), while `openloops2.f` divided by alpha/2pi. Three lines of `openloops2.f` changed:
- real and virtual are divided by `st_alpha/(2pi)`;
- `call set_parameter("mu", sqrt(st_muren2))` is active before the loop evaluation (OpenLoops used 100 GeV).

Run with `detailedNLO 1`: without it the fixed-order events carry no photon id and the lepton is not dressed
(W+ undressed: 1305.28 +- 0.76 fb).

    ./pwgrun.sh ew_nlo_wp <path>/NSC-lnuEW 31 chooseDY=1 ncall1=100000 itmx1=3 ncall2=600000 itmx2=6 detailedNLO=1
