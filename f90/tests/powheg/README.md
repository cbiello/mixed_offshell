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

## Results (2026-11-11)

| | POWHEG | mixed_ll | difference |
|---|---|---|---|
| LO W+ | 1417.2 +- 1.5 fb | 1416.41 +- 0.18 fb | -0.7 +- 1.5 |
| LO W- | 990.73 +- 0.47 fb | 990.88 +- 0.12 fb | +0.15 +- 0.49 |
| NLO QCD W+ (ns + gq + qg) | 1737.51 +- 0.61 fb | 1736.68 +- 0.31 fb | -0.83 +- 0.69 |
| NLO QCD W- (ns + gq + qg) | 1203.14 +- 0.44 fb | 1203.02 +- 0.18 fb | -0.12 +- 0.48 |
| NLO EW W+ (ns + aq + qa) | not available, see below | 1326.99 +- 0.21 fb | |
| NLO EW W- (ns + aq + qa) | not available, see below | 929.23 +- 0.13 fb | |

mixed_ll sectors, W+ / W- in fb. NLO QCD: ns v -283.75 / -198.50, ns s 560.18 / 389.76, ns r_is -28.54 / -17.22,
gq s -3.78 / -5.92, gq r_is 39.94 / 24.96, qg s -3.78 / -5.92, qg r_is 39.98 / 24.98.
NLO EW: ns v -74.30 / -50.38, ns s -5.54 / -3.94, ns r_is -1.94 / -0.97, ns r_fs_53 -2.66 / -5.85,
ns r_fs_54 -7.56 / -1.98, aq s 0.0056 / 0.0088, aq r_is 1.281 / 0.725, qa s 0.0056 / 0.0088, qa r_is 1.282 / 0.724.

## NLO EW: open

`NSC-lnuEW` as pushed (commit `eb7c4b3`) does not give a usable NLO result with the public POWHEG-BOX-RES
(`365bb20`) and OpenLoops 2.1.5: the total cross section comes out at about 38 pb against 4.06 pb at LO,
with large fluctuations between seeds, also with the original `testrun-lhc` input and with a PDF without photon.
Checked: the Born is right (1417.6 fb), the soft and collinear limit tests of POWHEG pass, the OpenLoops EW
virtual has a normal size (finite part / Born = -29.5 alpha/2pi at the test point). Not understood:
- `sigsoftvirt.f` of the public POWHEG-BOX-RES multiplies the virtual by alpha_s/2pi, while
  `openloops2.f` of `NSC-lnuEW` normalises it by alpha/2pi;
- the analysis used for the April 2026 benchmark (`pwhg_analysis-V_minimal.f`, named in the Makefile) is not in the repository.
The working copy used for that benchmark is needed to finish this comparison.
