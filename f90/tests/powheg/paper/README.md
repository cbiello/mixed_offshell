# POWHEG benchmarks with the paper setups (SETUP.md section 4)

| | `mixed_ll` | POWHEG |
|---|---|---|
| Setup A (mT > 200) | analysis `W_incl` | `nscsetup 1` |
| Setup B (ATLAS fiducial) | analysis `W_ATLASfid` | `nscsetup 2` |
| NLO QCD, NLO EW | `-require none` (default) | `nscrequire 0`, `NSC-lnuQCD`, `NSC-lnuEW` |
| l nu + photon | `-require photon` | `nscrequire 1`, `NSC-lnuaEW` |
| l nu + jet | `-require jet` | `nscrequire 2`, `NSC-lnujQCD` |

The cuts on the photon and the jet are in the headers of `W_leptons` (`f90/analysis/mod_cut_histo_aux.f90`) and of
`pwhg_analysis-NSCpaper.f`. Histograms have the same names and bins in the two codes.

## Build and run

    # mixed_ll: analysis=W_incl or W_ATLASfid in the cfg
    # POWHEG: copy pwhg_analysis-NSCpaper.f into the process directory, then (paths as in ../README.md)
    make ANALYSIS=W_ATLAS PWHGANAL="pwhg_bookhist-multi.o pwhg_analysis-NSCpaper.o" ... pwhg_main
    mv pwhg_main pwhg_main_paper
    PWG=pwhg_main_paper ../pwgrun.sh <tag> <path>/NSC-lnuQCD 4 chooseDY=1 nscsetup=1 nscrequire=0 ncall1=... 

`base.input` here has the PDF of the paper setup (336100), `detailedNLO 1` and `withdamp 0`.
