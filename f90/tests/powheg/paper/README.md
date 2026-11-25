# POWHEG benchmarks with the paper setups (SETUP.md section 4)

| | `mixed_ll` | POWHEG |
|---|---|---|
| Setup A (mT > 200) | analysis `W_incl` | `nscsetup 1` |
| Setup B (ATLAS fiducial) | analysis `W_ATLASfid` | `nscsetup 2` |
| A and B in one run | – | `nscsetup 3` (histograms of B with the prefix `B_`) |
| NLO QCD, NLO EW | `-require none` (default) | `nscrequire 0`, `NSC-lnuQCD`, `NSC-lnuEW` |
| l nu + photon | `-require photon` | `nscrequire 1`, `NSC-lnuaEW` |
| l nu + jet | `-require jet` | `nscrequire 2`, `NSC-lnujQCD` |

The cuts on the photon and the jet are in the headers of `W_leptons` (`f90/analysis/mod_cut_histo_aux.f90`) and of
`pwhg_analysis-NSCpaper.f`. Histograms have the same names and bins in the two codes.

## Scale

mu_R = mu_F = xi mT/2, mT of the dressed lepton and the neutrino, on the kinematics of each term.
- `mixed_ll`: `-dynscale .true.`; one run gives the three columns xi = 1, 1/2, 2.
- POWHEG: `fixedscale 0`, `btlscalereal 1`, `btlscalect 1` (real term on the real kinematics, counterterms on the
  Born), one run per xi with `renscfact = facscfact = xi`. Needs `Born_phsp_mT_scale.patch` in `NSC-lnuQCD` and
  `NSC-lnuEW` (the public codes only have mu = mZ). The EW code dresses with the radiated parton whatever its flavour:
  for the photon-induced channels this differs from `mixed_ll` when the final quark is within 0.1 of the lepton.

## Build

    # mixed_ll: analysis=W_incl or W_ATLASfid in the cfg
    # POWHEG: copy pwhg_analysis-NSCpaper.f into the process directory, then (paths as in ../README.md)
    make ANALYSIS=W_ATLAS PWHGANAL="pwhg_bookhist-multi.o pwhg_analysis-NSCpaper.o" ... pwhg_main
    mv pwhg_main pwhg_main_paper

## LO, NLO QCD, NLO EW production (2026-11-25)

Scripts (paths of this machine inside): `run_nsc.sh` (one sector, one seed), `run_pwg.sh` (one seed, LO / NLO QCD /
NLO EW at one xi, setups A and B together), `job.sh` (dispatcher for `xargs -P <cores> -L 1 ./job.sh < jobs.txt`).
`base.input` is the common POWHEG input (PDF 336100, `detailedNLO 1`, `withdamp 0`). `mixed_ll` runs with the default
technical mass rescaling, as in the NLO benchmarks of `../README.md`.

    python3 benchmark.py <run directory>      # averages the seeds, sums the sectors, compares: benchmark/rates.md,
                                              # benchmark/<setup>_<boson>.dat, plots mixed_ll / POWHEG per histogram

Notes: `-histodir` of `mixed_ll` holds 30 characters; `mixed_ll` returns exit code 0 also when it stops on an
unknown or repeated option (check the log for `# Result`); POWHEG stops on a repeated key in `powheg.input`.
