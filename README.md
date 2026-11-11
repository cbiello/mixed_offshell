# mixed_ll and mixed_nul

## External dependencies
In order to properly install and run `mixed_ll` the following packages must be available:

+ [LHAPDF][2]
+ [OpenLoops][3]
+ [HOPPET][4]
+ [handyG][5]

## OpenLoops installation and configuration
You can download `OpenLoops` from git via

    git clone https://gitlab.com/openloops/OpenLoops.git

or, if you are interested in the latest beta version, via

    git clone -b public_beta https://gitlab.com/openloops/OpenLoops.git

and then compile it using the following procedure:

    cd OpenLoops
    ./scons

In order to be able to access all the process libraries needed in this project you need to create an 'openloop.cfg' file adding the following:

    [OpenLoops]
    process_repositories = public, public_beta, nested_soft_coll

You can download and install all the relevant process libraries going into your `OpenLoops` folder and then typing

    ./openloops libinstall pplla
    ./openloops libinstall ppllj
    ./openloops libinstall ppllj_ew
    ./openloops libinstall pplljj_qcdew

## Build

The build is configured by a cfg file in `f90/` (one per user: `cfgcb.txt`, `cfgchiara.txt`, `cfgchristian.txt`, `cfgraoul.txt`). Copy one and adapt it.

| Key | Meaning |
|---|---|
| `FC`, `FCmpi`, `FFLAGS` | compilers and flags |
| `LHAPDF`, `OPENLOOPS`, `HOPPET`, `HANDYG` | paths of the dependencies |
| `analysis` | analysis module `analysis/mod_cut_histo_<name>.f90` (e.g. `dflt` for Z, `W_ATLAS` for W) |
| `MPI` | `true` / `false` |
| `CHECKS` | `true` enables the limit checks (`-checklim`) |
| `With2L` | `true` links the two-loop amplitudes |
| `OLDCODE` | `true` builds the legacy Z-only code, default `false` |

    cd f90
    ./configure --cfgfile=cfgcb.txt
    make remove
    make -j

## Run

The boson is chosen at run time with `-boson z|wm|wp` (default `z`); the analysis must match it.

    ./mixed_ll -boson wp -corr lo -ch ns -sec na
    ./mixed_ll -boson wm -corr nloewk -ch aq -sec r_is

`./mixed_ll -h` lists all options. The code layout is described in `STRUCTURE.md`, the status and open tasks in `PLAN.md`.

[2]: https://lhapdf.hepforge.org/ "LHAPDF"
[3]: https://openloops.hepforge.org/ "OpenLoops"
[4]: https://github.com/hoppet-code/hoppet "HOPPET"
[5]: https://gitlab.com/mule-tools/handyg "handyG"
