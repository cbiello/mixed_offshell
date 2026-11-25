# AI-plotting — paper plots with `NSCplots/`

`NSCplots/` turns `mixed_ll` (NSC) and POWHEG (PWG) histograms into gnuplot files and PDF plots:
main frame with the curves, lower frame with the ratio to the first curve and the scale band.
Setup of the runs: [SETUP.md](SETUP.md). How the histograms are produced and combined: [PLAN.md](PLAN.md) section 11.

## Chain

| Step | Input | Script | Output |
|---|---|---|---|
| 1 | `mixed_ll` histogram | `fromNSCtoPWG/fromHISTOtoTOP.py` (single run) or `fromavHISTOtoTOP.py` (summed seeds) | `output_CEN.top`, `output_MIN.top`, `output_MAX.top` |
| 2 | three top files (central and two scale variations) | `topfiles/convert_top_files_to_plots.py` | `<name>-run/distributions/<obs>__<name>.dat` |
| 3 | two sets of `.dat` files (`...-NSC`, `...-PWG`) | `LOscript-*.py`, `NLOQCDscript-*.py` | `gnuplot_minnlo/<obs>__<name>.gnu`, `.pdf`, `histograms/` |

POWHEG gives top files directly (`11.top` central, `HH.top` mu/2, `22.top` 2 mu): it enters at step 2.

## Layout

| Path | Content |
|---|---|
| `fromNSCtoPWG/` | Step 1 converters, with a sample `input.histo` |
| `topfiles/<order>-<boson>[-NSC or -PWG-<date>]/` | Top files of each run; `convert_top_files_to_plots.py` |
| `datfiles/<boson>-<order>-<code>-run/distributions/` | Step 2 output: columns x (left edge), central, error, min, error, max, error; fb per unit of x |
| `LOscript-Wm.py`, `LOscript-Z.py`, `NLOQCDscript-Wm.py`, `-Wp.py`, `-Z.py` | Step 3. Five copies of one 2500-line script (class `gnuplot`, then the main block) |
| `paper/<set>/` | Paper plots, one directory per set (see below) |

The five scripts differ only in the main block:

| Line | What |
|---|---|
| ~1116 | glob of the NSC `.dat` files; the PWG file is the same path with `-NSC` -> `-PWG` |
| ~1195 | `gnu.add_curve(path, {label, line_style})`; the first curve is the ratio reference |
| ~1227 | range of the ratio frame |
| ~1235-1280 | per observable (`mass`, `massT`, `pT-lep`, `pT-nu`, `y-lep`): labels, units, x range, rebinning |
| ~2441 | `process`, `collider`, `energy` labels |

## Paper plots

Work in a new directory `NSCplots/paper/<set>/` (e.g. `Wp-NLOQCD`), one per set of plots: the scripts read
`datfiles/` and write `gnuplot_minnlo/` relative to the current directory, and delete `gnuplot_minnlo/` at start.
Do not write paper plots into the top-level `datfiles/`, `topfiles/`: those are the validation runs.

    NSC=$(git rev-parse --show-toplevel)/NSCplots
    W=$NSC/paper/Wp-NLOQCD
    mkdir -p $W/topfiles/NSC $W/topfiles/PWG $W/datfiles

    # 1. mixed_ll histogram -> top files (input name is fixed: input.histo, or input.dat for fromavHISTOtoTOP.py)
    cd $W/topfiles/NSC
    cp <histogram> input.histo
    python3 $NSC/fromNSCtoPWG/fromHISTOtoTOP.py
    cp <POWHEG run>/{11,HH,22}.top $W/topfiles/PWG/

    # 2. top -> dat (arguments: central, min, max, name; writes <name>-run/ in the current directory)
    cd $W/datfiles
    python3 $NSC/topfiles/convert_top_files_to_plots.py ../topfiles/NSC/output_{CEN,MIN,MAX}.top Wp-NLOQCD-NSC
    python3 $NSC/topfiles/convert_top_files_to_plots.py ../topfiles/PWG/{11,HH,22}.top Wp-NLOQCD-PWG

    # 3. dat -> gnu and pdf
    cd $W
    cp $NSC/NLOQCDscript-Wp.py plot.py        # then edit the lines of the table above
    python3 plot.py
    ls gnuplot_minnlo/*.pdf

To change one plot by hand, edit its `.gnu` and run `gnuplot <file>.gnu`. Commit `plot.py`, `datfiles/`, and the
`.gnu` and `.pdf` files of the set.

For a plot without POWHEG (e.g. LO, NLO, NNLO from `mixed_ll`): make one `.dat` set per curve with step 1 and 2,
and add one `gnu.add_curve` per set in `plot.py`.

## Things to know

- Needs `gnuplot`, and the Python modules `pikepdf` and `scipy` (imported, not used).
- Step 1 divides by 1000 x bin width (fb per bin -> pb per unit of x), hard-coded for the `W_ATLAS` binning
  (10 GeV for `massT`, `mass`; 5 GeV for the pT; 0.2 for `y-lep`), and expects the observables in the order
  `rate, massT, pT-lep, pT-nu, mass, y-lep`. Update both when the analysis changes (the binning of
  [SETUP.md](SETUP.md) is different).
- MIN and MAX are the two scale variations in the order of the files (NSC: mu/2 then 2 mu), not a sorted envelope.
  The stored `*-PWG` files were made with the order `11 22 HH`.
- The observable names (`# <name> index n` in the top files) must be the same in the two codes.
- The scripts start gnuplot without waiting: the PDFs appear a moment after the script ends.
- The `.gnu` files contain absolute paths: regenerate them after moving a directory.
- The energy label is hard-coded to 14 TeV; the runs of [SETUP.md](SETUP.md) are at 13 TeV.
- `NLOQCDscript-Wp.py` reads `datfiles/Wp-NLOQCD-NSC-run`, which is not there yet.
