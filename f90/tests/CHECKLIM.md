# Limit checks (`-checklim`): how to run and read them

Build with `-D_withchecks=1` in `config.cfg` (the default here), then `make`. All commands from `f90/` (`tests/checklim_run.sh` also works from any directory).
The analysis compiled in (`-DANALYSIS`) sets the cuts: with `W_ATLAS` only W+ and W- can be checked; Z needs a build
with a Z analysis.

## 1. One point, all limits side by side

    yes '' | ./mixed_ll -boson wp -corr nnlo -ch gq -sec rr_5161c -checklim | grep -av 'overriding input\|limit depth'

Any NLO real (`-corr nloqcd|nloewk -sec r_*`), RV (`rv*`), RR (`rr_*`) or S sector with a resolved parton (`s_o*`).
For a random point that passes the cuts it prints one table per limit, then waits for Enter (Ctrl-C to stop).

- `limit S6TC (x(E6) -> 0, x(rho5) -> 0)`: the variables of that limit are set to 1e-1 ... 1e-10, the others are kept.
- Columns: H and every counterterm, with the sign of H, so the terms that must cancel show the same number. In a limit
  L every term T must approach T+L: soft S6, H = S6 and TC = S6TC; triple collinear TC, H = TC and S6 = S6TC; ...
  RR names: S5, S6 soft; TC triple collinear; C5, C6, C56 collinear. NLO/RV/ONLO: S, C1, C2, SC1, SC2, `/g` = kinematic group.
- `limits without own counterterm`: x^k |H + counterterms| at 1e-2 ... 1e-8 (k = variables at the edge). Must fall:
  `integrable`. `SINGULAR` means a counterterm is missing. `cut` means the point fails the cuts there.
  For W this covers the photon collinear to the neutrino and the gamma gamma -> l l triple-collinear limit.

## 2. Many points, graded (PLAN C1.5)

    tests/checklim_run.sh wp 200 csv            # all 94 sectors, 200 points each: about 5 minutes on 4 cores
    tests/checklim_run.sh wm 200 csv
    python3 tests/checklim_analyse.py csv ref   # a few minutes

`checklim_run.sh` runs every sector with `-checklim_csv -checklim_npts 200`: one line per point, limit and pair
(the ratios A/B at the ten depths), no pause. One sector alone:
`./mixed_ll -boson wp -corr nnlo -ch gq -sec rr_5161c -checklim_csv -checklim_npts 200 </dev/null | grep '^[PRU],' > x.csv`

For each point and pair, r(d) = |A/B - 1|. Its slope in log r against log d (between 1e-2 and 1e-6, or between 1e-4
and 1e-8) must be at least the expected one minus 0.2: expected 1 for a soft limit, 1/2 for any limit with a collinear
variable (faster is fine). A limit without own counterterm passes if x^k |H + counterterms| somewhere below 1e-3 is
3 times smaller than at 1e-2 ... 1e-4; a point where H is cut at 1e-10 is at a cut edge and not graded.

Output in `ref/`:

- `<boson>/summary.csv`: one row per sector, limit and pair: points, fraction passing, expected and median slope.
  Rows below 0.95 are printed at the end of the run.
- `<boson>/worst.csv`: the failing points (their number `n`) of each row.
- `<boson>/checklim_<boson>.pdf`: page 1 the heatmap sectors x limits (worst pair, green >= 95% pass, yellow 80-95%,
  red below); then one page per sector: median and 68% band of |A/B - 1| against depth for each pair, with the
  expected d^p line, and the histogram of the per-point slopes with the threshold.
- `dash/*.json`: the same, for the dashboard page.

## 3. Rerun a failing point as tables

    grep '^P,17,' csv/wp.nnlo.gq.rr_5161c.csv | cut -d, -f3- | tr ',' ' ' > p17.txt
    ./mixed_ll -boson wp -corr nnlo -ch gq -sec rr_5161c -checklim -checklim_point p17.txt | grep -av 'overriding input'

## Reading the results

- Good: soft rows slope ~1.0, collinear and mixed rows ~0.5, fraction passing ~1.
- A low fraction with a median slope at the expected value: a few points where A - B changes sign near 1e-3, or where
  double precision runs out (depths 1e-8 and below, earlier when three variables go to the edge together). Rerun them (3).
- A median slope near 0: the counterterm does not reproduce the limit (wrong charge, missing term). This is a bug.
- `SINGULAR` for a limit without counterterm: a missing counterterm. A single flag at a triple limit is usually precision.
