# NNLO QCD code for W+ (imported)

Imported unchanged on 2026-11-27, to compute the NNLO QCD corrections that are combined with the mixed QCD x EW
corrections of `f90/` for the comparison with the ATLAS data (PLAN.md).

- Source: `git@git.particle.kit.edu:roentsch/f90_wp.git`, branch `master`
- Commit: `e033944fa444ef62ab951c70d2331e416ce6d815` (2021-02-13, "removed seed setting in mod_chooser to allow input
  seed"), the last of 22 commits (2021-01-06 to 2021-02-13)
- Author: Raoul Röntsch; used with the authors' permission
- Imported with `git archive e033944`: every file of that commit, nothing added or removed except this file. Our
  changes follow as separate commits, so `git diff <import commit> -- nnlo_qcd/` shows all of them.

What it is (from a first reading): NNLO QCD for pp -> W+ -> e+ nu with the nested soft-collinear subtraction, the same
framework as `f90/` (LO, NLO, NNLO RR/RV/integrated sectors in the ns, qg, gq, gg channels, limit checks), aligned with
the mixed QCD x EW code of arXiv:2009.10386. As imported: on-shell W (`_PhaseSpaceType=6` in `config_wp.cfg`), W+ only,
M_W = 80.398 GeV, Gamma_W = 2.1054 GeV, G_F input, NNPDF3.1 luxQED, analysis `Histograms/cut_histo_w.f90`, `makefile`
set up for the authors' machines (`identity.mk`).
