#!/usr/bin/env python3
"""All observables of a combined run, LO and LO + dNLO with their scale bands (mu/2, 2 mu), and the ratio to LO, from f90/:
    python3 cluster_KIT/plot_obs.py [nloqcd|nloewk|nnlo] [--analysis W_incl] [--seeds 10]
Reads histo/comb_<analysis>/<analysis>_<b>_lo_sum10.histo and <b>_<corr>_total.histo (cluster_KIT/plot_lo.py,
combine_nlo.py) and writes <corr>_<b>_observables.png there."""
import os, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

opt = lambda k, d: type(d)(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
CORR = next((a for a in sys.argv[1:] if a in ('nloqcd', 'nloewk', 'nnlo')), 'nloqcd'); ANALYSIS = opt('--analysis', 'W_incl')
OUT = 'comb' if ANALYSIS == 'W_ATLAS' else f'comb_{ANALYSIS}'
OBS = {'W_incl': [('rate per $m_T$ window', [2], r'window $\Phi(i)$: 200–300, 300–500, 500–1500, >1500 GeV', True),
                  (r'$m_T$', [3, 4], r'$m_T$ [GeV]', True), (r'$m_{l\nu}$', [5, 6], r'$m_{l\nu}$ [GeV]', True),
                  (r'$p_{T,l}$', [7], r'$p_{T,l}$ [GeV]', True), (r'$y_l$', [8], r'$y_l$', False),
                  (r'$\eta_l$', [11], r'$\eta_l$', False), (r'$p_{T,l\nu}$', [9], r'$p_{T,l\nu}$ [GeV]', True),
                  (r'$\cos\theta^*_{CS}$', [10], r'$\cos\theta^*$', False)]}[ANALYSIS]
TITLE = {'nloqcd': 'NLO QCD', 'nloewk': 'NLO EW', 'nnlo': "NNLO q qbar'"}; BOS = {'wp': 'W$^+$', 'wm': 'W$^-$'}
C = {'lo': '#52514e', 'nlo': '#2a78d6'}; GRID = '#e4e3df'
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'histo', OUT))

def read(f, scale):                             #-- {histo: rows [x, v1, e1, v2, e2, v3, e3]} without the overflow row
    h = {}
    for l in open(f):
        x = l.split()
        if len(x) == 9 and not l.lstrip().startswith('#'): h.setdefault(int(x[0]), []).append([float(v) for v in x[2:]])
    return {k: (np.array(v) * np.array([1] + [scale] * 6))[:-1] for k, v in h.items()}

def obs(h, ids, perwidth):                      #-- edges and values per bin for the scale columns 1, 3, 5, and the MC error
    rows = np.vstack([np.hstack([h[k], np.full((len(h[k]), 1), h[k][1, 0] - h[k][0, 0])]) for k in ids])
    w = rows[:, 7] if perwidth else np.ones(len(rows))
    return np.append(rows[:, 0], rows[-1, 0] + rows[-1, 7]), rows[:, [1, 3, 5]] / w[:, None], rows[:, 2] / w

for b in ('wp', 'wm'):
    lo, dn = read(f'{ANALYSIS}_{b}_lo_sum10.histo', 1 / 10), read(f'{b}_{CORR}_total.histo', 1 / opt('--seeds', 10))
    fig = plt.figure(figsize=(15, 19)); sub = fig.subfigures(4, 2, hspace=0.04, wspace=0.04)
    for n, (name, ids, xl, logy) in enumerate(OBS):
        top, bot = sub[n // 2, n % 2].subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1.4], 'hspace': 0.06})
        perw = ids != [2]
        e, L, Le = obs(lo, ids, perw); _, D, De = obs(dn, ids, perw); N = L + D
        onebin = np.count_nonzero(L[:, 0] > 1e-3 * L[:, 0].max()) < 3          #-- p_T(l nu): LO is a single bin
        for y, col, lab in ((L, C['lo'], 'LO'), (N, C['nlo'], f'LO + δ{TITLE[CORR]}')):
            if logy and (y[:, 0] < 0).any():                                     #-- |y|, dashed where negative
                top.stairs(np.where(y[:, 0] > 0, y[:, 0], np.nan), e, baseline=None, color=col, lw=1.4, label=lab)
                top.stairs(np.where(y[:, 0] < 0, -y[:, 0], np.nan), e, baseline=None, color=col, lw=1.4, ls='--')
                continue
            top.stairs(y[:, 0], e, baseline=None, color=col, lw=1.4, label=lab)
            top.stairs(y.max(1), e, baseline=y.min(1), fill=True, color=col, alpha=0.25, lw=0)
        ref = N[:, 0] if onebin else np.where((L[:, 0] > 0) & (np.abs(N[:, 0] - L[:, 0]) < 2 * L[:, 0]), L[:, 0], np.nan)   #-- no ratio where LO ~ 0
        for y, col in (((N, C['nlo']),) if onebin else ((L, C['lo']), (N, C['nlo']))):
            bot.stairs(y[:, 0] / ref, e, baseline=None, color=col, lw=1.3)
            bot.stairs(y.max(1) / ref, e, baseline=y.min(1) / ref, fill=True, color=col, alpha=0.25, lw=0)
        bot.stairs((N[:, 0] + De) / ref, e, baseline=(N[:, 0] - De) / ref, fill=True, color='#104281', alpha=0.6, lw=0)
        if logy: top.set_yscale('log')
        top.set_title(name, loc='left', fontsize=11); top.set_ylabel('fb' if not perw else 'fb / unit')
        bot.set_ylabel(f'ratio to LO + δ{TITLE[CORR]}' if onebin else 'ratio to LO', fontsize=8); bot.set_xlabel(xl, fontsize=9)
        for a in (top, bot):
            a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True); a.spines[['top', 'right']].set_visible(False); a.set_xlim(e[0], e[-1])
        if n == 0: top.legend(frameon=False, fontsize=9)
    fig.suptitle(f'{BOS[b]}, {ANALYSIS}: LO and LO + δ{TITLE[CORR]}, bands: $\\mu_R = \\mu_F \\in [\\mu/2, 2\\mu]$, '
                 f'$\\mu = m_T/2$; dark band in the ratio: MC error', x=0.02, ha='left', y=1.005, fontsize=12)
    fig.savefig(f'{CORR}_{b}_observables.png', dpi=130, bbox_inches='tight', facecolor='#fcfcfb'); plt.close(fig)
    print(f'{OUT}/{CORR}_{b}_observables.png')
