#!/usr/bin/env python3
"""ATLAS d sigma/d m_T (arXiv:2502.21088, HEPData ins2895869, e and mu combined, Born level) against LO and LO + NLO QCD
of the W_ATLASfid run, from the f90 directory of that run:
    python3 cluster_KIT/plot_atlas.py [--seeds 10]
Reads histo/comb_W_ATLASfid/W_ATLASfid_<b>_lo_sum10.histo and <b>_nloqcd_total.histo (plot_lo.py, combine_nlo.py),
histogram 11 (mT_ATLAS, bin i = ATLAS bin i), and data/hepdata_ins2895869/lep_physical_{plus,minus}_mtw.csv of the
repository. Data error: statistical and all systematic sources in quadrature (larger of up/down) and the luminosity
(0.83%). Also shown: the ATLAS DYTurbo NNLO QCD prediction with NNPDF4.0 (PDF, alpha_s, scale, EW in quadrature).
Writes histo/comb_W_ATLASfid/atlas_mT.png and atlas_mT.md (table)."""
import os, sys, csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

NS = int(sys.argv[sys.argv.index('--seeds') + 1]) if '--seeds' in sys.argv else 10
here = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(here, '..', '..', 'data', 'hepdata_ins2895869')
OUT = os.path.join(here, '..', 'histo', 'comb_W_ATLASfid')
C = {'lo': '#52514e', 'nlo': '#2a78d6', 'dyt': '#eb6834', 'data': '#0b0b0b'}; GRID = '#e4e3df'
BOS = {'wp': ('plus', 'W$^+ \\to \\ell^+\\nu$'), 'wm': ('minus', 'W$^- \\to \\ell^-\\bar\\nu$')}

def hepdata(name):                              #-- blocks of the csv: [header, rows...]
    blocks, cur = [], []
    for r in csv.reader(l for l in open(os.path.join(DATA, name)) if not l.startswith('#')):
        if r: cur.append(r)
        elif cur: blocks.append(cur); cur = []
    return blocks + [cur] if cur else blocks

def band(rows, cols):                           #-- quadrature of the given (+, -) column pairs, larger side
    return np.sqrt(sum(np.maximum(np.abs(rows[:, c]), np.abs(rows[:, c + 1])) ** 2 for c in cols))

def histo(f, scale):                            #-- histogram 11 without the overflow row: values for mu, mu/2, 2mu and the error
    h = [l.split() for l in open(os.path.join(OUT, f)) if len(l.split()) == 9 and not l.lstrip().startswith('#')]
    a = np.array([[float(v) for v in x[2:]] for x in h if x[0] == '11'])[:-1] * scale
    return a[:, [1, 3, 5]], a[:, 2]

fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), sharex='col', gridspec_kw={'height_ratios': [3, 1.6], 'hspace': 0.06})
table = ['# ATLAS d sigma/d m_T against W_ATLASfid (fb/GeV)', '']
for j, b in enumerate(('wp', 'wm')):
    tag, label = BOS[b]
    blk = hepdata(f'lep_physical_{tag}_mtw.csv')
    d = np.array([[float(v) for v in r] for r in blk[0][1:]])
    lo_e, hi_e = d[:, 1], d[:, 2]; w = hi_e - lo_e; e = np.append(lo_e, hi_e[-1])
    data = d[:, 3] * 1e3; derr = np.hypot(band(d * 1e3, range(4, d.shape[1], 2)), 0.0083 * data)
    nn = next(x for x in blk if 'NNPDF40' in x[0][3]); t = np.array([[float(v) for v in r] for r in nn[1:]])
    dyt, dyterr = t[:, 3] * 1e3, band(t * 1e3, range(6, t.shape[1], 2))
    lo, _ = histo(f'W_ATLASfid_{b}_lo_sum10.histo', 0.1); nlo, nloerr = histo(f'{b}_nloqcd_total.histo', 1 / NS)
    lo, nlo, nloerr = lo / w[:, None], (lo + nlo) / w[:, None], nloerr / w
    top, bot = axes[0, j], axes[1, j]
    for y, col, lab in ((lo, C['lo'], 'LO'), (nlo, C['nlo'], 'LO + NLO QCD (this code)')):
        top.stairs(y[:, 0], e, baseline=None, color=col, lw=1.4, label=lab)
        top.stairs(y.max(1), e, baseline=y.min(1), fill=True, color=col, alpha=0.25, lw=0)
        bot.stairs(y[:, 0] / data, e, baseline=None, color=col, lw=1.3)
        bot.stairs(y.max(1) / data, e, baseline=y.min(1) / data, fill=True, color=col, alpha=0.25, lw=0)
    top.stairs(dyt, e, baseline=None, color=C['dyt'], lw=1.2, ls='--', label='DYTurbo NNLO QCD, NNPDF4.0 (ATLAS)')
    bot.stairs(dyt / data, e, baseline=None, color=C['dyt'], lw=1.1, ls='--')
    bot.stairs((dyt + dyterr) / data, e, baseline=(dyt - dyterr) / data, fill=True, color=C['dyt'], alpha=0.18, lw=0)
    xc = np.sqrt(lo_e * hi_e)                   #-- markers at the geometric bin centre (log x axis)
    top.errorbar(xc, data, yerr=derr, fmt='o', ms=4, color=C['data'], lw=1.1, label='ATLAS 140 fb$^{-1}$ (e + $\\mu$)')
    bot.errorbar(xc, np.ones_like(data), yerr=derr / data, fmt='o', ms=4, color=C['data'], lw=1.1)
    top.set_yscale('log'); top.set_title(label, loc='left', fontsize=12); top.set_ylabel('d$\\sigma$/d$m_T$ [fb/GeV]')
    bot.set_ylabel('ratio to data'); bot.set_xlabel('$m_T$ [GeV]'); bot.set_xscale('log'); bot.set_ylim(0.7, 1.3)
    for a in (top, bot):
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True); a.spines[['top', 'right']].set_visible(False); a.set_xlim(e[0], e[-1])
    if j == 0: top.legend(frameon=False, fontsize=9, loc='lower left')
    table += [f'## {label.replace("$", "")}', '', '| m_T [GeV] | ATLAS | LO | LO + NLO QCD (mu/2, 2mu) | NLO / data | DYTurbo NNLO / data |',
              '|---|---|---|---|---|---|']
    table += [f'| {l:g}–{h:g} | {x:.4g} ± {dx:.2g} | {a[0]:.4g} | {n[0]:.4g} ({n[1]:.4g}, {n[2]:.4g}) | {n[0] / x:.3f} | {y / x:.3f} |'
              for l, h, x, dx, a, n, y in zip(lo_e, hi_e, data, derr, lo, nlo, dyt)]
    sig = lambda v: float(np.sum(v * w))
    table += ['', f'Fiducial total: ATLAS {sig(data):.1f} fb, LO {sig(lo[:, 0]):.1f}, LO + NLO QCD {sig(nlo[:, 0]):.1f} '
              f'({sig(nlo[:, 1]):.1f}, {sig(nlo[:, 2]):.1f}), DYTurbo NNLO {sig(dyt):.1f}', '']
fig.suptitle('ATLAS high-$m_T$ W (arXiv:2502.21088) against LO and NLO QCD, fiducial region of `W_ATLASfid`; bands: '
             '$\\mu_R = \\mu_F \\in [\\mu/2, 2\\mu]$, $\\mu = m_T/2$', x=0.01, ha='left', fontsize=11)
fig.savefig(os.path.join(OUT, 'atlas_mT.png'), dpi=130, bbox_inches='tight', facecolor='#fcfcfb')
open(os.path.join(OUT, 'atlas_mT.md'), 'w').write('\n'.join(table) + '\n')
print(os.path.join(OUT, 'atlas_mT.png'))
