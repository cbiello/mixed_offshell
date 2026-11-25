#!/usr/bin/env python3
"""LO m_T for W+ and W-, from f90/: python3 cluster_KIT/plot_lo.py [--analysis W_ATLAS|W_incl|W_ATLASfid] [--seeds 10]
Sums the seeds 1..n of the cluster run with nsc-hsummer ad3 into histo/comb[_<analysis>]/<analysis>_<b>_lo_sum<n>.histo
(divided by n here; s00000 = warm-up, left out) and plots lo_mT.png there. Columns per bin: mu, mu/2, 2mu (muR = muF)."""
import os, subprocess, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
opt = lambda k, d: type(d)(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
NS, ANALYSIS = opt('--seeds', 10), opt('--analysis', 'W_ATLAS')
MT, XMAX = {'W_ATLAS': ([2], 3000), 'W_incl': ([3, 4], 3000), 'W_ATLASfid': ([2], 5000)}[ANALYSIS]   #-- m_T histograms
OUT = 'comb' if ANALYSIS == 'W_ATLAS' else f'comb_{ANALYSIS}'
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'histo')); os.makedirs(OUT, exist_ok=True)
for b in ('wp', 'wm'):
    subprocess.run(['../bin/nsc-hsummer', 'ad3', *[f'{ANALYSIS}_{b}_lo_ns_na_dynscale_xMuR_1.00_xMuF_1.00_s{k:05d}.histo'
                    for k in range(1, NS + 1)], f'{OUT}/{ANALYSIS}_{b}_lo_sum{NS}.histo'], check=True, stdout=subprocess.DEVNULL)
def load(f):                                    #-- rate, and m_T rows [x, v1, e1, v2, e2, v3, e3] in fb/GeV, bin edges
    r = [l.split() for l in open(f) if l.strip() and not l.lstrip().startswith('#')]
    h = {}
    for x in r: h.setdefault(int(x[0]), []).append([float(v) for v in x[2:9]])
    h = {k: np.array(v) * np.array([1] + [1 / NS] * 6) for k, v in h.items()}
    m = np.vstack([np.hstack([h[k][:-1], np.full((len(h[k]) - 1, 1), h[k][1, 0] - h[k][0, 0])]) for k in MT])
    m = m[(m[:, 0] >= 200) & (m[:, 0] < XMAX)]
    m[:, 1:7] /= m[:, 7:8]
    return h[1][0], m[:, :7], np.append(m[:, 0], m[-1, 0] + m[-1, 7])
d = {b: load(f'{OUT}/{ANALYSIS}_{b}_lo_sum{NS}.histo') for b in ('wp', 'wm')}
C = {'wp': '#2a78d6', 'wm': '#eb6834'}; L = {'wp': 'W$^+$', 'wm': 'W$^-$'}
ink, ink2, grid = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': ink2, 'axes.labelcolor': ink, 'xtick.color': ink2, 'ytick.color': ink2})
fig, ax = plt.subplots(4, 1, figsize=(7.5, 10), sharex=True, gridspec_kw={'height_ratios': [3, 1.4, 1.3, 1.3], 'hspace': 0.08})
for b in ('wp', 'wm'):
    rate, m, edges = d[b]
    y, e = m[:, 1], m[:, 2]; lo, hi = np.minimum(m[:, 3], m[:, 5]), np.maximum(m[:, 3], m[:, 5])
    rlo, rhi = min(rate[3], rate[5]) - rate[1], max(rate[3], rate[5]) - rate[1]
    ax[0].stairs(y, edges, color=C[b], lw=1.6, baseline=None, label=f"{L[b]}   $\\sigma$ = {rate[1]:.1f} $^{{+{rhi:.1f}}}_{{{rlo:.1f}}}$ (scale) $\\pm$ {rate[2]:.2f} (MC) fb")
    ax[0].stairs(hi, edges, baseline=lo, fill=True, color=C[b], alpha=0.25, lw=0)
    ax[1].stairs(m[:, 3] / m[:, 1], edges, color=C[b], lw=1.4, ls='--', baseline=None); ax[1].stairs(m[:, 5] / m[:, 1], edges, color=C[b], lw=1.4, baseline=None)
    ax[3].stairs(np.where(y > 0, e / y * 100, np.nan), edges, color=C[b], lw=1.4, baseline=None)
mp, mm, edges = d['wp'][1], d['wm'][1], d['wp'][2]
r = mm[:, 1] / mp[:, 1]; re = r * np.hypot(mm[:, 2] / mm[:, 1], mp[:, 2] / mp[:, 1])
ax[2].stairs(r, edges, color=ink2, lw=1.4, baseline=None); ax[2].stairs(r + re, edges, baseline=r - re, fill=True, color=ink2, alpha=0.2, lw=0)
ax[0].set_yscale('log'); ax[0].set_ylim(5e-6, 60); ax[0].set_ylabel(r'd$\sigma$/d$m_T$  [fb/GeV]')
ax[0].legend(frameon=False, loc='upper right', fontsize=9)
ax[0].set_title(rf'LO, {ANALYSIS}, {NS} seeds; band: $\mu_R = \mu_F \in [\mu/2, 2\mu]$', loc='left', color=ink, fontsize=10)
ax[1].set_ylabel('scale / central'); ax[1].text(0.01, 0.92, 'solid 2$\\mu$, dashed $\\mu$/2', transform=ax[1].transAxes, ha='left', va='top', color=ink2, fontsize=8)
ax[2].set_ylim(0, 0.85); ax[2].set_ylabel('W$^-$ / W$^+$')
ax[3].set_ylabel('MC error  [%]'); ax[3].set_yscale('log'); ax[3].set_xlabel(r'$m_T$  [GeV]')
for a in ax:
    a.grid(True, color=grid, lw=0.6); a.set_axisbelow(True); a.spines[['top', 'right']].set_visible(False)
ax[3].set_xlim(200, XMAX)
fig.savefig(f'{OUT}/lo_mT.png', dpi=200, bbox_inches='tight', facecolor='#fcfcfb')
for b in ('wp', 'wm'):
    m = d[b][1]; i = [np.argmin(abs(m[:, 0] - x)) for x in (200, 1000, XMAX - 10)]
    print(b, 'scale/central at 200, 1000, 2990 GeV:', [(round(m[k, 3] / m[k, 1], 3), round(m[k, 5] / m[k, 1], 3)) for k in i])
