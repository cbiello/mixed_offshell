#!/usr/bin/env python3
"""Combine the NLO seed histograms of the cluster runs (cluster_KIT/submit.sh), from f90/:
    python3 cluster_KIT/combine_nlo.py [nloqcd|nloewk ...] [--seeds 10] [--its 5] [--partial]
Per sector: nsc-hsummer ad3 of the seeds 1..n (s00000 = warm-up, left out), divided by n. Per correction: nsc-hsummer ad3
of the sectors. Writes to histo/comb/: <b>_<corr>_<ch>_<sec>_sum<n>.histo, <b>_<corr>_total.histo, <corr>_rates.md (total
and per-sector rates, error from the histograms and from the seed spread), <corr>_<b>_lo_nlo.png (LO, LO + dNLO, K factor),
<corr>_<b>_sectors.png (|d sigma/d m_T| per sector, sector / total). LO from histo/comb/W_ATLAS_<b>_lo_sum10.histo.
Without --partial a sector is used only if all its seeds finished --its iterations."""
import glob, os, re, subprocess, sys, statistics as st
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = lambda k, d: int(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
NS, ITS, PARTIAL = opt('--seeds', 10), opt('--its', 5), '--partial' in sys.argv
SECS = {'nloqcd': ['ns v', 'ns s', 'ns r_is', 'gq s', 'gq r_is', 'qg s', 'qg r_is'],
        'nloewk': ['ns v', 'ns s', 'ns r_is', 'ns r_fs_53', 'ns r_fs_54', 'aq s', 'aq r_is', 'qa s', 'qa r_is']}
TITLE = {'nloqcd': 'NLO QCD', 'nloewk': 'NLO EW'}; BOS = {'wp': 'W$^+$', 'wm': 'W$^-$'}
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948', '#52514e']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'histo')); os.makedirs('comb', exist_ok=True)
seedf = lambda b, c, ch, s, k: f'W_ATLAS_{b}_{c}_{ch}_{s}_dynscale_xMuR_1.00_xMuF_1.00_s{k:05d}.histo'

def read(f, scale=1.0):                         #-- {histo: array [x, v1, e1, v2, e2, v3, e3]}
    out = {}
    for l in open(f):
        x = l.split()
        if len(x) == 9 and not l.lstrip().startswith('#'): out.setdefault(int(x[0]), []).append([float(v) for v in x[2:]])
    return {h: np.array(v) * np.array([1] + [scale] * 6) for h, v in out.items()}

def iters(f):
    m = re.search(r'final run,\s+iteration\s+(\d+)', open(f, errors='replace').read()) if os.path.exists(f) else None
    return int(m.group(1)) if m else 0

def hsum(files, out):
    subprocess.run(['../bin/nsc-hsummer', 'ad3', *files, out], check=True, stdout=subprocess.DEVNULL)

def stairs(ax, x, y, **kw):                     #-- |y| on a log axis, dashed where y < 0
    e = np.append(x, x[-1] + (x[1] - x[0]))
    ax.stairs(np.where(y > 0, y, np.nan), e, baseline=None, **kw)
    ax.stairs(np.where(y < 0, -y, np.nan), e, baseline=None, ls='--', **{k: v for k, v in kw.items() if k != 'label'})

def style(ax):
    for a in ax:
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True); a.spines[['top', 'right']].set_visible(False)
    ax[-1].set_xlim(200, 3000); ax[-1].set_xlabel(r'$m_T$  [GeV]')

for corr in args or ['nloqcd', 'nloewk']:
    md = [f'# {TITLE[corr]} correction (W_ATLAS, {NS} seeds per sector, `nsc-hsummer ad3` / {NS})\n']
    for b in ('wp', 'wm'):
        secs, rows = [], []
        for cs in SECS[corr]:
            ch, s = cs.split(); fs = [seedf(b, corr, ch, s, k) for k in range(1, NS + 1)]
            done = [iters(f) for f in fs]
            if min(done) < ITS and not (PARTIAL and min(done) > 0):
                print(f'{b} {corr} {cs}: iterations per seed {done}, not complete: skipped'); continue
            out = f'comb/{b}_{corr}_{ch}_{s}_sum{NS}.histo'; hsum(fs, out)
            h = read(out, 1 / NS); r = [read(f)[1][0, 1] for f in fs]
            secs.append((cs, out, h)); rows.append((cs, h[1][0, 1], h[1][0, 2], st.stdev(r) / NS**0.5, min(done)))
        if not secs: continue
        hsum([o for _, o, _ in secs], f'comb/{b}_{corr}_total.histo'); tot = read(f'comb/{b}_{corr}_total.histo', 1 / NS)
        T, E = tot[1][0, 1], tot[1][0, 2]; Es = sum(r[3]**2 for r in rows)**0.5
        lo = read(f'comb/W_ATLAS_{b}_lo_sum10.histo', 1 / 10) if os.path.exists(f'comb/W_ATLAS_{b}_lo_sum10.histo') else None
        md += [f'\n## {BOS[b].replace("$", "").replace("^", "")}: total {TITLE[corr]} correction = {T:.4f} ± {E:.2g} fb '
               f'(seed spread {Es:.2g})' + (f', {100 * T / lo[1][0, 1]:+.2f}% of LO ({lo[1][0, 1]:.1f} fb)' if lo else '') + '\n',
               '| Sector | δσ [fb] | error (histo) | error (seed spread) | share of total | iterations |', '|---|---|---|---|---|---|']
        md += [f'| `{cs}` | {v:.4f} | {e:.2g} | {es:.2g} | {100 * v / T:+.1f}% | {it} |' for cs, v, e, es, it in rows]
        print('\n'.join(md[-len(rows) - 3:]))
        m = lambda h: h[2][(h[2][:, 0] >= 200) & (h[2][:, 0] < 3000)]
        if lo:                                  #-- plot 3: LO, LO + dNLO and K factor
            L, D = m(lo), m(tot); x = L[:, 0]; e = np.append(x, x[-1] + 10)
            fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={'height_ratios': [3, 1.5], 'hspace': 0.08})
            ax[0].stairs(L[:, 1] / 10, e, baseline=None, color=INK2, lw=1.4, label=f'LO   {lo[1][0, 1]:.1f} fb')
            ax[0].stairs((L[:, 1] + D[:, 1]) / 10, e, baseline=None, color=PAL[0], lw=1.6, label=f'LO + δ{TITLE[corr]}   {lo[1][0, 1] + T:.1f} fb')
            K, Ke = 1 + D[:, 1] / L[:, 1], D[:, 2] / L[:, 1]
            ax[1].stairs(K, e, baseline=None, color=PAL[0], lw=1.4); ax[1].stairs(K + Ke, e, baseline=K - Ke, fill=True, color=PAL[0], alpha=0.25, lw=0)
            ax[1].axhline(1, color=INK2, lw=0.8)
            ax[0].set_yscale('log'); ax[0].set_ylabel(r'd$\sigma$/d$m_T$  [fb/GeV]'); ax[0].legend(frameon=False)
            ax[0].set_title(f'{BOS[b]}, {TITLE[corr]}: central scale, band = MC error', loc='left', color=INK, fontsize=10)
            ax[1].set_ylabel(f'(LO + δ{TITLE[corr]}) / LO'); style(ax)
            fig.savefig(f'comb/{corr}_{b}_lo_nlo.png', dpi=200, bbox_inches='tight', facecolor='#fcfcfb'); plt.close(fig)
        D = m(tot)                              #-- plot 4: sectors and sector / total
        fig, ax = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True, gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.08})
        x = D[:, 0]; e = np.append(x, x[-1] + 10)
        stairs(ax[0], x, D[:, 1] / 10, color=INK, lw=1.8, label=f'total   {T:.2f} fb')
        for i, (cs, _, h) in enumerate(secs):
            S = m(h); y, ey = S[:, 1] / 10, S[:, 2] / 10
            stairs(ax[0], x, y, color=PAL[i % len(PAL)], lw=1.2, label=f'{cs}   {h[1][0, 1]:.2f} fb')
            ax[0].stairs(np.abs(y) + ey, e, baseline=np.clip(np.abs(y) - ey, 1e-12, None), fill=True, color=PAL[i % len(PAL)], alpha=0.2, lw=0)
            r = np.where(D[:, 1] != 0, S[:, 1] / D[:, 1], np.nan); er = np.where(D[:, 1] != 0, S[:, 2] / np.abs(D[:, 1]), np.nan)
            ax[1].stairs(r, e, baseline=None, color=PAL[i % len(PAL)], lw=1.2); ax[1].stairs(r + er, e, baseline=r - er, fill=True, color=PAL[i % len(PAL)], alpha=0.2, lw=0)
        ax[1].axhline(0, color=INK2, lw=0.8); ax[1].axhline(1, color=INK2, lw=0.8, ls=':')
        ax[0].set_yscale('log'); ax[0].set_ylabel(r'|d$\sigma$/d$m_T$|  [fb/GeV]'); ax[0].legend(frameon=False, fontsize=8, ncol=2)
        ax[0].set_title(f'{BOS[b]}, {TITLE[corr]} per sector (dashed: negative), band = MC error', loc='left', color=INK, fontsize=10)
        ax[1].set_ylabel(f'sector / total δ{TITLE[corr]}'); style(ax)
        fig.savefig(f'comb/{corr}_{b}_sectors.png', dpi=200, bbox_inches='tight', facecolor='#fcfcfb'); plt.close(fig)
    open(f'comb/{corr}_rates.md', 'w').write('\n'.join(md) + '\n')
