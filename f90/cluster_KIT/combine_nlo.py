#!/usr/bin/env python3
"""Combine the NLO seed histograms of the cluster runs (cluster_KIT/submit.sh), from f90/:
    python3 cluster_KIT/combine_nlo.py [nloqcd|nloewk|nnlo ...] [--analysis W_ATLAS|W_incl|W_ATLASfid] [--seeds 10] [--its 5] [--partial]
Per sector: nsc-hsummer ad3 of the seeds 1..n (s00000 = warm-up, left out), divided by n. Per correction: nsc-hsummer ad3
of the sectors. Writes to histo/comb/ (W_ATLAS) or histo/comb_<analysis>/: <b>_<corr>_<ch>_<sec>_sum<n>.histo,
<b>_<corr>_total.histo, <corr>_rates.md (total and per-sector rates, error from the histograms and from the seed spread),
<corr>_<b>_lo_nlo.png (LO, LO + dNLO, K factor), <corr>_<b>_sectors.png (|d sigma/d m_T| per sector, sector / total).
LO from <analysis>_<b>_lo_sum10.histo of the same directory (cluster_KIT/plot_lo.py).
Without --partial a sector is used only if all its seeds finished --its iterations."""
import os, re, subprocess, sys, statistics as st
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

args = [a for i, a in enumerate(sys.argv[1:]) if not a.startswith('--') and sys.argv[i] not in ('--seeds', '--its', '--analysis')]
opt = lambda k, d: type(d)(sys.argv[sys.argv.index(k) + 1]) if k in sys.argv else d
NS, ITS, PARTIAL, ANALYSIS = opt('--seeds', 10), opt('--its', 5), '--partial' in sys.argv, opt('--analysis', 'W_ATLAS')
MT = {'W_ATLAS': ([2], 3000), 'W_incl': ([3, 4], 3000), 'W_ATLASfid': ([2], 5000)}[ANALYSIS]   #-- m_T histograms, upper end
OUT = 'comb' if ANALYSIS == 'W_ATLAS' else f'comb_{ANALYSIS}'
SECS = {'nloqcd': ['ns v', 'ns s', 'ns r_is', 'gq s', 'gq r_is', 'qg s', 'qg r_is'],
        'nloewk': ['ns v', 'ns s', 'ns r_is', 'ns r_fs_53', 'ns r_fs_54', 'aq s', 'aq r_is', 'qa s', 'qa r_is'],
        'nnlo': [' '.join(l.split()[:2]) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nnlo_qqbar.txt'))
                 if l.strip() and not l.startswith('#')]}   #-- the sectors run so far: q qbar' channel
GROUPS = {'nnlo': {'s_vewk (EW virtual x int. QCD)': ['ns s_vewk'], 's_vqcd (QCD virtual x int. QED)': ['ns s_vqcd'],
                   's_12 (double integrated)': ['ns s_12'], 's_oqcd (gluon resolved)': ['ns s_oqcd'],
                   's_oewk (photon resolved)': ['ns s_oewk_is', 'ns s_oewk_fs_53', 'ns s_oewk_fs_54'],
                   'rvqcd (QCD loop + photon)': ['ns rvqcd_is', 'ns rvqcd_fs_53', 'ns rvqcd_fs_54'],
                   'rvewk (EW loop + gluon)': ['ns rvewk_is'], 'RR ns_ga (gluon + photon)': [f'ns_ga {x}' for x in
                   ['rr_5161a', 'rr_5161c', 'rr_5162', 'rr_5163', 'rr_5164', 'rr_5261', 'rr_5262a', 'rr_5262c', 'rr_5263', 'rr_5264']]}}
TITLE = {'nloqcd': 'NLO QCD', 'nloewk': 'NLO EW', 'nnlo': "NNLO q qbar'"}; BOS = {'wp': 'W$^+$', 'wm': 'W$^-$'}
PAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948', '#52514e']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'histo')); os.makedirs(OUT, exist_ok=True)
seedf = lambda b, c, ch, s, k: f'{ANALYSIS}_{b}_{c}_{ch}_{s}_dynscale_xMuR_1.00_xMuF_1.00_s{k:05d}.histo'

def read(f, scale=1.0):                         #-- {histo: array [x, v1, e1, v2, e2, v3, e3]}, last row = overflow
    out = {}
    for l in open(f):
        x = l.split()
        if len(x) == 9 and not l.lstrip().startswith('#'): out.setdefault(int(x[0]), []).append([float(v) for v in x[2:]])
    return {h: np.array(v) * np.array([1] + [scale] * 6) for h, v in out.items()}

def mt(h):                                      #-- m_T bins of all m_T histograms: edges, d sigma/d m_T and its error (central)
    rows = np.vstack([np.hstack([h[k][:-1], np.full((len(h[k]) - 1, 1), h[k][1, 0] - h[k][0, 0])]) for k in MT[0]])
    rows = rows[(rows[:, 0] >= 200) & (rows[:, 0] < MT[1])]
    return np.append(rows[:, 0], rows[-1, 0] + rows[-1, 7]), rows[:, 1] / rows[:, 7], rows[:, 2] / rows[:, 7]

def iters(f):
    m = re.search(r'final run,\s+iteration\s+(\d+)', open(f, errors='replace').read()) if os.path.exists(f) else None
    return int(m.group(1)) if m else 0

def hsum(files, out):
    subprocess.run(['../bin/nsc-hsummer', 'ad3', *files, out], check=True, stdout=subprocess.DEVNULL)

def stairs(ax, e, y, **kw):                     #-- |y| on a log axis, dashed where y < 0
    ax.stairs(np.where(y > 0, y, np.nan), e, baseline=None, **kw)
    ax.stairs(np.where(y < 0, -y, np.nan), e, baseline=None, ls='--', **{k: v for k, v in kw.items() if k != 'label'})

def style(ax):
    for a in ax:
        a.grid(True, color=GRID, lw=0.6); a.set_axisbelow(True); a.spines[['top', 'right']].set_visible(False)
    ax[-1].set_xlim(200, MT[1]); ax[-1].set_xlabel(r'$m_T$  [GeV]')

for corr in args or ['nloqcd', 'nloewk']:
    md = [f'# {TITLE[corr]} correction ({ANALYSIS}, {NS} seeds per sector, `nsc-hsummer ad3` / {NS})\n']
    for b in ('wp', 'wm'):
        secs, rows = [], []
        for cs in SECS[corr]:
            ch, s = cs.split(); fs = [seedf(b, corr, ch, s, k) for k in range(1, NS + 1)]
            done = [iters(f) for f in fs]
            if min(done) < ITS and not (PARTIAL and min(done) > 0):
                print(f'{b} {corr} {cs}: iterations per seed {done}, not complete: skipped'); continue
            out = f'{OUT}/{b}_{corr}_{ch}_{s}_sum{NS}.histo'; hsum(fs, out)
            h = read(out, 1 / NS); r = [read(f)[1][0, 1] for f in fs]
            secs.append((cs, out, h)); rows.append((cs, h[1][0, 1], h[1][0, 2], st.stdev(r) / NS**0.5, min(done)))
        if not secs: continue
        hsum([o for _, o, _ in secs], f'{OUT}/{b}_{corr}_total.histo'); tot = read(f'{OUT}/{b}_{corr}_total.histo', 1 / NS)
        T, E = tot[1][0, 1], tot[1][0, 2]; Es = sum(r[3]**2 for r in rows)**0.5
        lof = f'{OUT}/{ANALYSIS}_{b}_lo_sum10.histo'; lo = read(lof, 1 / 10) if os.path.exists(lof) else None
        md += [f'\n## {BOS[b].replace("$", "").replace("^", "")}: total {TITLE[corr]} correction = {T:.4f} ± {E:.2g} fb '
               f'(seed spread {Es:.2g})' + (f', {100 * T / lo[1][0, 1]:+.2f}% of LO ({lo[1][0, 1]:.1f} fb)' if lo else '') + '\n',
               '| Sector | δσ [fb] | error (histo) | error (seed spread) | share of total | iterations |', '|---|---|---|---|---|---|']
        md += [f'| `{cs}` | {v:.4f} | {e:.2g} | {es:.2g} | {100 * v / T:+.1f}% | {it} |' for cs, v, e, es, it in rows]
        print('\n'.join(md[-len(rows) - 3:]))
        e, D, De = mt(tot)
        if lo:                                  #-- plot 3: LO, LO + dNLO and K factor
            _, L, _ = mt(lo)
            fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={'height_ratios': [3, 1.5], 'hspace': 0.08})
            ax[0].stairs(L, e, baseline=None, color=INK2, lw=1.4, label=f'LO   {lo[1][0, 1]:.1f} fb')
            ax[0].stairs(L + D, e, baseline=None, color=PAL[0], lw=1.6, label=f'LO + δ{TITLE[corr]}   {lo[1][0, 1] + T:.1f} fb')
            K, Ke = 1 + D / L, De / L
            ax[1].stairs(K, e, baseline=None, color=PAL[0], lw=1.4); ax[1].stairs(K + Ke, e, baseline=K - Ke, fill=True, color=PAL[0], alpha=0.25, lw=0)
            ax[1].axhline(1, color=INK2, lw=0.8)
            ax[0].set_yscale('log'); ax[0].set_ylabel(r'd$\sigma$/d$m_T$  [fb/GeV]'); ax[0].legend(frameon=False)
            ax[0].set_title(f'{BOS[b]}, {TITLE[corr]} ({ANALYSIS}): central scale, band = MC error', loc='left', color=INK, fontsize=10)
            ax[1].set_ylabel(f'(LO + δ{TITLE[corr]}) / LO'); style(ax)
            fig.savefig(f'{OUT}/{corr}_{b}_lo_nlo.png', dpi=200, bbox_inches='tight', facecolor='#fcfcfb'); plt.close(fig)
        fig, ax = plt.subplots(2, 1, figsize=(7.5, 7.5), sharex=True, gridspec_kw={'height_ratios': [3, 2], 'hspace': 0.08})
        stairs(ax[0], e, D, color=INK, lw=1.8, label=f'total   {T:.2f} fb')     #-- plot 4: sectors and sector / total
        if corr in GROUPS:                      #-- many sectors: one curve per contribution (values add, errors in quadrature)
            hs = {cs: h for cs, _, h in secs}; grouped = []
            for g, members in GROUPS[corr].items():
                m = [hs[x] for x in members if x in hs]
                if not m: continue
                gh = {k: np.array(m[0][k]) for k in m[0]}
                for k in gh:
                    gh[k][:, 1::2] = sum(x[k][:, 1::2] for x in m); gh[k][:, 2::2] = np.sqrt(sum(x[k][:, 2::2]**2 for x in m))
                grouped.append((g, None, gh))
            secs = grouped
        for i, (cs, _, h) in enumerate(secs):
            _, y, ey = mt(h); c = PAL[i % len(PAL)]
            stairs(ax[0], e, y, color=c, lw=1.2, label=f'{cs}   {h[1][0, 1]:.2f} fb')
            ax[0].stairs(np.abs(y) + ey, e, baseline=np.clip(np.abs(y) - ey, 1e-12, None), fill=True, color=c, alpha=0.2, lw=0)
            r = np.where(D != 0, y / D, np.nan); er = np.where(D != 0, ey / np.abs(D), np.nan)
            ax[1].stairs(r, e, baseline=None, color=c, lw=1.2); ax[1].stairs(r + er, e, baseline=r - er, fill=True, color=c, alpha=0.2, lw=0)
        ax[1].axhline(0, color=INK2, lw=0.8); ax[1].axhline(1, color=INK2, lw=0.8, ls=':')
        ax[0].set_yscale('log'); ax[0].set_ylabel(r'|d$\sigma$/d$m_T$|  [fb/GeV]'); ax[0].legend(frameon=False, fontsize=8, ncol=2)
        ax[0].set_title(f'{BOS[b]}, {TITLE[corr]} ({ANALYSIS}) per {"contribution" if corr in GROUPS else "sector"} (dashed: negative), band = MC error', loc='left', color=INK, fontsize=10)
        ax[1].set_ylabel(f'{"contribution" if corr in GROUPS else "sector"} / total δ{TITLE[corr]}'); style(ax)
        fig.savefig(f'{OUT}/{corr}_{b}_sectors.png', dpi=200, bbox_inches='tight', facecolor='#fcfcfb'); plt.close(fig)
    open(f'{OUT}/{corr}_rates.md', 'w').write('\n'.join(md) + '\n')
