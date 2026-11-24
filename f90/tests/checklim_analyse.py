#!/usr/bin/env python3
"""Intensive limit checks (PLAN C1.5). Reads the -checklim_csv files <boson>.<corr>.<ch>.<sec>.csv of a directory (tests/checklim_run.sh):
    tests/checklim_analyse.py <csvdir> <outdir>
Per point and pair (A, B) that must agree in a limit: r(d) = |A/B - 1| at depths 1e-1..1e-10; slope p_fit of log r
against log d between 1e-2 and 1e-6 (median of pairwise slopes); pass if p_fit >= p - 0.2 (p = 1 soft limit, 1/2
otherwise), or the same between 1e-4 and 1e-8, or r(1e-6) < 1e-9 (a point cut at more than two of the depths 1e-2..1e-6
is not graded); coefficient c = r(1e-6) / 1e-6^p. Limits without own counterterm (U lines): pass if x^k |H + counterterms|
somewhere below 1e-3 is 3 times smaller than at 1e-2..1e-4; not graded if H is cut at 1e-10 (point at a cut edge).
Writes <outdir>/<boson>/summary.csv, worst.csv, checklim_<boson>.pdf (heatmap, bands, histograms) and checklim.json."""
import csv, glob, json, math, os, statistics as st, sys
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

csvdir, outdir = sys.argv[1:3]
D = np.array([10.0**-k for k in range(1, 11)])
ONLO = {'H': (), 'S': ('S',), 'C1': ('C1',), 'C2': ('C2',), 'SC1': ('S', 'C1'), 'SC2': ('S', 'C2')}
SOFT = {'S', 'S5', 'S6'}

def namer(sec):                                #-- physical names: X = TC or C5, Y = C6, C5 or C56; order soft, inner collinear, TC
    tc = sec[-1] in 'abcd'
    sub = {'X': 'TC' if tc else 'C5', 'Y': {'a': 'C6', 'c': 'C5', 'b': 'C56', 'd': 'C56'}.get(sec[-1], 'C6')}
    order = ['S5', 'S6', 'Y', 'X'] if tc else ['S5', 'S6', 'X', 'Y']
    def phys(lab):
        name, *g = lab.split('/')
        if name in ONLO: return lab
        parts = sorted(name.split('+'), key=lambda x: order.index(x) if x in order else -1)
        return ''.join(sub.get(x, x) for x in parts) + ('/' + g[0] if g else '')
    return phys

def parts(lim): return ONLO[lim] if lim in ONLO else tuple(lim.split('+'))

def slope(r, lo=2, hi=6):                      #-- Theil-Sen slope of log10 r vs log10 d over 1e-lo..1e-hi
    pts = [(-k, math.log10(r[k - 1])) for k in range(lo, hi + 1) if r[k - 1] > 0 and math.isfinite(r[k - 1])]
    s = [(y2 - y1) / (x2 - x1) for i, (x1, y1) in enumerate(pts) for x2, y2 in pts[i + 1:]]
    return st.median(s) if s else float('nan')

def falls(f):                                  #-- None: not graded (cut, or H cut at the deepest depth: point at a cut edge)
    v = [f[i] for i in range(1, 10) if f[i] > 0 and math.isfinite(f[i])]
    if len(v) < 3: return None
    if min(v[2:]) < 0.3 * st.median(v[:3]): return True     #-- min: double precision noise rises again at the deepest depths
    return False if f[9] != 0 else None

res = defaultdict(lambda: defaultdict(dict))   #-- res[boson][sector][(limit, pair)] = list of per-point records
unc = defaultdict(lambda: defaultdict(dict))
for path in sorted(glob.glob(os.path.join(csvdir, '*.csv'))):
    boson, corr, ch, sec = os.path.basename(path)[:-4].split('.')
    sector = f'{corr} {ch} {sec}'; phys = namer(sec)
    for row in csv.reader(open(path)):
        if row[0] == 'R':
            n, lim, a, b = int(row[1]), row[2], row[3], row[4]
            r = [abs(float(x) - 1) for x in row[5:15]]
            if sum(v > 0 and math.isfinite(v) for v in r[1:6]) < 3: continue   #-- cut at this point: not graded
            p = 1.0 if set(parts(lim)) <= SOFT else 0.5
            ps = slope(r)
            ok = ps >= p - 0.2 or slope(r, 4, 8) >= p - 0.2 or 0 <= r[5] < 1e-9   #-- 1e-4..1e-8: late onset, x log x
            c = r[5] / 1e-6**p if math.isfinite(r[5]) else float('nan')
            res[boson][sector].setdefault((phys(lim), f'{phys(a)} = {phys(b)}'), []).append(dict(n=n, r=r, p=p, ps=ps, ok=ok, c=c))
        elif row[0] == 'U':
            n, lim = int(row[1]), row[2]
            f = [float(x) for x in row[4:14]]
            unc[boson][sector].setdefault(phys(lim), []).append(dict(n=n, f=f, ok=falls(f)))

def q(v, x):
    v = [y for y in v if y > 0 and math.isfinite(y)]
    return float(np.quantile(v, x)) if v else None

js = {}
for boson in sorted(set(res) | set(unc)):
    od = os.path.join(outdir, boson); os.makedirs(od, exist_ok=True)
    out = {}
    with open(os.path.join(od, 'summary.csv'), 'w', newline='') as fs, open(os.path.join(od, 'worst.csv'), 'w', newline='') as fw:
        ws, ww = csv.writer(fs), csv.writer(fw)
        ws.writerow(['sector', 'limit', 'pair', 'points', 'pass', 'expected_slope', 'median_slope', 'median_log10_c'])
        ww.writerow(['sector', 'limit', 'pair', 'point', 'slope', 'r_1e-2', 'r_1e-4', 'r_1e-6'])
        for sector in sorted(set(res[boson]) | set(unc[boson])):
            S = out.setdefault(sector, {'pairs': [], 'unc': []})
            for (lim, pr), recs in res[boson][sector].items():
                npass = sum(x['ok'] for x in recs)
                ms = st.median([x['ps'] for x in recs if math.isfinite(x['ps'])] or [float('nan')])
                lc = [math.log10(x['c']) for x in recs if x['c'] > 0 and math.isfinite(x['c'])]
                ws.writerow([sector, lim, pr, len(recs), f'{npass/len(recs):.3f}', recs[0]['p'], f'{ms:.2f}', f'{st.median(lc):.2f}' if lc else ''])
                for x in sorted(recs, key=lambda x: x['ps'] if math.isfinite(x['ps']) else -9)[:3]:
                    if not x['ok']: ww.writerow([sector, lim, pr, x['n'], f"{x['ps']:.2f}", *(f'{x["r"][i]:.2e}' for i in (1, 3, 5))])
                S['pairs'].append({'limit': lim, 'pair': pr, 'n': len(recs), 'pass': npass / len(recs), 'p': recs[0]['p'], 'slope': ms,
                    'bands': [[q([x['r'][i] for x in recs], z) for z in (0.025, 0.16, 0.5, 0.84, 0.975)] for i in range(10)],
                    'slopes': [round(x['ps'], 3) for x in recs if math.isfinite(x['ps'])],
                    'logc': [round(v, 2) for v in lc], 'worst': [x['n'] for x in recs if not x['ok']][:10]})
            for lim, recs in unc[boson][sector].items():
                v = [x for x in recs if x['ok'] is not None]
                if not v: continue
                ws.writerow([sector, lim, 'H + counterterms integrable', len(v), f"{sum(x['ok'] for x in v)/len(v):.3f}", '', '', ''])
                for x in v:
                    if not x['ok']: ww.writerow([sector, lim, 'H + counterterms integrable', x['n'], '', *(f'{x["f"][i]:.2e}' for i in (1, 3, 5))])
                S['unc'].append({'limit': lim, 'n': len(v), 'pass': sum(x['ok'] for x in v) / len(v),
                    'bands': [[q([x['f'][i] for x in v], z) for z in (0.16, 0.5, 0.84)] for i in range(10)]})
    js[boson] = out

    #-- PDF: heatmap (worst pair per limit), then per sector the bands and slope histograms
    sectors = sorted(out); lims = sorted({p['limit'] for s in out.values() for p in s['pairs'] + s['unc']}, key=lambda l: (len(l), l))
    M = np.full((len(sectors), len(lims)), np.nan)
    for i, s in enumerate(sectors):
        for p in out[s]['pairs'] + out[s]['unc']:
            j = lims.index(p['limit']); M[i, j] = p['pass'] if np.isnan(M[i, j]) else min(M[i, j], p['pass'])
    with PdfPages(os.path.join(od, f'checklim_{boson}.pdf')) as pdf:
        fig, ax = plt.subplots(figsize=(max(8, 0.32 * len(lims) + 4), max(6, 0.16 * len(sectors) + 2)))
        im = ax.imshow(M, cmap='RdYlGn', vmin=0.5, vmax=1, aspect='auto', interpolation='nearest')
        ax.set_xticks(range(len(lims)), lims, rotation=90, fontsize=6); ax.set_yticks(range(len(sectors)), sectors, fontsize=5)
        ax.set_title(f'{boson}: fraction of points passing (worst pair per limit; limits without counterterm: integrable)', fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.02); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        for s in sectors:
            P = out[s]['pairs']
            if not P: continue
            bylim = defaultdict(list)
            for p in P: bylim[p['limit']].append(p)
            nl = len(bylim); fig, axs = plt.subplots(2, nl, figsize=(3.2 * nl, 6), squeeze=False)
            for j, (lim, ps) in enumerate(bylim.items()):
                a = axs[0][j]
                for p in ps:
                    b = np.array([[np.nan if v is None else v for v in row] for row in p['bands']])
                    l, = a.loglog(D, b[:, 2], lw=1, label=p['pair']); a.fill_between(D, b[:, 1], b[:, 3], alpha=0.25, color=l.get_color())
                p0 = ps[0]['p']; a.loglog(D, D**p0 * 10, 'k:', lw=0.8, label=f'd^{p0:g}')
                a.invert_xaxis(); a.set_title(f'{lim}', fontsize=8); a.legend(fontsize=5); a.set_xlabel('depth', fontsize=7)
                h = axs[1][j]
                for p in ps: h.hist(p['slopes'], bins=np.linspace(0, 1.5, 31), histtype='step', label=p['pair'])
                h.axvline(p0 - 0.2, color='r', ls='--', lw=0.8); h.set_xlabel('slope', fontsize=7)
            fig.suptitle(f'{boson} {s}: |A/B - 1| (median, 68% band) and per-point slopes', fontsize=9)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
json.dump(js, open(os.path.join(outdir, 'checklim.json'), 'w'), separators=(',', ':'))

#-- flat rows for the dashboard: pairs.json (one row per limit and pair), bands.json (median and 68% band per depth),
#-- hist.json (per-point slopes, bins of 0.05 from 0 to 1.5, values outside clamped to the end bins)
dd = os.path.join(outdir, 'dash'); os.makedirs(dd, exist_ok=True)
pairs, bands, hist = [], [], []
for boson, out in js.items():
    for sector, v in out.items():
        for p in v['pairs']:
            pairs.append(dict(boson=boson, sector=sector, limit=p['limit'], pair=p['pair'], kind='pair', n=p['n'], passed=round(p['pass'], 4),
                              expected=p['p'], slope=None if not math.isfinite(p['slope']) else round(p['slope'], 3),
                              logc=round(st.median(p['logc']), 2) if p['logc'] else None))
            for i, b in enumerate(p['bands']):
                if b[2] is not None: bands.append(dict(boson=boson, sector=sector, limit=p['limit'], pair=p['pair'], depth=10.0**-(i + 1),
                                                       lo=b[1], mid=b[2], hi=b[3]))
            cnt = defaultdict(int)
            for s in p['slopes']: cnt[min(29, max(0, int(s / 0.05)))] += 1
            for k, c in sorted(cnt.items()): hist.append(dict(boson=boson, sector=sector, limit=p['limit'], pair=p['pair'], bin=round(k * 0.05, 2), count=c))
        for u in v['unc']:
            pairs.append(dict(boson=boson, sector=sector, limit=u['limit'], pair='H + counterterms integrable', kind='integrable', n=u['n'],
                              passed=round(u['pass'], 4), expected=None, slope=None, logc=None))
for name, rows in (('pairs', pairs), ('bands', bands), ('hist', hist)):
    json.dump(rows, open(os.path.join(dd, f'{name}.json'), 'w'), separators=(',', ':'))
for boson, out in js.items():
    rows = [(s, p) for s, v in out.items() for p in v['pairs'] + v['unc']]
    bad = [(s, p['limit'], p.get('pair', 'integrable'), p['pass']) for s, p in rows if p['pass'] < 0.95]
    print(f"{boson}: {len(out)} sectors, {len(rows)} limit/pair rows, {len(bad)} with less than 95% of the points passing")
    for b in bad[:30]: print('   ', *b)
