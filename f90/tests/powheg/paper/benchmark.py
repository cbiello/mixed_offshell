#!/usr/bin/env python3
"""mixed_ll against POWHEG for the paper setups: LO, NLO QCD, NLO EW at mu = xi mT/2, xi = 1, 1/2, 2.

usage: benchmark.py [run directory]     (default /scratch2/cbiello/mixed_runs/prod)
Reads the finished seeds listed in logs/progress.log, averages the seeds, sums the mixed_ll sectors and
writes benchmark/rates.md, benchmark/<setup>_<boson>.dat (all histograms) and, with matplotlib, the plots.
"""
import sys, re, math, glob, os, collections

R = sys.argv[1] if len(sys.argv) > 1 else '/scratch2/cbiello/mixed_runs/prod'
SETUPS = {'A': ('W_incl', ''), 'B': ('W_ATLASfid', 'B_')}       # mixed_ll analysis, POWHEG prefix
SECTORS = {'lo': ['lo ns na'],
           'nloqcd': ['nloqcd ns v', 'nloqcd ns s', 'nloqcd ns r_is', 'nloqcd gq s', 'nloqcd gq r_is',
                      'nloqcd qg s', 'nloqcd qg r_is'],
           'nloew': ['nloewk ns v', 'nloewk ns s', 'nloewk ns r_is', 'nloewk ns r_fs_53', 'nloewk ns r_fs_54',
                     'nloewk aq s', 'nloewk aq r_is', 'nloewk qa s', 'nloewk qa r_is']}
SCALES = ['1', '0.5', '2']                                      # column order of the mixed_ll histograms


def finished():
    done = set()
    for l in open(f'{R}/logs/progress.log'):
        v = l.split()
        if 'exit' in v and v[v.index('exit') + 1] == '0':
            done.add(tuple(v[:v.index('exit')]))
    return done


def read_nsc(path):
    """{histogram: {left edge: [(value, error) per scale]}}, fb per bin; the last line of a histogram is the overflow"""
    H, name = {}, None
    for l in open(path, errors='replace'):
        m = re.match(r'\s*# -+Obs=(\S+?)-*\s', l)
        if m:
            name = m.group(1); H[name] = {}; continue
        v = l.split()
        if name is None or l.lstrip().startswith('#') or len(v) < 9: continue
        H[name][round(float(v[2]), 6)] = [(float(v[3 + 2 * i]), float(v[4 + 2 * i])) for i in range(3)]
    for h in H.values(): h.pop(max(h))
    return H


def read_pwg(path):
    """{histogram: {left edge: (value, error)}}, fb per bin"""
    H, name = {}, None
    for l in open(path):
        if l.startswith('#'):
            name = l.split()[1]; H[name] = {}; continue
        v = l.replace('D', 'E').split()
        if len(v) == 4:
            w = 1000 * (float(v[1]) - float(v[0]))
            H[name][round(float(v[0]), 6)] = (float(v[2]) * w, float(v[3]) * w)
    return H


def mean(rows):
    """mean over seeds of (value, error); error from the spread with at least 3 seeds"""
    n = len(rows); m = sum(x for x, _ in rows) / n
    if n >= 3: return m, math.sqrt(sum((x - m) ** 2 for x, _ in rows) / (n * (n - 1)))
    return m, math.sqrt(sum(e * e for _, e in rows)) / n


def combine(tables):
    """list of {hist: {edge: (v, e)}} over seeds -> the same with the mean"""
    out = collections.defaultdict(dict)
    for h in tables[0]:
        for x in tables[0][h]:
            out[h][x] = mean([t[h][x] for t in tables if x in t.get(h, {})])
    return out


def nsc(setup, boson, order, done):
    """sum of the sectors of an order (without the LO); returns {scale: {hist: {edge: (v, e)}}} and the seeds used"""
    ana = SETUPS[setup][0]; tot = {s: collections.defaultdict(dict) for s in SCALES}; nseeds = []
    for sec in SECTORS[order]:
        c, ch, s = sec.split()
        files = [f for f in glob.glob(f'{R}/nsc/{ana}/{boson}/{ana}_{boson}_{c}_{ch}_{s}_dynscale*_s*.histo')
                 if ('nsc', ana, boson, c, ch, s, str(int(f[-11:-6]))) in done]
        if not files: return None, 0
        nseeds.append(len(files)); tabs = [read_nsc(f) for f in files]
        for i, xi in enumerate(SCALES):
            m = combine([{h: {x: v[i] for x, v in t[h].items()} for h in t} for t in tabs])
            for h in m:
                for x, (v, e) in m[h].items():
                    a = tot[xi][h].get(x, (0, 0)); tot[xi][h][x] = (a[0] + v, math.hypot(a[1], e))
    return tot, min(nseeds)


def pwg(setup, boson, order, xi, done):
    pre = SETUPS[setup][1]
    files = [f for f in glob.glob(f'{R}/pwg/{boson}/{order}_mu{xi}/s*/pwg-NLO.top')
             if ('pwg', order, boson, xi, f.split('/')[-2][1:]) in done]
    if not files: return None, 0
    m = combine([read_pwg(f) for f in files])
    own = {'rate', 'rate_window', 'mT', 'mT_high', 'mlnu', 'mlnu_high', 'ptl', 'yl', 'ptlnu', 'costhCS'}
    return ({h[len(pre):]: v for h, v in m.items() if h.startswith('B_')} if pre else
            {h: v for h, v in m.items() if h in own}), len(files)


def add(a, b):
    return {h: {x: (a[h][x][0] + b[h][x][0], math.hypot(a[h][x][1], b[h][x][1])) for x in a[h]} for h in a}


def main():
    done = finished(); os.makedirs(f'{R}/benchmark', exist_ok=True)
    rows = ['| setup | boson | order | xi | mixed_ll [fb] | POWHEG [fb] | difference [fb] | ratio - 1 | pull | chi2/bin (worst histogram) | seeds |',
            '|---|---|---|---|---|---|---|---|---|---|---|']
    for setup in SETUPS:
        for boson in ('wp', 'wm'):
            out = open(f'{R}/benchmark/{setup}_{boson}.dat', 'w'); curves = {}
            lo, nlo_seeds = nsc(setup, boson, 'lo', done)
            for order in SECTORS:
                corr, ns = (lo, nlo_seeds) if order == 'lo' else nsc(setup, boson, order, done)
                if lo is None or corr is None: continue
                for xi in SCALES:
                    a = corr[xi] if order == 'lo' else add(lo[xi], corr[xi])
                    b, nb = pwg(setup, boson, order, xi, done)
                    if b is None: continue
                    worst = ('', 0)
                    for h in a:
                        if h not in b: continue
                        out.write(f'# {h} {order} xi={xi}: left edge, mixed_ll, error, POWHEG, error [fb per bin]\n')
                        chi = n = 0
                        for x in sorted(a[h]):
                            if x not in b[h]: continue
                            (v1, e1), (v2, e2) = a[h][x], b[h][x]
                            out.write(f'{x:12.5g} {v1:15.7g} {e1:12.4g} {v2:15.7g} {e2:12.4g}\n')
                            if e1 + e2 > 0 and v1 != 0: chi += (v1 - v2) ** 2 / (e1 ** 2 + e2 ** 2); n += 1
                        out.write('\n\n')
                        if n > 1 and h != 'rate' and chi / n > worst[1]: worst = (h, chi / n)
                        curves[(h, order, xi)] = (a[h], b[h])
                    (v1, e1), (v2, e2) = a['rate'][0.0], b['rate'][0.0]; e = math.hypot(e1, e2)
                    rows.append(f'| {setup} | {boson} | {order} | {xi} | {v1:.3f} +- {e1:.3f} | {v2:.3f} +- {e2:.3f} | '
                                f'{v1 - v2:+.3f} +- {e:.3f} | {(v1 / v2 - 1) * 100:+.3f}% | {(v1 - v2) / e:+.1f} | '
                                f'{worst[1]:.2f} ({worst[0]}) | {ns} / {nb} |')
            out.close(); plots(setup, boson, curves)
    open(f'{R}/benchmark/rates.md', 'w').write('\n'.join(rows) + '\n'); print('\n'.join(rows))


def plots(setup, boson, curves):
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    except ImportError:
        return
    for h in sorted({k[0] for k in curves} - {'rate'}):
        fig, ax = plt.subplots(2, 3, figsize=(15, 7), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
        for j, order in enumerate(SECTORS):
            for xi, col in zip(SCALES, ('k', 'tab:blue', 'tab:red')):
                if (h, order, xi) not in curves: continue
                a, b = curves[(h, order, xi)]; x = sorted(set(a) & set(b))
                ax[0, j].step(x, [a[i][0] for i in x], where='post', color=col, label=f'mixed_ll xi={xi}')
                ax[0, j].step(x, [b[i][0] for i in x], where='post', color=col, ls=':', label=f'POWHEG xi={xi}')
                r = [(a[i][0] / b[i][0] - 1, math.hypot(a[i][1], b[i][1]) / abs(b[i][0])) if b[i][0] else (0, 0) for i in x]
                ax[1, j].errorbar(x, [v for v, _ in r], [e for _, e in r], fmt='.', color=col, ms=3, lw=0.6)
            ax[0, j].set_title(f'{order}  setup {setup}  {boson}  {h}'); ax[0, j].set_yscale('symlog', linthresh=1e-3)
            ax[1, j].axhline(0, color='gray', lw=0.5); ax[1, j].set_ylim(-0.05, 0.05); ax[1, j].set_xlabel(h)
        ax[0, 0].set_ylabel('fb per bin'); ax[1, 0].set_ylabel('mixed_ll / POWHEG - 1'); ax[0, 0].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(f'{R}/benchmark/{setup}_{boson}_{h}.png', dpi=90); plt.close(fig)


main()
