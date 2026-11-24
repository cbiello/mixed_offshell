#!/usr/bin/env python3
"""Limit table of a sector (PLAN C1.5), from `-checklim` (check_limits_rr / check_limits_onlo), from f90/:
    tests/checklim_tables.py -boson wp -corr nnlo -ch gq -sec rr_5161c [--npts 20]
For each limit and each pair of terms that must agree there: median |a - b| / |b| over the points at depths 1e-2..1e-8,
slope of log(mismatch) against log(depth) between 1e-2 and 1e-6, grade ok if slope >= 0.3. Limits without own
counterterm: median x^k |H + counterterms| and the number of points graded integrable."""
import argparse, math, os, re, statistics as st, subprocess, sys, threading

ap = argparse.ArgumentParser()
ap.add_argument('--npts', type=int, default=20)
args, mixed = ap.parse_known_args()
DEPTHS = [1e-2, 1e-4, 1e-6, 1e-8]
ONLO = {'H': (), 'S': ('S',), 'C1': ('C1',), 'C2': ('C2',), 'SC1': ('S', 'C1'), 'SC2': ('S', 'C2')}

def key(lab):                                  #-- set of limits of a column: 'S5+X' (RR) or 'SC1/2' (ONLO, group 2)
    name, grp = (lab.split('/') + [''])[:2]
    s = ONLO[name] if grp else (() if name == 'H' else tuple(name.split('+')))
    return frozenset(s), grp

def run():
    cmd = ['./mixed_ll', *mixed, '-checklim']
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace')
    def keys():
        try:
            for _ in range(10**5): os.write(p.stdin.fileno(), b'\n' * 100)
        except OSError: pass
    threading.Thread(target=keys, daemon=True).start()
    pts, cur, tab, head, unc = [], None, None, None, False
    for line in p.stdout:
        s = line.strip()
        if s.startswith('point:'): cur = {'tabs': {}, 'unc': {}}; pts.append(cur); tab = None; unc = False; continue
        if cur is None: continue
        if s.startswith('limits without'): unc = True; continue
        if s.startswith('Press a key'):
            if len(pts) >= args.npts: break
            continue
        if unc:
            v = s.split()
            if len(v) == 6 and v[-1] in ('integrable', 'SINGULAR', 'cut'): cur['unc'][v[0]] = ([float(x) for x in v[1:5]], v[-1])
            continue
        if s.startswith('limit '): tab = s[6:]; cur['tabs'][tab] = []; continue
        if s.startswith('x('): head = s.split(); continue
        if tab and re.match(r'^[0-9.E+-]+ ', s):
            v = s.split(); nx = sum(h.startswith('x(') for h in head)
            cur['tabs'][tab].append(dict(zip(head[nx:], map(float, v[nx:]))))
    p.kill(); p.wait()
    return pts

def pairs(L, row):                             #-- column c and the column that must equal it in limit L
    l = frozenset(ONLO[L]) if L in ONLO else key(L)[0]
    for c in row:
        k, g = key(c)
        if k & l: continue
        for d in row:
            if key(d) == (k | l, g): yield c, d

pts = run()
sec = mixed[mixed.index('-sec') + 1]
tc = sec[-1] in 'abcd'                         #-- triple-collinear sector: X = TC, Y = collinear inside it
NAME = {'X': 'TC' if tc else 'C5', 'Y': {'a': 'C6', 'c': 'C5', 'b': 'C56', 'd': 'C56'}.get(sec[-1], 'C6')}
VAR = {'S5': 'x(E5)', 'S6': 'x(E6)', 'X': 'x(rho5)', 'Y': 'x(rho6)', 'S': 'x(E)', 'C1': 'x(rho)', 'C2': '1-x(rho)'}
def phys(lab):                                 #-- 'S6+X' -> 'S6TC'; ONLO names ('SC1/2') unchanged
    name, *g = lab.split('/')
    if name in ONLO: return lab
    order = ['S5', 'S6', 'Y', 'X'] if tc else ['S5', 'S6', 'X', 'Y']      #-- e.g. S6C5TC: soft, inner collinear, TC
    parts = sorted(name.split('+'), key=lambda x: order.index(x) if x in order else -1)
    return ''.join(NAME.get(x, x) for x in parts) + ('/' + g[0] if g else '')
def edge(L):
    parts = ONLO[L] if L in ONLO else L.split('+')
    return ', '.join(f'{VAR[x]} -> 0' for x in parts)

print(f"# {' '.join(mixed)}: {len(pts)} points")
rows = {}
for pt in pts:
    for L, tab in pt['tabs'].items():
        if len(tab) < 10: continue
        for c, d in pairs(L, tab[0]):
            r = [abs(tab[i][c] - tab[i][d]) / abs(tab[i][d]) if tab[i][d] else float('nan') for i in range(10)]
            rows.setdefault(L, {}).setdefault(f'{phys(c)} = {phys(d)}', []).append(r)
head = f"  {'pair':<22}{'n':>3}" + ''.join(f'{d:>10.0e}' for d in DEPTHS) + f"{'slope':>8}  grade"
for L, prs in rows.items():
    print(f"\nlimit {phys(L)}   ({edge(L)})\n{head}")
    for pr, series in prs.items():
        med = [st.median([s[i] for s in series if s[i] == s[i] and s[i] > 0] or [float('nan')]) for i in range(10)]
        xs = [-(i + 1) for i in range(1, 6) if med[i] == med[i] and med[i] > 0]
        ys = [math.log10(med[i]) for i in range(1, 6) if med[i] == med[i] and med[i] > 0]
        slope = (sum((x - st.mean(xs)) * (y - st.mean(ys)) for x, y in zip(xs, ys)) / sum((x - st.mean(xs))**2 for x in xs)) if len(xs) > 1 else float('nan')
        print(f"  {pr:<22}{len(series):>3}" + ''.join(f'{med[i]:>10.1e}' for i in (1, 3, 5, 7)) + f"{slope:>8.2f}  {'ok' if slope >= 0.3 else 'FAIL'}")
unc = {}
for pt in pts:
    for L, (f, v) in pt['unc'].items(): unc.setdefault(L, []).append((f, v))
unc = {L: xs for L, xs in unc.items() if any(v != 'cut' for f, v in xs)}
if unc:
    print(f"\nlimits without own counterterm: x^k |H + counterterms| (k = variables at the edge)\n  {'limit':<22}{'n':>3}" + ''.join(f'{d:>10.0e}' for d in DEPTHS) + '  integrable')
    for L, xs in unc.items():
        med = [st.median([f[i] for f, v in xs if v != 'cut']) for i in range(4)]
        ni = sum(v == 'integrable' for f, v in xs); nn = sum(v != 'cut' for f, v in xs)
        print(f"  {phys(L):<22}{nn:>3}" + ''.join(f'{m:>10.1e}' for m in med) + f"  {ni}/{nn}")
