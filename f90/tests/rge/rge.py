#!/usr/bin/env python3
# usage: rge.py <output dir of rgetest.sh>
# LHS: slope c1 of the mixed sectors, from the scale-variation weights (mu/2, 2 mu) of one run.
# RHS: in -rge runs the same weights hold PDFs -+ delta with all logs at the central scale, so (w3 - w2)/2 is the result with one PDF -> delta.
import sys,glob,math,re
ln4=math.log(4); c1=r0=0
rows={}
for f in sorted(glob.glob(sys.argv[1]+'/rg*')):
    tag,ch,sec=re.match(r'rg([A-Z])_([a-z]+)_(\w+?)__',f.split('/')[-1]).groups()
    L=open(f).read().split('\n'); w=L[[k for k,l in enumerate(L) if 'Obs=rate' in l][0]+2].split()
    cen,lo,hi,err=float(w[3]),float(w[5]),float(w[7]),float(w[4])
    x=(hi-lo)/(2*ln4) if tag=='L' else (hi-lo)/2
    rows.setdefault((tag,sec),[0,0,0]); r=rows[(tag,sec)]; r[0]+=x; r[1]+=cen; r[2]=math.hypot(r[2],err)
for (tag,sec),r in sorted(rows.items()):
    print('%s %-14s %10.5f    (value at mu0 %10.4f +- %.4f)'%(tag,sec,r[0],r[1],r[2]))
    if tag=='L': c1+=r[0]
    elif tag=='V': c1-=r[0]
    else: r0+=r[0]
print('slope of the mixed sectors (with -V for the two-loop sector) c1 = %.5f    prediction -r0 = %.5f    c1 + r0 = %.5f'%(c1,-r0,c1+r0))
