#!/usr/bin/env python3
# weighted average of a histogram of pwg-NLO.top over the seeds of a run
import sys,glob,re,math
tag,hist=sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'rate'
rows={}
for f in glob.glob(f'{tag}/s*/pwg-NLO.top'):
    on=False
    for l in open(f):
        if l.startswith('#'): on = l.split()[1]==hist; continue
        v=l.split()
        if on and len(v)==4:
            rows.setdefault((float(v[0]),float(v[1])),[]).append((float(v[2]),float(v[3])))
for k in sorted(rows):
    r=rows[k]; n=len(r)
    m=sum(x for x,_ in r)/n
    e=math.sqrt(sum((x-m)**2 for x,_ in r)/(n*(n-1))) if n>1 else r[0][1]
    print(f"{k[0]:8.1f} {k[1]:8.1f}  {m*1000:14.6f} +- {e*1000:10.6f} fb   ({n} seeds)")
