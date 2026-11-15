#!/usr/bin/env python3
# inverse-variance weighted average over seeds, using the error quoted in each pwg-NLO.top
import sys,glob,math
tag,hist=sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'rate'
v=[]
for f in glob.glob(f'{tag}/s*/pwg-NLO.top'):
    on=False
    for l in open(f):
        if l.startswith('#'): on=l.split()[1]==hist; continue
        w=l.split()
        if on and len(w)==4 and float(w[3])<100: v.append((float(w[2])*1000,float(w[3])*1000))
W=sum(1/e**2 for _,e in v); m=sum(x/e**2 for x,e in v)/W
chi=sum(((x-m)/e)**2 for x,e in v)/(len(v)-1)
print(f"{tag} {hist}: {m:.4f} +- {math.sqrt(1/W):.4f} fb  ({len(v)} seeds, chi2/dof {chi:.2f})")
