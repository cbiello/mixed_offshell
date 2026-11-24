#!/usr/bin/env python3
# usage: extract_kernels_fortran.py <f90/src/IntSub/mod_hoppet_nnlo.f90> > kernels.m
# Regular part of the gq, TCgq, TCaa, TCag kernels as Mathematica functions f["name"][x] (without the factor x of the HOPPET convention).
import re,sys
L=open(sys.argv[1]).read().split('\n')
def tr(e):
    e=re.sub(r'real\(dilog2\(([^()]*)\),kind=dp\)',r'PolyLog[2,\1]',e)
    e=re.sub(r'trilog\(([^()]*)\)',r'PolyLog[3,\1]',e)
    e=re.sub(r'\blog\(([^()]*)\)',r'Log[\1]',e)
    e=re.sub(r'(\d)\._dp',r'\1',e).replace('**','^')
    for a,b in [('one','1'),('two','2'),('three','3'),('four','4'),('zeta2','Zeta[2]'),('zeta3','Zeta[3]'),('pisq','Pi^2')]:
        e=re.sub(r'\b'+a+r'\b',b,e)
    return e
names=[l.split()[1].split('(')[0] for l in L if re.match(r'\s*function mySub_(gq|TC(gq|aa|ag))\w*\(',l)]
for name in names:
    i=[k for k,l in enumerate(L) if l.strip().startswith('function '+name+'(')][0]
    a=[k for k in range(i,min(i+40,len(L))) if 'case(cc_REAL,cc_REALVIRT)' in L[k]][0]
    b=[k for k in range(a,min(a+60,len(L))) if 'end select' in L[k]][0]
    body=' '.join(l.split('!')[0].strip() for l in L[a+1:b]).replace('& ','').replace('&','')
    st=[s.strip() for s in re.split(r'\s+(?=\w+\s*=[^=])',body) if s.strip()]
    loc=sorted(set(s.split('=')[0].strip() for s in st))
    print('f["%s"][x_]:=Module[{%s},%s;res];'%(name,','.join(loc),';'.join(tr(s) for s in st)))
