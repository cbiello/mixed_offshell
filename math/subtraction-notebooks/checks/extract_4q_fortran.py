#!/usr/bin/env python3
# usage: extract_4q_fortran.py <f90/src/IntSub/mod_hoppet_nnlo.f90> > kernels.m   -- four-quark kernels as Mathematica input
import re,sys
L=open(sys.argv[1]).read().split('\n')
for name in ['mySub_ns_qqb_func','mySub_ns_qqb_Lmu_func','mySub_ns_qq_func','mySub_ns_qq_Lmu_func']:
    i=[k for k,l in enumerate(L) if l.strip().startswith('function '+name)][0]
    j=[k for k in range(i,i+60) if L[k].strip().startswith('res = ') and 'zero' not in L[k] and 'res * x' not in L[k]][0]
    e=L[j].strip()[6:]
    while e.rstrip().endswith('&'):
        j+=1; e=e.rstrip()[:-1]+L[j].strip()
    print('f["%s"]=%s;'%(name,re.sub(r'(\d)\._dp',r'\1',e).replace('**','^')))
