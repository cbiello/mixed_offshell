#!/bin/bash
# usage: cp.sh "<corr> <ch> <sec>" [nc0 nc1]  -- W- in pp against W+ in pbar pbar (CP conjugate), 4 seeds each
mkdir -p tests/tmp_cp; set -- $1 ${2:-100000} ${3:-200000}; tag=$1_$2_$3
run() { ./mixed_ll -boson $1 $2 -corr $3 -ch $4 -sec $5 -vegasIt0 3 -vegasNc0 $6 -vegasIt1 3 -vegasNc1 $7 -seed $8 -pref cp$1$8_ -histodir tests/tmp_cp < /dev/null 2>&1 | grep -a "# Result\|# Error\|STOP" | tr -s ' ' | tr '\n' ' '; echo; }
export -f run
for s in 11 12 13 14; do run wm "" $1 $2 $3 $4 $5 $s > /tmp/cp_$$_wm_$s & run wp "-cptest" $1 $2 $3 $4 $5 $s > /tmp/cp_$$_wp_$s & done; wait
python3 - "$tag" /tmp/cp_$$_ <<'P'
import sys,glob,math
tag,pre=sys.argv[1:]
def avg(b):
    v=[]
    for f in glob.glob(pre+b+'_*'):
        t=open(f).read().split()
        if 'Result:' in t: v.append((float(t[t.index('Result:')+1]),float(t[t.index('Error:')+1])))
        else: print(tag,b,open(f).read().strip()[:100])
    if not v: return None
    w=[1/e**2 for _,e in v]; m=sum(x*wi for (x,_),wi in zip(v,w))/sum(w); return m,1/math.sqrt(sum(w))
a,b=avg('wm'),avg('wp')
if a and b:
    d=a[0]-b[0]; e=math.hypot(a[1],b[1])
    print(f"{tag:28s} W-: {a[0]:12.6f} +- {a[1]:9.6f}   W+ pbar pbar: {b[0]:12.6f} +- {b[1]:9.6f}   diff/err = {d/e:5.2f}   rel = {d/abs(a[0]) if a[0] else 0:8.1e}")
P
rm -rf /tmp/cp_$$_* grid/cpw* tests/tmp_cp
