"""Check literal row equivalence when absent degree labels duplicate objectives."""
import json
from itertools import combinations_with_replacement
from collections import Counter
from planecount5 import LP
from geometry14 import freeze
from lift11 import atomic
import lift14,nonlinear14,triple14
D=6
alltypes=[freeze(t) for t in json.load(open('out/types14_D6.json'))['types']]
types=[next(t for t in alltypes if t[:2]==(5,15) and sum(v>0 for v in t[2])==s) for s in (2,4)]
dms=list(combinations_with_replacement(range(1,D+1),5))
ys=[(t[:3],di) for t in types for di,dm in enumerate(dms) if all(c<=t[2][d-1] for d,c in Counter(dm).items())]
info={'_Y':ys,'_dms':dms}
def model():
    lp=LP()
    for i in range(len(types)):lp.var(f'N{i}')
    for i in range(len(ys)):lp.var(f'G{i}')
    return lp
def key(row):
    co,b,s,_=row
    return tuple(sorted(co.items())),b,s
results=[]
for mod in (lift14,nonlinear14,triple14):
    full=model();small=model()
    mod.add_cuts(full,info,types,D,'all',False)
    mod.add_cuts(small,info,types,D,'all',True)
    original={key(row) for row in full.rows};kept={key(row) for row in small.rows}
    assert kept<=original
    assert all(not row[0] and row[1]==0 or key(row) in kept for row in full.rows)
    results.append(dict(family=mod.FAMILY,full_rows=len(full.rows),compact_rows=len(small.rows),mismatches=0))
atomic('out/compact14_check.json',dict(status='EXACT_LITERAL_ROW_EQUIVALENCE',D=D,types=types,results=results))
print('EXACT_COMPACT14_ROW_EQUIVALENCE',results,flush=True)
