"""Independent Hungarian checks of selected-allocation point-crossing prices."""
import argparse,json,random
from functools import lru_cache
from collections import defaultdict
from math import comb
from sections import sections
from geometry14 import freeze
from lift11 import atomic,assignment
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--sample',type=int,default=3000);a=ap.parse_args();D=a.D
cache=json.load(open(f'out/pointtriple14_D{D}_prices.json'));assert cache['D']==D
profiles=defaultdict(set)
for k,r,u,c,sec,e in sections(5,5):
    mult=(0,)*(k-u)+tuple(b for b,count in enumerate(c,1) for _ in range(count))
    profiles[(r,k,u,tuple(sec),freeze(e))].add(mult)
core=[];others=[]
for tk in cache['prices']:
    t=freeze(json.loads(tk));(core if t[:2] in ((5,15),(6,17)) else others).append(tk)
rng=random.Random(140700+D);chosen=sorted(core)+rng.sample(sorted(others),min(a.sample,len(others)))
matrices=0
@lru_cache(None)
def price(mult,h,j,sign):
    global matrices
    ds=[d for d,c in enumerate(h,1) for _ in range(c)];INF=1000000
    matrix=[[sign*comb(b,3)*(d==j) if d>=max(1,b) else INF for d in ds] for b in mult]
    cost,perm,_,_=assignment(matrix);matrices+=1
    return None if any(matrix[i][p]==INF for i,p in enumerate(perm)) else cost
checks=0
for number,tk in enumerate(chosen,1):
    r,k,h,u,sec,e=freeze(json.loads(tk));values=[]
    for j,sign in cache['objectives']:
        candidates=[v for mult in profiles[(r,k,u,sec,e)] if (v:=price(mult,h,j,sign)) is not None]
        assert candidates;values.append(min(candidates))
    assert values==cache['prices'][tk],tk;checks+=len(values)
    if number%500==0:print('POINTTRIPLE14_HUNGARIAN_PROGRESS',D,number,len(chosen),flush=True)
record=dict(D=D,status='EXACT_POINTTRIPLE14_HUNGARIAN_AGREEMENT',types=len(chosen),core_types=len(core),sample_types=len(chosen)-len(core),prices=checks,assignment_matrices=matrices,mismatches=0,method='independent integer Hungarian optimization over every matching section profile; no selected-allocation pricer calls')
atomic(f'out/pointtriple14_D{D}_independent.json',record);print('EXACT_POINTTRIPLE14_HUNGARIAN_AGREEMENT',record,flush=True)
