"""Independent histogram-state exhaustive pricer; no Hungarian or lift13 pricing calls."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import json
from itertools import combinations
from functools import lru_cache
from plane10 import geometries,freeze
from lift11 import atomic,key
geos=geometries();cached=json.load(open('out/lift13_pair_prices.json'));answers={};states=0
for tk,expected in cached['prices'].items():
    r,k,h,u,sec,e=freeze(json.loads(tk));objs=[(i,j,s) for i,j in combinations(range(1,6),2) for s in (1,-1)]
    matching=set()
    for g in geos:
        sig=freeze(g['signature'])
        if sig[:2]+sig[3:]==(r,k,u,sec,e):matching.add(tuple(sorted(mask.bit_count() for mask in g['masks'])))
    minima=None
    for mult in matching:
        @lru_cache(None)
        def exhaustive(remaining):
            pos=k-sum(remaining)
            if pos==k:return (0,)*len(objs)
            b=mult[pos];best=None
            for d,count in enumerate(remaining,1):
                if not count or d<max(1,b):continue
                nxt=list(remaining);nxt[d-1]-=1;v=exhaustive(tuple(nxt))
                if v is None:continue
                v=tuple(x+s*b*(d==i or d==j) for x,(i,j,s) in zip(v,objs))
                best=v if best is None else tuple(min(x,y) for x,y in zip(best,v))
            return best
        got=exhaustive(h);states+=exhaustive.cache_info().currsize
        if got is not None:minima=got if minima is None else tuple(min(x,y) for x,y in zip(minima,got))
    assert minima is not None and list(minima)==expected,tk
    answers[tk]=list(minima)
atomic('out/lift13_dp_prices.json',dict(method='exhaustive histogram-state dynamic programming; no Hungarian',prices=answers))
atomic('out/lift13_dp_check.json',dict(types=len(answers),objectives=len(answers)*20,states=states,mismatches=0,status='EXACT_DP_PAIR_PRICES_VERIFIED'))
print('EXACT_DP_PAIR_PRICES_VERIFIED',len(answers),'objectives',len(answers)*20,'states',states,'mismatches 0',flush=True)
