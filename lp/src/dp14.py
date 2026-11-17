"""Independent exhaustive histogram-state checks; no matching-pricer calls."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import argparse,json,random
from functools import lru_cache
from collections import defaultdict
from geometry14 import geometries,freeze
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--sample',type=int,default=3000);a=ap.parse_args();D=a.D
cache=json.load(open(f'out/lift14_D{D}_prices.json'));assert cache['D']==D
objs=[(tuple(s),sgn) for s,sgn in cache['objectives']]
profiles=defaultdict(set)
for g in geometries(D):
    s=freeze(g['signature']);profiles[s[:2]+s[3:]].add(tuple(sorted((m.bit_count() for m in g['masks']),reverse=True)))
core=[];others=[]
for tk in cache['prices']:
    t=freeze(json.loads(tk));(core if t[:2] in ((5,15),(6,17)) else others).append(tk)
rng=random.Random(140000+D);chosen=sorted(core)+rng.sample(sorted(others),min(a.sample,len(others)));states=0;objectives=0
@lru_cache(None)
def solve(mult,h,num):
    global states
    n=sum(h);active=objs[:num];costs={(b,d):tuple(sgn*b*(d in subset) for subset,sgn in active) for b in set(mult) for d in range(1,D+1)}
    @lru_cache(None)
    def rec(rem):
        pos=n-sum(rem)
        if pos==n or mult[pos]==0:return (0,)*num
        b=mult[pos];best=None
        for d,count in enumerate(rem,1):
            if not count or d<b:continue
            nxt=list(rem);nxt[d-1]-=1;tail=rec(tuple(nxt))
            if tail is None:continue
            candidate=tuple(x+y for x,y in zip(tail,costs[b,d]))
            best=candidate if best is None else tuple(min(x,y) for x,y in zip(best,candidate))
        return best
    result=rec(h);states+=rec.cache_info().currsize;return result
for number,tk in enumerate(chosen,1):
    t=freeze(json.loads(tk));r,k,h,u,sec,e=t;expected=cache['prices'][tk];best=None
    for mult in profiles[(r,k,u,sec,e)]:
        v=solve(mult,h,len(expected))
        if v is not None:best=v if best is None else tuple(min(x,y) for x,y in zip(best,v))
    assert best is not None and list(best)==expected,(tk,best,expected)
    objectives+=len(expected)
    if number%500==0:print('DP14_PROGRESS',D,number,len(chosen),flush=True)
record=dict(D=D,status='EXACT_DP14_VERIFIED',core_types=len(core),sample_types=len(chosen)-len(core),types=len(chosen),objectives=objectives,states=states,mismatches=0,method='exhaustive histogram-state recursion, independent of Hungarian',sample_seed=140000+D)
atomic(f'out/dp14_D{D}.json',record);print('EXACT_DP14_VERIFIED',record,flush=True)
