"""Independent ordinary-colouring atlas and explicit sorted degree allocation."""
import json,random,argparse
from collections import Counter
from math import comb
from functools import lru_cache
from geometry14 import geometries,freeze
from lift11 import atomic,key
import extra14 as e
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--sample',type=int,default=3000);a=ap.parse_args()
gs=geometries(a.D);inside={tuple(sorted(m for m in g['masks'] if m)) for g in gs if g['signature'][0]==6 and g['signature'][3]==17}
atlases={}
for masks in inside:
    groups=[sum(1<<p for p,m in enumerate(masks) if m.bit_count()==b) for b in range(7)]
    lines=[sum(1<<p for p,m in enumerate(masks) if m>>i&1) for i in range(6)]
    direct={}
    for bits in range(1<<len(masks)):
        profile=tuple((bits&g).bit_count() for g in groups);counts=[(bits&l).bit_count() for l in lines]
        values=[sum(comb(c,power) for c in counts) for power in (2,3)]
        if profile not in direct:direct[profile]=[[v,v] for v in values]
        else:
            for p,v in enumerate(values):direct[profile][p]=[min(v,direct[profile][p][0]),max(v,direct[profile][p][1])]
    for power in (2,3):
        primary={p:[lo,hi] for p,lo,hi,_,_ in (e.pair if power==2 else e.triple).extrema(masks,6)}
        assert primary=={p:vals[power-2] for p,vals in direct.items()}
    atlases[masks]=direct;print('EXTRA14_INDEPENDENT_ATLAS',len(masks),len(direct),flush=True)
@lru_cache(None)
def independent(masks,h,s,power):
    extra=masks.count(0);inside=tuple(m for m in masks if m)
    ds=sorted(d for d,n in enumerate(h,1) for _ in range(n) if d in s)
    dc=sorted(d for d,n in enumerate(h,1) for _ in range(n) if d not in s)
    total=Counter(m.bit_count() for m in inside);low=None;high=None
    for profile,values in atlases[inside].items():
        v0=len(ds)-sum(profile)
        if not 0<=v0<=extra:continue
        selected=[1]*v0+[max(1,b) for b,n in enumerate(profile) for _ in range(n)]
        other=[1]*(extra-v0)+[max(1,b) for b in range(7) for _ in range(total[b]-profile[b])]
        selected.sort();other.sort()
        assert len(selected)==len(ds) and len(other)==len(dc)
        if any(b>d for b,d in zip(selected,ds)) or any(b>d for b,d in zip(other,dc)):continue
        lo,hi=values[power-2];low=lo if low is None else min(lo,low);high=hi if high is None else max(hi,high)
    assert low is not None
    return low,-high
cache=json.load(open(f'out/extra14_D{a.D}_prices.json'));types=list(cache['prices']);random.Random(140071).shuffle(types)
chosen=types[:a.sample];chosen+= [s for s in types if json.loads(s)[1] in (17,25) and s not in set(chosen)]
geometry_index={}
for g in gs:
    t=freeze(g['signature']);geometry_index.setdefault(t[:2]+t[3:],set()).add(tuple(sorted(g['masks'])))
checks=0
for n,ks in enumerate(chosen,1):
    t=freeze(json.loads(ks));masks=geometry_index[t[:2]+t[3:]]
    for z,(s,sg,p) in enumerate(e.objectives(a.D)):
        value=min(independent(m,t[2],s,p)[0 if sg==1 else 1] for m in masks)
        assert value==cache['prices'][ks][z],(t,s,sg,p,value,cache['prices'][ks][z]);checks+=1
    if n%100==0:print('EXTRA14_INDEPENDENT_PRICE',n,len(chosen),checks,flush=True)
atomic(f'out/extra14_D{a.D}_independent.json',dict(D=a.D,family=e.FAMILY,status='EXACT_EXTRA14_CHECKED',types=len(chosen),prices=checks,mismatches=0,method='ordinary subsets; explicit sorted selected/complement degree assignments including all extras'))
