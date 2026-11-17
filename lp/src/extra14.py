"""Exact r=6,u=17 within-line pair/triple costs, including extra points.

Conditional on the separately replayed maximal-r reduction. Full D-specific
u=17 types are priced before any prefix filtering. Zero-incidence extra points
are assigned by exact selected/complement Hall inequalities, not discarded.
"""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import json,hashlib
from pathlib import Path
from functools import lru_cache
from collections import defaultdict,Counter
from itertools import combinations,permutations
from math import comb
import numpy as np
import nonlinear14 as pair
import triple14 as triple
from geometry14 import geometries,freeze
from lift11 import atomic,key
FAMILY='r6-u17-extra-points-pair-triple-v1'
def objectives(D):
    return [(s,sign,power) for size in (1,2) for s in combinations(range(1,D+1),size) for power in (2,3) for sign in (1,-1)]
@lru_cache(None)
def atlas(masks,r,power):
    inside=tuple(m for m in masks if m)
    rows=(pair if power==2 else triple).extrema(inside,r)
    profiles=np.array([p for p,*_ in rows],dtype=np.int64)
    cumulative=np.cumsum(profiles[:,::-1],axis=1)[:,::-1]
    totals=np.array([sum(m.bit_count()==b for m in inside) for b in range(r+1)],dtype=np.int64)
    return profiles,cumulative,totals,np.array([[lo,hi] for _,lo,hi,_,_ in rows],dtype=np.int64)
@lru_cache(None)
def local(masks,r,h,subset,sign,power):
    # Strip zero masks for atlas reuse; retain their exact number separately.
    extra=masks.count(0);inside=tuple(m for m in masks if m)
    profiles,cum,totals,values=atlas(inside,r,power)
    selected=np.array([v if j in subset else 0 for j,v in enumerate(h,1)],dtype=np.int64)
    other=np.array(h,dtype=np.int64)-selected
    s=int(selected.sum());v0=s-cum[:,0]
    ok=(v0>=0)&(v0<=extra)
    for b in range(1,r+1):
        ok &= cum[:,b]<=selected[b-1:].sum()
        ok &= int(totals[b:].sum())-cum[:,b]<=other[b-1:].sum()
    if not np.any(ok):return None
    return int(values[ok,0].min()) if sign==1 else -int(values[ok,1].max())
@lru_cache(None)
def matching(D,t):
    wanted=t[:2]+t[3:]
    found={tuple(sorted(g['masks'])) for g in geometries(D) if (s:=freeze(g['signature']))[:2]+s[3:]==wanted}
    assert found
    return tuple(sorted(found))
def price(D,t,s,sign,power):
    vals=[local(m,t[0],t[2],tuple(s),sign,power) for m in matching(D,t)]
    vals=[v for v in vals if v is not None];assert vals,(t,s,sign,power)
    return min(vals)
def table(D,fresh=False):
    types=[freeze(t) for t in json.load(open(f'out/types14_D{D}.json'))['types'] if t[0]==6 and t[3]==17]
    path=Path(f'out/extra14_D{D}_prices.json');gh=hashlib.sha256(key(geometries(D)).encode()).hexdigest()
    if fresh:
        atlas.cache_clear();local.cache_clear();matching.cache_clear();pair.extrema.cache_clear();triple.extrema.cache_clear()
    obj=objectives(D);result=dict(D=D,family=FAMILY,geometry_sha256=gh,objectives=obj,prices={})
    for number,t in enumerate(types,1):
        result['prices'][key(t)]=[price(D,t,s,sg,p) for s,sg,p in obj]
        if number%500==0:print('EXTRA14_PRICE',D,number,len(types),flush=True)
    if fresh and path.exists():
        old=json.load(path.open());assert old['D']==D and old['family']==FAMILY and old['geometry_sha256']==gh and old['prices']==result['prices']
        atomic(f'out/extra14_D{D}_fresh_check.json',dict(D=D,types=len(types),prices=len(types)*len(obj),mismatches=0,method='fresh full inside-colouring atlas plus exact extra-point allocation'))
    atomic(path,result);return result

def validate(D):
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    records=[];obj=objectives(D)
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S);degree=Counter(p for l in W for p in ag.lines[l][1]);assert max(degree.values())<=D
        checks=extras=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=sorted(S&ps);r=len(ls)
            masks=tuple(sorted(sum(1<<i for i,l in enumerate(ls) if p in ag.lines[l][1]) for p in pts));h=tuple(sum(degree[p]==d for p in pts) for d in range(1,D+1))
            extras+=masks.count(0)
            for s,sign,power in obj:
                actual=sign*sum(comb(sum(degree[p] in s for p in ag.lines[l][1]),power) for l in ls)
                mu=local(masks,r,h,s,sign,power);assert mu is not None and actual>=mu
                checks+=1
        rec=dict(file=fn,checks=checks,extra_point_incidences=extras,violations=0);records.append(rec);print('EXTRA14_VALIDATE',rec,flush=True)
    brute=0
    for masks in ((0,0,1,2,3),(0,1,1,2,3),(0,0,1,3,3)):
        for ds in ((1,1,2,2,3),(1,2,2,3,3),(2,2,2,3,3)):
            h=tuple(ds.count(d) for d in range(1,4))
            for s,sign,power in objectives(3):
                vals=[sign*sum(comb(sum(d in s and bool(m>>i&1) for d,m in zip(perm,masks)),power) for i in range(2)) for perm in set(permutations(ds)) if all(d>=max(1,m.bit_count()) for d,m in zip(perm,masks))]
                assert local(tuple(sorted(masks)),2,h,s,sign,power)==(min(vals) if vals else None);brute+=1
    atomic(f'out/extra14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=records,brute_cases=brute))

def add_cuts(lp,info,types,D,compact=True):
    assert all(t[0]!=6 or t[3]==17 for t in types),'Requires proved maximal-r restriction'
    v=json.load(open(f'out/extra14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS')
    independent=json.load(open(f'out/extra14_D{D}_independent.json'));assert independent['status']=='EXACT_EXTRA14_CHECKED'
    cache=json.load(open(f'out/extra14_D{D}_prices.json'));assert cache['D']==D and cache['family']==FAMILY and cache['geometry_sha256']==hashlib.sha256(key(geometries(D)).encode()).hexdigest()
    obj=objectives(D);assert cache['objectives']==json.loads(key(obj));groups=defaultdict(list);gv=defaultdict(list)
    for i,t in enumerate(types):
        if t[0]==6 and t[1]>17:groups[t[:3]].append(i)
    for g,(a,di) in enumerate(info['_Y']):
        if a[0]==6 and a[1]>17:gv[a].append((lp.col[f'G{g}'],info['_dms'][di]))
    count=0
    for a,indices in groups.items():
        seen={};support={j for j,h in enumerate(a[2],1) if h}
        for z,(s,sign,power) in enumerate(obj):
            effective=tuple(j for j in s if j in support);signature=(effective,sign,power)
            if not effective:
                assert all(cache['prices'][key(types[i])][z]==0 for i in indices);continue
            if compact and signature in seen:
                old=seen[signature];assert all(cache['prices'][key(types[i])][z]==cache['prices'][key(types[i])][old] for i in indices);continue
            seen[signature]=z
            row=[(lp.col[f'N{i}'],cache['prices'][key(types[i])][z]) for i in indices]
            row += [(col,-sign*comb(sum(d in s for d in dm),power)) for col,dm in gv[a]]
            lp.add(row,'<',0,f'XIV extra {key(a)} {key(s)} {sign} {power}');count+=1
    print('EXTRA14_CUTS',D,count,flush=True)
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--fresh',action='store_true');a=ap.parse_args()
    validate(a.D);r=table(a.D,a.fresh);print('EXTRA14_COMPLETE',a.D,len(r['prices']),flush=True)
