"""Universal pair-of-degree incidence inequalities, with exact integer pricing."""
import os
for _v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_v]='1'
import json,hashlib
from collections import defaultdict,Counter
from functools import lru_cache
from itertools import combinations,permutations
from pathlib import Path
from lift11 import assignment,INF,atomic,key
from plane10 import geometries,freeze
FAMILY='signed-pair-subsets-v1'
PATH=Path('out/lift13_pair_prices.json')

def objectives(D):return [(i,j,s) for i,j in combinations(range(1,D+1),2) for s in (1,-1)]
@lru_cache(None)
def price_profile(mult,h,i,j,sgn):
    ds=[d for d,c in enumerate(h,1) for _ in range(c)]
    a=[[sgn*b*(d in (i,j)) if d>=max(1,b) else INF for d in ds] for b in mult]
    cost,perm,_,_=assignment(a)
    return None if any(a[k][perm[k]]==INF for k in range(len(a))) else cost

@lru_cache(None)
def profiles():
    out=defaultdict(set)
    for g in geometries():
        s=freeze(g['signature']);out[s[:2]+s[3:]].add(tuple(sorted(m.bit_count() for m in g['masks'])))
    return out

def prices(t):
    r,k,h,u,sec,e=t;pp=profiles()[(r,k,u,sec,e)]
    vals=[]
    for i,j,sgn in objectives(len(h)):
        v=[price_profile(p,h,i,j,sgn) for p in pp];v=[x for x in v if x is not None]
        assert v,('no assignment',t)
        vals.append(min(v))
    return vals

def admitted(t):return t[:2] in ((5,15),(6,17))
def table(types,fresh=False):
    gh=hashlib.sha256(key(geometries()).encode()).hexdigest()
    if fresh:price_profile.cache_clear();profiles.cache_clear()
    old=None if fresh or not PATH.exists() else json.load(PATH.open())
    if old:assert old['family']==FAMILY and old['geometry_sha256']==gh
    out=old or dict(family=FAMILY,geometry_sha256=gh,prices={})
    for t in types:
        if admitted(t) and key(t) not in out['prices']:out['prices'][key(t)]=prices(t)
    if fresh and PATH.exists():
        old=json.load(PATH.open());assert old['family']==FAMILY and old['geometry_sha256']==gh
        assert all(old['prices'].get(k)==v for k,v in out['prices'].items())
        atomic('out/lift13_pair_fresh_check.json',dict(types=len(out['prices']),mismatches=0,family=FAMILY,price_inputs='geometry only before reading old prices'))
        print('FRESH_PAIR13_PRICES_VERIFIED',len(out['prices']),flush=True)
    atomic(PATH,out)
    return out

def add_cuts(lp,info,types):
    val=json.load(open('out/lift13_validation.json'));assert val['status']=='ZERO_VIOLATIONS' and val['family']==FAMILY
    cache=table(types);groups=defaultdict(list);gvars=defaultdict(list)
    for i,t in enumerate(types):
        if admitted(t):groups[t[:3]].append(i)
    for g,(a,di) in enumerate(info['_Y']):
        if admitted(a):gvars[a].append((lp.col[f'G{g}'],info['_dms'][di]))
    count=0
    for a,inds in groups.items():
        for z,(j,k,s) in enumerate(objectives(5)):
            row=[(lp.col[f'N{i}'],cache['prices'][key(types[i])][z]) for i in inds]
            row += [(col,-s*(dm.count(j)+dm.count(k))) for col,dm in gvars[a]]
            lp.add(row,'<',0,f'XIII pair {key(a)} {j} {k} {s}');count+=1
    print('PAIR13_CUTS',count,'cached_types',len(cache['prices']),flush=True)

def validate():
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    results=[]
    for fn in ['solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json']:
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
        deg=Counter(p for l in W for p in ag.lines[l][1]);D=max(deg.values());checks=0;fullchecks=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=sorted(S&ps)
            mult=tuple(sorted(sum(p in ag.lines[l][1] for l in ls) for p in pts))
            h=tuple(sum(deg[p]==d for p in pts) for d in range(1,D+1));full_mu=None
            if ag.q==5 and len(ls)>=5:
                local=[(d,lp) for d,lp in ag.lines if lp<=ps];marked={ag.lines[l][0] for l in ls}
                sec=tuple(sum(len(lp&S)==d for _,lp in local) for d in range(1,6))
                e=tuple(sorted(tuple(sorted(len(lp&S) for dd,lp in local if dd==d)) for d in {d for d,_ in local}-marked)) if len(ls)==5 else None
                union=set().union(*(ag.lines[l][1] for l in ls))
                full_mu=prices((len(ls),len(pts),h,len(union),sec,e))
            for z,(i,j,s) in enumerate(objectives(D)):
                incidence=sum(sum(deg[p] in (i,j) for p in ag.lines[l][1]) for l in ls)
                mu=price_profile(mult,h,i,j,s);assert mu is not None and s*incidence>=mu
                if full_mu is not None:assert s*incidence>=full_mu[z];fullchecks+=1
                checks+=1
        rec=dict(file=fn,checks=checks,full_type_checks=fullchecks,violations=0);results.append(rec);print(rec,flush=True)
    brute=0
    for mult in [(0,1,1,2),(1,2,2,3),(0,0,1,3),(1,1,2,4)]:
        for ds in [(1,2,2,3),(1,1,2,3),(2,2,3,3),(1,2,3,4)]:
            h=tuple(ds.count(d) for d in range(1,5))
            for i,j,s in objectives(4):
                costs=[sum(s*b*(d in (i,j)) for b,d in zip(mult,p)) for p in set(permutations(ds)) if all(d>=max(1,b) for b,d in zip(mult,p))]
                assert price_profile(mult,h,i,j,s)==(min(costs) if costs else None);brute+=1
    atomic('out/lift13_validation.json',dict(family=FAMILY,status='ZERO_VIOLATIONS',results=results,brute_cases=brute))
    print('PAIR13_VALIDATION_OK',brute,flush=True)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--price',action='store_true');a=ap.parse_args();validate()
    if a.price:
        import planecount6 as p6
        t=p6.hist_types(5,5,46);out=table(t);print('PAIR13_UNIVERSAL_PRICES',len(out['prices']),flush=True)
