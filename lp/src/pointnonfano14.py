"""Exact m-free plane-point triples-minus-quadruples costs, split by global point degree.

Global costs are bounded by exact PG(2,q) direction-subset extrema. For q=5
and degree seven, the sum of C(b,3)-C(b,4) is between two and six.
Prices range over every admitted section profile and degree assignment.
"""
import argparse,json,hashlib
from pathlib import Path
from functools import lru_cache
from itertools import permutations
from collections import Counter,defaultdict
from math import comb
from sections import sections
from geometry14 import freeze
from lift11 import atomic,key
FAMILY='plane-point-triples-minus-quadruples-by-degree-v1'

def objectives(D):return [(j,sign) for j in range(3,D+1) for sign in (1,-1)]

@lru_cache(None)
def global_bounds(q,j):
    r=json.load(open(f'out/nonfano14_projective_q{q}.json'))
    assert r['q']==q and r['D']>=j and r['status']=='EXACT_NONFANO14_PROJECTIVE_ENUMERATED'
    v=r['records'][j];assert v['degree']==j
    if q==5:assert (v['minimum'],v['maximum'])=={3:(0,1),4:(0,3),5:(0,5),6:(0,5),7:(2,6)}[j]
    return v['minimum'],v['maximum']


@lru_cache(None)
def profiles(q):
    out=defaultdict(set)
    for k,r,u,c,t,e in sections(q,5):
        p=(k-u,)+tuple(c)
        assert sum(p)==k and sum(b*x for b,x in enumerate(p))==q*r
        assert sum(comb(b,2)*x for b,x in enumerate(p))==comb(r,2)
        out[(r,k,u,tuple(t),freeze(e))].add(p)
    return {k:tuple(sorted(v)) for k,v in out.items()}

@lru_cache(None)
def bounds(profile,h,j):
    """Enumerate selected allocations only at multiplicities b>=2.

    Multiplicities zero and one have zero cost and can be combined. The
    remaining degree multiset fits the remaining slots iff all tail (Hall)
    inequalities hold. Thus every enumerated allocation extends to an actual
    assignment, and every assignment is represented by one allocation.
    """
    r=len(profile)-1
    if sum(profile)!=sum(h) or any(sum(profile[b:])>sum(h[b-1:]) for b in range(2,r+1)):
        return None
    selected=h[j-1];low=sum(profile[:2]);x=[0]*(r+1);values=[]
    rest_tail={b:sum(h[b-1:])-(selected if j>=b else 0) for b in range(2,r+1)}
    groups=[b for b in range(2,min(j,r)+1) if profile[b]]
    def rec(pos,count,cost):
        if pos==len(groups):
            if not 0<=selected-count<=low:return
            if all(sum(profile[b:])-sum(x[b:])<=rest_tail[b] for b in range(2,r+1)):
                values.append(cost)
            return
        b=groups[pos]
        for v in range(min(profile[b],selected-count)+1):
            x[b]=v;rec(pos+1,count+v,cost+(comb(b,3)-comb(b,4))*v)
        x[b]=0
    rec(0,0,0)
    assert values,('compatible histogram failed allocation',profile,h,j)
    return min(values),max(values)

def prices(q,D,t):
    r,k,h,u,sec,e=t
    admitted=profiles(q)[(r,k,u,sec,e)]
    result=[]
    for j in range(3,D+1):
        vals=[v for p in admitted if (v:=bounds(p,h,j)) is not None]
        assert vals,('no admitted assignment',q,D,t,j)
        result.extend((min(v[0] for v in vals),-max(v[1] for v in vals)))
    return result

def validate(D):
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    results=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);q=ag.q;S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
        degree=Counter(p for l in W for p in ag.lines[l][1]);assert max(degree.values())<=D
        totals=Counter();point_totals=Counter();checks=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=sorted(S&ps);r=len(ls);k=len(pts)
            h=tuple(sum(degree[p]==j for p in pts) for j in range(1,D+1))
            mult={p:sum(p in ag.lines[l][1] for l in ls) for p in pts}
            u=sum(b>0 for b in mult.values())
            local=[(di,lp) for di,lp in ag.lines if lp<=ps];marked={ag.lines[l][0] for l in ls}
            sec=tuple(sum(len(lp&S)==j for _,lp in local) for j in range(1,q+1))
            e=tuple(sorted(tuple(sorted(len(lp&S) for dd,lp in local if dd==di)) for di in {di for di,_ in local}-marked)) if r>=5 and r<q+1 else None
            t=(r,k,h,u,sec,e);actual_profile=tuple(sum(b==i for b in mult.values()) for i in range(r+1))
            assert actual_profile in profiles(q)[(r,k,u,sec,e)]
            for p in pts:point_totals[p]+=(comb(mult[p],3)-comb(mult[p],4))
            mu=prices(q,D,t)
            for z,(j,sign) in enumerate(objectives(D)):
                actual=sum((comb(mult[p],3)-comb(mult[p],4)) for p in pts if degree[p]==j)
                assert sign*actual>=mu[z],(fn,t,j,sign,actual,mu[z])
                if sign==1:totals[j]+=actual
                checks+=1
        for j in range(3,D+1):
            lo,hi=global_bounds(q,j);m=sum(degree[p]==j for p in S)
            assert lo*m<=totals[j]<=hi*m
        for p in S:
            if degree[p]>=3:
                lo,hi=global_bounds(q,degree[p]);assert lo<=point_totals[p]<=hi
        results.append(dict(file=fn,checks=checks,global_bound_checks=D-2,point_bound_checks=sum(degree[p]>=3 for p in S),violations=0));print('POINTNONFANO14_VALIDATE',D,results[-1],flush=True)
    brute=0
    for mult in ((0,1,2,3),(1,1,2,2),(0,0,2,4),(0,1,1,2)):
        p=tuple(mult.count(b) for b in range(max(mult)+1))
        for ds in ((1,2,2,3),(1,1,2,3),(2,2,3,4),(1,2,3,4)):
            h=tuple(ds.count(j) for j in range(1,5))
            for j in range(3,5):
                vals=[sum((comb(b,3)-comb(b,4)) for b,d in zip(mult,perm) if d==j) for perm in set(permutations(ds)) if all(d>=max(1,b) for d,b in zip(perm,mult))]
                assert bounds(p,h,j)==((min(vals),max(vals)) if vals else None);brute+=1
    atomic(f'out/pointnonfano14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=results,brute_cases=brute))

def table(D,fresh=False):
    path=Path(f'out/pointnonfano14_D{D}_prices.json');partial=path.with_name(path.stem+'_partial.json')
    sh=hashlib.sha256(Path('out/sections3_q5_r5.json').read_bytes()).hexdigest()
    gh=hashlib.sha256(Path('out/nonfano14_projective_q5.json').read_bytes()).hexdigest()
    if fresh:profiles.cache_clear();bounds.cache_clear();global_bounds.cache_clear()
    for j in range(3,D+1):global_bounds(5,j)
    source=path if path.exists() else partial
    result=json.load(source.open()) if source.exists() and not fresh else dict(q=5,D=D,family=FAMILY,section_sha256=sh,projective_sha256=gh,objectives=objectives(D),prices={})
    assert (result['q'],result['D'],result['family'],result['section_sha256'],result['projective_sha256'])==(5,D,FAMILY,sh,gh)
    types=[freeze(t) for t in json.load(open(f'out/types14_D{D}.json'))['types']]
    for number,t in enumerate(types,1):
        tk=key(t)
        if tk not in result['prices']:result['prices'][tk]=prices(5,D,t)
        if number%10000==0:
            atomic(partial.with_name(partial.stem+'_fresh.json') if fresh else partial,result)
            print('POINTNONFANO14_PRICE',D,number,len(types),flush=True)
    if fresh and path.exists():
        old=json.load(path.open());assert old['prices']==result['prices'] and old['section_sha256']==sh and old['projective_sha256']==gh
        atomic(f'out/pointnonfano14_D{D}_fresh_check.json',dict(D=D,types=len(types),prices=sum(map(len,result['prices'].values())),mismatches=0,method='fresh exhaustive selected-multiplicity allocations before reading old prices'))
    atomic(path,result);return result

def add_cuts(lp,info,types,D):
    v=json.load(open(f'out/pointnonfano14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS')
    independent=json.load(open(f'out/pointnonfano14_D{D}_independent.json'))
    assert independent['status']=='EXACT_POINTNONFANO14_HUNGARIAN_AGREEMENT' and independent['mismatches']==0
    cache=json.load(open(f'out/pointnonfano14_D{D}_prices.json'));assert (cache['q'],cache['D'],cache['family'])==(5,D,FAMILY)
    assert cache['section_sha256']==hashlib.sha256(Path('out/sections3_q5_r5.json').read_bytes()).hexdigest()
    assert cache['projective_sha256']==hashlib.sha256(Path('out/nonfano14_projective_q5.json').read_bytes()).hexdigest()
    assert cache['objectives']==json.loads(key(objectives(D)))
    ordered_prices=[cache['prices'][key(t)] for t in types]
    for z,(j,sign) in enumerate(objectives(D)):
        row=[(lp.col[f'N{i}'],ordered_prices[i][z]) for i,t in enumerate(types)]
        lo,hi=global_bounds(5,j)
        row.append((lp.col[f'm{j}'],-hi if sign==1 else lo))
        lp.add(row,'<',0,f'XIV point triples-minus-quadruples degree {j} sign {sign}')
    print('POINTNONFANO14_CUTS',D,len(objectives(D)),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--fresh',action='store_true');a=ap.parse_args()
    validate(a.D);r=table(a.D,a.fresh);print('POINTNONFANO14_COMPLETE',a.D,len(r['prices']),flush=True)
