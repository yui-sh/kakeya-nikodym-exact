"""Joint prices for point crossing-pairs and rich secants at the same degree.

Global crossing-pair totals are exactly C(j,2)m_j. Rich-secant bounds are
conditional on the proved cap four, as in pointrich14. Pricing the sum before
minimization is stronger than adding separately minimized inequalities.
"""
import json,hashlib,argparse
from functools import lru_cache
from collections import Counter
from math import comb
from pathlib import Path
from itertools import permutations
from geometry14 import freeze
from lift11 import atomic,key
from pointrich14 import profiles,global_bounds
FAMILY='joint-point-pairs-rich-secants-v1'
LAMBDAS=(-2,-1,1,2)
def objectives(D):return [(j,lam,sign) for j in range(3,D+1) for lam in LAMBDAS for sign in (1,-1)]
@lru_cache(None)
def bounds(profile,h,j):
    r=len(profile)-1
    if sum(profile)!=sum(h) or any(sum(profile[b:])>sum(h[b-1:]) for b in range(2,r+1)):return None
    selected=h[j-1];low=sum(profile[:2]);x=[0]*(r+1);values=set()
    rest={b:sum(h[b-1:])-(selected if j>=b else 0) for b in range(2,r+1)}
    groups=[b for b in range(2,min(j,r)+1) if profile[b]]
    def rec(pos,count,pairs,rich):
        if pos==len(groups):
            if 0<=selected-count<=low and all(sum(profile[b:])-sum(x[b:])<=rest[b] for b in range(2,r+1)):values.add((pairs,rich))
            return
        b=groups[pos]
        for v in range(min(profile[b],selected-count)+1):
            x[b]=v;rec(pos+1,count+v,pairs+comb(b,2)*v,rich+int(b>=3)*v)
        x[b]=0
    rec(0,0,0,0);assert values
    return tuple((min(lam*p+r for p,r in values),max(lam*p+r for p,r in values)) for lam in LAMBDAS)
def prices(q,D,t):
    r,k,h,u,sec,e=t;admitted=profiles(q)[(r,k,u,sec,e)];result=[]
    for j in range(3,D+1):
        vals=[v for p in admitted if (v:=bounds(p,h,j)) is not None];assert vals
        for k in range(len(LAMBDAS)):result.extend((min(v[k][0] for v in vals),-max(v[k][1] for v in vals)))
    return result

def validate(D):
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    records=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);q=ag.q;S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S);degree=Counter(p for l in W for p in ag.lines[l][1]);tot=Counter();caps=Counter();checks=conditional=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=S&ps;r=len(ls);k=len(pts)
            mult={p:sum(p in ag.lines[l][1] for l in ls) for p in pts};u=sum(b>0 for b in mult.values());h=tuple(sum(degree[p]==j for p in pts) for j in range(1,D+1))
            lines=[(di,lp) for di,lp in ag.lines if lp<=ps];marked={ag.lines[l][0] for l in ls};sec=tuple(sum(len(lp&S)==j for _,lp in lines) for j in range(1,q+1));e=tuple(sorted(tuple(sorted(len(lp&S) for dd,lp in lines if dd==di)) for di in {di for di,_ in lines}-marked)) if r>=5 and r<q+1 else None
            t=(r,k,h,u,sec,e);mu=prices(q,D,t)
            for p,b in mult.items():
                caps[p]=max(caps[p],b)
                for lam in LAMBDAS:tot[p,lam]+=lam*comb(b,2)+int(b>=3)
            for z,(j,lam,sign) in enumerate(objectives(D)):
                actual=sign*sum(lam*comb(b,2)+int(b>=3) for p,b in mult.items() if degree[p]==j)
                assert actual>=mu[z];checks+=1
        projective=json.load(open(f'out/rich14_projective_q{q}.json'))['records']
        for p in S:
            j=degree[p]
            if j<3:continue
            for lam in LAMBDAS:
                lo,hi=global_bounds(q,j);shift=lam*comb(j,2);assert lo+shift<=tot[p,lam]<=hi+shift
                if caps[p]<=4:
                    assert projective[j]['cap_minimum']+shift<=tot[p,lam]<=projective[j]['cap_maximum']+shift;conditional+=1
        rec=dict(file=fn,checks=checks,conditional_cap_checks=conditional,violations=0);records.append(rec);print('POINTMIX14_VALIDATE',rec,flush=True)
    brute=0
    for mult in ((0,1,2,3),(0,0,2,4),(1,1,2,2)):
        profile=tuple(mult.count(b) for b in range(max(mult)+1))
        for ds in ((1,2,2,3),(1,2,3,4),(2,2,3,4)):
            h=tuple(ds.count(j) for j in range(1,5))
            for j in (3,4):
                perms=[p for p in set(permutations(ds)) if all(d>=max(1,b) for d,b in zip(p,mult))]
                expected=tuple((min(vals),max(vals)) for lam in LAMBDAS if (vals:=[sum(lam*comb(b,2)+int(b>=3) for b,d in zip(mult,p) if d==j) for p in perms])) if perms else None
                assert bounds(profile,h,j)==expected;brute+=1
    atomic(f'out/pointmix14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=records,brute_cases=brute))
def table(D,fresh=False):
    path=Path(f'out/pointmix14_D{D}_prices.json');sh=hashlib.sha256(Path('out/sections3_q5_r5.json').read_bytes()).hexdigest();gh=hashlib.sha256(Path('out/rich14_projective_q5.json').read_bytes()).hexdigest()
    if fresh:profiles.cache_clear();bounds.cache_clear();global_bounds.cache_clear()
    types=[freeze(t) for t in json.load(open(f'out/types14_D{D}.json'))['types']]
    result=dict(D=D,q=5,family=FAMILY,section_sha256=sh,projective_sha256=gh,objectives=objectives(D),prices={})
    for number,t in enumerate(types,1):
        result['prices'][key(t)]=prices(5,D,t)
        if number%10000==0:print('POINTMIX14_PRICE',number,len(types),flush=True)
    if fresh and path.exists():
        old=json.load(path.open());assert old['prices']==result['prices'] and old['section_sha256']==sh and old['projective_sha256']==gh
        atomic(f'out/pointmix14_D{D}_fresh_check.json',dict(D=D,types=len(types),prices=len(types)*len(objectives(D)),mismatches=0))
    atomic(path,result);return result
def add_cuts(lp,info,types,D):
    v=json.load(open(f'out/pointmix14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS') and all('conditional_cap_checks' in r for r in v['results'])
    chk=json.load(open(f'out/pointmix14_D{D}_independent.json'));assert chk['status']=='EXACT_POINTMIX14_HUNGARIAN_AGREEMENT'
    cache=json.load(open(f'out/pointmix14_D{D}_prices.json'));assert (cache['D'],cache['q'],cache['family'])==(D,5,FAMILY)
    assert cache['section_sha256']==hashlib.sha256(Path('out/sections3_q5_r5.json').read_bytes()).hexdigest() and cache['projective_sha256']==hashlib.sha256(Path('out/rich14_projective_q5.json').read_bytes()).hexdigest()
    obj=objectives(D);assert cache['objectives']==json.loads(key(obj));pg=json.load(open('out/rich14_projective_q5.json'))['records']
    ordered_prices=[cache['prices'][key(t)] for t in types]
    for z,(j,lam,sign) in enumerate(obj):
        lo=pg[j]['cap_minimum']+lam*comb(j,2);hi=pg[j]['cap_maximum']+lam*comb(j,2)
        row=[(lp.col[f'N{i}'],ordered_prices[i][z]) for i,t in enumerate(types)]+[(lp.col[f'm{j}'],-hi if sign==1 else lo)]
        lp.add(row,'<',0,f'XIV point mixed degree {j} lambda {lam} sign {sign}')
    print('POINTMIX14_CUTS',D,len(obj),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--fresh',action='store_true');a=ap.parse_args();validate(a.D);table(a.D,a.fresh);print('POINTMIX14_COMPLETE',flush=True)
