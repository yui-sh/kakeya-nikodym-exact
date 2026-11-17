"""Exact binary-colouring prices for within-line degree-subset pair costs."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import json,hashlib
from pathlib import Path
from functools import lru_cache
from itertools import combinations,permutations
from collections import defaultdict,Counter
from geometry14 import geometries,freeze
from lift11 import atomic,key
CORE=((5,15),(6,17))
FAMILY='binary-subset-within-line-pairs-v1'
def objectives(D):
    return [(s,sign) for size in range(1,D) for s in combinations(range(1,D+1),size) for sign in (1,-1)]
@lru_cache(None)
def extrema(masks,r):
    """All 2^k colourings, Gray traversal, extrema indexed by counts of inside multiplicity."""
    counts=[0]*(r+1);line=[0]*r;cost=0;previous=0;out={}
    incident=[tuple(i for i in range(r) if m>>i&1) for m in masks]
    for number in range(1<<len(masks)):
        bits=number^(number>>1)
        if number:
            changed=bits^previous;p=changed.bit_length()-1;add=bool(bits&changed);b=len(incident[p]);counts[b]+=1 if add else -1
            for i in incident[p]:
                if add:cost+=line[i];line[i]+=1
                else:line[i]-=1;cost-=line[i]
        profile=tuple(counts)
        if profile not in out:out[profile]=[cost,cost,bits,bits]
        else:
            if cost<out[profile][0]:out[profile][0]=cost;out[profile][2]=bits
            if cost>out[profile][1]:out[profile][1]=cost;out[profile][3]=bits
        previous=bits
    # Independently check every recorded extremal witness directly.
    for profile,(lo,hi,lb,hb) in out.items():
        for value,bits in ((lo,lb),(hi,hb)):
            assert tuple(sum(bool(bits>>p&1) and m.bit_count()==b for p,m in enumerate(masks)) for b in range(r+1))==profile
            assert sum((c:=sum(bool(bits>>p&1) and bool(m>>i&1) for p,m in enumerate(masks)))*(c-1)//2 for i in range(r))==value
    return tuple((p,*v) for p,v in sorted(out.items()))
def compatible(profile,h):
    return sum(profile)==sum(h) and all(sum(profile[b:])<=sum(h[b-1:]) for b in range(1,len(profile)))
@lru_cache(None)
def local_price(masks,r,h,subset,sign):
    total=tuple(sum(m.bit_count()==b for m in masks) for b in range(r+1))
    hs=tuple(c if d in subset else 0 for d,c in enumerate(h,1));hc=tuple(c-a for c,a in zip(h,hs));best=None
    for p,lo,hi,_,_ in extrema(masks,r):
        if compatible(p,hs) and compatible(tuple(t-a for t,a in zip(total,p)),hc):
            value=lo if sign==1 else -hi
            best=value if best is None else min(best,value)
    return best
@lru_cache(None)
def matching(D,t):
    wanted=t[:2]+t[3:]
    return tuple(tuple(g['masks']) for g in geometries(D) if (s:=freeze(g['signature']))[:2]+s[3:]==wanted)
def price(D,t,subset,sign):
    vals=[local_price(m,t[0],t[2],tuple(subset),sign) for m in matching(D,t)]
    vals=[v for v in vals if v is not None];assert vals,('missing full type',t,subset)
    return min(vals)
def table(D,fresh=False):
    path=Path(f'out/nonlinear14_D{D}_prices.json');types=[freeze(t) for t in json.load(open(f'out/types14_D{D}.json'))['types'] if tuple(t[:2]) in CORE]
    gh=hashlib.sha256(key(geometries(D)).encode()).hexdigest();obj=objectives(D)
    if fresh:extrema.cache_clear();local_price.cache_clear();matching.cache_clear()
    result=dict(D=D,family=FAMILY,geometry_sha256=gh,objectives=obj,prices={})
    for number,t in enumerate(types,1):
        result['prices'][key(t)]=[price(D,t,s,sign) for s,sign in obj]
        if number%100==0:print('NONLINEAR14_PRICE',D,number,len(types),flush=True)
    if fresh and path.exists():
        old=json.load(path.open());assert old['D']==D and old['family']==FAMILY and old['geometry_sha256']==gh and old['prices']==result['prices']
        atomic(f'out/nonlinear14_D{D}_fresh_check.json',dict(D=D,types=len(types),prices=sum(map(len,result['prices'].values())),mismatches=0,method='fresh exhaustive binary colourings before loading cached prices'))
    atomic(path,result);return result

def validate(D):
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    results=[];obj=objectives(D)
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S);degree=Counter(p for l in W for p in ag.lines[l][1]);assert max(degree.values())<=D
        checks=full=planes=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=sorted(S&ps);r=len(ls);k=len(pts)
            if ag.q==5 and (r,k) not in CORE:continue
            masks=tuple(sum(1<<i for i,l in enumerate(ls) if p in ag.lines[l][1]) for p in pts);h=tuple(sum(degree[p]==d for p in pts) for d in range(1,D+1));planes+=1;t=None
            if ag.q==5:
                local=[(d,lp) for d,lp in ag.lines if lp<=ps];marked={ag.lines[l][0] for l in ls};sec=tuple(sum(len(lp&S)==j for _,lp in local) for j in range(1,6))
                e=tuple(sorted(tuple(sorted(len(lp&S) for dd,lp in local if dd==d)) for d in {d for d,_ in local}-marked)) if r==5 else None
                union=set().union(*(ag.lines[l][1] for l in ls));t=(r,k,h,len(union),sec,e)
            for subset,sign in obj:
                actual=sign*sum((c:=sum(degree[p] in subset for p in ag.lines[l][1]))*(c-1)//2 for l in ls)
                mu=local_price(masks,r,h,subset,sign);assert mu is not None and actual>=mu
                if t is not None:assert actual>=price(D,t,subset,sign);full+=1
                checks+=1
        rec=dict(file=fn,planes=planes,checks=checks,full_type_checks=full,violations=0);results.append(rec);print('NONLINEAR14_VALIDATE',D,rec,flush=True)
    brute=0
    for masks in ((1,1,2,3),(1,2,3,3),(0,1,2,3)):
        for ds in ((1,2,2,3),(1,1,2,3),(2,2,3,3)):
            h=tuple(ds.count(d) for d in range(1,4))
            for subset,sign in objectives(3):
                vals=[]
                for perm in set(permutations(ds)):
                    if not all(d>=max(1,m.bit_count()) for d,m in zip(perm,masks)):continue
                    vals.append(sign*sum((c:=sum(d in subset and bool(m>>i&1) for d,m in zip(perm,masks)))*(c-1)//2 for i in range(2)))
                assert local_price(masks,2,h,subset,sign)==(min(vals) if vals else None);brute+=1
    atomic(f'out/nonlinear14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=results,brute_cases=brute))

def add_cuts(lp,info,types,D,scope='all',compact=False):
    assert scope in ('all','pair','singletons')
    v=json.load(open(f'out/nonlinear14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS')
    cache=json.load(open(f'out/nonlinear14_D{D}_prices.json'));assert cache['D']==D and cache['family']==FAMILY and cache['geometry_sha256']==hashlib.sha256(key(geometries(D)).encode()).hexdigest()
    obj=objectives(D);assert cache['objectives']==json.loads(key(obj));groups=defaultdict(list);gv=defaultdict(list)
    for i,t in enumerate(types):
        if t[:2] in CORE:groups[t[:3]].append(i)
    for g,(a,di) in enumerate(info['_Y']):
        if a[:2] in CORE:gv[a].append((lp.col[f'G{g}'],info['_dms'][di]))
    count=0
    for a,indices in groups.items():
        seen={};support={j for j,h in enumerate(a[2],1) if h}
        if compact:assert all(set(dm)<=support for _,dm in gv[a])
        for z,(subset,sign) in enumerate(obj):
            if scope=='pair' and len(subset)>2:continue
            if scope=='singletons' and len(subset)>1:continue
            if compact:
                effective=tuple(d for d in subset if d in support)
                if not effective:
                    assert all(cache['prices'][key(types[i])][z]==0 for i in indices)
                    continue
                signature=(effective,sign)
                if signature in seen:
                    old=seen[signature]
                    assert all(cache['prices'][key(types[i])][z]==cache['prices'][key(types[i])][old] for i in indices)
                    continue
                seen[signature]=z
            row=[(lp.col[f'N{i}'],cache['prices'][key(types[i])][z]) for i in indices]
            row += [(col,-sign*(c:=sum(d in subset for d in dm))*(c-1)//2) for col,dm in gv[a]]
            lp.add(row,'<',0,f'XIV nonlinear {key(a)} {key(subset)} {sign}');count+=1
    print('NONLINEAR14_CUTS',D,count,flush=True)
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--fresh',action='store_true');a=ap.parse_args()
    validate(a.D);r=table(a.D,a.fresh);print('NONLINEAR14_COMPLETE',a.D,len(r['prices']),flush=True)
