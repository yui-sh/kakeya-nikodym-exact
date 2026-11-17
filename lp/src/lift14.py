"""D-parameterized universal additive subset cuts; isolated XIV caches."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import json,hashlib
from pathlib import Path
from itertools import combinations,permutations
from functools import lru_cache
from collections import defaultdict,Counter
from geometry14 import geometries,freeze,signature
from lift11 import assignment,INF,atomic,key
FAMILY='all-subsets-core-plus-all-coordinate-v1'
CORE=((5,15),(6,17))

def subsets(D):
    out=[(i,) for i in range(1,D+1)]
    for k in range(2,D//2+1):
        for s in combinations(range(1,D+1),k):
            complement=tuple(i for i in range(1,D+1) if i not in s)
            if 2*k==D and s>complement:continue
            out.append(s)
    return out

def objectives(D):return [(s,sign) for s in subsets(D) for sign in (1,-1)]
@lru_cache(None)
def profile_price(mult,h,subset,sign):
    ds=[d for d,c in enumerate(h,1) for _ in range(c)]
    matrix=[[sign*b*(d in subset) if d>=max(1,b) else INF for d in ds] for b in mult]
    value,perm,_,_=assignment(matrix)
    if any(matrix[i][perm[i]]==INF for i in range(len(matrix))):return None
    return value

@lru_cache(None)
def profile_index(D):
    idx=defaultdict(set)
    for g in geometries(D):
        s=freeze(g['signature']);idx[s[:2]+s[3:]].add(tuple(sorted(m.bit_count() for m in g['masks'])))
    return idx

def price(t,D,objs):
    r,k,h,u,sec,e=t;profiles=profile_index(D)[(r,k,u,sec,e)];values=[]
    for subset,sign in objs:
        vals=[profile_price(p,h,tuple(subset),sign) for p in profiles];vals=[v for v in vals if v is not None]
        assert vals,('missing type assignments',t,subset)
        values.append(min(vals))
    return values

def table(D,types,fresh=False):
    path=Path(f'out/lift14_D{D}_prices.json');gh=hashlib.sha256(key(geometries(D)).encode()).hexdigest();obj=objectives(D)
    if fresh:profile_price.cache_clear();profile_index.cache_clear()
    partial=path.with_name(path.stem+'_partial.json')
    source=path if path.exists() else partial
    cache=None if fresh or not source.exists() else json.load(source.open())
    if cache:assert cache['D']==D and cache['family']==FAMILY and cache['geometry_sha256']==gh and cache['objectives']==json.loads(key(obj))
    cache=cache or dict(D=D,family=FAMILY,geometry_sha256=gh,objectives=obj,prices={})
    count=0
    for t in types:
        if t[0]<5:continue
        kk=key(t)
        if kk not in cache['prices']:
            selected=obj if t[:2] in CORE else obj[:2*D]
            cache['prices'][kk]=price(t,D,selected);count+=1
            if count%10000==0:
                atomic(path.with_name(path.stem+('_fresh_partial' if fresh else '_partial')+'.json'),cache)
                print('PRICE14_PROGRESS',D,count,flush=True)
    if fresh and path.exists():
        old=json.load(path.open());assert old['D']==D and old['geometry_sha256']==gh and old['family']==FAMILY
        assert old['prices']==cache['prices']
        atomic(f'out/lift14_D{D}_fresh_check.json',dict(D=D,types=len(cache['prices']),mismatches=0,inputs='new geometry-based pricing before reading old table'))
        print('FRESH_PRICES14_VERIFIED',D,len(cache['prices']),'mismatches 0',flush=True)
    atomic(path,cache);return cache

def add_cuts(lp,info,types,D,scope='all',compact=False):
    assert scope in ('all','core','pair','coordinate')
    validation=json.load(open(f'out/validation14_D{D}.json'));assert validation['status']=='ZERO_VIOLATIONS' and validation['D']==D and validation['family']==FAMILY
    cache=json.load(open(f'out/lift14_D{D}_prices.json'));assert cache['D']==D and cache['family']==FAMILY
    assert cache['geometry_sha256']==hashlib.sha256(key(geometries(D)).encode()).hexdigest()
    obj=objectives(D);assert cache['objectives']==json.loads(key(obj))
    groups=defaultdict(list);gv=defaultdict(list)
    for i,t in enumerate(types):
        if t[0]>=5:groups[t[:3]].append(i)
    for g,(a,di) in enumerate(info['_Y']):
        if a[0]>=5:gv[a].append((lp.col[f'G{g}'],info['_dms'][di]))
    count=0
    for a,indices in groups.items():
        seen={};support={j for j,h in enumerate(a[2],1) if h}
        if compact:assert all(set(dm)<=support for _,dm in gv[a])
        if scope!='all' and a[:2] not in CORE:continue
        selected=obj if a[:2] in CORE else obj[:2*D]
        for z,(subset,sign) in enumerate(selected):
            if scope=='pair' and len(subset)>2:continue
            if scope=='coordinate' and len(subset)>1:continue
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
            row += [(col,-sign*sum(d in subset for d in dm)) for col,dm in gv[a]]
            lp.add(row,'<',0,f'XIV subset {key(a)} {key(subset)} {sign}');count+=1
    print('CUTS14',D,count,'types',sum(map(len,groups.values())),flush=True)

def validate(D):
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    ag2=AG(2,5);geos=geometries(D)
    for g in geos:
        ls=g['lines'];pts=g['points'];assert signature(ag2,ls,set(pts))==freeze(g['signature'])
        assert g['masks']==[sum(1<<i for i,l in enumerate(ls) if p in ag2.lines[l][1]) for p in pts]
    results=[]
    for fn in ['solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json']:
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
        degree=Counter(p for l in W for p in ag.lines[l][1]);actualD=max(degree.values());slotD=max(D,actualD)
        # Exercise coordinate, pair and all higher campaign subset sizes, using true example degrees.
        objs=[(ss,sgn) for k in range(1,min(max(2,D//2),slotD)+1) for ss in combinations(range(1,slotD+1),k) for sgn in (1,-1)]
        checks=full=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=sorted(S&ps)
            mult=tuple(sorted(sum(p in ag.lines[l][1] for l in ls) for p in pts));h=tuple(sum(degree[p]==d for p in pts) for d in range(1,slotD+1));mu=None
            if ag.q==5 and len(ls)>=5:
                local=[(d,lp) for d,lp in ag.lines if lp<=ps];marked={ag.lines[l][0] for l in ls}
                sec=tuple(sum(len(lp&S)==j for _,lp in local) for j in range(1,6))
                e=tuple(sorted(tuple(sorted(len(lp&S) for dd,lp in local if dd==d)) for d in {d for d,_ in local}-marked)) if len(ls)==5 else None
                union=set().union(*(ag.lines[l][1] for l in ls));mu=price((len(ls),len(pts),h,len(union),sec,e),D,objs)
            for z,(subset,sign) in enumerate(objs):
                actual=sign*sum(sum(degree[p] in subset for p in ag.lines[l][1]) for l in ls)
                p=profile_price(mult,h,subset,sign);assert p is not None and actual>=p
                if mu is not None:assert actual>=mu[z];full+=1
                checks+=1
        rec=dict(file=fn,actual_D=actualD,slot_D=slotD,checks=checks,full_type_checks=full,violations=0);results.append(rec);print('VALIDATION14',D,rec,flush=True)
    brute=0
    for mult in ((0,1,1,2),(1,2,2,3),(0,0,1,3),(1,1,2,4)):
        for ds in ((1,2,2,3),(1,1,2,3),(2,2,3,3),(1,2,3,4)):
            h=tuple(ds.count(d) for d in range(1,5))
            for k in (1,2,3):
                for subset in combinations(range(1,5),k):
                    for sign in (1,-1):
                        vals=[sum(sign*b*(d in subset) for b,d in zip(mult,p)) for p in set(permutations(ds)) if all(d>=max(1,b) for b,d in zip(mult,p))]
                        assert profile_price(mult,h,subset,sign)==(min(vals) if vals else None);brute+=1
    atomic(f'out/validation14_D{D}.json',dict(D=D,status='ZERO_VIOLATIONS',family=FAMILY,geometries=len(geos),results=results,brute_cases=brute))

if __name__=='__main__':
    import argparse,planecount6 as p6
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--fresh',action='store_true');a=ap.parse_args()
    validate(a.D);types=p6.hist_types(5,a.D,46)
    atomic(f'out/types14_D{a.D}.json',dict(D=a.D,q=5,n=46,types=types))
    print('TYPES14',a.D,len(types),'high_r',sum(t[0]>=5 for t in types),flush=True)
    result=table(a.D,types,a.fresh);print('PREPARE14_COMPLETE',a.D,len(result['prices']),flush=True)
