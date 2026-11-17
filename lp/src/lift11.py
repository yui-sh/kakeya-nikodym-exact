"""M-free incidence cuts: exact assignment pricing, canonical caches, geometry checks."""
import os
for _v in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[_v]="1"
import json, hashlib
from pathlib import Path
from collections import defaultdict, Counter
from functools import lru_cache
from plane10 import geometries, freeze
INF = 10**6

def key(x): return json.dumps(x, separators=(',', ':'))
def atomic(path, obj):
    p=Path(path); tmp=p.with_suffix(p.suffix+'.tmp.'+str(os.getpid())); tmp.write_text(json.dumps(obj,indent=1)); tmp.replace(p)

def assignment(a):
    """Integer Hungarian algorithm; dual feasibility and equality independently checked."""
    n=len(a); u=[0]*(n+1);v=[0]*(n+1);p=[0]*(n+1);way=[0]*(n+1)
    for i in range(1,n+1):
        p[0]=i;j0=0;mv=[INF*100]*(n+1);used=[False]*(n+1)
        while True:
            used[j0]=True;i0=p[j0];delta=INF*100;j1=0
            for j in range(1,n+1):
                if not used[j]:
                    cur=a[i0-1][j-1]-u[i0]-v[j]
                    if cur<mv[j]:mv[j]=cur;way[j]=j0
                    if mv[j]<delta:delta=mv[j];j1=j
            for j in range(n+1):
                if used[j]:u[p[j]]+=delta;v[j]-=delta
                else:mv[j]-=delta
            j0=j1
            if p[j0]==0:break
        while True:
            j1=way[j0];p[j0]=p[j1];j0=j1
            if j0==0:break
    perm=[0]*n
    for j in range(1,n+1):perm[p[j]-1]=j-1
    cost=sum(a[i][perm[i]] for i in range(n))
    assert cost==sum(u[1:])+sum(v[1:])
    assert all(u[i+1]+v[j+1]<=a[i][j] for i in range(n) for j in range(n))
    return cost,perm,u[1:],v[1:]

@lru_cache(None)
def price_profile(mult,h,j,sgn):
    ds=[d for d,c in enumerate(h,1) for _ in range(c)]
    a=[[sgn*b*(d==j) if d>=max(1,b) else INF for d in ds] for b in mult]
    cost,perm,u,v=assignment(a)
    if any(a[i][perm[i]]==INF for i in range(len(a))):return None
    return cost

def prices(t,geos):
    r,k,h,u,sec,e=t;wanted=(r,k,u,sec,e)
    profiles=set()
    for g in geos:
        s=freeze(g['signature'])
        if s[:2]+s[3:]==wanted:profiles.add(tuple(sorted(m.bit_count() for m in g['masks'])))
    vals=[]
    for j in range(1,len(h)+1):
        for sign in (1,-1):
            pp=[price_profile(p,h,j,sign) for p in profiles]
            pp=[v for v in pp if v is not None]
            vals.append(min(pp) if pp else None)
    assert vals and all(v is not None for v in vals), ('missing full type',t)
    return vals

def add_cuts(lp,info,types,scope="core"):
    assert scope in ("core","all")
    geos=geometries();groups=defaultdict(list);gvars=defaultdict(list)
    path=Path('out/lift11_prices.json')
    gh=hashlib.sha256(key(geos).encode()).hexdigest()
    cache=json.loads(path.read_text()) if path.exists() else dict(geometry_sha256=gh,prices={})
    assert cache['geometry_sha256']==gh
    for i,t in enumerate(types):
        if t[0]>=5:
            groups[t[:3]].append(i)
            kk=key(t)
            if kk not in cache['prices']:cache['prices'][kk]=prices(t,geos)
    atomic(path,cache)
    for g,(a,di) in enumerate(info['_Y']):
        if a[0]>=5:gvars[a].append((lp.col[f'G{g}'],info['_dms'][di]))
    count=0
    for a,inds in groups.items():
        if scope=='core' and a[:2] not in ((5,15),(6,17)):continue
        count+=10
        for j in range(1,6):
            for sign in (1,-1):
                z=2*(j-1)+(sign==-1)
                row=[(lp.col[f'N{i}'],cache['prices'][key(types[i])][z]) for i in inds]
                row += [(col,-sign*dm.count(j)) for col,dm in gvars[a]]
                lp.add(row,'<',0,f'XI incidence {key(a)} {j} {sign}')
    print('LIFT11',count,'cuts',sum(len(v) for v in groups.values()),'full types',len(cache['prices']),'cached',flush=True)
    return cache

def validate():
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    results=[]
    for fn in ['solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json']:
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
        deg=Counter(p for l in W for p in ag.lines[l][1]);D=max(deg.values());checks=0;high=0
        # All planes exercise the assignment family, including q=3 where r>=5 is vacuous.
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];pts=sorted(S&ps)
            mult=tuple(sorted(sum(p in ag.lines[l][1] for l in ls) for p in pts))
            h=tuple(sum(deg[p]==j for p in pts) for j in range(1,D+1))
            high+=len(ls)>=5
            full_mu=None
            if ag.q==5 and len(ls)>=5:
                local=[(d,lp) for d,lp in ag.lines if lp<=ps]
                marked={ag.lines[l][0] for l in ls}
                sec=tuple(sum(len(lp&S)==j for d,lp in local) for j in range(1,6))
                e=tuple(sorted(tuple(sorted(len(lp&S) for dd,lp in local if dd==d)) for d in {d for d,lp in local}-marked)) if len(ls)==5 else None
                union=set().union(*(ag.lines[l][1] for l in ls))
                full_mu=prices((len(ls),len(pts),h,len(union),sec,e),geometries())
            for j in range(1,D+1):
                incidence=sum(sum(deg[p]==j for p in ag.lines[l][1]) for l in ls)
                for sign in (1,-1):
                    mu=price_profile(mult,h,j,sign)
                    assert mu is not None and sign*incidence>=mu
                    if full_mu is not None:assert sign*incidence>=full_mu[2*(j-1)+(sign==-1)]
                    checks+=1
        rec=dict(file=fn,checks=checks,high_planes=high,violations=0);results.append(rec);print(rec,flush=True)
    # Independent exhaustive permutation regression, signed costs / forbidden assignments.
    from itertools import permutations
    for mult in [(0,1,1,2),(1,2,2,3),(0,0,1,3)]:
        for ds in [(1,2,2,3),(1,1,2,3),(2,2,3,3)]:
            h=tuple(ds.count(j) for j in range(1,4))
            for j in range(1,4):
                for sign in (1,-1):
                    costs=[sum(sign*b*(d==j) for b,d in zip(mult,p)) for p in set(permutations(ds)) if all(d>=max(1,b) for b,d in zip(mult,p))]
                    assert price_profile(mult,h,j,sign)==(min(costs) if costs else None)
    atomic('out/lift11_validation.json',dict(results=results,brute_cases=54,status='ZERO_VIOLATIONS'))
    print('LIFT11_VALIDATE_OK',flush=True)
if __name__=='__main__':validate()
