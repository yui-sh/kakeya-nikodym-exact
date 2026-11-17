"""Geometry validation and independent brute-force checks, before LP use."""
from collections import Counter
from itertools import combinations
import json, random, subprocess
from pathlib import Path
from geom import AG
from structure import planes_of
from planecount_check import witness
from plane10 import geometries,signature,freeze

def geometry(fn):
    data=json.load(open(fn));q=data['q'];ag=AG(3,q);a2=AG(2,q)
    S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
    deg=Counter(p for l in W for p in ag.lines[l][1]);D=max(deg.values())
    expected=Counter();actual=Counter();num=0;cols=Counter()
    for _,_,ps in planes_of(ag):
        ls=[l for l in W if ag.lines[l][1]<=ps]
        if len(ls)<5:continue
        num+=1;o=ag.points[min(ps)];vecs=[ag.sub(ag.points[p],o) for p in sorted(ps) if ag.points[p]!=o]
        u=vecs[0];span={ag.scale(t,u) for t in range(q)};v=next(v for v in vecs if v not in span)
        lift=[ag.pidx[ag.add(o,ag.add(ag.scale(x,u),ag.scale(y,v)))] for x,y in a2.points];inv={p:i for i,p in enumerate(lift)}
        localS={inv[p] for p in ps&S};localL=[a2.line_of[frozenset(inv[p] for p in ag.lines[l][1])] for l in ls]
        ds={inv[p]:deg[p] for p in ps&S};sig=signature(a2,localL,localS,ds);h=sig[2]+(0,)*(D-len(sig[2]));sig=sig[:2]+(h,)+sig[3:]
        U=set().union(*(a2.lines[l][1] for l in localL))
        assert all(ds[p]>=sum(p in a2.lines[l][1] for l in localL) for p in localS)
        dms=tuple(sorted(tuple(sorted(ds[p] for p in a2.lines[l][1])) for l in localL))
        cols[(sig,dms)]+=1
        for dm in dms:actual[(sig[:3],dm)]+=1
        if q==5:
            # Independently find an affine image of a catalogue line representative.
            candidates=[g for g in geometries() if freeze(g['signature'])[:2]+freeze(g['signature'])[3:]==sig[:2]+sig[3:]]
            assert candidates
            wanted={tuple(g['lines']) for g in candidates};target=tuple(sorted(localL))
            maps=[lambda x,y:((x+1)%5,y),lambda x,y:(x,(y+1)%5),lambda x,y:(y,x),lambda x,y:((x+y)%5,y),lambda x,y:((2*x)%5,y)]
            perms=[]
            for f in maps:
                pp=[a2.pidx[f(*p)] for p in a2.points];perms.append([a2.line_of[frozenset(pp[p] for p in lp)] for _,lp in a2.lines])
            seen={target};stack=[target];found=False
            while stack:
                cur=stack.pop()
                if cur in wanted:found=True;break
                for perm in perms:
                    nxt=tuple(sorted(perm[l] for l in cur))
                    if nxt not in seen:seen.add(nxt);stack.append(nxt)
            assert found
    # independent line->plane incidence count (opposite traversal)
    for l in W:
        for _,_,ps in planes_of(ag):
            ls=[ll for ll in W if ag.lines[ll][1]<=ps]
            if len(ls)>=5 and ag.lines[l][1]<=ps:
                key=(len(ls),len(ps&S),tuple(sum(deg[p]==j for p in ps&S) for j in range(1,D+1)))
                expected[(key,tuple(sorted(deg[p] for p in ag.lines[l][1])))]+=1
    assert actual==expected and sum(cols.values())==num
    print(fn,'high_planes',num,'incidences',sum(actual.values()),'violations 0',flush=True)

def brute():
    ag=AG(2,5);g=next(g for g in geometries() if g['signature'][0:2]==[5,15]);pts=g['points'];masks=g['masks']
    assignments=[]
    for four in combinations(range(15),3):
        ds=[4 if p in four else 3 for p in range(15)]
        dms=[tuple(sorted(ds[p] for p in range(15) if masks[p]>>j&1)) for j in range(5)]
        assignments.append((ds,dms))
    all_dm=sorted({dm for _,dms in assignments for dm in dms});rng=random.Random(1010)
    for case in range(24):
        costs={dm:rng.randrange(-100,101) for dm in all_dm if rng.randrange(5) or case==0}
        if case==20:costs={}
        if case>=21:costs={dm:0 for dm in all_dm}
        valid=[(sum(costs[dm] for dm in dms),ds,dms) for ds,dms in assignments if all(dm in costs for dm in dms)]
        with open('out/plane10_test.in','w') as f:
            print('1 0\n0 0 5 15 5 0\n0 0 12 3 0',file=f);print(*masks,file=f);print(len(costs),file=f)
            for i,(dm,c) in enumerate(costs.items()):print(i,c,*dm,file=f)
        subprocess.run(['out/plane10_enum','out/plane10_test.in','out/plane10_test.out'],check=True,stderr=subprocess.DEVNULL)
        out=Path('out/plane10_test.out').read_text().split()
        assert bool(out)==bool(valid)
        if out:assert int(out[2])==min(c for c,_,_ in valid)
    print('24 exhaustive pricing vs independent brute-force tests: zero discrepancies',flush=True)

if __name__=='__main__':
    geometries();brute()
    for fn in ['solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json']:geometry(fn)
    print('P10_VALIDATE_OK',flush=True)
