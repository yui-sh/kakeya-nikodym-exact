"""Exact local high-r plane column catalogue and implicit assignment enumeration."""
from collections import Counter, defaultdict
from itertools import combinations, product
import json, subprocess, hashlib
from pathlib import Path
from geom import AG
import numpy as np
from scipy import sparse as sp
import pencil9 as p9

def freeze(x):
    return tuple(map(freeze,x)) if isinstance(x,list) else x

def signature(ag, lines, S, degree=None):
    r=len(lines); U=set().union(*(ag.lines[l][1] for l in lines))
    t=tuple(sum(len(ps&S)==j for _,ps in ag.lines) for j in range(1,ag.q+1))
    dirs={ag.lines[l][0] for l in lines}
    e=tuple(sorted(tuple(sorted(len(ag.lines[l][1]&S) for l in ag.lines_by_dir[d])) for d in range(ag.q+1) if d not in dirs)) if r>=5 and r<ag.q+1 else None
    h=tuple(Counter(degree[p] for p in S)[j] for j in range(1,max(degree.values())+1)) if degree else None
    return (r,len(S),h,len(U),t,e)

def geometries(D,q=5):
    """Orbit traversal using generators of AGL(2,q), q=5. All extras enumerated.
    One representative per (line orbit,k,t,e); extras have no line incidences.
    """
    path=Path(f'out/geometry14_D{D}.json')
    if path.exists():return json.load(path.open())
    assert q==5
    ag=AG(2,q)
    maps=[lambda x,y:((x+1)%5,y),lambda x,y:(x,(y+1)%5),lambda x,y:(y,x),lambda x,y:((x+y)%5,y),lambda x,y:((2*x)%5,y)]
    perms=[]
    for f in maps:
        pp=[ag.pidx[f(*p)] for p in ag.points]
        perms.append([ag.line_of[frozenset(pp[p] for p in ps)] for _,ps in ag.lines])
    out=[]
    for r in (5,6):
        seen=set(); orbits=[]
        for ds in combinations(range(6),r):
            for ls in product(*(ag.lines_by_dir[d] for d in ds)):
                ls=tuple(sorted(ls))
                if ls in seen:continue
                orbit={ls};stack=[ls]
                while stack:
                    cur=stack.pop()
                    for perm in perms:
                        nxt=tuple(sorted(perm[l] for l in cur))
                        if nxt not in orbit:orbit.add(nxt);stack.append(nxt)
                seen.update(orbit);orbits.append((ls,len(orbit)))
        assert len(seen)==(6 if r==5 else 1)*5**r
        for ls,sz in orbits:
            U=set().union(*(ag.lines[l][1] for l in ls));outside=sorted(set(range(25))-U);sgseen=set()
            for bits in range(1<<len(outside)):
                S=U|{p for i,p in enumerate(outside) if bits>>i&1}
                sig=signature(ag,ls,S)
                if sig in sgseen:continue
                sgseen.add(sig)
                masks={p:sum(1<<i for i,l in enumerate(ls) if p in ag.lines[l][1]) for p in S}
                pts=sorted(S,key=lambda p:(masks[p],p))
                out.append(dict(lines=ls,points=pts,masks=[masks[p] for p in pts],signature=sig,orbit_size=sz))
        print('GEOMETRY',r,'orbits',len(orbits),'covered',len(seen),flush=True)
    json.dump(out,path.open('w'));return out
