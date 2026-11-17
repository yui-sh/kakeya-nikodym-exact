"""Independent geometry checks for new P9 rows and exhaustive oracle tests."""
from collections import Counter, defaultdict
import itertools
import json
import random
import subprocess
from pathlib import Path
from geom import AG
from structure import planes_of
from planecount_check import witness
from pencil9 import column_counts, compatible
import planecount6 as p6


def geometry(fn):
    data=json.load(open(fn));q=data['q'];ag=AG(3,q)
    S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
    degree=Counter(p for l in W for p in ag.lines[l][1]);assert set(degree)==S
    D=max(degree.values());m=tuple(sum(degree[p]==j for p in S) for j in range(1,D+1))
    planes=planes_of(ag)
    labels=[(sum(ag.lines[l][1]<=ps for l in W),len(ps&S),tuple(sum(degree[p]==j for p in ps&S) for j in range(1,D+1))) for _,_,ps in planes]
    allowed={t[:3] for t in p6.hist_types(q,D,len(S))};assert all(key in allowed for key in labels)
    gs={};zs={};columns=[];expected_g=Counter();expected_z=Counter();high=0
    for l in W:
        pts=ag.lines[l][1];dm=tuple(sorted(degree[p] for p in pts));dv=tuple(dm.count(j) for j in range(1,D+1))
        ids=[i for i,(_,_,ps) in enumerate(planes) if pts<=ps];labs=[labels[i] for i in ids]
        assert compatible(m,dv,labs,q)
        high+=sum(key[0]>=5 for key in labs)
        # Test high-r subtraction followed by exact low-r residual completion.
        residual=[m[j]+q*dv[j]-sum(key[2][j] for key in labs if key[0]>=5) for j in range(D)]
        assert residual==[sum(key[2][j] for key in labs if key[0]<5) for j in range(D)]
        gids=[]
        for key in labs:
            pair=(key,dm);g=gs.setdefault(pair,len(gs));gids.append(g)
        z=zs.setdefault((dm,tuple(sorted(key[:2] for key in labs))),len(zs));columns.append((gids,z))
    # Independently count G from plane -> contained witness line incidences.
    for i,(_,_,ps) in enumerate(planes):
        for l in W:
            if ag.lines[l][1]<=ps:
                dm=tuple(sorted(degree[p] for p in ag.lines[l][1]));expected_g[gs[(labels[i],dm)]]+=1
    for gids,z in columns:expected_z[z]+=1
    grow={g:g for g in gs.values()};zrow={z:len(gs)+z for z in zs.values()};actual=Counter()
    for gids,z in columns:actual.update(column_counts(gids,z,grow,zrow))
    for g,v in expected_g.items():actual[grow[g]]-=v
    for z,v in expected_z.items():actual[zrow[z]]-=v
    assert all(v==0 for v in actual.values())
    print(fn,'q',q,'m',m,'pencils',len(columns),'high_r_incidences',high,'new_link_rows',len(actual),'violations',0,flush=True)


def brute_tests():
    rng=random.Random(913);cases=[]
    for case in range(100):
        D=3;groups=[];gid=0
        for ii in range(4):
            group=[]
            for a in range(rng.randrange(1,5)):
                group.append((gid,rng.randrange(-100,101),tuple(rng.randrange(4) for _ in range(D))));gid+=1
            groups.append(group)
        seq=[rng.randrange(4) for _ in range(rng.randrange(1,7))]
        if case%2:
            pick=[rng.choice(groups[ii]) for ii in seq];target=tuple(sum(c[2][j] for c in pick) for j in range(D))
        else:target=tuple(rng.randrange(12) for _ in range(D))
        zcost=rng.randrange(-100,101);best=None
        for choice in itertools.product(*(groups[ii] for ii in seq)):
            if tuple(sum(c[2][j] for c in choice) for j in range(D))==target:
                cost=zcost+sum(c[1] for c in choice);best=cost if best is None else min(best,cost)
        cases.append((groups,seq,target,zcost,best))
    with open('out/pencil9_test.in','w') as f:
        print(3,len(cases),file=f)
        for ci,(groups,seq,target,zcost,best) in enumerate(cases):
            print(len(groups),1,file=f)
            for group in groups:
                print(len(group),file=f)
                for g,c,h in group:print(g,c,*h,file=f)
            print(ci,zcost,len(seq),*seq,*target,file=f)
    with open('out/pencil9_test_price.log','w') as log:subprocess.run(['out/pencil9_price','out/pencil9_test.in','out/pencil9_test.out'],stderr=log,check=True)
    lines=Path('out/pencil9_test.out').read_text().splitlines();assert len(lines)==len(cases)
    for line in lines:
        fs=line.split();i=int(fs[0]);groups,seq,target,zcost,best=cases[i]
        assert (None if fs[1]=='NONE' else int(fs[1]))==best,(i,line,best)
        if best is not None:
            gids=list(map(int,fs[2:]));assert len(gids)==len(seq)
            choice=[next(c for c in groups[ii] if c[0]==g) for ii,g in zip(seq,gids)]
            assert tuple(sum(c[2][j] for c in choice) for j in range(3))==target
            assert zcost+sum(c[1] for c in choice)==best
    print('100 exact oracle vs brute-force cases: zero violations',flush=True)

if __name__=='__main__':
    brute_tests()
    for fn in ['solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json']:geometry(fn)
    print('P9_VALIDATE_OK',flush=True)
