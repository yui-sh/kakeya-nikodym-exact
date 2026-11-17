"""Exact product envelopes for global degree-class pair counts."""
import json
from pathlib import Path
from itertools import combinations
from collections import Counter
from lift11 import atomic
FAMILY='degree-count-product-box-v1'
def rows(D,lower,upper):
    for j,(L,U) in enumerate(zip(lower,upper),1):
        assert 0<=L<=U
        yield [(f'm{j}',1)],'<',U,f'PB14 upper {j}'
        yield [(f'm{j}',-1)],'<',-L,f'PB14 lower {j}'
        yield [(f'M{j},{j}',2),(f'm{j}',-(L+U-1))],'<',-L*U,f'PB14 diagonal chord {j}'
    for i,j in combinations(range(1,D+1),2):
        Li,Ui=lower[i-1],upper[i-1];Lj,Uj=lower[j-1],upper[j-1];M=f'M{i},{j}'
        yield [(M,-1),(f'm{j}',Li),(f'm{i}',Lj)],'<',Li*Lj,f'PB14 product lower LL {i},{j}'
        yield [(M,-1),(f'm{j}',Ui),(f'm{i}',Uj)],'<',Ui*Uj,f'PB14 product lower UU {i},{j}'
        yield [(M,1),(f'm{j}',Ui*-1),(f'm{i}',-Lj)],'<',-Ui*Lj,f'PB14 product upper UL {i},{j}'
        yield [(M,1),(f'm{j}',-Li),(f'm{i}',-Uj)],'<',-Li*Uj,f'PB14 product upper LU {i},{j}'
def add_cuts(lp,D,lower,upper):
    v=json.load(open(f'out/pairbounds14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS')
    count=0
    for coef,sense,rhs,name in rows(D,lower,upper):
        lp.add([(lp.col[n],c) for n,c in coef],sense,rhs,name);count+=1
    print('PAIR_BOUNDS14',D,count,'lower',lower,'upper',upper,flush=True)
def validate(D):
    from geom import AG
    from planecount_check import witness
    results=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S);degree=Counter(p for l in W for p in ag.lines[l][1]);assert max(degree.values())<=D
        m=tuple(sum(degree[p]==j for p in S) for j in range(1,D+1));assert sum(m)==len(S)
        pairs=Counter(tuple(sorted((degree[a],degree[b]))) for a,b in combinations(sorted(S),2))
        val={f'm{j}':m[j-1] for j in range(1,D+1)}
        for i in range(1,D+1):
            for j in range(i,D+1):val[f'M{i},{j}']=pairs[i,j]
        boxes=[]
        for prefix in range(D+1):boxes.append((tuple(m[j] if j<prefix else 0 for j in range(D)),tuple(m[j] if j<prefix else len(S) for j in range(D))))
        boxes.append((tuple(max(0,x-2) for x in m),tuple(min(len(S),x+3) for x in m)))
        checks=0
        for lo,hi in boxes:
            for coef,_,rhs,name in rows(D,lo,hi):
                lhs=sum(c*val[n] for n,c in coef);assert lhs<=rhs,(fn,name,lhs,rhs);checks+=1
        results.append(dict(file=fn,degree_counts=m,boxes=len(boxes),checks=checks,violations=0));print('PAIR_BOUNDS14_VALIDATE',D,results[-1],flush=True)
    # All integer pairs in small boxes, independently of geometric examples.
    brute=0
    for L in range(4):
        for U in range(L,5):
            for l in range(4):
                for u in range(l,5):
                    for x in range(L,U+1):
                        for y in range(l,u+1):
                            val={'m1':x,'m2':y,'M1,1':x*(x-1)//2,'M2,2':y*(y-1)//2,'M1,2':x*y}
                            for coef,_,rhs,_ in rows(2,(L,l),(U,u)):
                                assert sum(c*val[n] for n,c in coef)<=rhs;brute+=1
    atomic(f'out/pairbounds14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=results,brute_checks=brute))
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);a=ap.parse_args();validate(a.D)
