"""Exact convex-hull lift of degree counts and their pair moments."""
import json
from itertools import combinations
from collections import Counter
from lift11 import atomic
FAMILY='raw-degree-vector-pair-moment-hull-v1'
def add_cuts(lp,D,vectors,check=True):
    if check:
        v=json.load(open(f'out/degreehull14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS')
    vectors=tuple(tuple(m) for m in vectors);assert vectors and len(vectors)==len(set(vectors)) and all(len(m)==D and all(type(x)==int and x>=0 for x in m) for m in vectors)
    assert len({sum(m) for m in vectors})==1 and len({sum(j*x for j,x in enumerate(m,1)) for m in vectors})==1
    Q=[lp.var('Q14:'+','.join(map(str,m))) for m in vectors]
    lp.add([(q,1) for q in Q],'=',1,'XIV degree hull normalization')
    for j in range(1,D+1):
        lp.add([(lp.col[f'm{j}'],1)]+[(q,-m[j-1]) for q,m in zip(Q,vectors) if m[j-1]],'=',0,f'XIV degree hull m{j}')
    for i in range(1,D+1):
        for j in range(i,D+1):
            def moment(m):return m[i-1]*(m[i-1]-1)//2 if i==j else m[i-1]*m[j-1]
            lp.add([(lp.col[f'M{i},{j}'],1)]+[(q,-moment(m)) for q,m in zip(Q,vectors) if moment(m)],'=',0,f'XIV degree hull M{i},{j}')
    print('DEGREE_HULL14',D,'columns',len(Q),'rows',1+D+D*(D+1)//2,flush=True)
    return Q

def validate(D):
    from geom import AG
    from planecount_check import witness
    from planecount5 import LP
    results=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S);degree=Counter(p for l in W for p in ag.lines[l][1]);assert max(degree.values())<=D
        m=tuple(sum(degree[p]==j for p in S) for j in range(1,D+1));assert sum(m)==len(S)
        pairs=Counter(tuple(sorted((degree[a],degree[b]))) for a,b in combinations(sorted(S),2))
        lp=LP()
        for j in range(1,D+1):lp.var(f'm{j}')
        for i in range(1,D+1):
            for j in range(i,D+1):lp.var(f'M{i},{j}')
        candidates=[m]
        for j in range(1,D-1):
            if m[j]>=2:
                v=list(m);v[j]-=2;v[j-1]+=1;v[j+1]+=1;candidates.append(tuple(v))
        Q=add_cuts(lp,D,candidates,False);values=[0]*len(lp.names);values[Q[0]]=1
        for j in range(1,D+1):values[lp.col[f'm{j}']]=m[j-1]
        for i in range(1,D+1):
            for j in range(i,D+1):values[lp.col[f'M{i},{j}']]=pairs[i,j]
        for row,rhs,sense,name in lp.rows:assert sense=='=' and sum(c*values[k] for k,c in row.items())==rhs,(fn,name)
        results.append(dict(file=fn,degree_counts=m,columns=len(Q),checks=len(lp.rows),violations=0));print('DEGREE_HULL14_VALIDATE',D,results[-1],flush=True)
    atomic(f'out/degreehull14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=results,method='delta column for actual degree vector, with pair counts independently counted over geometric point pairs'))
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);a=ap.parse_args();validate(a.D)
