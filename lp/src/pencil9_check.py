"""Solver-free replay of saved pencil9 primals and all completion witnesses.
Also attempts exact rational reconstruction; failure never certifies anything.
"""
import argparse
from collections import Counter
from fractions import Fraction
import json
import numpy as np
from pencil9 import build, compatible, column_counts


def check(tag, rational_only=False):
    result=json.load(open(f'out/{tag}.json'));m=tuple(result['m']);config=tuple(result['config'])
    lp,info,types,grow,zrow=build(config,m)
    columns=json.load(open(f'out/{tag}_columns.json'));x=np.load(f'out/{tag}_primal.npz')['x']
    assert len(x)==len(lp.names)+len(columns)
    assert len(set((z,tuple(gids)) for z,gids in columns))==len(columns)
    sparse=[]
    for z,gids in columns:
        di,wi=info['_Z'][z];dm=info['_dms'][di];dv=tuple(dm.count(j) for j in range(1,len(m)+1))
        keys=[info['_Y'][g][0] for g in gids]
        assert all(info['_Y'][g][1]==di for g in gids)
        assert compatible(m,dv,keys,5)
        items=sorted(next(ii for ii,it in enumerate(info['_items']) if key[0]==it[0] and key[1] in it[3]) for key in keys)
        assert items==list(info['_wt'][wi])
        sparse.append(column_counts(gids,z,grow,zrow))
    if rational_only:
        saved=json.load(open(f'out/{tag}_rational_primal.json'))
        rational=[Fraction(0)]*len(x)
        for i,v in saved.items():rational[int(i)]=Fraction(v)
        assert min(rational)>=0
        left=[sum((coef*rational[v] for v,coef in row[0].items()),Fraction(0)) for row in lp.rows]
        for a,col in zip(rational[len(lp.names):],sparse):
            if a:
                for r,c in col.items():left[r]+=c*a
        for i,row in enumerate(lp.rows):
            assert (left[i]==row[1] if row[2]=='=' else left[i]<=row[1]),(i,row[3],left[i],row[1])
        print(tag,'SAVED_RATIONAL_VERIFIED',len(lp.rows),'rows',len(columns),'completion columns',sum(v>0 for v in rational),'positive variables',flush=True)
        return
    # Recompute the scalar rows directly from dictionaries, not solver matrices.
    lhs=np.array([sum(coef*x[v] for v,coef in row[0].items()) for row in lp.rows])
    for a,col in zip(x[len(lp.names):],sparse):
        for r,c in col.items():lhs[r]+=c*a
    eq=max(abs(lhs[i]-row[1]) for i,row in enumerate(lp.rows) if row[2]=='=')
    ub=max(0,max(lhs[i]-row[1] for i,row in enumerate(lp.rows) if row[2]=='<'))
    assert eq<1e-7 and ub<1e-7 and min(x)>-1e-7
    report=dict(tag=tag,columns_checked=len(columns),positive_columns=int(sum(x[len(lp.names):]>1e-8)),max_eq_residual=eq,max_ineq_violation=ub,min_variable=float(min(x)))
    for den in [10000,1000000,100000000]:
        rational=[Fraction(0) if abs(v)<1e-8 else Fraction(float(v)).limit_denominator(den) for v in x]
        if min(rational)<0:continue
        left=[sum((coef*rational[v] for v,coef in row[0].items()),Fraction(0)) for row in lp.rows]
        for a,col in zip(rational[len(lp.names):],sparse):
            if not a:continue
            for r,c in col.items():left[r]+=c*a
        bad=[i for i,row in enumerate(lp.rows) if (left[i]!=row[1] if row[2]=='=' else left[i]>row[1])]
        print(tag,'denominator limit',den,'violated exact rows',len(bad),flush=True)
        if not bad:
            report['exact_rational_feasible']=True;report['max_denominator']=max(v.denominator for v in rational)
            json.dump({str(i):str(v) for i,v in enumerate(rational) if v},open(f'out/{tag}_rational_primal.json','w'),indent=1)
            break
    else:
        # Reconstruct on numerical support by exact sparse Gaussian elimination.
        support=[i for i,v in enumerate(x) if v>1e-8]; support_set=set(support)
        equations=[{i:Fraction(c) for i,c in row[0].items() if i in support_set} for row in lp.rows]
        for j,col in enumerate(sparse,len(lp.names)):
            if j in support_set:
                for r,c in col.items():equations[r][j]=Fraction(c)
        basis={}
        for r,row in enumerate(lp.rows):
            if row[2]!='=' and abs(lhs[r]-row[1])>1e-8:continue
            d=equations[r].copy();rhs=Fraction(row[1])
            while d:
                pivot=min(d)
                if pivot not in basis:
                    scale=d[pivot];basis[pivot]=({i:c/scale for i,c in d.items()},rhs/scale);break
                factor=d[pivot];bd,br=basis[pivot];rhs-=factor*br
                for i,c in bd.items():
                    d[i]=d.get(i,Fraction(0))-factor*c
                    if not d[i]:del d[i]
            else:assert rhs==0, ('support inconsistent',r,rhs)
        rational=[Fraction(0)]*len(x)
        for i in support:
            if i not in basis:rational[i]=Fraction(float(x[i])).limit_denominator(1000000)
        for pivot in sorted(basis,reverse=True):
            d,rhs=basis[pivot];rational[pivot]=rhs-sum(c*rational[i] for i,c in d.items() if i!=pivot)
        left=[sum((coef*rational[v] for v,coef in row[0].items()),Fraction(0)) for row in lp.rows]
        for a,col in zip(rational[len(lp.names):],sparse):
            if a:
                for r,c in col.items():left[r]+=c*a
        bad=[i for i,row in enumerate(lp.rows) if (left[i]!=row[1] if row[2]=='=' else left[i]>row[1])]
        exact=not bad and min(rational)>=0
        report.update(exact_rational_feasible=exact,support=len(support),rank=len(basis))
        if exact:
            report['max_denominator']=max(v.denominator for v in rational)
            json.dump({str(i):str(v) for i,v in enumerate(rational) if v},open(f'out/{tag}_rational_primal.json','w'),indent=1)
        print(tag,'exact reconstruction: rank',len(basis),'support',len(support),'bad rows',len(bad),'nonnegative',min(rational)>=0,flush=True)
    json.dump(report,open(f'out/{tag}_check.json','w'),indent=2)
    print(json.dumps(report),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('tags',nargs='+');ap.add_argument('--rational-only',action='store_true');a=ap.parse_args()
    for tag in a.tags:check(tag,a.rational_only)
    print('P9_PRIMALS_CHECKED',flush=True)
