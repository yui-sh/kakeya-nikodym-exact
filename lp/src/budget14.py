"""Exact plane-cost caps from the inherited quadratic counting certificate."""
import json,hashlib
from pathlib import Path
from functools import lru_cache
from itertools import product
from math import comb
from geom import AG
from lift11 import atomic
FAMILY='nonnegative-quadratic-plane-budget-v1'

def polynomial(r,k):return 10*k*k-(150+17*r)*k+560+143*r

@lru_cache(None)
def minima(q):
    ag=AG(2,q);f=[q*q+1]*(q+2);counts=[0]*(q+2);witnesses=[None]*(q+2)
    masks=[sum(1<<p for p in pts) for _,pts in ag.lines]
    choices=[[-1]+list(ag.lines_by_dir[d]) for d in range(q+1)]
    for selected in product(*choices):
        ls=[l for l in selected if l>=0];mask=0
        for l in ls:mask|=masks[l]
        r=len(ls);u=mask.bit_count();counts[r]+=1
        if u<f[r]:f[r]=u;witnesses[r]=ls
    assert counts==[comb(q+1,r)*q**r for r in range(q+2)]
    return dict(q=q,minima=f,configurations=counts,witnesses=witnesses)

def shift(q):
    f=minima(q)['minima']
    return -min(0,min(polynomial(r,k) for r in range(q+2) for k in range(f[r],q*q+1)))

def total(q,n):
    nd=q*q+q+1
    return 10*((q+1)*n*n+q*q*n)-150*nd*n-17*nd*(n+q*q)+(560+shift(q))*q*nd+143*nd*(q+1)

def caps(q,n):
    f=minima(q)['minima'];B=total(q,n)
    return {str(r):[k for k in range(f[r],q*q+1) if polynomial(r,k)+shift(q)<=B] for r in range(q+2)}

def prepare():
    from structure import planes_of
    from planecount_check import witness
    results=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));q=data['q'];ag=AG(3,q);S={ag.pidx[tuple(p)] for p in data['points']};n=len(S);W=witness(ag,S)
        rk=[(sum(ag.lines[l][1]<=ps for l in W),len(S&ps)) for _,_,ps in planes_of(ag)]
        nd=q*q+q+1
        assert len(rk)==q*nd and sum(k for r,k in rk)==nd*n
        assert sum(k*k for r,k in rk)==(q+1)*n*n+q*q*n
        assert sum(r for r,k in rk)==nd*(q+1)
        assert sum(r*k for r,k in rk)==nd*(n+q*q)
        costs=[polynomial(r,k)+shift(q) for r,k in rk];B=total(q,n)
        assert min(costs)>=0 and sum(costs)==B and max(costs)<=B
        results.append(dict(file=fn,q=q,n=n,shift=shift(q),planes=len(rk),total_cost=B,max_plane_cost=max(costs),violations=0))
        print('BUDGET14_VALIDATE',results[-1],flush=True)
    assert shift(5)==0 and total(5,46)==541
    record=dict(family=FAMILY,status='EXACT_PLANE_BUDGET14_VERIFIED',q=5,n=46,budget=541,shift=0,plane_cost='10*k^2-(150+17*r)*k+560+143*r',total_cost='60*n^2-4927*n+100223',enumerations=[minima(q) for q in (3,4,5)],admitted_k=caps(5,46),results=results,argument='All plane costs are nonnegative; their sum is the stated fixed budget. A plane occurs as one actual plane, so its cost cannot exceed the total. The type filter is a necessary condition for actual witness systems, using integer plane populations.')
    atomic('out/budget14_validation.json',record);print('EXACT_PLANE_BUDGET14_VERIFIED',record['admitted_k'],flush=True)

def guard():
    p=Path('out/budget14_validation.json');r=json.loads(p.read_text())
    assert (r['family'],r['status'],r['q'],r['n'],r['budget'],r['shift'])==(FAMILY,'EXACT_PLANE_BUDGET14_VERIFIED',5,46,541,0)
    assert r['enumerations'][2]==minima(5) and r['admitted_k']==caps(5,46)
    assert total(5,46)==541 and shift(5)==0
    return hashlib.sha256(p.read_bytes()).hexdigest()

if __name__=='__main__':prepare();guard()
