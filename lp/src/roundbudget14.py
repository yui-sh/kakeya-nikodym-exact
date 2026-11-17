"""Integer-plane rounding of the proven nonnegative quadratic budget."""
import json,hashlib
from pathlib import Path
from budget14 import polynomial,total,shift,guard as budget_guard
from lift11 import atomic

FAMILY='rounded-nonnegative-plane-budget-v1'


def denominators(B):
    # With fixed floor(B/t), the smallest t dominates all larger t in that
    # interval, coefficient by coefficient, for every nonnegative plane cost.
    out=[];last=None
    for t in range(2,B+1):
        bound=B//t
        if bound!=last:out.append(t);last=bound
    return out


def prepare():
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    ds=denominators(541);records=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json',
               'solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));q=data['q'];ag=AG(3,q)
        S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S)
        costs=[]
        for _,_,ps in planes_of(ag):
            r=sum(ag.lines[l][1]<=ps for l in W);k=len(S&ps)
            costs.append(polynomial(r,k)+shift(q))
        B=total(q,len(S));assert min(costs)>=0 and sum(costs)==B
        for t in ds:assert sum(c//t for c in costs)<=B//t
        records.append(dict(file=fn,total_cost=B,shift=shift(q),checks=len(ds),violations=0))
    # Check coefficientwise domination for every omitted denominator and
    # every possible integer plane cost under the n=46 cap.
    witnesses={541//t:t for t in ds};checks=0
    for t in range(2,542):
        representative=witnesses[541//t];assert representative<=t
        for cost in range(542):assert cost//representative>=cost//t;checks+=1
    record=dict(family=FAMILY,status='ZERO_VIOLATIONS',q=5,n=46,budget=541,
                budget_sha256=budget_guard(),denominators=ds,results=records,
                domination_checks=checks,
                argument='Actual planes occur with integer multiplicities. Sum_pi floor(c_pi/t) is an integer at most (sum_pi c_pi)/t=541/t, hence at most floor(541/t). Grouping identical full types gives the stated linear inequality. All costs are constant on (r,k), independent of degree assignments and of global m.',
                validation_scope='Known sets use their actual q,n budget and proven nonnegative cost shift; the n=46 RHS is not imposed on n=53 examples.')
    atomic('out/roundbudget14_validation.json',record)
    print('ROUNDBUDGET14_VALIDATED',len(ds),'rows',records,flush=True)


def guard():
    p=Path('out/roundbudget14_validation.json');r=json.loads(p.read_text())
    assert (r['family'],r['status'],r['q'],r['n'],r['budget'])==(FAMILY,'ZERO_VIOLATIONS',5,46,541)
    assert r['budget_sha256']==budget_guard() and r['denominators']==denominators(541)
    return hashlib.sha256(p.read_bytes()).hexdigest()


def add_cuts(lp,types):
    guard();costs=[polynomial(t[0],t[1]) for t in types]
    assert all(0<=v<=541 for v in costs)
    for t in denominators(541):
        row=[(lp.col[f'N{i}'],cost//t) for i,cost in enumerate(costs) if cost>=t]
        lp.add(row,'<',541//t,f'XIV rounded plane budget denominator {t}')
    print('ROUNDBUDGET14_CUTS',len(denominators(541)),flush=True)


if __name__=='__main__':prepare();guard()
