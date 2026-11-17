"""Synthetic regression: an omitted hull vertex must never give false credit."""
import os,random,json
os.environ['XIV_Q_SAMPLE']='1'
os.environ['XIV_PRESOLVE']='on'
from types import MethodType
from planecount5 import LP
import coverage14 as c
from thin14 import ProjectedMaster
from lift11 import atomic

picked=random.Random(140007).sample(list(range(4)),1)[0]
target=next(i for i in (1,2) if i!=picked)
records=[]
for full_feasible in (True,False):
    lp=LP();m=lp.var('m1')
    qs=[lp.var(f'Q14:{j},{3-j}') for j in range(4)]
    lp.add([(m,1)],'=',1,'population')
    lp.add([(q,1) for q in qs],'=',1,'normalization')
    lp.add([(q,j) for j,q in enumerate(qs)],'=',target,'first moment')
    lp.add([(q,j*j) for j,q in enumerate(qs)],'=',target*target,'second moment')
    if not full_feasible:lp.add([(m,1)],'<',0,'incompatible population')
    for i in range(12):lp.add([(m,1)],'<',100,'XIV upper synthetic '+str(i))
    if full_feasible:
        x=[0]*len(lp.names);x[m]=1;x[qs[target]]=1
        for row,b,s,name in lp.rows:
            a=sum(v*x[j] for j,v in row.items())
            assert a==b if s=='=' else a<=b
    lp.matrices=MethodType(c.fast_matrices,lp)
    c.CERT_SUFFIX='_restricted_check'
    master=ProjectedMaster(lp,'direct')
    assert qs[target] in master.omitted_q
    status,data,detail=master.solve(lp,10)
    if full_feasible:
        assert status=='open' and detail['negative_Q_columns']>=1
    else:
        assert status=='LP-certified';c.verify(lp,data)
    records.append(dict(full_feasible=full_feasible,status=status,detail=detail))
    print('RESTRICTED14_CHECK',records[-1],flush=True)
atomic('out/restricted14_check.json',dict(status='RESTRICTED14_CHECK_PASSED',records=records))
print('RESTRICTED14_CHECK_PASSED',flush=True)
