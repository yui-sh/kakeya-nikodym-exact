"""Exact feasible/infeasible regression for scaled numerical rows and restored duals."""
import os
os.environ['XIV_PRESOLVE']='on'
from types import MethodType
from planecount5 import LP
import coverage14 as c
from scaled14 import ScaledMaster
from lift11 import atomic
records=[]
for method in ('direct','phase1'):
 for feasible in (True,False):
    lp=LP();m=lp.var('m1');q=lp.var('Q14:1');r=lp.var('RD14:1:1')
    lp.add([(m,2**20)],'=',2**20,'population')
    lp.add([(q,1)],'=',1,'normalization')
    lp.add([(r,1)],'=',1,'marginal normalization')
    lp.add([(m,2**12),(q,-2**12)],'=',0,'link')
    if not feasible:lp.add([(m,2**18)],'<',0,'incompatible')
    for i in range(12):lp.add([(m,2**i)],'<',100*2**i,'XIV upper synthetic')
    lp.matrices=MethodType(c.fast_matrices,lp)
    master=ScaledMaster(lp,method);status,data,detail=master.solve(lp,30)
    if feasible:assert status=='numerical-feasible' and c.numerical_feasibility(lp,data)[0]
    else:assert status=='LP-certified' and c.verify(lp,data)['rhs']<0
    records.append(dict(method=method,feasible=feasible,status=status,detail=detail))
atomic('out/scaled14_check.json',dict(status='SCALED14_CHECK_PASSED',records=records));print('SCALED14_CHECK_PASSED',flush=True)
