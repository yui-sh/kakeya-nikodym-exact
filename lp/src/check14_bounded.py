"""Feasible/infeasible controls for homogeneous numerical search and exact restoration."""
import os
os.environ['XIV_HOMOG_BOUNDS']='1'
from types import MethodType
from planecount5 import LP
from homogeneous14 import HomogeneousMaster
from lift11 import atomic
import coverage14 as c
records=[]
for solver in ('simplex','ipm'):
 for feasible in (True,False):
    lp=LP();m=lp.var('m1');q=lp.var('Q14:1');rd=lp.var('RD14:1:1')
    lp.add([(m,1024)],'=',1024,'population');lp.add([(q,1)],'=',1,'normalization');lp.add([(rd,1)],'=',1,'marginal normalization')
    lp.add([(m,16),(q,-16)],'=',0,'link')
    if not feasible:lp.add([(m,1)],'<',0,'incompatible')
    for j in range(12):lp.add([(m,2**j)],'<',100*2**j,'XIV upper synthetic')
    lp.matrices=MethodType(c.fast_matrices,lp);master=HomogeneousMaster(lp);master.solver.setOptionValue('solver',solver)
    if solver=='ipm':master.solver.setOptionValue('run_crossover','off')
    status,data,detail=master.solve(lp,30)
    if feasible:assert status=='numerical-feasible' and c.numerical_feasibility(lp,data)[0]
    else:assert status=='LP-certified' and c.verify(lp,data)['rhs']<0
    records.append(dict(solver=solver,feasible=feasible,status=status,detail=detail))
atomic('out/bounded14_check.json',dict(status='BOUNDED14_CHECK_PASSED',records=records));print('BOUNDED14_CHECK_PASSED',flush=True)
