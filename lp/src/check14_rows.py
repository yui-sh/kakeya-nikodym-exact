"""Regression: omitted cuts cannot turn a false full-model point into a survivor."""
from types import MethodType
from pathlib import Path
from planecount5 import LP
import coverage14 as c
from row14 import RowMaster
from lift11 import atomic
records=[]
for feasible in (True,False):
 c.CERT_SUFFIX='_rows_check_'+str(feasible);Path('out/rows14'+c.CERT_SUFFIX+'_active.json').unlink(missing_ok=True)
 lp=LP();m=lp.var('m1');lp.add([(m,1)],'=',1,'population')
 lp.add([(m,1)],'<',2 if feasible else 0,'XIV nonlinear [6,17,[1]] [1,2] 1')
 for i in range(12):lp.add([(m,1)],'<',100,'XIV upper synthetic')
 lp.matrices=MethodType(c.fast_matrices,lp);master=RowMaster(lp);status,data,detail=master.solve(lp,30)
 if feasible:assert status=='numerical-feasible' and c.numerical_feasibility(lp,data)[0]
 else:assert status=='LP-certified' and master.rounds>=1 and c.verify(lp,data)['rhs']<0
 records.append(dict(feasible=feasible,status=status,detail=detail))
atomic('out/rows14_check.json',dict(status='ROW14_CHECK_PASSED',records=records));print('ROW14_CHECK_PASSED',flush=True)
