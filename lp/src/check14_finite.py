"""Regression: finite column pricing must add a missing feasible vertex."""
from types import MethodType
from planecount5 import LP
import coverage14 as c
from finite14 import FiniteMaster
from lift11 import atomic
from pathlib import Path
Path('out/finite14_finite_check_seeds.json').unlink(missing_ok=True)
records=[]
for feasible in (True,False):
    lp=LP();m=lp.var('m1');qs=[lp.var(f'Q14:{j},{3-j}') for j in range(4)]
    lp.add([(m,1)],'=',1,'population')
    lp.add([(q,1) for q in qs],'=',1,'normalization')
    lp.add([(q,j) for j,q in enumerate(qs)],'=',1,'first moment')
    lp.add([(q,j*j) for j,q in enumerate(qs)],'=',1,'second moment')
    if not feasible:lp.add([(m,1)],'<',0,'incompatible population')
    for i in range(12):lp.add([(m,1)],'<',100,'XIV upper synthetic '+str(i))
    lp.matrices=MethodType(c.fast_matrices,lp);c.CERT_SUFFIX='_finite_check'
    master=FiniteMaster(lp);status,data,detail=master.solve(lp,30)
    if feasible:
        assert status=='numerical-feasible' and c.numerical_feasibility(lp,data)[0]
        assert abs(data[qs[1]]-1)<1e-8 and master.rounds>=1
    else:assert status=='LP-certified' and c.verify(lp,data)['rhs']<0
    records.append(dict(feasible=feasible,status=status,detail=detail))
    print('FINITE14_CHECK',records[-1],flush=True)
atomic('out/finite14_check.json',dict(status='FINITE14_CHECK_PASSED',records=records))
print('FINITE14_CHECK_PASSED',flush=True)
