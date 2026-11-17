from planecount5 import LP
from coverage14 import Master14,fast_matrices,verify
from types import MethodType
for feasible in (True,False):
    lp=LP();x=lp.var('m1');lp.add([(x,1)],'=',1,'population');lp.add([(x,-1 if feasible else 1)],'<',0,'homogeneous');lp.add([(x,-1)],'<',-1,'configuration')
    for i in range(12):lp.add([(x,1)],'<',100,'bound'+str(i))
    lp.matrices=MethodType(fast_matrices,lp);m=Master14(lp);status,data,detail=m.solve(lp,10)
    if feasible:assert status=='numerical-feasible' and abs(data[0]-1)<1e-8
    else:assert status=='LP-certified' and verify(lp,data)['rhs']<0
    print('NONZERO_PHASE14_CHECK',feasible,status,detail)
print('NONZERO_PHASE14_CHECK_PASSED')
