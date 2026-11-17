"""Numerical power-of-two row scaling; certificates stay in original integers."""
import runpy,sys
from types import SimpleNamespace
import numpy as np
from scipy import sparse
import coverage14 as c
Base=c.Master14
class Proxy:
    def __init__(self,solver,factors,lp):self.inner=solver;self.factors=factors;self.lp=lp
    def __getattr__(self,n):return getattr(self.inner,n)
    def changeRowBounds(self,i,lo,hi):return self.inner.changeRowBounds(i,lo*self.factors[i],hi*self.factors[i])
    def getDualRay(self):
        # HiGHS can silently solve the whole LP again to obtain a named ray.
        # Avoid that only when its current original-row candidate ALREADY has
        # an exact all-column Farkas repair. The outer checker repeats it.
        if self.inner.getOptionValue('solver')[1]=='ipm':
            dual=np.asarray(self.inner.getSolution().row_dual)
            if len(dual)==len(self.lp.rows) and np.all(np.isfinite(dual)) and np.any(dual):
                from plane10_cert import bump,combo
                from scipy.optimize._highspy._core import HighsStatus
                dual=dual*self.factors;by,bc,_=bump(self.lp,{}, {})
                for scale in (10**7,10**9,10**11):
                    y=[int(round(-v*scale)) for v in dual]
                    for i,row in enumerate(self.lp.rows):
                        if row[2]=='<':y[i]=max(0,y[i])
                    cc,_=combo(self.lp,y);repair=max(0,max((-v+bc[j]-1)//bc[j] for j,v in enumerate(cc)))
                    yy=[a+repair*b for a,b in zip(y,by)];cc,rhs=combo(self.lp,yy)
                    if rhs<0 and min(cc)>=0:
                        print('SCALED14_PRE_RAY_EXACT_CANDIDATE',rhs,min(cc),flush=True)
                        return HighsStatus.kOk,False,np.zeros(len(self.lp.rows))
        status,exists,y=self.inner.getDualRay();return status,exists,np.asarray(y)*self.factors if exists else y
    def getSolution(self):
        x=self.inner.getSolution()
        dual=np.asarray(x.row_dual);values=np.asarray(x.row_value)
        return SimpleNamespace(col_value=x.col_value,col_dual=x.col_dual,row_dual=dual*self.factors if len(dual)==len(self.factors) else dual,row_value=values/self.factors if len(values)==len(self.factors) else values)
class ScaledMaster(Base):
    def __init__(self,lp,method='phase1'):
        original=lp.matrices;A,b,sense=original()
        maximum=np.asarray(abs(A).max(axis=1).toarray()).reshape(-1)
        maximum=np.maximum(1,np.maximum(maximum,abs(b)))
        factors=np.exp2(-np.ceil(np.log2(maximum)))
        assert np.all(factors>0)
        scaled=sparse.diags(factors)@A
        lp.matrices=lambda:(scaled,b*factors,sense)
        try:super().__init__(lp,method)
        finally:lp.matrices=original
        self.solver=Proxy(self.solver,factors,lp)
        print('SCALED14_ORIGINAL_INTEGER_REPLAY',len(factors),float(min(factors)),float(max(factors)),flush=True)
    def solve(self,lp,seconds):
        status,data,detail=super().solve(lp,seconds)
        detail['numerical_row_scaling']='positive powers of two; duals restored to original integer rows before all exact checks'
        return status,data,detail
if __name__=='__main__':
    c.Master14=ScaledMaster
    if '--projected' in sys.argv:
        sys.argv.remove('--projected');runpy.run_path('thin14.py',run_name='__main__')
    else:runpy.run_path('probe14.py',run_name='__main__')
