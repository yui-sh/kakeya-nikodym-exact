"""Bounded homogeneous feasibility search; exact checks use original integer rows.

Replace A x <=/= b by A x-b*t <=/= 0, with 0<=t<=1, and minimize -t.
The zero vector is feasible. The optimum is -1 for a feasible original LP and
0 otherwise. At the latter optimum the original-row dual gives a Farkas
candidate without requesting a separate infeasibility-ray solve.
"""
import sys,runpy,os
from types import SimpleNamespace
import numpy as np
from scipy import sparse
from scipy.optimize._highspy._core import HighsStatus
from coverage11 import _Highs,HighsLp,MatrixFormat
from plane10_cert import bump,combo
import coverage14 as c
class View:
    def __init__(self,inner,factors):self.inner=inner;self.factors=factors
    def __getattr__(self,k):return getattr(self.inner,k)
    def getSolution(self):
        x=self.inner.getSolution();d=np.asarray(x.row_dual)
        return SimpleNamespace(col_value=x.col_value,row_dual=d*self.factors if len(d)==len(self.factors) else d)
class HomogeneousMaster:
    def __init__(self,lp,method='phase1'):
        self.lp=lp;self.base=len(lp.names);self.bump=None
        A,b,sense=lp.matrices();maximum=np.asarray(abs(A).max(axis=1).toarray()).reshape(-1)
        maximum=np.maximum(1,np.maximum(maximum,abs(b)));factors=np.exp2(-np.ceil(np.log2(maximum)))
        matrix=sparse.hstack([sparse.diags(factors)@A,sparse.csc_matrix((-b*factors).reshape(-1,1))],format='csc')
        hp=HighsLp();hp.num_col_=self.base+1;hp.num_row_=len(b)
        hp.col_cost_=np.r_[np.zeros(self.base),-1.0];hp.col_lower_=np.zeros(self.base+1);hp.col_upper_=np.r_[np.full(self.base,np.inf),1.0]
        hp.row_lower_=np.array([-np.inf if s=='<' else 0.0 for s in sense]);hp.row_upper_=np.zeros(len(b))
        hp.a_matrix_.format_=MatrixFormat.kColwise;hp.a_matrix_.start_=matrix.indptr;hp.a_matrix_.index_=matrix.indices;hp.a_matrix_.value_=matrix.data
        strategy=int(os.environ.get('XIV_HOMOG_STRATEGY','4'));assert strategy in (1,4)
        solver=_Highs();solver.setOptionValue('threads',1);solver.setOptionValue('output_flag',False);solver.setOptionValue('solver','simplex');solver.setOptionValue('simplex_strategy',strategy)
        assert solver.passModel(hp)==HighsStatus.kOk
        self.bound_rows=None
        if os.environ.get('XIV_HOMOG_BOUNDS')=='1':
            from bound14 import upper_bounds
            upper,detail,self.bound_rows=upper_bounds(lp)
            assert solver.changeColsBounds(self.base,np.arange(self.base,dtype=np.int32),np.zeros(self.base),np.asarray(upper,dtype=float))==HighsStatus.kOk
            print('HOMOGENEOUS14_REDUNDANT_BOUNDS',detail,flush=True)
        self.solver=View(solver,factors)
        print('HOMOGENEOUS14_FEASIBLE_ZERO',self.base,len(b),'objective bounded in [-1,0]',flush=True)
    def solve(self,lp,seconds):
        assert lp is self.lp
        self.solver.setOptionValue('time_limit',self.solver.getRunTime()+seconds)
        self.solver.run();status=self.solver.getModelStatus();solution=self.solver.getSolution()
        detail=dict(redundant_original_row_bounds=self.bound_rows is not None,simplex_strategy=self.solver.getOptionValue('simplex_strategy')[1],method='homogeneous bounded feasibility',solver_status=str(status),objective=float(self.solver.getObjectiveValue()),numerical_row_scaling='positive powers of two; original-row dual restored')
        values=np.asarray(solution.col_value,dtype=float)
        if len(values)==self.base+1 and np.all(np.isfinite(values)):
            t=float(values[-1]);detail['homogenizing_t']=t
            if t>1e-8:
                x=values[:self.base]/t;valid,residual=c.numerical_feasibility(lp,x);detail['residual_check']=residual
                if valid:
                    detail['degree_values']=[float(x[lp.col[f'm{j}']]) for j in range(1,8) if f'm{j}' in lp.col]
                    detail['positive_degree_vectors']=sorted([dict(m=list(map(int,n[4:].split(','))),weight=float(x[j])) for j,n in enumerate(lp.names) if n.startswith('Q14:') and x[j]>1e-9],key=lambda v:-v['weight'])
                    return 'numerical-feasible',x,detail
        dual=np.asarray(solution.row_dual,dtype=float)
        if len(dual)==len(lp.rows) and np.all(np.isfinite(dual)):
            by,bc,_=bump(lp,{}, {})
            for scale in (10**7,10**9,10**11):
                y=[int(round(-v*scale)) for v in dual]
                for i,row in enumerate(lp.rows):
                    if row[2]=='<':y[i]=max(0,y[i])
                if self.bound_rows is not None:
                    from bound14 import repair_bounds
                    y=repair_bounds(lp,y,self.bound_rows)
                cc,_=combo(lp,y);repair=max(0,max((-v+bc[j]-1)//bc[j] for j,v in enumerate(cc)))
                yy=[a+repair*b for a,b in zip(y,by)];cc,rhs=combo(lp,yy)
                if rhs<0 and min(cc)>=0:return 'LP-certified',yy,dict(detail,**c.verify(lp,yy),scale=scale,bump=repair)
        return 'open',None,dict(detail,reason='No verified original-model primal or exact Farkas candidate')
if __name__=='__main__':
    c.Master14=HomogeneousMaster
    if '--projected' in sys.argv:
        sys.argv.remove('--projected');runpy.run_path('thin14.py',run_name='__main__')
    else:runpy.run_path('probe14.py',run_name='__main__')
