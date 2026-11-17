"""Numerical column generation over an already-complete finite LP catalogue.

No geometric columns or types are omitted from exact verification. Restricted
masters only guide numerical search; every accepted ray is checked on ALL input
columns, then thin14 (if used) checks all original columns again. Progress files
do not cache leaf-local geometric indices or claim certificates.
"""
import os,time,runpy,sys,json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize._highspy._core import HighsBasis,HighsBasisStatus,HighsStatus
import coverage14 as c
from coverage11 import _Highs,HighsLp,MatrixFormat,HighsModelStatus
from plane10_cert import bump,combo
from lift11 import atomic,key


class FiniteMaster:
    def __init__(self,lp,method='phase1'):
        assert method=='phase1'
        self.lp=lp;self.base=len(lp.names);self.method='finite-phase1'
        A,b,sense=lp.matrices();self.A=A.tocsc();self.b=b;self.sense=sense
        self.ub=np.array([i for i,s in enumerate(sense) if s=='<'],dtype=np.int32)
        eq=np.array([i for i,s in enumerate(sense) if s=='=' and b[i]!=0],dtype=np.int32)
        au=np.array([i for i in self.ub if b[i]!=0],dtype=np.int32)
        self.canonical=getattr(lp,'canonical14',{n:('Q',tuple(map(int,n[4:].split(',')))) for n in lp.names if n.startswith('Q14:')})
        self.metadata=getattr(lp,'model14',dict(q=5,D=2,config='synthetic',prefix=[]))
        self.seedpath=Path('out/finite14'+c.CERT_SUFFIX+'_seeds.json')
        lookup={key(self.canonical[n]):j for j,n in enumerate(lp.names) if n in self.canonical}
        assert len(lookup)==len(self.canonical),'Canonical finite-column keys must be unique'
        # Other count variables are reconstructed freshly, never cached by a
        # leaf-local index. Geometric/count-hull seeds use semantic keys only.
        seeds=[j for j,n in enumerate(lp.names) if n not in self.canonical and not n.startswith(('X','Y'))]
        source=Path(os.environ.get('XIV_FINITE_SEEDS',str(self.seedpath)))
        if source.exists():
            old=json.loads(source.read_text())
            assert old['family']=='finite14-canonical-active-columns-v1'
            assert (old['model']['q'],old['model']['D'])==(self.metadata['q'],self.metadata['D'])
            seeds+= [lookup[key(v)] for v in old['columns'] if key(v) in lookup]
        seeds=sorted(set(seeds))
        self.active=np.zeros(self.base,dtype=bool);self.active[seeds]=True
        self.positions={j:k for k,j in enumerate(seeds)}
        na=2*len(eq)+len(au)
        arts=sparse.csc_matrix((np.r_[np.ones(len(eq)),-np.ones(len(eq)),-np.ones(len(au))],
              (np.r_[eq,eq,au],np.arange(na))),shape=(len(b),na))
        mat=sparse.hstack([self.A[:,seeds],arts],format='csc');nv=mat.shape[1]
        hp=HighsLp();hp.num_col_=nv;hp.num_row_=len(b)
        hp.col_cost_=np.r_[np.zeros(len(seeds)),np.ones(na)]
        hp.col_lower_=np.zeros(nv);hp.col_upper_=np.full(nv,np.inf)
        lower=b.copy();lower[self.ub]=-np.inf;hp.row_lower_=lower;hp.row_upper_=b
        hp.a_matrix_.format_=MatrixFormat.kColwise;hp.a_matrix_.start_=mat.indptr
        hp.a_matrix_.index_=mat.indices;hp.a_matrix_.value_=mat.data
        self.solver=_Highs();self.solver.setOptionValue('threads',1)
        self.solver.setOptionValue('output_flag',False)
        assert self.solver.passModel(hp)==HighsStatus.kOk
        basis=HighsBasis();cs=[HighsBasisStatus.kLower]*nv;rs=[HighsBasisStatus.kBasic]*len(b)
        values=np.zeros(nv)
        for k,i in enumerate(eq):
            j=len(seeds)+k+(0 if b[i]>0 else len(eq))
            cs[j]=HighsBasisStatus.kBasic;rs[i]=HighsBasisStatus.kUpper;values[j]=abs(b[i])
        for k,i in enumerate(au):
            if b[i]<0:
                j=len(seeds)+2*len(eq)+k
                cs[j]=HighsBasisStatus.kBasic;rs[i]=HighsBasisStatus.kUpper;values[j]=-b[i]
        activity=np.asarray(mat@values).reshape(-1)
        assert np.all(activity>=lower) and np.all(activity<=b)
        basis.col_status=cs;basis.row_status=rs
        assert self.solver.setBasis(basis)==HighsStatus.kOk
        self.solver.setOptionValue('simplex_strategy',4)
        self.solver.setOptionValue('presolve','off')
        self.bump=None;self.rounds=0
        print('FINITE14_INITIAL',len(seeds),'of',self.base,'rows',len(b),'artificials',na,flush=True)

    def checkpoint(self,**detail):
        selected=[self.canonical[self.lp.names[j]] for j in np.flatnonzero(self.active)
                  if self.lp.names[j] in self.canonical]
        atomic(self.seedpath,dict(family='finite14-canonical-active-columns-v1',model=self.metadata,
                                 columns=sorted(selected,key=key)))
        atomic('out/finite14'+c.CERT_SUFFIX+'.json',dict(
            method='numerical column generation over all existing finite columns',
            rows=len(self.lp.rows),full_columns=self.base,active_columns=int(self.active.sum()),
            rounds=self.rounds,resume=str(self.seedpath),
            **detail))

    def exact_candidate(self,dual):
        if self.bump is None:self.bump=bump(self.lp,{}, {})[:2]
        by,bc=self.bump
        for scale in (10**7,10**9,10**11):
            y=[int(round(-v*scale)) for v in dual]
            for i in self.ub:y[i]=max(0,y[i])
            cc,_=combo(self.lp,y)
            repair=max(0,max((-v+bc[j]-1)//bc[j] for j,v in enumerate(cc)))
            yy=[a+repair*b for a,b in zip(y,by)]
            cc,rhs=combo(self.lp,yy)
            if rhs<0 and min(cc)>=0:
                return yy,dict(c.verify(self.lp,yy),scale=scale,bump=repair)
        return None

    def solve(self,lp,seconds):
        assert lp is self.lp
        deadline=time.monotonic()+seconds;self.solver.setOptionValue('solver','simplex')
        while time.monotonic()<deadline and not c.stop:
            self.solver.setOptionValue('time_limit',self.solver.getRunTime()+min(60,deadline-time.monotonic()))
            self.solver.run();status=self.solver.getModelStatus();sol=self.solver.getSolution()
            objective=self.solver.getObjectiveValue()
            detail=dict(method='finite-column-phase1',solver_status=str(status),phase1=objective,
                        active_columns=int(self.active.sum()),full_columns=self.base,rounds=self.rounds)
            if status==HighsModelStatus.kOptimal and objective<1e-7:
                x=np.zeros(self.base)
                for j,k in self.positions.items():x[j]=sol.col_value[k]
                valid,residual=c.numerical_feasibility(lp,x);detail['residual_check']=residual
                if valid:
                    positive=[j for j,n in enumerate(lp.names) if n.startswith('Q14:') and x[j]>1e-8]
                    detail['positive_degree_vectors']=[dict(m=list(map(int,lp.names[j][4:].split(','))),weight=float(x[j]))
                        for j in sorted(positive,key=lambda j:-x[j])[:32]]
                    self.checkpoint(status='numerical-feasible',detail=detail)
                    return 'numerical-feasible',x,detail
                self.checkpoint(status='open',reason='Numerical point failed full input-row check',detail=detail)
                return 'open',None,detail
            dual=np.asarray(sol.row_dual,dtype=float)
            if len(dual)!=len(lp.rows) or not np.all(np.isfinite(dual)):
                self.checkpoint(status='open',reason='No finite dual',detail=detail)
                return 'open',None,detail
            costs=np.asarray(self.A.T@(-dual)).reshape(-1)
            detail['minimum_numeric_column']=float(min(costs))
            if status==HighsModelStatus.kOptimal and min(costs)>=-1e-6:
                exact=self.exact_candidate(dual)
                if exact:
                    y,verified=exact;detail.update(verified)
                    self.checkpoint(status='exact Farkas on every input column; outer projection replay still required',detail=detail)
                    print('FINITE14_ALL_INPUT_COLUMNS_VERIFIED',detail,flush=True)
                    return 'LP-certified',y,detail
            if status!=HighsModelStatus.kOptimal:
                self.checkpoint(status='running' if status==HighsModelStatus.kTimeLimit else 'open',detail=detail)
                print('FINITE14_SOLVER_INTERVAL',detail,flush=True)
                if status==HighsModelStatus.kTimeLimit:continue
                return 'open',None,detail
            missing=np.flatnonzero((costs < -1e-7)&(~self.active))
            if not len(missing):
                self.checkpoint(status='open',reason='No negative missing columns; exact repair failed',detail=detail)
                return 'open',None,detail
            chosen=set(missing[np.argsort(costs[missing])[:512]].tolist())
            for prefix in ('N','G','Z','Q','X','Y','C','W','V'):
                group=np.array([j for j in missing if lp.names[j].startswith(prefix)],dtype=np.int32)
                if len(group):chosen.update(group[np.argsort(costs[group])[:128]].tolist())
            chosen=sorted(chosen);new=self.A[:,chosen].tocsc();start=self.solver.getNumCol()
            nn=len(chosen)
            assert self.solver.addCols(nn,np.zeros(nn),np.zeros(nn),np.full(nn,np.inf),
                                      new.nnz,new.indptr,new.indices,new.data)==HighsStatus.kOk
            for k,j in enumerate(chosen):self.positions[j]=start+k
            self.active[chosen]=True;self.rounds+=1
            self.checkpoint(status='running',added_columns=nn,detail=detail)
            print('FINITE14_PRICE_ROUND',self.rounds,'added',nn,'active',int(self.active.sum()),
                  'objective',objective,'minimum',float(min(costs)),flush=True)
        detail=dict(method='finite-column-phase1',reason='time limit or requested stop',rounds=self.rounds,
                    active_columns=int(self.active.sum()),full_columns=self.base)
        self.checkpoint(status='open',detail=detail)
        return 'open',None,detail


if __name__=='__main__':
    os.environ['XIV_METHOD']='phase1'
    os.environ['XIV_FINITE_COLUMNS']='1'
    c.Master14=FiniteMaster
    if '--full-model' in sys.argv:
        sys.argv.remove('--full-model');runpy.run_path('probe14.py',run_name='__main__')
    else:
        runpy.run_path('thin14.py',run_name='__main__')
