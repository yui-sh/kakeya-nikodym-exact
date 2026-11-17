"""Adaptive selection from the already priced full local-cut family.

A smaller row system is a relaxation. Its Farkas vector is padded by zero and
checked on every original column. A numerical point is accepted only after ALL
original rows pass; otherwise violated original cuts are added and checkpointed
by canonical row name, never by leaf-local geometric indices.
"""
import json,time,runpy
from pathlib import Path
from types import MethodType,SimpleNamespace
import numpy as np
from planecount5 import LP
from homogeneous14 import HomogeneousMaster
import coverage14 as c
from lift11 import atomic
FAMILY='adaptive-full-priced-row-family-v1'
def local_family(name):
    return next((p for p in ('XIV subset ','XIV nonlinear ','XIV triple ','XIV extra ') if name.startswith(p)),None)
def initial(name):
    family=local_family(name)
    if not family:return True
    subset=json.loads(name.split()[3])
    return len(subset)<=(2 if family=='XIV subset ' else 1)
class OuterView:
    def __init__(self,master):self.master=master
    def __getattr__(self,k):return getattr(self.master.inner.solver,k)
    def getSolution(self):
        s=self.master.inner.solver.getSolution();dual=np.zeros(len(self.master.full.rows))
        if len(s.row_dual)==len(self.master.keep):dual[self.master.keep]=s.row_dual
        return SimpleNamespace(col_value=s.col_value,row_dual=dual)
class RowMaster:
    def __init__(self,lp,method='phase1'):
        self.full=lp;self.method=method;self.rounds=0
        self.D=max(int(n[1:]) for n in lp.names if n.startswith('m') and n[1:].isdigit())
        self.seed=Path('out/rows14'+c.CERT_SUFFIX+'_active.json')
        names=set()
        if self.seed.exists():
            old=json.load(self.seed.open());assert old['family']==FAMILY and old['D']==self.D;names=set(old['rows'])
        self.keep=sorted(i for i,row in enumerate(lp.rows) if initial(row[3]) or row[3] in names)
        self.rebuild();self.solver=OuterView(self);self.checkpoint('initialized')
        print('ROW14_INITIAL',len(self.keep),'of',len(lp.rows),'all columns',len(lp.names),flush=True)
    def rebuild(self):
        small=LP();small.names=self.full.names;small.col=self.full.col;small.rows=[self.full.rows[i] for i in self.keep]
        small.matrices=MethodType(c.fast_matrices,small);self.small=small;self.inner=HomogeneousMaster(small,self.method)
    def checkpoint(self,status,**detail):
        atomic(self.seed,dict(family=FAMILY,D=self.D,rows=[self.full.rows[i][3] for i in self.keep if local_family(self.full.rows[i][3])],status=status,rounds=self.rounds,full_rows=len(self.full.rows),active_rows=len(self.keep),**detail))
    def solve(self,lp,seconds):
        assert lp is self.full
        options={name:self.inner.solver.getOptionValue(name)[1] for name in ('solver','run_crossover','output_flag','log_to_console','log_file')}
        deadline=time.monotonic()+seconds
        while not c.stop and time.monotonic()<deadline:
            status,data,detail=self.inner.solve(self.small,max(1,deadline-time.monotonic()))
            detail.update(numerical_row_selection='adaptive subset of full validated local rows; all original columns retained',active_rows=len(self.keep),full_rows=len(lp.rows),row_rounds=self.rounds)
            if status=='LP-certified':
                y=[0]*len(lp.rows)
                for i,v in zip(self.keep,data):y[i]=v
                detail.update(c.verify(lp,y));self.checkpoint('full original-column Farkas verified',detail=detail)
                print('ROW14_FULL_COLUMN_VERIFIED',detail,flush=True);return status,y,detail
            if status!='numerical-feasible':
                self.checkpoint('open',detail=detail);return 'open',None,detail
            valid,residual=c.numerical_feasibility(lp,data);detail['full_original_rows']=residual
            if valid:
                self.checkpoint('numerical feasible on all original rows; exact joint check still required',detail=detail)
                return status,data,detail
            A,b,sense=lp.matrices();activity=np.asarray(A@data).reshape(-1);delta=activity-b
            violation=np.array([max(0,v) if s=='<' else abs(v) for v,s in zip(delta,sense)])/(1+abs(b));del A
            omitted=set(range(len(lp.rows)))-set(self.keep)
            bad=[i for i in omitted if violation[i]>1e-7]
            if not bad:
                self.checkpoint('open',reason='original-row residual failure without a violated omitted cut',detail=detail)
                return 'open',None,detail
            added=sorted(bad,key=lambda i:-violation[i])[:1000]
            assert all(local_family(lp.rows[i][3]) for i in added)
            self.keep=sorted(set(self.keep)|set(added));self.rounds+=1
            self.checkpoint('adding violated original cuts',added=len(added),violated=len(bad),detail=detail)
            print('ROW14_ADD_VIOLATED',len(added),'of',len(bad),'active',len(self.keep),flush=True)
            if time.monotonic()>=deadline:break
            self.rebuild()
            for name,value in options.items():self.inner.solver.setOptionValue(name,value)
        self.checkpoint('open',reason='time limit or requested stop')
        return 'open',None,dict(reason='row generation time limit or requested stop',row_rounds=self.rounds,active_rows=len(self.keep),full_rows=len(lp.rows))
if __name__=='__main__':
    c.Master14=RowMaster;runpy.run_path('probe14.py',run_name='__main__')
