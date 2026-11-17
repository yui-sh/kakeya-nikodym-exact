"""Joint plane/pencil column generation; no uncertified exclusions."""
import json,sys,time,os
from pathlib import Path
import numpy as np
from scipy import sparse as sp
from scipy.optimize._highspy._core import _Highs,HighsLp,MatrixFormat,HighsModelStatus,HighsStatus
from plane10 import build,counts,PlaneOracle
from pencil9 import Oracle

def solve(short):
    full=os.environ.get('PLANE10_FULL')=='1'
    tag='plane10_'+('full' if full else '')+short;old=json.load(open('out/pencil9_'+short+'.json'));config=tuple(old['config']);m=tuple(old['m']);t0=time.monotonic()
    lp,info,types,grow,zrow,nr,gr=build(config,m);A,b,sense=lp.matrices();base=A.shape[1]
    eq=np.array([i for i,s in enumerate(sense) if s=='=']);ub=np.array([i for i,s in enumerate(sense) if s=='<'])
    po=PlaneOracle(types,info,nr,gr,tag);wo=Oracle(info,m,grow,zrow,tag)
    columns=[];seen=set()
    def add(priced,seed=False):
        rr=[];cc=[];vv=[];num=0
        for cost,col in priced:
            if cost>=0 and not seed:continue
            key=(col[0],col[1],tuple(col[2]))
            if key in seen:continue
            seen.add(key);columns.append(col)
            for r,c in counts(col,grow,zrow,nr,gr).items():rr.append(r);cc.append(num);vv.append(c)
            num+=1
        return sp.csc_matrix((vv,(rr,cc)),shape=(len(b),num))
    seeds=[(0,['pencil',z,gs]) for z,gs in json.load(open('out/pencil9_'+short+'_columns.json'))]
    seeds+=po.price([0]*len(b),first=True)
    if full:
        dtype=np.dtype([('type','<u4'),('geo','<u4'),('gids','<u2',(6,)),('degrees','<u8')])
        records=np.memmap(f'out/plane10_catalog_{short}.bin',dtype=dtype,mode='r')
        assert 'CATALOG_COMPLETE' in Path(f'out/plane10_catalog_{short}.log').read_text()
        for rec in records:
            i=int(rec['type']);gi=int(rec['geo']);packed=int(rec['degrees']);k=types[i][1]
            seeds.append((0,['plane',i,[int(g) for g in rec['gids'] if g!=65535],gi,[(packed>>(3*p))&7 for p in range(k)]]))
        print('FULL_PLANE_SEEDS',len(records),flush=True)
    checkpoint=Path(f'out/{tag}_checkpoint_columns.json')
    if checkpoint.exists():
        resumed=json.load(checkpoint.open())
        for col in resumed:
            if col[0]=='plane':po.verify(col)
        seeds.extend((0,col) for col in resumed)
        print('RESUME_COLUMNS',len(resumed),flush=True)
    C=add(seeds,True)
    print('BUILD',A.shape,'plane types',len(nr),'seed columns',len(columns),'seconds',time.monotonic()-t0,flush=True)
    # Keep the master and its artificial variables in place as columns are added.
    arts=sp.csc_matrix((np.r_[np.ones(len(eq)),-np.ones(len(eq)),-np.ones(len(ub))],(np.r_[eq,eq,ub],np.arange(2*len(eq)+len(ub)))),shape=(len(b),2*len(eq)+len(ub)))
    na=arts.shape[1];mat=sp.hstack([A,arts,C],format='csc');nv=mat.shape[1]
    hp=HighsLp();hp.num_col_=nv;hp.num_row_=len(b)
    hp.col_cost_=np.r_[np.zeros(base),np.ones(na),np.zeros(C.shape[1])]
    hp.col_lower_=np.zeros(nv);hp.col_upper_=np.full(nv,np.inf)
    lower=b.copy();lower[ub]=-np.inf;hp.row_lower_=lower;hp.row_upper_=b
    hp.a_matrix_.format_=MatrixFormat.kColwise;hp.a_matrix_.start_=mat.indptr;hp.a_matrix_.index_=mat.indices;hp.a_matrix_.value_=mat.data
    solver=_Highs();solver.setOptionValue('threads',1);solver.setOptionValue('output_flag',False)
    solver.setOptionValue('time_limit',120.0);solver.setOptionValue('solver','ipm');solver.passModel(hp)
    recs=[];status='iteration_limit';x=None
    for iteration in range(300):
        json.dump(columns,open(f'out/{tag}_checkpoint_columns.json','w'))
        ts=time.monotonic();solver.run();st=solver.getModelStatus()
        if st==HighsModelStatus.kTimeLimit:
            # Retry this same master with interior point; no timeout is a conclusion.
            print('MASTER_TIMEOUT_RETRY_IPM',iteration,flush=True)
            solver.clearSolver();solver.setOptionValue('solver','ipm');solver.setOptionValue('time_limit',solver.getRunTime()+300.0);solver.run();st=solver.getModelStatus()
        rec=dict(iteration=iteration,columns=len(columns),status=str(st),solve_s=time.monotonic()-ts);recs.append(rec)
        if st!=HighsModelStatus.kOptimal:status='solver_failure';break
        sol=solver.getSolution();rec['phase1']=float(solver.getObjectiveValue());print('MASTER',rec,flush=True)
        if rec['phase1']<1e-7:
            xx=np.array(sol.col_value);x=np.r_[xx[:base],xx[base+na:]];status='computed_feasible_pending_exact';break
        y=-np.array(sol.row_dual);yi=[int(round(v*10**9)) for v in y]
        np.savez_compressed(f'out/{tag}_dual.npz',y=y)
        ts=time.monotonic();priced=[(c,['pencil',z,gs]) for c,z,gs in wo.price(yi)]
        pp=po.price(yi);priced+=pp;new=add(priced)
        if full:assert all(c>=-10000 for c,col in pp), 'full-plane dual has a large negative plane cost'
        rec.update(pricing_s=time.monotonic()-ts,new_columns=new.shape[1],min_price=min(c for c,_ in priced));print('PRICING',rec,flush=True)
        if not new.shape[1] or (rec['phase1']>1e-5 and rec['min_price']>-rec['phase1']*10**9/100):
            from plane10_cert import run as certify
            try:
                certify(('full' if full else '')+short)
                certify(('full' if full else '')+short,replay=True)
                status='LP-certified infeasible';break
            except RuntimeError as error:
                print('CERTIFICATE_NOT_OBTAINED',str(error),flush=True)
        if not new.shape[1]:status='uncertified_stall';break
        nn=new.shape[1]
        assert solver.addCols(nn,np.zeros(nn),np.zeros(nn),np.full(nn,np.inf),new.nnz,new.indptr,new.indices,new.data)==HighsStatus.kOk
        solver.setOptionValue('solver','ipm');solver.setOptionValue('time_limit',solver.getRunTime()+120.0)
    result=dict(config=config,m=m,status=status,elapsed_s=time.monotonic()-t0,base_variables=base,rows=len(b),columns=len(columns),plane_jobs=len(po.jobs),rounds=recs)
    json.dump(result,open(f'out/{tag}.json','w'),indent=2);json.dump(columns,open(f'out/{tag}_columns.json','w'))
    if x is not None:np.savez_compressed(f'out/{tag}_primal.npz',x=x)
    print('DONE',json.dumps(result),flush=True)
if __name__=='__main__':solve(sys.argv[1])
