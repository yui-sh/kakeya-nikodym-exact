"""Separate D=6/7 exact prefix coverage; no writes to the D=5 ledger."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import argparse,json,time,copy,signal,datetime,fcntl,hashlib
from pathlib import Path
from collections import Counter
from functools import lru_cache
import numpy as np
import planecount6 as p6
from coverage11 import Master,verify as _verify
from geometry14 import freeze
from lift14 import add_cuts,table,FAMILY
from lift11 import atomic,key
HARD_STOP=datetime.datetime(2026,11,4,11,15,tzinfo=datetime.timezone.utc).timestamp()
# Leave time to rebuild certificates for the final audit.
CAMPAIGN_STOP=datetime.datetime(2026,11,4,9,30,tzinfo=datetime.timezone.utc).timestamp()
stop=False
CERT_SUFFIX=""

def verify(lp,y):
    assert all(type(b)==int and all(type(c)==int for c in row.values()) for row,b,_,_ in lp.rows)
    return _verify(lp,y)

def numerical_feasibility(lp,x):
    if len(x)!=len(lp.names) or not np.all(np.isfinite(x)):
        return False,dict(reason='nonfinite or wrong-sized numerical point')
    A,b,sense=lp.matrices();activity=np.asarray(A@x).reshape(-1)
    residual=activity-b
    violation=np.array([max(0,v) if s=='<' else abs(v) for v,s in zip(residual,sense)])/(1+np.abs(b))
    detail=dict(min_variable=float(min(x)),max_scaled_row_violation=float(max(violation,default=0)))
    return detail['min_variable']>=-1e-7 and detail['max_scaled_row_violation']<=1e-7,detail

class Master14(Master):
    def __init__(self,lp,method='phase1'):
        if method!='phase1':
            super().__init__(lp,method);return
        # x=0 satisfies every homogeneous row. Artificial variables are needed
        # only at nonzero RHS rows to give a feasible phase-I starting problem.
        from scipy import sparse as sp
        from coverage11 import _Highs,HighsLp,MatrixFormat
        self.method=method;self.lp=lp;A,b,sense=lp.matrices();self.base=A.shape[1]
        self.ub=np.array([i for i,s in enumerate(sense) if s=='<'],dtype=np.int32)
        eq=np.array([i for i,s in enumerate(sense) if s=='=' and b[i]!=0],dtype=np.int32)
        au=np.array([i for i in self.ub if b[i]!=0],dtype=np.int32)
        na=2*len(eq)+len(au)
        arts=sp.csc_matrix((np.r_[np.ones(len(eq)),-np.ones(len(eq)),-np.ones(len(au))],(np.r_[eq,eq,au],np.arange(na))),shape=(len(b),na))
        mat=sp.hstack([A,arts],format='csc');nv=mat.shape[1]
        hp=HighsLp();hp.num_col_=nv;hp.num_row_=len(b);hp.col_cost_=np.r_[np.zeros(self.base),np.ones(na)];hp.col_lower_=np.zeros(nv);hp.col_upper_=np.full(nv,np.inf)
        lower=b.copy();lower[self.ub]=-np.inf;hp.row_lower_=lower;hp.row_upper_=b
        hp.a_matrix_.format_=MatrixFormat.kColwise;hp.a_matrix_.start_=mat.indptr;hp.a_matrix_.index_=mat.indices;hp.a_matrix_.value_=mat.data
        self.solver=_Highs();self.solver.setOptionValue('threads',1);self.solver.setOptionValue('output_flag',False);self.solver.setOptionValue('solver','simplex');self.solver.passModel(hp);self.bump=None
        self.feasible_start=os.environ.get('XIV_PRIMAL_START')=='1'
        if self.feasible_start:
            from scipy.optimize._highspy._core import HighsBasis,HighsBasisStatus,HighsStatus
            basis=HighsBasis()
            cs=[HighsBasisStatus.kLower]*nv
            rs=[HighsBasisStatus.kBasic]*len(b)
            selected=[]
            for k,i in enumerate(eq):
                j=self.base+k+(0 if b[i]>0 else len(eq))
                cs[j]=HighsBasisStatus.kBasic;rs[i]=HighsBasisStatus.kUpper
                selected.append((int(i),int(j),abs(float(b[i]))))
            for k,i in enumerate(au):
                if b[i]<0:
                    j=self.base+2*len(eq)+k
                    cs[j]=HighsBasisStatus.kBasic;rs[i]=HighsBasisStatus.kUpper
                    selected.append((int(i),int(j),-float(b[i])))
            # This diagonal basis has original variables zero. Its selected
            # artificial values satisfy all nonzero equalities and negative
            # upper bounds; remaining basic row slacks satisfy the other rows.
            start=np.zeros(nv)
            for i,j,value in selected:start[j]=value
            activity=np.asarray(mat@start).reshape(-1)
            assert np.all(activity>=lower) and np.all(activity<=b)
            basis.col_status=cs;basis.row_status=rs
            assert self.solver.setBasis(basis)==HighsStatus.kOk
            self.solver.setOptionValue('simplex_strategy',4)
            print('PHASE14_EXPLICIT_FEASIBLE_BASIS',len(selected),'objective',sum(start[self.base:]),flush=True)
        print('PHASE14_NONZERO_RHS_ARTIFICIALS',na,'of',2*sum(s=='=' for s in sense)+len(self.ub),flush=True)
    def solve(self,lp,seconds):
        if self.method=='direct' and os.environ.get('XIV_PRESOLVE')=='on':
            self.solver.setOptionValue('presolve','on')
            return self._solve_once(lp,seconds)
        if self.solver.getOptionValue('solver')[1]!='simplex':return self._solve_once(lp,seconds)
        # HiGHS may discard the reduced-model basis/dual at a time limit.
        # Solve original rows so partial duals are available for exact checks.
        self.solver.setOptionValue('presolve','off')
        deadline=time.monotonic()+seconds
        while True:
            result=self._solve_once(lp,min(120,max(1,deadline-time.monotonic())))
            status,_,detail=result
            detail['simplex_presolve']='off'
            print('PHASE14_SIMPLEX_CHUNK',self.solver.getRunTime(),status,detail,flush=True)
            if status!='open' or stop or time.monotonic()>=deadline or 'kTimeLimit' not in detail.get('solver_status',''):
                return result
    def _solve_once(self,lp,seconds):
        status,data,detail=super().solve(lp,seconds)
        detail['phase1_mode']='nonzero_rhs_artificials' if self.method=='phase1' else 'none_direct'
        detail['explicit_feasible_start']=getattr(self,'feasible_start',False)
        if status=='numerical-feasible':
            valid,residual=numerical_feasibility(lp,data);detail['residual_check']=residual
            if valid:
                detail['degree_values']=[float(data[lp.col[f'm{j}']]) for j in range(1,8) if f'm{j}' in lp.col]
                selected=[dict(m=list(map(int,name[4:].split(','))),weight=float(data[j])) for j,name in enumerate(lp.names) if name.startswith('Q14:') and data[j]>1e-9]
                selected.sort(key=lambda x:-x['weight']);detail['positive_degree_vectors']=selected
            else:status,data='open',None;detail['reason']='numerical point failed original-row residual check'
        if status!='open' or detail.get('reason')=='exact repair failed':return status,data,detail
        # A candidate need not be solver-certified: only the exact final combination matters.
        dual=self.solver.getSolution().row_dual
        if len(dual)!=len(lp.rows) or not np.all(np.isfinite(dual)):
            return status,data,detail
        from plane10_cert import bump,combo
        if self.bump is None:self.bump=bump(lp,{}, {})[:2]
        by,bc=self.bump
        for scale in (10**7,10**9,10**11):
            y=[int(round(-v*scale)) for v in dual]
            for i,row in enumerate(lp.rows):
                if row[2]=='<':y[i]=max(0,y[i])
            cc,_=combo(lp,y)
            repair=max(0,max((-v+bc[j]-1)//bc[j] for j,v in enumerate(cc)))
            yy=[a+repair*b for a,b in zip(y,by)];cc,rhs=combo(lp,yy)
            if rhs<0 and min(cc)>=0:
                return 'LP-certified',yy,dict(detail,**verify(lp,yy),scale=scale,bump=repair,candidate_after_stop=True)
        return status,data,detail

def children(prefix,D):
    j=len(prefix)+1;n=46-sum(prefix);s=155-sum((i+1)*v for i,v in enumerate(prefix))
    if j==D:return [n] if n>=1 and s==D*n else []
    return [h for h in range(n+1) if (j+1)*(n-h)<=s-j*h<=D*(n-h)]
@lru_cache(None)
def raw(D):
    def rec(prefix):
        if len(prefix)==D:yield prefix;return
        for h in children(prefix,D):yield from rec(prefix+(h,))
    ans=tuple(rec(()))
    # Independent enumeration: choose the first D-2 counts, solve the final two linear equations.
    independent=set()
    def check(prefix,left):
        if len(prefix)==D-2:
            s=155-sum((i+1)*v for i,v in enumerate(prefix));last=s-(D-1)*left;penultimate=left-last
            if last>=1 and penultimate>=0:independent.add(prefix+(penultimate,last))
            return
        for v in range(left+1):check(prefix+(v,),left-v)
    check((),46);assert independent==set(ans) and len(ans)==len(set(ans))
    return ans

def nid(cfg,p):return cfg+':'+','.join(map(str,p))
def path(D):return Path(f'out/coverage14_D{D}.json')
@lru_cache(None)
def descendants(D,p):return tuple(m for m in raw(D) if m[:len(p)]==p)
def newnode(cfg,p):return dict(config=cfg,prefix=list(p),status='open',children=[])
def initialize(D):
    universe=raw(D)
    atomic(f'out/raw14_D{D}.json',dict(q=5,n=46,D=D,count=len(universe),vectors=universe))
    return dict(schema=14,q=5,n=46,D=D,family=FAMILY,raw_vectors=len(universe),total=2*len(universe),roots=['515:','617:'],nodes={nid(c,()):newnode(c,()) for c in ('515','617')})
def frontier(ledger):
    def rec(idx):
        n=ledger['nodes'][idx]
        if n['status']=='split':
            for child in n['children']:yield from rec(child)
        else:yield idx,n
    for r in ledger['roots']:yield from rec(r)
def split(ledger,idx):
    n=ledger['nodes'][idx];p=tuple(n['prefix']);D=ledger['D'];assert len(p)<D
    values=sorted({m[len(p)] for m in descendants(D,p)})
    n['children']=[nid(n['config'],p+(v,)) for v in values]
    for v,child in zip(values,n['children']):ledger['nodes'].setdefault(child,newnode(n['config'],p+(v,)))
    n['status']='split'
def refresh(ledger):
    c=Counter();opened=[]
    for idx,n in frontier(ledger):
        category='pruned_certified' if n['status'] in ('LP-certified','joint-LP-certified') else 'surviving' if n['status']=='exact-primal' else 'open'
        c[category]+=len(descendants(ledger['D'],tuple(n['prefix'])))
        if category=='open':opened.append(idx)
    assert sum(c.values())==ledger['total']
    ledger['counts']={k:c[k] for k in ('pruned_certified','surviving','open')};ledger['open_nodes']=opened
    ledger['coverage_fraction']=c['pruned_certified']/ledger['total'];ledger['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    atomic(path(ledger['D']),ledger)
def check_tree(ledger):
    assert (ledger['schema'],ledger['q'],ledger['n'],ledger['family'])==(14,5,46,FAMILY)
    D=ledger['D'];assert ledger['raw_vectors']==len(raw(D)) and ledger['total']==2*len(raw(D)) and ledger['roots']==['515:','617:']
    seen=set()
    for idx,n in ledger['nodes'].items():
        p=tuple(n['prefix']);assert idx==nid(n['config'],p) and descendants(D,p)
        assert n['config'] in ('515','617')
        if n['status'] in ('joint-LP-certified','exact-primal'):assert len(p)==D
        if n['children']:
            expected={nid(n['config'],p+(m[len(p)],)) for m in descendants(D,p)}
            assert set(n['children'])==expected and len(n['children'])==len(expected)
    for idx,n in frontier(ledger):
        for m in descendants(D,tuple(n['prefix'])):
            k=(n['config'],m);assert k not in seen;seen.add(k)
    assert len(seen)==ledger['total'];print('COVERAGE14_TREE_VERIFIED',D,len(seen),flush=True)
def fast_matrices(lp):
    # Same explicit rows, assembled directly as CSR to avoid large LIL insertion overhead.
    from scipy import sparse
    lengths=np.fromiter((len(row) for row,_,_,_ in lp.rows),dtype=np.int64,count=len(lp.rows))
    starts=np.r_[0,np.cumsum(lengths)];nnz=int(starts[-1])
    indices=np.fromiter((j for row,_,_,_ in lp.rows for j in row),dtype=np.int32,count=nnz)
    data=np.fromiter((v for row,_,_,_ in lp.rows for v in row.values()),dtype=np.float64,count=nnz)
    A=sparse.csr_matrix((data,indices,starts),shape=(len(lp.rows),len(lp.names)));A.sort_indices()
    b=np.fromiter((b for _,b,_,_ in lp.rows),dtype=np.float64,count=len(lp.rows))
    return A,b,[s for _,_,s,_ in lp.rows]

def build(D,cfg,prefix,nonlinear=False,cut_scope='all',nonlinear_scope='all',maximal_only=False,pair_bounds=False,degree_hull=False,triple_scope='none',compact_objectives=False,point_pairs=False,plane_budget=False,point_triples=False,star_bounds=False,rounded_budget=False,nonfano=False,extra=False,rich=False,point_mix=False,round_degree=False,point_subset=False):
    assert p6.SECTIONS and p6.DMH and p6.UNMARKED and p6.EXACT_K==4 and p6.F6_RMIN==5 and not os.environ.get('PC6_MFIX')
    v=json.load(open('out/validation14_base_D7.json' if D==7 else 'out/validation14_base.json'));assert v['status']=='ZERO_VIOLATIONS'
    if D==7:assert v['degree_cap']==7
    assert v['projection_table_sha256']==hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest()
    t=time.monotonic();candidates=descendants(D,tuple(prefix));assert candidates
    caps=tuple(max(m[j] for m in candidates) for j in range(D));rmax=5 if cfg=='515' else 6
    stored=json.load(open(f'out/types14_D{D}.json'));assert stored['D']==D and stored['q']==5 and stored['n']==46
    if maximal_only:
        from maximal14 import guard
        guard(D,cfg)
        assert {t[3] for t in stored['types'] if t[0]==rmax}==({15,16,17,18,21} if rmax==5 else {17,18,19,21,25})
    if plane_budget:
        from budget14 import guard as budget_guard,polynomial
        budget_guard()
    allowed_union=15 if cfg=='515' else 17
    types=[freeze(t) for t in stored['types'] if t[0]<=rmax and (not plane_budget or polynomial(t[0],t[1])<=541) and (not maximal_only or t[0]!=rmax or t[3]==allowed_union) and all(h<=c for h,c in zip(t[2],caps))];del stored
    original=p6.hist_types
    def filtered(q,d,n,*args,**kwargs):assert (q,d,n)==(5,D,46);return types
    p6.hist_types=filtered
    try:lp,info,types=p6.build(5,46,D,max_types=2000000)
    finally:p6.hist_types=original
    add_cuts(lp,info,types,D,cut_scope,compact_objectives)
    if nonlinear:
        from nonlinear14 import add_cuts as add_nonlinear
        add_nonlinear(lp,info,types,D,nonlinear_scope,compact_objectives)
    if triple_scope!='none':
        from triple14 import add_cuts as add_triple
        add_triple(lp,info,types,D,triple_scope,compact_objectives)
    if point_pairs:
        from point14 import add_cuts as add_point_pairs
        add_point_pairs(lp,info,types,D)
    if point_triples:
        from pointtriple14 import add_cuts as add_point_triples
        add_point_triples(lp,info,types,D)
    if star_bounds:
        assert maximal_only and plane_budget and point_triples
        from star14 import add_cuts as add_star_bounds
        add_star_bounds(lp,types,D,cfg)
    if rounded_budget:
        assert plane_budget
        from roundbudget14 import add_cuts as add_rounded_budget
        add_rounded_budget(lp,types)
    if nonfano:
        from pointnonfano14 import add_cuts as add_nonfano
        add_nonfano(lp,info,types,D)
    if point_subset:
        from pointsubset14 import add_cuts as add_point_subset
        add_point_subset(lp,info,types,D)
    if point_mix:
        assert star_bounds and maximal_only and plane_budget
        from pointmix14 import add_cuts as add_point_mix
        add_point_mix(lp,info,types,D)
    if rich:
        assert star_bounds and maximal_only and plane_budget
        from pointrich14 import add_cuts as add_rich
        add_rich(lp,info,types,D)
    if extra:
        assert maximal_only and cfg=='617'
        from extra14 import add_cuts as add_extra
        add_extra(lp,info,types,D,compact_objectives)
    if pair_bounds:
        from pairbounds14 import add_cuts as add_pair_bounds
        floors=tuple(min(m[j] for m in candidates) for j in range(D))
        add_pair_bounds(lp,D,floors,caps)
    if degree_hull:
        from degreehull14 import add_cuts as add_degree_hull
        add_degree_hull(lp,D,candidates)
    if round_degree:
        assert degree_hull
        from rounddegree14 import add_cuts as add_round_degree
        add_round_degree(lp,info,types,D,candidates)
    r,u=(5,15) if cfg=='515' else (6,17)
    lp.add([(lp.col[f'N{i}'],-1) for i,t in enumerate(types) if t[0]==r and t[3]==u],'<',-1,'XIV configuration')
    lp.add([(lp.col[f'N{i}'],1) for i,t in enumerate(types) if t[0]>r],'<',0,'XIV rmax')
    for j in range(D):
        lo=hi=prefix[j] if j<len(prefix) else None
        if j>=len(prefix):lo,hi=0,46
        lp.add([(lp.col[f'm{j+1}'],1)],'<',hi,f'XIV upper m{j+1}')
        lp.add([(lp.col[f'm{j+1}'],-1)],'<',-lo,f'XIV lower m{j+1}')
    print('BASE14_BUILT',D,cfg,prefix,len(lp.names),len(lp.rows),'seconds',time.monotonic()-t,flush=True)
    from types import MethodType
    lp.matrices=MethodType(fast_matrices,lp)
    if os.environ.get('XIV_FINITE_COLUMNS')=='1':
        canonical={f'N{i}':('N',t) for i,t in enumerate(types)}
        canonical.update({f'G{g}':('G',a,tuple(info['_dms'][di])) for g,(a,di) in enumerate(info['_Y'])})
        canonical.update({f'Z{z}':('Z',tuple(info['_dms'][di]),tuple(sorted(
            (freeze(info['_items'][ii]) for ii in info['_wt'][wi]),key=key)))
            for z,(di,wi) in enumerate(info['_Z'])})
        canonical.update({name:('Q',tuple(map(int,name[4:].split(',')))) for name in lp.names if name.startswith('Q14:')})
        lp.canonical14=canonical;lp.model14=dict(q=5,D=D,config=cfg,prefix=list(prefix))
    return lp,info,types

def solve_node(ledger,idx,seconds,solver,nonlinear=False,cut_scope='all',nonlinear_scope='all',maximal_only=False,pair_bounds=False,degree_hull=False,triple_scope='none',compact_objectives=False,point_pairs=False,plane_budget=False,point_triples=False,star_bounds=False,rounded_budget=False,nonfano=False,extra=False,rich=False,point_mix=False,round_degree=False,point_subset=False):
    if nonlinear:ledger['nonlinear_used']=True;refresh(ledger)
    if triple_scope!='none':ledger['triple_used']=True;refresh(ledger)
    if point_pairs:ledger['point_pairs_used']=True;refresh(ledger)
    if point_triples:ledger['point_triples_used']=True;refresh(ledger)
    if nonfano:ledger['nonfano_used']=True;refresh(ledger)
    if extra:ledger['extra_used']=True;refresh(ledger)
    if rich:ledger['rich_used']=True;refresh(ledger)
    if point_mix:ledger['point_mix_used']=True;refresh(ledger)
    if round_degree:ledger['round_degree_used']=True;refresh(ledger)
    if point_subset:ledger['point_subset_used']=True;refresh(ledger)
    n=ledger['nodes'][idx];D=ledger['D'];p=tuple(n['prefix']);lp,info,types=build(D,n['config'],p,nonlinear,cut_scope,nonlinear_scope,maximal_only,pair_bounds,degree_hull,triple_scope,compact_objectives,point_pairs,plane_budget,point_triples,star_bounds,rounded_budget,nonfano,extra,rich,point_mix,round_degree,point_subset)
    if stop:return False
    del info,types
    # A new Master is built for this exact LP; no inter-node bound changes or reused bases.
    numerical_method=os.environ.get('XIV_METHOD','phase1')
    assert numerical_method in ('phase1','direct')
    master=Master14(lp,numerical_method);master.solver.setOptionValue('solver',solver)
    if solver=='ipm':master.solver.setOptionValue('run_crossover','off')
    master.solver.setOptionValue('output_flag',True)
    master.solver.setOptionValue('log_to_console',False)
    master.solver.setOptionValue('log_file',f"out/solver14_D{D}_{n['config']}_"+('_'.join(map(str,p)) or 'root')+('_nonlinear' if nonlinear else '')+'_'+cut_scope+'_'+nonlinear_scope+('_maximal' if maximal_only else '')+('_products' if pair_bounds else '')+('_hull' if degree_hull else '')+CERT_SUFFIX+'.log')
    print('NODE14',D,idx,solver,seconds,flush=True);t=time.monotonic()
    status,data,detail=master.solve(lp,min(seconds,max(1,CAMPAIGN_STOP-time.time())))
    if status=='open':
        candidate=master.solver.getSolution()
        diagnostic=f"out/dual14_D{D}_{n['config']}_"+('_'.join(map(str,p)) or 'root')+f'_{time.time_ns()}.npz'
        np.savez_compressed(diagnostic,y=np.asarray(candidate.row_dual))
        detail['candidate_dual']=diagnostic
        if len(candidate.col_value)>=len(lp.names):
            detail['candidate_degree_values']=[float(candidate.col_value[lp.col[f'm{j}']]) for j in range(1,D+1)]
    detail.update(solver=solver,point_subset=point_subset,round_degree=round_degree,point_mix=point_mix,rich=rich,extra=extra,nonfano=nonfano,rounded_budget=rounded_budget,star_bounds=star_bounds,point_triples=point_triples,plane_budget=plane_budget,point_pairs=point_pairs,compact_objectives=compact_objectives,triple_scope=triple_scope,degree_hull=degree_hull,pair_bounds=pair_bounds,maximal_only=maximal_only,nonlinear=nonlinear,cut_scope=cut_scope,nonlinear_scope=nonlinear_scope,seconds=time.monotonic()-t,allowance=seconds,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    n.setdefault('probes',[]).append(dict(detail,status=status));n['prefix_probe_done']=True
    if status=='LP-certified':
        fn=f"out/cert14_D{D}_{n['config']}_"+('_'.join(map(str,p)) or 'root')+CERT_SUFFIX+'.json'
        from maximal14 import guard
        reduction_digest=guard(D,n['config']) if maximal_only else None
        from budget14 import guard as budget_guard
        budget_digest=budget_guard() if plane_budget else None
        from star14 import guard as star_guard
        star_digest=star_guard(D,n['config']) if star_bounds else None
        from roundbudget14 import guard as rounded_guard
        rounded_digest=rounded_guard() if rounded_budget else None
        atomic(fn,dict(point_subset=point_subset,round_degree=round_degree,point_mix=point_mix,rich=rich,extra=extra,nonfano=nonfano,rounded_budget=rounded_budget,rounded_budget_sha256=rounded_digest,star_bounds=star_bounds,star_sha256=star_digest,D=D,q=5,n=46,config=n['config'],prefix=p,family=FAMILY,point_triples=point_triples,plane_budget=plane_budget,budget_sha256=budget_digest,point_pairs=point_pairs,compact_objectives=compact_objectives,triple_scope=triple_scope,degree_hull=degree_hull,pair_bounds=pair_bounds,maximal_only=maximal_only,reduction_sha256=reduction_digest,nonlinear=nonlinear,cut_scope=cut_scope,nonlinear_scope=nonlinear_scope,projection_table_sha256=hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest(),multipliers=data,detail=detail))
        verify(lp,json.load(open(fn))['multipliers']);n.update(status=status,certificate=fn)
    elif len(p)==D:
        n['status']='open';n['reason']='joint IX/X pending' if status=='numerical-feasible' else 'unfinished prefix solve'
        if status=='numerical-feasible':
            n['joint_pending']=True;np.savez_compressed(f"out/prefix14_D{D}_{n['config']}_"+'_'.join(map(str,p))+'.npz',x=data)
    else:
        n['reason']='numerical prefix feasible' if status=='numerical-feasible' else 'unfinished prefix solve'
        split(ledger,idx)
    refresh(ledger);print('RESULT14',D,idx,status,detail,ledger['counts'],flush=True)
    return status=='numerical-feasible' and len(p)==D

def release_price_memory():
    # Fresh prices have already been saved and checked. Certificate construction
    # reads those tables; temporary optimization memoization is not proof data.
    import sys,gc
    for name in ('lift14','nonlinear14','triple14','point14','pointtriple14','pointnonfano14','pointrich14','pointmix14','extra14','rounddegree14','pointsubset14'):
        module=sys.modules.get(name)
        if module:
            for value in list(vars(module).values()):
                clear=getattr(value,'cache_clear',None)
                if clear:clear()
    gc.collect()

def audit(ledger,fresh=False):
    replayed=[];domain_replayed={}
    D=ledger['D'];types=p6.hist_types(5,D,46);stored=json.load(open(f'out/types14_D{D}.json'))
    assert freeze(stored['types'])==tuple(types)
    if fresh:table(D,types,True);release_price_memory()
    else:
        from lift14 import price,objectives,CORE
        cache=json.load(open(f'out/lift14_D{D}_prices.json'));obj=objectives(D)
        for t in types:
            if t[0]>=5:assert cache['prices'][key(t)]==price(t,D,obj if t[:2] in CORE else obj[:2*D])
    if ledger.get('nonlinear_used'):
        from nonlinear14 import table as nonlinear_table
        nonlinear_table(D,True);release_price_memory()
    if ledger.get('triple_used'):
        from triple14 import table as triple_table
        triple_table(D,True);release_price_memory()
    if ledger.get('point_pairs_used'):
        from point14 import table as point_table
        point_table(D,True);release_price_memory()
    if ledger.get('point_triples_used'):
        from pointtriple14 import table as point_triple_table
        point_triple_table(D,True);release_price_memory()
    if ledger.get('nonfano_used'):
        from pointnonfano14 import table as nonfano_table
        nonfano_table(D,True);release_price_memory()
    if ledger.get('point_subset_used'):
        from pointsubset14 import table as point_subset_table
        point_subset_table(D,True);release_price_memory()
    if ledger.get('round_degree_used'):
        from rounddegree14 import table as round_degree_table
        round_degree_table(D,True);release_price_memory()
    if ledger.get('point_mix_used'):
        from pointmix14 import table as point_mix_table
        point_mix_table(D,True);release_price_memory()
    if ledger.get('rich_used'):
        from pointrich14 import table as rich_table
        rich_table(D,True);release_price_memory()
    if ledger.get('extra_used'):
        from extra14 import table as extra_table
        extra_table(D,True);release_price_memory()
    del types,stored
    release_price_memory()
    for idx,n in frontier(ledger):
        if n['status']=='LP-certified':
            r=json.load(open(n['certificate']))
            if r.get('domain'):
                from domain14 import membership,replay_record
                assert (r['D'],r['config'])==(D,n['config'])
                cases=membership(D,n['prefix'],r['domain'])
                sha=hashlib.sha256(Path(n['certificate']).read_bytes()).hexdigest();cachekey=(n['certificate'],sha)
                if cachekey not in domain_replayed:
                    domain_replayed[cachekey]=replay_record(r)
                    import gc;gc.collect()
                detail=domain_replayed[cachekey]
                replayed.append(dict(node=idx,certificate=n['certificate'],sha256=sha,cases=cases,domain=r['domain'],**detail))
                print('EXACT_DOMAIN14_PREFIX_COVERED',idx,cases,detail,flush=True)
                continue
            assert (r['D'],r['q'],r['n'],r['config'],r['prefix'],r['family'])==(D,5,46,n['config'],n['prefix'],FAMILY)
            assert r['projection_table_sha256']==hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest()
            if r.get('rounded_budget'):
                from roundbudget14 import guard as rounded_guard
                assert rounded_guard()==r['rounded_budget_sha256']
            if r.get('star_bounds'):
                from star14 import guard as star_guard
                assert star_guard(D,n['config'])==r['star_sha256']
            if r.get('plane_budget'):
                from budget14 import guard as budget_guard
                assert budget_guard()==r['budget_sha256']
            if r.get('maximal_only'):
                from maximal14 import guard
                assert guard(D,n['config'])==r['reduction_sha256']
            lp,_,_=build(D,n['config'],n['prefix'],r.get('nonlinear',False),r.get('cut_scope','all'),r.get('nonlinear_scope','all'),r.get('maximal_only',False),r.get('pair_bounds',False),r.get('degree_hull',False),r.get('triple_scope','none'),r.get('compact_objectives',False),r.get('point_pairs',False),r.get('plane_budget',False),r.get('point_triples',False),r.get('star_bounds',False),r.get('rounded_budget',False),r.get('nonfano',False),r.get('extra',False),r.get('rich',False),r.get('point_mix',False),r.get('round_degree',False),r.get('point_subset',False));detail=verify(lp,r['multipliers']);print('EXACT_PREFIX14_REPLAYED',D,idx,detail,flush=True)
            replayed.append(dict(node=idx,certificate=n['certificate'],sha256=hashlib.sha256(Path(n['certificate']).read_bytes()).hexdigest(),rows=len(lp.rows),columns=len(lp.names),cases=len(descendants(D,tuple(n['prefix']))),**detail));del lp
            import gc;gc.collect()
        elif n['status'] in ('joint-LP-certified','exact-primal'):
            from joint14 import replay_node
            replay_node(D,n)
    check_tree(ledger);refresh(ledger)
    atomic(f'out/replay14_D{D}_final.json',dict(status='EXACT_COVERAGE14_REPLAYED',D=D,fresh_prices=fresh,counts=ledger['counts'],tree_cases=ledger['total'],certificates=replayed,ledger_sha256=hashlib.sha256(path(D).read_bytes()).hexdigest(),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    print('EXACT_COVERAGE14_REPLAYED',D,ledger['counts'],flush=True)

def main():
    global stop
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True,choices=[6,7]);ap.add_argument('--init',action='store_true');ap.add_argument('--replay',action='store_true');ap.add_argument('--fresh-prices',action='store_true');ap.add_argument('--nodes',type=int,default=50);ap.add_argument('--seconds',type=int,default=900);ap.add_argument('--solver',choices=['ipm','simplex'],default='ipm');ap.add_argument('--config',choices=['515','617']);ap.add_argument('--roots-only',action='store_true');ap.add_argument('--target');ap.add_argument('--nonlinear',action='store_true');ap.add_argument('--cut-scope',choices=['all','core','pair','coordinate'],default='all');ap.add_argument('--nonlinear-scope',choices=['all','pair','singletons'],default='all');ap.add_argument('--maximal-only',action='store_true');ap.add_argument('--pair-bounds',action='store_true');ap.add_argument('--degree-hull',action='store_true');ap.add_argument('--triple-scope',choices=['none','all','pair','singletons'],default='none');ap.add_argument('--compact-objectives',action='store_true');ap.add_argument('--point-pairs',action='store_true');ap.add_argument('--plane-budget',action='store_true');ap.add_argument('--point-triples',action='store_true');ap.add_argument('--star-bounds',action='store_true');ap.add_argument('--rounded-budget',action='store_true');ap.add_argument('--nonfano',action='store_true');ap.add_argument('--extra',action='store_true');ap.add_argument('--rich',action='store_true');ap.add_argument('--point-mix',action='store_true');ap.add_argument('--round-degree',action='store_true');ap.add_argument('--point-subset',action='store_true');a=ap.parse_args()
    def handler(*_):
        global stop
        stop=True
    signal.signal(signal.SIGTERM,handler);signal.signal(signal.SIGINT,handler)
    lock=open(f'out/coverage14_D{a.D}.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    ledger=json.load(path(a.D).open()) if path(a.D).exists() else initialize(a.D)
    refresh(ledger);check_tree(ledger)
    if a.init:return
    if a.replay:audit(ledger,a.fresh_prices);return
    attempts=0;tried=set()
    def attempt(idx):
        nonlocal attempts
        if stop or attempts>=a.nodes or time.time()+a.seconds+180>=CAMPAIGN_STOP:return True
        n=ledger['nodes'][idx]
        if n['status'] in ('LP-certified','joint-LP-certified','exact-primal'):return False
        if idx in tried:return False
        tried.add(idx)
        attempts+=1;return solve_node(ledger,idx,a.seconds,a.solver,a.nonlinear,a.cut_scope,a.nonlinear_scope,a.maximal_only,a.pair_bounds,a.degree_hull,a.triple_scope,a.compact_objectives,a.point_pairs,a.plane_budget,a.point_triples,a.star_bounds,a.rounded_budget,a.nonfano,a.extra,a.rich,a.point_mix,a.round_degree,a.point_subset)
    if a.target:attempt(a.target)
    else:
        # Every configuration root is tried before any descendant is considered.
        for idx in ledger['roots']:
            n=ledger['nodes'][idx]
            if a.config and n['config']!=a.config:continue
            if not n.get('prefix_probe_done') and attempt(idx):break
        if not a.roots_only:
            while not stop and attempts<a.nodes:
                todo=[idx for idx,n in frontier(ledger) if n['status']=='open' and idx not in tried and (not a.config or n['config']==a.config) and not n.get('joint_pending')]
                if not todo or attempt(todo[0]):break
    check_tree(ledger);refresh(ledger);print('COVERAGE14_CHECKPOINT',a.D,ledger['counts'],'attempts',attempts,'stop',stop,flush=True)
if __name__=='__main__':main()
