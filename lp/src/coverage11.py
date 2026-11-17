"""Resumable D=5 prefix coverage. Only exact replayed certificates close nodes."""
import argparse,copy,json,time,hashlib,os,signal
for _v in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[_v]="1"
from pathlib import Path
from collections import Counter
import numpy as np
from scipy import sparse as sp
from scipy.optimize._highspy._core import _Highs,HighsLp,MatrixFormat,HighsModelStatus
import planecount6 as p6
from degvec_bb import children
from lift11 import atomic,add_cuts,key,prices
from plane10 import geometries
from plane10_cert import combo,bump
LEDGER=Path('out/coverage11.json')
HARD_STOP=None # assigned from the explicit UTC deadline in main
stop=False

def nid(cfg,p):return cfg+':'+','.join(map(str,p))
def leaves(p=()):
    if len(p)==5:yield p;return
    for h in children(p,46,155,5):yield from leaves(p+(h,))
RAW=list(leaves());assert len(RAW)==2599

def initialize():
    nodes={}
    for cfg in ('515','617'):
        for m in RAW:
            for d in range(6):
                p=m[:d];idx=nid(cfg,p)
                if idx not in nodes:nodes[idx]=dict(config=cfg,prefix=list(p),status='open',children=[])
                if d:
                    parent=nodes[nid(cfg,m[:d-1])]
                    if idx not in parent['children']:parent['children'].append(idx)
    return dict(schema=1,q=5,n=46,D=5,raw_vectors=2599,total=5198,roots=['515:','617:'],nodes=nodes)

def frontier(ledger):
    def rec(idx):
        node=ledger['nodes'][idx]
        if node['status']=='split':
            for child in node['children']:yield from rec(child)
        else:yield idx,node
    for root in ledger['roots']:yield from rec(root)

def refresh(ledger):
    counts=Counter();opened=[]
    for idx,node in frontier(ledger):
        count=sum(m[:len(node['prefix'])]==tuple(node['prefix']) for m in RAW)
        st=node['status'];cat='pruned_certified' if st in ('LP-certified','inherited-LP-certified') else 'surviving' if st=='exact-primal' else 'open'
        counts[cat]+=count
        if cat=='open':opened.append(idx)
    ledger['counts']={k:counts[k] for k in ('pruned_certified','surviving','open')}
    assert sum(counts.values())==5198
    ledger['open_nodes']=opened;ledger['coverage_fraction']=counts['pruned_certified']/5198
    ledger['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    atomic(LEDGER,ledger)

def build(scope="core",config=None,prefix=()):
    assert p6.SECTIONS and p6.DMH and p6.UNMARKED and p6.EXACT_K==4 and p6.F6_RMIN==5
    assert not os.environ.get('PC6_MFIX')
    v=json.load(open('out/lift11_validation.json'));assert v['status']=='ZERO_VIOLATIONS'
    t=time.monotonic()
    original=p6.hist_types
    if config is not None:
        candidates=[m for m in RAW if m[:len(prefix)]==tuple(prefix)]
        assert candidates
        cap=tuple(max(m[j] for m in candidates) for j in range(5))
        rmax=5 if config=='515' else 6
        def filtered(*args,**kwargs):
            return [t for t in original(*args,**kwargs) if t[0]<=rmax and all(h<=c for h,c in zip(t[2],cap))]
        p6.hist_types=filtered
    try:lp,info,types=p6.build(5,46,5)
    finally:p6.hist_types=original
    add_cuts(lp,info,types,'core' if scope=='pair' else scope)
    if scope=='pair':
        from lift13 import add_cuts as pair_cuts
        pair_cuts(lp,info,types)
    print('BASE_BUILT',len(lp.names),len(lp.rows),'seconds',time.monotonic()-t,flush=True)
    return lp,info,types

def node_lp(base,types,cfg,prefix):
    lp=copy.copy(base);lp.rows=list(base.rows);r,u=(5,15) if cfg=='515' else (6,17)
    lp.add([(lp.col[f'N{i}'],-1) for i,t in enumerate(types) if t[0]==r and t[3]==u],'<',-1,'XI config')
    lp.add([(lp.col[f'N{i}'],1) for i,t in enumerate(types) if t[0]>r],'<',0,'XI rmax')
    for j in range(5):
        lo=hi=prefix[j] if j<len(prefix) else None
        if j>=len(prefix):lo,hi=0,46
        lp.add([(lp.col[f'm{j+1}'],1)],'<',hi,f'XI upper m{j+1}')
        lp.add([(lp.col[f'm{j+1}'],-1)],'<',-lo,f'XI lower m{j+1}')
    return lp

def verify(lp,y):
    assert len(y)==len(lp.rows) and all(type(v)==int for v in y)
    assert all(v>=0 for v,row in zip(y,lp.rows) if row[2]=='<')
    c,b=combo(lp,y);assert min(c)>=0 and b<0
    return dict(rhs=b,min_column=min(c))

class Master:
    def __init__(self,lp,method="direct"):
        self.method=method;self.lp=lp;A,b,sense=lp.matrices();self.ub=np.array([i for i,s in enumerate(sense) if s=='<']);eq=np.array([i for i,s in enumerate(sense) if s=='=']);self.base=A.shape[1]
        arts=sp.csc_matrix((np.r_[np.ones(len(eq)),-np.ones(len(eq)),-np.ones(len(self.ub))],(np.r_[eq,eq,self.ub],np.arange(2*len(eq)+len(self.ub)))),shape=(len(b),2*len(eq)+len(self.ub)))
        mat=sp.hstack([A,arts],format='csc') if method=='phase1' else A.tocsc();nv=mat.shape[1]
        hp=HighsLp();hp.num_col_=nv;hp.num_row_=len(b);hp.col_cost_=np.r_[np.zeros(self.base),np.ones(arts.shape[1])] if method=='phase1' else np.zeros(nv);hp.col_lower_=np.zeros(nv);hp.col_upper_=np.full(nv,np.inf)
        lower=b.copy();lower[self.ub]=-np.inf;hp.row_lower_=lower;hp.row_upper_=b
        hp.a_matrix_.format_=MatrixFormat.kColwise;hp.a_matrix_.start_=mat.indptr;hp.a_matrix_.index_=mat.indices;hp.a_matrix_.value_=mat.data
        self.solver=_Highs();self.solver.setOptionValue('threads',1);self.solver.setOptionValue('output_flag',False);self.solver.setOptionValue('solver','simplex');self.solver.passModel(hp)
        self.bump=None
    def solve(self,lp,seconds):
        for i in range(len(lp.rows)-12,len(lp.rows)):
            row=lp.rows[i];self.solver.changeRowBounds(i,-np.inf,row[1])
        self.solver.setOptionValue('time_limit',self.solver.getRunTime()+seconds)
        self.solver.run();st=self.solver.getModelStatus()
        if self.method=='direct':
            if st==HighsModelStatus.kOptimal:return 'numerical-feasible',np.array(self.solver.getSolution().col_value[:self.base]),dict(method='direct')
            if st!=HighsModelStatus.kInfeasible:return 'open',None,dict(solver_status=str(st),method='direct')
            _,exists,dual=self.solver.getDualRay()
            if not exists:return 'open',None,dict(reason='no numerical dual ray',solver_status=str(st),method='direct')
            detail=dict(method='direct')
        else:
            if st!=HighsModelStatus.kOptimal:return 'open',None,dict(solver_status=str(st))
            sol=self.solver.getSolution();obj=self.solver.getObjectiveValue();detail=dict(phase1=obj)
            if obj<1e-7:return 'numerical-feasible',np.array(sol.col_value[:self.base]),detail
            dual=sol.row_dual
        if self.bump is None:self.bump=bump(lp,{}, {})[:2]
        by,bc=self.bump
        for scale in (10**7,10**9,10**11):
            y=[int(round(-v*scale)) for v in dual]
            for i,row in enumerate(lp.rows):
                if row[2]=='<':y[i]=max(0,y[i])
            cc,rhs=combo(lp,y);T=max(0,max((-v+bc[j]-1)//bc[j] for j,v in enumerate(cc)))
            yy=[a+T*b for a,b in zip(y,by)];cc,rhs=combo(lp,yy)
            if rhs<0 and min(cc)>=0:
                detail.update(verify(lp,yy));detail.update(scale=scale,bump=T)
                return 'LP-certified',yy,detail
        return 'open',None,dict(detail,reason='exact repair failed')

def check_tree(ledger):
    independent=set()
    for a in range(47):
        for b in range(47-a):
            for c in range(47-a-b):
                remaining=46-a-b-c;d5=155-a-2*b-3*c-4*remaining;d4=remaining-d5
                if d4>=0 and d5>=1:independent.add((a,b,c,d4,d5))
    assert independent==set(RAW)
    # Verify exhaustive disjoint prefix coverage, independently from the displayed count.
    covered={}
    assert ledger['roots']==['515:','617:']
    for idx,node in ledger['nodes'].items():
        assert idx==nid(node['config'],node['prefix'])
        p=tuple(node['prefix'])
        expected={nid(node['config'],m[:len(p)+1]) for m in RAW if len(p)<5 and m[:len(p)]==p}
        assert set(node['children'])==expected and len(node['children'])==len(expected)
    for idx,node in frontier(ledger):
        for m in RAW:
            if m[:len(node['prefix'])]==tuple(node['prefix']):
                k=(node['config'],m);assert k not in covered;covered[k]=idx
    assert len(covered)==5198
    print('COVERAGE11_TREE_VERIFIED',len(covered),flush=True)

def fresh_prices(types):
    """Recompute from geometry with no disk-price inputs or prior memoized prices."""
    from lift11 import price_profile
    price_profile.cache_clear()
    geos=geometries()
    fresh=dict(geometry_sha256=hashlib.sha256(key(geos).encode()).hexdigest(),prices={})
    for t in types:
        if t[0]>=5:fresh['prices'][key(t)]=prices(t,geos)
    # Old values are read only AFTER the fresh computation, for comparison.
    old=json.load(open('out/lift11_prices.json'))
    assert old['geometry_sha256']==fresh['geometry_sha256']
    assert all(old['prices'].get(k)==v for k,v in fresh['prices'].items())
    session=os.environ.get('COVERAGE_SESSION','12')
    atomic(f'out/lift{session}_fresh_prices.json',fresh)
    atomic('out/lift11_prices.json',fresh)
    atomic(f'out/lift{session}_fresh_check.json',dict(types=len(fresh['prices']),mismatches=0,
           price_inputs='none; integer assignment pricing from geometry',geometry_sha256=fresh['geometry_sha256']))
    print('FRESH_PRICES_VERIFIED',len(fresh['prices']),'mismatches 0',flush=True)

def audit(ledger,base,types,history=False):
    cache=json.load(open('out/lift11_prices.json'));geos=geometries();num=0
    for t in types:
        if t[0]>=5:assert cache['prices'][key(t)]==prices(t,geos);num+=t[0]>=5
    if Path('out/lift13_pair_prices.json').exists():
        from lift13 import admitted,prices as pair_prices
        pair_cache=json.load(open('out/lift13_pair_prices.json'))
        for t in types:
            if admitted(t):assert pair_cache['prices'][key(t)]==pair_prices(t)
        print('PAIR13_PRICES_REPLAYED',len(pair_cache['prices']),flush=True)
    for idx,node in (ledger['nodes'].items() if history else frontier(ledger)):
        if node['status']=='LP-certified':
            rec=json.load(open(node['certificate']));assert rec['prefix']==node['prefix'] and rec['config']==node['config']
            if rec.get('model_kind')=='prefix-filtered':
                nb,_,nt=build(rec.get('cut_scope',ledger['cut_scope']),node['config'],node['prefix'])
            else:nb,_,nt=build(rec.get('cut_scope',ledger['cut_scope']))
            detail=verify(node_lp(nb,nt,node['config'],node['prefix']),rec['multipliers'])
            print('EXACT_PREFIX_REPLAYED',idx,detail,flush=True)
    check_tree(ledger)
    refresh(ledger);print('EXACT_COVERAGE11_REPLAYED',ledger['counts'],'priced_types',num,flush=True)

def main():
    import datetime
    global HARD_STOP
    HARD_STOP=float(os.environ.get('COVERAGE_DEADLINE',datetime.datetime(2026,11,4,11,15,tzinfo=datetime.timezone.utc).timestamp()))
    ap=argparse.ArgumentParser();ap.add_argument('--stop-on-leaf',action='store_true');ap.add_argument('--nodes',type=int,default=40);ap.add_argument('--seconds',type=float,default=180);ap.add_argument('--replay',action='store_true');ap.add_argument('--replay-history',action='store_true');ap.add_argument('--fresh-prices',action='store_true');ap.add_argument('--config',choices=['515','617']);ap.add_argument('--target',help='retry one prefix, e.g. 617:0; combine with --nodes 1');ap.add_argument('--import-x',action='store_true');ap.add_argument('--import-roots',action='store_true');ap.add_argument('--scope',choices=['core','all','pair'],default='core');ap.add_argument('--method',choices=['direct','phase1'],default='phase1');a=ap.parse_args()
    def handler(*_):
        global stop
        stop=True
    signal.signal(signal.SIGTERM,handler);signal.signal(signal.SIGINT,handler)
    import fcntl
    lock=open('out/coverage11.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    ledger=json.load(LEDGER.open()) if LEDGER.exists() else initialize();refresh(ledger)
    if a.import_x:
        from plane10_cert import run
        for cfg,m in [('515',(0,4,24,15,3)),('617',(0,3,27,12,4))]:
            run(cfg,replay=True)
            # Record this proven leaf separately until a prefix is opened around it.
            node=ledger['nodes'][nid(cfg,m)];node.update(status='inherited-LP-certified',certificate=f'out/plane10_{cfg}_farkas.json')
            for d in range(5):
                ancestor=ledger['nodes'][nid(cfg,m[:d])]
                if ancestor['status']!='LP-certified':ancestor['status']='split'
        refresh(ledger)
    if 'cut_scope' not in ledger:
        assert not any(n['status']=='LP-certified' for n in ledger['nodes'].values())
        ledger['cut_scope']=a.scope;refresh(ledger)
    ledger['model_kind']='prefix-filtered';refresh(ledger)
    if a.import_roots:
        assert a.scope=='core'
        for cfg in ('515','617'):
            result_path=Path('out/coverage11_root_'+cfg+'_result.json')
            if result_path.exists():ledger['nodes'][cfg+':']['root_probe']=json.load(result_path.open())
            path=Path('out/coverage11_cert_'+cfg+'_root.json')
            if not path.exists():continue
            rec=json.load(path.open());assert rec['config']==cfg and rec['prefix']==[] and rec['model_kind']=='prefix-filtered'
            nb,_,nt=build(a.scope,cfg,[]);verify(node_lp(nb,nt,cfg,[]),rec['multipliers'])
            ledger['nodes'][cfg+':'].update(status='LP-certified',certificate=str(path),prefix_probe_done=True,reason='exact root certificate replayed')
            refresh(ledger);print('ROOT_IMPORTED',cfg,ledger['counts'],flush=True)
    base=None;info=None;types=p6.hist_types(5,5,46)
    if a.fresh_prices:
        assert a.replay or a.replay_history, '--fresh-prices requires a replay'
        fresh_prices(types)
        if Path('out/lift13_pair_prices.json').exists():
            from lift13 import table
            table(types,fresh=True)
    if a.replay or a.replay_history:
        # Imported X certificates are rechecked as well; never inferred from log text.
        for n in ledger['nodes'].values():
            if n['status'] in ('inherited-LP-certified','exact-primal'):
                assert len(n['prefix'])==5
                tag=n.get('joint_tag',n['config'])
                model=json.load(open('out/pencil9_'+tag+'.json'))
                assert model['m']==n['prefix']
                assert model['config']==([5,15] if n['config']=='515' else [6,17])
        if any(n['status']=='inherited-LP-certified' for n in ledger['nodes'].values()):
            from plane10_cert import run
            for tag in sorted({n.get('joint_tag',n['config']) for n in ledger['nodes'].values() if n['status']=='inherited-LP-certified'}):run(tag,replay=True)
        from plane10_check import check
        for n in ledger['nodes'].values():
            if n['status']=='exact-primal':check(n['joint_tag'],replay=True)
        audit(ledger,base,types,a.replay_history);return
    # Internal split nodes without a prefix probe are still eligible for coarse pruning.
    masters={};attempts=0
    def visit(idx):
        nonlocal attempts
        global stop
        node=ledger['nodes'][idx]
        if stop or attempts>=a.nodes or time.time()+a.seconds+30>=HARD_STOP:return
        if a.config and node['config']!=a.config:return
        if node['status'] in ('LP-certified','inherited-LP-certified','exact-primal'):return
        if not node.get('prefix_probe_done'):
            cfg=node['config'];masters.clear()
            nb,_,nt=build(a.scope,cfg,node['prefix']);lp=node_lp(nb,nt,cfg,node['prefix'])
            if stop:return
            masters[cfg]=Master(lp,a.method)
            attempts+=1;print('NODE',idx,'attempt',attempts,flush=True)
            st,data,detail=masters[cfg].solve(lp,a.seconds);detail['cut_scope']=a.scope;node['last_probe']=detail
            node.setdefault('probes',[]).append(dict(detail,scope=a.scope,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
            if st=='LP-certified':
                path=f'out/coverage11_cert_{cfg}_'+('_'.join(map(str,node['prefix'])) or 'root')+'.json'
                atomic(path,dict(config=cfg,prefix=node['prefix'],model_kind='prefix-filtered',cut_scope=a.scope,multipliers=data,detail=detail))
                # Reload from disk; only this exact replay permits pruning.
                verify(lp,json.load(open(path))['multipliers']);node.update(status=st,certificate=path,prefix_probe_done=True,reason='exact certificate replayed')
            elif st=='numerical-feasible':
                node['prefix_probe_done']=True
                if node['children']:node['status']='split'
                else:
                    node.update(status='open',reason='prefix LP feasible; joint IX/X pricing pending')
                    if a.stop_on_leaf:stop=True
                    np.savez_compressed(f'out/coverage11_prefix_{cfg}_'+ '_'.join(map(str,node['prefix']))+'.npz',x=data)
            else:
                node['reason']='prefix solve unfinished; no exclusion'
                if node['children']:
                    # Splitting needs no feasibility claim: children partition the parent.
                    node.update(status='split',prefix_probe_done=True,split_after_unfinished=True)
            refresh(ledger);print('RESULT',idx,st,detail,ledger['counts'],flush=True)
            if st=='open' and node['status']!='split':return
        if node['status']=='split':
            for child in node['children']:visit(child)
    if a.target:
        node=ledger['nodes'][a.target];p=node['prefix'];cfg=node['config']
        assert node['status'] not in ('LP-certified','inherited-LP-certified','exact-primal')
        for d in range(len(p)):
            assert ledger['nodes'][nid(cfg,p[:d])]['status']=='split', 'target lies beneath a closed or unsplit ancestor'
        node['prefix_probe_done']=False
        visit(a.target)
    else:
        for root in ledger['roots']:visit(root)
    check_tree(ledger);refresh(ledger);print('COVERAGE11_CHECKPOINT',ledger['counts'],'attempts',attempts,'stop',stop,flush=True)
if __name__=='__main__':main()
