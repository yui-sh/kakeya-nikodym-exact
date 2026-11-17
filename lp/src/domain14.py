"""Exact range certificates for unions of explicit m1-prefix branches.

Every raw degree vector in m1>=bound is retained, regardless of prior pruning.
An imported certificate is shared by explicit tree nodes whose entire prefix
lies in that domain. No range probe ever claims the unrestricted root.
"""
import argparse,json,time,hashlib,fcntl,os,signal,re
from pathlib import Path
import coverage14 as c
from lift11 import atomic,key
SCHEMA='XIV-domain-m1-lower-v1'
OPTS=dict(nonlinear=False,cut_scope='coordinate',nonlinear_scope='all',maximal_only=True,
          pair_bounds=True,degree_hull=True,triple_scope='none',compact_objectives=True,
          point_pairs=True,plane_budget=True,point_triples=True,star_bounds=True,
          rounded_budget=True,nonfano=True,extra=False,rich=True,point_mix=True,round_degree=True)
def vectors(D,domain):
    assert set(domain)=={'m1_min'} and type(domain['m1_min'])==int and 1<=domain['m1_min']<=27
    return tuple(m for m in c.raw(D) if m[0]>=domain['m1_min'])
def domain_digest(D,domain):return hashlib.sha256(key(vectors(D,domain)).encode()).hexdigest()
def membership(D,prefix,domain):
    allowed=vectors(D,domain);assert allowed
    assert prefix and prefix[0]>=domain['m1_min'],'A certificate must cover every vector in the node'
    assert all(m[0]>=domain['m1_min'] for m in c.descendants(D,tuple(prefix)))
    return len(c.descendants(D,tuple(prefix)))
def build(D,cfg,domain,options):
    assert D==7 and cfg=='617'
    candidates=vectors(D,domain);original=c.descendants
    def selected(d,p):
        assert d==D and tuple(p)==()
        return candidates
    c.descendants=selected
    try:lp,info,types=c.build(D,cfg,(),**options)
    finally:c.descendants=original
    # Preserve the final degree-bound rows used by the numerical adapter.
    insertion=len(lp.rows)-2*D
    lp.add([(lp.col['m1'],-1)],'<',-domain['m1_min'],'XIV domain m1 lower')
    lp.rows.insert(insertion,lp.rows.pop())
    return lp,info,types

def check_record(r):
    assert (r['schema'],r['D'],r['q'],r['n'],r['config'],r['prefix'],r['family'])==(SCHEMA,7,5,46,'617',[],c.FAMILY)
    assert r['raw_domain_count']==len(vectors(r['D'],r['domain']))
    assert r['raw_domain_sha256']==domain_digest(r['D'],r['domain'])
    from maximal14 import guard,digest
    from budget14 import guard as bg
    from star14 import guard as sg
    from roundbudget14 import guard as rg
    assert r['projection_table_sha256']==digest('out/union_tables_q5.json')
    assert r['reduction_sha256']==guard(r['D'],r['config'])
    assert r['budget_sha256']==bg() and r['star_sha256']==sg(r['D'],r['config']) and r['rounded_budget_sha256']==rg()
    options={k:r[k] for k in OPTS};assert options==OPTS
    if r.get('point_subset',False):options['point_subset']=True
    return options

def replay_record(r):
    options=check_record(r);lp,_,_=build(r['D'],r['config'],r['domain'],options)
    result=dict(c.verify(lp,r['multipliers']),rows=len(lp.rows),columns=len(lp.names))
    print('EXACT_DOMAIN14_REPLAYED',r['domain'],result,flush=True)
    return result

def register_attempt(path):
    with open('out/coverage14_D7.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        ledger=json.load(c.path(7).open());attempts=ledger.setdefault('domain_attempts',[])
        if str(path) not in attempts:attempts.append(str(path))
        record=json.load(open(path))
        for flag in ('point_pairs','point_triples','nonfano','rich','point_mix','round_degree','point_subset'):
            if record['options'].get(flag,False):ledger[flag+'_used']=True
        c.check_tree(ledger);c.refresh(ledger)

def solve(a):
    if a.engine=='homogeneous':
        from homogeneous14 import HomogeneousMaster
        c.Master14=HomogeneousMaster
    else:
        from scaled14 import ScaledMaster
        c.Master14=ScaledMaster
    from thin14 import ProjectedMaster
    assert re.fullmatch('[A-Za-z0-9_]+',a.tag)
    lock=open('out/domain14_'+a.tag+'.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    D=7;cfg='617';domain=dict(m1_min=a.m1_min);c.CERT_SUFFIX='_'+a.tag
    options=dict(OPTS,point_subset=a.point_subset)
    path=Path(f'out/domain14_D7_617_{a.tag}.json');cert=Path(f'out/cert14_D7_617_domain_{a.tag}.json')
    if path.exists():assert json.load(path.open())['domain']==domain
    status=dict(schema=SCHEMA,D=D,config=cfg,domain=domain,status='open',raw_domain_count=len(vectors(D,domain)),raw_domain_sha256=domain_digest(D,domain),options=options,reason='building or solving; no exclusion yet')
    atomic(path,status);register_attempt(path)
    if time.time()>=c.CAMPAIGN_STOP:return
    def stop(*_):c.stop=True
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    lp,_,_=build(D,cfg,domain,options)
    if c.stop:return
    master=ProjectedMaster(lp,a.method);master.solver.setOptionValue('solver',a.solver)
    if a.solver=='ipm':master.solver.setOptionValue('run_crossover','off')
    master.solver.setOptionValue('output_flag',True);master.solver.setOptionValue('log_to_console',False)
    master.solver.setOptionValue('log_file',f'out/solver14_domain_{a.tag}.log')
    started=time.monotonic();print('DOMAIN14_SOLVING',domain,flush=True)
    outcome,data,detail=master.solve(lp,min(a.seconds,max(1,c.CAMPAIGN_STOP-time.time())))
    detail.update(seconds=time.monotonic()-started,solver=a.solver,requested_method=a.method,engine=a.engine)
    status.update(status=outcome,detail=detail,reason='exact certificate pending independent import' if outcome=='LP-certified' else 'unfinished domain; no coverage credit')
    if outcome=='LP-certified':
        from maximal14 import guard,digest
        from budget14 import guard as bg
        from star14 import guard as sg
        from roundbudget14 import guard as rg
        verified=c.verify(lp,data)
        record=dict(schema=SCHEMA,D=D,q=5,n=46,config=cfg,prefix=[],family=c.FAMILY,domain=domain,
                    raw_domain_count=len(vectors(D,domain)),raw_domain_sha256=domain_digest(D,domain),
                    projection_table_sha256=digest('out/union_tables_q5.json'),reduction_sha256=guard(D,cfg),
                    budget_sha256=bg(),star_sha256=sg(D,cfg),rounded_budget_sha256=rg(),
                    multipliers=data,detail=dict(detail,**verified),**options)
        atomic(cert,record);status['certificate']=str(cert)
    atomic(path,status);print('DOMAIN14_RESULT',outcome,detail,flush=True)

def import_certificate(fn):
    lock=open('out/coverage14_D7.lock','a');fcntl.flock(lock,fcntl.LOCK_EX)
    r=json.load(open(fn));detail=replay_record(r)
    ledger=json.load(c.path(7).open());changed=[]
    assert ledger['nodes']['617:']['children']
    for idx,node in ledger['nodes'].items():
        if node['config']!='617' or len(node['prefix'])!=1 or node['prefix'][0]<r['domain']['m1_min']:continue
        membership(7,node['prefix'],r['domain'])
        if node['status'] in ('LP-certified','joint-LP-certified','exact-primal'):continue
        node.update(status='LP-certified',certificate=fn,certificate_domain=r['domain'],prefix_probe_done=True)
        changed.append(idx)
    for flag in ('point_pairs','point_triples','nonfano','rich','point_mix','round_degree','point_subset'):
        if r.get(flag,False):ledger[flag+'_used']=True
    c.check_tree(ledger);c.refresh(ledger)
    atomic('out/import14_domain_'+Path(fn).stem+'.json',dict(status='EXACT_DOMAIN14_IMPORTED',certificate=fn,nodes=changed,detail=detail,counts=ledger['counts']))
    print('EXACT_DOMAIN14_IMPORTED',changed,ledger['counts'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--solve',action='store_true');ap.add_argument('--import-certificate');ap.add_argument('--replay-certificate');ap.add_argument('--m1-min',type=int,default=1);ap.add_argument('--seconds',type=int,default=1800);ap.add_argument('--method',choices=['phase1','direct'],default='phase1');ap.add_argument('--engine',choices=['scaled','homogeneous'],default='scaled');ap.add_argument('--solver',choices=['ipm','simplex'],default='ipm');ap.add_argument('--tag',default='m1ge1');ap.add_argument('--point-subset',action='store_true');a=ap.parse_args()
    assert sum((a.solve,bool(a.import_certificate),bool(a.replay_certificate)))==1
    if a.solve:solve(a)
    elif a.import_certificate:import_certificate(a.import_certificate)
    else:replay_record(json.load(open(a.replay_certificate)))
