"""One isolated configuration worker, with serialized exact node imports.

Two workers may solve disjoint configurations. The shared ledger is only locked
for an import; each solver reads a snapshot and cannot overwrite another node.
"""
import argparse, fcntl, json, subprocess, sys, time
from pathlib import Path
import coverage14 as c
from lift11 import atomic

ap=argparse.ArgumentParser()
ap.add_argument('--D',type=int,choices=[6,7],required=True)
ap.add_argument('--config',choices=['515','617'],required=True)
ap.add_argument('--seconds',type=int,default=600)
ap.add_argument('--engine',choices=['plain','scaled','homogeneous','rows'],default='plain')
ap.add_argument('--solver',choices=['simplex','ipm'],default='simplex')
ap.add_argument('--degree-hull',action='store_true')
ap.add_argument('--pair-bounds',action='store_true')
ap.add_argument('--level',type=int,default=0,help='Probe existing ancestor nodes at this depth; zero uses the frontier')
ap.add_argument('--cut-scope',choices=['all','core','pair','coordinate'],default='pair')
ap.add_argument('--triple-scope',choices=['none','all','pair','singletons'],default='none')
ap.add_argument('--nonlinear-scope',choices=['all','pair','singletons'],default='singletons')
ap.add_argument('--compact-objectives',action='store_true')
ap.add_argument('--point-pairs',action='store_true')
ap.add_argument('--point-triples',action='store_true')
ap.add_argument('--plane-budget',action='store_true')
ap.add_argument('--star-bounds',action='store_true')
ap.add_argument('--rounded-budget',action='store_true')
ap.add_argument('--nonfano',action='store_true')
ap.add_argument('--extra',action='store_true')
ap.add_argument('--rich',action='store_true')
ap.add_argument('--point-mix',action='store_true')
ap.add_argument('--round-degree',action='store_true');ap.add_argument('--point-subset',action='store_true')
ap.add_argument('--m1-min',type=int,default=0)
ap.add_argument('--m1-max',type=int,default=46)
ap.add_argument('--reverse',action='store_true')
ap.add_argument('--parallel-worker',action='store_true')
ap.add_argument('--tag',required=True)
a=ap.parse_args()
assert 0<=a.m1_min<=a.m1_max<=46
lock=open(f'out/campaign14_D{a.D}_{a.config}'+('_'+a.tag if a.parallel_worker else '')+'.lock','a')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
statepath=Path(f'out/campaign14_D{a.D}_{a.config}_{a.tag}.json')
state=json.load(statepath.open()) if statepath.exists() else dict(D=a.D,config=a.config,tag=a.tag,attempted=[],results=[],status='running')
assert (state['D'],state['config'],state['tag'])==(a.D,a.config,a.tag)
state.update(status='running',options=vars(a),degree_hull=a.degree_hull,pair_bounds=a.pair_bounds,seconds=a.seconds)
stalls=0
while time.time()+a.seconds+240<c.CAMPAIGN_STOP:
    ledger=json.load(c.path(a.D).open());c.check_tree(ledger)
    pending=[idx for idx,n in c.frontier(ledger) if n['config']==a.config and n.get('joint_pending') and n['status']=='open']
    if pending:
        state.update(status='joint_pending',pending=pending);break
    if ledger['nodes'][a.config+':']['status']=='LP-certified':
        state['status']='root_certified';break
    candidates=ledger['nodes'].items() if a.level else c.frontier(ledger)
    todo=[(idx,n) for idx,n in candidates if n['config']==a.config and n['status'] in ('open','split') and idx not in state['attempted'] and (not a.level or len(n['prefix'])==a.level) and (not n['prefix'] or a.m1_min<=n['prefix'][0]<=a.m1_max)]
    todo.sort(key=lambda item:tuple(item[1]['prefix']),reverse=a.reverse)
    if not todo:
        state['status']='no_unattempted_open_nodes';break
    idx,node=todo[0];prefix=','.join(map(str,node['prefix']))
    tag=a.tag+'_'+('_'.join(map(str,node['prefix'])) or 'root')
    base=[sys.executable,'-u',{'plain':'probe14.py','scaled':'scaled14.py','homogeneous':'homogeneous14.py','rows':'row14.py'}[a.engine],'--D',str(a.D),'--config',a.config,'--prefix',prefix,'--tag',tag,'--seconds',str(a.seconds),'--solver',a.solver,'--nonlinear','--maximal-only','--cut-scope',a.cut_scope,'--nonlinear-scope',a.nonlinear_scope,'--triple-scope',a.triple_scope]
    if a.degree_hull:base.append('--degree-hull')
    if a.pair_bounds:base.append('--pair-bounds')
    for flag in ('compact_objectives','point_pairs','point_triples','plane_budget','star_bounds','rounded_budget','nonfano','extra','rich','point_mix','round_degree','point_subset'):
        if getattr(a,flag):base.append('--'+flag.replace('_','-'))
    if a.parallel_worker:base.append('--parallel-probe')
    logfile=f'out/campaign14_D{a.D}_{a.config}_{tag}.log'
    state.update(current=idx,current_log=logfile);atomic(statepath,state)
    print('CAMPAIGN14_START',idx,logfile,flush=True)
    with open(logfile,'ab',buffering=0) as out:
        result=subprocess.run(base,stdout=out,stderr=subprocess.STDOUT)
        if result.returncode:
            state.update(status='probe_failed',returncode=result.returncode);break
        result=subprocess.run([sys.executable,'-u','import14.py','--D',str(a.D),'--config',a.config,'--prefix',prefix,'--tag',tag,'--wait-lock'],stdout=out,stderr=subprocess.STDOUT)
        if result.returncode:
            state.update(status='import_failed',returncode=result.returncode);break
    ledger=json.load(c.path(a.D).open());n=ledger['nodes'][idx]
    state['attempted'].append(idx)
    state['results'].append(dict(node=idx,status=n['status'],counts=ledger['counts'],utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    state.pop('current',None);atomic(statepath,state)
    print('CAMPAIGN14_IMPORTED',idx,n['status'],ledger['counts'],flush=True)
    stalls=0 if n['status']=='LP-certified' else stalls+1
    if stalls>=3:
        state['status']='three_stalls_review_stronger_model';break
else:state['status']='campaign_time_reserve'
atomic(statepath,state)
print('CAMPAIGN14_STOP',state['status'],flush=True)
