"""Import an isolated root probe into the ledger after full exact reconstruction."""
import argparse,json,fcntl,hashlib
from pathlib import Path
import coverage14 as c
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True,choices=[6,7]);ap.add_argument('--config',required=True,choices=['515','617']);ap.add_argument('--tag',default='');ap.add_argument('--prefix',default='');ap.add_argument('--wait-lock',action='store_true');a=ap.parse_args();prefix=tuple(map(int,a.prefix.split(','))) if a.prefix else ();suffix=('_'+a.tag) if a.tag else ('_p'+'_'.join(map(str,prefix)) if prefix else '')
lock=open(f'out/coverage14_D{a.D}.lock','a');fcntl.flock(lock,fcntl.LOCK_EX if a.wait_lock else fcntl.LOCK_EX|fcntl.LOCK_NB)
ledger=json.load(c.path(a.D).open());probe=json.load(open(f'out/probe14_D{a.D}_{a.config}{suffix}.json'));idx=c.nid(a.config,prefix);n=probe['nodes'][idx];old=ledger['nodes'][idx]
assert probe['D']==ledger['D'] and n['prefix']==list(prefix) and n['config']==a.config
if n['status']=='LP-certified':
    record=json.load(open(n['certificate']));assert (record['D'],record['q'],record['n'],record['config'],record['prefix'],record['family'])==(a.D,5,46,a.config,list(prefix),c.FAMILY)
    assert record['projection_table_sha256']==hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest()
    if record.get('rounded_budget'):
        from roundbudget14 import guard as rounded_guard
        assert rounded_guard()==record['rounded_budget_sha256']
    if record.get('star_bounds'):
        from star14 import guard as star_guard
        assert star_guard(a.D,a.config)==record['star_sha256']
    if record.get('plane_budget'):
        from budget14 import guard as budget_guard
        assert budget_guard()==record['budget_sha256']
    if record.get('maximal_only'):
        from maximal14 import guard
        assert guard(a.D,a.config)==record['reduction_sha256']
    lp,_,_=c.build(a.D,a.config,list(prefix),record.get('nonlinear',False),record.get('cut_scope','all'),record.get('nonlinear_scope','all'),record.get('maximal_only',False),record.get('pair_bounds',False),record.get('degree_hull',False),record.get('triple_scope','none'),record.get('compact_objectives',False),record.get('point_pairs',False),record.get('plane_budget',False),record.get('point_triples',False),record.get('star_bounds',False),record.get('rounded_budget',False),record.get('nonfano',False),record.get('extra',False),record.get('rich',False),record.get('point_mix',False),record.get('round_degree',False),record.get('point_subset',False));detail=c.verify(lp,record['multipliers'])
    old.update(status='LP-certified',certificate=n['certificate'],prefix_probe_done=True)
    print('EXACT_NODE14_IMPORTED',a.D,idx,detail,flush=True)
elif old['status'] not in ('LP-certified','joint-LP-certified','exact-primal'):
    old['prefix_probe_done']=True
    old['reason']=n.get('reason','unfinished isolated prefix probe')
    if n.get('joint_pending'):old['joint_pending']=True
    if old['status']=='open' and n['status']=='split':c.split(ledger,idx)
for p in n.get('probes',[]):
    if p not in old.setdefault('probes',[]):old['probes'].append(p)
if probe.get('nonlinear_used'):ledger['nonlinear_used']=True
if probe.get('triple_used'):ledger['triple_used']=True
if probe.get('point_pairs_used'):ledger['point_pairs_used']=True
if probe.get('point_triples_used'):ledger['point_triples_used']=True
if probe.get('nonfano_used'):ledger['nonfano_used']=True
if probe.get('extra_used'):ledger['extra_used']=True
if probe.get('rich_used'):ledger['rich_used']=True
if probe.get('point_mix_used'):ledger['point_mix_used']=True
if probe.get('round_degree_used'):ledger['round_degree_used']=True
if probe.get('point_subset_used'):ledger['point_subset_used']=True
c.check_tree(ledger);c.refresh(ledger);print('IMPORT14_CHECKPOINT',a.D,ledger['counts'],flush=True)
