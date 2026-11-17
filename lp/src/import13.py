"""Replay and import one independent Session XIII prefix probe."""
import argparse,fcntl,json
from coverage11 import LEDGER,build,node_lp,verify,refresh,check_tree,nid
ap=argparse.ArgumentParser();ap.add_argument('tag');a=ap.parse_args()
lock=open('out/coverage11.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
ledger=json.load(LEDGER.open());result=json.load(open('out/coverage13_probe_'+a.tag+'.json'))
idx=nid(result['config'],result['prefix']);node=ledger['nodes'][idx]
node.setdefault('independent_probes',{})[a.tag]=result
if result['status']=='LP-certified':
    path='out/coverage13_cert_'+a.tag+'.json';rec=json.load(open(path))
    assert rec['config']==node['config'] and rec['prefix']==node['prefix'] and rec['cut_scope']==result['cut_scope'] and rec['model_kind']=='prefix-filtered'
    base,_,types=build(rec['cut_scope'],rec['config'],rec['prefix'])
    detail=verify(node_lp(base,types,rec['config'],rec['prefix']),rec['multipliers'])
    node.update(status='LP-certified',certificate=path,prefix_probe_done=True,reason='XIII independent prefix certificate replayed')
    print('EXACT_PREFIX13_IMPORTED',idx,detail,flush=True)
check_tree(ledger);refresh(ledger)
print('COVERAGE13_IMPORTED',ledger['counts'],flush=True)
