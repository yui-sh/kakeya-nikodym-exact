"""Import an independent all-scope root result only after exact reconstruction."""
import fcntl,json
from coverage11 import LEDGER,build,node_lp,verify,refresh,check_tree
lock=open('out/coverage11.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
ledger=json.load(LEDGER.open())
result=json.load(open('out/coverage12_all_root_617_result.json'))
node=ledger['nodes']['617:'];node['session12_all_root_probe']=result
if result['status']=='LP-certified':
    path='out/coverage12_all_cert_617_root.json'
    rec=json.load(open(path))
    assert rec['config']=='617' and rec['prefix']==[] and rec['cut_scope']=='all' and rec['model_kind']=='prefix-filtered'
    base,_,types=build('all','617',[])
    detail=verify(node_lp(base,types,'617',[]),rec['multipliers'])
    node.update(status='LP-certified',certificate=path,prefix_probe_done=True,reason='XII all-scope root certificate independently replayed')
    print('EXACT_ROOT12_IMPORTED',detail,flush=True)
check_tree(ledger);refresh(ledger)
print('COVERAGE12_IMPORTED',ledger['counts'],flush=True)
