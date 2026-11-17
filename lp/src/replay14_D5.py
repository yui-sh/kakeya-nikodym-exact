"""Read-only D=5 ledger check under the corrected degree-independent projection table."""
import json,hashlib,argparse,os
from pathlib import Path
from coverage11 import build,node_lp,verify,check_tree
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--fresh-prices',action='store_true');a=ap.parse_args()
p=Path('out/coverage11.json');before=hashlib.sha256(p.read_bytes()).hexdigest();ledger=json.load(p.open());results=[]
if a.fresh_prices:
 import planecount6 as p6
 from coverage11 import fresh_prices
 from lift13 import table as pair_table
 os.environ['COVERAGE_SESSION']='14_D5'
 universe=p6.hist_types(5,5,46)
 fresh_prices(universe);pair_table(universe,True)
 print('D5_ALL_MFREE_PRICES_FRESHLY_RECOMPUTED',flush=True)
for cfg in ('515','617'):
 node=ledger['nodes'][cfg+':'];assert node['status']=='LP-certified'
 rec=json.load(open(node['certificate']));scope=rec.get('cut_scope',ledger['cut_scope'])
 base,_,types=build(scope,cfg,[]);lp=node_lp(base,types,cfg,[])
 detail=verify(lp,rec['multipliers']);results.append(dict(config=cfg,certificate=node['certificate'],certificate_sha256=hashlib.sha256(Path(node['certificate']).read_bytes()).hexdigest(),detail=detail))
 print('EXACT_D5_CORRECTED_TABLE_REPLAYED',cfg,detail,flush=True);del lp,base,types
check_tree(ledger);after=hashlib.sha256(p.read_bytes()).hexdigest();assert after==before
atomic('out/replay14_D5_corrected.json',dict(status='EXACT_D5_CORRECTED_TABLE_REPLAYED',fresh_prices=a.fresh_prices,results=results,ledger_sha256=before,ledger_untouched=True,projection_table_sha256=hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest()))
print('D5_LEDGER_UNTOUCHED',before,flush=True)
