import argparse,json,hashlib
from pathlib import Path
from coverage14 import build,verify,descendants,FAMILY
from maximal14 import guard
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('certificate');a=ap.parse_args();r=json.load(open(a.certificate));D=r['D'];cfg=r['config'];prefix=r['prefix']
assert (r['q'],r['n'],r['family'])==(5,46,FAMILY)
assert r['projection_table_sha256']==hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest()
if r.get('rounded_budget'):
    from roundbudget14 import guard as rounded_guard
    assert rounded_guard()==r['rounded_budget_sha256']
if r.get('star_bounds'):
    from star14 import guard as star_guard
    assert star_guard(D,cfg)==r['star_sha256']
if r.get('plane_budget'):
    from budget14 import guard as budget_guard
    assert budget_guard()==r['budget_sha256']
if r.get('maximal_only'):assert guard(D,cfg)==r['reduction_sha256']
lp,_,_=build(D,cfg,prefix,r.get('nonlinear',False),r.get('cut_scope','all'),r.get('nonlinear_scope','all'),r.get('maximal_only',False),r.get('pair_bounds',False),r.get('degree_hull',False),r.get('triple_scope','none'),r.get('compact_objectives',False),r.get('point_pairs',False),r.get('plane_budget',False),r.get('point_triples',False),r.get('star_bounds',False),r.get('rounded_budget',False),r.get('nonfano',False),r.get('extra',False),r.get('rich',False),r.get('point_mix',False),r.get('round_degree',False),r.get('point_subset',False))
assert all(type(b)==int and all(type(c)==int for c in row.values()) for row,b,_,_ in lp.rows)
detail=verify(lp,r['multipliers']);record=dict(status='EXACT_PREFIX14_INDEPENDENT_REPLAY',certificate=a.certificate,certificate_sha256=hashlib.sha256(Path(a.certificate).read_bytes()).hexdigest(),all_coefficients_python_integers=True,cases=len(descendants(D,tuple(prefix))),variables=len(lp.names),rows=len(lp.rows),detail=detail)
atomic('out/replay14_'+Path(a.certificate).stem+'.json',record);print('EXACT_PREFIX14_INDEPENDENT_REPLAY',record,flush=True)
