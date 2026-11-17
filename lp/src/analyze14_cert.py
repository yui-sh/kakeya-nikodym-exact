"""Exact support analysis and optional whole-family zeroing experiments."""
import argparse,json,hashlib
from pathlib import Path
from collections import defaultdict
import coverage14 as c
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('certificate');ap.add_argument('--support-only',action='store_true');a=ap.parse_args()
r=json.load(open(a.certificate));D=r['D']
if r.get('domain'):
    from domain14 import check_record,build as domain_build
    lp,_,_=domain_build(D,r['config'],r['domain'],check_record(r))
else:
    lp,_,_=c.build(D,r['config'],r['prefix'],r.get('nonlinear',False),r.get('cut_scope','all'),r.get('nonlinear_scope','all'),r.get('maximal_only',False),r.get('pair_bounds',False),r.get('degree_hull',False),r.get('triple_scope','none'),r.get('compact_objectives',False),r.get('point_pairs',False),r.get('plane_budget',False),r.get('point_triples',False),r.get('star_bounds',False),r.get('rounded_budget',False),r.get('nonfano',False),r.get('extra',False),r.get('rich',False),r.get('point_mix',False),r.get('round_degree',False),r.get('point_subset',False))
y=r['multipliers'];detail=c.verify(lp,y)
def family(name):
    for p in ('XIV domain','XIV point subset crossing','XIV rounded degree','XIV extra','XIV point rich-secants','XIV point mixed','XIV point triples-minus-quadruples','XIV star-four','XIV rounded plane budget','XIV point crossing','XIV point coplanar','XIV triple','XIV subset','XIV nonlinear','XIV degree hull','PB14','F8 Y1','F8 Y2','F8 Z1','F8 Z2','F6'):
        if name.startswith(p):return p
    return name.split(':')[0].split(' ')[0]
groups=defaultdict(lambda:dict(rows=0,nonzero=0,rhs_contribution=0,max_abs_multiplier=0))
for v,(_,b,_,name) in zip(y,lp.rows):
    z=groups[family(name)];z['rows']+=1;z['nonzero']+=bool(v);z['rhs_contribution']+=v*b;z['max_abs_multiplier']=max(z['max_abs_multiplier'],abs(v))
ablations=[]
for prefix in (() if a.support_only else ('XIV point triples-minus-quadruples','XIV star-four','XIV rounded plane budget','XIV point crossing','XIV point coplanar','XIV triple','XIV degree hull','XIV nonlinear','XIV subset','F8')):
    candidate=[0 if row[3].startswith(prefix) else v for v,row in zip(y,lp.rows)]
    try:result=dict(valid=True,**c.verify(lp,candidate))
    except AssertionError:result=dict(valid=False)
    ablations.append(dict(removed=prefix,**result))
record=dict(certificate=a.certificate,sha256=hashlib.sha256(Path(a.certificate).read_bytes()).hexdigest(),detail=detail,groups=dict(groups),ablations=ablations)
atomic('out/analysis14_'+Path(a.certificate).stem+'.json',record)
print('CERTIFICATE14_ANALYSIS',json.dumps(record),flush=True)
