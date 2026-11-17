"""Check original-integer certificate lifting and positive rounding repair with RD variables."""
import json
from plane10_cert import bump
from lift11 import atomic
import coverage14 as c
r=json.load(open('out/cert14_D7_617_0_2_extra02.json'))
opts=dict(nonlinear=True,cut_scope='all',nonlinear_scope='all',maximal_only=True,pair_bounds=True,degree_hull=True,triple_scope='all',compact_objectives=True,point_pairs=True,plane_budget=True,point_triples=True,star_bounds=True,rounded_budget=True,nonfano=True,extra=True,rich=True)
old,_,_=c.build(7,'617',(0,2),**opts);new,_,_=c.build(7,'617',(0,2),point_mix=True,round_degree=True,**opts)
assert old.names==new.names[:len(old.names)]
y=[0]*len(new.rows);i=0
for j,row in enumerate(new.rows):
 if i<len(old.rows) and row==old.rows[i]:y[j]=r['multipliers'][i];i+=1
assert i==len(old.rows)
result=c.verify(new,y);assert result['rhs']==r['detail']['rhs']
by,bc,bb=bump(new,{},{});assert min(bc)>0 and bb>0
atomic('out/extended14_check.json',dict(status='EXACT_EXTENDED14_CHECKED',old_columns=len(old.names),new_columns=len(new.names),old_rows=len(old.rows),new_rows=len(new.rows),lifted_certificate=result,bump_min=min(bc),bump_rhs=bb))
print('EXACT_EXTENDED14_CHECKED',result,min(bc),bb,flush=True)
