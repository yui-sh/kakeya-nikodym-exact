"""Exhaustive projective direction-subset triples-minus-quadruples bounds."""
import subprocess,json
from math import comb
from geom import AG
from pc6_validate import dot
from lift11 import atomic
for q in (3,4,5):
    ag=AG(3,q);pts=ag.dirs;n=len(pts);D=7
    lines=[sum(1<<i for i,p in enumerate(pts) if dot(ag.F,normal,p)==0) for normal in pts]
    assert len(set(lines))==n and all(l.bit_count()==q+1 for l in lines)
    for i in range(n):
        for j in range(i):assert sum(bool(l>>i&1) and bool(l>>j&1) for l in lines)==1
    inp=f'{n} {D} {n}\n'+' '.join(map(str,lines))+'\n'
    records=json.loads(subprocess.run(['out/pg14_nonfano_enum'],input=inp,text=True,capture_output=True,check=True).stdout)
    for rec in records:
        j=rec['degree'];assert rec['subsets']==comb(n,j)
        for prefix in ('min','max'):
            mask=rec[prefix+'_mask'];assert mask.bit_count()==j and mask>>n==0
            value=sum(comb((mask&l).bit_count(),3)-comb((mask&l).bit_count(),4) for l in lines)
            assert value==rec['minimum' if prefix=='min' else 'maximum']
            rec[prefix+'_directions']=[p for i,p in enumerate(pts) if mask>>i&1]
    if q==5:assert (records[7]['minimum'],records[7]['maximum'])==(2,6)
    atomic(f'out/nonfano14_projective_q{q}.json',dict(q=q,D=D,status='EXACT_NONFANO14_PROJECTIVE_ENUMERATED',directions=pts,line_masks=lines,records=records))
    print('EXACT_NONFANO14_PROJECTIVE_ENUMERATED',q,[(r['degree'],r['subsets'],r['minimum'],r['maximum']) for r in records],flush=True)
