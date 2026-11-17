"""Independent ordinary subset enumeration of core pair-cost profiles."""
import argparse
from collections import defaultdict
from geometry14 import geometries
from nonlinear14 import extrema,CORE
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);a=ap.parse_args();records=[]
for g in geometries(a.D):
    if tuple(g['signature'][:2]) not in CORE:continue
    masks=tuple(g['masks']);r=g['signature'][0];n=len(masks)
    groups=[sum(1<<p for p,m in enumerate(masks) if m.bit_count()==b) for b in range(r+1)]
    lines=[sum(1<<p for p,m in enumerate(masks) if m>>i&1) for i in range(r)]
    direct={}
    for bits in range(1<<n):
        profile=tuple((bits&group).bit_count() for group in groups)
        counts=[(bits&line).bit_count() for line in lines];value=sum(c*(c-1)//2 for c in counts)
        if profile not in direct:direct[profile]=[value,value]
        else:direct[profile]=[min(direct[profile][0],value),max(direct[profile][1],value)]
    assert direct=={p:[lo,hi] for p,lo,hi,_,_ in extrema(masks,r)}
    records.append(dict(signature=g['signature'],subsets=1<<n,profiles=len(direct)));print('EXACT_BINARY14_ENUMERATOR_AGREEMENT',a.D,r,n,1<<n,flush=True)
atomic(f'out/nonlinear14_D{a.D}_independent.json',dict(D=a.D,status='EXACT_BINARY14_ENUMERATOR_AGREEMENT',records=records,mismatches=0))
