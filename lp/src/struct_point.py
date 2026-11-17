import sys, numpy as np
sys.path.insert(0,'.')
from bm_family_probe import overlap_for
for q in [11,19,23,31,43,47,59,101,103,107,109]:
    S=sorted({a*a%q for a in range(q)}); sq=np.zeros(q,bool); sq[S]=True
    G=np.array([(a,b) for a in S for b in S]); x,y=np.meshgrid(np.arange(q),np.arange(q),indexing='ij')
    base=sq[(x*x-y)%q].astype(np.int32); nb=1-base
    vx,vy=[a.ravel() for a in np.meshgrid(np.arange(q),np.arange(q),indexing='ij')]
    g=((q+1)//2)**2; target=(g+2*q-(q%4==1))//2
    b,e=overlap_for(q,(1,q-1),(0,1),base,nb,G,vx,vy); c=b+e
    w=(-pow(4,-1,q))%q; i=w*q
    print(q,q%4,'target',target,'struct(k=-1,l=1,v=(-1/4,0))',int(c[i]),'b',int(b[i]),'e',int(e[i]),'fam max',int(c.max()),flush=True)
