import sys, numpy as np
sys.path.insert(0,'.')
from bm_family_probe import overlap_for
for q in map(int,sys.argv[1:]):
    S=sorted({a*a%q for a in range(q)}); sq=np.zeros(q,bool); sq[S]=True
    G=np.array([(a,b) for a in S for b in S]); x,y=np.meshgrid(np.arange(q),np.arange(q),indexing='ij')
    base=sq[(x*x-y)%q].astype(np.int32); nb=1-base
    vx,vy=[a.ravel() for a in np.meshgrid(np.arange(q),np.arange(q),indexing='ij')]
    tab={}
    for k in range(q):
        for l in range(1,q):
            b,e=overlap_for(q,(1,k),(0,l),base,nb,G,vx,vy); tab[(k,l)]=int((b+e).max())
    m=max(tab.values()); g=((q+1)//2)**2
    win=[kl for kl,v in tab.items() if v==m]
    print(q,q%8,'max',m,'G/2+q',g/2+q,'#win',len(win),'win sample',win[:6],
          'k-sq?',sorted({(bool(sq[k]),bool(sq[l])) for k,l in win}),flush=True)
