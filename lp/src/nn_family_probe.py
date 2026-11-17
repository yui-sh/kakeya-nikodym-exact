import sys, numpy as np
sys.path.insert(0,'.')
from bm_family_probe import overlap_for
from sympy import primerange
for q in primerange(7,140):
    S=sorted({a*a%q for a in range(q)}); sq=np.zeros(q,bool); sq[S]=True
    G=np.array([(a,b) for a in S for b in S]); x,y=np.meshgrid(np.arange(q),np.arange(q),indexing='ij')
    base=sq[(x*x-y)%q].astype(np.int32); nb=1-base
    vx,vy=[a.ravel() for a in np.meshgrid(np.arange(q),np.arange(q),indexing='ij')]
    n=min(a for a in range(2,q) if not sq[a])
    k,l=(n,n) if q%4==3 else (1,1)
    b,e=overlap_for(q,(1,k),(0,l),base,nb,G,vx,vy); c=b+e
    g=((q+1)//2)**2; target=(g+2*q-(q%4==1))//2
    i=int(c.argmax()); top=np.sort(c)[::-1][:4]
    print(q,q%8,'fam',(k,l),'max',int(c.max()),'target',target,'diff',int(c.max())-target,'argmax v',(int(vx[i]),int(vy[i])),'b,e',int(b[i]),int(e[i]),'#at max',int((c==c.max()).sum()),flush=True)
