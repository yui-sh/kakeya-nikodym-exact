"""Exact local high-r plane column catalogue and implicit assignment enumeration."""
from collections import Counter, defaultdict
from itertools import combinations, product
import json, subprocess, hashlib
from pathlib import Path
from geom import AG
import numpy as np
from scipy import sparse as sp
import pencil9 as p9

def freeze(x):
    return tuple(map(freeze,x)) if isinstance(x,list) else x

def signature(ag, lines, S, degree=None):
    r=len(lines); U=set().union(*(ag.lines[l][1] for l in lines))
    t=tuple(sum(len(ps&S)==j for _,ps in ag.lines) for j in range(1,ag.q+1))
    dirs={ag.lines[l][0] for l in lines}
    e=tuple(sorted(tuple(sorted(len(ag.lines[l][1]&S) for l in ag.lines_by_dir[d])) for d in range(ag.q+1) if d not in dirs)) if r>=5 and r<ag.q+1 else None
    h=tuple(Counter(degree[p] for p in S)[j] for j in range(1,max(degree.values())+1)) if degree else None
    return (r,len(S),h,len(U),t,e)

def geometries(q=5):
    """Orbit traversal using generators of AGL(2,q), q=5. All extras enumerated.
    One representative per (line orbit,k,t,e); extras have no line incidences.
    """
    path=Path('out/plane10_geometries.json')
    if path.exists():return json.load(path.open())
    assert q==5
    ag=AG(2,q)
    maps=[lambda x,y:((x+1)%5,y),lambda x,y:(x,(y+1)%5),lambda x,y:(y,x),lambda x,y:((x+y)%5,y),lambda x,y:((2*x)%5,y)]
    perms=[]
    for f in maps:
        pp=[ag.pidx[f(*p)] for p in ag.points]
        perms.append([ag.line_of[frozenset(pp[p] for p in ps)] for _,ps in ag.lines])
    out=[]
    for r in (5,6):
        seen=set(); orbits=[]
        for ds in combinations(range(6),r):
            for ls in product(*(ag.lines_by_dir[d] for d in ds)):
                ls=tuple(sorted(ls))
                if ls in seen:continue
                orbit={ls};stack=[ls]
                while stack:
                    cur=stack.pop()
                    for perm in perms:
                        nxt=tuple(sorted(perm[l] for l in cur))
                        if nxt not in orbit:orbit.add(nxt);stack.append(nxt)
                seen.update(orbit);orbits.append((ls,len(orbit)))
        assert len(seen)==(6 if r==5 else 1)*5**r
        for ls,sz in orbits:
            U=set().union(*(ag.lines[l][1] for l in ls));outside=sorted(set(range(25))-U);sgseen=set()
            for bits in range(1<<len(outside)):
                S=U|{p for i,p in enumerate(outside) if bits>>i&1}
                sig=signature(ag,ls,S)
                if sig in sgseen:continue
                sgseen.add(sig)
                masks={p:sum(1<<i for i,l in enumerate(ls) if p in ag.lines[l][1]) for p in S}
                pts=sorted(S,key=lambda p:(masks[p],p))
                out.append(dict(lines=ls,points=pts,masks=[masks[p] for p in pts],signature=sig,orbit_size=sz))
        print('GEOMETRY',r,'orbits',len(orbits),'covered',len(seen),flush=True)
    json.dump(out,path.open('w'));return out

def build(config,m):
    lp,info,types,grow,zrow=p9.build(config,m)
    nr={};gr={}
    for i,t in enumerate(types):
        if t[0]>=5:
            nr[i]=len(lp.rows);lp.add([(i,-1)],'=',0,f'P10 N{i}')
    for g,(key,di) in enumerate(info['_Y']):
        if key[0]>=5:
            gr[g]=len(lp.rows);lp.add([(lp.col[f'G{g}'],-1)],'=',0,f'P10 G{g}')
    return lp,info,types,grow,zrow,nr,gr

def counts(col,grow,zrow,nr,gr):
    if col[0]=='pencil':return p9.column_counts(col[2],col[1],grow,zrow)
    _,i,gids,geo,deg=col
    c=Counter(gr[g] for g in gids);c[nr[i]]+=1;return dict(c)

class PlaneOracle:
    def __init__(self,types,info,nr,gr,tag):
        self.types,self.info,self.nr,self.gr,self.tag=types,info,nr,gr,tag
        self.geos=geometries();self.ag=AG(2,5);self.jobs=[]
        self.catalog_signature=hashlib.sha256(json.dumps([types,info['_Y'],info['_dms']],separators=(',',':')).encode()).hexdigest()
        glookup={(key,info['_dms'][di]):g for g,(key,di) in enumerate(info['_Y'])}
        for i in nr:
            r,k,h,u,t,e=types[i];key=(r,k,h)
            allowed=[(dm,g) for (kk,dm),g in glookup.items() if kk==key]
            for gi,geo in enumerate(self.geos):
                sg=freeze(geo['signature'])
                if sg[:2]+sg[3:]!=(r,k,u,t,e):continue
                self.jobs.append((i,gi,h,allowed))
        print('PLANE_JOBS',len(self.jobs),'types',len(nr),flush=True)
    def price(self,y,first=False):
        # The fully enumerated catalogue is an equivalent exact pricing oracle.
        short='617' if any(t[0]==6 for t in self.types) else '515'
        path=Path(f'out/plane10_catalog_{short}.bin')
        done=Path(f'out/plane10_catalog_{short}.log')
        check=Path(f'out/plane10_catalog_{short}_check.json')
        if not first and path.exists() and done.exists() and check.exists() and 'CATALOG_COMPLETE' in done.read_text():
            meta=json.load(check.open())
            if meta.get('model_signature')==self.catalog_signature and meta.get('sha256')==hashlib.sha256(path.read_bytes()).hexdigest():
                return self.catalog_price(y,path)
        return self.dfs_price(y,first)

    def catalog_price(self,y,path):
        assert max(map(abs,y))<10**15  # seven terms fit signed int64 exactly
        dtype=np.dtype([('type','<u4'),('geo','<u4'),('gids','<u2',(6,)),('degrees','<u8')])
        records=np.memmap(path,dtype=dtype,mode='r');assert path.stat().st_size%dtype.itemsize==0
        gc=np.zeros(65536,dtype=np.int64)
        for g,row in self.gr.items():gc[g]=y[row]
        nc=np.zeros(len(self.types),dtype=np.int64)
        for i,row in self.nr.items():nc[i]=y[row]
        costs=nc[records['type']]+gc[records['gids']].sum(axis=1)
        results=[];rng=np.random.default_rng(1010)
        for i in self.nr:
            ids=np.flatnonzero(records['type']==i)
            if not len(ids):continue
            best=ids[np.argmin(costs[ids])];chosen={int(best)}
            neg=ids[costs[ids]<0]
            if len(neg):
                chosen.update(map(int,neg[np.argsort(costs[neg],kind='stable')[:8]]))
                chosen.update(map(int,rng.choice(neg,size=min(24,len(neg)),replace=False)))
            for j in sorted(chosen):
                rec=records[j];gi=int(rec['geo']);k=self.types[i][1];packed=int(rec['degrees'])
                ds=[(packed>>(3*p))&7 for p in range(k)]
                col=['plane',i,[int(g) for g in rec['gids'] if g!=65535],gi,ds]
                self.verify(col);results.append((int(costs[j]),col))
        return results

    def dfs_price(self,y,first=False):
        inp=f'out/{self.tag}_plane.in';out=f'out/{self.tag}_plane.out'
        with open(inp,'w') as f:
            print(len(self.jobs),int(first),file=f)
            for i,gi,h,allowed in self.jobs:
                geo=self.geos[gi];print(i,gi,len(h),len(geo['points']),self.types[i][0],int(y[self.nr[i]]),file=f)
                print(*h,file=f);print(*geo['masks'],file=f);print(len(allowed),file=f)
                for dm,g in allowed:print(g,int(y[self.gr[g]]),*dm,file=f)
        with open(f'out/{self.tag}_plane_price.log','a') as log:subprocess.run(['out/plane10_enum',inp,out],stderr=log,check=True)
        results=[]
        for line in Path(out).read_text().splitlines():
            a=list(map(int,line.split()));i,gi,cost=a[:3];r=self.types[i][0];gids=a[3:3+r];degrees=a[3+r:]
            col=['plane',i,sorted(gids),gi,degrees];self.verify(col)
            assert cost==y[self.nr[i]]+sum(y[self.gr[g]] for g in gids)
            results.append((cost,col))
        return results
    def verify(self,col):
        _,i,gids,gi,ds=col;geo=self.geos[gi];t=self.types[i];r,k,h,u,sec,e=t
        # Recompute the geometric witness, independently of cached signature/masks.
        points=geo['points'];lines=geo['lines'];assert len(set(points))==k
        assert len({self.ag.lines[l][0] for l in lines})==r
        assert all(self.ag.lines[l][1]<=set(points) for l in lines)
        assert geo['masks']==[sum(1<<j for j,l in enumerate(lines) if p in self.ag.lines[l][1]) for p in points]
        check=signature(self.ag,lines,set(points),dict(zip(points,ds)))
        check=check[:2]+(check[2]+(0,)*(len(h)-len(check[2])),)+check[3:]
        assert check==t
        assert len(ds)==k and tuple(ds.count(j) for j in range(1,len(h)+1))==h
        assert all(d>=mask.bit_count() and 1<=d<=len(h) for d,mask in zip(ds,geo['masks']))
        sg=freeze(geo['signature']);assert sg[:2]+sg[3:]==t[:2]+t[3:]
        actual=sorted(tuple(sorted(d for d,mask in zip(ds,geo['masks']) if mask>>j&1)) for j in range(r))
        assert actual==sorted(self.info['_dms'][self.info['_Y'][g][1]] for g in gids)
        assert all(self.info['_Y'][g][0]==t[:3] for g in gids)
