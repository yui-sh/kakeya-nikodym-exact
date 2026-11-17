"""Exact all-column Farkas repair/replay; success only after BOTH exhaustive oracles."""
import argparse,json,math
import numpy as np
from plane10 import build,PlaneOracle
from pencil9 import Oracle

def combo(lp,y):
    col=[0]*len(lp.names);rhs=0
    for a,(row,b,s,name) in zip(y,lp.rows):
        if a:
            for j,c in row.items():col[j]+=a*c
            rhs+=a*b
    return col,rhs

def bump(lp,nr,zr):
    y=[int(s=='=' and all(c>0 for c in row.values())) for row,b,s,name in lp.rows]
    for i,(_,_,_,name) in enumerate(lp.rows):
        if name.startswith('F8 Y1'):y[i]+=1
        if name.startswith('F8 Z1'):y[i]+=2
    for i in nr.values():y[i]+=1
    for i in zr.values():y[i]+=1
    def adjust(prefix,rowprefix,sign=1):
        c,_=combo(lp,y);w=min((c[j] for j,n in enumerate(lp.names) if n.startswith(prefix)),default=1)
        if w<=0:
            for i,row in enumerate(lp.rows):
                if row[3].startswith(rowprefix):y[i]+=sign*(1-w)
    adjust('M','H6');adjust('m','deg: sum m_j');adjust('Y','F6 link',-1);adjust('X','F6 R1');adjust('N','I2a');adjust('W','wlines: sum W');adjust('V','wdeg: sum V')
    c,b=combo(lp,y);assert min(c)>0 and b>0
    # All implicit plane/pencil columns have bump cost exactly one.
    assert all(y[i]==1 for i in nr.values()) and all(y[i]==1 for i in zr.values())
    return y,c,b

def run(short,replay=False):
    tag='plane10_'+short;old=json.load(open('out/pencil9_'+short.replace('full','')+'.json'));config=tuple(old['config']);m=tuple(old['m'])
    lp,info,types,gr,zr,nr,pgr=build(config,m);po=PlaneOracle(types,info,nr,pgr,tag+'_cert');wo=Oracle(info,m,gr,zr,tag+'_cert')
    def prices(y):
        assert max(map(abs,y))<10**15
        pp=po.price(y);wp=wo.price(y)
        return min((c for c,col in pp),default=0),min((c for c,z,gs in wp),default=0)
    if replay:
        record=json.load(open(f'out/{tag}_farkas.json'));y=record['multipliers']
        assert len(y)==len(lp.rows) and all(isinstance(v,int) for v in y)
        assert all(y[i]>=0 for i,row in enumerate(lp.rows) if row[2]=='<')
        c,b=combo(lp,y);p,w=prices(y);assert min(c)>=0 and p>=0 and w>=0 and b<0
        print(tag,'EXACT_ALL_COLUMN_FARKAS_VERIFIED',b,'base_min',min(c),'plane_min',p,'pencil_min',w,flush=True);return
    dual=np.load(f'out/{tag}_dual.npz')['y'];by,bc,bb=bump(lp,nr,zr)
    for scale in (10**9,10**11):
        y=[int(round(v*scale)) for v in dual]
        for i,row in enumerate(lp.rows):
            if row[2]=='<':y[i]=max(0,y[i])
        c,b=combo(lp,y);pp=po.price(y);wp=wo.price(y)
        p=min((cost for cost,col in pp),default=0);w=min((cost for cost,z,gs in wp),default=0)
        # Repair each implicit family using its population, before the small generic rounding repair.
        pd={i:0 for i in nr};zd={z:0 for z in zr}
        for cost,col in pp:pd[col[1]]=max(pd[col[1]],-cost)
        for cost,z,gs in wp:zd[z]=max(zd[z],-cost)
        yy=y.copy()
        for i,d in pd.items():yy[nr[i]]+=d
        for z,d in zd.items():yy[zr[z]]+=d
        alpha=max((d+types[i][0]-1)//types[i][0] for i,d in pd.items())
        wd={}
        for z,d in zd.items():
            wi=info['_Z'][z][1];wd[wi]=max(wd.get(wi,0),d)
        beta=max(wd.values(),default=0)
        for i,row in enumerate(lp.rows):
            name=row[3]
            if name.startswith('I2b:'):yy[i]+=alpha
            if name.startswith('F8 Z1 w'):yy[i]+=wd.get(int(name.split('w')[-1]),0)
            if name=='wlines: sum W = nd':yy[i]+=beta
        cc,rhs=combo(lp,yy)
        T=max(0,max((-v+bc[j]-1)//bc[j] for j,v in enumerate(cc)))
        yy=[a+T*bbv for a,bbv in zip(yy,by)];cc,rhs=combo(lp,yy)
        print(tag,'repair',scale,'plane_alpha',alpha,'pencil_beta',beta,'T',T,'rhs',rhs,flush=True)
        if rhs>=0:continue
        assert min(cc)>=0;pp,ww=prices(yy);assert pp>=0 and ww>=0
        record=dict(multipliers=yy,rhs=rhs,base_min=min(cc),plane_min=pp,pencil_min=ww,scale=scale,bump=T,plane_alpha=alpha,pencil_beta=beta)
        json.dump(record,open(f'out/{tag}_farkas.json','w'),indent=1)
        print(tag,'EXACT_ALL_COLUMN_FARKAS_SAVED',rhs,flush=True);return
    raise RuntimeError('No certificate: current dual does not separate all columns')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('short');ap.add_argument('--replay',action='store_true');a=ap.parse_args();run(a.short,a.replay)
