"""Solver-free exact replay of joint plane and completed-pencil primal witnesses."""
import argparse,json
from fractions import Fraction as F
import numpy as np
from plane10 import build,counts,PlaneOracle
from pencil9 import compatible

def check(short,replay=False):
    tag='plane10_'+short;r=json.load(open(f'out/{tag}.json'));m=tuple(r['m']);config=tuple(r['config'])
    lp,info,types,grow,zrow,nr,gr=build(config,m);po=PlaneOracle(types,info,nr,gr,tag)
    cols=json.load(open(f'out/{tag}_columns.json'));x=np.load(f'out/{tag}_primal.npz')['x'];base=len(lp.names)
    assert len(x)==base+len(cols);sparse=[]
    for col in cols:
        if col[0]=='plane':po.verify(col)
        else:
            _,z,gs=col;di,wi=info['_Z'][z];dm=info['_dms'][di];dv=tuple(dm.count(j) for j in range(1,len(m)+1));keys=[info['_Y'][g][0] for g in gs]
            assert all(info['_Y'][g][1]==di for g in gs) and compatible(m,dv,keys,5)
            assert sorted(next(ii for ii,it in enumerate(info['_items']) if key[0]==it[0] and key[1] in it[3]) for key in keys)==list(info['_wt'][wi])
        sparse.append(counts(col,grow,zrow,nr,gr))
    def verify(v):
        assert min(v)>=0
        left=[sum((c*v[j] for j,c in row[0].items()),F(0)) for row in lp.rows]
        for a,col in zip(v[base:],sparse):
            if a:
                for i,c in col.items():left[i]+=a*c
        bad=[i for i,row in enumerate(lp.rows) if (left[i]!=row[1] if row[2]=='=' else left[i]>row[1])]
        return bad
    if replay:
        v=[F(0)]*len(x)
        for i,a in json.load(open(f'out/{tag}_rational_primal.json')).items():v[int(i)]=F(a)
        assert not verify(v)
        print(tag,'SAVED_RATIONAL_VERIFIED',len(lp.rows),'rows',len(cols),'columns',flush=True);return
    support={i for i,a in enumerate(x) if a>1e-8};eqs=[{j:F(c) for j,c in row[0].items() if j in support} for row in lp.rows]
    lhs=np.array([sum(c*x[j] for j,c in row[0].items()) for row in lp.rows])
    for j,col in enumerate(sparse,base):
        for i,c in col.items():
            lhs[i]+=c*x[j]
            if j in support:eqs[i][j]=F(c)
    basis={}
    for i,row in enumerate(lp.rows):
        if row[2]!='=' and abs(lhs[i]-row[1])>1e-8:continue
        d=eqs[i].copy();b=F(row[1])
        while d:
            p=min(d)
            if p not in basis:
                s=d[p];basis[p]=({j:c/s for j,c in d.items()},b/s);break
            a=d[p];bd,bb=basis[p];b-=a*bb
            for j,c in bd.items():
                d[j]=d.get(j,F(0))-a*c
                if not d[j]:del d[j]
        else:assert b==0,(i,b)
    v=[F(0)]*len(x)
    for j in support:
        if j not in basis:v[j]=F(float(x[j])).limit_denominator(1000000)
    for p in sorted(basis,reverse=True):
        d,b=basis[p];v[p]=b-sum(c*v[j] for j,c in d.items() if j!=p)
    assert not verify(v)
    json.dump({str(i):str(a) for i,a in enumerate(v) if a},open(f'out/{tag}_rational_primal.json','w'),indent=1)
    report=dict(exact_rational_feasible=True,rows=len(lp.rows),columns=len(cols),support=len(support),rank=len(basis),max_denominator=max(a.denominator for a in v),positive_planes=sum(v[base+j]>0 for j,c in enumerate(cols) if c[0]=='plane'),positive_pencils=sum(v[base+j]>0 for j,c in enumerate(cols) if c[0]=='pencil'))
    json.dump(report,open(f'out/{tag}_check.json','w'),indent=2);print(tag,report,flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('tags',nargs='+');ap.add_argument('--rational-only',action='store_true');a=ap.parse_args()
    for tag in a.tags:check(tag,a.rational_only)
    print('P10_PRIMALS_CHECKED',flush=True)
