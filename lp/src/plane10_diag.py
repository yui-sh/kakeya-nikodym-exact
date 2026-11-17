"""Exact diagnostic summaries, no geometric nonexistence inference."""
from collections import defaultdict
from fractions import Fraction as F
import json
from pencil9 import build as oldbuild
from plane10 import build

def old(short):
    tag='pencil9_'+short;r=json.load(open(f'out/{tag}.json'));lp,info,ts,_,_=oldbuild(tuple(r['config']),tuple(r['m']))
    x={int(i):F(a) for i,a in json.load(open(f'out/{tag}_rational_primal.json')).items()};bad=[]
    for i,t in enumerate(ts):
        if t[0]!=5 or t[1]!=15 or not x.get(i):continue
        for j in range(1,6):
            mass=sum(x.get(lp.col[f'G{g}'],F(0))*info['_dms'][di].count(j) for g,(key,di) in enumerate(info['_Y']) if key==t[:3]);lo=t[2][j-1]*x[i];hi=2*lo
            if mass<lo or mass>hi:bad.append(dict(type=i,degree=j,mass=str(mass),lo=str(lo),hi=str(hi)))
    json.dump(bad,open('out/plane10_old_diag_'+short+'.json','w'),indent=2)
    print(short,'IX violation count',len(bad),flush=True)

def new(short):
    tag='plane10_'+short;r=json.load(open(f'out/{tag}.json'));lp,info,ts,gr,zr,nr,pgr=build(tuple(r['config']),tuple(r['m']))
    x={int(i):F(a) for i,a in json.load(open(f'out/{tag}_rational_primal.json')).items()};cols=json.load(open(f'out/{tag}_columns.json'));base=len(lp.names)
    planes=defaultdict(F);high=[]
    for i,t in enumerate(ts):
        if x.get(i):
            planes[str(t[:2])]+=x[i]
            if t[0]>=5:high.append(dict(type=i,label=t,value=str(x[i])))
    pc=[];wc=[]
    for j,c in enumerate(cols):
        a=x.get(base+j,F(0))
        if a:(pc if c[0]=='plane' else wc).append(dict(index=j,value=str(a),weight=float(a),column=c))
    report=dict(plane_population={k:str(v) for k,v in planes.items()},positive_high_types=high,plane_columns=pc,pencil_columns=wc,nonintegral_variables=sum(a.denominator!=1 for a in x.values()),positive_variables=len(x),plane_columns_below_one=sum(z['weight']<1 for z in pc))
    json.dump(report,open(f'out/{tag}_diag.json','w'),indent=2)
    print(short,'final diag', {k:v for k,v in report.items() if k not in ('positive_high_types','plane_columns','pencil_columns')},flush=True)
if __name__=='__main__':
    for short in ('617','515'):
        old(short)
        if __import__('pathlib').Path(f'out/plane10_{short}_rational_primal.json').exists():new(short)
