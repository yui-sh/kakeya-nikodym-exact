"""Validate every materialized column and compare full-catalogue minima with DFS."""
import json,hashlib
from pathlib import Path
import numpy as np
from plane10 import build,PlaneOracle

def check(short,config,m):
    lp,info,types,gr,zr,nr,pgr=build(config,m);po=PlaneOracle(types,info,nr,pgr,'plane10_catalogcheck_'+short)
    dtype=np.dtype([('type','<u4'),('geo','<u4'),('gids','<u2',(6,)),('degrees','<u8')]);path=Path(f'out/plane10_catalog_{short}.bin');recs=np.memmap(path,dtype=dtype,mode='r')
    assert path.stat().st_size%dtype.itemsize==0
    # Every witness checked via fast independently derived point/line masks and histograms.
    ag=po.ag;unique=set()
    for i in nr:
        rs=recs[recs['type']==i];t=types[i];r,k,h,u,sec,e=t
        for gi in set(map(int,rs['geo'])):
            sub=rs[rs['geo']==gi];geo=po.geos[gi];pts=geo['points'];lines=geo['lines']
            masks=np.array([sum(1<<j for j,l in enumerate(lines) if p in ag.lines[l][1]) for p in pts])
            ds=((sub['degrees'][:,None]>>(3*np.arange(k,dtype=np.uint64)))&7).astype(np.int64)
            assert all(np.all(np.sum(ds==j,axis=1)==h[j-1]) for j in range(1,6))
            assert np.all(ds>=np.array([int(z).bit_count() for z in masks])) and np.all((ds>=1)&(ds<=5))
            # Encode sorted line degrees as base-6 integers; compare joint multisets.
            actual=np.stack([np.sort(ds[:,masks>>j&1!=0],axis=1)@ (6**np.arange(5)) for j in range(r)],axis=1);actual.sort(axis=1)
            expected=np.array([sum(d*6**j for j,d in enumerate(info['_dms'][info['_Y'][int(g)][1]])) for g in sub['gids'][:,:r].flat]).reshape(-1,r);expected.sort(axis=1)
            assert np.array_equal(actual,expected)
            assert all(info['_Y'][int(g)][0]==t[:3] for g in set(map(int,sub['gids'][:,:r].flat)))
            # The scalar verifier also recomputes the fixed geometry signature.
            if len(sub):po.verify(['plane',i,list(map(int,sub[0]['gids'][:r])),gi,list(map(int,ds[0]))])
        unique.update((i,tuple(map(int,row['gids']))) for row in rs)
    rng=np.random.default_rng(110)
    for trial in range(3):
        y=rng.integers(-1000,1001,size=len(lp.rows)).tolist()
        dfs=po.dfs_price(y);cat=po.catalog_price(y,path)
        def mins(priced):
            ans={}
            for c,col in priced:ans[col[1]]=min(c,ans.get(col[1],c))
            return ans
        assert mins(dfs)==mins(cat)
    result=dict(records=len(recs),distinct_type_joint_multisets=len(unique),types=len(nr),jobs=len(po.jobs),violations=0,pricing_comparisons=3,model_signature=po.catalog_signature,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    json.dump(result,open(f'out/plane10_catalog_{short}_check.json','w'),indent=2);print(short,result,flush=True)
if __name__=='__main__':
    check('617',(6,17),(0,3,27,12,4));check('515',(5,15),(0,4,24,15,3));print('CATALOGS_CHECKED',flush=True)
