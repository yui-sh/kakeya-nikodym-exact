"""IX/X escalation for a prefix-LP survivor, with canonical plane-column cache.
Run under run_coverage11.sh's environment. Unfinished/failing leaves remain open.
"""
import argparse,json,hashlib
from pathlib import Path
from lift11 import key,atomic
from plane10 import PlaneOracle,freeze
from coverage11 import LEDGER,nid,refresh,frontier

class CachedPlaneOracle(PlaneOracle):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.cache_dir=Path('out/coverage11_columns');self.cache_dir.mkdir(exist_ok=True)
        self.geo_lookup={key(g):i for i,g in enumerate(self.geos)}
    def save_columns(self,priced):
        grouped={}
        for _,col in priced:
            _,i,gs,gi,ds=col;t=self.types[i];tk=key(t)
            dm=sorted(self.info['_dms'][self.info['_Y'][g][1]] for g in gs)
            rec=dict(type=t,dms=dm,geometry=self.geos[gi],degrees=ds)
            grouped.setdefault(tk,{})[key(dm)]=rec
        for tk,new in grouped.items():
            path=self.cache_dir/(hashlib.sha256(tk.encode()).hexdigest()+'.json')
            old=json.load(path.open()) if path.exists() else {}
            old.update(new);atomic(path,old)
    def price(self,y,first=False):
        # The cache only seeds the master. Exhaustive pricing ALWAYS runs for a dual.
        result=super().price(y,first)
        self.save_columns(result)
        if first:
            lookup={(key(a),key(self.info['_dms'][di])):g for g,(a,di) in enumerate(self.info['_Y'])}
            for i in self.nr:
                t=self.types[i];path=self.cache_dir/(hashlib.sha256(key(t).encode()).hexdigest()+'.json')
                if not path.exists():continue
                for rec in json.load(path.open()).values():
                    assert freeze(rec['type'])==t
                    gs=[lookup.get((key(t[:3]),key(dm))) for dm in rec['dms']]
                    if any(g is None for g in gs):continue
                    gi=self.geo_lookup[key(rec['geometry'])]
                    col=['plane',i,gs,gi,rec['degrees']];self.verify(col);result.append((0,col))
        return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('node',nargs='?');ap.add_argument('--cache-x',action='store_true');ap.add_argument('--replay',action='store_true');a=ap.parse_args()
    if a.cache_x:
        from plane10 import build
        for cfg in ('515','617'):
            record=json.load(open('out/plane10_'+cfg+'.json'))
            lp,info,types,grow,zrow,nr,gr=build(tuple(record['config']),tuple(record['m']))
            oracle=CachedPlaneOracle(types,info,nr,gr,'cache11_'+cfg)
            path=Path('out/plane10_'+cfg+'_columns.json')
            if not path.exists():path=Path('out/plane10_'+cfg+'_checkpoint_columns.json')
            cols=[c for c in json.load(path.open()) if c[0]=='plane']
            for c in cols:oracle.verify(c)
            oracle.save_columns([(0,c) for c in cols])
            print('CANONICAL_CACHE_IMPORTED',cfg,len(cols),flush=True)
        return
    assert a.node
    import fcntl
    lock=open('out/coverage11.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    ledger=json.load(LEDGER.open());node=ledger['nodes'][a.node]
    assert len(node['prefix'])==5
    if not a.replay:assert node['status']=='open', 'leaf is already certified'
    assert any(idx==a.node for idx,_ in frontier(ledger)), 'not an active leaf'
    assert node.get('prefix_probe_done'), 'prefix cuts must be attempted before joint pricing'
    m=tuple(node['prefix']);cfg=(5,15) if node['config']=='515' else (6,17)
    short='xi_'+node['config']+'_'+'_'.join(map(str,m))
    import pencil9,plane10_solve,plane10_cert,plane10_check
    for module in (plane10_solve,plane10_cert,plane10_check):module.PlaneOracle=CachedPlaneOracle
    if not a.replay:
        pencil9.solve(cfg,m,'pencil9_'+short,200)
        plane10_solve.solve(short)
    result=json.load(open('out/plane10_'+short+'.json'))
    if result['status']=='LP-certified infeasible':
        plane10_cert.run(short,replay=True)
        node.update(status='inherited-LP-certified',joint_tag=short,certificate='out/plane10_'+short+'_farkas.json')
    elif result['status']=='computed_feasible_pending_exact':
        if not a.replay:plane10_check.check(short)
        plane10_check.check(short,replay=True)
        node.update(status='exact-primal',joint_tag=short,primal='out/plane10_'+short+'_rational_primal.json')
    else:node.update(status='open',reason='joint '+result['status'],joint_tag=short)
    refresh(ledger);print('JOINT11_CHECKPOINT',ledger['counts'],flush=True)
if __name__=='__main__':main()
