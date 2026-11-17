"""D-specific IX/X escalation and exact replay; isolated ledger and canonical seeds."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import argparse,json,hashlib,fcntl,time
from pathlib import Path
from lift11 import atomic,key
from geometry14 import geometries,freeze
from coverage14 import path,frontier,refresh
import plane10,pencil9,plane10_solve,plane10_cert,plane10_check
_BASE_PENCIL_BUILD=pencil9.build

def configure(D,strong=False):
    from maximal14 import guard
    from pairbounds14 import add_cuts as add_pair_bounds
    def build_with_exact_global_pairs(config,m):
        cfg='515' if tuple(config)==(5,15) else '617'
        assert tuple(config) in ((5,15),(6,17)) and len(m)==D
        guard(D,cfg)
        original=pencil9.p6.hist_types;original_build=pencil9.p6.build
        def maximal_types(*args,**kwargs):
            from budget14 import polynomial
            return [t for t in original(*args,**kwargs) if (t[0]!=config[0] or t[3]==config[1])
                    and (not strong or polynomial(t[0],t[1])<=541)]
        def large_build(*args,**kwargs):
            kwargs['max_types']=2000000
            return original_build(*args,**kwargs)
        pencil9.p6.hist_types=maximal_types;pencil9.p6.build=large_build
        try:result=_BASE_PENCIL_BUILD(config,m)
        finally:pencil9.p6.hist_types=original;pencil9.p6.build=original_build
        add_pair_bounds(result[0],D,tuple(m),tuple(m))
        if strong:
            assert D==7
            lp,info,types=result[:3]
            from lift14 import add_cuts as additive
            from nonlinear14 import add_cuts as pairs
            from triple14 import add_cuts as triples
            from point14 import add_cuts as point_pairs
            from pointtriple14 import add_cuts as point_triples
            from pointnonfano14 import add_cuts as nonfano
            from star14 import add_cuts as star
            from roundbudget14 import add_cuts as rounded
            from degreehull14 import add_cuts as hull
            additive(lp,info,types,D,'all',True)
            pairs(lp,info,types,D,'all',True);triples(lp,info,types,D,'all',True)
            point_pairs(lp,info,types,D);point_triples(lp,info,types,D)
            star(lp,types,D,cfg);rounded(lp,types);nonfano(lp,info,types,D)
            from pointrich14 import add_cuts as rich
            from extra14 import add_cuts as extra
            rich(lp,info,types,D)
            if cfg=='617':extra(lp,info,types,D,True)
            from pointmix14 import add_cuts as point_mix
            from rounddegree14 import add_cuts as round_degree
            point_mix(lp,info,types,D)
            hull(lp,D,[tuple(m)])
            round_degree(lp,info,types,D,[tuple(m)])
        return result
    pencil9.build=build_with_exact_global_pairs
    plane10.geometries=lambda q=5:geometries(D,q)
    class Oracle14(plane10.PlaneOracle):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            assert all(len(t[2])==D for t in self.types)
            self.cache_dir=Path(f'out/columns14_D{D}');self.cache_dir.mkdir(exist_ok=True)
            self.geo_lookup={key(g):i for i,g in enumerate(self.geos)}
        def price(self,y,first=False):
            # No D=5 materialized catalogue is consulted. Exhaustive DFS always prices the current universe.
            result=self.dfs_price(y,first)
            grouped={}
            for _,col in result:
                _,i,gs,gi,ds=col;t=self.types[i];tk=key(t)
                dms=sorted(self.info['_dms'][self.info['_Y'][g][1]] for g in gs)
                grouped.setdefault(tk,{})[key(dms)]=dict(type=t,dms=dms,geometry=self.geos[gi],degrees=ds)
            for tk,new in grouped.items():
                p=self.cache_dir/(hashlib.sha256(tk.encode()).hexdigest()+'.json');old=json.load(p.open()) if p.exists() else {};old.update(new);atomic(p,old)
            if first:
                lookup={(key(a),key(self.info['_dms'][di])):g for g,(a,di) in enumerate(self.info['_Y'])}
                for i in self.nr:
                    t=self.types[i];p=self.cache_dir/(hashlib.sha256(key(t).encode()).hexdigest()+'.json')
                    if not p.exists():continue
                    for rec in json.load(p.open()).values():
                        assert freeze(rec['type'])==t
                        gs=[lookup.get((key(t[:3]),key(dm))) for dm in rec['dms']]
                        if any(g is None for g in gs):continue
                        gi=self.geo_lookup[key(rec['geometry'])];col=['plane',i,gs,gi,rec['degrees']]
                        self.verify(col);result.append((0,col))
            return result
    for module in (plane10_solve,plane10_cert,plane10_check):module.PlaneOracle=Oracle14

def strong_manifest(D,cfg):
    assert D==7
    from maximal14 import digest
    from budget14 import guard as budget_guard
    from star14 import guard as star_guard
    from roundbudget14 import guard as rounded_guard
    files=['out/validation14_D7.json','out/validation14_base_D7.json']
    files += [f'out/{name}_D7_validation.json' for name in
              ('nonlinear14','triple14','point14','pointtriple14','pointnonfano14','pointrich14','extra14','pointmix14','rounddegree14','pairbounds14','degreehull14')]
    return dict(version='XIV_D7_all_validated_v3',budget=budget_guard(),star=star_guard(D,cfg),
                rounded=rounded_guard(),validations={f:digest(f) for f in files})

def replay_node(D,node):
    from maximal14 import guard,digest
    model=node['joint_model'];assert model['maximal_only'] and model['pair_bounds']
    assert model['reduction_sha256']==guard(D,node['config']) and model['projection_table_sha256']==digest('out/union_tables_q5.json')
    strong=model.get('strong_D7',False)
    if strong:assert model['strong_manifest']==strong_manifest(D,node['config'])
    configure(D,strong);tag=node['joint_tag'];m=node['prefix'];cfg=[5,15] if node['config']=='515' else [6,17]
    for fn in ('out/pencil9_'+tag+'.json','out/plane10_'+tag+'.json'):
        data=json.load(open(fn));assert data['m']==m and data['config']==cfg and len(m)==D
    if node['status']=='joint-LP-certified':plane10_cert.run(tag,replay=True)
    elif node['status']=='exact-primal':plane10_check.check(tag,replay=True)
    else:raise ValueError(node['status'])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True,choices=[6,7]);ap.add_argument('node');a=ap.parse_args();D=a.D
    lock=open(f'out/coverage14_D{D}.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    ledger=json.load(path(D).open());n=ledger['nodes'][a.node]
    assert len(n['prefix'])==D and n['status']=='open' and n.get('joint_pending')
    assert any(idx==a.node for idx,_ in frontier(ledger))
    strong=(D==7)
    tag=f"xiv_D{D}_{n['config']}_"+'_'.join(map(str,n['prefix']))+('_strongv3' if strong else '');n['joint_tag']=tag
    from maximal14 import guard,digest
    n['joint_model']=dict(maximal_only=True,pair_bounds=True,reduction_sha256=guard(D,n['config']),projection_table_sha256=digest('out/union_tables_q5.json'))
    if strong:n['joint_model'].update(strong_D7=True,strong_manifest=strong_manifest(D,n['config']))
    n.setdefault('joint_attempts',[]).append(dict(started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
    if strong:
        ledger['rich_used']=True;ledger['point_mix_used']=True;ledger['round_degree_used']=True
        if n['config']=='617':ledger['extra_used']=True
    refresh(ledger);configure(D,strong)
    cfg=(5,15) if n['config']=='515' else (6,17);m=tuple(n['prefix'])
    pencil9.solve(cfg,m,'pencil9_'+tag,200);plane10_solve.solve(tag)
    result=json.load(open('out/plane10_'+tag+'.json'))
    if result['status']=='LP-certified infeasible':
        plane10_cert.run(tag,replay=True)
        n.update(status='joint-LP-certified',certificate='out/plane10_'+tag+'_farkas.json',joint_pending=False)
    elif result['status']=='computed_feasible_pending_exact':
        plane10_check.check(tag);plane10_check.check(tag,replay=True)
        n.update(status='exact-primal',primal='out/plane10_'+tag+'_rational_primal.json',joint_pending=False)
    else:n.update(status='open',reason='joint '+result['status'])
    refresh(ledger);print('JOINT14_CHECKPOINT',D,ledger['counts'],flush=True)
if __name__=='__main__':main()
