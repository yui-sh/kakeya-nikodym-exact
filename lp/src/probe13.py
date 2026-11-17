"""Independent exact-certifying prefix probe; never writes the coverage ledger."""
import argparse,json,time,datetime
from coverage11 import build,node_lp,Master,verify
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--prefix',default='');ap.add_argument('--scope',default='pair',choices=['core','all','pair']);ap.add_argument('--seconds',type=int,default=900);ap.add_argument('--solver',default='ipm',choices=['ipm','simplex']);ap.add_argument('--tag',required=True);a=ap.parse_args()
prefix=tuple(map(int,a.prefix.split(','))) if a.prefix else ()
assert time.time()+a.seconds+180<datetime.datetime(2026,11,4,10,0,tzinfo=datetime.timezone.utc).timestamp()
lp,info,types=build(a.scope,'617',prefix);lp=node_lp(lp,types,'617',prefix)
master=Master(lp,'phase1');master.solver.setOptionValue('solver',a.solver)
if a.solver=='ipm':master.solver.setOptionValue('run_crossover','off')
print('PROBE13_START',a.tag,prefix,a.scope,a.solver,flush=True)
st,y,detail=master.solve(lp,a.seconds);detail.update(scope=a.scope,solver=a.solver,seconds=a.seconds)
if st=='LP-certified':
    path='out/coverage13_cert_'+a.tag+'.json'
    atomic(path,dict(config='617',prefix=list(prefix),model_kind='prefix-filtered',cut_scope=a.scope,multipliers=y,detail=detail))
    verify(lp,json.load(open(path))['multipliers'])
result=dict(status=st,detail=detail,prefix=list(prefix),config='617',cut_scope=a.scope)
if st=='numerical-feasible':
    import numpy as np
    np.savez_compressed('out/coverage13_probe_'+a.tag+'.npz',x=y)
    result['note']='numerical prefix solution only; no exact-primal survival claim'
atomic('out/coverage13_probe_'+a.tag+'.json',result)
print('PROBE13_RESULT',a.tag,result,flush=True)
