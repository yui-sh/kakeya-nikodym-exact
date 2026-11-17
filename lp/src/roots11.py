"""Independent root probes; no coverage ledger writes. Import only after exact replay."""
import json,time,datetime
from coverage11 import build,node_lp,Master,verify
from lift11 import atomic
for cfg in ('515','617'):
    if time.time()+650>=datetime.datetime(2026,11,4,11,30,tzinfo=datetime.timezone.utc).timestamp():break
    lp,info,types=build('core',cfg,[]);lp=node_lp(lp,types,cfg,[])
    master=Master(lp,'phase1');print('ROOT_START',cfg,flush=True)
    st,y,detail=master.solve(lp,600)
    if st=='LP-certified':
        path='out/coverage11_cert_'+cfg+'_root.json'
        atomic(path,dict(config=cfg,prefix=[],model_kind='prefix-filtered',multipliers=y,detail=detail))
        verify(lp,json.load(open(path))['multipliers'])
    atomic('out/coverage11_root_'+cfg+'_result.json',dict(status=st,detail=detail))
    print('ROOT_RESULT',cfg,st,detail,flush=True)
    del master,lp,info,types
