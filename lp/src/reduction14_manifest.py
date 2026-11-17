"""Check that inherited certificates tile every excluded degree/configuration case."""
import json,hashlib
from pathlib import Path
from lift11 import atomic
from geometry14 import geometries
required={'out/lp_probe_rmax4.json':list(range(4,12))}
for r,u in [(5,21),(5,18),(5,17),(5,16),(5,15),(6,25),(6,21),(6,19),(6,18),(6,17)]:
 required[f'out/lp_probe_cfg_{r}_{u}.json']=list(range(8,12)) if (r,u) in [(5,16),(5,15),(6,17)] else list(range(4,12))
for r,u in [(5,15),(5,16),(6,17)]:required[f'out/lp_probe_f7_cfg_{r}_{u}.json']=[4,5,6,7] if (r,u)==(5,16) else [4]
records=[]
for fn,Ds in required.items():
 data=json.load(open(fn))
 for D in Ds:assert data[f'46,{D}']['status']=='INFEASIBLE' and data[f'46,{D}'].get('cert'),(fn,D)
 records.append(dict(file=fn,required_D=Ds,sha256=hashlib.sha256(Path(fn).read_bytes()).hexdigest()))
orbits={}
for g in geometries(6):orbits[tuple(g['lines'])]=(len(g['lines']),g['signature'][3])
assert sorted(u for r,u in orbits.values() if r==5)==[15,16,17,17,18,21]
assert sorted(u for r,u in orbits.values() if r==6)==[17,18,19,19,21,25]
atomic('out/reduction14_manifest.json',dict(required_certificates=sum(map(len,required.values())),files=records,affine_orbits=[list(v) for v in orbits.values()],elementary_degree_range=[4,11],remaining_after_V_VI=dict(degrees=[5,6,7],configurations=[[5,15],[6,17]])))
print('REDUCTION14_REQUIRED_CASES_VERIFIED',sum(map(len,required.values())),'affine_orbits',len(orbits),flush=True)
