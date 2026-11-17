"""Independent coordinate-bitmask replay of degree-independent q=5 inputs."""
import json,hashlib,itertools
from pathlib import Path
from lift11 import atomic
q=5;allpoints=(1<<25)-1
lines=[]
for slope in range(5):
    lines.append([sum(1<<(5*x+(slope*x+b)%5) for x in range(5)) for b in range(5)])
lines.append([sum(1<<(5*b+y) for y in range(5)) for b in range(5)])
seen={r:set() for r in range(7)};systems=0
for choice in itertools.product(range(6),repeat=6):
    selected=[lines[d][b] for d,b in enumerate(choice) if b<5]
    levels=[allpoints]
    for mask in selected:
        levels.append(0)
        for b in range(len(levels)-1,0,-1):levels[b]=(levels[b]&~mask)|(levels[b-1]&mask)
        levels[0]&=~mask
    assert sum(v.bit_count() for v in levels)==25
    profile=(25-levels[0].bit_count(),)+tuple(v.bit_count() for v in levels[1:])
    seen[len(selected)].add(profile);systems+=1
assert systems==46656
cached=json.load(open('out/inside_profiles_q5.json'))
assert all(seen[r]=={tuple(p) for p in cached[str(r)]} for r in seen)
minima={r:min(p[0] for p in values) for r,values in seen.items()}
f=json.load(open('out/plane_f.json'))['5']
assert minima=={int(r):value for r,value in f.items()}
section=json.load(open('out/degree_independent14_tables.json'))
assert hashlib.sha256(Path(section['file']).read_bytes()).hexdigest()==section['sha256']
inputs=['out/inside_profiles_q5.json','out/plane_f.json',section['file']]
r=dict(status='EXACT_GEOMETRIC_INPUTS14_REPLAYED',q=5,line_systems=systems,
       profile_counts={r:len(v) for r,v in seen.items()},planar_minima=minima,
       method='direct slope/intercept coordinates and bitmask occupancy levels; no inside_profiles generator calls',
       degree_independence='These tables record only geometric unions, multiplicities, and section signatures keyed by q/r; no ambient D, m, h, or degree assignment enters them.',
       input_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs})
atomic('out/geometric_inputs14_check.json',r);print('EXACT_GEOMETRIC_INPUTS14_REPLAYED',r,flush=True)
