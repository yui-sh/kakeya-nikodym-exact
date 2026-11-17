"""Certified n=46 case restriction: every maximal-r plane has the surviving union."""
import json,hashlib
from pathlib import Path
from lift11 import atomic
from geom import AG
from structure import planes_of
from planecount_check import witness

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def guard(D,cfg):
    assert D in (6,7) and cfg in ('515','617')
    r=json.load(open('out/reduction14.json'));assert r['status']=='EXACT_REDUCTION14_VERIFIED'
    assert sum(x['verified'] for x in r['results'])==82
    assert r['projection_table_sha256']==digest('out/union_tables_q5.json')
    assert r['manifest_sha256']==digest('out/reduction14_manifest.json')
    m=json.load(open('out/reduction14_manifest.json'));assert m['required_certificates']==82
    for entry in m['files']:assert digest(entry['file'])==entry['sha256']
    v=json.load(open('out/maximal14_validation.json'));assert v['status']=='EXACT_MAXIMAL14_REDUCTION'
    assert v['reduction_sha256']==digest('out/reduction14.json')
    return digest('out/reduction14.json')

def prepare():
    results=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        d=json.load(open(fn));ag=AG(3,d['q']);S={ag.pidx[tuple(p)] for p in d['points']};W=witness(ag,S)
        degrees={p:sum(p in ag.lines[l][1] for l in W) for p in S};D=max(degrees.values());unions=[];rmax=0
        for _,_,ps in planes_of(ag):
            ls=[l for l in W if ag.lines[l][1]<=ps];r=len(ls)
            if r>rmax:rmax=r;unions=[]
            if r==rmax:unions.append(len(set().union(*(ag.lines[l][1] for l in ls))))
        active=(ag.q==5 and len(S)==46 and D in (6,7) and rmax in (5,6))
        violations=sum(u!={5:15,6:17}[rmax] for u in unions) if active else 0
        assert violations==0
        results.append(dict(file=fn,q=ag.q,n=len(S),D=D,rmax=rmax,premises_active=active,violations=violations))
    atomic('out/maximal14_validation.json',dict(status='EXACT_MAXIMAL14_REDUCTION',reduction_sha256=digest('out/reduction14.json'),results=results,known_set_warning='All four examples are outside n=46 premises; these conditional checks are vacuous. Validity is proved by the 82 replayed exclusion certificates, not empirically by these examples.',excluded_unions={'515':[16,17,18,21],'617':[18,19,21,25]},argument='Existence of any other maximal-r plane is one of the globally excluded configurations. Lower-r unions are not restricted.'))
    for D in (6,7):
        for cfg in ('515','617'):guard(D,cfg)
    print('EXACT_MAXIMAL14_REDUCTION',flush=True)
if __name__=='__main__':prepare()
