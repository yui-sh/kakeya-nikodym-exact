"""Revalidate all m-free base rows on four sets; cache the common n=53,D=6 model."""
import os,hashlib,argparse
from pathlib import Path
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import planecount6 as p6
import pc6_validate,f9_validate
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--D-cap',type=int);a=ap.parse_args();suffix=f'_D{a.D_cap}' if a.D_cap else ''
files=['solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json']
original=p6.build;cache={}
def build(q,n,D,*args,**kwargs):
    kwargs.pop('mfix',None)
    if a.D_cap:kwargs['max_types']=2000000
    key=(q,n,D,args,tuple(sorted(kwargs.items())))
    if key not in cache:
        cache.clear();print('VALIDATION14_BUILD',q,n,D,flush=True)
        cache[key]=original(q,n,D,*args,mfix=None,**kwargs)
    return cache[key]
p6.build=build
results=[]
for fn in files:
    assert pc6_validate.validate(fn,a.D_cap),fn
    f9_validate.validate(fn,a.D_cap)
    results.append(dict(file=fn,violations=0))
    atomic(f'out/validation14_base{suffix}_progress.json',dict(status='RUNNING',results=results))
atomic(f'out/validation14_base{suffix}.json',dict(status='ZERO_VIOLATIONS',mfree=True,degree_cap=a.D_cap,inactive_premise='m_D>=1 is skipped only when the known example has smaller actual maximum degree',projection_table_sha256=hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest(),families='F1-F9 including sections, exact pencil items, R4, F7, F8 and F9; Y3 not used',results=results))
print('VALIDATION14_BASE_ZERO_VIOLATIONS',flush=True)
