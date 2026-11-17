"""Replay inherited reduction certificates with their original model environments."""
import json,os,subprocess,hashlib
from pathlib import Path
from lift11 import atomic
base=dict(os.environ)
for k in ('PC6_SECTIONS','PC6_DMH','PC6_UNMARKED','PC6_EXACT_K'):base[k]='0'
base['PC6_F6_RMIN']='2';base.pop('PC6_MFIX',None)
jobs=[('out/lp_probe_rmax4.json',['--rmax','4'],base)]
for cfg in ('5,21','5,18','5,17','5,16','5,15','6,25','6,21','6,19','6,18','6,17'):
 jobs.append(('out/lp_probe_cfg_'+cfg.replace(',','_')+'.json',['--config',cfg],base))
sec=dict(base,PC6_SECTIONS='1',PC6_EXACT_K='4',PC6_F6_RMIN='5')
for cfg in ('5,15','5,16','6,17'):jobs.append(('out/lp_probe_f7_cfg_'+cfg.replace(',','_')+'.json',['--config',cfg],sec))
results=[]
for i,(fn,args,env) in enumerate(jobs):
 print('REDUCTION14_START',fn,flush=True)
 log=Path(f'out/reduction14_{i:02}.log')
 with log.open('w') as out:r=subprocess.run(['.venv/bin/python','-u','cert5_check.py',fn]+args,env=env,stdout=out,stderr=subprocess.STDOUT)
 text=log.read_text();assert r.returncode==0 and 'ALL CERTIFICATES VERIFIED' in text and 'FAILED' not in text,(fn,text[-2000:])
 rec=dict(file=fn,log=str(log),verified=text.count('-> VERIFIED'));results.append(rec)
 atomic('out/reduction14.json',dict(status='RUNNING',results=results));print('REDUCTION14_VERIFIED',rec,flush=True)
atomic('out/reduction14.json',dict(status='EXACT_REDUCTION14_VERIFIED',results=results,degree_range_before_LP=[4,11],degree_range_after_V_VI=[5,7],configurations=[[5,15],[6,17]],projection_table_sha256=hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest(),manifest_sha256=hashlib.sha256(Path('out/reduction14_manifest.json').read_bytes()).hexdigest()))
with open('out/reduction14.json','a') as f:f.write('\n')
print('EXACT_REDUCTION14_VERIFIED',flush=True)
