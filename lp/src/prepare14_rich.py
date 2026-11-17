"""Resumable field-specific point-cost preparation, without LP solvers."""
import json,os,subprocess,sys,time
from pathlib import Path
from lift11 import atomic
os.environ['TMPDIR']=str(Path('out').resolve())
p=Path('out/prepare14_rich.json')
state=json.loads(p.read_text()) if p.exists() else dict(completed=[])
stages=[('compile',['g++','-O3','-std=c++20','-pipe','pg14_rich_enum.cpp','-o','out/pg14_rich_enum']),
        ('projective',[sys.executable,'-u','pg14_rich.py']),
        ('prices',[sys.executable,'-u','pointrich14.py','--D','7']),
        ('conditional-validation',[sys.executable,'-u','-c','import pointrich14; pointrich14.validate(7)']),
        ('independent',[sys.executable,'-u','pointrich14_check.py','--D','7','--sample','3000'])]
for name,cmd in stages:
    if name in state['completed']:continue
    state.update(status='running',current=name,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));atomic(p,state)
    print('RICH14_PREPARE',name,flush=True)
    subprocess.run(cmd,check=True)
    state['completed'].append(name);atomic(p,state)
state.update(status='ready',current=None);atomic(p,state)
print('RICH14_READY',flush=True)
