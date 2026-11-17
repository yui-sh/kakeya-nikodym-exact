"""Resumable D-specific catalogue, pricing, and validation preparation."""
import argparse,json,subprocess,sys,time
from pathlib import Path
from lift11 import atomic
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,choices=[6,7],required=True);a=ap.parse_args();D=a.D
if D==7:
    lower=json.load(open('out/coverage14_D6.json'))
    assert lower['counts']==dict(pruned_certified=lower['total'],surviving=0,open=0)
steps=[('lift14.py',[]),('nonlinear14.py',[]),('triple14.py',[]),('pairbounds14.py',[]),('degreehull14.py',[]),('nonlinear14_check.py',[]),('triple14_check.py',[]),('dp14.py',['--sample','1000000000'])]
path=Path(f'out/prepare14_D{D}.json')
state=json.load(path.open()) if path.exists() else dict(D=D,completed=[],status='running')
assert state['D']==D
for script,extra in steps:
    if script in state['completed']:continue
    state.update(status='running',current=script,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));atomic(path,state)
    print('PREPARE14_STAGE_START',D,script,flush=True)
    result=subprocess.run([sys.executable,'-u',script,'--D',str(D)]+extra)
    if result.returncode:
        state.update(status='stage_failed',returncode=result.returncode);atomic(path,state);sys.exit(result.returncode)
    state['completed'].append(script);state.pop('current',None);atomic(path,state)
    print('PREPARE14_STAGE_COMPLETE',D,script,flush=True)
state.update(status='COMPLETE',utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()));atomic(path,state)
print('PREPARE14_ALL_COMPLETE',D,flush=True)
