import subprocess,sys,json
from pathlib import Path
from lift11 import atomic
p=Path('out/prepare14_mix.json');state=json.load(p.open()) if p.exists() else dict(completed=[])
for name,cmd in [('prices',[sys.executable,'-u','pointmix14.py','--D','7']),('independent',[sys.executable,'-u','pointmix14_check.py','--D','7','--sample','3000'])]:
 if name in state['completed']:continue
 state.update(status='running',current=name);atomic(p,state);subprocess.run(cmd,check=True);state['completed'].append(name);atomic(p,state)
state.update(status='ready',current=None);atomic(p,state)
