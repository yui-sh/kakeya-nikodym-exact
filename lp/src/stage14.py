"""Bounded preparation/validation/root sequence; all stage logs and status are resumable."""
import os,time,subprocess,json,datetime
from pathlib import Path
from lift11 import atomic
D=6
while 'PREPARE14_COMPLETE' not in Path('out/prepare14_D6.log').read_text():
    assert time.time()<datetime.datetime(2026,11,4,8,0,tzinfo=datetime.timezone.utc).timestamp()
    time.sleep(10)
steps=[('base_validation',['validate14_base.py'],'out/validation14_base.log'),('D6_roots',['coverage14.py','--D','6','--roots-only','--nodes','2','--seconds','1800'],'out/coverage14_D6_roots.log')]
for name,args,fn in steps:
    atomic('out/stage14.json',dict(stage=name,status='running',log=fn));print('STAGE14_START',name,flush=True)
    with open(fn,'w') as f:r=subprocess.run(['bash','run14.sh']+args,stdout=f,stderr=subprocess.STDOUT)
    atomic('out/stage14.json',dict(stage=name,status='complete' if r.returncode==0 else 'failed',returncode=r.returncode,log=fn))
    assert r.returncode==0,(name,r.returncode)
    print('STAGE14_COMPLETE',name,flush=True)
