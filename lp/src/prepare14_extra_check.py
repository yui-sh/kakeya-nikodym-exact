"""Wait for primary prices, then independently check; explicit resumable status."""
import time,json,subprocess,sys
from pathlib import Path
from lift11 import atomic
p=Path('out/extra14_D7_prices.json')
while not p.exists():
    atomic('out/prepare14_extra_check.json',dict(status='waiting for primary full-type prices'))
    time.sleep(20)
atomic('out/prepare14_extra_check.json',dict(status='checking'))
subprocess.run([sys.executable,'-u','extra14_check.py','--D','7','--sample','3000'],check=True)
atomic('out/prepare14_extra_check.json',dict(status='ready'))
