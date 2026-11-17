import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import json,subprocess,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
from geom import AG
from lift11 import atomic
ag=AG(2,5);inp=Path('out/union14_enum.in')
inp.write_text('\n'.join(' '.join(str(sum(1<<p for p in ag.lines[l][1])) for l in ag.lines_by_dir[d]) for d in range(6)))
subprocess.run(['g++','-O3','-std=c++17','union14_enum.cpp','-o','out/union14_enum'],check=True)
subprocess.run(['out/union14_enum',str(inp),'out/union14_enum.tsv'],check=True)
best={}
for line in Path('out/union14_enum.tsv').read_text().splitlines():
 idx,k,mask=map(int,line.split());b=[]
 for _ in range(6):b.append(idx%6);idx//=6
 b=tuple(reversed(b));S={p for p in range(25) if mask>>p&1}
 counts=tuple(sorted(sum(ag.lines[l][1]<=S for l in ag.lines_by_dir[d]) for d in range(6)))
 assert len(S)==k and counts==b;best[b]=k
F={a:min(k for b,k in best.items() if sum(b)>=a) for a in range(6,31)}
Ft={a:min(k for b,k in best.items() if all(x<=y for x,y in zip(a,b))) for a in combinations_with_replacement(range(1,6),6)}
old=json.load(open('out/union_tables_q5.json'));mismatches=[]
for a,v in F.items():
 if old['F'][str(a)]!=v:mismatches.append(dict(kind='F',key=a,old=old['F'][str(a)],exact=v))
for a,v in Ft.items():
 key=','.join(map(str,a))
 if old['Ft'][key]!=v:mismatches.append(dict(kind='Ft',key=key,old=old['Ft'][key],exact=v))
unsafe=[m for m in mismatches if m['old']>m['exact']]
atomic('out/union14_audit.json',dict(status='EXACT_UNION14_SAFE' if not unsafe else 'UNSAFE_OLD_TABLE',enumerated=2**25,F_entries=len(F),Ft_entries=len(Ft),mismatches=mismatches,unsafe=unsafe,old_sha256=hashlib.sha256(Path('out/union_tables_q5.json').read_bytes()).hexdigest()))
atomic('out/union14_exact.json',dict(F=F,Ft={','.join(map(str,k)):v for k,v in Ft.items()}))
print('UNION14_AUDIT',len(mismatches),'mismatches',len(unsafe),'unsafe',flush=True)
assert not unsafe,unsafe
