import json
from pathlib import Path
from geom import AG
from lift11 import atomic
first={int(a):int(b) for a,b,_ in (line.split() for line in Path('out/union14_enum.tsv').read_text().splitlines())}
second={int(a):int(b) for a,b,_ in (line.split() for line in Path('out/union14_gray.tsv').read_text().splitlines())}
assert first==second
ag=AG(2,5);target=(1,1,2,2,3,3);found=False
for line in Path('out/union14_enum.tsv').read_text().splitlines():
 idx,k,mask=map(int,line.split())
 if k!=21:continue
 S={p for p in range(25) if mask>>p&1};available=[[l for l in ag.lines_by_dir[d] if ag.lines[l][1]<=S] for d in range(6)]
 order=sorted(range(6),key=lambda d:len(available[d]))
 if any(len(available[d])<a for d,a in zip(order,target)):continue
 chosen=[l for d,a in zip(order,target) for l in available[d][:a]];U=set().union(*(ag.lines[l][1] for l in chosen));assert len(U)==21
 counts=[sum(ag.lines[l][0]==d for l in chosen) for d in range(6)];assert tuple(sorted(counts))==target
 atomic('out/union14_counterexample.json',dict(q=5,n=2,points=[ag.points[p] for p in sorted(U)],cardinality=21,selected_lines=chosen,direction_counts=counts,sorted_counts=target,invalid_old_lower_bound=22))
 found=True;print('EXACT_21_POINT_FT_COUNTEREXAMPLE',counts,flush=True);break
assert found
atomic('out/union14_double_check.json',dict(status='EXACT_TWO_ENUMERATORS_AGREE',subsets_each=2**25,admitted_each=123276,attained_sorted_profiles=len(first),fixed_order_minimum=22,unordered_minimum=21))
print('EXACT_TWO_ENUMERATORS_AGREE',len(first),flush=True)
