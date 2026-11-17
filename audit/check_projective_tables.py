"""AUDIT: compare the projective bound files the D=7 rows read (kakeya-fq/out/{point14,rich14,nonfano14}_projective_q5.json)
with the fresh exhaustive enumeration in out/projective_result.json (check_projective.py)."""
import json
mine = json.load(open('out/projective_result.json')); O = '../../kakeya-fq/out/'
p = json.load(open(O + 'point14_projective_q5.json')); r = json.load(open(O + 'rich14_projective_q5.json')); nf = json.load(open(O + 'nonfano14_projective_q5.json'))
ok = True
for j in range(3, 8):
    m = mine[str(j)]
    pj = p['records'][j]; rj = r['records'][j]; nj = nf['records'][j]
    row = dict(j=j, T_theirs=(pj['minimum'], pj['maximum']), T_mine=m['T'], rich_cap_theirs=(rj['cap_minimum'], rj['cap_maximum']), rich_cap_mine=m['rich_cap4'], g_theirs=(nj['minimum'], nj['maximum']), g_mine=m['g'])
    good = tuple(row['T_theirs']) == tuple(m['T']) and tuple(row['rich_cap_theirs']) == tuple(m['rich_cap4']) and tuple(row['g_theirs']) == tuple(m['g'])
    ok &= good; print(row, 'OK' if good else 'MISMATCH')
print('projective tables agree with fresh enumeration:', ok)
