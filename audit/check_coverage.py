"""AUDIT: the global case split. (1) raw degree-vector universes for D=5,6,7 (fresh enumeration) vs the ledgers;
(2) every raw vector is under a certified node of the ledger (tree check, fresh); (3) the 82 reduction certificates:
which (configuration, D) each covers and which side rows (config present / rmax / forbid) carry nonzero multipliers;
(4) the eight root certificates' options. Reads lp/certs/*.json only."""
import json, itertools, math
from collections import Counter
C = '../lp/certs/'
def raw(D):
    out = []
    def rec(pre, left_n, left_s):
        j = len(pre) + 1
        if j == D:
            if left_s == D * left_n and (left_n >= 1): out.append(pre + (left_n,))
            return
        for h in range(left_n + 1):
            rn, rs = left_n - h, left_s - j * h
            if rs < (j + 1) * rn or rs > D * rn: continue
            rec(pre + (h,), rn, rs)
    rec((), 46, 155); return out
counts = {D: len(raw(D)) for D in (4, 5, 6, 7, 8)}
print('raw degree vectors (sum m=46, sum j m_j=155, m_D>=1):', counts, ' expected 5:2599 6:25622 7:137209')
for D, fn in ((6, 'raw14_D6.json'), (7, 'raw14_D7.json')):
    their = json.load(open(C + fn)); their = their['vectors'] if isinstance(their, dict) else their
    print(f'D={D} ledger raw file: {len(their)} vectors; equal to fresh set:', set(map(tuple, their)) == set(raw(D)))
# ledgers: tree coverage
def tree_check(led, D, cfgs):
    nodes = led['nodes']; R = set(raw(D)); covered = Counter(); ok = True
    def desc(p): return [m for m in R if m[:len(p)] == tuple(p)]
    for cfg in cfgs:
        root = nodes[cfg + ':']
        def walk(idx):
            node = nodes[idx]
            st = node['status']
            if st in ('LP-certified', 'joint-LP-certified'):
                for m in desc(node['prefix']): covered[(cfg, m)] += 1
                return
            if st == 'exact-primal': raise AssertionError('survivor ' + idx)
            if not node['children']: raise AssertionError('open leaf ' + idx + ' status ' + st)
            kids = [nodes[k]['prefix'] for k in node['children']]
            # children must partition the parent's raw vectors
            assert all(len(k) == len(node['prefix']) + 1 for k in kids)
            assert set(map(tuple, (m[:len(node['prefix']) + 1] for m in desc(node['prefix'])))) == set(map(tuple, kids)), idx
            for k in node['children']: walk(k)
        walk(cfg + ':')
    allc = {(cfg, m) for cfg in cfgs for m in R}
    return len(allc), sum(1 for k in allc if covered[k] >= 1), [k for k in allc if covered[k] == 0][:3], max(covered.values())
for D, fn in ((5, 'coverage11.json'), (6, 'coverage14_D6.json'), (7, 'coverage14_D7.json')):
    led = json.load(open(C + fn))
    print(f'D={D} ledger {fn}: counts={led.get("counts")}; fresh tree check (cases, covered, uncovered sample, max multiplicity):', tree_check(led, D, ['515', '617']))
    certs = Counter(); 
    for idx, node in led['nodes'].items():
        if node['status'] in ('LP-certified', 'joint-LP-certified'): certs[node.get('certificate')] += 1
    print('   certificates referenced by certified nodes:', dict(certs))
# reduction certificates
print('\n== reduction certificates (session V/VI) ==')
summary = {}
for fn, expect in [('lp_probe_rmax4.json', 'rmax<=4'), ('lp_probe_cfg_5_15.json', 'config (5,15) present'), ('lp_probe_cfg_5_16.json', 'config (5,16) present'), ('lp_probe_cfg_5_17.json', None), ('lp_probe_cfg_5_18.json', None), ('lp_probe_cfg_5_21.json', None), ('lp_probe_cfg_6_17.json', None), ('lp_probe_cfg_6_18.json', None), ('lp_probe_cfg_6_19.json', None), ('lp_probe_cfg_6_21.json', None), ('lp_probe_cfg_6_25.json', None), ('lp_probe_f7_cfg_5_15.json', None), ('lp_probe_f7_cfg_5_16.json', None), ('lp_probe_f7_cfg_6_17.json', None)]:
    d = json.load(open(C + fn)); row = {}
    for key, rec in d.items():
        side = [r for r in rec.get('cert', {}) if r.startswith(('config', 'rmax', 'forbid', '#planes', 'mfix'))]
        row[key] = (rec['status'][:10], side)
    summary[fn] = row
    print(fn, {k: v[0] for k, v in row.items()})
    sides = {k: v[1] for k, v in row.items() if v[1]}
    print('    side rows with nonzero multiplier:', sides)
json.dump({'raw_counts': counts, 'reduction': summary}, open('out/coverage_result.json', 'w'), indent=1)
