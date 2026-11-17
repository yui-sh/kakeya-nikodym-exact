"""Exact geometric validation of F9 rows, without rebuilding the large F8 model.
Checks every unmarked direction of EVERY plane, including r < F6_RMIN.
The old rows were validated separately in out/f8_validate.log.
"""
import json
import sys
from collections import defaultdict
from geom import AG
from structure import planes_of
from planecount_check import witness
from pc6_validate import dot
from planecount6 import unmarked_deficit


def validate(fn, degree_cap=None):
    data = json.load(open(fn))
    q = data['q']; ag = AG(3, q)
    S = {ag.pidx[tuple(p)] for p in data['points']}
    W = witness(ag, S)
    degree = {p: 0 for p in S}
    for l in W:
        for p in ag.lines[l][1]: degree[p] += 1
    D = max(max(degree.values()), degree_cap or 0)
    planes = planes_of(ag)
    stats = [(sum(ag.lines[l][1] <= pts for l in W), len(S & pts)) for _, _, pts in planes]
    wdir = {ag.lines[l][0]: l for l in W}
    rows = defaultdict(int); count = selected = bad = 0
    for i, (normal, _, pts) in enumerate(planes):
        r, k = stats[i]
        for e, l in wdir.items():
            lpts = ag.lines[l][1]
            if dot(ag.F, normal, ag.dirs[e]) or lpts <= pts: continue
            mates = [j for j, (nm, _, ps) in enumerate(planes) if nm == normal and lpts <= ps]
            assert len(mates) == 1 and mates[0] != i
            mate = stats[mates[0]]
            pencil = [j for j, (_, _, ps) in enumerate(planes) if lpts <= ps]
            assert len(pencil) == q + 1
            s = sum(stats[j][0] for j in pencil) - 1
            assert s == sum(degree[p] for p in lpts)
            assert sum(degree[p] for p in planes[mates[0]][2] & S) == q*q + q*mate[0]
            profile = tuple(sorted(len(ag.lines[a][1] & S) for a in ag.lines_by_dir[e] if ag.lines[a][1] <= pts))
            # Independently evaluate the inequality, then compare the implemented coefficient.
            deficit = max(0, q*q + q*mate[0] - D*(mate[1]-q) - sum(degree[p] for p in lpts))
            assert deficit == unmarked_deficit(q, D, mate, s)
            rows[(profile, mate)] += deficit
            bad += deficit != 0; count += 1; selected += r >= 5
    assert all(lhs == 0 for lhs in rows.values()) and not bad
    print(f'{fn}: n={len(S)} D={D} flags={count} r>=5_flags={selected} aggregated_rows={len(rows)} violated={bad}', flush=True)


if __name__ == '__main__':
    for fn in sys.argv[1:]: validate(fn)
    print('ALL OK', flush=True)
