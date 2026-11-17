"""Conditional coplanar-triple upper bounds from the n=46 plane restrictions.

Local objectives are IDENTICAL to pointtriple14.py; no new assignment prices.
"""
import json, hashlib
from math import comb
from pathlib import Path
from collections import Counter
from lift11 import atomic, key

FAMILY = 'maximal-plane-star-four-triple-upper-v1'


def upper(j, B):
    """Upper bound when no projective line contains more than B selected points."""
    bounds = [0]
    if min(j, B) >= 3:
        bounds.append(comb(j, 2) // 3)
    for a in range(4, min(j, B) + 1):
        b = j - a
        bounds.append(comb(a, 3) + comb(b, 2) + comb(b, 3))
    return max(bounds)


def prepare():
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    from budget14 import polynomial
    from sections import sections
    records = []
    for fn in ('solutions/best_n3_q3_13.json', 'solutions/best_n3_q4_28_dvir.json',
               'solutions/best_n3_q5_53.json', 'solutions/bc_n3_q5_53.json'):
        data = json.load(open(fn)); ag = AG(3, data['q'])
        S = {ag.pidx[tuple(p)] for p in data['points']}; W = witness(ag, S)
        degree = Counter(p for l in W for p in ag.lines[l][1])
        maxima = Counter(); triples = Counter()
        for _, _, ps in planes_of(ag):
            lines = [l for l in W if ag.lines[l][1] <= ps]
            for p in S & ps:
                b = sum(p in ag.lines[l][1] for l in lines)
                maxima[p] = max(maxima[p], b)
                triples[p] += comb(b, 3)
        for p in S:
            assert triples[p] <= upper(degree[p], maxima[p])
        records.append(dict(file=fn, points=len(S), max_inside_multiplicity=max(maxima.values()),
                            checks=len(S), violations=0))
    sections_checked = []
    for cfg, rmax, u0 in [('515', 5, 15), ('617', 6, 17)]:
        admitted = [(k,r,u,c,t,e) for k,r,u,c,t,e in sections(5,5)
                    if r <= rmax and (r < rmax or u == u0) and polynomial(r,k) <= 541]
        B = max(max((b for b,n in enumerate(c,1) if n), default=0)
                for k,r,u,c,t,e in admitted)
        assert B == 4
        sections_checked.append(dict(config=cfg, profiles=len(admitted), max_multiplicity=B))
    assert [upper(j,4) for j in range(3,8)] == [1,4,4,5,8]
    atomic('out/star14_validation.json', dict(family=FAMILY, status='ZERO_VIOLATIONS',
           results=records, section_checks=sections_checked, upper_bounds=[1,4,4,5,8],
           argument='Five concurrent inside lines force union at least 21. At r=5, k>=21 has quadratic cost at least 750>541. At r=6 it contradicts maximal union 17. Thus inside multiplicity is at most four. For a maximal secant of size a>=4, triples are at most C(a,3)+C(j-a,2)+C(j-a,3); for a<=3 use 3T<=C(j,2).',
           validation_scope='Known sets test the general bound using their actual maximum secant; n=46 premises are separately proved, not assumed for n=53.',
           local_prices='Identical full-type objective and assignments to pointtriple14 positive-sign prices; no new local pricing problem.'))
    print('STAR14_VALIDATED', records, sections_checked, flush=True)


def guard(D, cfg):
    from maximal14 import guard as reduction_guard
    from budget14 import guard as budget_guard
    reduction_guard(D,cfg); budget_guard()
    data = json.load(open('out/star14_validation.json'))
    assert (data['family'], data['status']) == (FAMILY, 'ZERO_VIOLATIONS')
    assert data['upper_bounds'] == [upper(j,4) for j in range(3,8)]
    return hashlib.sha256(Path('out/star14_validation.json').read_bytes()).hexdigest()


def add_cuts(lp, types, D, cfg):
    from pointtriple14 import FAMILY as PRICE_FAMILY, objectives
    guard(D,cfg)
    cache = json.load(open(f'out/pointtriple14_D{D}_prices.json'))
    assert (cache['q'],cache['D'],cache['family']) == (5,D,PRICE_FAMILY)
    assert cache['section_sha256'] == hashlib.sha256(Path('out/sections3_q5_r5.json').read_bytes()).hexdigest()
    assert cache['objectives'] == json.loads(key(objectives(D)))
    for z,(j,sign) in enumerate(objectives(D)):
        if sign != 1 or j < 5: continue
        row = [(lp.col[f'N{i}'],cache['prices'][key(t)][z]) for i,t in enumerate(types)]
        row.append((lp.col[f'm{j}'],-upper(j,4)))
        lp.add(row,'<',0,f'XIV star-four triples degree {j}')
    print('STAR14_CUTS',D,D-4,flush=True)


if __name__ == '__main__':
    prepare()
    for cfg in ('515','617'): guard(7,cfg)
