"""Session VIII: F8 + F9 probes and a size estimate for high-r histogram pencils.
Run with run_f9.sh; uses one HiGHS thread. Certificates use the existing phase-1 machinery.
"""
import copy
import functools
import json
import math
import time
from collections import Counter
import numpy as np
import planecount5 as p5
import planecount6 as p6

p5.linprog = functools.partial(p5.linprog, options={'threads': 1})
assert p6.SECTIONS and p6.DMH and p6.UNMARKED and p6.EXACT_K == 4 and p6.F6_RMIN == 5


def estimate(info):
    # Count high-r histogram multiset refinements of existing (dm,w) columns.
    # Identical item occurrences are indistinguishable. Low-r items stay bucketed.
    # This is before Y3 completion and exact-k-sum pruning, hence an upper bound.
    items, wt = info['_items'], info['_wt']
    hist = {}
    for (r, k, h), di in info['_Y']:
        for ii, (rr, lo, hi, kk) in enumerate(items):
            if r == rr and r >= 5 and k in kk:
                hist.setdefault((ii, di), set()).add((k, h))
    counts = Counter(); total = maximum = 0
    for di, wi in info['_Z']:
        factor = 1
        for ii, multiplicity in Counter(wt[wi]).items():
            if items[ii][0] >= 5:
                choices = len(hist.get((ii, di), ()))
                factor *= math.comb(choices + multiplicity - 1, multiplicity) if choices else 0
        total += factor; maximum = max(maximum, factor)
        counts[sum(items[ii][0] >= 5 for ii in wt[wi])] += factor
    return dict(base_Z=len(info['_Z']), refined_Z_upper_bound=total,
                max_refinements_per_Z=maximum, by_high_item_count=dict(counts),
                compatible_histograms_per_high_item=[dict(item=items[ii][:3], min=min(vals), max=max(vals))
                    for ii in range(len(items)) if items[ii][0] >= 5
                    for vals in [[len(v) for (jj, di), v in hist.items() if jj == ii]] if vals])


def main():
    D = 5
    while D <= 7:
        t0 = time.monotonic()
        lp0, info, types = p6.build(5, 46, D)
        public = {k: v for k, v in info.items() if not k.startswith('_')}
        print(f'D={D} built in {time.monotonic()-t0:.1f}s: {public}', flush=True)
        if D == 5:
            sizing = estimate(info)
            json.dump(sizing, open('out/f9_disaggregation_size.json', 'w'), indent=2)
            print('SIZE', sizing, flush=True)
        killed = []
        for r, u in [(5, 15), (6, 17)]:
            lp = copy.copy(lp0); lp.rows = list(lp0.rows)
            lp.add([(lp.col[f'N{i}'], -1) for i, t in enumerate(types) if t[0] == r and t[3] == u], '<', -1, f'config ({r},{u}) present')
            idx = [i for i, t in enumerate(types) if t[0] > r]
            if idx: lp.add([(lp.col[f'N{i}'], 1) for i in idx], '<', 0, f'rmax<={r}')
            A, b, sense = lp.matrices()
            eq = [i for i, s in enumerate(sense) if s == '=']; ub = [i for i, s in enumerate(sense) if s == '<']
            ts = time.monotonic()
            yr, msg, x = p5.farkas_phase1(lp, A, b, sense, eq, ub)
            result = dict(status='feasible' if x is not None else ('INFEASIBLE' if yr is not None else 'UNDECIDED'),
                          cert_msg=msg, solve_s=round(time.monotonic()-ts, 2), info=public, nvars=len(lp.names), nrows=len(lp.rows))
            if yr is not None:
                result['cert'] = {lp.rows[i][3]: str(y) for i, y in enumerate(yr) if y}
            if x is not None:
                residual = A @ x - b
                result['max_eq_residual'] = float(np.max(np.abs(residual[eq])))
                result['max_ineq_violation'] = float(max(0, np.max(residual[ub])))
                result['min_variable'] = float(np.min(x))
                np.savez_compressed(f'out/f9_primal_{r}_{u}_D{D}.npz', x=x)
            fn = f'out/lp_probe_f9_cfg_{r}_{u}.json'
            try: saved = json.load(open(fn))
            except FileNotFoundError: saved = {}
            saved[f'46,{D}'] = result
            json.dump(saved, open(fn, 'w'), indent=2)
            print(f'D={D} config={r},{u}: {result}', flush=True)
            if yr is not None:
                import subprocess
                subprocess.run(['.venv/bin/python', 'cert5_check.py', fn, '--config', f'{r},{u}'], check=True)
            killed.append(yr is not None)
        if not all(killed): break
        D += 1
    print('F9_DONE', flush=True)


if __name__ == '__main__': main()
