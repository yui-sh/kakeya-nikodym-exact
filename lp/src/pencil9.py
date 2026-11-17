"""Fixed-leaf histogram pencil column generation; see WRITEUP9.md.
New columns link to existing F8 G and Z variables. High-r choices are explicit;
low-r items keep their old item labels, with an integer histogram completion
carried as a column witness. All pricing arithmetic is integral.
"""
import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy import sparse as sp
from scipy.optimize import linprog
import planecount6 as p6


def column_counts(gids, zid, grow, zrow):
    ans = Counter(grow[g] for g in gids)
    ans[zrow[zid]] += 1
    return dict(ans)


def compatible(m, dv, labels, q):
    return (len(labels) == q + 1 and all(all(d <= h <= cap for d, h, cap in zip(dv, key[2], m)) for key in labels)
            and all(sum(key[2][j] for key in labels) == m[j] + q * dv[j] for j in range(len(m))))


def build(config, m):
    q, n, D = 5, 46, len(m)
    original = p6.hist_types
    # Necessary fixed-leaf bounds: each plane's points are a subset of K.
    def filtered(*args, **kwargs):
        return [t for t in original(*args, **kwargs) if t[0] <= config[0] and all(h <= cap for h, cap in zip(t[2], m))]
    p6.hist_types = filtered
    try: lp, info, types = p6.build(q, n, D, mfix=m)
    finally: p6.hist_types = original
    for j, val in enumerate(m, 1): lp.add([(lp.col[f'm{j}'], 1)], '=', val, f'leaf m{j}')
    lp.add([(i, -1) for i, t in enumerate(types) if t[0] == config[0] and t[3] == config[1]], '<', -1, f'config {config} present')
    grow = {}; zrow = {}
    for g in range(len(info['_Y'])):
        grow[g] = len(lp.rows); lp.add([(lp.col[f'G{g}'], -1)], '=', 0, f'P9 G{g}')
    for z in range(len(info['_Z'])):
        zrow[z] = len(lp.rows); lp.add([(lp.col[f'Z{z}'], -1)], '=', 0, f'P9 Z{z}')
    return lp, info, types, grow, zrow


class Oracle:
    def __init__(self, info, m, grow, zrow, tag):
        self.info, self.m, self.grow, self.zrow, self.tag = info, m, grow, zrow, tag
        items = info['_items']; D = len(m)
        groups = defaultdict(lambda: defaultdict(list))
        for g, (key, di) in enumerate(info['_Y']):
            dv = tuple(info['_dms'][di].count(j) for j in range(1, D + 1))
            if any(d > cap for d, cap in zip(dv, m)): continue
            ii = next(i for i, item in enumerate(items) if key[0] == item[0] and key[1] in item[3])
            groups[di][ii].append((g, tuple(h-d for h,d in zip(key[2], dv))))
        patterns = defaultdict(list)
        for z, (di, wi) in enumerate(info['_Z']):
            dv = tuple(info['_dms'][di].count(j) for j in range(1, D + 1))
            if any(d > cap for d, cap in zip(dv, m)): continue
            w = info['_wt'][wi]
            if any(ii not in groups[di] for ii in w): continue
            # Enumerate high-r histograms first; share low-r completion suffix DP.
            seq = sorted(w, key=lambda ii: (items[ii][0] < 5, ii))
            patterns[di].append((z, seq, tuple(cap-d for cap,d in zip(m,dv))))
        self.batches = [(di, dict(groups[di]), patterns[di]) for di in sorted(patterns)]

    def price(self, y):
        inp = Path(f'out/{self.tag}_price.in'); output = Path(f'out/{self.tag}_price.out')
        with inp.open('w') as f:
            print(len(self.m), len(self.batches), file=f)
            for di, groups, patterns in self.batches:
                ids = {ii:i for i,ii in enumerate(sorted(groups))}
                print(len(groups), len(patterns), file=f)
                for ii in sorted(groups):
                    print(len(groups[ii]), file=f)
                    for g, residual in groups[ii]: print(g, int(y[self.grow[g]]), *residual, file=f)
                for z, seq, target in patterns: print(z, int(y[self.zrow[z]]), len(seq), *(ids[ii] for ii in seq), *target, file=f)
        with open(f'out/{self.tag}_price.log', 'a') as log:
            subprocess.run(['out/pencil9_price', str(inp), str(output)], stderr=log, check=True)
        result = []
        for line in output.read_text().splitlines():
            fields = line.split(); z = int(fields[0])
            if fields[1] == 'NONE': continue
            cost = int(fields[1]); gids = tuple(sorted(map(int, fields[2:])))
            di, wi = self.info['_Z'][z]; dm = self.info['_dms'][di]
            dv = tuple(dm.count(j) for j in range(1, len(self.m)+1))
            labels = [self.info['_Y'][g][0] for g in gids]
            assert all(self.info['_Y'][g][1] == di for g in gids)
            assert compatible(self.m, dv, labels, 5)
            actual_items = sorted(next(ii for ii, it in enumerate(self.info['_items']) if k[0] == it[0] and k[1] in it[3]) for k in labels)
            assert actual_items == list(self.info['_wt'][wi])
            assert cost == int(y[self.zrow[z]]) + sum(int(y[self.grow[g]]) for g in gids)
            result.append((cost, z, gids))
        return result


def solve(config, m, tag, rounds):
    t0 = time.monotonic()
    lp, info, types, grow, zrow = build(config, m)
    A, b, sense = lp.matrices(); base_n = A.shape[1]
    eq = np.array([i for i,s in enumerate(sense) if s == '=']); ub = np.array([i for i,s in enumerate(sense) if s == '<'])
    print('BUILD', A.shape, 'Y', len(grow), 'Z', len(zrow), 'seconds', time.monotonic()-t0, flush=True)
    oracle = Oracle(info, m, grow, zrow, tag)
    columns = []; seen = set(); rounds_log = []
    def add(priced, seed=False):
        rr=[]; cc=[]; vv=[]; count=0
        for cost,z,gids in priced:
            if not seed and cost >= 0: continue
            key = (z,gids)
            if key in seen: continue
            seen.add(key); columns.append(key)
            for row,value in column_counts(gids,z,grow,zrow).items(): rr.append(row); cc.append(count); vv.append(value)
            count += 1
        return sp.csc_matrix((vv,(rr,cc)), shape=(A.shape[0],count))
    priced = oracle.price([0]*A.shape[0])
    C = add(priced, True)
    print('SEED', C.shape[1], 'completed patterns out of', len(zrow), 'seconds', time.monotonic()-t0, flush=True)
    status = 'iteration_limit'; x = None
    for iteration in range(rounds):
        ts=time.monotonic(); master=sp.hstack([A,C],format='csr'); nv=master.shape[1]
        Ae=sp.hstack([master[eq],sp.eye(len(eq)), -sp.eye(len(eq)),sp.csr_matrix((len(eq),len(ub)))],format='csr')
        Au=sp.hstack([master[ub],sp.csr_matrix((len(ub),2*len(eq))),-sp.eye(len(ub))],format='csr')
        res=linprog(np.r_[np.zeros(nv),np.ones(2*len(eq)+len(ub))], A_eq=Ae,b_eq=b[eq],A_ub=Au,b_ub=b[ub],bounds=(0,None),method='highs',options={'threads':1})
        rec={'iteration':iteration,'columns':C.shape[1],'solver_status':res.status,'solve_s':time.monotonic()-ts}
        if res.status != 0: status='solver_failure'; rounds_log.append(rec); break
        rec['phase1']=float(res.fun)
        if res.fun < 1e-7:
            x=res.x[:nv]; residual=master@x-b
            rec.update(max_eq_residual=float(np.max(np.abs(residual[eq]))),max_ineq_violation=float(max(0,np.max(residual[ub]))),min_variable=float(np.min(x)))
            status='computed_feasible'; rounds_log.append(rec); print(rec,flush=True); break
        y=np.zeros(len(b));y[eq]=-res.eqlin.marginals;y[ub]=-res.ineqlin.marginals
        # Integer objective is the exact rounded dual, not a floating lower bound.
        yi=[int(round(v*10**9)) for v in y]
        ts=time.monotonic(); priced=oracle.price(yi); new=add(priced)
        rec.update(pricing_s=time.monotonic()-ts,new_columns=new.shape[1],min_integer_price=min((p[0] for p in priced),default=None))
        rounds_log.append(rec); print(rec,flush=True)
        np.savez_compressed(f'out/{tag}_dual.npz', y=y)
        if not new.shape[1]:
            # No exclusion is claimed without a separate exact all-column check.
            status='uncertified_stall'; break
        C=sp.hstack([C,new],format='csc')
    result=dict(config=config,m=m,status=status,elapsed_s=time.monotonic()-t0,base_variables=base_n,rows=len(b),columns=C.shape[1],rounds=rounds_log)
    json.dump(result,open(f'out/{tag}.json','w'),indent=2)
    json.dump(columns,open(f'out/{tag}_columns.json','w'))
    if x is not None: np.savez_compressed(f'out/{tag}_primal.npz',x=x)
    print('DONE',json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--m',required=True);ap.add_argument('--tag',required=True);ap.add_argument('--rounds',type=int,default=200)
    a=ap.parse_args()
    assert p6.SECTIONS and p6.DMH and p6.UNMARKED and p6.EXACT_K==4 and p6.F6_RMIN==5
    solve(tuple(map(int,a.config.split(','))),tuple(map(int,a.m.split(','))),a.tag,a.rounds)
