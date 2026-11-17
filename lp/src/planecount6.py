"""Session V: merged aggregated LP (WRITEUP5 §2-3). Plane types carry the degree histogram h = (h_1..h_D) of K∩π and
must admit an inside configuration (inside_profiles.py) whose crossing numbers are dominated by h (a point on i inside
lines has degree >= i). Degree-pair counts M_ij linearise the second-order degree identities without fixing (m_j).
Families: plane types N[r,k,h]; degrees m_j and pairs M_ij; line sections L_j; parallel-class types C; direction
types A + per-direction line vectors P (planecount5); bucketed witness-line types W.
Usage: planecount6.py --q 5 --n 46 --D 4 [--no-classes --no-dir --no-lines --no-hist] [--cert] [--show]"""
import argparse, math, json, sys, time, itertools
from fractions import Fraction
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog
from plane_f import get
from union_tables import tables
from inside_profiles import profiles
import planecount5 as p5
import os
SECTIONS = os.environ.get("PC6_SECTIONS", "0") == "1"   # session VI family F5: exact per-plane line-section profiles (sections.py)
F6_RMIN = int(os.environ.get("PC6_F6_RMIN", "2"))        # flags (pi, e) for planes with r >= F6_RMIN and every unmarked direction e
DMH = os.environ.get("PC6_DMH", "0") == "1"             # session VII family F8: witness-line degree multiset vs. plane histograms
UNMARKED = os.environ.get("PC6_UNMARKED", "0") == "1"   # session VIII: exact mate version of Lemma C
EXACT_K = int(os.environ.get("PC6_EXACT_K", "0"))       # pencil items: k exact for k <= f(r)+EXACT_K-1 (0 = session V buckets {lo},{lo+1,lo+2},[lo+3,..])
from sections import sections as _sections

def unmarked_deficit(q, D, mate, pencil_s):
    """Positive part of the exact mate-plane Lemma C deficit (WRITEUP8 §1)."""
    r, k = mate
    return max(0, q * q + q * r - D * (k - q) - pencil_s)

def hist_types(q, D, n, use_hist=True, u_for_r=None, use_sections=None):
    """plane types (r, k, h, u) or, with sections (F5), (r, k, h, u, t): t = exact line-section profile of the plane
    (t_j = # lines of pi with exactly j points of K), u the union size of the r witness lines, h dominated by the exact
    crossing profile c of that (u, c, t)-realisation. h = None when use_hist is False (then only (r,k))."""
    if use_sections is None: use_sections = SECTIONS
    nd = q * q + q + 1; base = nd - q - 1; f = get(q); prof = profiles(q)
    types = []
    if use_sections:
        def hists(k, S):
            def rec(j, left_k, left_S, cur):
                if j == D:
                    if left_S == D * left_k: yield tuple(cur + [left_k])
                    return
                for h in range(left_k + 1):
                    rk, rS = left_k - h, left_S - j * h
                    if rS < (j + 1) * rk or rS > D * rk: continue
                    yield from rec(j + 1, rk, rS, cur + [h])
            return list(rec(1, k, S, []))
        hcache = {}
        seen = set()
        for (k, r, u, c, t, e) in _sections(q, F6_RMIN):
            if k > min(q * q, n) or r > q + 1 or k < max(f[r], 1): continue
            S = base + q * r
            if S < k or S > D * k: continue
            if not use_hist:
                key = (r, k, None, u, t, e)
                if key not in seen: seen.add(key); types.append(key)
                continue
            if (k, S) not in hcache: hcache[(k, S)] = hists(k, S)
            for h in hcache[(k, S)]:
                tail_h = [sum(h[j - 1:]) for j in range(1, D + 2)]
                if all(sum(c[i - 1:]) <= (tail_h[i - 1] if i <= D else 0) for i in range(1, r + 1)):
                    key = (r, k, h, u, t, e)
                    if key not in seen: seen.add(key); types.append(key)
        return types
    for r in range(q + 2):
        S = base + q * r
        for k in range(max(f[r], 1), min(q * q, n) + 1):
            if S < k or S > D * k: continue
            if not use_hist: types.append((r, k, None, None)); continue
            profs = [p for p in prof[r] if p[0] <= k]
            if not profs: continue
            # multisets of k values in [1,D] with sum S, as histogram tuples
            def rec(j, left_k, left_S, cur):
                if j == D:
                    if left_S == D * left_k: yield tuple(cur + [left_k])
                    return
                for h in range(left_k + 1):
                    rk, rS = left_k - h, left_S - j * h
                    if rS < (j + 1) * rk or rS > D * rk: continue
                    yield from rec(j + 1, rk, rS, cur + [h])
            for h in rec(1, k, S, []):
                tail_h = [sum(h[j - 1:]) for j in range(1, D + 2)]  # tail_h[i-1] = #{d_p >= i}, i = 1..D+1
                us = set()
                for p in profs:
                    c = p[1:]  # c_1..c_r
                    if all(sum(c[i - 1:]) <= (tail_h[i - 1] if i <= D else 0) for i in range(1, r + 1)): us.add(p[0])
                if u_for_r is not None and r != u_for_r: us = {None} if us else set()
                for u in sorted(us, key=lambda v: -1 if v is None else v): types.append((r, k, h, u))   # type = (r, k, degree histogram, union size u)
    return types

def build(q, n, D, use_classes=True, use_lines=True, use_dir=True, use_hist=True, max_types=400000, u_for_r=None, use_sections=None, mfix=None):
    if use_sections is None: use_sections = SECTIONS
    nd = q * q + q + 1; npl = q * nd; f = get(q); bm = f[q + 1]; base = nd - q - 1
    types = hist_types(q, D, n, use_hist, u_for_r, use_sections); info = {"plane_types": len(types), "sections": use_sections}
    if len(types) > max_types: raise RuntimeError(f"too many plane types: {len(types)}")
    lp = p5.LP()
    N = [lp.var(f"N{i}") for i in range(len(types))]
    m = {j: lp.var(f"m{j}") for j in range(1, D + 1)}
    L = {j: lp.var(f"L{j}") for j in range(1, q + 1)}
    R = np.array([t[0] for t in types]); K = np.array([t[1] for t in types])
    def tsum(coef): return [(N[i], int(coef[i])) for i in range(len(types)) if coef[i]]
    lp.add(tsum(np.ones(len(types), int)), '=', npl, "I2a: sum_pi 1")
    lp.add(tsum(R), '=', (q + 1) * nd, "I2b: sum_pi r")
    lp.add(tsum(K), '=', nd * n, "I2c: sum_pi k")
    lp.add([(m[j], 1) for j in m], '=', n, "deg: sum m_j")
    lp.add([(m[j], j) for j in m], '=', nd * q, "deg: sum j m_j")
    lp.add([(m[D], -1)], '<', -1, "m_D >= 1")
    lp.add(tsum(R * (R - 1) // 2) + [(m[j], -math.comb(j, 2)) for j in m], '=', 0, "I3: sum C(r,2) = sum C(d,2)")
    lp.add(tsum(K * K), '=', (q + 1) * n * n + base * n, "I4: sum k^2")
    lp.add(tsum(K * R), '=', nd * n + q * q * nd, "I5: sum k r")
    lp.add([(L[j], j) for j in L], '=', nd * n, "I6a: sum j L_j")
    lp.add([(L[j], math.comb(j, 2)) for j in L], '=', math.comb(n, 2), "I6b: sum C(j,2) L_j")
    lp.add([(L[j], -1) for j in L], '<', -nd * bm, "I6c: sum L_j >= nd*BM")
    lp.add([(L[q], -1)], '<', -nd, "I6d: L_q >= nd")
    lp.add(tsum(K * (K - 1) * (K - 2) // 6) + [(L[j], -q * math.comb(j, 3)) for j in L], '=', math.comb(n, 3), "I6e: sum C(k,3) - q sum C(j,3) L_j")
    if use_sections:   # F5: every line of AG(3,q) lies in q+1 planes, so sum_pi t_j(pi) = (q+1) L_j
        T = np.array([t[4] for t in types])
        for j in range(1, q + 1):
            lp.add(tsum(T[:, j - 1]) + [(L[j], -(q + 1))], '=', 0, f"F5: sum_pi t{j} = (q+1) L{j}")
    if use_hist:
        H = np.array([t[2] for t in types])  # T x D
        M = {}
        for i in range(1, D + 1):
            for j in range(i, D + 1): M[(i, j)] = lp.var(f"M{i},{j}")
        for j in range(1, D + 1):
            hj = H[:, j - 1]
            lp.add(tsum(hj) + [(m[j], -nd)], '=', 0, f"H1: sum_pi h{j} = nd m{j}")
            lp.add(tsum(hj * K) + [(m[j], -((q + 1) * n + base))], '=', 0, f"H2: sum_pi h{j} k = ((q+1)n+q^2-q) m{j}")
            lp.add(tsum(hj * R) + [(m[j], -(nd + q * j))], '=', 0, f"H3: sum_pi h{j} r = (nd+q j) m{j}")
            lp.add(tsum(hj * (hj - 1) // 2) + [(M[(j, j)], -(q + 1))], '=', 0, f"H4: sum_pi C(h{j},2) = (q+1) M{j}{j}")
            for i in range(1, j):
                lp.add(tsum(H[:, i - 1] * hj) + [(M[(i, j)], -(q + 1))], '=', 0, f"H5: sum_pi h{i} h{j} = (q+1) M{i}{j}")
            lp.add([(M[(min(i, j), max(i, j))], 1) for i in range(1, D + 1) if i != j] + [(M[(j, j)], 2), (m[j], -(n - 1))], '=', 0, f"H6: sum_i M_i{j} + 2 M{j}{j} = (n-1) m{j}")
            lp.add([(M[(j, j)], 2), (m[j], -(n - 1))], '<', 0, f"H7: 2 M{j}{j} <= (n-1) m{j}")
            for a in range(1, n):
                lp.add([(M[(j, j)], -1), (m[j], a)], '<', math.comb(a + 1, 2), f"H8: M{j}{j} >= {a} m{j} - C({a+1},2)")
    rk_index = {}
    for i, t in enumerate(types): rk_index.setdefault((t[0], t[1]), []).append(i)
    pairs = sorted(rk_index)
    ctypes = []; C = []
    if use_classes:
        ctypes = p5.multisets(pairs, q, n, sum_r_exact=q + 1); info["class_types"] = len(ctypes)
        C = [lp.var(f"C{i}") for i in range(len(ctypes))]
        lp.add([(c, 1) for c in C], '=', nd, "classes: sum C = nd")
        for rk in pairs:
            lp.add([(C[i], t.count(rk)) for i, t in enumerate(ctypes) if rk in t] + [(N[i], -1) for i in rk_index[rk]], '=', 0, f"class link N{rk}")
    if use_dir:
        assert use_classes
        F, Ft = tables(q); utypes = sorted(Ft); info["dir_types"] = len(utypes)
        A = [lp.var(f"A{i}") for i in range(len(utypes))]
        lp.add([(a, 1) for a in A], '=', nd, "dirs: sum A = nd")
        def phi(t, v):
            rho = sum(1 for r, k in t if r >= 1); n1 = sum(1 for r, k in t if r == 1)
            return (q + 1 - n1) * (rho == v) + n1 * (rho - 1 == v)
        for v in range(1, q + 1):
            lp.add([(A[i], u.count(v)) for i, u in enumerate(utypes)] + [(C[i], -phi(t, v)) for i, t in enumerate(ctypes) if phi(t, v)], '=', 0, f"dir flags a=={v}")
        wvecs = p5.line_vectors(q, n); info["line_vectors"] = len(wvecs)
        P = [lp.var(f"P{i}") for i in range(len(wvecs))]
        lp.add([(p, 1) for p in P], '=', nd, "dirs: sum P = nd")
        for j in range(1, q + 1):
            lp.add([(P[i], w[j - 1]) for i, w in enumerate(wvecs) if w[j - 1]] + [(L[j], -1)], '=', 0, f"dir link L{j}")
        lam = [sum(w) for w in wvecs]
        for theta in range(bm + 1, q * q + 1):
            lp.add([(P[i], -1) for i in range(len(wvecs)) if lam[i] >= theta] + [(A[i], 1) for i, u in enumerate(utypes) if Ft[u] >= theta], '<', 0, f"dominance lambda>={theta}")
    if use_lines:
        # buckets of k per r: {f(r)}, {f(r)+1, f(r)+2}, [f(r)+3, kmax]
        items = []
        for r in range(1, q + 2):
            ks = sorted(k for (rr, k) in pairs if rr == r)
            if not ks: continue
            lo = ks[0]
            buckets = [(lo, lo), (lo + 1, lo + 2), (lo + 3, ks[-1])] if EXACT_K == 0 else [(lo + i, lo + i) for i in range(EXACT_K)] + [(lo + EXACT_K, ks[-1])]
            for (a, b) in buckets:
                kk = [k for k in ks if a <= k <= b]
                if kk: items.append((r, min(kk), max(kk), tuple(kk)))
        wt = []
        def rec(start, left, sr, slo, shi, cur):
            if left == 0:
                if slo <= n + q * q <= shi and q + 1 <= sr <= 1 + q * D:
                    # F7 (session VI): s = sum_{p in l} d_p = sr - 1; for every plane pi of the pencil, (I1) gives
                    # sum_{p in pi\l} d_p = base + q r_pi - s <= D (k_pi - q), so s >= base + q r_pi - D (k_pi - q).
                    s = sr - 1
                    if all(s >= base + q * items[i][0] - D * (items[i][2] - q) for i in cur): wt.append(tuple(cur))
                return
            for i in range(start, len(items)):
                r, lo, hi, kk = items[i]
                if slo + lo > n + q * q or sr + r > 1 + q * D: continue
                cur.append(i); rec(i, left - 1, sr + r, slo + lo, shi + hi, cur); cur.pop()
        rec(0, q + 1, 0, 0, 0, []); info["wline_types"] = len(wt)
        W = [lp.var(f"W{i}") for i in range(len(wt))]
        lp.add([(w, 1) for w in W], '=', nd, "wlines: sum W = nd")
        for ii, (r, lo, hi, kk) in enumerate(items):
            lp.add([(W[i], t.count(ii)) for i, t in enumerate(wt) if ii in t] + [(N[j], -r) for k in kk for j in rk_index.get((r, k), [])], '=', 0, f"wline link r*N[{r},{lo}..{hi}]")
        lp.add([(W[i], sum(items[ii][0] for ii in t) - 1) for i, t in enumerate(wt)] + [(m[j], -j * j) for j in m], '=', 0, "wlines: sum s_l = sum d_p^2")
        # degree multisets on witness lines: V[dm], dm = q-multiset of degrees in [1,D]; s = sum dm = (sum r over pencil) - 1
        dms = list(itertools.combinations_with_replacement(range(1, D + 1), q)); info["wline_degree_types"] = len(dms)
        V = [lp.var(f"V{i}") for i in range(len(dms))]
        lp.add([(v, 1) for v in V], '=', nd, "wdeg: sum V = nd")
        for j in range(1, D + 1):
            lp.add([(V[i], dm.count(j)) for i, dm in enumerate(dms) if j in dm] + [(m[j], -j)], '=', 0, f"wdeg: sum_l #(deg {j} on l) = {j} m{j}")
        for s in range(q, q * D + 1):
            lp.add([(V[i], 1) for i, dm in enumerate(dms) if sum(dm) == s] + [(W[i], -1) for i, t in enumerate(wt) if sum(items[ii][0] for ii in t) - 1 == s], '=', 0, f"wdeg/wline link s={s}")
        if use_hist:
            for i in range(1, D + 1):
                for j in range(i, D + 1):
                    coef = []
                    for ii, dm in enumerate(dms):
                        ci, cj = dm.count(i), dm.count(j)
                        c = math.comb(ci, 2) if i == j else ci * cj
                        if c: coef.append((V[ii], c))
                    if coef: lp.add(coef + [(M[(i, j)], -1)], '<', 0, f"wdeg: witness pairs ({i},{j}) <= M{i}{j}")
    if use_sections and use_lines and use_classes:
        # F6 (session VI, general form): flags (pi, e) for every plane pi with r_pi >= F6_RMIN and every unmarked
        # direction e of its line at infinity (profile s = sections of the q e-lines of pi, from the enumeration).
        # The plane pi' of pi's class containing l_e has type (r',k'), r' >= 1; the other q planes through l_e cut pi in
        # its q e-lines, so each must accommodate one s_j: k_sigma >= f(r_sigma) + max(0, s_j - r_sigma + 1) (Lemma B).
        # Each plane of type (r',k') in a class hosts exactly r' of the unmarked directions of the other planes.
        rmin = F6_RMIN
        def ok_item(ii, sv): r, lo, hi, kk = items[ii]; return hi >= f[r] + max(0, sv - r + 1)
        def matchable(rest, e):
            masks = [sum(1 << j for j in range(q) if ok_item(ii, e[j])) for ii in rest]
            reach = {0}
            for mk in masks:
                reach = {rm | (1 << j) for rm in reach for j in range(q) if (mk >> j) & 1 and not (rm >> j) & 1}
                if not reach: return False
            return True
        # class coefficients g[c][(rk, rk')] = # (plane of type rk, plane of type rk' != it) pairs in c times r'
        prof_of = {}   # (r,k) -> set of profiles s occurring
        for t in types:
            if t[0] >= rmin and t[5]:
                for sp in t[5]: prof_of.setdefault((t[0], t[1]), set()).add(sp)
        kinds = set(); gco = {}
        for ci, ct in enumerate(ctypes):
            cnt = {}
            for rk in ct: cnt[rk] = cnt.get(rk, 0) + 1
            for rk, m1 in cnt.items():
                if rk[0] < rmin or rk not in prof_of: continue
                for rk2, m2 in cnt.items():
                    if rk2[0] < 1: continue
                    g = m1 * (m2 - (1 if rk2 == rk else 0)) * rk2[0]
                    if g: gco.setdefault((rk, rk2), []).append((C[ci], -g)); kinds.update((rk[0], rk[1], sp, rk2) for sp in prof_of[rk])
        kinds = sorted(kinds)
        X1 = {}
        for kd in kinds: X1[kd] = lp.var(f"X{len(X1)}")
        srk2 = sorted({(kd[2], kd[3]) for kd in kinds}); srk2set = set(srk2)
        X2 = {}; cache = {}
        for i, w in enumerate(wt):
            for pos in {w[p]: p for p in range(len(w))}.values():
                r2, lo, hi, kk = items[w[pos]]
                if r2 < 1: continue
                rest = tuple(sorted(w[:pos] + w[pos + 1:]))
                for sp in {kd[2] for kd in kinds}:
                    key = (rest, sp)
                    if key not in cache: cache[key] = matchable(rest, sp)
                    if cache[key]:
                        for k2 in kk:
                            if (sp, (r2, k2)) in srk2set and ((sp, (r2, k2)), i) not in X2: X2[((sp, (r2, k2)), i)] = lp.var(f"Y{len(X2)}")
        info["flag_kinds"] = len(kinds); info["flag_vars"] = len(X1) + len(X2); info["_fkinds"] = kinds; info["_X1"] = list(X1); info["_X2"] = list(X2); info["_f6_rmin"] = rmin
        by_rks = {}; by_rkrk2 = {}; by_srk2 = {}; by_w = {}
        for kd, v in X1.items():
            by_rks.setdefault((kd[0], kd[1], kd[2]), []).append(v); by_rkrk2.setdefault(((kd[0], kd[1]), kd[3]), []).append(v); by_srk2.setdefault((kd[2], kd[3]), []).append((v, 1))
        for ((sp, rk2), i), v in X2.items(): by_srk2.setdefault((sp, rk2), []).append((v, -1)); by_w.setdefault(i, []).append(v)
        types_rks = {}
        for j, t in enumerate(types):
            if t[0] >= rmin and t[5]:
                for sp in set(t[5]): types_rks.setdefault((t[0], t[1], sp), []).append((N[j], -t[5].count(sp)))
        for key, coef in types_rks.items():
            lp.add([(v, 1) for v in by_rks.get(key, [])] + coef, '=', 0, f"F6 R1 flags (r={key[0]},k={key[1]},s={key[2]})")
        for (rk, rk2), coef in gco.items():
            lp.add([(v, 1) for v in by_rkrk2.get((rk, rk2), [])] + coef, '=', 0, f"F6 R2 flags ({rk} hosted by {rk2})")
        for key, coef in sorted(by_srk2.items()):
            lp.add(coef, '=', 0, f"F6 link flags (s={key[0]},mate={key[1]})")
        if UNMARKED:
            cuts = {key: [] for key in srk2}
            for (key, i), v in X2.items():
                deficit = unmarked_deficit(q, D, key[1], sum(items[ii][0] for ii in wt[i]) - 1)
                if deficit: cuts[key].append((v, deficit))
            for key, coef in cuts.items():
                lp.add(coef, '<', 0, f"F9 unmarked C (s={key[0]},mate={key[1]})")
            info["unmarked_rows"] = len(cuts)
            info["unmarked_excluded_pairs"] = sum(len(coef) for coef in cuts.values())
        for i, vs in by_w.items():
            cap = sum((q + 1 - items[ii][0]) // rmin for ii in wt[i] if items[ii][0] >= 1)
            lp.add([(v, 1) for v in vs] + [(W[i], -cap)], '<', 0, f"F6 R3 flags on pencil {i} <= {cap} W{i}")
        # F6 R4 (capacity cut for (5,15)-flags, q = 5 only): x planes (5,15) with the same unmarked direction e give 5x
        # 3-point e-lines hosted by the six planes through l_e; two of them coincide only if they are the intersection
        # of two of the x planes (at most C(x,2) coincidences), and a pencil plane of type (r,k) hosts at most cap3(r,k)
        # distinct 3-point e-lines: (5,15): 2 (Lemma A), (1,k): (k-5)//3, (2,k): (k-9)//2, (3,k): min(4,k-12), else 4.
        if q == 5:
            def cap3(ii):
                r, lo, hi, kk = items[ii]
                if r == 5 and hi == 15: return 2
                if r == 1: return min(4, (hi - 5) // 3)
                if r == 2: return min(4, (hi - 9) // 2)
                if r == 3: return min(4, hi - 12)
                return 4
            full3 = (3,) * q
            # Preserve X2 insertion order while avoiding a full X2 scan per pencil.
            full3_by_w = {}
            for ((sp, rk2), j), v in X2.items():
                if sp == full3: full3_by_w.setdefault(j, []).append(v)
            for i, vs in by_w.items():
                vs3 = full3_by_w.get(i, [])
                if not vs3: continue
                tot = sum(cap3(ii) for ii in wt[i]); xmax = 0
                while xmax + 1 <= 6 and 5 * (xmax + 1) - math.comb(xmax + 1, 2) <= tot: xmax += 1
                lp.add([(v, 1) for v in vs3] + [(W[i], -xmax)], '<', 0, f"F6 R4 (3,3,3,3,3)-flags on pencil {i} <= {xmax} W{i}")
    if mfix is None and os.environ.get("PC6_MFIX"):   # degree vector fixed from the environment (lp_probe / milp_probe leaves)
        mfix = [int(v) for v in os.environ["PC6_MFIX"].split(",")]; assert len(mfix) == D
        for j in range(1, D + 1): lp.add([(m[j], 1)], '=', mfix[j - 1], f"mfix m{j}")
    if DMH and use_lines and use_hist:
        # F8 (session VII): degree multiset of a witness line vs. the histograms of the planes through it.
        # For a witness line l with degree multiset dm and any plane pi ⊃ l: the points of l are points of K∩pi, so
        # dm ⊆ h_pi (as multisets). Variables Y[(r,k,h), dm] = # incidences (plane of histogram type (r,k,h), witness
        # line inside it with multiset dm); Z[dm, w] = # witness lines with multiset dm and pencil type w.
        #   (Y1) sum_dm Y[(r,k,h),dm] = r * N[(r,k,h)]                 (each plane has r witness lines)
        #   (Y2) sum_{(r,k,h): (r,k) in item ii} Y[.,dm] = sum_w #ii(w) Z[dm,w]   (pencils carry the same planes)
        #   (Z1) sum_dm Z[dm,w] = W[w],  (Z2) sum_w Z[dm,w] = V[dm]          (replaces the s-link)
        # With the degree vector fixed (mfix), the planes through l partition K \ l, so sum_{pi ⊃ l} h_pi = m + q*dm:
        #   (Y3) sum_{(r,k,h)} h_j Y[(r,k,h),dm] = (m_j + q*#j(dm)) V[dm]  for every dm, j.
        rkh_idx = {}
        for i, t in enumerate(types):
            if t[0] >= 1: rkh_idx.setdefault((t[0], t[1], tuple(t[2])), []).append(i)
        dmv = [tuple(dm.count(j) for j in range(1, D + 1)) for dm in dms]
        Y = {}
        for key in rkh_idx:
            h = key[2]
            for di, dv in enumerate(dmv):
                if all(h[j] >= dv[j] for j in range(D)): Y[(key, di)] = lp.var(f"G{len(Y)}")
        info["dmh_Y"] = len(Y)
        y_by_key = {}; y_by_rkdm = {}
        for (key, di), v in Y.items():
            y_by_key.setdefault(key, []).append(v); y_by_rkdm.setdefault((key[0], key[1], di), []).append(v)
        for key, idxs in rkh_idx.items():
            lp.add([(v, 1) for v in y_by_key.get(key, [])] + [(N[i], -key[0]) for i in idxs], '=', 0, f"F8 Y1 {key}")
        # compatible (dm, w): every item of w has some (r,k in item, h) with a Y-pair for dm
        rk_dm_ok = {(r, k, di) for (r, k, di) in y_by_rkdm}
        Z = {}
        s_of = [sum(items[ii][0] for ii in t) - 1 for t in wt]
        for i, t in enumerate(wt):
            for di, dm in enumerate(dms):
                if sum(dm) != s_of[i]: continue
                if all(any((items[ii][0], k, di) in rk_dm_ok for k in items[ii][3]) for ii in set(t)): Z[(di, i)] = lp.var(f"Z{len(Z)}")
        info["dmh_Z"] = len(Z); info["_Y"] = list(Y); info["_Z"] = list(Z)
        z_by_w = {}; z_by_dm = {}; z_by_dm_item = {}
        for (di, i), v in Z.items():
            z_by_w.setdefault(i, []).append(v); z_by_dm.setdefault(di, []).append(v)
            for ii in set(wt[i]): z_by_dm_item.setdefault((di, ii), []).append((v, wt[i].count(ii)))
        for i in range(len(wt)):
            lp.add([(v, 1) for v in z_by_w.get(i, [])] + [(W[i], -1)], '=', 0, f"F8 Z1 w{i}")
        for di in range(len(dms)):
            lp.add([(v, 1) for v in z_by_dm.get(di, [])] + [(V[di], -1)], '=', 0, f"F8 Z2 dm{dms[di]}")
        for ii, (r, lo, hi, kk) in enumerate(items):
            for di in range(len(dms)):
                ys = [(v, -1) for k in kk for v in y_by_rkdm.get((r, k, di), [])]
                zs = z_by_dm_item.get((di, ii), [])
                if ys or zs: lp.add(zs + ys, '=', 0, f"F8 Y2 item{ii} dm{dms[di]}")
        if mfix is not None:
            y_by_dm = {}
            for (key, d2), v in Y.items(): y_by_dm.setdefault(d2, []).append((key, v))
            for di, dv in enumerate(dmv):
                for j in range(D):
                    coef = [(v, key[2][j]) for key, v in y_by_dm.get(di, []) if key[2][j]]
                    lp.add(coef + [(V[di], -(mfix[j] + q * dv[j]))], '=', 0, f"F8 Y3 dm{dms[di]} j{j+1}")
    info["_ctypes"] = ctypes
    if use_dir: info["_utypes"] = utypes; info["_wvecs"] = wvecs
    if use_lines: info["_items"] = items; info["_wt"] = wt; info["_dms"] = dms
    return lp, info, types

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--q", type=int, default=5); ap.add_argument("--n", type=int, nargs="+", default=[46])
    ap.add_argument("--D", type=int, nargs="+", default=None)
    for fl in ("classes", "lines", "dir", "hist"): ap.add_argument(f"--no-{fl}", action="store_true")
    ap.add_argument("--cert", action="store_true"); ap.add_argument("--show", action="store_true"); ap.add_argument("--tag", default="")
    ap.add_argument("--count-only", action="store_true")
    a = ap.parse_args(); q = a.q; nd = q * q + q + 1
    results = {}
    for n in a.n:
        Ds = a.D or list(range(math.ceil(nd * q / n), (n - 1) // (q - 1) + 1))
        row = []
        for D in Ds:
            t0 = time.time()
            if a.count_only:
                print(f"n={n} D={D}: plane types {len(hist_types(q, D, n, not a.no_hist))} ({time.time()-t0:.1f}s)", flush=True); continue
            lp, info, types = build(q, n, D, not a.no_classes, not a.no_lines, not a.no_dir, not a.no_hist)
            tb = time.time() - t0
            res, A, b, sense, eq, ub = p5.solve(lp)
            st = {0: "feasible", 2: "INFEASIBLE"}.get(res.status, res.message)
            rec = dict(status=st, info=info, nvars=len(lp.names), nrows=len(lp.rows), build_s=round(tb, 1), time=round(time.time() - t0, 1))
            if res.status == 2 and a.cert:
                yr, msg = p5.farkas_exact(lp, A, b, sense, eq, ub); rec["cert_msg"] = msg
                if yr is not None: rec["cert"] = {lp.rows[i][3]: str(yr[i]) for i in range(len(yr)) if yr[i] != 0}
            if res.status == 0 and a.show:
                sol = p5.show_solution(lp, res.x); rec["solution"] = sol
                print("   planes:", [(types[int(k[1:])][:2], types[int(k[1:])][2], v) for k, v in sol.items() if k[0] == 'N'])
                print("   m:", {k: v for k, v in sol.items() if k[0] == 'm'}, " L:", {k: v for k, v in sol.items() if k[0] == 'L'})
            results[f"{n},{D}"] = rec
            info = {k: v for k, v in info.items() if not k.startswith("_")}; rec["info"] = info
            print(f"n={n} D={D}: {st}  vars={rec['nvars']} rows={rec['nrows']} {info} build {rec['build_s']}s total {rec['time']}s" + (f"  cert: {rec.get('cert_msg','')}" if a.cert else ""), flush=True)
            row.append(f"D={D}:{'X' if st == 'INFEASIBLE' else 'ok'}")
        if not a.count_only: print(f"n={n}: " + "  ".join(row), flush=True)
    if not a.count_only:
        tag = a.tag or "full"
        json.dump(results, open(f"out/planecount6_q{q}_n{'_'.join(map(str,a.n))}_{tag}.json", "w"), indent=1)
