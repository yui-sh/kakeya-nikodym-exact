"""Session V: validate every row of the merged LP (planecount6.build) on actual witness systems.
For each solution file: build the witness system, compute the exact value of every LP variable from the geometry
(plane types with degree histograms, degree pairs, line sections, parallel-class types, direction types a_{d,L}
computed directly from the geometry, per-direction line vectors, bucketed witness-line types, degree multisets on
witness lines) and check all equality / inequality rows exactly.  Usage: pc6_validate.py solutions/*.json"""
import sys, json, math, itertools
from fractions import Fraction
from geom import AG
from structure import planes_of
from planecount_check import witness
import planecount6 as p6

def dot(F, u, v):
    s = 0
    for a, b in zip(u, v): s = F.add[s][F.mul[a][b]]
    return s

def validate(fn, degree_cap=None):
    d = json.load(open(fn)); q = d["q"]; ag = AG(3, q); F = ag.F; nd = len(ag.dirs)
    S = {ag.pidx[tuple(p)] for p in d["points"]}; n = len(S); W = witness(ag, S); Wset = set(W)
    deg = {p: 0 for p in S}
    for l in W:
        for p in ag.lines[l][1]: deg[p] += 1
    actual_D = max(deg.values())
    D = max(actual_D, degree_cap or 0)
    mvec = [sum(1 for p in S if deg[p] == j) for j in range(1, D + 1)]
    lp, info, types = p6.build(q, n, D, mfix=mvec)
    val = {name: 0 for name in lp.names}
    tindex = {t: i for i, t in enumerate(types)}
    planes = planes_of(ag)  # (normal, c, pts)
    # plane statistics
    prk = {}; flag_planes = []
    for (nrm, c, pts) in planes:
        r = sum(1 for l in W if ag.lines[l][1] <= pts); k = len(pts & S)
        h = [0] * D
        for p in pts & S: h[deg[p] - 1] += 1
        u = len(set().union(*[ag.lines[l][1] for l in W if ag.lines[l][1] <= pts])) if r else 0
        t = (r, k, tuple(h), u)
        if p6.SECTIONS:   # session VI: exact section profile of the 30 lines of the plane
            tt = [0] * q
            for di, lpts in ag.lines:
                if lpts <= pts:
                    jj = len(lpts & S)
                    if jj: tt[jj - 1] += 1
            ee = None
            if r >= p6.F6_RMIN:   # profiles of all unmarked directions (directions at infinity whose witness line is not in the plane)
                marked = {ag.lines[l][0] for l in W if ag.lines[l][1] <= pts}
                unm = [di for di in range(nd) if dot(F, nrm, ag.dirs[di]) == 0 and di not in marked]; assert len(unm) == q + 1 - r
                profs = []
                for di in unm:
                    pr = tuple(sorted(len(ag.lines[l][1] & S) for l in ag.lines_by_dir[di] if ag.lines[l][1] <= pts))
                    profs.append(pr); flag_planes.append((nrm, c, pts, r, k, pr, di))
                ee = tuple(sorted(profs)) if profs else None
            t = t + (tuple(tt), ee)
        assert t in tindex, ("plane type not generated", t)
        val[f"N{tindex[t]}"] += 1; prk[(nrm, c)] = (r, k)
    for j in range(1, D + 1): val[f"m{j}"] = sum(1 for p in S if deg[p] == j)
    Sl = sorted(S)
    for a_, b_ in itertools.combinations(Sl, 2):
        i, j = sorted((deg[a_], deg[b_])); val[f"M{i},{j}"] += 1
    for di, pts in ag.lines:
        jj = len(pts & S)
        if jj: val[f"L{jj}"] += 1
    # classes
    ctypes = info["_ctypes"]; cindex = {t: i for i, t in enumerate(ctypes)}
    cls = {}
    for (nrm, c), rk in prk.items(): cls.setdefault(nrm, []).append(rk)
    for nrm, lst in cls.items():
        t = tuple(sorted(lst)); assert t in cindex, ("class type missing", t); val[f"C{cindex[t]}"] += 1
    # directions: a_{d,L} from geometry. class L (normal nrm) contains direction d iff dot(nrm, d) = 0
    utypes = info["_utypes"]; uindex = {u: i for i, u in enumerate(utypes)}
    wvecs = info["_wvecs"]; windex = {w: i for i, w in enumerate(wvecs)}
    wdir = {ag.lines[l][0]: l for l in W}
    for di, dvec in enumerate(ag.dirs):
        ld = wdir[di]; a_list = []
        for nrm in ag.dirs:
            if dot(F, nrm, dvec) != 0: continue
            cnt = 0
            for (nrm2, c, pts) in planes:
                if nrm2 != nrm: continue
                if any(ag.lines[l][1] <= pts for l in W if l != ld): cnt += 1
            a_list.append(cnt)
        u = tuple(sorted(a_list)); assert len(u) == q + 1 and u in uindex, u; val[f"A{uindex[u]}"] += 1
        Lv = [0] * q
        for l in ag.lines_by_dir[di]:
            jj = len(ag.lines[l][1] & S)
            if jj: Lv[jj - 1] += 1
        w = tuple(Lv); assert w in windex, w; val[f"P{windex[w]}"] += 1
    # witness lines
    items = info["_items"]; wt = info["_wt"]; wtindex = {t: i for i, t in enumerate(wt)}; dms = info["_dms"]; dmindex = {t: i for i, t in enumerate(dms)}
    item_of = {}; wtype_of_line = {}
    for ii, (r, lo, hi, kk) in enumerate(items):
        for k in kk: item_of[(r, k)] = ii
    for l in W:
        lpts = ag.lines[l][1]; its = []
        for (nrm, c, pts) in planes:
            if lpts <= pts: its.append(item_of[prk[(nrm, c)]])
        t = tuple(sorted(its)); assert t in wtindex, ("wline type missing", t); val[f"W{wtindex[t]}"] += 1
        dm = tuple(sorted(deg[p] for p in lpts)); val[f"V{dmindex[dm]}"] += 1
        wtype_of_line[l] = wtindex[t]
        if "_Y" in info:   # F8 (session VII): (plane histogram type, dm) incidences and (dm, pencil) pairs
            if "_Yidx" not in info: info["_Yidx"] = {k: i for i, k in enumerate(info["_Y"])}; info["_Zidx"] = {k: i for i, k in enumerate(info["_Z"])}
            di = dmindex[dm]
            zk = (di, wtindex[t]); assert zk in info["_Zidx"], ("dm/pencil pair missing", dm, t); val[f"Z{info['_Zidx'][zk]}"] += 1
            for (nrm, c, pts) in planes:
                if lpts <= pts:
                    h = [0] * D
                    for p in pts & S: h[deg[p] - 1] += 1
                    yk = ((prk[(nrm, c)][0], prk[(nrm, c)][1], tuple(h)), di); assert yk in info["_Yidx"], ("Y pair missing", yk); val[f"G{info['_Yidx'][yk]}"] += 1
    if p6.SECTIONS and "_X1" in info:   # F6 flags from geometry
        X1index = {key: i for i, key in enumerate(info["_X1"])}; X2index = {key: i for i, key in enumerate(info["_X2"])}
        for (nrm, c, pts, r, k, pr, ed) in flag_planes:
            lpts = ag.lines[wdir[ed]][1]
            mate = [prk[(nrm2, c2)] for (nrm2, c2, pts2) in planes if nrm2 == nrm and lpts <= pts2]
            assert len(mate) == 1 and mate[0][0] >= 1 and not lpts <= pts, mate
            rk2 = mate[0]; wi = wtype_of_line[wdir[ed]]
            kd = (r, k, pr, rk2); assert kd in X1index, ("flag kind missing", kd); val[f"X{X1index[kd]}"] += 1
            key = ((pr, rk2), wi); assert key in X2index, ("flag/pencil incompatible", key); val[f"Y{X2index[key]}"] += 1
    # check rows
    bad = 0
    for (coefs, rhs, sense, name) in lp.rows:
        if name == 'm_D >= 1' and actual_D < D:continue  # Exact-max premise, not a universal row.
        lhs = sum(c * val[lp.names[v]] for v, c in coefs.items())
        ok = (lhs == rhs) if sense == '=' else (lhs <= rhs)
        if not ok: bad += 1; print("  VIOLATED:", name, "lhs", lhs, sense, rhs)
    print(f"{fn}: n={n} D={D} rows={len(lp.rows)} violated={bad}  (plane types used {sum(1 for k,v in val.items() if k[0]=='N' and v)})")
    if actual_D < D:print(f'  padded cap {D}: universal rows checked; exact-max premise m_D>=1 inactive (actual D={actual_D})',flush=True)
    return bad == 0

if __name__ == "__main__":
    allok = all([validate(fn) for fn in sys.argv[1:]])
    print("ALL OK" if allok else "FAILURES")
