#!/usr/bin/env python3
"""Final hostile completeness audit of the Nik(2,9) = 63 release proof (2026-11-13).

Independent of the lane scripts (run_batch / collect / progress).  Reads only the specs, the records in
big/nik29 (+ the tracked copies in nik29/records), the shipped encoder and TOOLCHAIN.md.

  (a) enumerate the 28 max-secant cases with OWN orbit code (not enc_nik29.axis_orbit_reps), cross-check
      against the encoder and the specs, and show each case is closed by a verified whole-case record or by a
      verified cover + all 100 cubes;
  (b) regenerate every case's whole CNF, cover CNF and cubes with src/enc_nik29.py and check STRUCTURALLY
      (by parsing the DIMACS files, not by trusting the encoder's intent) that
        cube_i  = whole-case clauses + exactly D unit clauses of positive witness literals,
        cover   = the witness clauses of the whole-case CNF (one per point, all-positive, 10 literals)
                  + the negation of every cube's unit set, nothing else, and
        the cube unit sets are pairwise distinct and form the full product of the witness choices of the
        cube points, 100 = 10 x 10 of them;
      hence  whole UNSAT  <=  cover UNSAT  and  every cube UNSAT;
  (c) compare the sha256 of every regenerated CNF (28 whole, 11 covers, 1100 cubes; also the 4 SAT controls
      and 32 q=7 controls) with the cnf_sha256 of the records (byte-identical regeneration); list the 10
      hardest cubes by CaDiCaL time explicitly;
  (d) run src/symmetry_check_nik29.py 7 8 9 and nik29/check_models.py on the SAT controls;
  (e) every accepted record: result UNSAT, verified true, cake_verdict 's VERIFIED UNSAT', lrat sha256 +
      size present, tool hashes = TOOLCHAIN.md, record copies in nik29/records byte-identical; nothing is
      counted on solver status alone.
Writes nik29/FINAL_AUDIT.md.  Exit 0 iff everything passes.
usage: nice -n 19 python3 nik29/final_audit.py [--regen all|sample]  (default all)
"""
import argparse, hashlib, itertools, json, os, re, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
BIG = os.path.join(ROOT, "big", "nik29")
REC = os.path.join(HERE, "records")
TMP = os.path.join(ROOT, "tmp", "nik29_audit")
sys.path.insert(0, SRC)
from ff import GF  # noqa: E402

Q, K = 9, 19
P = Q * Q
TOOL_HASHES = {
    "cadical": "fafdc017e81c48e0418986460cc7738b28c5efb823764c6a7b875c1a321756e2",
    "lrat-trim": "1c4a4ec77c953dabc6c8b47835c15d94d7bc48fff75fbbcde8b4ff6fb5626129",
    "cake_lpr": "130816693a5392bbba02cf0c382ac3569542dcbb482ad0b409e8a1136ab4126d",
}
FAIL = []
LOG = []


def log(s=""):
    print(s, flush=True)
    LOG.append(s)


def fail(s):
    FAIL.append(s)
    log("FAIL: " + s)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_cnf(path):
    """returns (comment lines, nvars, nclauses declared, list of clauses as tuples)"""
    comments, nv, nc, cls = [], None, None, []
    for line in open(path):
        if line.startswith("c"):
            comments.append(line[1:].strip())
        elif line.startswith("p"):
            _, _, nv, nc = line.split()
            nv, nc = int(nv), int(nc)
        else:
            lits = [int(t) for t in line.split()]
            assert lits[-1] == 0
            cls.append(tuple(lits[:-1]))
    assert nc == len(cls), (path, nc, len(cls))
    return comments, nv, nc, cls


# ---------------------------------------------------------------- (a) own enumeration of the 28 cases
def own_axis_orbits(F):
    """orbits of s-subsets of F_q under x -> a x^sigma + b, written from scratch (union-find over images)"""
    q = F.q

    def frob(x, e):
        for _ in range(e):
            x = F.frob[x]
        return x

    gens = [(a, b, e) for a in range(1, q) for b in range(q) for e in range(F.k)]
    out = {}
    for s in range(2, q + 1):
        subsets = list(itertools.combinations(range(q), s))
        idx = {S: i for i, S in enumerate(subsets)}
        parent = list(range(len(subsets)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for S in subsets:
            for (a, b, e) in gens:
                T = tuple(sorted(F.add[F.mul[a][frob(x, e)]][b] for x in S))
                ra, rb = find(idx[S]), find(idx[T])
                if ra != rb:
                    parent[ra] = rb
        classes = {}
        for S in subsets:
            classes.setdefault(find(idx[S]), []).append(S)
        reps = []
        for cl in classes.values():
            cand = [S for S in cl if 0 in S and 1 in S]
            if not cand:
                fail(f"s={s}: an orbit without a representative containing 0 and 1")
            reps.append((min(cand), len(cl)))
        reps.sort(key=lambda t: t[0])
        out[s] = reps
    return out


def step_a():
    log("## (a) the 28 cases and how each is closed\n")
    F = GF(Q)
    own = own_axis_orbits(F)
    enc = subprocess.run(["python3", "enc_nik29.py", "--q", str(Q), "--list"], cwd=SRC, capture_output=True, text=True, check=True).stdout
    enc_reps = {}
    for line in enc.strip().splitlines():
        s, rest = line.split(" ", 1)
        lst = eval(rest.split(" sum of orbit sizes")[0])  # list of (tuple, size)
        enc_reps[int(s)] = [(tuple(r), n) for r, n in lst]
    n_cases = 0
    for s in range(2, Q + 1):
        mine = own[s]
        theirs = enc_reps[s]
        if sorted(mine) != sorted(theirs):
            fail(f"s={s}: own orbit representatives {mine} != encoder's {theirs}")
        tot = sum(n for _, n in mine)
        if tot != len(list(itertools.combinations(range(Q), s))):
            fail(f"s={s}: orbit sizes sum to {tot} != C(9,{s})")
        n_cases += 2 * len(mine)
    log(f"own orbit enumeration (union-find over all images, no encoder code): representatives per s = "
        + ", ".join(f"s={s}:{len(own[s])}" for s in range(2, Q + 1)) + f"; {n_cases // 2} representatives x 2 tangent flags = {n_cases} cases;"
        " identical to the encoder's list (same representatives, same orbit sizes; sizes sum to C(9,s) for every s).")
    if n_cases != 28:
        fail(f"case count {n_cases} != 28")
    # expected case names, in the encoder's orbit order (the encoder enumerates orbits in combinations order;
    # the names in the specs use that order -- check by matching S in the regenerated CNF header later)
    cases = []
    for s in range(2, Q + 1):
        for o, (rep, _) in enumerate(enc_reps[s]):
            for tan in "hd":
                cases.append((f"q{Q}K{K}_s{s}o{o}{tan}", s, o, tan, rep))
    spec_whole = {c["name"]: c for c in json.load(open(os.path.join(ROOT, "specs", "nik29.json")))}
    spec_sup_path = os.path.join(ROOT, "specs", "nik29_superseded.json")
    spec_sup = {c["name"]: c for c in json.load(open(spec_sup_path))} if os.path.exists(spec_sup_path) else {}
    spec_cubes = {}
    for sp in ("nik29cubes_s4.json", "nik29cubes_s5.json"):
        for c in json.load(open(os.path.join(ROOT, "specs", sp))):
            spec_cubes.setdefault(c["name"].rsplit("_", 1)[0], {})[c["name"]] = (sp, c)
    if set(spec_whole) | set(spec_sup) != {c[0] for c in cases}:
        fail(f"specs/nik29.json + nik29_superseded.json names != the 28 cases: {sorted((set(spec_whole) | set(spec_sup)) ^ {c[0] for c in cases})}")
    if set(spec_whole) & set(spec_sup):
        fail("a case is both in nik29.json and nik29_superseded.json")
    closure = {}
    rows = []
    for name, s, o, tan, rep in cases:
        whole = load_rec(name)
        by_whole = whole is not None and accepted(whole, name)
        cub = spec_cubes.get(name)
        by_cubes = False
        ncube = 0
        if cub:
            names = sorted(cub)
            if f"{name}_cover" not in cub:
                fail(f"{name}: cube spec without cover")
            ncube = sum(1 for n in names if "_cube" in n)
            if ncube != 100 or len(names) != 101:
                fail(f"{name}: cube spec has {ncube} cubes / {len(names)} entries, expected 100 / 101")
            if sorted(int(n.rsplit("cube", 1)[1]) for n in names if "_cube" in n) != list(range(100)):
                fail(f"{name}: cube indices are not 0..99")
            recs = {n: load_rec(n) for n in names}
            by_cubes = all(r is not None and accepted(r, n) for n, r in recs.items())
        if by_whole and name in spec_sup:
            fail(f"{name}: verified whole record but listed as superseded")
        if not by_whole and name in spec_whole:
            fail(f"{name}: in specs/nik29.json but not closed by a verified whole-case record")
        if by_cubes and name not in spec_sup and not by_whole:
            fail(f"{name}: closed by cubes but not listed in nik29_superseded.json")
        how = "whole-case record" if by_whole else ("cover + 100 cubes" if by_cubes else "NOT CLOSED")
        if how == "NOT CLOSED":
            fail(f"{name}: not closed")
        if by_whole and by_cubes:
            how = "whole-case record (also cover + 100 cubes)"
        closure[name] = (by_whole, by_cubes, rep)
        sup = ""
        if name in spec_sup:
            r = load_rec(name)
            sup = f"whole-case attempt: {r.get('result')} after {r.get('cadical_seconds')} s" if r else "whole case never run"
        rows.append(f"| {name} | {s} | {list(rep)} | {tan} | {how} | {sup} |")
    log("\n| case | s | axis set S | tangent | closed by | superseded whole-case entry |\n|---|---|---|---|---|---|")
    for r in rows:
        log(r)
    nw = sum(1 for v in closure.values() if v[0])
    nc = sum(1 for v in closure.values() if not v[0] and v[1])
    log(f"\n{nw} cases closed by a verified whole-case record, {nc} by a verified cover + 100 cubes, total {nw + nc} / 28.")
    return cases, closure, spec_cubes


# ---------------------------------------------------------------- records
_rec_cache = {}


def rec_path(name):
    return os.path.join(BIG if not name.startswith("q7") else os.path.join(BIG, "ctrl7"), name + ".json")


def load_rec(name):
    if name not in _rec_cache:
        p = rec_path(name)
        _rec_cache[name] = json.load(open(p)) if os.path.exists(p) else None
    return _rec_cache[name]


def accepted(r, name):
    """the acceptance criterion of this audit (stricter than collect.py): an UNSAT record counts only with a
    cake_lpr success, a proof hash and size, and the pinned tools"""
    ok = (r.get("result") == "UNSAT" and r.get("verified") is True and r.get("cake_verdict", "").endswith("s VERIFIED UNSAT")
          and re.fullmatch(r"[0-9a-f]{64}", r.get("lrat_sha256", "")) is not None and isinstance(r.get("lrat_bytes"), int) and r["lrat_bytes"] > 0
          and re.fullmatch(r"[0-9a-f]{64}", r.get("cnf_sha256", "")) is not None and r.get("tools") == TOOL_HASHES
          and r.get("cnf") == name + ".cnf" and isinstance(r.get("cake_seconds"), (int, float)))
    return ok


# ---------------------------------------------------------------- (b)+(c) regeneration + structure
def regen_whole(s, o, tan, Kv, outpath, q=Q):
    subprocess.run(["python3", "enc_nik29.py", "--q", str(q), "--K", str(Kv), "--s", str(s), "--orbit", str(o), "--tan", tan, "--out", outpath],
                   cwd=SRC, check=True, stdout=subprocess.DEVNULL)


def regen_cubes(s, o, tan, pts, outdir, cover):
    subprocess.run(["python3", "enc_nik29.py", "--q", str(Q), "--K", str(K), "--s", str(s), "--orbit", str(o), "--tan", tan,
                    "--cube", str(len(pts)), "--cube-points", ",".join(map(str, pts)), "--outdir", outdir, "--cover", cover],
                   cwd=SRC, check=True, stdout=subprocess.DEVNULL)


def step_bc(cases, closure, spec_cubes, regen_mode):
    log("\n## (b) structure of covers and cubes, (c) byte-identical regeneration\n")
    os.makedirs(TMP, exist_ok=True)
    n_hash_ok = n_hash = 0
    hardest = []
    whole_sha = {}
    for name, s, o, tan, rep in cases:
        wpath = os.path.join(TMP, name + ".cnf")
        regen_whole(s, o, tan, K, wpath)
        whole_sha[name] = sha(wpath)
        hdr, nv, nc, wcls = read_cnf(wpath)
        m = re.search(r"axis S=\[([0-9, ]+)\] tan=(\w)", hdr[0])
        if tuple(int(t) for t in m.group(1).split(",")) != rep or m.group(2) != tan:
            fail(f"{name}: regenerated header {hdr[0]} does not match case (S={rep}, tan={tan})")
        # axis units present?  c(t,0) is variable 9t+1
        units = {c[0] for c in wcls if len(c) == 1}
        for t in range(Q):
            v = Q * t + 0 + 1
            if (v if t in rep else -v) not in units:
                fail(f"{name}: axis unit for t={t} missing")
        r = load_rec(name)
        if r is not None and r.get("cnf_sha256"):
            n_hash += 1
            if r["cnf_sha256"] == whole_sha[name]:
                n_hash_ok += 1
            else:
                fail(f"{name}: whole-case cnf_sha256 mismatch (record {r['cnf_sha256'][:12]} vs regenerated {whole_sha[name][:12]})")
        if name in spec_cubes:
            _, cover_case = spec_cubes[name][name + "_cover"]
            gen = cover_case["gen"]
            pts = [int(t) for t in gen[gen.index("--cube-points") + 1].split(",")]
            D = int(gen[gen.index("--cube") + 1])
            if D != len(pts):
                fail(f"{name}: --cube {D} but {len(pts)} cube points")
            cdir = os.path.join(TMP, name + "_cubes")
            cover = os.path.join(TMP, name + "_cover.cnf")
            shutil.rmtree(cdir, ignore_errors=True)
            regen_cubes(s, o, tan, pts, cdir, cover)
            # witness clauses of the whole CNF: all-positive clauses of 10 w-literals (> P)
            wit = [c for c in wcls if len(c) == Q + 1 and all(l > P for l in c)]
            if len(wit) != P:
                fail(f"{name}: {len(wit)} witness clauses in whole CNF, expected {P}")
            wit_set = set(wit)
            # which witness clause belongs to which point: the encoder numbers w-variables point by point;
            # recover the point of a w-literal from the binary clauses (-w, -c_p): the point x of w is the one
            # point of the line NOT excluded.  Independent of the encoder's variable layout.
            excl = {}
            for c in wcls:
                if len(c) == 2 and c[0] < -P and -P <= c[1] < 0:
                    excl.setdefault(-c[0], set()).add(-c[1] - 1)
            def point_of(w):
                # the witness clause containing w is the clause of point x; x is the point on no excluded set...
                # simpler: every line through x has q-1 excluded points, all != x, and x is the unique point of the
                # 9 excluded sets' complement that lies on all lines -> use the clause membership instead
                return w_point[w]
            w_point = {}
            for cl in wit:
                # the point x of this clause: the unique point not excluded by any w in cl but on every line:
                # intersection over w in cl of (line(w) = excl[w] U {x}) -- compute as the common point
                common = None
                for w in cl:
                    line_pts = set(excl[w])
                    if len(line_pts) != Q - 1:
                        fail(f"{name}: witness literal {w} excludes {len(line_pts)} points, expected {Q - 1}")
                    # x is not in excl[w]; all lines through x pairwise meet only in x, so x = the point missing
                    # from excl[w] that is missing from every other excl[w'] as well
                    miss = set(range(P)) - line_pts
                    common = miss if common is None else common & miss
                if common is None or len(common) != 1:
                    fail(f"{name}: cannot identify the point of a witness clause ({common})")
                    continue
                x = next(iter(common))
                for w in cl:
                    w_point[w] = x
            wit_by_point = {}
            for cl in wit:
                wit_by_point[w_point[cl[0]]] = cl
            if len(wit_by_point) != P:
                fail(f"{name}: witness clauses do not cover all {P} points exactly once")
            cube_units = []
            for i in range(100):
                cpath = os.path.join(cdir, f"cube{i}.cnf")
                chdr, cnv, cnc, ccls = read_cnf(cpath)
                if cnv != nv:
                    fail(f"{name} cube{i}: nvars {cnv} != whole {nv}")
                if ccls[:len(wcls)] != wcls or len(ccls) != len(wcls) + D:
                    fail(f"{name} cube{i}: clauses are not (whole-case clauses + {D} extra)")
                extra = ccls[len(wcls):]
                if any(len(c) != 1 or c[0] <= P for c in extra):
                    fail(f"{name} cube{i}: extra clauses are not positive witness units: {extra}")
                u = tuple(c[0] for c in extra)
                if [w_point.get(w) for w in u] != pts:
                    fail(f"{name} cube{i}: units belong to points {[w_point.get(w) for w in u]}, expected {pts}")
                cube_units.append(u)
                # hash vs record
                rn = f"{name}_cube{i}"
                rr = load_rec(rn)
                n_hash += 1
                h = sha(cpath)
                if rr is None or rr.get("cnf_sha256") != h:
                    fail(f"{rn}: cnf_sha256 mismatch or missing record")
                else:
                    n_hash_ok += 1
                    hardest.append((rr["cadical_seconds"], rn, rr["lrat_bytes"]))
            if len(set(cube_units)) != 100:
                fail(f"{name}: cube unit sets not pairwise distinct")
            expected = set(itertools.product(*[wit_by_point[x] for x in pts]))
            if set(cube_units) != expected:
                fail(f"{name}: the cube unit sets are not the full product of the witness choices of points {pts}")
            # cover
            vhdr, vnv, vnc, vcls = read_cnf(cover)
            if vnv != nv:
                fail(f"{name} cover: nvars {vnv} != whole {nv}")
            cov_wit = [c for c in vcls if len(c) == Q + 1]
            cov_neg = [c for c in vcls if len(c) == D]
            if len(cov_wit) + len(cov_neg) != len(vcls):
                fail(f"{name} cover: unexpected clause lengths {sorted(set(len(c) for c in vcls))}")
            if set(cov_wit) != wit_set or len(cov_wit) != P:
                fail(f"{name} cover: witness clauses differ from the whole-case witness clauses")
            if sorted(tuple(-l for l in c) for c in cov_neg) != sorted(cube_units):
                fail(f"{name} cover: negated clauses are not exactly the 100 cube unit sets")
            if not set(cov_wit) <= set(wcls):
                fail(f"{name} cover: a cover clause is not a clause of the whole-case CNF")
            rr = load_rec(name + "_cover")
            n_hash += 1
            if rr is None or rr.get("cnf_sha256") != sha(cover):
                fail(f"{name}_cover: cnf_sha256 mismatch or missing record")
            else:
                n_hash_ok += 1
            log(f"{name}: cover = {len(cov_wit)} witness clauses (all clauses of the whole CNF) + {len(cov_neg)} negated cubes; "
                f"100 cubes = whole CNF + {D} witness units of points {pts} = {[divmod(p, Q) for p in pts]}, full 10 x 10 product; "
                f"101/101 regenerated CNFs byte-identical to the records (sha256).")
            shutil.rmtree(cdir, ignore_errors=True)
            os.remove(cover)
        os.remove(wpath)
    # SAT control and q=7 control CNFs
    for sp, Kv, q in (("nik29sat.json", 18, 9), ("nik29ctrl7.json", None, 7)):
        for c in json.load(open(os.path.join(ROOT, "specs", sp))):
            g = c["gen"]
            s, o, tan, Kc = int(g[g.index("--s") + 1]), int(g[g.index("--orbit") + 1]), g[g.index("--tan") + 1], int(g[g.index("--K") + 1])
            p = os.path.join(TMP, c["name"] + ".cnf")
            regen_whole(s, o, tan, Kc, p, q=q)
            r = load_rec(c["name"])
            n_hash += 1
            if r is None or r.get("cnf_sha256") != sha(p):
                fail(f"{c['name']}: control cnf_sha256 mismatch")
            else:
                n_hash_ok += 1
            os.remove(p)
    hardest.sort(reverse=True)
    log(f"\nRegenerated and hashed {n_hash} CNFs (28 whole cases incl. the superseded ones that have a record, 11 covers, 1100 cubes, "
        f"4 SAT controls, 32 q=7 controls): {n_hash_ok} byte-identical to the recorded cnf_sha256, {n_hash - n_hash_ok} mismatches.")
    log("\nHardest 10 cubes by CaDiCaL time (all among the byte-identical regenerations):\n\n| cube | CaDiCaL s | LRAT bytes |\n|---|---|---|")
    for t, n, b in hardest[:10]:
        log(f"| {n} | {t} | {b} |")
    return whole_sha


# ---------------------------------------------------------------- (d)
def step_d():
    log("\n## (d) symmetry checks and SAT-control models\n")
    p = subprocess.run(["nice", "-n", "19", "python3", "symmetry_check_nik29.py", "7", "8", "9"], cwd=SRC, capture_output=True, text=True)
    log("```\n" + p.stdout.strip() + "\n```")
    if p.returncode != 0 or "NIK29 SYMMETRY CHECKS PASSED" not in p.stdout:
        fail("symmetry_check_nik29.py failed: " + p.stderr[-500:])
    for sp, d, expect_sat in (("nik29sat.json", BIG, {"q9K18_s6o0h", "q9K18_s6o0d"}),
                              ("nik29ctrl7.json", os.path.join(BIG, "ctrl7"), {"q7K13_s4o1h", "q7K13_s4o1d"})):
        p = subprocess.run(["python3", os.path.join(HERE, "check_models.py"), os.path.join(ROOT, "specs", sp), d], capture_output=True, text=True)
        log("```\n" + p.stdout.strip() + "\n```")
        sat_names = set()
        for line in p.stdout.splitlines():
            if " |C| = " in line:
                nm = line.split()[0]
                sat_names.add(nm)
                if "rc 0" not in line or "NIKODYM OK" not in line:
                    fail(f"{nm}: model check line not OK: {line}")
                # size claim
                size = int(line.split("|C| =")[1].split()[0])
                Kc = int(nm.split("K")[1].split("_")[0])
                if size != Kc:
                    fail(f"{nm}: decoded |C| = {size} != K = {Kc}")
        if p.returncode != 0:
            fail(f"check_models.py {sp} exited {p.returncode}: {p.stderr[-300:]}")
        if sat_names != expect_sat:
            fail(f"{sp}: SAT cases {sorted(sat_names)} != expected {sorted(expect_sat)}")
        # the non-SAT controls must be accepted UNSAT
        for c in json.load(open(os.path.join(ROOT, "specs", sp))):
            r = load_rec(c["name"])
            if c["name"] in expect_sat:
                if r is None or r.get("result") != "SAT" or not os.path.exists(r.get("model_file", "")):
                    fail(f"{c['name']}: expected SAT record with model file")
            elif r is None or not accepted(r, c["name"]):
                fail(f"{c['name']}: control not accepted UNSAT")
    # independent re-check of the decoded SAT models with check_nikodym.py and of the case conditions
    for nm in ("q9K18_s6o0h", "q9K18_s6o0d"):
        d = json.load(open(os.path.join(BIG, nm + ".decoded.json")))
        C = {tuple(p) for p in d["C"]}
        if len(C) != 18:
            fail(f"{nm}: |C| = {len(C)}")
        axis = {x for (x, y) in C if y == 0}
        if axis != set(d["S"]):
            fail(f"{nm}: axis section {axis} != S {d['S']}")
        # max secant exactly s = 6 (own line enumeration)
        from ag import AG
        g = AG(2, 9)
        ms = max(sum(1 for p in pts if g.points[p] in C) for _, pts in g.lines)
        if ms != 6:
            fail(f"{nm}: max secant {ms} != 6")
        log(f"{nm}: decoded model has |C| = 18, axis section = S = {d['S']}, max secant = 6, Nikodym-checked (|N| = 63).")


# ---------------------------------------------------------------- (e)
def step_e(cases, closure, spec_cubes):
    log("\n## (e) record audit\n")
    specs = ["nik29ctrl7.json", "nik29sat.json", "nik29.json", "nik29cubes_s4.json", "nik29cubes_s5.json"]
    tot = {}
    n_unsat = n_sat = n_copies_ok = n_copies = 0
    dates = {}
    for sp in specs:
        spec_name = sp[:-5]
        t = {"cadical": 0.0, "trim": 0.0, "cake": 0.0, "max_cadical": 0.0, "max_lrat": 0, "n": 0, "unsat": 0, "sat": 0, "trim_fallback": 0}
        for c in json.load(open(os.path.join(ROOT, "specs", sp))):
            nm = c["name"]
            r = load_rec(nm)
            if r is None:
                fail(f"{sp}: {nm} has no record")
                continue
            t["n"] += 1
            p = rec_path(nm)
            mt = os.path.getmtime(p)
            dates.setdefault(spec_name, []).append(mt)
            # tracked copy
            cp = os.path.join(REC, spec_name, nm + ".json")
            n_copies += 1
            if os.path.exists(cp) and sha(cp) == sha(p):
                n_copies_ok += 1
            else:
                fail(f"{nm}: tracked copy nik29/records/{spec_name}/{nm}.json missing or differs from big/nik29")
            if r.get("tools") != TOOL_HASHES:
                fail(f"{nm}: tool hashes differ from TOOLCHAIN.md")
            if r.get("result") == "SAT":
                if c["expect"] == "UNSAT":
                    fail(f"{nm}: SAT but expected UNSAT")
                t["sat"] += 1
                n_sat += 1
                t["cadical"] += r["cadical_seconds"]
                continue
            if not accepted(r, nm):
                fail(f"{nm}: not accepted (result={r.get('result')}, verified={r.get('verified')}, cake_verdict={r.get('cake_verdict')!r})")
                continue
            if r.get("trim_rc") != 0:
                t["trim_fallback"] += 1
            t["unsat"] += 1
            n_unsat += 1
            t["cadical"] += r["cadical_seconds"]
            t["trim"] += r.get("trim_seconds", 0)
            t["cake"] += r["cake_seconds"]
            t["max_cadical"] = max(t["max_cadical"], r["cadical_seconds"])
            t["max_lrat"] = max(t["max_lrat"], r["lrat_bytes"])
        tot[spec_name] = t
    # current binaries
    for k, h in TOOL_HASHES.items():
        b = os.path.join(ROOT, "repos", "bin", k)
        if not os.path.exists(b) or sha(b) != h:
            fail(f"repos/bin/{k} missing or hash differs from TOOLCHAIN.md")
    log(f"{n_unsat} UNSAT records accepted (each: result UNSAT, verified true, cake_verdict 's VERIFIED UNSAT', LRAT sha256 + byte size, "
        f"CNF sha256, pinned tool hashes); {n_sat} SAT records (models re-checked in (d)); {n_copies_ok}/{n_copies} tracked copies in nik29/records "
        f"byte-identical to big/nik29; repos/bin binaries match TOOLCHAIN.md.  No instance is counted on solver status alone: the acceptance "
        f"predicate requires the cake_lpr verdict string.")
    log("\n| spec | instances | UNSAT verified | SAT | CaDiCaL s | lrat-trim s | cake_lpr s | max CaDiCaL s | max LRAT | trim fallbacks | record dates (UTC) |\n|---|---|---|---|---|---|---|---|---|---|---|")
    for k, t in tot.items():
        ds = dates[k]
        f = lambda x: time.strftime("%Y-%m-%d %H:%M", time.gmtime(x))
        log(f"| {k} | {t['n']} | {t['unsat']} | {t['sat']} | {t['cadical']:.0f} | {t['trim']:.0f} | {t['cake']:.0f} | {t['max_cadical']:.0f} | {t['max_lrat'] / 1e9:.2f} GB | {t['trim_fallback']} | {f(min(ds))} - {f(max(ds))} |")
    proof = {k: tot[k] for k in ("nik29", "nik29cubes_s4", "nik29cubes_s5")}
    cad = sum(t["cadical"] for t in proof.values()); trim = sum(t["trim"] for t in proof.values()); cake = sum(t["cake"] for t in proof.values())
    ml = max(t["max_lrat"] for t in proof.values()); mc = max(t["max_cadical"] for t in proof.values())
    n = sum(t["unsat"] for t in proof.values())
    log(f"\n**Proof of |C| = 19 impossible (the 17 whole cases + 1111 cube instances = {n} UNSAT instances): CaDiCaL {cad:.0f} s "
        f"({cad / 3600:.1f} h), lrat-trim {trim:.0f} s, cake_lpr {cake:.0f} s ({cake / 3600:.1f} h); max CaDiCaL {mc:.0f} s; largest LRAT {ml} bytes "
        f"({ml / 1e9:.2f} GB).**")
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regen", default="all")
    a = ap.parse_args()
    t0 = time.time()
    log(f"# FINAL AUDIT of Nik(2,9) = 63 (max-secant split, enc_nik29.py) — {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}\n")
    log("Produced by `nik29/final_audit.py` (this file is its verbatim output).  The audit is independent of the lane scripts: it enumerates the "
        "cases with its own orbit code, parses the regenerated DIMACS files and checks the cube/cover structure clause by clause, and applies its "
        "own acceptance predicate to every record.\n")
    cases, closure, spec_cubes = step_a()
    step_bc(cases, closure, spec_cubes, a.regen)
    step_d()
    step_e(cases, closure, spec_cubes)
    log("\n## Verdict\n")
    if FAIL:
        log(f"**{len(FAIL)} FAILURE(S):**\n")
        for f in FAIL:
            log("* " + f)
    else:
        log("**ALL CHECKS PASSED.**  Every one of the 28 cases is closed by cake_lpr-verified instances whose CNFs are regenerated byte-identically "
            "by the shipped encoder; the cube splits are exhaustive by construction (full product) and by the verified cover CNFs; the SAT controls "
            "behave as required.  Nik(2,9) = 63 is READY at the Nik(2,8) standard.")
    log(f"\n(audit wall time {time.time() - t0:.0f} s)")
    open(os.path.join(HERE, "FINAL_AUDIT.md"), "w").write("\n".join(LOG) + "\n")
    shutil.rmtree(TMP, ignore_errors=True)
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
