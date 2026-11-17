"""Cardinality constraints as CNF clauses (fresh implementation, no library).

atmost(lits, k, top) encodes sum(lits) <= k with the sequential counter of
Sinz (CP 2005): auxiliary variables s[i][j] = "at least j of the first i literals are true".
Returns (clauses, new_top).  Literals are nonzero ints (DIMACS convention).
Correctness: every assignment of the lits with <= k true extends to the s variables
(set s[i][j] = [count of first i >= j]); conversely the clauses force s[i][j] whenever
j of the first i are true, and the last clause family forbids a (k+1)-st true literal.
"""


def atmost(lits, k, top):
    lits = list(lits)
    n = len(lits)
    clauses = []
    if k >= n:
        return clauses, top
    if k == 0:
        return [[-l] for l in lits], top
    s = [[0] * (k + 1) for _ in range(n)]  # s[i][j], 1<=j<=k
    for i in range(n):
        for j in range(1, k + 1):
            top += 1
            s[i][j] = top
    for i in range(n):
        x = lits[i]
        clauses.append([-x, s[i][1]])                       # x_i -> s_i1
        if i > 0:
            for j in range(1, k + 1):
                clauses.append([-s[i - 1][j], s[i][j]])     # s_{i-1,j} -> s_ij
            for j in range(2, k + 1):
                clauses.append([-x, -s[i - 1][j - 1], s[i][j]])  # x_i & s_{i-1,j-1} -> s_ij
            clauses.append([-x, -s[i - 1][k]])              # x_i & s_{i-1,k} -> false
    return clauses, top


def atleast(lits, k, top):
    lits = list(lits)
    return atmost([-l for l in lits], len(lits) - k, top)


def equals(lits, k, top):
    c1, top = atmost(lits, k, top)
    c2, top = atleast(lits, k, top)
    return c1 + c2, top


def write_dimacs(path, nvars, clauses, comments=()):
    with open(path, "w") as f:
        for c in comments:
            f.write(f"c {c}\n")
        f.write(f"p cnf {nvars} {len(clauses)}\n")
        for cl in clauses:
            f.write(" ".join(map(str, cl)) + " 0\n")


if __name__ == "__main__":
    # exhaustive self-test on small sizes with a tiny DPLL
    import itertools

    def sat(nv, clauses, fixed):
        # fixed: dict var->bool for the first variables; brute force the rest? use simple backtracking
        assign = dict(fixed)

        def rec(cls):
            cls2 = []
            for c in cls:
                if any(assign.get(abs(l), None) == (l > 0) for l in c):
                    continue
                c2 = [l for l in c if abs(l) not in assign]
                if not c2:
                    return False
                cls2.append(c2)
            if not cls2:
                return True
            # unit propagation / branch on first literal
            c = min(cls2, key=len)
            l = c[0]
            for val in (l > 0, l < 0):
                assign[abs(l)] = val
                if rec(cls2):
                    return True
                del assign[abs(l)]
            return False
        return rec(clauses)

    for n in range(1, 6):
        for k in range(0, n + 1):
            cls, top = atmost(list(range(1, n + 1)), k, n)
            cls2, top2 = atleast(list(range(1, n + 1)), k, n)
            for bits in itertools.product((False, True), repeat=n):
                fixed = {i + 1: b for i, b in enumerate(bits)}
                assert sat(top, cls, fixed) == (sum(bits) <= k), (n, k, bits)
                assert sat(top2, cls2, fixed) == (sum(bits) >= k), (n, k, bits)
    print("card.py: sequential counter exhaustively verified for n<=5")
