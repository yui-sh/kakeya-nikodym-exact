"""Test Singer-cycle / norm-type candidates in F_{q^3} = F_q^3 as Kakeya sets.
Hypothesis from the q=3 optimum (|Stab|=39=13*3): the set is a coset of the order-(q^2+q+1)
subgroup of F_{q^3}^* (norm-1 elements), invariant under Frobenius.
Lines: {x + t d : t in F_q}. Directions: F_{q^3}^* / F_q^*."""
import sys, itertools
from gf import GF, _IRRED
_IRRED[64] = (2, [1, 1, 0, 0, 0, 0])   # x^6+x+1
_IRRED[125] = (5, [1, 1, 0])           # x^3+x+1 (no roots mod 5)
_IRRED[343] = (7, [2, 0, 0])           # x^3+2: 2 is a non-cube mod 7 -> irreducible

def test(q):
    Q = q ** 3
    F = GF(Q)
    pw = lambda a, e: _pow(F, a, e)
    sub = [a for a in range(Q) if pw(a, q) == a]           # subfield F_q
    assert len(sub) == q
    units = [a for a in range(1, Q)]
    N = q * q + q + 1
    norm = {a: pw(a, N) for a in units}                    # x^{q^2+q+1} in F_q^*
    cosets = {}
    for a in units:
        cosets.setdefault(norm[a], set()).add(a)
    # directions: representatives of F^*/F_q^*
    seen = set(); dirs = []
    for a in units:
        if a in seen: continue
        orb = {F.mul[a][t] for t in sub if t}
        seen |= orb; dirs.append(a)
    assert len(dirs) == N
    def is_kakeya(S):
        S = set(S)
        for d in dirs:
            if not any(all(F.add[x][F.mul[t][d]] in S for t in sub) for x in S):
                return False
        return True
    results = {}
    for c, S in sorted(cosets.items()):
        results[f"norm={c} (size {len(S)})"] = is_kakeya(S)
        results[f"norm={c}+0 (size {len(S)+1})"] = is_kakeya(S | {0})
    return results

def _pow(F, a, e):
    r = 1
    while e:
        if e & 1: r = F.mul[r][a]
        a = F.mul[a][a]; e >>= 1
    return r

if __name__ == "__main__":
    for q in map(int, sys.argv[1:]):
        print(q, test(q))
