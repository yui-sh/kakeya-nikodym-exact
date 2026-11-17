"""Known constructions / bounds as functions of q, n. Also builds Dvir's set explicitly and checks it."""
import itertools, json, sys
from math import comb
from gf import GF
from check import check

def plane_min(q):      # Blokhuis-Mazzocca
    return q * (q + 1) // 2 + ((q - 1) // 2 if q % 2 else 0)

def dvir_set(n, q):
    F = GF(q)
    sq = {F.mul[a][a] for a in range(q)}
    S = set()
    if q % 2:
        for v in itertools.product(range(q), repeat=n):
            *al, b = v
            b2 = F.mul[b][b]
            if all(F.add[a][b2] in sq for a in al):
                S.add(v)
        for v in itertools.product(range(q), repeat=n - 1):
            S.add(v + (0,))
    else:
        for v in itertools.product(range(q), repeat=n):
            *al, b = v
            if all(any(F.add[F.mul[g][g]][F.mul[g][b]] == a for g in range(q)) for a in al):
                S.add(v)
    return S

def dvir_bound(n, q):  # Dvir: |K| >= C(q+n-1, n)
    return comb(q + n - 1, n)

def dkss_bound(n, q):  # Dvir-Kopparty-Saraf-Sudan: |K| >= q^n / 2^n  (asymptotic form; exact: (q/2)^n ... keep ceil)
    return -(-q**n // 2**n)

def union_of_planes_construction(n, q):
    """Trivial n=3 construction: not used; placeholder for comparison (q^2 * ... )."""
    return None

if __name__ == "__main__":
    rows = []
    for n, qs in [(2, [2, 3, 4, 5, 7, 8, 9]), (3, [2, 3, 4, 5, 7])]:
        for q in qs:
            S = dvir_set(n, q)
            miss = check(n, q, S)
            assert not miss, (n, q, miss)
            rows.append({"n": n, "q": q, "dvir_size": len(S), "dvir_lb": dvir_bound(n, q),
                         "plane_min": plane_min(q) if n == 2 else None})
            print(rows[-1])
    json.dump(rows, open("out/constructions.json", "w"), indent=1)
