"""Finite field GF(q) for q in {2,3,4,5,7,8,9} (fresh implementation for the release).

Elements are the integers 0..q-1.  For a prime power q = p^k the integer encodes the
coefficient vector (c_0, ..., c_{k-1}) of a polynomial c_0 + c_1 X + ... over GF(p) in base p
(c_0 = least significant digit).  Multiplication reduces modulo the irreducible polynomial
listed in IRRED.  This is the same convention as the research-session code, so the saved
witness files are read with the same meaning; the tables are recomputed here from scratch
and every field axiom is checked exhaustively in __init__.
"""
IRRED = {  # q: (p, monic irreducible polynomial X^k + a_{k-1} X^{k-1} + ... + a_0 given as [a_0, ..., a_{k-1}])
    4: (2, [1, 1]),       # X^2 + X + 1
    8: (2, [1, 1, 0]),    # X^3 + X + 1
    9: (3, [1, 0]),       # X^2 + 1
}


def _prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            k, m = 0, q
            while m % p == 0:
                m //= p
                k += 1
            if m != 1:
                raise ValueError(f"{q} is not a prime power")
            return p, k
    raise ValueError(q)


class GF:
    def __init__(self, q):
        self.q = q
        self.p, self.k = p, k = _prime_power(q)
        if k == 1:
            self.add = [[(a + b) % p for b in range(q)] for a in range(q)]
            self.mul = [[(a * b) % p for b in range(q)] for a in range(q)]
        else:
            poly = IRRED[q][1]
            assert len(poly) == k

            def digits(a):
                return [(a // p ** i) % p for i in range(k)]

            def value(v):
                return sum(c * p ** i for i, c in enumerate(v))

            def mul(a, b):
                u, v = digits(a), digits(b)
                prod = [0] * (2 * k - 1)
                for i in range(k):
                    for j in range(k):
                        prod[i + j] = (prod[i + j] + u[i] * v[j]) % p
                for d in range(2 * k - 2, k - 1, -1):  # X^d = X^{d-k} * X^k, X^k = -(a_0 + ... )
                    c = prod[d]
                    if c:
                        prod[d] = 0
                        for i in range(k):
                            prod[d - k + i] = (prod[d - k + i] - c * poly[i]) % p
                return value(prod[:k])

            self.add = [[value([(x + y) % p for x, y in zip(digits(a), digits(b))]) for b in range(q)] for a in range(q)]
            self.mul = [[mul(a, b) for b in range(q)] for a in range(q)]
        self.neg = [next(b for b in range(q) if self.add[a][b] == 0) for a in range(q)]
        self.inv = [None] + [next(b for b in range(1, q) if self.mul[a][b] == 1) for a in range(1, q)]
        self.sub = [[self.add[a][self.neg[b]] for b in range(q)] for a in range(q)]
        self.frob = [self.power(a, p) for a in range(q)]  # x -> x^p, a field automorphism
        self._check_axioms()

    def power(self, a, e):
        r = 1
        for _ in range(e):
            r = self.mul[r][a]
        return r

    def _check_axioms(self):
        q, add, mul = self.q, self.add, self.mul
        for a in range(q):
            assert add[a][0] == a and mul[a][1] == a and mul[a][0] == 0
            for b in range(q):
                assert add[a][b] == add[b][a] and mul[a][b] == mul[b][a]
                if a and b:
                    assert mul[a][b] != 0  # no zero divisors
                for c in range(q):
                    assert add[add[a][b]][c] == add[a][add[b][c]]
                    assert mul[mul[a][b]][c] == mul[a][mul[b][c]]
                    assert mul[a][add[b][c]] == add[mul[a][b]][mul[a][c]]
        assert all(mul[a][self.inv[a]] == 1 for a in range(1, q))
        # Frobenius is additive and multiplicative
        for a in range(q):
            for b in range(q):
                assert self.frob[add[a][b]] == add[self.frob[a]][self.frob[b]]
                assert self.frob[mul[a][b]] == mul[self.frob[a]][self.frob[b]]


if __name__ == "__main__":
    for q in (2, 3, 4, 5, 7, 8, 9):
        F = GF(q)
        print(f"GF({q}): p={F.p} k={F.k} axioms OK")
