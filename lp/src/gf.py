"""Finite field GF(q) arithmetic for prime powers q = p^k via lookup tables.
Elements are integers 0..q-1; for prime powers an integer encodes the
coefficient vector of a polynomial over GF(p) in base p."""
import itertools

# irreducible monic polynomials over GF(p), as coefficient lists low->high (excluding leading 1)
_IRRED = {
    4: (2, [1, 1]),        # x^2 + x + 1 over GF(2)
    8: (2, [1, 1, 0]),     # x^3 + x + 1
    9: (3, [1, 0]),        # x^2 + 1 over GF(3) (irreducible since -1 is a non-square mod 3)
    16: (2, [1, 1, 0, 0]), # x^4 + x + 1
    25: (5, [2, 0]),       # x^2 + 2 over GF(5): 2 is a non-square mod 5
    27: (3, [1, 2, 0]),    # x^3 + 2x + 1 over GF(3)
}

def _factor_prime_power(q):
    for p in range(2, q + 1):
        if q % p == 0:
            k = 0
            m = q
            while m % p == 0:
                m //= p; k += 1
            assert m == 1, f"{q} is not a prime power"
            return p, k
    raise ValueError

class GF:
    def __init__(self, q):
        self.q = q
        self.p, self.k = _factor_prime_power(q)
        p, k = self.p, self.k
        if k == 1:
            self.add = [[(a + b) % p for b in range(q)] for a in range(q)]
            self.mul = [[(a * b) % p for b in range(q)] for a in range(q)]
        else:
            _, poly = _IRRED[q]
            assert len(poly) == k
            def to_vec(a):
                return [(a // p**i) % p for i in range(k)]
            def to_int(v):
                return sum(c * p**i for i, c in enumerate(v))
            def mulv(u, v):
                prod = [0] * (2 * k - 1)
                for i, a in enumerate(u):
                    for j, b in enumerate(v):
                        prod[i + j] = (prod[i + j] + a * b) % p
                # reduce: x^k = -(poly)
                for d in range(2 * k - 2, k - 1, -1):
                    c = prod[d]
                    if c:
                        prod[d] = 0
                        for i in range(k):
                            prod[d - k + i] = (prod[d - k + i] - c * poly[i]) % p
                return prod[:k]
            self.add = [[to_int([(x + y) % p for x, y in zip(to_vec(a), to_vec(b))]) for b in range(q)] for a in range(q)]
            self.mul = [[to_int(mulv(to_vec(a), to_vec(b))) for b in range(q)] for a in range(q)]
        self.neg = [next(b for b in range(q) if self.add[a][b] == 0) for a in range(q)]
        self.inv = [None] + [next(b for b in range(1, q) if self.mul[a][b] == 1) for a in range(1, q)]
        self._check()

    def _check(self):
        q = self.q
        # field axioms spot-check: distributivity + no zero divisors + inverse
        for a in range(q):
            for b in range(q):
                assert self.mul[a][b] == self.mul[b][a]
                assert self.add[a][b] == self.add[b][a]
                for c in range(q):
                    assert self.mul[a][self.add[b][c]] == self.add[self.mul[a][b]][self.mul[a][c]]
                    assert self.mul[self.mul[a][b]][c] == self.mul[a][self.mul[b][c]]
                if a and b:
                    assert self.mul[a][b] != 0
        assert all(self.mul[a][self.inv[a]] == 1 for a in range(1, q))

if __name__ == "__main__":
    for q in [2, 3, 4, 5, 7, 8, 9, 11, 16, 25, 27]:
        F = GF(q)
        print(q, "ok; p=%d k=%d" % (F.p, F.k))
