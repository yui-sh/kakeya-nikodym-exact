"""Hand-checkable certificate for K(3,5) >= 46 (session IV, WRITEUP4 §2). For every plane pi of AG(3,5) with r = # witness
lines inside and k = |K ∩ pi| define c(r,k) = 10k^2 - (150+17r)k + 560 + 143r.
(1) c(r,k) >= 0 for every admissible (r,k): k any integer for r <= 2; k >= 12, 14, 15, 17 for r = 3, 4, 5, 6.
(2) sum_pi c(r_pi,k_pi) = 10*(6n^2+25n) - 150*31n - 17*(31n+775) + 560*155 + 143*186 = 60n^2 - 4927n + 100223,
    which is < 0 for 38 <= n <= 44 (roots 37.13, 44.98; n <= 38 is excluded separately by the integer model / Bukh-Chao) and = 8 for n = 45, while c takes no value in {1..5, 8} -> n >= 46."""
import sympy as sp
f = {0: 0, 1: 5, 2: 9, 3: 12, 4: 14, 5: 15, 6: 17}
c = lambda r, k: 10 * k * k - (150 + 17 * r) * k + 560 + 143 * r
vals = set()
for r in range(7):
    for k in range(max(f[r], 1), 26):
        assert c(r, k) >= 0, (r, k, c(r, k)); vals.add(c(r, k))
    if r <= 2:
        assert all(c(r, k) >= 0 for k in range(0, 200)), r
print("c(r,k) >= 0 on all admissible pairs; smallest values:", sorted(vals)[:8])
n = sp.symbols("n")
total = sp.expand(10 * (6 * n**2 + 25 * n) - 150 * 31 * n - 17 * (31 * n + 775) + 560 * 155 + 143 * 186)
print("sum_pi c =", total)
assert total == 60 * n**2 - 4927 * n + 100223
for v in range(39, 47): print(f"  n={v}: {total.subs(n, v)}")
assert all(total.subs(n, v) < 0 for v in range(38, 45)) and total.subs(n, 45) == 8 and 8 not in vals and all(v not in vals for v in range(1, 6))
print("=> every witness system in AG(3,5) has n >= 46: certificate OK")
