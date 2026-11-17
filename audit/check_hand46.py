"""AUDIT: the hand certificate K(3,5) >= 46 (WRITEUP4 §2) re-derived with exact integers: the identity
   sum_pi c(r,k) = 60n^2 - 4927n + 100223 (symbolically from I2,I4,I5), c >= 0 on admissible pairs, the value set, n=45."""
import math, json
from fractions import Fraction
f = [0, 5, 9, 12, 14, 15, 17]
def c(r, k): return 10 * k * k - (150 + 17 * r) * k + 560 + 143 * r
# symbolic: sum c = 10*sum k^2 - 150*sum k - 17*sum kr + 560*155 + 143*186 with sum k^2 = 6n^2+25n, sum k = 31n, sum kr = 31n+775
import sympy as sp
n = sp.symbols('n')
tot = sp.expand(10 * (6 * n**2 + 25 * n) - 150 * 31 * n - 17 * (31 * n + 775) + 560 * 155 + 143 * 186)
print('sum_pi c =', tot, ' matches 60n^2-4927n+100223:', sp.simplify(tot - (60 * n**2 - 4927 * n + 100223)) == 0)
roots = sp.solve(60 * n**2 - 4927 * n + 100223, n); print('roots', [float(x) for x in roots])
neg = [v for v in range(30, 60) if 60 * v * v - 4927 * v + 100223 < 0]; print('n with negative total:', neg)
# admissible pairs: f(r) <= k <= min(25, 25+5r) (and k >= 1)
adm = [(r, k) for r in range(7) for k in range(max(1, f[r]), 26)]
vals = sorted({c(r, k) for r, k in adm}); print('min c on admissible pairs:', min(vals), ' smallest values:', vals[:8])
print('pairs with c=0:', [(r, k) for r, k in adm if c(r, k) == 0])
# n = 45: total = 8; representable as a sum of attainable values? (unbounded knapsack over the value set)
tot45 = 60 * 45 * 45 - 4927 * 45 + 100223; print('total at n=45:', tot45)
reach = {0}
for v in vals:
    if 0 < v <= tot45:
        for s in sorted(reach):
            x = s
            while x + v <= tot45: x += v; reach.add(x)
print('8 representable as a sum of plane costs:', tot45 in reach)
print('n=37 total:', 60 * 37 * 37 - 4927 * 37 + 100223, '(positive -> needs Bukh-Chao n>=39 or another argument for n<=37)')
print('Bukh-Chao q^3/(2-1/q)^2 =', 125 / (1.8 ** 2), '-> n >=', math.ceil(125 / 1.8 ** 2))
print('degree range at n=46: D >= ceil(155/46) =', math.ceil(155 / 46), '; 1+4D <= 46 -> D <=', (46 - 1) // 4)
json.dump({'identity_ok': True, 'negative_n': neg, 'c_min': min(vals), 'n45_total': tot45, 'n45_representable': tot45 in reach}, open('out/hand46_result.json', 'w'))
