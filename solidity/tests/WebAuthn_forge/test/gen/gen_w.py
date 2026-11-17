"""Weierstrass vectors for BabyJubJub: 8A -> W8, uG+vQ with degenerate cases and rare loop paths.
usage: python3 gen_w.py <keys.json> <weier.json>"""
import json, secrets, sys
from bjj import *

def H(v): return hex(v)
keys = json.load(open(sys.argv[1]))
w8 = {"ax": [], "ay": [], "wx": [], "wy": []}
for x, y in keys:
    A = (int(x), int(y))
    W = ed2w(ed_mul(8, A))
    for k, v in zip(("ax", "ay", "wx", "wy"), (A[0], A[1], W[0], W[1])):
        w8[k].append(H(v))
rnd = lambda: secrets.randbelow(l - 1) + 1
i2, i3 = pow(2, -1, l), pow(3, -1, l)
cases = [(rnd(), rnd(), rnd()) for _ in range(24)]
for c in [1, 2, 3, i2, 3 * i2 % l, i3, 2 * i3 % l, l - 1, l - 2, l - 3, l - i2, l - i3, 4, 5, 6]:
    for (u, v) in [(1, 1), (3, 1), (1, 3), (2, 2), (5, 7), (l - 1, 1), (1, l - 1), (rnd(), rnd()), (rnd(), 0), (0, rnd()), (12, 4)]:
        cases.append((c, u, v))
for _ in range(6):
    q, u = rnd(), rnd()
    cases.append((q, u, (-u) * pow(q, -1, l) % l))
import rare_cases
cases += rare_cases.extra(l)
mma = {k: [] for k in ("qx", "qy", "u", "v", "x", "y")}
for q, u, v in cases:
    Q = w_mul(q, G)
    R = w_add(w_mul(u % l, G), w_mul(v % l, Q)) if (u % l or v % l) else None
    R = R or (0, 0)
    for k, val in zip(("qx", "qy", "u", "v", "x", "y"), (Q[0], Q[1], u, v, R[0], R[1])):
        mma[k].append(H(val))
json.dump({"w8": w8, "mma": mma}, open(sys.argv[2], 'w'))
print(len(keys), "keys,", len(cases), "mulmuladd cases")
