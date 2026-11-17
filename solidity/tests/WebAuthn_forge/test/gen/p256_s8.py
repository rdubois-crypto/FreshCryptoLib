"""8-dim precomputation table (FCL Precalc_8dim / sage format) and vectors for S8_extcode / S8_hackmem.
usage: python3 p256_s8.py <out.json> [2G]"""
import json, secrets, sys
from p256 import *

def table(Q):
    base = [mul(1 << (64 * j), G) for j in range(4)] + [mul(1 << (64 * j), Q) if Q else None for j in range(4)]
    out = b''
    for i in range(256):
        P = None
        for j in range(8):
            if i >> j & 1:
                P = add(P, base[j])
        x, y = aff(P)
        out += x.to_bytes(32, 'big') + y.to_bytes(32, 'big')
    return out

if len(sys.argv) > 2 and sys.argv[2] == "2G":
    # Q = 2G: (u,v) = (2,1) goes through acc = E (doubling), (n-2,1) through acc = -E (result O)
    Q = mul(2, G)
    cases = [(2, 1), (n - 2, 1), (n - 4, 2), (6, 1), (1, 0), (0, 1)] + [(rnd(), rnd()) for _ in range(4)]
else:
    d = rnd()
    Q = mul(d, G)
    cases = []
    for _ in range(16):
        u, v = rnd(), rnd()
        cases.append((u, v))
    cases += [(1, 0), (0, 1), (n - 1, 1), (1 << 255 >> 1, 3), (5, 7)]
tab = table(Q)
res = {"qx": hex(Q[0]), "qy": hex(Q[1]), "table": "0x" + tab.hex(),
       "u": [hex(u) for u, v in cases], "v": [hex(v) for u, v in cases],
       "x": [hex(aff(mulmuladd(Q, u, v))[0]) for u, v in cases]}
json.dump(res, open(sys.argv[1], 'w'))
print(len(tab), "table bytes,", len(cases), "cases")
