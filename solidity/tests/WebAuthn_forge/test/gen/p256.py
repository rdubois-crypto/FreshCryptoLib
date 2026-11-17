"""Reference P-256 arithmetic (affine, None = point at infinity) and vector generation.
usage: python3 p256.py <p256.json>"""
import hashlib
import json
import secrets
import sys

p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
a = p - 3
b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
G = (0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
     0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5)


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0:
            return None
        lam = (3 * P[0] * P[0] + a) * pow(2 * P[1], -1, p) % p
    else:
        lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x = (lam * lam - P[0] - Q[0]) % p
    return (x, (lam * (P[0] - x) - P[1]) % p)


def neg(P):
    return None if P is None else (P[0], (-P[1]) % p)


def mul(k, P):
    R = None
    k %= n
    while k:
        if k & 1:
            R = add(R, P)
        P = add(P, P)
        k >>= 1
    return R


def mulmuladd(Q, u, v):
    return add(mul(u, G), mul(v % n, Q))


def aff(P):
    return (0, 0) if P is None else P


def sign(d, h, k):
    R = mul(k, G)
    r = R[0] % n
    s = pow(k, -1, n) * (h + r * d) % n
    return r, s


def rnd():
    return secrets.randbelow(n - 1) + 1


def gen(path):
    out = {"ecdsa": [], "mma": []}
    # random valid and invalid ECDSA signatures
    for i in range(48):
        d = rnd()
        Q = mul(d, G)
        h = int.from_bytes(hashlib.sha256(secrets.token_bytes(32)).digest(), "big")
        r, s = sign(d, h, rnd())
        ok = True
        if i % 4 == 3:  # bit flip
            h ^= 1 << secrets.randbelow(256)
            ok = False
        out["ecdsa"].append({"h": h, "r": r, "s": s, "qx": Q[0], "qy": Q[1], "ok": ok})
    # mulmuladd: random + edge cases (Q = c.G for small c, extreme scalars)
    cases = []
    for _ in range(24):
        q = rnd()
        cases.append((q, rnd(), rnd()))
    inv2, inv3 = pow(2, -1, n), pow(3, -1, n)
    special = [1, 2, 3, 4, n - 1, n - 2, n - 3, n - inv2, n - 3 * inv2 % n, n - inv3,
               n - 2 * inv3 % n, inv2, 3 * inv2 % n, inv3, 2 * inv3 % n, 5, n - 5, 6, 9]
    for c in special:
        for (u, v) in [(1, 1), (3, 1), (1, 3), (2, 2), (5, 7), (n - 1, 1), (1, n - 1),
                       (n - 1, n - 1), (rnd(), rnd()), (rnd(), 0), (0, rnd()),
                       (12, 4), (6, 6), (1 << 255 >> 1, 1 << 200)]:
            cases.append((c, u % n, v % n))
    for _ in range(8):  # point at infinity result
        q, u = rnd(), rnd()
        v = (-u) * pow(q, -1, n) % n
        cases.append((q, u, v))
    import rare_cases  # rare loop paths: acc = E, acc = -E, acc = O mid loop
    cases += rare_cases.extra(n)
    for (q, u, v) in cases:
        Q = mul(q, G)
        R = aff(mulmuladd(Q, u, v))
        out["mma"].append({"qx": Q[0], "qy": Q[1], "u": u, "v": v, "x": R[0], "y": R[1]})
    flat = {}
    for k in ("ecdsa", "mma"):
        for key in out[k][0]:
            flat[k + "_" + key] = [hex(int(e[key])) for e in out[k]]
    with open(path, "w") as f:
        json.dump(flat, f)
    print(len(out["ecdsa"]), "ecdsa,", len(out["mma"]), "mulmuladd")


if __name__ == "__main__":
    gen(sys.argv[1])
