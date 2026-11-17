"""(c, u, v) driving the joint window loop through each rare path (Q = c.G, c outside the small multiples)."""
def events(n, c, u, v):
    ev = []
    if u == 0 and v == 0:
        return ev
    i = 252
    while not ((u >> i) | (v >> i)):
        i -= 4
    b, first = 0, True
    while i >= 0:
        for sh in (2, 0):
            du, dv = (u >> (i + sh)) & 3, (v >> (i + sh)) & 3
            b = 4 * b % n
            if du or dv:
                t = (du + dv * c) % n
                if b == 0 and not first:
                    ev.append("O_mid")
                elif b != 0 and b == t:
                    ev.append("dbl")
                elif b != 0 and (b + t) % n == 0:
                    ev.append("inf")
                first = False
                b = (b + t) % n
        i -= 4
    return ev

def find(n, cs=(4, 5, 6, 7, 9, 10, 11, 12)):
    want = {"dbl": [], "inf": [], "O_mid": []}
    for c in cs:
        for u in range(0, 256):
            for v in range(0, 64):
                for e in set(events(n, c, u, v)):
                    if len(want[e]) < 4:
                        want[e].append((c, u, v))
        if all(len(x) >= 4 for x in want.values()):
            break
    return want


def balanced(n, c):
    """positive (U, V) ~ sqrt(n) with U + c.V = 0 mod n (truncated extended Euclid)."""
    r0, r1, t0, t1 = n, c, 0, 1
    while r1 * r1 >= n:
        q = r0 // r1
        r0, r1, t0, t1 = r1, r0 - q * r1, t1, t0 - q * t1
    # r1 = t1.c mod n -> r1 - t1.c = 0 ; on veut V = -t1 > 0
    U, V = r1, -t1
    if V < 0:
        U, V = r0, -t0
    assert (U + c * V) % n == 0 and U > 0 and V > 0 and U < 1 << 130 and V < 1 << 130
    return U, V


def extra(n, cs=(0x5f1d6b2a9c3e8d4f7a0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6f708192a3b4c5,
                  0x3a7c9e1b5d2f4a6c8e0b1d3f5a7c9e1b2d4f6a8c0e2b4d6f8a1c3e5b7d9f0a2c)):
    """rare cases: dbl (small scalars), inf at the end, O in the middle of the loop."""
    out = []
    for c in cs:
        c %= n
        U, V = balanced(n, c)
        out += [(c, U, V)]                                    # result O, last addition -> inf
        out += [(c, (U << 124) | 0x9, (V << 124) | 0x6)]       # O mid loop, then more digits
        out += [(c, (U << 124), (V << 124) | 0x3)]
    out += [(4, 4, 1), (4, 16, 4), (4, 12, 3), (5, 20, 4), (7, 28, 4)]
    return out
