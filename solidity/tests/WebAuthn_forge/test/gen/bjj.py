"""BabyJubJub (circomlib): twisted Edwards <-> Montgomery <-> Weierstrass a=-3, Solidity constants.

Edwards   : 168700 x^2 + y^2 = 1 + 168696 x^2 y^2 over F_r (r = BN254 group order)
Montgomery: v^2 = u^3 + 168698 u^2 + u              u = (1+y)/(1-y), v = u/x
Weierstr. : Y^2 = X^3 + a_w X + b_w                  X = u + A/3, Y = v
Iso a=-3  : y'^2 = x'^3 - 3 x' + b'                 x' = mu^2 X, y' = mu^3 Y, mu^4 = -3/a_w
"""
r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
l = 2736030358979909402780800718157159386076813972158567259200215660948447373041
AE, DE = 168700, 168696
AM = 168698
B8 = (5299619240641551281634865583518297030282874472190772894086521144482721001553,
      16950150798460657717958625567821834550301663161624707787222815936182638968203)

inv = lambda x: pow(x, -1, r)
a_w = (3 - AM * AM) * inv(3) % r
b_w = (2 * AM ** 3 - 9 * AM) * inv(27) % r


def sqrt(n):
    """Tonelli-Shanks, None for a non residue."""
    n %= r
    if n == 0:
        return 0
    if pow(n, (r - 1) // 2, r) != 1:
        return None
    q, s = r - 1, 0
    while q % 2 == 0:
        q, s = q // 2, s + 1
    z = 5
    while pow(z, (r - 1) // 2, r) != r - 1:
        z += 1
    m, c, t, R = s, pow(z, q, r), pow(n, q, r), pow(n, (q + 1) // 2, r)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % r, i + 1
        b = pow(c, 1 << (m - i - 1), r)
        m, c, t, R = i, b * b % r, t * b * b % r, R * b % r
    return R


# mu^4 = -3/a_w: square root of a square root which is itself a square
t = (-3) * inv(a_w) % r
s1 = sqrt(t)
if pow(s1, (r - 1) // 2, r) != 1:
    s1 = r - s1
mu = sqrt(s1)
assert pow(mu, 4, r) == t
MU2, MU3 = mu * mu % r, pow(mu, 3, r)
b_iso = b_w * pow(mu, 6, r) % r
A3MU2 = AM * inv(3) % r * MU2 % r   # mu^2 * A/3


def ed_on(P):
    x, y = P
    return (AE * x * x + y * y - 1 - DE * x * x * y * y) % r == 0


def ed_add(P, Q):
    (x1, y1), (x2, y2) = P, Q
    k = DE * x1 * x2 * y1 * y2 % r
    return ((x1 * y2 + y1 * x2) * inv(1 + k) % r, (y1 * y2 - AE * x1 * x2) * inv(1 - k) % r)


def ed_mul(k, P):
    R = (0, 1)
    while k:
        if k & 1:
            R = ed_add(R, P)
        P = ed_add(P, P)
        k >>= 1
    return R


def ed2w(P):
    """Edwards -> Weierstrass a=-3, None for the neutral point (0,1)."""
    x, y = P
    if (x, y) == (0, 1):
        return None
    u = (1 + y) * inv(1 - y) % r
    if x == 0:  # point d'ordre 2 (0,-1) -> (u,v) = (0,0)
        return (A3MU2, 0)
    v = u * inv(x) % r
    return ((MU2 * u + A3MU2) % r, MU3 * v % r)


def w2ed(W):
    if W is None:
        return (0, 1)
    X, Y = W
    u = (X - A3MU2) * inv(MU2) % r
    v = Y * inv(MU3) % r
    if v == 0:
        return (0, r - 1)
    return (u * inv(v) % r, (u - 1) * inv(u + 1) % r)


def w_on(W):
    X, Y = W
    return (Y * Y - X ** 3 + 3 * X - b_iso) % r == 0


def w_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % r == 0:
            return None
        lam = (3 * P[0] * P[0] - 3) * inv(2 * P[1]) % r
    else:
        lam = (Q[1] - P[1]) * inv(Q[0] - P[0]) % r
    x = (lam * lam - P[0] - Q[0]) % r
    return (x, (lam * (P[0] - x) - P[1]) % r)


def w_mul(k, P):
    R = None
    while k:
        if k & 1:
            R = w_add(R, P)
        P = w_add(P, P)
        k >>= 1
    return R


G = ed2w(B8)

if __name__ == "__main__":
    assert ed_on(B8) and w_on(G) and w2ed(G) == B8
    # group morphism
    for k in (2, 3, 5, 12345):
        assert ed2w(ed_mul(k, B8)) == w_mul(k, G)
    assert w_mul(l, G) is None
    print("a_w   ", hex(a_w))
    print("b_iso ", hex(b_iso))
    print("MU2   ", hex(MU2))
    print("MU3   ", hex(MU3))
    print("A3MU2 ", hex(A3MU2))
    print("C3H   ", hex(3 * inv(2) % r))
    print("gx    ", hex(G[0]))
    print("gy    ", hex(G[1]))
    for k in (2, 3):
        P = w_mul(k, G)
        print(f"g{k}x   ", hex(P[0]))
        print(f"g{k}y   ", hex(P[1]))
    i2, i3 = pow(2, -1, l), pow(3, -1, l)
    for c, nm in zip([1, 2, 3, i2, 3 * i2 % l, i3, 2 * i3 % l], ["1", "2", "3", "1/2", "3/2", "1/3", "2/3"]):
        P = w_mul(c, G)
        print("small", nm, hex(c), hex(P[0]), hex(P[1]))
    print("2^254 > 4r ?", 4 * r < 2 ** 256, " 3r < 2^256 ?", 3 * r < 2 ** 256)
