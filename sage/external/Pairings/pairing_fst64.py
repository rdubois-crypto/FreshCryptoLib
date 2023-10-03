from pairing import *

def miller_loop_opt_ate_fst64(Q, P, a, u):
    """Miller loop (f_{u,Q}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fq4) of order r
    - `P`: point on E(Fp)
    - `a`: curve coefficient in short Weierstrass y^2 = x^3 + a*x
    - `u`: seed for curve parameters

    E has type FST 6.4, 4 | k, D=1, t=x+1
    """
    m, uQ = miller_function_ate(Q, P, a, u, m0=1)
    return m

def miller_loop_opt_ate_fst64_cln_b0(Q, P, a, u):
    """Miller loop f_{u,Q}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fq4) of order r
    - `P`: point on E(Fp)
    - `a`: curve coefficient in short Weierstrass y^2 = x^3 + a*x
    - `u`: seed for curve parameters

    E has type FST 6.4, 4 | k, D=1, t=x+1
    """
    m, uQ = miller_function_ate_cln_b0(Q, P, a, u)
    return m

def final_exp_hard_fst64_k20(m, u):
    """
    m^e where
    e = ((u-1)//2)^2 (u^2+1) (q + u) ((u^2 - 1 + q^2) (q^4 + (u^4 - u^2 + 1)) + (u^4 - u^2)) + 1

    total cost 2*exp((u-1)//2) + 9*exp(|u|) + 9 M + frob + frob(2) + frob(4) + 3 conj + 1 conj if u < 0
    """
    assert (u % 2) == 1
    u1 = abs((u-1)//2)
    u_ = abs(u)
    m1 = m**u1
    m1 = m1**u1                                         # m^((u-1)/2)^2
    m1 = (m1**u_)**u_ * m1                              # m1^(u^2 + 1)
    m2 = m1**u_
    if u < 0:
        m2 = m2.conjugate()                              # m1^u
    m1 = m1.frobenius() * m2                             # m1^(q + u)
    m2 = (m1**u_)**u_ * m1.conjugate()                   # m1^(u^2 - 1)
    m2  = (m2**u_)**u_                                   # m1^(u^4 - u^2)
    m1 = m1.frobenius(4) * m2 * m1                       # m1^(q^4 + (u^4 - u^2 + 1))
    m1 = (m1**u_)**u_ * m1.conjugate() * m1.frobenius(2) # m1^(u^2 - 1 + q^2)
    m1 = m1 * m2
    m1 = m1 * m
    return m1

def final_exp_hard_fst64_k20_v1(m, u):
    """
    ((u-1)//2)^2 * (u^2+1) * (q + u) * ((q^2 + u^2 - 1)*(q^4 + u^4 + 1) -u^2*q^2) + 1;

    total cost 2*exp((u-1)//2) + 9*exp(|u|) + 8 M + frob + 2*frob(2) + frob(4) + 2 conj + 1 conj if u < 0
    """
    assert (u % 2) == 1
    u1 = abs((u-1)//2)
    u_ = abs(u)
    m1 = m**u1
    m1 = m1**u1                                         # m^((u-1)/2)^2
    m1 = (m1**u_)**u_ * m1                              # m1^(u^2 + 1)
    m2 = m1**u_
    if u < 0:
        m2 = m2.conjugate()                              # m1^u
    m1 = m1.frobenius() * m2                             # m1^(q + u)
    m3 = (m1**u_)**u_                                    # m1^(u^2)
    m1 = m1.frobenius(2) * m3 * m1.conjugate()           # m1^(q^2 + u^2 - 1)
    m1 = m1.frobenius(4) * (((m1**u_)**u_)**u_)**u_ * m1 # m1^(q^4 + u^4 + 1)
    m1 = m1 * m3.frobenius(2).conjugate()                # m1^((q^2 + u^2 - 1)*(q^4 + u^4 + 1) -u^2*q^2)
    m1 = m1 * m
    return m1

def final_exp_hard_fst64_k28(m, u):
    """
    Formula from
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya

    The generic formula with Th.1 is not very nice, though we obtain
    k = 28
    Phi_k = cyclotomic_polynomial(k)
    tx = x+1
    px = (x^16 - 2*x^15 + x^14 + x^2 + 2*x + 1)/4
    rx = Phi_k
    Tx = tx-1
    h1x = (px+1-tx) // rx # = (x - 1)^2/4 * (x^2 + 1)
    h2x = Phi_k(tx-1) // rx # = 1
    ex = Phi_k(px)//rx
    C = Phi_k.list()
    d = Phi_k.degree() # d == len(C)-1
    L = [0]*d
    L[d-1] = C[d]
    for i in range(d-2, -1, -1):
        print(i)
        L[i] = Tx*L[i+1] + C[i+1]
    L
    l0 = x^11 - x^9 + x^7 - x^5 + x^3 - x
    l1 = x^10 - x^8 + x^6 - x^4 + x^2 - 1
    l2 = x^9 - x^7 + x^5 - x^3 + x
    l3 = x^8 - x^6 + x^4 - x^2 + 1
    l4 = x^7 - x^5 + x^3 - x
    l5 = x^6 - x^4 + x^2 - 1
    l6 = x^5 - x^3 + x
    l7 = x^4 - x^2 + 1
    l8 = x^3 - x
    l9 = x^2 - 1
    l10 = x
    l11 = 1
    QQqu.<q,u> = QQ[]
    e1_uq = sum([QQx(L[i])(u)*q^i for i in range(len(L))])
    e1_uq.factor()
    (e1_uq // (q+u))
    e1_uq == (q+u)*(q^10 + (u^2-1)*q^8 + (u^4-u^2+1)*q^6 + (u^6-u^4+u^2-1)*q^4 + (u^8-u^6+u^4-u^2+1)*q^2 + (u^10-u^8+u^6-u^4+u^2-1))
    assert h1x * e1_uq([px, x]) + h2x == ex
    i8 = u^2-1
    i6 = u^4-u^2+1;              i6 == i8*u^2 + 1
    i4 = u^6-u^4+u^2-1;          i4 == i6*u^2 - 1
    i2 = u^8-u^6+u^4-u^2+1;      i2 == i4*u^2 + 1
    i0 = u^10-u^8+u^6-u^4+u^2-1; i0 == i2*u^2 - 1
    Phi_28(q)/r = (u-1)^2/4 * (u^2+1) * (q+u) * (q^10 + i8*q^8 + i6*q^6 + i4*q^4 + i2*q^2 + i0) + 1

    total cost 2*exp((u-1)//2) + 13*exp(|u|) + 13 M + frob + 5*frob(2) + conj (+ 1 conj if u < 0)
    """
    u_ = abs(u)
    u1 = abs(u-1)//2
    m1 = (m**u1)**u1
    m1 = (m1**u_)**u_ * m1
    mu = m1**u_
    if u < 0:
        mu = mu.conjugate()
    m1 = m1.frobenius() * mu           # m**((u-1)^2/4 * (u^2+1) * (q+u))
    # 2 exp(|u-1|/2) + 3 exp(|u|) + 2 M + f
    m1_inv = m1.conjugate()
    r = m1.frobenius(2)                # q10
    m2 = (m1**u_)**u_ * m1_inv         # m8 = m1**i8 = m1**u^2 * m1^(-1)
    r = (r * m2).frobenius(2)          # q8
    m2 = (m2**u_)**u_ * m1             # m6 = m1**i6 = m8**u^2 * m1
    r = (r * m2).frobenius(2)          # q6
    m2 = (m2**u_)**u_ * m1_inv         # m4 = m1**i4 = m6**u^2 * m1^(-1)
    r = (r * m2).frobenius(2)          # q4
    m2 = (m2**u_)**u_ * m1             # m2 = m1**i2 = m4**u^2 * m1
    r = (r * m2).frobenius(2)          # q2
    m2 = (m2**u_)**u_ * m1_inv         # m0 = m1**i0 = m2**u^2 * m1^(-1)
    r = (r * m2)                       # q0
    # 10 exp(|u|) + 10 M + 4f
    return r * m                       # M

def final_exp_hard_fst64_k28_v1(m, u):
    """
    Formula from
    https://eprint.iacr.org/2020/875
    Efficient Final Exponentiation via Cyclotomic Structure for
    Pairings over Families of Elliptic Curves
    Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya

    e = (u-1)^2/4 (u^2+1) (q+u) e1_uq + 1
    e1_uq = (u^2(u^2+q^2) + q^4)(u^2(u^2(u^2-1) + 1) + q^6-1) + (1-q^6)((u^2-1) + q^2)

    total cost 2*exp((u-1)//2) + 13*exp(|u|) + 12 M + frob + 2*frob(2)+frob(4)+2*frob(6) + 2conj (+ 1 conj if u < 0)
    """
    u_ = abs(u)
    u1 = abs(u-1)//2
    m1 = (m**u1)**u1
    m1 = (m1**u_)**u_ * m1
    mu = m1**u_
    if u < 0:
        mu = mu.conjugate()
    m1 = m1.frobenius() * mu           # m**((u-1)^2/4 * (u^2+1) * (q+u))
    # 2 exp(|u-1|/2) + 3 exp(|u|) + 2 M + f
    # second part: e1_uq
    m1_inv = m1.conjugate()
    m2 = (m1**u_)**u_ * m1_inv         # ^(u^2-1)
    m3 = m2 * m1.frobenius(2)          # ^((u^2-1) + q^2)
    m3 = m3.frobenius(6).conjugate() * m3
                                                 # m2 is m1^(u^2-1)
    m2 = (m2**u_)**u_ * m1                       #       u^2(u^2-1) + 1
    m2 = (m2**u_)**u_ * m1.frobenius(6) * m1_inv # ^(u^2(u^2(u^2-1) + 1) + q^6-1) done
    m4 = (m2**u_)**u_ * m2.frobenius(2)          # (u^2+q^2)
    m4 = (m2**u_)**u_ * m2.frobenius(4)
    m2 = m4 * m3
    # 10 exp(|u|) + 9 M + 5 f
    return m2 * m
