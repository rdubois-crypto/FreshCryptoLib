from pairing import *
from sage.misc.functional import cyclotomic_polynomial

"""
Optimal ate pairing on Brezing-Weng curves of Clarisse, Duquesne and Sanders
CANS'2020
eprint 2020/760 https://eprint.iacr.org/2020/760
"""

def add_line_j_denom(S, Q, P):
    """
    Compute S+Q and l_{S,Q}(P) the line through S and Q evaluated at P

    Jacobian coordinates (X,Y,Z,Z^2) correspond to the affine coordinates
    (x,y) = (X/Z^2,Y/Z^3)

    INPUT:
    - `S`: point in extended Jacobian coordinates (X, Y, Z, Z^2)
    - `Q`: point in affine coordinates (xQ, yQ)
    - `P`: point in affine coordinates (xP, yP)

    RETURN:
    (ln, ld) line ln/ld through S, Q evaluated at P, point S+Q=(X',Y',Z',Z'^2)
    numerator and denominator of the line without simplification
    """
    ln, SQ = add_line_j(S, Q, P)
    ld = SQ[2]
    return (ln, ld), SQ

def double_line_j_denom(S, P, a):
    """
    Compute 2*S and l_{S,S}(P) the line tangent at S evaluated at P
    Extended Jacobian coordinates (X,Y,Z,Z^2) correspond to the affine
    coordinates (x,y) = (X/Z^2,Y/Z^3)

    INPUT:
    - `S`: elliptic curve point in Jacobian coordinates (X,Y,Z,Z^2)
    - `P`: elliptic curve point in affine coordinates (x,y)
    - `a`: curve coefficient y^2=x^3+a*x+b

    RETURN:
    (ln, ld) line ln/ld tangent at S evaluated at P, point 2*S = (X',Y',Z',Z'^2)
    numerator and denominator of the line without simplification
    """
    ln, S2 = double_line_j(S,P,a)
    ld = S2[2]*S[3]
    return (ln, ld), S2

def vertical_j_denom(S, P):
    """
    The vertical v_{S}(P) at S (in extended Jacobian coordinates) evaluated at P
    Extended Jacobian coordinates (X,Y,Z,Z^2) correspond to the affine
    coordinates (x,y) = (X/Z^2,Y/Z^3)

    INPUT:
    - `S`: elliptic curve point in Jacobian coordinates (X,Y,Z,Z^2)
    - `P`: elliptic curve point in affine coordinates (x,y)

    RETURN:
    (vn, vd) vertical vn/vd vertical at S evaluated at P
    numerator and denominator of the vertical without simplification
    """
    Z2 = S[3]
    X = S[0]
    xP = P[0]
    vd = Z2
    vn = Z2 * xP - X
    return (vn, vd)

def miller_function_tate_with_denom(P, Q, a, T, m0=1, m0_inv=1):
    """
    Miller function f_{T, P}(Q) with denominators and verticals
    Simplify only the factors in Fp

    INPUT:
    - `P`: r-torsion point on E(Fp) in affine coordinates
    - `Q`: r-torsion point on E(Fpk) of trace 0 in affine coordinates
           that is, Q = [k]Q0 - trace(Q0) for some Q0, so that pi(Q) = [p]Q
    - `T`: scalar, can be the order of the points
    - `a`: curve coefficient in y^2=x^3+a*x+b (short Weierstrass)
    - `m0`: optional parameter, for multi-exponentiation optimization,
            this is not needed for simple ate pairing
    - `m0_inv`: optional parameter, for multi-exponentiation optimization,
            this is not needed for simple tate pairing

    Problem of negative T:
    If T < 0, then f_{-|T|, P}(Q) is computed thanks to the formula
    f_{uv,P} = f_{u,P}^v*f_{v,[u]P} and with u=-1, v=|T|:
    f_{-|T|,P} = f_{-1,P}^|T|*f_{|T|,-P} and f_{-1,P} is a vectical line,
    then f_{-|T|,P} = f_{|T|,-P} * (v_P(Q))^|T|
    So for negative T, the miller loop is initialized at v_P(Q) instead of 1.

    The denominators are in Fp and can be removed.
    """
    xP, yP = P[0], P[1]
    xQ, yQ = Q[0], Q[1]
    PP = (xP, yP)
    QQ = (xQ, yQ)
    S = (xP, yP, 1, 1) # extended Jacobian coordinates
    m = m0 # for accumulating the numerators
    md = m0_inv # for accumulating the denominators
    with_m0 = m0 != 1
    with_m0_inv = m0_inv != 1
    negative_T = False
    if T < 0:
        negative_T = True
        T = -T
        PP = (xP, -yP)
        S = (xP, -yP, 1, 1)
        vP_Q = xQ - xP
        if with_m0_inv:
            md = md * vQ_P
        else:
            md = vP_Q
    loop = Integer(T).digits(2)
    # do not compute the final vertical if it is zero
    i = len(loop)-2
    while i >= 0 and S[2] != 0:
        bi = loop[i]
        ln, S = double_line_j(S, QQ, a)
        m = m**2 * ln
        # the denominator ld = S[2] is in Fp, can ignore it
        if S[2] != 0:
            vn, vd = vertical_j_denom(S, QQ)
            # vd is in Fp, can ignore it
            md = md**2 * vn
        if bi == 1 and S[2] != 0:
            ln, S = add_line_j(S, PP, QQ)
            if with_m0:
                m = m * m0
            if with_m0_inv:
                md = md * m0_inv
            if S[2] == 0: # S and PP are opposite, so actually the line is a vertical, do it manually
                m = m * (xQ - xP)
            else:
                m = m * ln
                # divide the line by the vertical at (S+P) evaluated at Q:
                # the new updated S is S+P with the former S
                vn, vd = vertical_j_denom(S, QQ)
                md = md * vn
            if negative_T:
                md = md * vP_Q
        i = i-1
    return m, md, S

def miller_loop_tate_with_denom(P, Q, a, r):
    m, md, S = miller_function_tate_with_denom(P, Q, a, r)
    return m/md

def miller_function_ate_with_denom(Q, P, a, T, m0=1, m0_inv=1):
    """
    Miller function f_{T, Q}(P) with denominators and verticals
    Simplify only the factors in Fp

    INPUT:
    - `Q`: r-torsion point on E(Fpk) of trace 0 in affine coordinates
           that is, Q = [k]Q0 - trace(Q0) for some Q0, so that pi(Q) = [p]Q
    - `P`: r-torsion point on E(Fp) in affine coordinates
    - `T`: scalar, T=(t-1) for ate pairing for example, can be negative
    - `a`: curve coefficient in y^2=x^3+a*x+b (short Weierstrass)
    - `m0`: optional parameter, for multi-exponentiation optimization,
            this is not needed for simple ate pairing
    - `m0_inv`: optional parameter, for multi-exponentiation optimization,
            this is not needed for simple ate pairing

    Problem of negative T:
    If T < 0, then f_{-|T|, Q}(P) is computed thanks to the formula
    f_{uv,Q} = f_{u,Q}^v*f_{v,[u]Q} and with u=-1, v=|T|:
    f_{-|T|,Q} = f_{-1,Q}^|T|*f_{|T|,-Q} and f_{-1,Q} is a vectical line,
    then f_{-|T|,Q} = f_{|T|,-Q} * (v_Q(P))^|T|
    So for negative T, the miller loop is initialized at v_Q(P) instead of 1.
    """
    xQ, yQ = Q[0], Q[1]
    xP, yP = P[0], P[1]
    QQ = (xQ, yQ)
    PP = (xP, yP)
    S = (xQ, yQ, 1, 1) # extended Jacobian coordinates
    m = m0 # for accumulating the numerators
    md = m0_inv # for accumulating the denominators
    with_m0 = m0 != 1
    with_m0_inv = m0_inv != 1
    negative_T = False
    if T < 0:
        negative_T = True
        T = -T
        QQ = (xQ, -yQ)
        S = (xQ, -yQ, 1, 1)
        vQ_P = xP - xQ
        if with_m0_inv:
            md = md * vQ_P
        else:
            md = vQ_P
    loop = Integer(T).digits(2)
    # do not compute the final vertical if it is zero
    i = len(loop)-2
    while i >= 0 and S[2] != 0:
        bi = loop[i]
        (ln, ld), S = double_line_j_denom(S, PP, a)
        m = m**2 * ln
        md = md**2 * ld
        if S[2] != 0:
            vn, vd = vertical_j_denom(S, PP)
            m = m * vd
            md = md * vn
        if bi == 1 and S[2] != 0:
            (ln, ld), S = add_line_j_denom(S, QQ, PP)
            if with_m0:
                m = m * m0
            if with_m0_inv:
                md = md * m0_inv
            if S[2] == 0: # S and QQ are opposite, so actually the line is a vertical, do it manually
                m = m * (xP - xQ)
            else:
                m = m * ln
                md = md * ld
                # divide the line by the vertical at (S+Q) evaluated at P:
                # the new updated S is S+Q with the former S
                vn, vd = vertical_j_denom(S, PP)
                m = m * vd
                md = md * vn
            if negative_T:
                md = md * vQ_P
        i = i-1
    return m, md, S

def miller_loop_ate_with_denom(Q, P, a, T):
    m, md, S = miller_function_ate_with_denom(Q, P, a, T)
    return m/md

def cyclotomic_inverse_k5(m):
    """
    Inverse in the cyclotomic subgroup of GF(p^5)
    m^(p^5-1) = 1
    cyclotomic subgroup means m = m0^(p-1) for some m0
    and then m^(1+p+p^2+p^3+p^4) = m0^(p^5-1) = 1
    finally, 1/m = m^(p+p^2+p^3+p^4)
    
    (p+p^2+p^3+p^4) = p * (p+1) * (p^2 + 1)

    cost: 3 Frobenius + 2 M
    """
    if m == 0 or m == 1:
        return m
    m = m.frobenius()
    m = m.frobenius() * m
    m = m.frobenius(2) * m
    return m

def cyclotomic_inverse_k7(m):
    """
    Inverse in the cyclotomic subgroup of GF(p^7)
    m^(p^7-1) = 1
    cyclotomic subgroup means m = m0^(p-1) for some m0
    and then m^(1+p+p^2+...+p^6) = m0^(p^7-1) = 1
    finally, 1/m = m^(p+p^2+...+p^6)
    
    (p+p^2+...+p^6) = p * (p + 1) * (p^4 + p^2 + 1)
                    = p * (p^3 + 1) * (p^2 + p + 1)

    cost: 4 Frobenius + 3 M
    """
    if m == 0 or m == 1:
        return m
    m = m.frobenius()
    m = m.frobenius() * m
    m = m.frobenius(4) * m.frobenius(2) * m
    return m

def cyclotomic_inverse_k11(m):
    """
    Inverse in the cyclotomic subgroup of GF(p^11)
    m^(p^11-1) = 1
    cyclotomic subgroup means m = m0^(p-1) for some m0
    and then m^(1+p+p^2+...+p^10) = m0^(p^11-1) = 1
    finally, 1/m = m^(p+p^2+...+p^10)
    
    (p+p^2+...+p^10) = p * (p + 1) * (p^8 + p^6 + p^4 + p^2 + 1)
                     = p * (p^5 + 1) * (p^4 + p^3 + p^2 + p + 1)

    cost: 6 Frobenius + 5 M
    """
    if m == 0 or m == 1:
        return m
    m = m.frobenius()
    m = m.frobenius(5) * m
    m = m.frobenius(4) * m.frobenius(3) * m.frobenius(2) * m.frobenius() * m
    return m

def cyclotomic_inverse_k13(m):
    """
    Inverse in the cyclotomic subgroup of GF(p^13)
    m^(p^13-1) = 1
    cyclotomic subgroup means m = m0^(p-1) for some m0
    and then m^(1+p+p^2+...+p^12) = m0^(p^13-1) = 1
    finally, 1/m = m^(p+p^2+...+p^12)
    
    (p+p^2+...+p^12) = p * (p + 1) * (p^2 + 1) * (p^8 + p^4 + 1)
                     = p * (p + 1) * (p^6 + 1) * (p^4 + p^2 + 1)
                     = p * (p^3 + 1) * (p^6 + 1) * (p^2 + p + 1)

    cost: 5 Frobenius + 4 M
    """
    if m == 0 or m == 1:
        return m
    m = m.frobenius()
    m = m.frobenius() * m
    m = m.frobenius(2) * m
    m = m.frobenius(8) * m.frobenius(4) * m
    return m
 
def cyclotomic_inverse_k17(m):
    """
    Inverse in the cyclotomic subgroup of GF(p^17)
    m^(p^17-1) = 1
    cyclotomic subgroup means m = m0^(p-1) for some m0
    and then m^(1+p+p^2+...+p^16) = m0^(p^17-1) = 1
    finally, 1/m = m^(p+p^2+...+p^16)
    
    (p+p^2+...+p^16) = p * (p + 1) * (p^2 + 1) * (p^4 + 1) * (p^8 + 1)

    cost: 5 Frobenius + 4 M
    """
    if m == 0 or m == 1:
        return m
    m = m.frobenius()
    m = m.frobenius() * m
    m = m.frobenius(2) * m
    m = m.frobenius(4) * m
    m = m.frobenius(8) * m
    return m

def cyclotomic_inverse_k19(m):
    """
    Inverse in the cyclotomic subgroup of GF(p^19)
    m^(p^19-1) = 1
    cyclotomic subgroup means m = m0^(p-1) for some m0
    and then m^(1+p+p^2+...+p^18) = m0^(p^19-1) = 1
    finally, 1/m = m^(p+p^2+...+p^18)
    
    (p+p^2+...+p^18) = p * (p+1) * (p^4 + p^2 + 1) * (p^12 + p^6 + 1)
                     = p * (p^3 + 1) * (p^2 + p + 1) * (p^12 + p^6 + 1)
                     = p * (p^9 + 1) * (p^2 + p + 1) * (p^6 + p^3 + 1)

    cost: 6 Frobenius + 5 M
    """
    if m == 0 or m == 1:
        return m
    m = m.frobenius()
    m = m.frobenius(9) * m
    m = m.frobenius(2) * m.frobenius() * m
    m = m.frobenius(6) * m.frobenius(3) * m
    return m

def cyclotomic_inverse(m, k):
    if not k in [5,7,11,13,17,19]:
        return 1/m
    if k == 5:
        return cyclotomic_inverse_k5(m)
    elif k == 7:
        return cyclotomic_inverse_k7(m)
    elif k == 11:
        return cyclotomic_inverse_k11(m)
    elif k == 13:
        return cyclotomic_inverse_k13(m)
    elif k == 17:
        return cyclotomic_inverse_k17(m)
    elif k == 19:
        return cyclotomic_inverse_k19(m)
    
####### optimal ate pairing

def miller_loop_opt_ate_bw13_bw19(Q, P, u, k):
    """
    Optimal ate Miller loop for BW13, BW19
    u^2 -u*p + p^2 = 0 mod r
    f_{u^2, Q}(P) * f_{-u, pQ}(P) * l_{[-u*p]Q,[p^2]Q}(P)
    where
    f_{-up, Q}(P) = f_{-u,Q}(P)^p f_{p, [-u]Q}(P)
                  = f_{p, Q}(P)^(-u) f_{-u, [p]Q}(P)
    after simplifying by the f_{p, iQ}(P) terms which are bilinear pairings,
    f_{-u, pQ}(P) = f_{-u,Q}(P)^p

    moreover,
    f_{u^2, Q}(P) = f_{u, Q}(P)^u * f_{u, [u]Q}(P)
                  = f_{-u, Q}(P)^(-u) * f_{-u, [-u]Q}(P)

    f_{-u, Q}(P) = f_{-1, Q}(P)^u f_{u, -Q}(P)
                 = f_{u, Q}(P)^(-1) f_{-1, uQ}(P)
                 = 1/(f_{u, Q}(P) * v_{uQ}(P))

    f_{-u, pQ}(P) = f_{-u,Q}(P)^p
                  = 1/(f_{u, Q}(P) * v_{uQ}(P))^p

    finally, if u > 0:
        f_{u^2, Q}(P) = f_{u, Q}(P)^u f_{u, uQ}(P)
    if u < 0:
        f_{u^2, Q}(P) = f_{-u, Q}(P)^(-u) f_{-u, -uQ}(P)
    """
    p = P.curve().base_field().characteristic()
    r = cyclotomic_polynomial(k)(u)
    ee = (p**k-1)//r
    if u < 0: # compute f_{-u, Q}(P)
        m_u, md_u, uQ = miller_function_ate_with_denom(Q, P, 0, -u)
        ##uQ = (uQ[0]/uQ[3], uQ[1]/(uQ[2]*uQ[3]))
        invZ3 = 1/(uQ[2]*uQ[3]) # 1/Z^3
        invZ2 = invZ3 * uQ[2]   # 1/Z^2
        uQ = (uQ[0]*invZ2, uQ[1]*invZ3)
        #assert Q.curve()(uQ) == -u*Q
        m_u2, md_u2, u2Q = miller_function_ate_with_denom(uQ, P, 0, -u, m0=m_u, m0_inv=md_u)
    else:
        m_u, md_u, uQ = miller_function_ate_with_denom(Q, P, 0, u)
        ##uQ = (uQ[0]/uQ[3], uQ[1]/(uQ[2]*uQ[3]))
        invZ3 = 1/(uQ[2]*uQ[3]) # 1/Z^3
        invZ2 = invZ3 * uQ[2]   # 1/Z^2
        uQ = (uQ[0]*invZ2, uQ[1]*invZ3)
        #assert Q.curve()(uQ) == u*Q
        m_u2, md_u2, u2Q = miller_function_ate_with_denom(uQ, P, 0, u, m0=m_u, m0_inv=md_u)
        #compute    f_{-u, Q}(P) = 1/(f_{u, Q}(P) * v_{uQ}(P))
        #then later f_{-u, pQ}(P) = 1/(f_{u, Q}(P) * v_{uQ}(P))^p
        m_u, md_u = md_u, m_u
        #vn, vd = vertical_j_denom(uQ, P) but uQ is now in affine coordinates
        vn = P[0] - uQ[0]
        md_u = md_u * vn
        uQ = (uQ[0], -uQ[1])
    #compute f_{-u, pQ}(P) = f_{-u,Q}(P)^p
    upQ = (uQ[0].frobenius(), uQ[1].frobenius())
    m_u = m_u.frobenius()
    md_u = md_u.frobenius()
    # add -upQ and p2Q
    p2Q = (Q[0].frobenius(2), Q[1].frobenius(2))
    # could do it in affine coordinates
    #(ln, ld), upp2Q = add_line_j_denom(upQ, p2Q, (P[0], P[1]))
    ln, upp2Q = add_line_affine_j(upQ, p2Q, (P[0], P[1]))
    ld = upp2Q[2]
    m = m_u * m_u2 * ln
    md = md_u * md_u2 * ld
    # divide by the vertical at upp2Q = -u^2Q
    # but next, multiply by the line through upp2Q and -u^2Q: this is the same vertical
    # the line through u2Q and [-u*P + p^2]Q should be a vertical
    #assert u2Q[0]*upp2Q[3] == upp2Q[0]*u2Q[3]
    #assert u2Q[1]*upp2Q[2]*upp2Q[3] == -upp2Q[1]*u2Q[2]*u2Q[3]
    return m / md

def miller_loop_opt_ate_bw13_bw19_seq(Q, P, u, k):
    """
    Optimal ate Miller loop for BW13, BW19
    with two independent sequential miller functions
    u^2 -u*p + p^2 = 0 mod r
    f_{u^2, Q}(P) * f_{-u, pQ}(P) * l_{[-u*p]Q,[p^2]Q}(P)
    """
    p = P.curve().base_field().characteristic()
    m_u2, md_u2, u2Q = miller_function_ate_with_denom(Q, P, 0, u**2)
    pQ = (Q[0].frobenius(), Q[1].frobenius())
    m_u, md_u, upQ = miller_function_ate_with_denom(pQ, P, 0, -u)
    # add -upQ and p2Q
    p2Q = (Q[0].frobenius(2), Q[1].frobenius(2))
    (ln, ld), upp2Q = add_line_j_denom(upQ, p2Q, (P[0], P[1]))
    m = m_u * m_u2 * ln
    md = md_u * md_u2 * ld
    # divide by the vertical at upp2Q = -u^2Q
    # but next, multiply by the line through upp2Q and -u^2Q: this is the same vertical
    # the line through u2Q and [-u*P + p^2]Q should be a vertical
    #assert u2Q[0]*upp2Q[3] == upp2Q[0]*u2Q[3]
    #assert u2Q[1]*upp2Q[2]*upp2Q[3] == -upp2Q[1]*u2Q[2]*u2Q[3]
    return m / md

####### Final exponentiation
def final_exp_easy_k_prime(m):
    """
    Easy part m^(p-1) of the final exponentiation in GF(p^k) prime k
    
    INPUT:
    - `m`: element of GF(p^k) absolute extension

    RETURN: m^(p-1) = m^p / m = frobenius(m) * (1/m)
    """
    return m.frobenius() / m

def final_exp_hard_k13(m, u):
    """
    exponent Phi_k(p)/r multiplied by 3*(u^12+u^11-u^9-u^8+u^6+u^5-u^3-u^2+1) has sparse form
    u > 0:
    (p^3 - u^3)*(p^6 + u^6)*u*(u + p)*((u^2-u+1)*(p*u^13 - p + u) - 3*p*u) + 3 + (u^2-u+1) * (p^3+1)*(p^6+1)*(p^2+p+1)
    u < 0:
    (p^3 - u^3)*(p^6 + u^6)*(-u)*(-u-p)*((u^2-u+1)*(-(-p*u^13 + p - u)) - 3*p*u) + 3 + (u^2-u+1) * (p^3+1)*(p^6+1)*(p^2+p+1)

    """
    # (u^2-u+1)
    if u < 0:
        mu = m**(-u)
        m1 = mu**(-u) * mu * m
    else:
        mu = m**u
        m1 = mu**u
        mu = cyclotomic_inverse_k13(mu) # because m**(-u) is needed later
        m1 = m1 * mu * m
    #assert m1 == m**(u**2-u+1)
    n0 = m1 * m1.frobenius(3)
    n0 = n0 * n0.frobenius(6)
    n0 = n0 * n0.frobenius() * n0.frobenius(2)

    #p = m.base_ring().characteristic()
    # p*u^13 - p + u
    mp = m1.frobenius()
    if u > 0: # (p*u^12 + 1)*u - p
        m2 = (mp**(u**12) * m1)**u * cyclotomic_inverse_k13(mp)
    else:     # -((p*u^12 + 1)*(-u) + p)
        m2 = cyclotomic_inverse_k13((mp**(u**12) * m1)**(-u) * mp)
    #assert m2 == m1**(p*u**13 - p + u)

    m2 = m2 * (mu**2 * mu).frobenius()
    #assert m2 == m**((u**2-u+1)*(p*u**13 - p + u) - 3*p*u)
    # u*(u+p) or -u*(-u-p)
    if u > 0:
        m2 = m2**u * m2.frobenius()
        m2 = m2**u
    else:
        m2 = m2**(-u) * cyclotomic_inverse_k13(m2.frobenius())
        m2 = m2**(-u)
    # (p^3-u^3)*(p^6+u^6)
    m2 = m2.frobenius(6) * m2**(u**6)
    if u < 0:
        m2 = m2.frobenius(3) * m2**(-u**3)
    else:
        m2 = m2.frobenius(3) * cyclotomic_inverse_k13(m2**(u**3))
    return m2 * n0 * m**2 * m

def final_exp_hard_k19(m, u):
    """
    exponent is Phi_19(p)/r*3*(u^18+u^17-u^15-u^14+u^12+u^11-u^9-u^8+u^6+u^5-u^3-u^2+1)
    (-u^15+p^3*u^12-p^6*u^9+p^9*u^6-p^12*u^3+p^15)*u*(p+u)*((u^2-u+1)*(p*u^19-p+u) -3*p*u) + 3 + (u^2-u+1)*((p^9+1)*(p^2+p+1)*(p^6+p^3+1))
    (p^3-u^3)*(p^12+p^6*u^6+u^12)*u*(p+u)*((u^2-u+1)*(p*u^19-p+u) -3*p*u) + 3 + (u^2-u+1)*((p^9+1)*(p^2+p+1)*(p^6+p^3+1))
    5 * exp(u) + exp(u**3) + 2*exp(u**6) + exp(u**18)
    """
    # (u^2-u+1)
    if u < 0:
        mu = m**(-u)                 # exp(u)
        m1 = mu**(-u) * mu * m       # exp(u)
    else:
        mu = m**u
        m1 = mu**u
        mu = cyclotomic_inverse_k19(mu) # because m**(-u) is needed later
        m1 = m1 * mu * m
    assert m1 == m**(u**2-u+1)
    n0 = m1 * m1.frobenius(9)
    n0 = n0 * n0.frobenius(3) * n0.frobenius(6)
    n0 = n0 * n0.frobenius() * n0.frobenius(2)

    #p = m.base_ring().characteristic()
    # p*u^19 - p + u
    mp = m1.frobenius()
    if u > 0: # (p*u^18 + 1)*u - p
        m2 = (mp**(u**18) * m1)**u * cyclotomic_inverse_k19(mp)     # exp(u**18) + exp(u)
    else:     # -((p*u^18 + 1)*(-u) + p)
        m2 = cyclotomic_inverse_k19((mp**(u**18) * m1)**(-u) * mp)
    #assert m2 == m1**(p*u**19 - p + u)

    m2 = m2 * (mu**2 * mu).frobenius()
    #assert m2 == m**((u**2-u+1)*(p*u**19 - p + u) - 3*p*u)
    # u*(u+p) or -u*(-u-p)
    if u > 0:
        m2 = m2**u * m2.frobenius()                                 # exp(u)
        m2 = m2**u                                                  # exp(u)
    else:
        m2 = m2**(-u) * cyclotomic_inverse_k19(m2.frobenius())
        m2 = m2**(-u)
    # = (p^3-u^3)*(p^12 + p^6*u^6 + u^12)
    # = (p^3-u^3)*(p^12 + (p^6 + u^6)*u^6)
    m6 = m2.frobenius(6)
    m2 = m6.frobenius(6) * (m6 * m2**(u**6))**(u**6)                # exp(u**6) + exp(u**6)
    if u < 0:
        m2 = m2.frobenius(3) * m2**(-u**3)                          # exp(u**3)
    else:
        m2 = m2.frobenius(3) * cyclotomic_inverse_k19(m2**(u**3))
    return m2 * n0 * m**2 * m
