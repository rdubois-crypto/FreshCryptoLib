from pairing import *

def miller_loop_opt_ate_kss18(Q, P, u):
    """Miller loop f_{u,Q}(P) f_{3,Q}(P)^p l_{[u]Q,3\pi(Q)}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fqd) of order r
    - `P`: point on E(Fp)
    - `u`: seed for curve parameters

    The curve coefficient a is zero.
    """
    m, uQ = miller_function_ate(Q, P, 0, u, m0=1)
    # it could be possible the first time of doubling to take into account that both P and Q are in affine coordinates, and also to store the first result to avoid re-computing 2*Q later.
    PP = (P[0], P[1])
    QQ = (Q[0], Q[1], 1, 1)
    l2, Q2 = double_line_j(QQ, PP, 0) # this is already computed in miller_function_ate(...) above
    QQ = (Q[0],Q[1])
    l3, Q3 = add_line_j(Q2, QQ, PP)

    piQ3 = [(Q3[0]).frobenius(), (Q3[1]).frobenius(), (Q3[2]).frobenius(), (Q3[3]).frobenius()]
    li, uQ_3piQ = add_line_j_with_z(uQ, piQ3, PP)

    return m * (l2*l3).frobenius() * li

def miller_loop_opt_ate_aklgl_kss18(Q, P, b_t, u, Fq6, map_Fq6_Fpk, D_twist=False, xi=None, check=False):
    """KSS18 optimal ate Miller loop f_{u,Q}(P) f_{3,Q}(P)^p l_{[u]Q,3\pi(Q)}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fqd) of order r
    - `P`: point on E(Fp)
    - `b_t`: the curve coefficient of the twist curve 
    - `u`: seed for curve parameters
    - `Fq6`: the extension compatible with the twist
    - `map_Fq6_Fpk`: map from Fq6 (with explicit degree 6) to absolute extension Fpk
    - `D_twist`: is it a D-twist of a M-twist

    The curve coefficient a is zero.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    m, uQ = miller_function_ate_aklgl(Q, P, b_t, u, Fq6, D_twist=D_twist)
    if check:
        uQ_test = u*Q
        uQ_aff = (uQ[0]/uQ[2], uQ[1]/uQ[2])
        assert uQ_aff[0] == uQ_test[0]
        assert uQ_aff[1] == uQ_test[1]
    # it could be possible the first time of doubling to take into account that both P and Q are in affine coordinates, and also to store the first result to avoid re-computing 2*Q later.
    PP = (P[0], P[1])
    QQ = (Q[0], Q[1], 1)
    l2, Q2 = double_line_h_a0_twist6_aklgl(QQ,PP,b_t,D_twist=D_twist)
    # this is already computed in miller_function_ate_aklgl(...) above
    QQ = (Q[0],Q[1])
    l3, Q3 = add_line_h_a0_twist6_aklgl(Q2,QQ,PP,D_twist=D_twist)

    # how to compute a frobenius 'before' the twist map?
    if check:
        p = Fq6.characteristic()
        E2 = Q.parent()
    if D_twist:
        # psi_sextic_d_twist(x,y) = (x*xi^2, y*xi^3)
        # in homogeneous coordinates: (X,Y,Z) ~ (X/Z,Y/Z,1) = (x,y)
        # a^6 = xi
        # pi(x*a^2, y*a^3) = (x.frobenius()*a^2.frobenius(), y.frobenius()*a^3.frobenius())
        # a^p = a^((p-1)/6 * 6)*a = xi^((p-1)/6) * a
        # a^2^p = xi^((p-1)/3) * a^2
        # a^3^p = xi^((p-1)/2) * a^3 = -a^3
        w = xi**((p-1)//6)
        piQ3 = [(Q3[0]).frobenius() * w**2, (Q3[1]).frobenius() * w**3, (Q3[2]).frobenius()]
    else:
        w = 1/xi**((p-1)//6)
        piQ3 = [(Q3[0]).frobenius() * w**2, (Q3[1]).frobenius() * w**3, (Q3[2]).frobenius()]
        
    li, uQ_3piQ = add_line_h_a0_twist6_aklgl_with_z(uQ, piQ3, PP, D_twist=D_twist)
    if D_twist:
        l2_l3 = sparse_sparse_mult_d6_twist(l2[0], l2[1], l2[3], l3[0], l3[1], l3[3], xi, Fq6)
        m_li = sparse_mult_d6_twist(li[0], li[1], li[3], m, xi, Fq6)
    else:
        l2_l3 = sparse_sparse_mult_m6_twist(l2[0], l2[2], l2[3], l3[0], l3[2], l3[3], xi, Fq6)
        m_li = sparse_mult_m6_twist(li[0], li[2], li[3], m, xi, Fq6)

    return map_Fq6_Fpk(m_li) * map_Fq6_Fpk(l2_l3).frobenius()

def final_exp_hard_kss18_v0(f, u):
    """
    formulas from Guide to pairing-based crypto book
    Algorithm 7.16 in chapter 16
    Cost: 6 Sq + 53 M + 29 Frobenius powers + 7 exponentiations to u + 8 Inversions with Frobenius power p^9
    """
    f_ = f.frobenius(9)
    A = f_.frobenius()    # 1
    t0 = A**2             # 2
    B = f_.frobenius(4)   # 3
    t0 = t0 * B           # 4
    fx = f**u
    fx_ = fx.frobenius(9)
    C = fx_.frobenius()   # 5
    t1 = t0 * C           # 6
    B = fx.frobenius(3)   # 7
    t0 = t1 * B           # 8
    fx2 = fx**u
    fx3 = fx2**u
    B = (fx2.frobenius() * fx3.frobenius(2)).frobenius(9) # 9
    t1 = t1 * B           # 10
    A = f.frobenius(2)    # 11
    B = fx_.frobenius(4)  # 12
    t6 = A * B            # 13
    t0 = t0 * B           # 14
    fx4 = fx3**u
    fx5 = fx4**u
    fx5_ = fx5.frobenius(9)
    A = fx5_.frobenius(4) * f.frobenius(5) # 15
    t4 = A * B            # 16
    B = fx                # 17
    t2 = t0 * B           # 18
    t0 = t0 * t1          # 19
    A = fx2.frobenius(3)  # 20
    t1 = A * C            # 21
    fx4_ = fx4.frobenius(9)
    B = fx4_.frobenius(2) # 22
    t3 = A * B            # 23
    t2 = t1 * t2          # 24
    B = fx4_.frobenius(4) # 25
    t5 = t1 * B           # 26
    fx2_ = fx2.frobenius(9)
    B = fx2_.frobenius(2) # 27
    t1 = t2 * B           # 28
    B = fx4.frobenius(3)  # 29
    t8 = t2 * B           # 30
    B = fx2               # 31
    t2 = t1 * B           # 32
    B = fx4_.frobenius()  # 33
    t1 = t1 * B           # 34
    t0 = t2 * t0          # 35
    B = fx5.frobenius(3)  # 36
    t7 = t2 * B           # 37
    t0 = t0**2            # 38
    B = fx2_.frobenius(4) # 39
    t2 = t0 * B           # 40
    fx3_ = fx3.frobenius(9)
    B = fx3_.frobenius() * fx3.frobenius(3) # 41
    t0 = t2 * B           # 42
    t2 = t2 * t8          # 43
    t1 = t0 * t1          # 44
    t0 = t0 * t7          # 45
    t3 = t1 * t3          # 46
    t1 = t1 * t6          # 47
    B = fx3 * fx3_.frobenius(4) # 48
    t6 = t3 * B           # 49
    fx6 = fx5**u
    A = fx6.frobenius(3)  # 50
    t3 = A * B            # 51
    t2 = t6 * t2          # 52
    fx6_ = fx6.frobenius(9)
    B = fx5 * fx5_.frobenius() * fx6_.frobenius(2) * fx3.frobenius(5) # 53
    t6 = t6 * B           # 54
    t2 = t2 * t5          # 55
    fx7_ = fx6_**u
    B = fx6 * fx7_.frobenius(2) * fx4.frobenius(5) # 56
    t5 = t5 * B           # 57
    t2 = t2**2            # 58
    B = fx4 * fx5_.frobenius(2) * fx2.frobenius(5) # 59
    t2 = t2 * B           # 60
    t0 = t0**2            # 61
    t0 = t0 * t6          # 62
    t1 = t2 * t1          # 63
    t2 = t2 * t5          # 64
    t1 = t1**2            # 65
    t1 = t1 * t4          # 66
    t1 = t1 * t0          # 67
    t0 = t0 * t3          # 68
    t0 = t0 * t1          # 69
    t1 = t1 * t2          # 70
    t0 = t0**2            # 71
    t0 = t0 * t1          # 72
    return t0             # 73

def final_exp_hard_kss18(m, u):
    """
    (p^6 - p^3 + 1)/r
    1. replace p by (r*c + t-1)
    2. expand the powers, simplify by r
    3. re-factor r*c into (p+1-t)
    There is a more systematic way to do it: https://eprint.iacr.org/2020/875

    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    px = (x**8 + 5*x**7 + 7*x**6 + 37*x**5 + 188*x**4 + 259*x**3 + 343*x**2 + 1763*x + 2401)/21
    rx = (x**6 + 37*x**3 + 343)/343 # 343 = 7^3
    cx = (x**2 + 5*x + 7)*49/3
    tx = (x**4 + 16*x + 7)/7

    Formula can be deduced from https://eprint.iacr.org/2020/875:
    B = px^5 + px^4*(tx-1) + px^3*(tx-1)^2 + px^2*((tx-1)^3-1) + px*((tx-1)^4 - tx + 1) + (tx-1)^5  - (tx-1)^2
    D = ((tx-1)^6 - (tx-1)^3 + 1) // rx
    cx*B + D == (px^6 - px^3 + 1) // rx

    Magma:
    px := (x^8 + 5*x^7 + 7*x^6 + 37*x^5 + 188*x^4 + 259*x^3 + 343*x^2 + 1763*x + 2401)/21;
    rx := (x^6 + 37*x^3 + 343)/343; // 343 = 7^3
    cx := (x^2 + 5*x + 7)*49/3;
    tx := (x^4 + 16*x + 7)/7;
    B := px^5 + px^4*(tx-1) + px^3*(tx-1)^2 + px^2*((tx-1)^3-1) + px*((tx-1)^4 - tx + 1) + (tx-1)^5  - (tx-1)^2;
    D := ((tx-1)^6 - (tx-1)^3 + 1) div rx;
    cx*B + D eq (px^6 - px^3 + 1) div rx;
    M := Matrix(QQx, 6, 6, [B,0,0,0,0,0, -px,1,0,0,0,0,  0,-px,1,0,0,0,  0,0,-px,1,0,0,  0,0,0,-px,1,0,  0,0,0,0,-px,1]);
    R := LLL(M);

    The 4th row of R corresponds to the formulas in https://link.springer.com/content/pdf/10.1007%2F978-3-642-28496-0_25.pdf
    Fuentes-Castaneda  -- Rodriguez-Henriquez SAC 2011 page 420.
    Claimed cost: 8 Sq + 52 M + Frobenius powers + 7 exponentiations to u
    But then it is possible to factor the leading term (of p^5)
    -> writing Euclidean division by (x^4 + 5*x^3 + 7*x^2 + 3) for each term in p^i, we obtain:

    l0 =      -u**6 - 5*u**5 - 7*u**4 - 21*u**3 - 108*u**2 - 147*u #= (-x^2) * (x^4 + 5*x^3 + 7*x^2 + 3) + 21 * (-x^3 - 5*x^2 - 7*x)
    l1 =        5*u**5 + 25*u**4 + 35*u**3 + 98*u**2 + 505*u + 686 #= (5*x) * (x^4 + 5*x^3 + 7*x^2 + 3) + 98 * (x^2 + 5*x + 7)
    l2 = u**7 + 5*u**6 + 7*u**5 + 19*u**4 + 98*u**3 + 133*u**2 - 6 #= (x^3 + 19) * (x^4 + 5*x^3 + 7*x^2 + 3) + 1 * (-63)
    l3 =  -2*u**6 - 10*u**5 - 14*u**4 - 35*u**3 - 181*u**2 - 245*u #= (-2*x^2) * (x^4 + 5*x^3 + 7*x^2 + 3) + 35 * (-x^3 - 5*x^2 - 7*x)
    l4 =        3*u**5 + 15*u**4 + 21*u**3 + 49*u**2 + 254*u + 343 #= (3*x) * (x^4 + 5*x^3 + 7*x^2 + 3) + 49 * (x^2 + 5*x + 7)
    l5 =                               -u**4 - 5*u**3 - 7*u**2 - 3 #= (-1) * (x^4 + 5*x^3 + 7*x^2 + 3) + 1 * (0)

    Then we write the Euclidean division by (x^2 + 5*x + 7) for each remainder. We obtain
    Our own formula: assuming easy inversion f^(-1) = f^(p^9)
    e4 eq (x^4 + 5*x^3 + 7*x^2 + 3)*(-px^5 + 3*x*px^4 -2*x^2*px^3 + (x^3+16 + 3)*px^2 + 5*x*px - x^2) + 7*((x^2 + 5*x + 7)*(7*(px^4 + 2*px) - x*(5*px^3 + 3)) - 9*px^2);
    e4 eq (x^4 + 5*x^3 + 7*x^2 + 3)*(-px^5 + x*(3*px^4 + 5*px) -x^2*(2*px^3 + 1) + x^3*px^2 + (16 + 3)*px^2) + 7*((x^2 + 5*x + 7)*(7*(px^4 + 2*px) - x*(5*px^3 + 3)) - 9*px^2);

    Cost: 7 exp(u) + 19 S + 26 M + 10 Frobenius powers + 6 Inversions with Frobenius power p^9

    """
    # 1. x^2 + 5*x + 7
    m1 = m**u
    m2 = m1**u
    m01 = m * m1 # m^(x+1)           # M
    f = m2 * ((m01**2 * m)**2 * m01) # 3M + 2S
    #assert f == m**(u**2 + 5*u + 7)
    f1 = f.frobenius(4) * f.frobenius()**2 # M + S
    fu = f**u #
    fu_ = fu.frobenius(9) # inv
    f2 = fu_.frobenius(3)
    f12 = f1 * f2                               # M
    f3 = (f12**2 * (f1 * fu_))**2 * (f12 * fu_) # 4M + 2S
    #f3 = f1**7 * f2**5 * fu_**3
    f9 = (m.frobenius(2)).frobenius(9)
    f9 = ((f9**2)**2)**2 * f9                  # 3 S + M
    f4 = f3 * f9                               # M
    f5 = f4.frobenius(9) * ((f4**2)**2)**2     # 3 S + M
    #i4 = 49 * (u**2 + 5*u + 7)
    #i3 = 35 * (-u**3 - 5*u**2 - 7*u)
    #i2 = -63
    #i1 = 98 * (u**2 + 5*u + 7)
    #i0 = 21 * (-u**3 - 5*u**2 - 7*u)
    # assert 7*((x^2 + 5*x + 7)*(7*(px^4 + 2*px) - x*(5*px^3 + 3)) - 9*px^2) == i0 + i1*px + i2 * px^2 + i3*px^3 + i4*px^4
    #assert f5 == m**i0 * m.frobenius()**i1 * m.frobenius(2)**i2 * m.frobenius(3)**i3 * m.frobenius(4)**i4
    g = fu**u * m**2 * m # m^(x^4 + 5*x^3 + 7*x^2 + 3) # 2M + S + exp(u)
    #assert g == m**(u**4 + 5*u**3 + 7*u**2 + 3)
    gu = g**u
    gu2 = gu**u
    gu3 = gu2**u
    g5 = (g.frobenius(5)).frobenius(9)
    #gx = gu.frobenius()**5 * gu.frobenius(4)**3
    gup = gu.frobenius()
    gup4 = gu.frobenius(4)
    gx = (gup**2 * gup4)**2 * gup * gup4 # 2 S + 3M
    gx2 = (gu2 * gu2.frobenius(3)**2).frobenius(9) #S + M
    gx3 = gu3.frobenius(2)
    #g2 = (g**19).frobenius(2)
    g4 = (g**2)**2                            # 2S
    g2 = ((g4)**2)**2 * g4 * g.frobenius(9)   # 2S + 2M
    #assert g2 == g**19
    g2 = g2.frobenius(2)
    h = (g5 * gx * gx2 * gx3 * g2) * f5
    #assert h == m.frobenius(5)**l5 * m.frobenius(4)**l4 * m.frobenius(3)**l3 * m.frobenius(2)**l2 * m.frobenius()**l1 * m**l0
    return h

def final_exp_hard_kss18_v1(m, u):
    """
    Final exponentiation, hard part 3/49*u**2 * (p^6 - p^3 + 1)/r

    From Faster Final Exponentiation on the KSS18 Curve,
    Shiping Cai, Zhi Hu, and Chang-An Zhao
    https://eprint.iacr.org/2021/1309
    7 exp(u) + 24 M + 11 S + 7 Frob(9) + 5 Frobenius
    """
    t0 = m**u
    t1 = m**2
    t4 = m * t1                     # m**3
    t2 = t1 * t4                    # m**5
    t1 = t1 * t2                    # m**7
    t2 = (t0 * t2)**u               # m**(5*u + u^2)
    c = t1 * t2                     # m**(7 + 5*u + u^2) = m**l6
    #l6 = (7 + 5*u + u**2); assert c == m**l6
    t0 = (c**2 * c)**2 * c          # m**(7*l6)
    t1 = t0**2                      # c**14 = m**(14*l6)
    t3 = (c**u).frobenius(9)        # c**(-u) = m**(-u*l6)
    c = (t3**u).frobenius(9) * t4   # c**(u^2) * m**3 = m**(u^2*l6 + 3) = m**l5
    #l5 = (u**2 * l6 + 3); assert c == m**l5
    t4 = c**u                       # m**(u*l5)
    t2 = t4.frobenius(9)            # m**(-u*l5)
    t41 = t4 * t1                   # m**(u*l5+14*l6)
    t1 = t41**2 * t41 * t0          # m**(3*u*l5 + 42*l6 + 7*l6) = m**(3*u*l5 + 49*l6)
    t2 = (t2**u).frobenius(9)       # m**(u**2*l5)
    t0 = t1.frobenius(9)            # m**(-3*u*l5 - 49*l6)
    #l4 = -3*u*l5 - 49*l6; assert t0 == m**l4
    t1 = t0**2 * t4   # m**(2*l4 + u*l5)
    #l1 = 2*l4 + u*l5; assert t1 == m**l1
    t1 = t1.frobenius()
    t4 = t1 * (c.frobenius() * t0).frobenius(4) # m**(l1*p + l4*p^4 + l5*p^5)
    #assert t4 == (m**l1).frobenius() * (m**l4).frobenius(4) * (m**l5).frobenius(5)
    t3_ = t3.frobenius(9)
    t3 = (t3_**2 * t3_)**2 * t3_  # m**(7*u*l6)
    t1 = t3**2                    # m**(14*u*l6)
    t0 = (t1 * t2)**2 * t3        # m**(2*u^2*l5 + 35*u*l6)
    #l3 = 2*u**2*l5 + 35*u*l6; assert t0 == m**l3
    t3 = t1 * t3                  # m**(21*u*l6)
    t1 = t2 * t3                  # m**(u**2*l5 + 21*u*l6)
    t4 = t1 * t4                  # m**(u**2*l5 + 21*u*l6 + l1*p + l4*p4 + l5*p5)
    t1 = (t1**u).frobenius(9)     # m**(-u**3*l5 - 21*u**2*l6)
    t2 = c**2                     # m**(2*l5)
    t1 = t1 * t2                  # m**(-u**3*l5 + 2*l5 - 21*u**2*l6)
    t0 = (t0.frobenius() * t1).frobenius(2)
    c = t4 * t0
    #l0 = 2*l3 + u*l4
    #l2 =-u*l0 + 2*l5
    return c
