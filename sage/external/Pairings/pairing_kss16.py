from pairing import *

def miller_loop_opt_ate_kss16(Q, P, a, u):
    """Miller loop (f_{u,Q}(P) l_{[u]Q,\pi(Q)}(P))^{p^3} l_{Q,Q}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fq4) of order r
    - `P`: point on E(Fp)
    - `a`: curve coefficient in short Weierstrass y^2 = x^3 + a*x + b
    - `u`: seed for curve parameters
    """
    m, uQ = miller_function_ate(Q, P, a, u, m0=1)
    piQ = [(Q[0]).frobenius(), (Q[1]).frobenius()]
    PP = (P[0], P[1])
    l1, uQ_piQ = add_line_j(uQ, piQ, PP)
    QQ = (Q[0], Q[1])
    l2, Q2 = double_line_affine_j(QQ, PP, a)

    return (l1 *  m).frobenius(3) * l2

def miller_loop_opt_ate_kss16_v2(Q, P, a, u):
    """Miller loop l_{Q,Q}(P) * f_{u, pi3(Q)} * l_{pi4Q, u*pi3(Q)}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fq4) of order r
    - `P`: point on E(Fp)
    - `a`: curve coefficient in short Weierstrass y^2 = x^3 + a*x + b
    - `u`: seed for curve parameters
    """
    pi3Q = [(Q[0]).frobenius(3), (Q[1]).frobenius(3)]
    pi4Q = [(Q[0]).frobenius(4), (Q[1]).frobenius(4)]
    m, upi3Q = miller_function_ate(pi3Q, P, a, u, m0=1)
    QQ = (Q[0],Q[1])
    PP = (P[0],P[1])
    l2, Q2 = double_line_affine_j(QQ, PP, a)
    l1, pi4Q_u_pi3Q = add_line_j(upi3Q, pi4Q, PP)

    return l2 *  m * l1

def miller_loop_opt_ate_kss16_cln_b0(Q, P, a, u):
    """Miller loop (f_{u,Q}(P) l_{[u]Q,\pi(Q)}(P))^{p^3} l_{Q,Q}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fq4) of order r
    - `P`: point on E(Fp)
    - `a`: curve coefficient in short Weierstrass y^2 = x^3 + a*x
    - `u`: seed for curve parameters
    """
    m, uQ = miller_function_ate_cln_b0(Q, P, a, u)
    piQ = [(Q[0]).frobenius(), (Q[1]).frobenius()]
    PP = (P[0], P[1])
    l1, uQ_piQ = add_line_cln_b0(uQ, piQ, PP)
    QQ = (Q[0], Q[1], 1)
    l2, Q2 = double_line_cln_b0(QQ, PP, a)

    return (l1 *  m).frobenius(3) * l2

def miller_loop_opt_ate_kss16_cln_b0_v2(Q, P, a, u):
    """Miller loop l_{Q,Q}(P) * f_{u, pi3(Q)} * l_{pi4Q, u*pi3(Q)}(P)

    INPUT:
    - `Q`: point on E(Fpk) or E(Fq4) of order r
    - `P`: point on E(Fp)
    - `a`: curve coefficient in short Weierstrass y^2 = x^3 + a*x + b
    - `u`: seed for curve parameters
    """
    pi3Q = [(Q[0]).frobenius(3), (Q[1]).frobenius(3)]
    pi4Q = [(Q[0]).frobenius(4), (Q[1]).frobenius(4)]
    m, upi3Q = miller_function_ate_cln_b0(pi3Q, P, a, u, m0=1)
    QQ = (Q[0],Q[1], 1)
    PP = (P[0],P[1])
    l2, Q2 = double_line_cln_b0(QQ, PP, a)
    l1, pi4Q_u_pi3Q = add_line_cln_b0(upi3Q, pi4Q, PP)

    return l2 *  m * l1

def gfp16_frobenius(s, i):
    return s.frobenius(i)

def inversion_easy_gfp16(s):
    """1/s = s^(p^8) in the cyclotomic subgroup """
    return gfp16_frobenius(s, 8)

def final_exp_hard_kss16(m, u):
    """
    k := 16;
    D := 1;
    px := (x^10 + 2*x^9 + 5*x^8 + 48*x^6 + 152*x^5 + 240*x^4 + 625*x^2 + 2398*x + 3125)/980;
    rx := (x^8 + 48*x^4 + 625)/61250; // 625 = 5^4, 61250 = 2*5^4*7^2
    tx := (2*x^5 + 41*x + 35)/35;
    cx := 125 * (x^2 + 2*x + 5)/2; // C such that P+1-T = C*R
    exponent := (px^8+1) div rx;

    C5[1] =   -2*x^8 - 4*x^7 - 10*x^6 - 55*x^4 - 222*x^3 - 275*x^2 = (  -2*x^3) * (x^5 + 2*x^4 + 5*x^3 + 56) +  55 * (-x^4 - 2*x^3 - 5*x^2) = ... + ((-x^2) * (x^2 + 2*x + 5) + 0
    C5[2] =      4*x^7 + 8*x^6 + 20*x^5 + 75*x^3 + 374*x^2 + 375*x = (   4*x^2) * (x^5 + 2*x^4 + 5*x^3 + 56) +  75 * (x^3 + 2*x^2 + 5*x)    = ... + ((x)    * (x^2 + 2*x + 5) + 0
    C5[3] =         2*x^6 + 4*x^5 + 10*x^4 + 125*x^2 + 362*x + 625 = (     2*x) * (x^5 + 2*x^4 + 5*x^3 + 56) + 125 * (x^2 + 2*x + 5)        = ... + ((1)    * (x^2 + 2*x + 5) + 0
    C5[4] = x^9 + 2*x^8 + 5*x^7 + 24*x^5 + 104*x^4 + 120*x^3 - 196 = (x^4 + 24) * (x^5 + 2*x^4 + 5*x^3 + 56) +   1 * (-1540)                = ... + ((0)    * (x^2 + 2*x + 5) + -1540
    C5[5] =        -x^8 - 2*x^7 - 5*x^6 - 10*x^4 - 76*x^3 - 50*x^2 = (    -x^3) * (x^5 + 2*x^4 + 5*x^3 + 56) +  10 * (-x^4 - 2*x^3 - 5*x^2) = ... + ((-x^2) * (x^2 + 2*x + 5) + 0
    C5[6] =    -3*x^7 - 6*x^6 - 15*x^5 - 100*x^3 - 368*x^2 - 500*x = (  -3*x^2) * (x^5 + 2*x^4 + 5*x^3 + 56) + 100 * (-x^3 - 2*x^2 - 5*x)   = ... + ((-x)   * (x^2 + 2*x + 5) + 0
    C5[7] =     11*x^6 + 22*x^5 + 55*x^4 + 250*x^2 + 1116*x + 1250 = (    11*x) * (x^5 + 2*x^4 + 5*x^3 + 56) + 250 * (x^2 + 2*x + 5)        = ... + ((1)    * (x^2 + 2*x + 5) + 0
    C5[8] =                         -7*x^5 - 14*x^4 - 35*x^3 - 392 = (      -7) * (x^5 + 2*x^4 + 5*x^3 + 56) +   1 * (0)                    = ... + ((0)    * (x^2 + 2*x + 5) + 0

    e5 := (x^3*(x^2+2*x+5)+56)*(-(8-1)*px^7 + (16+8)*px^3 + 4*px + x*((8+2+1)*px^6 + 2*px^2) + x^2*(-(4-1)*px^5) + x^4*px^3 -x^3*(px^4+2)) + 5*((x^2+2*x+5)*(25*(px^2 + 2*px^6) -5*x*(4*px^5 - 3*px) -x^2*(2*px^4 + 11)) - 308*px^3);
    e5 eq -14/125*x^3 * exponent;

    # total cost 33 S + 31 M + 2 exp(u+1) + 7 exp(u) + 14*f16 + 8 inv_frob_8
    """
    # 0. compute m^(308*p^3)
    m4 = (m**2)**2
    m8 = m4**2                                   # 3 S
    f0 = (((((m8 * m)**2 * m)**2)**2 * m)**2)**2 # 5 S + 3 M
    #assert f0 == m**308
    f0 = f0.frobenius(3)
    f0 = f0.frobenius(8)                         # cyclo inv
    #1. compute m^(x^2+2*x+5) = m^((x+1)^2 + 4)
    f1 = (m**(u+1))**(u+1) * m4                  # M + 2*exp(u+1)
    #assert f1 == m**(u**2 + 2*u + 5)
    # compute the term in (x^2+2*x+5)
    f1a = f1.frobenius(2) * f1.frobenius(6)**2   # (px^2 + 2*px^6)
    f1a = (f1a**2)**2 * f1a                      # power 5 -> 5*(px^2 + 2*px^6)
    f1x = f1**u
    # -5*x*(4*px^5 - 3*px) = -5*x*(4*px^5 - (4-1)*px) = -5*x*(4*(px^5-px) + px)
    f1xp = f1x.frobenius()
    f1b = ((f1x.frobenius(5) * f1xp.frobenius(8))**2)**2 * f1xp
    #assert f1b == f1x.frobenius(5)**4 * f1x.frobenius()**(-3)
    f1ab = f1a * f1b.frobenius(8)
    f1ab = (f1ab**2)**2 * f1ab
    # f1ab = f1^(25*(px^2 + 2*px^6) -5*x*(4*px^5 - 3*px))
    f1x2 = f1x**u
    # -x^2*(2*px^4 + 11) = (2*px^4 + 8+2+1)
    f1c = (f1x2.frobenius(4) * ((f1x2**2)**2 * f1x2))**2 * f1x2 # 3 S + 3 M
    #assert f1c == (f1.frobenius(4)**2 * f1**11)**(u**2)
    f1c = f1c.frobenius(8)
    f1abc = f1ab * f1c
    f01 = f1abc * f0
    f01 = (f01**2)**2 * f01
    # all the term 5*((x^2+2*x+5)*(25*(px^2 + 2*px^6) -5*x*(4*px^5 - 3*px) -x^2*(2*px^4 + 11)) - 308*px^3)
    # so far 20 S + 16 M + 2 exp(u+2) + 2 exp(u) + 6 Frob + 4 cyclo-inv with Frob-8
    # (x^2+2*x+5)*(50*px^6 -20*x*px^5 -2*x^2*px^4 + 25*px^2 + 15*x*px - 11*x^2)
    #assert f1abc == f1.frobenius(6)**50 * f1.frobenius(5)**(-20*u) * f1.frobenius(4)**(-2*u**2) * f1.frobenius(2)**25 * f1.frobenius()**(15*u) * f1**(-11*u**2)

    # 2. compute m^(x^3*(x^2+2*x+5)+56) from f1 = m^(x^2+2*x+5)
    f1x3 = f1x2**u
    m56 = (((m8 * m.frobenius(8))**2)**2)**2 # M + 3 S
    assert m56 == m**56
    f2 = f1x3 * m56                           # M
    # exp(u) + 2M + 3 S

    # 3. compute ((-8+1)*px^7 + (16+8)*px^3 + x*((8+2+1)*px^6 + 2*px^2) + x^2*(-(4-1)*px^5 + 4*px) + x^4*px^3 -x^3*(px^4+2))
    f2p3 = f2.frobenius(3)
    f2p7 = f2.frobenius(7)
    f2a = (((f2p3**2 * f2p3 * f2p7.frobenius(8))**2)**2)**2 * f2p7
    #assert f2a == f2.frobenius(3)**24 * f2.frobenius(7)**(-7)
    f2x = f2**u
    f2xp6 = f2x.frobenius(6)
    f2b = ((f2xp6**2)**2 * f2xp6 * f2x.frobenius(2))**2 * f2xp6 # 3 S + 3 M
    #assert f2b == f2.frobenius(6)**(11*u) * f2.frobenius(2)**(2*u)
    f2x2 = f2x**u
    f2x2p5 = f2x2.frobenius(5)
    f2c = ((f2x2p5.frobenius(8) * f2x2.frobenius())**2)**2 * f2x2p5
    #assert f2c == f2.frobenius(5)**(-3*u**2) * f2.frobenius()**(4*u**2)
    f2x3 = f2x2**u
    f2e = (f2x3.frobenius(4) * f2x3**2).frobenius(8)
    #assert f2e == f2.frobenius(4)**(-u**3) * f2**(-2*u**3)
    f2x4 = f2x3**u
    f2d = f2x4.frobenius(3)
    #assert f2d == f2.frobenius(3)**(u**4)
    # more 10 S + 4 exp(u) + 9 M
    fC = f2a * f2b * f2c * f2d * f2e
    # -7*px^7 + 11*x*px^6 -3*x^2*px^5 -x^3*px^4 + (x^4+24)*px^3 + 2*x*px^2 + 4*x^2*px -2*x^3
    # p = m.base_ring().characteristic()
    #assert fC == f2.frobenius(7)**(-7) * f2.frobenius(6)**(11*u) * f2.frobenius(5)**(-3*u**2) * f2.frobenius(4)**(-u**3) * f2.frobenius(3)**(u**4+24) * f2.frobenius(2)**(2*u) * f2.frobenius()**(4*u**2) * f2**(-2*u**3)
    # total cost 33 S + 31 M + 2 exp(u+1) + 7 exp(u)
    return (fC) * f01

def final_exp_hard_kss16_ghammam(m, u):
    """Hard part of final exponentiation Phi_8(p)/r

    Loubna Ghammam
    https://tel.archives-ouvertes.fr/tel-01469981v1
    Utilisation des Couplages en Cryptographie asymétrique pour la micro-électronique.
    PhD thesis, Universite de Rennes I, France
    Section 4.3 p. 107
    Eq (4.9) p. 114 Section 4.3.

    Loubna Ghammam and Emmanuel Fouotsa. Adequate Elliptic Curves for Computing the Product of n Pairings.
    International Workshop on the Arithmetic of Finite Fields, WAIFI 2016
    https://link.springer.com/chapter/10.1007/978-3-319-55227-9_3
    https://eprint.iacr.org/2016/472

    s = 14*(u//5)**3
    # for KSS16 curve we have u = 25,45 mod 70 <=> +/- 25 mod 70 --> 5 | u
    m_0 = 2*u**8 + 4*u**7 + 10*u**6 + 55*u**4 + 222*u**3 + 275*u**2
    m_1 = -4*u**7 - 8*u**6 - 20*u**5 - 75*u**3 - 374*u**2 - 375*u
    m_2 = -2*u**6 - 4*u**5 - 10*u**4 - 125*u**2 - 362*u - 625
    m_3 = -*u**9 - 2*u**8 - 5*u**7 - 24*u**5 - 104*u**4 - 120*u**3 + 196
    m_4 = *u**8 + 2*u**7 + 5*u**6 + 10*u**4 + 76*u**3 + 50*u**2
    m_5 = 3*u**7 + 6*u**6 + 15*u**5 + 100*u**3 + 368*u**2 + 500*u
    m_6 = -11*u**6 - 22*u**5 - 55*u**4 - 250*u**2 - 1116*u - 1250
    m_7 = 7*u**5 + 14*u**4 + 35*u**3 + 392
    m_0 + m_1*p + m_2*p**2 + m_3*p**3 + m_4*p**4 + m_5*p**5 + m_6*p**6 + m_7*p**7 == s*d1

    B = (u+1)**2 + 4
    A = u**3*B + 56
    m0 == 2*u**3*A + 55*u**2*B
    m1 == -4*u**2*A - 75*u*B
    m2 == -2*u*A - 125*B
    m3 == -u**4*A -24*u**3*B + 196
    m4 == u**3*A + 10*u**2*B
    m5 == 3*u**2*A + 100*u*B
    m6 == -11*u*A-250*B
    m7 == 7*A

    Cost: 37 S + 35 M + 2*exp(u+1) 7*exp(u)
    expected cost: 7*exp_u + 2*exp_u1 + 34*s16_cyclo + 32*m16 + 7*f16 + 3*c16_cyclo
    where c16_cyclo is a cube, hence 3*c16_cyclo = 3*(S_cyclo + M) = 3*S_cyclo + 3*M
    """
    t0 = m**2                               # t0  = m^2
    t1 = t0**2                              # t1  = m^4
    t2 = m**(u+1)                           # t2  = m^(u+1)
    t3 = t2**(u+1)                          # t3  = m^(u+1)^2
    t4 = t3 * t1                            # t4  = m^(4 + (u+1)^2)                = m^B
    t5 = t4**u                              # t5  = m^(u*B)
    t6 = (t4**2)**2 * t4                    # t6  = m^(5*B)
    t7 = ((t1**2)**2)**2                    # t7  = m^32
    t8 = t7**2                              # t8  = m^64
    inv_t1 = inversion_easy_gfp16(t1)       # m^(-4)
    t9 = t7 * inv_t1                        # t9  = m^(28)
    t10 = t9**2                             # t10 = m^56
    t11 = t5**u                             # t11 = m^(u^2*B)
    t12 = t11**u                            # t12 = m^(u^3*B)
    t01 = t12 * t10                         # t01 = m^(u^3*B + 56)                 = m^A
    t14 = t01**u                            # t14 = m^(u*A)
    inv_t14 = inversion_easy_gfp16(t14)     # m^(-u*A)
    t13 = inv_t14**2                        # t13 = m^(-2*u*A)
    t00 = (t6**2)**2 * t6                   # t00 = m^(25*B)
    t15 = (t00**2)**2 * t00                 # t15 = m^(125*B)
    inv_t15 = inversion_easy_gfp16(t15)     # m^(-125*B)
    t0 = t13 * inv_t15                      # t0 = m^(-2*u*A - 125*B)              = m^m2
    t16 = t0**2                             # t16 = m^(-4*u*A - 250*B)
    t17 = (t13**2)**2                       # t17 = m^(-8*u*A)
    t18 = t17 * t14                         # t18 = m^(-7*u*A)
    t2 = t16 * t18                          # t2  = m^(-11*u*A - 250*B)            = m^m6
    t19 = t14**u                            # t19 = m^(u^2*A)
    t20 = t19**u                            # t20 = m^(u^3*A)
    t21 = t20**u                            # t21 = m^(u^4*A)
    t22 = t19**2                            # t22 = m^(2*u^2*A)
    t23 = (t5**2)**2 * t5                   # t23 = m^(5*u*B)
    t24 = (t23**2)**2 * t23                 # t24 = m^(25*u*B)
    t25 = t24**2 * t24                      # t25 = m^(75*u*B)
    t26 = t24 * t25                         # t26 = m^(100*u*B)
    t27 = t22**2                            # t27 = m^(4*u^2*A)
    t37 = inversion_easy_gfp16(t27 * t25)   # t37 = m^(-4*u^2*A - 75*u*B) = m^m1
    t28 = t27 * inversion_easy_gfp16(t19)   # t28 = m^(3*u^2*A)
    t3 = t28 * t26                          # t3  = m^(3*u^2*A + 100*u*B)          = m^m5
    t29 = (t11**2)**2 * t11                 # t29 = m^(5*u^2*B)
    t30 = t29**2                            # t30 = m^(10*u^2*B)
    t4 = t20 * t30                          # t4  = m^(u^3*A + 10*u^2*B)           = m^m4
    s0 = t20**2                             # s0  = m^(2*u^3*A)
    s1 = (t30**2)**2 * t30                  # s1  = m^(50*u^2*B)
    s2 = s1 * t29                           # s2  = m^(55*u^2*B)
    s3 = s0 * s2                            # s3  = m^(2*u^3*A + 55*u^2*B)         = m^m0
    #t31 = t12**24    # t31 = m^(24*u^3*B)
    t31 = (((t12**2 * t12)**2)**2)**2       # t31 = m^(24*u^3*B)
    t5 = inversion_easy_gfp16(t21 * t31)    # m^(-u^4*A-24*u^3*B)
    t6 = (t8**2 * t8) * t1                  # t6  = m^196
    t7 = t5 * t6                            # t7  = m^(-u^4*A-24*u^3*B+196)        = m^m3
    t8 = (t01**2 * t01)**2 * t01            # m^(7*A)                             = m^m7
    t32 = t37.frobenius() * t7.frobenius(3) * t3.frobenius(5) * t8.frobenius(7)
    # t32 = m^(m1^p + m3*p^3 + m5*p^5 + m7*p^7)
    t33 = t0.frobenius(2) * t2.frobenius(6) # t33 = m^(m2*p^2 + m6*p^6)
    t = s3 * t32 * t33 * t4.frobenius(4)
    # m^(m0 + m1^p + m2*p^2 + m3*p^3 + m4*p^4 + m5*p^5 + m6*p^6 + m7*p^7)
    # B = (u+1)**2 + 4
    # A = u**3*B + 56
    # m0 = 2*u**3*A + 55*u**2*B
    # m1 = -4*u**2*A - 75*u*B
    # m2 = -2*u*A - 125*B
    # m3 = -u**4*A -24*u**3*B + 196
    # m4 = u**3*A + 10*u**2*B
    # m5 = 3*u**2*A + 100*u*B
    # m6 = -11*u*A-250*B
    # m7 = 7*A
    # #assert t4 == m**B
    # assert t01 == m**A
    # assert t0 == m**m2
    # assert t2 == m**m6
    # assert t37 == m**m1
    # assert t3 == m**m5
    # assert t4 == m**m4
    # assert s3 == m**m0
    # assert t7 == m**m3
    # assert t8 == m**m7
    return t

def final_exp_hard_kss16_v2(m, u):
    """
c=5
C5[1] =     -2*x^8 - 4*x^7 - 10*x^6 - 55*x^4 - 222*x^3 - 275*x^2 = (    2*x^3) * (-x^5 - 2*x^4 - 5*x^3 - 56) +  55 * (-x^4 - 2*x^3 - 5*x^2) = ... + (( -55*x^2) * (x^2 + 2*x + 5) + 0
C5[2] =        4*x^7 + 8*x^6 + 20*x^5 + 75*x^3 + 374*x^2 + 375*x = (   -4*x^2) * (-x^5 - 2*x^4 - 5*x^3 - 56) +  75 * (   x^3 + 2*x^2 + 5*x) = ... + ((    75*x) * (x^2 + 2*x + 5) + 0
C5[3] =           2*x^6 + 4*x^5 + 10*x^4 + 125*x^2 + 362*x + 625 = (     -2*x) * (-x^5 - 2*x^4 - 5*x^3 - 56) + 125 * (       x^2 + 2*x + 5) = ... + ((     125) * (x^2 + 2*x + 5) + 0
C5[4] =   x^9 + 2*x^8 + 5*x^7 + 24*x^5 + 104*x^4 + 120*x^3 - 196 = (-x^4 - 24) * (-x^5 - 2*x^4 - 5*x^3 - 56) + 125 * (             -308/25) = ... + ((       0) * (x^2 + 2*x + 5) + -1540
C5[5] =          -x^8 - 2*x^7 - 5*x^6 - 10*x^4 - 76*x^3 - 50*x^2 = (      x^3) * (-x^5 - 2*x^4 - 5*x^3 - 56) +  10 * (-x^4 - 2*x^3 - 5*x^2) = ... + (( -10*x^2) * (x^2 + 2*x + 5) + 0
C5[6] =      -3*x^7 - 6*x^6 - 15*x^5 - 100*x^3 - 368*x^2 - 500*x = (    3*x^2) * (-x^5 - 2*x^4 - 5*x^3 - 56) + 100 * (  -x^3 - 2*x^2 - 5*x) = ... + ((  -100*x) * (x^2 + 2*x + 5) + 0
C5[7] =       11*x^6 + 22*x^5 + 55*x^4 + 250*x^2 + 1116*x + 1250 = (    -11*x) * (-x^5 - 2*x^4 - 5*x^3 - 56) + 250 * (       x^2 + 2*x + 5) = ... + ((     250) * (x^2 + 2*x + 5) + 0
C5[8] =                           -7*x^5 - 14*x^4 - 35*x^3 - 392 = (        7) * (-x^5 - 2*x^4 - 5*x^3 - 56) + 250 * (                   0) = ... + ((       0) * (x^2 + 2*x + 5) + 0

(7*p^7 - 11*p^6*z + 3*p^5*z^2 + p^4*z^3 - p^3*z^4 - 24*p^3 - 2*p^2*z - 4*p*z^2 + 2*z^3) * (-x^5 - 2*x^4 - 5*x^3 - 56) + (250*p^6 - 100*p^5*z - 10*p^4*z^2 + 125*p^2 + 75*p*z - 55*z^2) * (x^2 + 2*x + 5) + (-1540*p^3)

    9*exp(u) + 34*S + 34*M + 13 Frobenius + 10 conjugations
    """
    p = m.parent().characteristic()
    # 1. (x^2 + 2*x + 5) = x^2 + 2*(x+2) + 1
    mu = m**u                              # exp(u)
    mu2 = mu**u                            # exp(u)
    m2 = m**2
    m4 = m2**2
    mc = mu2 * mu**2 * m4 * m              # 3M + 3S

    # 2.  -(x^5 + 2*x^4 + 5*x^3 + 56) = -(x^3*(x^2 + 2*x + 5) + 56)
    # 56 = 2*28 and 28 = 32 - 4
    m8 = m4**2
    m32 = (m8**2)**2                       # 2S
    m28 = m32 * inversion_easy_gfp16(m4)   # M
    m56 = m28**2                           # S

    mcu = mc**u
    mcu2 = mcu**u
    md = mcu2**u * m56                     # 3*exp(u) + M
    md = inversion_easy_gfp16(md)
    #assert md == m**(-(u**5 + 2*u**4 + 5*u**3 + 56))

    # 3. 7*p^7 - 11*p^6*x + 3*p^5*x^2 + p^4*x^3 - p^3*x^4 - 24*p^3 - 2*p^2*x - 4*p*x^2 + 2*x^3
    #    p^4 * (7*p^3 - 11*p^2*x + 3*p*x^2 + x^3) - p^3*x^4 - 24*p^3 - 2*p^2*x - 4*p*x^2 + 2*x^3
    mdu = md**u
    mdu2 = mdu**u
    mdu3 = mdu2**u
    mdu4 = mdu3**u
    mdu2_2= mdu2**2
    # 11  = 8 + 4 - 1 = 8+2+1 = 1011
    mdu_11= ((mdu**2)**2 * mdu)**2 * mdu
    md4 = (md**2)**2
    md8 = md4**2
    mm = (mdu3 * inversion_easy_gfp16((mdu2_2).frobenius() * mdu.frobenius(2) * (md8 * md4).frobenius(3)))**2 * inversion_easy_gfp16(mdu4.frobenius(3))
    mn = mdu3 * (mdu2_2 * mdu2).frobenius() * inversion_easy_gfp16((mdu_11).frobenius(2)) * (md8 * inversion_easy_gfp16(md)).frobenius(3)
    mf = mm * mn.frobenius(4)
    #assert mf == md ** (p**4 * (7*p**3 - 11*p**2*u + 3*p*u**2 + u**3) - p**3*u**4 - 24*p**3 - 2*p**2*u - 4*p*u**2 + 2*u**3)

    # 4. (250*p^6 - 100*p^5*x - 10*p^4*x^2 + 125*p^2 + 75*p*x - 55*x^2)
    #= 5*(50*p^6 - 20*p^5*x - 2*p^4*x^2 + 25*p^2 + 15*p*x - 11*x^2)
    #= 5*(2*(25*p^2 - 10*p*x - x^2)*p^4 + 25*p^2 + 15*p*x - 11*x^2)
    #(250*p^6 - 100*p^5*z - 10*p^4*z^2 + 125*p^2 + 75*p*z - 55*z^2) eq 5*(2*(25*p^2 - 10*p*z - z^2)*p^4 + 25*p^2 + 15*p*z - 11*z^2)

    # 11 = 1011 in binary
    mcu2_11 = ((mcu2**2)**2 * mcu2)**2 * mcu2 # 3 S + 2M
    # 5*(5*p^2 + 3*p*x)
    # 5 = 101 in bits 5 = 4+1 + 4+2-1
    # 25 = 11001
    mc_25 = (((mc**2 * mc)**2)**2)**2 * mc
    # 15 = 16-1 = 1111 = 10000 - 1
    # 10 = 8+2 = 1010
    mcu_2 = mcu**2
    mcu_8 = (mcu_2**2)**2
    mcu_10 = mcu_8 * mcu_2
    mcu_15 = mcu_8**2 * inversion_easy_gfp16(mcu)

    mc_25_p2 = mc_25.frobenius(2)
    me0 = (mc_25_p2 * mcu_15.frobenius() * inversion_easy_gfp16(mcu2_11))
    me1 = ((mc_25_p2 * inversion_easy_gfp16(mcu_10.frobenius() * mcu2))**2).frobenius(4)
    #assert me0 == mc**(25*p**2 + 15*p*u - 11*u**2)
    #assert me1 == mc**(2*(25*p**2 - 10*p*u - u**2)*p**4)
    me = me1 * me0
    # me^5
    me = (me**2)**2 * me
    #assert me == mc**(250*p**6 - 100*p**5*u - 10*p**4*u**2 + 125*p**2 + 75*p*u - 55*u**2)

    # 1540 = 28 * 55
    # 56 = 2*28
    # 55 = 64 - 8 - 1
    m28_8 = (m56**2)**2
    m28_64 = ((m28_8**2)**2)**2
    m1540 = inversion_easy_gfp16(m28_64) * m28_8 * m28

    return me * mf * m1540.frobenius(3)
