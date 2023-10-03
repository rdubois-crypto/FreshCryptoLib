from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.arith.misc import XGCD, xgcd
from sage.arith.misc import valuation

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from external.Pairings.pairing import *
from external.Pairings.pairing_bn import *
from external.Pairings.tests.test_pairing import *
from external.Pairings.tests.test_scalar_mult import test_glv_scalar_mult_g1
from external.Pairings.tests.test_pairing_bls12 import test_G2_endomorphism
from external.Pairings.tests.test_pairing_bw6_761 import test_g2_frobenius_eigenvalue_bw6
from external.Pairings.tests.test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761, test_bilinear_miller_loop_opt_ate_bw6_761_aklgl
from external.Pairings.tests.test_pairing_cp12_bls12 import test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar

from cost_pairing import cost_pairing_bn

def test_untwist_Frobenius_twist(E2, Fqd, map_Fqd_Fpk, map_Fpk_Fqd, r, c2, D_twist=False):
    """Test pi(Q) on G2 with the compressed form over Fq of Q

    INPUT:
    -`E2`: the d-twist over Fq = Fp^(k/d), Fp2 for BN, BLS12 curves, Fp4 for BLS24 curves
    -`Fqd`: with the explicit extension of degree d
    -`map_Fqd_Fpk`: map for the change of variables from Fqd to Fpk (isomorphism)
    -`map_Fpk_Fqd`: map for the change of variables from Fpk to Fqd (isomorphism)
    -`r`: prime integer
    -`c2`: twist cofactor (for G2)
    -`D_twist`: is E2 a D-twist or M-twist of E
    """
    p = Fqd.characteristic()
    w = Fqd.gen(0)
    Fq = E2.base_field()
    E20 = E2(0)
    Q2 = c2*E2.random_element()
    while Q2 == E20 or r*Q2 != E20:
        Q2 = c2 * E2.random_element()
    # now apply the twist, then the map from Fqd to Fpk
    Q2_ = (Fqd(Q2[0]), Fqd(Q2[1]))
    if D_twist:
        Q = psi_sextic_d_twist(Q2_, w)
    else:
        Q = psi_sextic_m_twist(Q2_, w)
    Q = (map_Fqd_Fpk(Q[0]), map_Fqd_Fpk(Q[1]))
    pQ = (map_Fpk_Fqd(Q[0].frobenius()), map_Fpk_Fqd(Q[1].frobenius()))
    # untwist: reverse map of the twist
    if D_twist:
        piQ = psi_sextic_m_twist(pQ, w)
    else:
        piQ = psi_sextic_d_twist(pQ, w)
    try:
        piQ = (Fq(piQ[0].list()[0]), Fq(piQ[1].list()[0]))
    except AssertionError:
        print("piQ[0] = {}\npiQ[1] = {}\n".format(piQ[0], piQ[1]))
        return False
    piQ2 = E2(piQ)
    p_mod_r_Q2 = (p % r)*Q2
    ok = piQ2 == p_mod_r_Q2
    print("test_untwist_Frobenius_twist(Q) = (p mod r)*Q on G2: {}".format(ok))
    # now test if directly on G2 is working too
    xi = -Fqd.modulus().constant_coefficient() # works with absolute extension and tower extension
    assert w**6 == xi
    xi_p16 = xi**((p-1)//6)
    if D_twist:
        xi_p13 = xi_p16**2
        xi_p12 = xi_p13 * xi_p16
        pQ2 = ((Q2[0]).frobenius() * xi_p13, (Q2[1]).frobenius() * xi_p12)
    else:
        xip16 = 1/xi_p16
        xip13 = xip16**2
        xip12 = xip13 * xip16
        pQ2 = ((Q2[0]).frobenius() * xip13, (Q2[1]).frobenius() * xip12)
    pQ2 = E2(pQ2)
    ok = pQ2 == p_mod_r_Q2
    print("test (xQ^p*xi^(p-1)/3, yQ^p*xi^(p-1)/2) = (p mod r)*Q on G2: {}".format(ok))
    if not ok:
        print("(p mod r)*Q2 = {}".format(p_mod_r_Q2))
        print("map(Q2) = {}".format(pQ2))
    return ok

def test_frobenius_in_Fq6(Fq6, D_twist=False):
    """
    Testing manual Frobenius computations in Fq6
    """
    p = Fq6.characteristic()
    p2 = p**2
    p3 = p*p2
    p4 = p2**2
    p5 = p*p4
    p6 = p3**2
    w0 = Fq6.gen(0)
    Fq = Fq6.base_field()
    xi = -Fq6.modulus().constant_coefficient()
    omega1 = xi**((p-1)//6)
    omega2 = omega1**2
    omega3 = omega1 * omega2
    omega4 = omega2**2
    omega5 = omega4 * omega1

    xip1 = xi.frobenius() * xi # xi^(p+1)
    omega1_p2 = xip1**((p-1)//6)
    omega2_p2 = omega1_p2 - 1
    omega3_p2 = Fq(-1)
    omega4_p2 = -omega1_p2
    omega5_p2 = -omega1_p2 + 1

    # omegai_p3 = omegai_p2 * omegai
    omega1_p3 = omega1_p2 * omega1
    omega2_p3 = omega2_p2 * omega2
    omega3_p3 = - omega3
    omega4_p3 = omega4_p2 * omega4
    omega5_p3 = omega5_p2 * omega5

    # omegai_p4 = omegai_p2^2
    omega1_p4 = omega2_p2
    omega2_p4 = omega4_p2
    omega3_p4 = 1
    omega4_p4 = omega2_p2
    omega5_p4 = omega4_p2

    # omegai_p5 = omegai_p4 * omegai_p
    omega1_p5 = omega1_p4 * omega1
    omega2_p5 = omega2_p4 * omega2
    omega3_p5 = omega3
    omega4_p5 = omega4_p4 * omega4
    omega5_p5 = omega5_p4 * omega5

    assert omega1_p2 == xi**((p2-1)//6)
    assert omega2_p2 == omega1_p2**2
    assert omega3_p2 == omega1_p2**3
    assert omega4_p2 == omega1_p2**4
    assert omega5_p2 == omega1_p2**5

    assert omega1_p3 == xi**((p3-1)//6)
    assert omega2_p3 == omega1_p3**2
    assert omega3_p3 == omega1_p3**3
    assert omega4_p3 == omega1_p3**4
    assert omega5_p3 == omega1_p3**5

    assert omega1_p4 == xi**((p4-1)//6)
    assert omega2_p4 == omega1_p4**2
    assert omega3_p4 == omega1_p4**3
    assert omega4_p4 == omega1_p4**4
    assert omega5_p4 == omega1_p4**5

    assert omega1_p5 == xi**((p5-1)//6)
    assert omega2_p5 == omega1_p5**2
    assert omega3_p5 == omega1_p5**3
    assert omega4_p5 == omega1_p5**4
    assert omega5_p5 == omega1_p5**5

    w1 = w0**p
    w2 = w0**p2
    w3 = w0**p3
    w4 = w0**p4
    w5 = w0**p5
    w6 = w0**p6
    print("w^p     == w*omega1:    {}".format(w1 == w0 * omega1))
    print("w^(p^2) == w*omega1_p2: {}".format(w2 == w0 * omega1_p2))
    print("w^(p^3) == w*omega1_p3: {}".format(w3 == w0 * omega1_p3))
    print("w^(p^4) == w*omega1_p4: {}".format(w4 == w0 * omega1_p4))
    print("w^(p^5) == w*omega1_p5: {}".format(w5 == w0 * omega1_p5))
    print("w^(p^6) == -w:          {}".format(w6 == -w0))
    i = 0
    ok = True
    while ok and i < 10:
        m = Fq6.random_element()
        mp1 = m**p
        mp2 = m**p2
        mp3 = m**p3
        mp4 = m**p4
        mp5 = m**p5
        mp6 = m**p6

        cm = m.list()
        m1 = cm[:]
        m1[0] = m1[0].frobenius()
        m1[1] = m1[1].frobenius() * omega1
        m1[2] = m1[2].frobenius() * omega2
        m1[3] = m1[3].frobenius() * omega3
        m1[4] = m1[4].frobenius() * omega4
        m1[5] = m1[5].frobenius() * omega5
        m1 = Fq6(m1)
        ok1 = mp1 == m1

        m2 = cm[:]
        m2[1] = m2[1]*omega1_p2
        m2[2] = m2[2]*omega2_p2
        m2[3] = -m2[3] #m2[3] = m2[3]*omega3_p2
        m2[4] = m2[4]*omega4_p2
        m2[5] = m2[5]*omega5_p2
        m2 = Fq6(m2)
        ok2 = mp2 == m2

        m3 = cm[:]
        m3[0] = m3[0].frobenius()
        m3[1] = m3[1].frobenius() * omega1_p3
        m3[2] = m3[2].frobenius() * omega2_p3
        m3[3] = m3[3].frobenius() * omega3_p3
        m3[4] = m3[4].frobenius() * omega4_p3
        m3[5] = m3[5].frobenius() * omega5_p3
        m3 = Fq6(m3)
        ok3 = mp3 == m3

        m4 = cm[:]
        m4[1] = m4[1]*omega1_p4
        m4[2] = m4[2]*omega2_p4
        # m4[3] does not change
        m4[4] = m4[4]*omega4_p4
        m4[5] = m4[5]*omega5_p4
        m4 = Fq6(m4)
        ok4 = mp4 == m4

        m5 = cm[:]
        m5[0] = m5[0].frobenius()
        m5[1] = m5[1].frobenius() * omega1_p5
        m5[2] = m5[2].frobenius() * omega2_p5
        m5[3] = m5[3].frobenius() * omega3_p5
        m5[4] = m5[4].frobenius() * omega4_p5
        m5[5] = m5[5].frobenius() * omega5_p5
        m5 = Fq6(m5)
        ok5 = mp5 == m5

        m6 = cm[:]
        m6[1] = -m6[1]
        m6[3] = -m6[3]
        m6[5] = -m6[5]
        m6 = Fq6(m6)
        ok6 = mp6 == m6
        ok = ok1 and ok2 and ok3 and ok4 and ok5 and ok6
        i = i+1
    if not D_twist:
        print("test_frobenius_in_Fq6M: {}".format(ok))
    else:
        print("test_frobenius_in_Fq6D: {}".format(ok))
    return ok

def test_formulas_tate_pairing(E, E2, r, c, c2, u):
    """
    Test that [a0]P + [a1][p^2 mod r]P = 0 on G1 = E(Fp)
    With a0 = 2*u+1, a1 = 6*u^2 + 2*u
    And  a0 = 6*u^2+4*u+1, a1 = -2*u-1
    """
    p = E.base_field().characteristic()
    T = (p % r)
    T2 = T**2 % r
    ok = True
    okG1a = True
    okG2a = True
    okG1b = True
    okG2b = True
    a0 = 2*u+1; a1 = 6*u**2 + 2*u
    b0 = 6*u**2+4*u+1; b1 = -2*u-1
    i = 0
    while ok and i < 10:
        P = c*E.random_element()
        while P == E(0) or r*P != E(0):
            P = c*E.random_element()
        Q = c2*E2.random_element()
        while  Q == E2(0) or r*Q != E2(0):
            Q = c2*E2.random_element()
        T2P = T2*P
        T2Q = T2*Q
        okG1a = a0*P + a1*T2P == E(0)
        okG2a = a0*Q + a1*T2Q == E2(0)
        okG1b = b0*P + b1*T2P == E(0)
        okG2b = b0*Q + b1*T2Q == E2(0)
        ok = okG1a and okG2a and okG1b and okG2b
        i = i+1
    print("test_formulas_tate_pairing: {} (2*u+1)+(6*u^2+2*u)*(p^2%r) = 0, (6*u^2+4*u+1)-(2*u+1)*(p^2%r) = 0".format(ok))
    return ok

def test_bilinear_miller_loop_with_frobenius(miller_loop_function, E, E_Fpk, E2, Fqd, map_Fqd_Fpk, r, c, c2, u0, D_twist=False, Tate=False):
    """Test the bilinearity of the arg function

    INPUT:
    -`miller_loop_function`: a function
    -`E`: elliptic curve over prime field Fp of order c*r (subgroup of order r is G1)
    -`E_Fpk`: E(Fpk) (just a change in the coordinate domain)
    -`E2`: the d-twist over Fq = Fp^(k/d), Fp2 for BN, BLS12 curves, Fp4 for BLS24 curves
    -`Fqd`: with the explicit extension of degree d
    -`map_Fqd_Fpk`: map for the change of variables from Fqd to Fpk (isomorphism)
    -`r`: prime integer
    -`c`: curve cofactor (for G1)
    -`c2`: twist cofactor (for G2)
    -`u0`: seed
    -`D_twist`: is E2 a D-twist or M-twist of E

    Functions:
    miller_loop_opt_ate_bn
    miller_loop_opt_ate_bn_2naf
    miller_loop_opt_tate_bn (Tate=True)
    """
    E0 = E(0)
    E20 = E2(0)
    Fpk = E_Fpk.base_field()
    w = Fqd.gen(0)
    P = c * E.random_element()
    while P == E0 or r*P != E0:
        P = c * E.random_element()
    Q2 = c2*E2.random_element()
    while Q2 == E20 or r*Q2 != E20:
        Q2 = c2 * E2.random_element()
    # now apply the twist, then the map from Fqd to Fpk
    Q2 = (Fqd(Q2[0]), Fqd(Q2[1]))
    if D_twist:
        Q = psi_sextic_d_twist(Q2, w)
    else:
        Q = psi_sextic_m_twist(Q2, w)
    Q = E_Fpk((map_Fqd_Fpk(Q[0]), map_Fqd_Fpk(Q[1])))
    if Tate:
        f = miller_loop_function(P, Q, u0)
    else:
        f = miller_loop_function(Q, P, u0)
    assert (Fpk.cardinality()-1) % r == 0
    exponent = (Fpk.cardinality()-1) // r
    assert exponent % r != 0
    g = f**exponent
    assert g != 1 and g**r == 1
    ok = True
    aa = 1
    while ok and aa < 4:
        bb = 1
        while ok and bb < 4:
            if Tate:
                fij = miller_loop_function(aa*P, bb*Q, u0)
            else:
                fij = miller_loop_function(bb*Q, aa*P, u0)
            gij = fij**exponent
            ok1 = gij != 1 and gij**r == 1
            if not ok1:
                print("error gij == 1: {}, gij^r == 1: {}".format(gij != 1, gij**r == 1))
            ok = gij == g**(aa*bb)
            if not ok:
                print("error gij ! g^(a*b), a={}, b={}".format(aa,bb))
                print("gij = {}\n1/gij = {}".format(gij, 1/gij))
                print("g = {}\ng^(a*b) = {}".format(g, g**(aa*bb)))
            bb += 1
        aa += 1
    print("test bilinearity {}: {}".format(miller_loop_function.__name__, ok))
    return ok

def test_miller_loop_opt_tate_a0_twist6_aklgl_2_multi_scalar(E, E2, E_Fqd, E_Fpk, u, p, r, tr, c, c2, a0, a1, omega, map_Fqd_Fpk, D_twist=True, function_name=miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar):
    # test miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar
    P = c*E.random_element()
    while P == E(0):
        P = c * E.random_element()
    Q = c2*E2.random_element()
    while Q == E2(0):
        Q = c2 * E2.random_element()
    Fqd = E_Fqd.base_field() # assume that the degree d extension is explicit
    Fpk = E_Fpk.base_field()
    Fp = Fpk.base_ring()
    w = Fqd.gen(0)
    if Fqd.base_ring().degree() == 1:
        xi = -Fqd.polynomial().constant_coefficient()
    else:
        xi = -Fqd.modulus().constant_coefficient()
    exponent = (Fqd.cardinality()-1) // r # (q^d-1)//r = (p^k-1)//r
    omega = Fp(omega)
    assert omega**2 + omega + 1 == 0
    psiP = E([omega * P[0], P[1]])
    k = Fpk.degree()
    Fq = E2.base_field()
    e = Fq.degree()
    d = k//e
    tre_mod_r = (tr-1)**e % r
    #assert ((r - a0) % a1) == 0
    #lambd = (r - a0) // a1
    # there is a co-factor...
    if ((a0 + a1*tre_mod_r) % r) == 0:
        lambd = tre_mod_r
    elif ((a0 - a1*tre_mod_r) % r) == 0:
        lambd = -tre_mod_r
    elif ((a0 + a1*(tre_mod_r-1)) % r) == 0:
        lambd = tre_mod_r-1
    elif ((a0 - a1*(tre_mod_r-1)) % r) == 0:
        lambd = -tre_mod_r+1
    else:
        print("eigenvalue lambda such that a0 + a1*lambda = 0 mod r not found")
    neg_psi = False
    if ((lambd**2 + lambd + 1) % r) != 0 and ((lambd**2 - lambd + 1) % r) == 0:
        psiP = -psiP
        neg_psi = True
    if not (psiP == lambd*P or psiP == -lambd*P):
        print("adjusting omega -> -omega-1")
        omega = -omega - 1
        psiP = E([omega * P[0], P[1]])
        assert (psiP == lambd*P or psiP == -lambd*P)

    ok = True
    if function_name is miller_loop_opt_tate_bn_aklgl_a0_2_parallel:
        f, scalarsP = function_name(P, Q, u, E.a6(), Fqd, map_Fqd_Fpk, D_twist, omega, xi)
    elif function_name is miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar:
        f, scalarsP = function_name(P, Q, E.a6(), u, omega, xi, Fqd, D_twist)
    else:
        f, scalarsP = function_name(P, Q, E.a6(), a0, a1, omega, xi, Fqd, D_twist)

    ok = f != Fpk(1)
    if not ok:
        print("test bilinearity {}: {} pairing result is 1".format(function_name.__name__, ok))
        return False
    #print("type(scalarsP) = {}".format(type(scalarsP)))
    check_ai_P = ((type(scalarsP) is list) or (type(scalarsP) is tuple)) and len(scalarsP) == 2
    if check_ai_P:
        a0P, a1P = scalarsP
        a0P = (a0P[0]/a0P[2], a0P[1]/a0P[2])
        a1P = (a1P[0]/a1P[2], a1P[1]/a1P[2])
        ok_a0P = E(a0P) == a0*P
        ok_a1P = E(a1P) == a1*P
    else:
        ok_ai_P = scalarsP[2] == 0
    g = f**exponent
    aa = 1
    bb = 1
    while ok and aa < 4:
        Pa = aa*P
        psiPa = E([Pa[0]*omega, Pa[1]])
        if neg_psi:
            psiPa = -psiPa
        bb = 1
        while ok and bb < 4:
            Qb = bb*Q
            if function_name is miller_loop_opt_tate_bn_aklgl_a0_2_parallel:
                fab, scalars_abP = function_name(Pa, Qb, u, E.a6(), Fqd, map_Fqd_Fpk, D_twist, omega, xi)
            elif function_name is miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar:
                fab, scalars_abP = function_name(Pa, Qb, E.a6(), u, omega, xi, Fqd, D_twist)
            else:
                fab, scalars_abP = function_name(Pa, Qb, E.a6(), a0, a1, omega, xi, Fqd, D_twist)
            gab = fab**exponent
            ok = ok and gab != Fpk(1) and gab == g**(aa*bb)
            if check_ai_P:
                a0P, a1P = scalars_abP
                a0P = (a0P[0]/a0P[2], a0P[1]/a0P[2])
                a1P = (a1P[0]/a1P[2], a1P[1]/a1P[2])
                ok_a0P = ok_a0P and E(a0P) == a0*Pa
                ok_a1P = ok_a1P and E(a1P) == a1*Pa
            else:
                ok_ai_P = ok_ai_P and scalars_abP[2] == 0
            bb += 1
        aa += 1
    print("test bilinearity {} scalars {},{}: {} ({} tests)".format(function_name.__name__, a0, a1, ok, (aa-1)*(bb-1)))
    if check_ai_P:
        print("a0P == [a0]P, a1P == [a1](P): {}".format(ok_a0P and ok_a1P))
    else:
        print("a0*P + a1*psiP == 0: {}".format(ok_ai_P))
    return ok

def test_final_exp_hard_bn(Fpk, u, r, exponent_hard=None):
    """
    test the hard part of the final exponentiation
    """
    p = Fpk.characteristic()
    ee = (p**4 - p**2 + 1)//r
    if exponent_hard is None:
        hh = 12*u**3 + 6*u**2 + 2*u
        exponent_hard = hh * ee
    else:
        assert (exponent_hard % ee) == 0
    ok = True
    i = 0
    while ok and i < 10:
        f = Fpk.random_element()
        g = final_exp_easy_k12(f)
        h = final_exp_hard_bn(g, u)
        ok = g**exponent_hard == h
        i += 1
    print("test final_exp_hard_bn: {}".format(ok))
    return ok

def test_curve(test_vector_value):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # BN polynomials
    px = 36*x**4 + 36*x**3 + 24*x**2 + 6*x + 1
    tx = 6*x**2 + 1
    rx = 36*x**4 + 36*x**3 + 18*x**2 + 6*x + 1
    cx = QQx(1)
    yx = 6*x**2 + 4*x + 1
    assert px == (tx**2 + 3*yx**2)//4
    assert rx == ((tx-2)**2 + 3*yx**2)//4
    # lambx
    # 36*x^3 + 18*x^2 + 6*x + 1
    # -36*x^3 - 18*x^2 - 6*x - 2
    lambx = 36*x**3 + 18*x**2 + 6*x + 1 # eigenvalue mod r, s.t. lambx^2 + lambx + 1 = 0
    lambx_ = -36*x**3 - 18*x**2 - 6*x - 2
    assert ((lambx**2 + lambx + 1) % rx) == 0
    assert ((lambx_**2 + lambx_ + 1) % rx) == 0
    # -18*x^3 - 18*x^2 - 9*x - 2
    # 18*x^3 + 18*x^2 + 9*x + 1
    betax = 18*x**3 + 18*x**2 + 9*x + 1 # mod p, s.t. betax^2 + betax + 1 = 0
    betax_ = -18*x**3 - 18*x**2 - 9*x - 2
    assert ((betax**2 + betax + 1) % px) == 0
    assert ((betax_**2 + betax_ + 1) % px) == 0
    # cofactor for G2
    t2x = tx**2 - 2*px
    y2x = tx*yx
    # (t2x**2 - 4*p2x) == -3*(tx*yx)**2
    tw2x = (t2x + 3*y2x)//2
    #tw2x_ = (t2x - 3*y2x)//2
    #assert ((px2 + 1 - tw2x) % rx) == 0
    #c2x = (px2 + 1 - tw2x) // rx
    c2x = 36*x**4 + 36*x**3 + 30*x**2 + 6*x + 1
    k = 12
    # hard part of the final exponentiation
    exponent_hard = (px**4 - px**2 + 1) // rx
    ex = exponent_hard
    l0 = 1+6*x + 12*x**2 + 12*x**3
    l1 = 4*x + 6*x**2 + 12*x**3
    l2 = 6*x + 6*x**2 + 12*x**3
    l3 = -1 + 4*x + 6*x**2 + 12*x**3
    assert ((l0 + l1*px + l2*px**2 + l3*px**3) % ex) == 0
    v = test_vector_value
    u0 = ZZ(v['u'])
    u = u0
    # check formula for final exponentiation
    assert ((x**12-1) // cyclotomic_polynomial(12)) == (x**6-1)*(x**2+1)
    exponent_x = (px**4-px**2+1)//rx
    exponent = ZZ(exponent_x(u0))
    exponent_easy = (px**6-1)*(px**2+1)
    exponent_hard = (px**4-px**2+1)//rx
    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    c2 = ZZ(c2x(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    b, E = find_curve_parameter_b(Fp, r, c)
    print("BN-{} curve u={:#x}".format(p.nbits(), u0))
    print("curve parameter b = {}".format(b))
    print("p = {} mod 4".format(p % 4))
    print("r is prime: {}".format(r.is_prime()))
    print("r-1 has 2-adicity {}".format((r-1).valuation(2)))
    print("p-1 has 2-adicity {}".format((p-1).valuation(2)))
    print("exponent_easy = (px**6-1)*(px**2+1)")
    print("exponent_hard = (px**4-px**2+1)//rx = {}".format(exponent_hard))
    print("(p^4-p^2+1)//r:")
    print("exponent = {}".format(exponent))
    print("exponent is prime: {}\n".format(exponent.is_prime()))

    cost_pairing_bn(u)

    lambda_mod_r = ZZ(lambx(u0))
    beta_mod_p = Fp(betax(u0))
    omega = beta_mod_p
    omega_ = -omega-1
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    if (p % 4) == 3:
        Fp2 = Fp.extension(z**2 + 1, names=('i',)); (i,) = Fp2._first_ngens(1)
        a = -1
        print("Fp2 = Fp[x]/(x^2 + 1)")
    else:
        a = 2
        while not (z**2 - a).is_irreducible():
            a = a+1
        print("Fp2 = Fp[x]/(x^2 - {})".format(a))
        Fp2 = Fp.extension(z**2 - a, names=('i',)); (i,) = Fp2._first_ngens(1)
    Fq = Fp2
    Fp2s = Fp2['s']; (s,) = Fp2s._first_ngens(1)
    xiM, btwM = find_twist_curve_parameter_xi_ab(b, Fp2, r, c2, d=6, D_twist=False)
    EM = EllipticCurve([Fp2(0), Fp2(btwM)])
    Fq6M = Fp2.extension(s**6 - xiM, names=('wM',)); (wM,) = Fq6M._first_ngens(1)
    E_Fq6M = EllipticCurve([Fq6M(0), Fq6M(b)])

    # note: s^6 = i+1 => (s^6-1)^2 = i^2 = -1 so s^12 - 2*s^6 + 2 = 0
    try:
        coeffs_xiM = xiM.polynomial().list()
    except AttributeError as err:
        coeffs_xiM = xiM.list()
    i0M = coeffs_xiM[0]
    i1M = coeffs_xiM[1]
    i0m = ZZ(i0M)
    if abs(i0m - p) < abs(i0m):
        i0m = i0m - p
    i1m = ZZ(i1M)
    if abs(i1m - p) < abs(i1m):
        i1m = i1m - p
    if i0m == 0:
        str_xiM = ""
    else:
        str_xiM = "{}".format(i0m)
    if i1m == 1:
        if len(str_xiM) == 0:
            str_xiM = "i"
        else:
            str_xiM += "+i"
    elif i1m == -1:
        str_xiM += "-i"
    elif i1m != 0:
        if len(str_xiM) == 0:
            str_xiM = "{}*i".format(i1m)
        else:
            str_xiM += "{:+}*i".format(i1m)
    print("M-twist xiM = {}".format(str_xiM))
    # xi = i0 + i*i1
    # s^6 - xi = 0 <=> s^6 - i0 = i1*i <=> (s^6 - i0)^2 = i1^2*a
    # resultant(s^6-xi, z^2 - a)
    Fp12M = Fp.extension((z**6 - i0M)**2 - i1M**2*a, names=('SM',)); (SM,) = Fp12M._first_ngens(1)
    E_Fp12M = EllipticCurve([Fp12M(0), Fp12M(b)])
    # (z^6 - i0M)^2 - i1M^2*a = 0
    # z^12 -2*i0M*z^6 + i0M^2 - i1M^2*a = 0

    #def map_Fq6M_Fp12M(x):
        # evaluate elements of Fq6M = Fp[i]/(i^2+1)[s]/(s^6-(i+1)) at i=S^6-1 and s=S
        # i <-> wM^6-1 = SM^6-1 and wM <-> SM
        #return sum([sum([yj*(SM**6-1)**j for j,yj in enumerate(xi.polynomial())]) * SM**i for i,xi in enumerate(x.list())])
        #return sum([xi.polynomial()((SM**6-i0M)/i1M) * SM**e for e,xi in enumerate(x.list())])
    def map_Fq6M_Fp12M(x, aM=None):
        if aM is None:
            aM = SM
        return sum([xi.polynomial()((aM**6-i0M)/i1M) * aM**e for e,xi in enumerate(x.list())])

    def map_Fp12M_Fq6M(x, w_Fq6=None):
        if w_Fq6 is None:
            w_Fq6 = wM
        return x.polynomial()(w_Fq6)

    def map_Fp2_Fp12M(x):
        # evaluate elements of Fq=Fp[i] at i=wM^6-1 = SM^6-1
        return x.polynomial()((SM**6-i0M)/i1M)

    xiD, btwD = find_twist_curve_parameter_xi_ab(b, Fp2, r, c2, d=6, D_twist=True)
    ED = EllipticCurve([Fp2(0), Fp2(btwD)])
    Fq6D = Fp2.extension(s**6 - xiD, names=('wD',)); (wD,) = Fq6D._first_ngens(1)
    E_Fq6D = EllipticCurve([Fq6D(0), Fq6D(b)]) # but un Sage, an elliptic curve over an extension field in two layers is not handled easily
    # note: s^6 = i+2 <=> s^6-2 = i => s^12 -4*s^6 + 4 = i^2 = -1 so s^12 - 4*s^6 + 5 = 0
    try:
        coeffs_xiD = xiD.polynomial().list()
    except AttributeError as err:
        coeffs_xiD = xiD.list()
    i0D = coeffs_xiD[0]
    i1D = coeffs_xiD[1]
    i0d = ZZ(i0D)
    if abs(i0d - p) < abs(i0d):
        i0d = i0d - p
    i1d = ZZ(i1D)
    if abs(i1d - p) < abs(i1d):
        i1d = i1d - p
    if i0d == 0:
        str_xiD = ""
    else:
        str_xiD = "{}".format(i0d)
    if i1d == 1:
        if len(str_xiD) == 0:
            str_xiD = "i"
        else:
            str_xiD += "+i"
    elif i1d == -1:
        str_xiD += "-i"
    elif i1d != 0:
        if len(str_xiD) == 0:
            str_xiD = "{}*i".format(i1d)
        else:
            str_xiD += "{:+}*i".format(i1d)
    print("D-twist xiD = {}".format(str_xiD))
    # xiD = i0D + i*i1D
    # s^6 - xiD = 0 <=> s^6 - i0D = i1D*i <=> (s^6 - i0D)^2 = i1D^2*a
    # resultant(s^6-xiD, z^2 - a)
    Fp12D = Fp.extension((z**6 - i0D)**2 - i1D**2*a, names=('SD',)); (SD,) = Fp12D._first_ngens(1)
    E_Fp12D = EllipticCurve([Fp12D(0), Fp12D(b)])

    #def map_Fq6D_Fp12D(x):
        # evaluate elements of Fq6D = Fp[i]/(i^2+1)[s]/(s^6-(i+2)) at i=S^6-2 and s=S
        # i <-> s^6-2 = SD^6-2 and s <-> SD
        #return sum([sum([yj*(SD**6-2)**j for j,yj in enumerate(xi.polynomial())]) * SD**e for e,xi in enumerate(x.list())])
    #    return sum([xi.polynomial()((SD**6-i0D)/i1D) * SD**e for e,xi in enumerate(x.list())])
    def map_Fq6D_Fp12D(x, aD=None):
        if aD is None:
            aD = SD
        return sum([xi.polynomial()((aD**6-i0D)/i1D) * aD**e for e,xi in enumerate(x.list())])
    def map_Fp12D_Fq6D(x, w_Fq6=None):
        if w_Fq6 is None:
            w_Fq6 = wD
        return x.polynomial()(w_Fq6)

    def map_Fp2_Fp12D(x):
        # evaluate elements of Fq=Fp[i] at i=s^6-1 = S^6-1
        return x.polynomial()((SD**6-i0D)/i1D)

    print("test E (G1)")
    test_order(E,r*c)
    print("test E' (G2) M-twist")
    test_order(EM,r*c2)

    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(E_Fp12M,EM,Fq6M,map_Fq6M_Fp12M,r,c2,D_twist=False)
    test_g2_frobenius_eigenvalue_alt(E_Fp12M,EM,map_Fp2_Fp12M,r,c2,D_twist=False)
    print("test Frobenius map on G2 with D-twist")
    test_g2_frobenius_eigenvalue(E_Fp12D,ED,Fq6D,map_Fq6D_Fp12D,r,c2,D_twist=True)
    test_g2_frobenius_eigenvalue_alt(E_Fp12D,ED,map_Fp2_Fp12D,r,c2,D_twist=True)

    print("Test endomorphism on G2")
    test_G2_endomorphism(ED, E_Fq6D, D_twist=True)
    test_G2_endomorphism(EM, E_Fq6M, D_twist=False)

    print("test GLV on G1")
    test_glv_scalar_mult_g1(E, lambda_mod_r, beta_mod_p, r, c)

    test_frobenius_in_Fq6(Fq6M, D_twist=False)
    test_frobenius_in_Fq6(Fq6D, D_twist=True)

    print("test sparse-dense and sparse-sparse multiplications")
    test_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_d6_twist(Fq6D)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_sparse_mult_d6_twist(Fq6D)

    print("\ntest pairings with M-twist")

    test_final_exp_hard_bn(Fp12M, u, r)
    test_double_line_j(E,EM,Fq6M,D_twist=False)
    test_double_line_affine_j(E,EM,Fq6M,D_twist=False)
    test_add_line_j(E,EM,Fq6M,D_twist=False)
    test_add_line_affine_j(E,EM,Fq6M,D_twist=False)
    test_double_line_j_csb(E,EM,Fq6M,D_twist=False)
    test_add_line_j_csb(E,EM,Fq6M,D_twist=False)

    test_miller_function_ate(E, E_Fq6M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_2naf(E, E_Fq6M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_csb(E, E_Fq6M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_tate(E, E_Fq6M, EM, r, c, c2, D_twist=False)
    test_miller_function_tate_2naf(E, E_Fq6M, EM, r, c, c2, D_twist=False)

    test_untwist_Frobenius_twist(EM, Fq6M, map_Fq6M_Fp12M, map_Fp12M_Fq6M, r, c2, D_twist=False)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_ate_bn, E, E_Fp12M, EM, Fq6M, map_Fq6M_Fp12M, r, c, c2, u0, D_twist=False)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_ate_bn_2naf, E, E_Fp12M, EM, Fq6M, map_Fq6M_Fp12M, r, c, c2, u0, D_twist=False)

    test_formulas_tate_pairing(E, EM, r, c, c2, u0)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_tate_bn, E, E_Fp12M, EM, Fq6M, map_Fq6M_Fp12M, r, c, c2, u0, D_twist=False, Tate=True)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_tate_bn_alt, E, E_Fp12M, EM, Fq6M, map_Fq6M_Fp12M, r, c, c2, u0, D_twist=False, Tate=True)

    test_double_line_h_a0_twist6_aklgl(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_add_line_h_a0_twist6_aklgl_test(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl_with_z(E,EM,wM,D_twist=False)
    test_sparse_sparse_mult_d6_twist(Fq6M)
    test_sparse_mult_d6_twist(Fq6M)

    test_miller_function_ate_aklgl(E,EM,Fq6M,xiM,r,c,c2,t-1,D_twist=False)
    test_miller_function_ate_2naf_aklgl(E, EM, Fq6M, xiM, r, c, c2, t-1, D_twist=False)

    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bn_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bn_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)

    test_miller_loop_opt_tate_aklgl(miller_loop_opt_tate_bn_aklgl_a0, E, EM, Fq6M, map_Fq6M_Fp12M, r, c, c2, u, xiM, D_twist=False)
    a0 = 2*u+1
    a1 = 6*u**2+2*u
    test_miller_loop_opt_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fp12M, u, p, r, t, c, c2, a0, a1, omega, map_Fq6M_Fp12M, D_twist=False, function_name=miller_loop_opt_tate_bn_aklgl_a0_2_parallel)
    test_miller_loop_opt_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fp12M, u, p, r, t, c, c2, a0, a1, omega, map_Fq6M_Fp12M, D_twist=False, function_name=miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar)

    print("\ntest pairings with D-twist")

    test_final_exp_hard_bn(Fp12D, u, r)
    test_double_line_j(E,ED,Fq6D,D_twist=True)
    test_double_line_affine_j(E,ED,Fq6D,D_twist=True)
    test_add_line_j(E,ED,Fq6D,D_twist=True)
    test_add_line_affine_j(E,ED,Fq6D,D_twist=True)
    test_double_line_j_csb(E,ED,Fq6D,D_twist=True)
    test_add_line_j_csb(E,ED,Fq6D,D_twist=True)

    test_miller_function_ate(E, E_Fq6D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_2naf(E, E_Fq6D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_csb(E, E_Fq6D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_tate(E, E_Fq6D, ED, r, c, c2, D_twist=True)
    test_miller_function_tate_2naf(E, E_Fq6D, ED, r, c, c2, D_twist=True)

    test_untwist_Frobenius_twist(ED, Fq6D, map_Fq6D_Fp12D, map_Fp12D_Fq6D, r, c2, D_twist=True)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_ate_bn, E, E_Fp12D, ED, Fq6D, map_Fq6D_Fp12D, r, c, c2, u0, D_twist=True)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_ate_bn_2naf, E, E_Fp12D, ED, Fq6D, map_Fq6D_Fp12D, r, c, c2, u0, D_twist=True)

    test_formulas_tate_pairing(E, ED, r, c, c2, u0)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_tate_bn, E, E_Fp12D, ED, Fq6D, map_Fq6D_Fp12D, r, c, c2, u0, D_twist=True, Tate=True)
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_tate_bn_alt, E, E_Fp12D, ED, Fq6D, map_Fq6D_Fp12D, r, c, c2, u0, D_twist=True, Tate=True)

    test_double_line_h_a0_twist6_aklgl(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_add_line_h_a0_twist6_aklgl_test(E, ED, wD, D_twist=True)
    test_add_line_h_a0_twist6_aklgl_with_z(E,ED,wD,D_twist=True)
    test_sparse_sparse_mult_d6_twist(Fq6D)
    test_sparse_mult_d6_twist(Fq6D)

    test_miller_function_ate_aklgl(E,ED,Fq6D,xiD,r,c,c2,t-1,D_twist=True)
    
    print(" E, ED, Fq6D, xiD, r, c, c2, t-1",E, ED, Fq6D, xiD, r, c, c2, t-1);
    
    test_miller_function_ate_2naf_aklgl(E, ED, Fq6D, xiD, r, c, c2, t-1, D_twist=True)

    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bn_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bn_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)

    test_miller_loop_opt_tate_aklgl(miller_loop_opt_tate_bn_aklgl_a0, E, ED, Fq6D, map_Fq6D_Fp12D, r, c, c2, u, xiD, D_twist=True)
    test_miller_loop_opt_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fp12D, u, p, r, t, c, c2, a0, a1,  omega, map_Fq6D_Fp12D, D_twist=True, function_name=miller_loop_opt_tate_bn_aklgl_a0_2_parallel)
    test_miller_loop_opt_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fp12D, u, p, r, t, c, c2, a0, a1, omega, map_Fq6D_Fp12D, D_twist=True, function_name=miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar)

if __name__ == "__main__":
    arithmetic(False)
    test_vector_bn = [
 #       {'u':             0x36ab400000000000, 'b':  11, 'pnbits': 253,  'rnbits':253, 'label':"u=+2^62-2^59-2^56-2^54-2^52-2^50-2^48+2^46 Hw2NAF=8"},#twist-secure, BW6-509, BW6-510, BW6-512
 #       {'u':             0x381edc0000000000, 'b':  15, 'pnbits': 253,  'rnbits':253, 'label':"u=+2^62-2^59+2^53-2^48-2^45-2^42     Hw2NAF=6"},#twist-secure, BW6-512
 #       {'u':              2**62-2**54+2**44, 'b':   5, 'pnbits': 254,  'rnbits':254, 'label':"Nogami et al PAIRING'2008, TEPLA"},
  #      {'u':               -(2**62+2**55+1), 'b':   2, 'pnbits': 254,  'rnbits':254, 'label':"Pereira et al eprint 2010/429"},
        {'u':             0x44e992b44a6909f1, 'b':   3, 'pnbits': 254,  'rnbits':254, 'label':"Ethereum BN254, 2^28 | r-1"},
   #     {'u':2**62+2**59+2**55+2**15+2**10-1, 'b':   5, 'pnbits': 254,  'rnbits':254, 'label':"CBMNPZ15 eprint 2015/247"},
    #    {'u':                     1868033**3, 'b':   3, 'pnbits': 256,  'rnbits':256, 'label':"Naehrig-Niederhagen-Schwabe LATINCRYPT 2010 e2010/186"},
     #   {'u':     0x3e59c4000000000000000000, 'b':  10, 'pnbits': 382,  'rnbits':382, 'label':"u=+2^94-2^89+2^87-2^85-2^83+2^81-2^78+2^74 Hw2NAF=8"},#twist-secure, BW6-767
      #  {'u':     0x467eae000000000000000000, 'b':   7, 'pnbits': 382,  'rnbits':382, 'label':"u=+2^94+2^91-2^89+2^87-2^80-2^78-2^76-2^73 Hw2NAF=8"},#twist-secure, BW6-768
       # {'u':         -(2**94+2**76+2**72+1), 'b':   2, 'pnbits': 382,  'rnbits':382, 'label':"Pereira et al eprint 2010/429"},
        #{'u':                 2**110+2**36+1, 'b': 257, 'pnbits': 446,  'rnbits':446, 'label':"Pereira et al eprint 2010/429"},
        #{'u':-0x4000000000001000008780000000, 'b':  57, 'pnbits': 446,  'rnbits':446, 'label':"Pluto curve https://github.com/daira/pluto-eris/"},
    ]
    for tv in test_vector_bn:
        print(tv['label'])
        test_curve(tv)
