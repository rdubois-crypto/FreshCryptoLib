from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from external.Pairings.pairing import *
from external.Pairings.pairing_bw6_bls12 import *
from external.Pairings.tests.test_pairing import *

# test subgroup membership testing and cofactor multiplication
def test_bw6_761_phi(E,r,c):
    #print("test bw6_761_phi(P) (P in G1)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        assert r*c*P == E0
        phiP = bw6_761_phi(P)
        phi2P = bw6_761_phi(phiP)
        ok = phi2P + phiP + P == E0
        i = i+1
    print("test bw6_761_phi(P) (P in G1): {}".format(ok))
    return ok

def test_bw6_761_g1_mult_by_cofactor(E,r,c):
    #print("test bw6_761_g1_mult_by_cofactor(P)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        rP = r*P
        assert rP != E0 and c*rP == E0 and c*P != E0
        R = bw6_761_g1_mult_by_cofactor(P)
        ok = r*R == E0
        if ok:
            R = bw6_761_g1_mult_by_cofactor(rP)
            ok = R == E0
        i = i+1
    print("test bw6_761_g1_mult_by_cofactor(P): {}".format(ok))
    return ok

def test_bw6_761_g1_mult_by_cofactor_alt(E,r,c):
    #print("test bw6_761_g1_mult_by_cofactor_alt(P)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        rP = r*P
        assert rP != E0 and c*rP == E0 and c*P != E0
        R = bw6_761_g1_mult_by_cofactor_alt(P)
        ok = r*R == E0
        if ok:
            R = bw6_761_g1_mult_by_cofactor_alt(rP)
            ok = R == E0
        i = i+1
    print("test bw6_761_g1_mult_by_cofactor_alt(P): {}".format(ok))
    return ok

def test_bw6_761_g1_mult_by_r(E,r,c):
    #print("test bw6_761_g1_mult_by_r(P)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        cP = c*P
        assert cP != E0 and r*cP == E0 and r*P != E0
        R = bw6_761_g1_mult_by_r(P)
        ok = c*R == E0
        if ok:
            R = bw6_761_g1_mult_by_r(cP)
            ok = R == E0
        i = i+1
    print("test bw6_761_g1_mult_by_r(P): {}".format(ok))
    return ok

def test_bw6_761_g1_mult_by_r_alt(E,r,c):
    #print("test bw6_761_g1_mult_by_r_alt(P)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        cP = c*P
        assert cP != E0 and r*cP == E0 and r*P != E0
        R = bw6_761_g1_mult_by_r_alt(P)
        ok = c*R == E0
        if ok:
            R = bw6_761_g1_mult_by_r_alt(cP)
            ok = R == E0
        i = i+1
    print("test bw6_761_g1_mult_by_r_alt(P): {}".format(ok))
    return ok

def test_bw6_761_g1_check_membership(E,r,c):
    #print("test bw6_761_g1_check_membership(P)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        rP = r*P
        cP = c*P
        assert rP != E0 and cP != E0 and c*rP == E0
        val = bw6_761_g1_check_membership(P)
        ok = not val
        if ok:
            val = bw6_761_g1_check_membership(rP)
            ok = not val
            if ok:
                val = bw6_761_g1_check_membership(cP)
                ok = val
                if ok:
                    C = bw6_761_g1_mult_by_cofactor(P)
                    val = bw6_761_g1_check_membership(C)
                    ok = val
        i = i+1
    print("test bw6_761_g1_check_membership(P): {}".format(ok))
    return ok

def test_bw6_761_g1_check_membership_alt(E,r,c):
    #print("test bw6_761_g1_check_membership_alt(P)")
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        cP = c*P
        rP = r*P
        assert rP != E0 and cP != E0 and c*rP == E0
        val = bw6_761_g1_check_membership_alt(P)
        ok = not val
        if ok:
            val = bw6_761_g1_check_membership_alt(rP)
            ok = not val
            if ok:
                val = bw6_761_g1_check_membership_alt(cP)
                ok = val
                if ok:
                    C = bw6_761_g1_mult_by_cofactor_alt(P)
                    val = bw6_761_g1_check_membership_alt(C)
                    ok = val
        i = i+1
    print("test bw6_761_g1_check_membership_alt(P): {}".format(ok))
    return ok

# test subgroup membership testing and cofactor multiplication
def test_bw6_761_phi_g2(E2,r,c2):
    #print("test bw6_761_phi(Q) (Q in G2)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        assert r*c2*Q == E20
        phiQ = bw6_761_phi(Q)
        phi2Q = bw6_761_phi(phiQ)
        ok = phi2Q + phiQ + Q == E20
        i = i+1
    print("test bw6_761_phi(Q) (Q in G2): {}".format(ok))
    return ok

def test_bw6_761_g2_mult_by_cofactor(E2,r,c2):
    #print("test bw6_761_g2_mult_by_cofactor(Q)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        C = r*Q
        assert C != E20 and c2*Q != E20 and c2*C == E20
        R = bw6_761_g2_mult_by_cofactor(Q)
        ok = r*R == E20
        if ok:
            R = bw6_761_g2_mult_by_cofactor(C)
            ok = R == E20
        i = i+1
    print("test bw6_761_g2_mult_by_cofactor(Q): {}".format(ok))
    return ok

def test_bw6_761_g2_mult_by_cofactor_alt(E2,r,c2):
    #print("test bw6_761_g2_mult_by_cofactor_alt(Q)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        C = r*Q
        assert C != E20 and c2*Q != E20 and c2*C == E20
        R = bw6_761_g2_mult_by_cofactor_alt(Q)
        ok = r*R == E20
        if ok:
            R = bw6_761_g2_mult_by_cofactor(C)
            ok = R == E20
        i = i+1
    print("test bw6_761_g2_mult_by_cofactor_alt(Q): {}".format(ok))
    return ok

def test_bw6_761_g2_mult_by_r(E2,r,c2):
    #print("test bw6_761_g2_mult_by_r(Q)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        C = c2*Q
        assert C != E20 and r*Q != E20 and r*C == E20
        R = bw6_761_g2_mult_by_r(Q)
        ok = c2*R == E20
        if ok:
            R = bw6_761_g2_mult_by_r(C)
            ok = R == E20
        i = i+1
    print("test bw6_761_g2_mult_by_r(Q): {}".format(ok))
    return ok

def test_bw6_761_g2_mult_by_r_alt(E2,r,c2):
    #print("test bw6_761_g2_mult_by_r_alt(Q)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        C = c2*Q
        assert C != E20 and r*Q != E20 and r*C == E20
        R = bw6_761_g2_mult_by_r_alt(Q)
        ok = c2*R == E20
        if ok:
            R = bw6_761_g2_mult_by_r(C)
            ok = R == E20
        i = i+1
    print("test bw6_761_g2_mult_by_r_alt(Q): {}".format(ok))
    return ok

def test_bw6_761_g2_check_membership(E2,r,c2):
    #print("test bw6_761_g2_check_membership(Q)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        rQ = r*Q
        c2Q = c2*Q
        assert rQ != E20 and c2Q != E20 and c2*rQ == E20
        val = bw6_761_g2_check_membership(Q)
        ok = not val
        if ok:
            val = bw6_761_g2_check_membership(rQ)
            ok = not val
            if ok:
                val = bw6_761_g2_check_membership(c2Q)
                ok = val
                if ok:
                    C = bw6_761_g2_mult_by_cofactor(Q)
                    val = bw6_761_g2_check_membership(C)
                    ok = val
        i = i+1
    print("test bw6_761_g2_check_membership(Q): {}".format(ok))
    return ok

def test_bw6_761_g2_check_membership_alt(E2,r,c2):
    #print("test bw6_761_g2_check_membership_alt(Q)")
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        Q = E2.random_element()
        rQ = r*Q
        c2Q = c2*Q
        assert rQ != E20 and c2Q != E20 and c2*rQ == E20
        val = bw6_761_g2_check_membership_alt(Q)
        ok = not val
        if ok:
            val = bw6_761_g2_check_membership_alt(rQ)
            ok = not val
            if ok:
                val = bw6_761_g2_check_membership_alt(c2Q)
                ok = val
                if ok:
                    C = bw6_761_g2_mult_by_cofactor_alt(Q)
                    val = bw6_761_g2_check_membership_alt(C)
                    ok = val
        i = i+1
    print("test bw6_761_g2_check_membership(Q): {}".format(ok))
    return ok

def test_g2_frobenius_eigenvalue_bw6(E_Fqk,E2,r,c2,D_twist=False):
    """E_Fqk: elliptic curve defined over Fqk = Fq[w]/(w^d-xi), d = twist degree, Fq a prime finite field
    note that we cannot have Fq an extension: implementation is not yet available in Sage
    documentation from Fp.extension? says:
    Extensions of non-prime finite fields by polynomials are not yet supported:
    we fall back to generic code:
      sage: k.extension(x^5 + x^2 + x - 1)
      Univariate Quotient Polynomial Ring in x over Finite Field in z4 of size 3^4 with modulus x^5 + x^2 + x + 2
    """
    E20 = E2(0)
    ok = True
    Fqk = E_Fqk.base_field()
    w = Fqk.gen(0)
    q = Fqk.characteristic()
    i = 0
    while ok and i < 10:
        Q2 = c2 * E2.random_element()
        while Q2 == E20:
            Q2 = c2 * E2.random_element()
        assert r*Q2 == E20
        if D_twist:
            Q = E_Fqk(psi_sextic_d_twist(Q2, w))
        else:
            Q = E_Fqk(psi_sextic_m_twist(Q2, w))
        piQ = E_Fqk([(Q[0]).frobenius(), (Q[1]).frobenius()])
        ok = piQ == q*Q
        i = i+1
    print("test Frobenius(Q) == q*Q: {}".format(ok))
    return ok

def test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761_function, E, E_Fqk, E2, r, c, c2, u0, D_twist=False):
    """Test the bilinearity of the arg function

    INPUT:
    -`miller_loop_bw6_function`: a function
    -`E`: elliptic curve over prime field Fq of order c*r (subgroup of order r is G1)
    -`E_Fqk`: E(Fqk) (just a change in the coordinate domain)
    -`E2`: the d-twist over Fq^(k/d), actually Fq for BW6 curves
    -`r`: prime integer
    -`c`: curve cofactor
    -`c2`: twist cofactor
    -`u0`: seed
    -`D_twist`: is E2 a D-twist or M-twist of E

    Functions:
    miller_loop_opt_ate_bw6_761
    miller_loop_opt_ate_bw6_761_2naf
    miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u
    miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_2naf

    """
    # define two base-points P and Q for pairing tests
    E0 = E(0)
    E20 = E2(0)
    Fqk = E_Fqk.base_field()
    w = Fqk.gen(0)
    P = c * E.random_element()
    while P == E0 or r*P != E0:
        P = c * E.random_element()
    Q2 = c2*E2.random_element()
    while Q2 == E20 or r*Q2 != E20:
        Q2 = c2 * E2.random_element()
    if D_twist:
        Q = E_Fqk(psi_sextic_d_twist(Q2, w))
    else:
        Q = E_Fqk(psi_sextic_m_twist(Q2, w))
    f = miller_loop_opt_ate_bw6_761_function(Q,P,u0)
    assert (Fqk.cardinality()-1) % r == 0
    exponent = (Fqk.cardinality()-1) // r
    assert exponent % r != 0
    g = f**exponent
    assert g != 1 and g**r == 1

    ok = True
    aa = 1
    while ok and aa < 4:
        bb = 1
        while ok and bb < 4:
            fij = miller_loop_opt_ate_bw6_761_function(bb*Q,aa*P,u0)
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
    print("test bilinearity {}: {}".format(miller_loop_opt_ate_bw6_761_function.__name__, ok))
    return ok

def test_miller_loop_opt_ate_bw6_761_all(E,E_Fqk,E2,r,c,c2,u0,D_twist=False):
    #Fqk.cardinality() == q^k
    Fqk = E_Fqk.base_field()
    w = Fqk.gen(0)
    exponent = (Fqk.cardinality()-1) // r
    E0 = E(0)
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = c*E.random_element()
        while P == E0:
            P = c*E.random_element()
        Q2 = c2*E2.random_element()
        while Q2 == E20:
            Q2 = c2 * E2.random_element()
        if D_twist:
            Q = E_Fqk(psi_sextic_d_twist(Q2, w))
        else:
            Q = E_Fqk(psi_sextic_m_twist(Q2, w))
        
        m1 = miller_loop_opt_ate_bw6_761(Q,P,u0)
        m2 = miller_loop_opt_ate_bw6_761_2naf(Q,P,u0)
        f1 = m1**exponent
        f2 = m2**exponent
        ok = f1 == f2
        i = i+1
    print("test all miller_loop_opt_ate_bw6_761 (default,2naf): {}".format(ok))
    return ok

def test_final_exp_bw6_761(Fqk,r,u0,expected_exponent=None):
    #k = Fqk.degree()
    #q = Fqk.characteristic()
    #Fqk.cardinality() == q^k
    #exponent = (Fqk.cardinality()-1) // r
    # actually this function uses a multiple of this exponent
    Fqk_1 = Fqk(1)
    ok = True
    i = 0
    while ok and i < 10:
        f = Fqk.random_element()
        g = final_exp_bw6_761(f,u0)
        ok = g**r == Fqk_1
        if ok and expected_exponent is not None and expected_exponent != 0:
            ok = g == f**expected_exponent
        i += 1
    print("test final_exp_bw6_761: {}".format(ok))
    return ok

def test_optimal_ate_pairing_bw6_761(E,E_Fqk,E2,r,c,c2,u0,D_twist=False):
    E0 = E(0)
    E20 = E2(0)
    Fqk = E_Fqk.base_field()
    w = Fqk.gen(0)
    P = c * E.random_element()
    while P == E0:
        P = c * E.random_element()
    Q2 = c2*E2.random_element()
    while Q2 == E20:
        Q2 = c2 * E2.random_element()
    if D_twist:
        Q = E_Fqk(psi_sextic_d_twist(Q2, w))
    else:
        Q = E_Fqk(psi_sextic_m_twist(Q2, w))
    f = optimal_ate_pairing_bw6_761(Q,P,u0)

    ok = True
    aa = 1
    while ok and aa < 4:
        bb = 1
        while ok and bb < 4:
            fij = optimal_ate_pairing_bw6_761(bb*Q,aa*P,u0)
            ok = fij == f**(aa*bb)
            bb = bb+1
        aa = aa+1
    print("test bilinear optimal_ate_pairing_bw6_761: {}".format(ok))
    return ok

def test_formula_miller_loop_opt_ate_bw6_761_aklgl(E,E2,Fqk,r,c,c2,u0,D_twist=False):
    Q1 = c2*E2.random_element()
    assert r*Q1 == E2(0)
    P = c*E.random_element()
    assert r*P == E(0)
    Pxy = (P[0],P[1])
    Q1xy = (Q1[0],Q1[1])
    u = u0+1
    mu, Su = miller_function_ate_aklgl(Q1xy,Pxy,E2.a6(),u,Fqk,m0=1,D_twist=D_twist)
    Qu1 = (u0+1)*E2(Q1)
    ok1 = Su[0]/Su[2] == Qu1[0] and Su[1]/Su[2] == Qu1[1]
    #print("S == (u+1)*Q: {}".format(ok1))
    v = u0*(u0**2-u0-1)
    mv, Sv = miller_function_ate_aklgl(Q1xy,Pxy,E2.a6(),v,Fqk,m0=1,D_twist=D_twist)
    Qv = v*E2(Q1)
    ok2 = Sv[0]/Sv[2] == Qv[0] and Sv[1]/Sv[2] == Qv[1]
    #print("S == v*Q: {}".format(ok2))
    ok = ok1 and ok2
    print("test miller_function_ate_aklgl: {}".format(ok))
    return ok

def test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_761_aklgl_function,E,E2,Fqk,r,c,c2,u0,D_twist=False):
    """Test with functions:
    miller_loop_opt_ate_bw6_761_aklgl
    miller_loop_opt_ate_bw6_761_aklgl_2naf
    """
    P = c*E.random_element()
    while P == E(0):
        P = c*E.random_element()
    Q2 = c2*E2.random_element()
    while Q2 == E2(0):
        Q2 = c2 * E2.random_element()
    exponent = (Fqk.cardinality()-1) // r
    m = miller_loop_opt_ate_bw6_761_aklgl_function(Q2,P,E2.a6(),u0,Fqk,D_twist=D_twist)
    f = m**exponent
    ok = True
    aa = 1
    while ok and aa < 4:
        bb = 1
        while ok and bb < 4:
            mij = miller_loop_opt_ate_bw6_761_aklgl_function(bb*Q2,aa*P,E2.a6(),u0,Fqk,D_twist=D_twist)
            fij = mij**exponent
            ok = fij == f**(aa*bb)
            bb += 1
        aa += 1
    print("test bilinear {}: {}".format(miller_loop_opt_ate_bw6_761_aklgl_function.__name__, ok))
    return ok

def test_miller_loop_opt_ate_bw6_761_aklgl_all(E,E2,Fqk,r,c,c2,u0,D_twist=False):
    #Fqk.cardinality() == q^k
    exponent = (Fqk.cardinality()-1) // r
    E0 = E(0)
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = c*E.random_element()
        while P == E0:
            P = c*E.random_element()
        Q2 = c2*E2.random_element()
        while Q2 == E20:
            Q2 = c2 * E2.random_element()
        m1 = miller_loop_opt_ate_bw6_761_aklgl(Q2,P,E2.a6(),u0,Fqk,D_twist=D_twist)
        m2 = miller_loop_opt_ate_bw6_761_aklgl_2naf(Q2,P,E2.a6(),u0,Fqk,D_twist=D_twist)
        f1 = m1**exponent
        f2 = m2**exponent
        ok = f1 == f2
        i = i+1
    print("test all miller_loop_opt_ate_bw6_761_aklgl (default,2naf): {}".format(ok))
    return ok

def test_miller_loop_opt_ate_all_bw6_761_and_aklgl(E,E_Fqk,E2,r,c,c2,u0,D_twist=False):
    Fqk = E_Fqk.base_field()
    w = Fqk.gen(0)
    #Fqk.cardinality() == q^k
    exponent = (Fqk.cardinality()-1) // r
    E0 = E(0)
    E20 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = c*E.random_element()
        while P == E0:
            P = c*E.random_element()
        Q2 = c2*E2.random_element()
        while Q2 == E20:
            Q2 = c2 * E2.random_element()
        if D_twist:
            Q = E_Fqk(psi_sextic_d_twist(Q2, w))
        else:
            Q = E_Fqk(psi_sextic_m_twist(Q2, w))
        
        m1 = miller_loop_opt_ate_bw6_761(Q,P,u0)
        m2 = miller_loop_opt_ate_bw6_761_2naf(Q,P,u0)
        m4 = miller_loop_opt_ate_bw6_761_aklgl(Q2,P,E2.a6(),u0,Fqk,D_twist=D_twist)
        m5 = miller_loop_opt_ate_bw6_761_aklgl_2naf(Q2,P,E2.a6(),u0,Fqk,D_twist=D_twist)
        f1 = m1**exponent
        f2 = m2**exponent
        f4 = m4**exponent
        f5 = m5**exponent
        ok = f1 == f2 and f2 == f4 and f4 == f5
        i = i+1
    print("test all miller_loop_opt_ate_bw6_761_aklgl and bw6_761 (default,2naf): {}".format(ok))
    return ok

def cost_pairing_bw6_761():
    # cost
    # a=0
    # AKLGL'11
    u0 = ZZ(0x8508C00000000001)
    m = 1
    s = 1
    m3 = 6*m
    s3 = 5*m
    m6 = 18*m
    s6 = 12*m
    s2 = 2*m
    s6_cyclo = 3*s2
    f6 = 4*m
    f3 = 0 # frobenius to the power p^3 costs 3 negations
    i = 25*m
    i6 = 34*m+i
    k = 6
    DoubleLine_ate = 3*m+6*s+(k//3)*m
    AddLine_ate    = 11*m+2*s+(k//3)*m
    Update1        = 13*m+s6
    Update2        = 13*m
    
    HW2naf_u = sum([1 for bi in bits_2naf(u0) if bi != 0])
    cost_ate1 = (u0.nbits()-1)*DoubleLine_ate + (u0.nbits()-2)*Update1 + (HW2naf_u-1)*(AddLine_ate+Update2)
    v = u0**2-u0-1
    HW2naf_v = sum([1 for bi in bits_2naf(v) if bi != 0])
    cost_ate2 = (v.nbits()-1)*(DoubleLine_ate + Update1) + (HW2naf_v-1)*(AddLine_ate+Update2+m6)
    cost_ate = cost_ate1 + i+2*m+i6+AddLine_ate+Update2 + cost_ate2 + f6 + m6
    
    print("cost ate Miller = {}m".format(cost_ate))
    print("({0}-1)*(3*m+6*s+(k//3)*m) + ({0}-2)*(13*m+s6) + ({1}-1)*(11*m+2*s+(k//3)*m+13*m)".format(u0.nbits(),HW2naf_u))
    print("+ ({}-1)*(3*m+6*s+(k//3)*m+13*m+s6) + ({}-1)*(11*m+2*s+(k//3)*m+13*m+m6)".format(v.nbits(),HW2naf_v))
    print("+ i+2*m+i6+13*m+f6+m6")
    cost_easy_part = f3 + f6 + 2*m6 + i6 # f3 + i6 + 2*m6 + f
    # final exp hard:51 m6 + 9 s6_cyclo + 9 exponentiations to u0
    bits_2naf_u0 = bits_2naf(abs(u0))
    Hw2naf_u0 = sum([1 for bi in bits_2naf_u0 if bi != 0])
    cost_exp_u0 = (len(bits_2naf_u0)-1) * s6_cyclo + (Hw2naf_u0-1) * m6
    cost_hard_part = 9*cost_exp_u0 + 8*f6 + 10*f3 + 51*m6 + 9*s6_cyclo
    print("cost final exp easy: {}m".format(cost_easy_part))
    print("cost final exp hard: {}m (with cyclotomic squaring)".format(cost_hard_part))
    print("cost final exp:      {}m".format(cost_easy_part + cost_hard_part))

# even faster exponentiation
def script_fast_final_exponentiation():
    ZZx = ZZ['x']; (x,) = ZZx._first_ngens(1)
    R0 = -103*x**7+70*x**6+269*x**5-197*x**4-314*x**3-73*x**2-263*x-220
    R1 = 103*x**9-276*x**8+77*x**7+492*x**6-445*x**5-65*x**4+452*x**3-181*x**2+34*x+229
    
    for ri in (R0.list() + R1.list()):
        if ri < 0:
            print("{}".format([-ei for ei in bits_2naf(-ri)]))
        else:
            print(bits_2naf(ri))
    
    j=0; j0=True
    HW = 0
    for ri in (R0.list() + R1.list()):
        ri = int(ri)
        if ri < 0:
            C = [-ei for ei in bits_2naf(-ri)]
            HW += sum([1 for ei in C if ei != 0])
            C += [0 for i in range(len(C)-1,9)]
            C.reverse()
            s = ""
            for ci in C:
                s += "{: 2d}".format(ci)
            print("{} {} {: 3d}".format(j,s,ri))
        else:
            C = bits_2naf(ri)
            HW += sum([1 for ei in C if ei != 0])
            C += [0 for i in range(len(C)-1,9)]
            C.reverse()
            s = ""
            for ci in C:
                s += "{: 2d}".format(ci)
            print("{} {} {: 3d}".format(j,s,ri))
        j += 1
        if j == len(R0.list()) and j0:
            j = 0
            j0 = False
    print("HW = {}".format(HW))

    j=0; j0=True
    HW = 0
    for ri in (R0.list() + R1.list()):
        ri = ZZ(ri)
        if ri < 0:
            C = [-ei for ei in (-ri).digits(2)]
            HW += sum([1 for ei in C if ei != 0])
            C += [0 for i in range(len(C)-1,9)]
            C.reverse()
            s = ""
            for ci in C:
                s += "{: 2d}".format(ci)
            print("{} {} {: 3d}".format(j,s,ri))
        else:
            C = (ri).digits(2)
            HW += sum([1 for ei in C if ei != 0])
            C += [0 for i in range(len(C)-1,9)]
            C.reverse()
            s = ""
            for ci in C:
                s += "{: 2d}".format(ci)
            print("{} {} {: 3d}".format(j,s,ri))
        j += 1
        if j == len(R0.list()) and j0:
            j = 0
            j0 = False
    print("HW = {}".format(HW))


"""
0  0-1 0 0 1 0 0 1 0 0 -220
1  0-1 0 0 0 0-1 0 0 1 -263
2  0 0 0-1 0 0-1 0 0-1 -73
3  0-1 0-1 0 0 1 0-1 0 -314
4  0-1 0 1 0 0 0-1 0-1 -197
5  0 1 0 0 0 1 0-1 0 1  269
6  0 0 0 1 0 0 1 0-1 0  70
7  0 0-1 0 1 0-1 0 0 1 -103
0  0 1 0 0-1 0 0 1 0 1  229
1  0 0 0 0 1 0 0 0 1 0  34
2  0-1 0 1 0 1 0-1 0-1 -181
3  1 0 0-1 0 0 0 1 0 0  452
4  0 0 0-1 0 0 0 0 0-1 -65
5 -1 0 0 1 0 0 0 1 0-1 -445
6  1 0 0 0 0-1 0-1 0 0  492
7  0 0 0 1 0 1 0-1 0 1  77
8  0-1 0 0 0-1 0-1 0 0 -276
9  0 0 1 0-1 0 1 0 0-1  103
HW = 62
0  0 0-1-1 0-1-1-1 0 0 -220
1  0-1 0 0 0 0 0-1-1-1 -263
2  0 0 0-1 0 0-1 0 0-1 -73
3  0-1 0 0-1-1-1 0-1 0 -314
4  0 0-1-1 0 0 0-1 0-1 -197
5  0 1 0 0 0 0 1 1 0 1  269
6  0 0 0 1 0 0 0 1 1 0  70
7  0 0 0-1-1 0 0-1-1-1 -103
0  0 0 1 1 1 0 0 1 0 1  229
1  0 0 0 0 1 0 0 0 1 0  34
2  0 0-1 0-1-1 0-1 0-1 -181
3  0 1 1 1 0 0 0 1 0 0  452
4  0 0 0-1 0 0 0 0 0-1 -65
5  0-1-1 0-1-1-1-1 0-1 -445
6  0 1 1 1 1 0 1 1 0 0  492
7  0 0 0 1 0 0 1 1 0 1  77
8  0-1 0 0 0-1 0-1 0 0 -276
9  0 0 0 1 1 0 0 1 1 1  103
HW = 76
"""

if __name__ == "__main__":
    arithmetic(False)
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    u0 = ZZ(0x8508C00000000001)
    #proof.arithmetic(False)
    #BW6_761 parameters
    q = ZZ(6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786613451161759120554247759349511699125301598951605099378508850372543631423596795951899700429969112842764913119068299)
    Fq = GF(q, proof=False)
    #Fqz.<z> = Fq[]
    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1)
    r = ZZ(258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177)
    c = ZZ(26642435879335816683987677701488073867751118270052650655942102502312977592501693353047140953112195348280268661194876)
    t = ZZ(3362637538168598222219435186298528655381674028954528064283340709388076588006567983337308081752755143497537638367248)
    y = ZZ(2327979834116721846122857819342346041630394402507777770613906795574054381627779834062290838568927395079900712927242)
    k = 6
    rx = (x**6 - 2*x**5 + 2*x**3 + x + 1)/3
    qx = (103*x**12 - 379*x**11 + 250*x**10 + 691*x**9 - 911*x**8 - 79*x**7 + 623*x**6 - 640*x**5 + 274*x**4 + 763*x**3 + 73*x**2 + 254*x + 229)/9
    tx = (13*x**6 - 23*x**5 - 9*x**4 + 35*x**3 + 10*x + 22)/3
    yx = (9*x**6 - 17*x**5 - 3*x**4 + 21*x**3 + 8*x + 12)/3
    cx = (103*x**6 - 173*x**5 - 96*x**4 + 293*x**3 + 21*x**2 + 52*x + 172)/3
    c2x = (103*x**6 - 173*x**5 - 96*x**4 + 293*x**3 + 21*x**2 + 52*x + 151)/3
    
    q_ = ZZ(qx(u0))
    r_ = ZZ(rx(u0))
    c_ = ZZ(cx(u0))
    t_ = ZZ(tx(u0))
    y_ = ZZ(yx(u0))
    c2_ = ZZ(c2x(u0))
    c2 = c2_
    
    assert t**2 + 3*y**2 == 4*q
    assert t_**2 + 3*y_**2 == 4*q_
    assert q == q_
    assert r == r_
    assert c == c_
    assert t == t_
    assert y == y_
    
    E = EllipticCurve([Fq(0),Fq(-1)])
    b = Fq(-1)
    print(E)
    P = E.random_element()
    assert r*c*P == E(0)
    
    # M-twist
    E_M = EllipticCurve([Fq(0),Fq(4)])
    print(E_M)
    QM = E_M.random_element()
    assert c2*r *QM == E_M(0)
    #Fq6M.<wM> = Fq.extension(z**6 + 4)
    #preparse("Fq6M.<wM> = Fq.extension(z**6 + 4)")
    Fq6M = Fq.extension(z**6 + 4, names=('wM',),proof=False); (wM,) = Fq6M._first_ngens(1)
    E6M = EllipticCurve([Fq6M(0),Fq6M(-1)])

    E_D = EllipticCurve([Fq(0), Fq(-1)/Fq(2)])
    print(E_D)
    QD = E_D.random_element()
    assert c2*r *QD == E_D(0)
    #Fq6D.<wD> = Fq.extension(z**6 - 2)
    Fq6D = Fq.extension(z**6 - 2, names=('wD',),proof=False); (wD,) = Fq6D._first_ngens(1)
    E6D = EllipticCurve([Fq6D(0),Fq6D(-1)])

    xiD = Fq(2)
    xiM = Fq(-4)
    try:
        test_xiD = -Fq6D.modulus().constant_coefficient()
        print("xiD == -Fq6D.modulus().constant_coefficient(): {}".format(xiD == test_xiD))
    except AttributeError as err:
        print("xiD = -Fq6D.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiD = -Fq6D.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiD == -Fq6D.polynomial().constant_coefficient(): {}".format(xiD == test_xiD))
    except AttributeError as err:
        print("xiD = -Fq6D.polynomial().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiM = -Fq6M.modulus().constant_coefficient()
        print("xiM == -Fq6M.modulus().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fq6M.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiM = -Fq6M.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiM == -Fq6M.polynomial().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fq6M.polynomial().constant_coefficient() raised an error:\n{}".format(err))

    def map_ED_E_Fq6D(Q2):
        return E6D(psi_sextic_d_twist(Q2, wD))
    
    def map_EM_E_Fq6M(Q2):
        return E6M(psi_sextic_m_twist(Q2, wM))
    
    # final exponentiation
    R0 = -103*x**7+70*x**6+269*x**5-197*x**4-314*x**3-73*x**2-263*x-220;
    R1 = 103*x**9-276*x**8+77*x**7+492*x**6-445*x**5-65*x**4+452*x**3-181*x**2+34*x+229;
    exponent = cyclotomic_polynomial(6)(qx) // rx
    assert (R0 + R1*qx) == 3*(x**3 - x**2 + 1)*exponent
    e0 = ZZ((((x**k-1) // cyclotomic_polynomial(k))(qx))(u0))
    e1 = ZZ(exponent(u0))
    ee = ZZ((R0 + R1*qx)(u0))

    print("tests G1")
    test_bw6_761_phi(E,r,c)
    test_bw6_761_g1_mult_by_cofactor(E,r,c)
    test_bw6_761_g1_mult_by_cofactor_alt(E,r,c)
    test_bw6_761_g1_mult_by_r(E,r,c)
    test_bw6_761_g1_mult_by_r_alt(E,r,c)
    test_bw6_761_g1_check_membership(E,r,c)
    test_bw6_761_g1_check_membership_alt(E,r,c)
    
    print("\ntests G2 with M-twist")
    test_bw6_761_phi_g2(E_M,r,c2)
    test_bw6_761_g2_mult_by_cofactor(E_M,r,c2)
    test_bw6_761_g2_mult_by_cofactor_alt(E_M,r,c2)
    test_bw6_761_g2_mult_by_r(E_M,r,c2)
    test_bw6_761_g2_mult_by_r_alt(E_M,r,c2)
    test_bw6_761_g2_check_membership(E_M,r,c2)
    test_bw6_761_g2_check_membership_alt(E_M,r,c2)
    
    print("\ntests G2 with D-twist")
    test_bw6_761_phi_g2(E_D,r,c2)
    test_bw6_761_g2_mult_by_cofactor(E_D,r,c2)
    test_bw6_761_g2_mult_by_cofactor_alt(E_D,r,c2)
    test_bw6_761_g2_mult_by_r(E_D,r,c2)
    test_bw6_761_g2_mult_by_r_alt(E_D,r,c2)
    test_bw6_761_g2_check_membership(E_D,r,c2)
    test_bw6_761_g2_check_membership_alt(E_D,r,c2)

    print("\ntest pairings with M-twist")
    test_g2_frobenius_eigenvalue_bw6(E6M,E_M,r,c2,D_twist=False)
    test_double_line_j(E,E_M,Fq6M,D_twist=False)
    test_add_line_j(E,E_M,Fq6M,D_twist=False)
    test_double_line_j_csb(E,E_M,Fq6M,D_twist=False)
    test_add_line_j_csb(E,E_M,Fq6M,D_twist=False)

    test_miller_function_tate(E, E6M, E_M, r, c, c2, D_twist=False)
    test_miller_function_tate_2naf(E, E6M, E_M, r, c, c2, D_twist=False)
    test_miller_function_ate(E, E6M, E_M, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_2naf(E, E6M, E_M, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_csb(E, E6M, E_M, r, c, c2, t-1, D_twist=False)

    test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761,E,E6M,E_M,r,c,c2,u0,D_twist=False)
    test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761_2naf,E,E6M,E_M,r,c,c2,u0,D_twist=False)
    test_miller_loop_opt_ate_bw6_761_all(E,E6M,E_M,r,c,c2,u0,D_twist=False)
    
    test_final_exp_bw6_761(Fq6M,r,u0)
    test_final_exp_bw6_761(Fq6M,r,u0,expected_exponent=e0*ee)
    test_optimal_ate_pairing_bw6_761(E,E6M,E_M,r,c,c2,u0,D_twist=False)

    test_double_line_h_a0_twist6_aklgl(E,E_M,Fq6M,r,c,c2,D_twist=False,verbose=False)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,E_M,Fq6M,r,c,c2,D_twist=False,verbose=False)
    test_add_line_h_a0_twist6_aklgl(E,E_M,xiM,D_twist=False)
    test_add_line_h_a0_twist6_aklgl_test(E,E_M,xiM,D_twist=False)
    test_add_line_h_a0_twist6_aklgl_with_z(E,E_M,xiM,D_twist=False)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_m6_twist(Fq6M)
    
    test_formula_miller_loop_opt_ate_bw6_761_aklgl(E,E_M,Fq6M,r,c,c2,u0,D_twist=False)

    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_761_aklgl,E,E_M,Fq6M,r,c,c2,u0,D_twist=False)
    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_761_aklgl_2naf,E,E_M,Fq6M,r,c,c2,u0,D_twist=False)
    test_miller_loop_opt_ate_bw6_761_aklgl_all(E,E_M,Fq6M,r,c,c2,u0,D_twist=False)
    test_miller_loop_opt_ate_all_bw6_761_and_aklgl(E,E6M,E_M,r,c,c2,u0,D_twist=False)

    print("\ntest pairings with D-twist")
    
    test_g2_frobenius_eigenvalue_bw6(E6D,E_D,r,c2,D_twist=True)
    test_double_line_j(E,E_D,Fq6D,D_twist=True)
    test_add_line_j(E,E_D,Fq6D,D_twist=True)
    test_double_line_j_csb(E,E_D,Fq6D,D_twist=True)
    test_add_line_j_csb(E,E_D,Fq6D,D_twist=True)

    test_miller_function_tate(E, E6D, E_D, r, c, c2, D_twist=True)
    test_miller_function_tate_2naf(E, E6D, E_D, r, c, c2, D_twist=True)
    test_miller_function_ate(E, E6D, E_D, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_2naf(E, E6D, E_D, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_csb(E, E6D, E_D, r, c, c2, t-1, D_twist=True)

    test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761,E,E6D,E_D,r,c,c2,u0,D_twist=True)
    test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761_2naf,E,E6D,E_D,r,c,c2,u0,D_twist=True)
    test_miller_loop_opt_ate_bw6_761_all(E,E6D,E_D,r,c,c2,u0,D_twist=True)
    
    test_final_exp_bw6_761(Fq6D,r,u0)
    test_final_exp_bw6_761(Fq6D,r,u0,expected_exponent=e0*ee)
    test_optimal_ate_pairing_bw6_761(E,E6D,E_D,r,c,c2,u0,D_twist=True)
    
    test_double_line_h_a0_twist6_aklgl(E,E_D,Fq6D,r,c,c2,D_twist=True,verbose=False)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,E_D,Fq6D,r,c,c2,D_twist=True,verbose=False)
    test_add_line_h_a0_twist6_aklgl(E,E_D,xiD,D_twist=True)
    test_add_line_h_a0_twist6_aklgl_test(E,E_D,xiD,D_twist=True)
    test_add_line_h_a0_twist6_aklgl_with_z(E,E_D,xiD,D_twist=True)
    test_sparse_sparse_mult_d6_twist(Fq6D)
    test_sparse_mult_d6_twist(Fq6D)
    
    test_formula_miller_loop_opt_ate_bw6_761_aklgl(E,E_D,Fq6D,r,c,c2,u0,D_twist=True)
    
    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_761_aklgl,E,E_D,Fq6D,r,c,c2,u0,D_twist=True)
    test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_761_aklgl_2naf,E,E_D,Fq6D,r,c,c2,u0,D_twist=True)
    test_miller_loop_opt_ate_bw6_761_aklgl_all(E,E_D,Fq6D,r,c,c2,u0,D_twist=True)
    test_miller_loop_opt_ate_all_bw6_761_and_aklgl(E,E6D,E_D,r,c,c2,u0,D_twist=True)
    
    script_fast_final_exponentiation()
    cost_pairing_bw6_761()
