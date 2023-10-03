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

from pairing import *
from pairing_cocks_pinch import *
from test_pairing import *

#from test_scalar_mult import test_glv_scalar_mult_g1

from testvector_bls12_377_cp8d1_768 import testvector_bls12_377_cp8d1_768

from cost_pairing import cost_pairing_cp8_cp12

def test_miller_function_ate_cln_b0_4_multi_scalar(E, E2, E_Fqd, E_Fpk, p, r, tr, c, c2, miller_scalars, map_Fqd_Fpk, D_twist=True, function_name=miller_function_ate_cln_b0_4_multi_scalar):
    # test miller_function_ate_cln_b0_4_multi_scalar

    P = c*E.random_element()
    while P == E(0):
        P = c * E.random_element()
    Q = c2*E2.random_element()
    while Q == E2(0):
        Q = c2 * E2.random_element()
    Fqd = E_Fqd.base_field() # assume that the degree d extension is explicit
    Fpk = E_Fpk.base_field()
    w = Fqd.gen(0)
    exponent = (Fqd.cardinality()-1) // r # (q^d-1)//r = (p^k-1)//r
    if not D_twist:
        Q2 = psi_sextic_m_twist(Q, w)
    else:
        Q2 = psi_sextic_d_twist(Q, w)
    Q2 = E_Fpk([map_Fqd_Fpk(Q2[0]), map_Fqd_Fpk(Q2[1])])

    pQ2 = E_Fpk([Q2[0].frobenius(), Q2[1].frobenius()])
    p2Q2 = E_Fpk([Q2[0].frobenius(2), Q2[1].frobenius(2)])
    p3Q2 = E_Fpk([Q2[0].frobenius(3), Q2[1].frobenius(3)])

    a0, a1, a2, a3 = miller_scalars
    
    # note that E.a4() == E_Fpk.a4()
    f, scalarsQ = function_name(Q2, P, E.a4(), miller_scalars)
    ok = f != Fpk(1)
    if not ok:
        print("test bilinearity {}: {} pairing result is 1".format(function_name.__name__, ok))
        return False
    #print("type(scalarsQ) = {}".format(type(scalarsQ)))
    check_ai_piQ = ((type(scalarsQ) is list) or (type(scalarsQ) is tuple)) and len(scalarsQ) == 4
    if check_ai_piQ:
        a0Q, a1Q, a2Q, a3Q = scalarsQ
        a0Q = (a0Q[0]/a0Q[2], a0Q[1]/a0Q[2]**2)
        a1Q = (a1Q[0]/a1Q[2], a1Q[1]/a1Q[2]**2)
        a2Q = (a2Q[0]/a2Q[2], a2Q[1]/a2Q[2]**2)
        a3Q = (a3Q[0]/a3Q[2], a3Q[1]/a3Q[2]**2)
        ok_a0Q = E_Fpk(a0Q) == a0*Q2
        ok_a1Q = E_Fpk(a1Q) == a1*pQ2
        ok_a2Q = E_Fpk(a2Q) == a2*p2Q2
        ok_a3Q = E_Fpk(a3Q) == a3*p3Q2
    else:
        # should be zero?
        ok_ai_piQ = scalarsQ[2] == 0
    g = f**exponent
    aa = 1
    bb = 1
    while ok and aa < 4:
        Pa = aa*P
        bb = 1
        while ok and bb < 4:
            Qb = bb*Q
            if not D_twist:
                Q2b = psi_sextic_m_twist(Qb, w)
            else:
                Q2b = psi_sextic_d_twist(Qb, w)
            Q2b = E_Fpk([map_Fqd_Fpk(Q2b[0]), map_Fqd_Fpk(Q2b[1])])
            pQ2b = E_Fpk([Q2b[0].frobenius(), Q2b[1].frobenius()])
            p2Q2b = E_Fpk([Q2b[0].frobenius(2), Q2b[1].frobenius(2)])
            p3Q2b = E_Fpk([Q2b[0].frobenius(3), Q2b[1].frobenius(3)])
            fab, scalars_abQ = function_name(Q2b, Pa, E.a4(), miller_scalars)
            gab = fab**exponent
            ok = ok and gab != Fpk(1) and gab == g**(aa*bb)
            if check_ai_piQ:
                a0Q, a1Q, a2Q, a3Q = scalars_abQ
                a0Q = (a0Q[0]/a0Q[2], a0Q[1]/a0Q[2]**2)
                a1Q = (a1Q[0]/a1Q[2], a1Q[1]/a1Q[2]**2)
                a2Q = (a2Q[0]/a2Q[2], a2Q[1]/a2Q[2]**2)
                a3Q = (a3Q[0]/a3Q[2], a3Q[1]/a3Q[2]**2)
                ok_a0Q = ok_a0Q and E_Fpk(a0Q) == a0*Q2b
                ok_a1Q = ok_a1Q and E_Fpk(a1Q) == a1*pQ2b
                ok_a2Q = ok_a2Q and E_Fpk(a2Q) == a2*p2Q2b
                ok_a3Q = ok_a3Q and E_Fpk(a3Q) == a3*p3Q2b
            else:
                ok_ai_piQ = ok_ai_piQ and scalars_abQ[2] == 0
            bb += 1
        aa += 1
    print("test bilinearity {}: {} ({} tests)".format(function_name.__name__, ok, (aa-1)*(bb-1)))
    if check_ai_piQ:
        print("aiQ == [ai]*pi_p^i(Q): {}".format(ok_a0Q and ok_a1Q and ok_a2Q and ok_a3Q))
    else:
        print("sum ai*piQ == 0: {}".format(ok_ai_piQ))
    return ok

def test_miller_function_tate_cln_b0_2_multi_scalar(E, E2, E_Fqd, E_Fpk, p, r, tr, c, c2, miller_scalars, map_Fqd_Fpk, D_twist=True, function_name=miller_function_tate_cln_b0_2_multi_scalar):
    # test miller_function_tate_cln_b0_2_multi_scalar

    P = c*E.random_element()
    while P == E(0):
        P = c * E.random_element()
    Q = c2*E2.random_element()
    while Q == E2(0):
        Q = c2 * E2.random_element()
    Fqd = E_Fqd.base_field() # assume that the degree d extension is explicit
    Fpk = E_Fpk.base_field()
    w = Fqd.gen(0)
    exponent = (Fqd.cardinality()-1) // r # (q^d-1)//r = (p^k-1)//r
    if not D_twist:
        Q2 = psi_sextic_m_twist(Q, w)
    else:
        Q2 = psi_sextic_d_twist(Q, w)
    Q2 = E_Fpk([map_Fqd_Fpk(Q2[0]), map_Fqd_Fpk(Q2[1])])

    a0, a1, omega = miller_scalars
    psiP = E([-P[0], P[1] * omega])

    # note that E.a4() == E_Fpk.a4()
    f, scalarsP = function_name(P, Q2, E.a4(), a0, a1, omega)
    ok = f != Fpk(1)
    if not ok:
        print("test bilinearity {}: {} pairing result is 1".format(function_name.__name__, ok))
        return False
    #print("type(scalarsP) = {}".format(type(scalarsP)))
    check_ai_P = ((type(scalarsP) is list) or (type(scalarsP) is tuple)) and len(scalarsP) == 2
    if check_ai_P:
        a0P, a1P = scalarsP
        a0P = (a0P[0]/a0P[2], a0P[1]/a0P[2]**2)
        a1P = (a1P[0]/a1P[2], a1P[1]/a1P[2]**2)
        ok_a0P = E(a0P) == a0*P
        ok_a1P = E(a1P) == a1*psiP
    else:
        # should be zero?
        ok_ai_P = scalarsP[2] == 0
    g = f**exponent
    aa = 1
    bb = 1
    while ok and aa < 4:
        Pa = aa*P
        psiPa = E([-Pa[0], Pa[1] * omega])
        bb = 1
        while ok and bb < 4:
            Qb = bb*Q
            if not D_twist:
                Q2b = psi_sextic_m_twist(Qb, w)
            else:
                Q2b = psi_sextic_d_twist(Qb, w)
            Q2b = E_Fpk([map_Fqd_Fpk(Q2b[0]), map_Fqd_Fpk(Q2b[1])])

            fab, scalars_abP = function_name(Pa, Q2b, E.a4(), a0, a1, omega)
            gab = fab**exponent
            ok = ok and gab != Fpk(1) and gab == g**(aa*bb)
            if check_ai_P:
                a0P, a1P = scalars_abP
                a0P = (a0P[0]/a0P[2], a0P[1]/a0P[2]**2)
                a1P = (a1P[0]/a1P[2], a1P[1]/a1P[2]**2)
                ok_a0P = ok_a0P and E(a0P) == a0*Pa
                ok_a1P = ok_a1P and E(a1P) == a1*psiPa
            else:
                ok_ai_P = ok_ai_P and scalars_abP[2] == 0
            bb += 1
        aa += 1
    print("test bilinearity {}: {} ({} tests)".format(function_name.__name__, ok, (aa-1)*(bb-1)))
    if check_ai_P:
        print("a0P == [a0]P, a1P == [a1]psi(P): {}".format(ok_a0P and ok_a1P))
    else:
        print("a0*P + a1*psiP == 0: {}".format(ok_ai_P))
    return ok

def test_curve(test_vector_value):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    v = test_vector_value
    k = v['k']
    u0 = ZZ(v['u'])
    u = u0
    D = v['D']
    ht = v['ht']
    hy = v['hy']
    sign_ht = v['sign_ht']
    sign_y0 = v['sign_y0']
    tr0 = ZZ(v['tr0'])
    y0 = ZZ(v['y0'])
    tr = ZZ(v['tr'])
    t = tr
    y = ZZ(v['y'])
    p = ZZ(v['p'])
    r = ZZ(v['r'])
    miller_scalars = [ZZ(ai) for ai in v['miller']]
    a = v['a']
    pnbits = v['pnbits']
    rnbits = v['rnbits']
    rx = v['rx']
    rx_denom = v['rx_denom']
    label = v['label']
    
    rx = QQx(rx)/QQ(rx_denom)
    assert ((p+1-t) % r) == 0
    c = (p+1-t)//r
    # compute c2 (cofactor of twist E2) and E2 (curve coefficient)
    print("{}".format(label))
    print("ht= {} sign_ht = {}".format(ht, sign_ht))
    print("hy= {} sign_y0 = {}".format(hy, sign_y0))
    print("u = {:#x}".format(u))
    print("p = {:#x} # {} bits".format(p, p.nbits()))
    print("r = {:#x}".format(r))
    print("c = {:#x}".format(c))
    print("y = {:#x}".format(y))
    print("t = {:#x}".format(t))

    # find smallest integer a so that given x mod p, x^a is a permutation
    # will work for gcd(a, p-1) == 1
    # since 2 | (p-1) by construction, start at 3.
    a_perm = 3
    while gcd(a_perm, p-1) > 1:
        a_perm = a_perm + 1
    print("Smallest a such that x^a mod p is a permutation: a={}".format(a_perm))
    
    print("Elligator2 available: {} (<=> the curve has a 2-torsion point <=> cofactor is even)".format((c % 2) == 0))
    print("small factors of y that are also factors of (p+1-t) or (p+1+t):")
    print("gcd(y, (p+1-t)*(p+1+t)) = {}".format( (gcd(y, (p+1-t)*(p+1+t))).factor()))
    # try with quartic twists
    print("gcd(y, (p+1-y)*(p+1+y)) = {}".format( (gcd(y, (p+1-y)*(p+1+y))).factor()))
    print("--> isogenies of these degrees are defined over Fp.")
    
    curve_order = p+1-t
    assert curve_order == r*c
    Fp = GF(p, proof=False)
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    E = EllipticCurve([Fp(a), Fp(0)])
    a = Fp(a)
    print(E)
    P = E.random_element()
    assert curve_order*P == E(0)
    # define a quadratic extension
    a0 = 1
    while not (z**2 + a0).is_irreducible():
        if a0 < 0:
            a0 = -a0
        else:
            a0 = -(a0 + 1)
    Fq = Fp.extension(z**2 + a0, names=('v',)); (v,) = Fq._first_ngens(1)
    FqZ = Fq['Z']; (Z,) = FqZ._first_ngens(1)
    # compute g2c
    p2 = p**2
    tr2 = tr**2 - 2*p
    y2 = tr * y
    if ((p2+1-y2) % r) == 0:
        twist_order = (p2+1-y2)
        g2c = twist_order // r
    elif ((p2+1+y2) % r) == 0:
        twist_order = (p2+1+y2)
        g2c = twist_order // r
    else:
        raise ValueError("Error no appropriate quartic twist")
    
    xiD, aD = find_twist_curve_parameter_xi_ab(a, Fq, r, g2c, d=4, D_twist=True)
    
    c2 = g2c

    print("xiD = {}".format(xiD))
    ED = EllipticCurve([aD, Fq(0)])
    Q = ED.random_element()
    orderD = g2c * r
    assert orderD * Q == ED(0)
    print("ED ok")
    Fq4D = Fq.extension(Z**4 - xiD, names=('wD',)); (wD,) = Fq4D._first_ngens(1)
    E_Fq4D = EllipticCurve([Fq4D(a), Fq4D(0)])
    # absolute extension
    # Z**4 == xiD
    # xiM has the form v0 + v1*v in Fq
    v0D, v1D = xiD.polynomial().list()
    assert xiD == v0D + v1D*v
    # (Z**4 - v0) == v1*v
    # (Z**4 - v0)**2 == v1**2 * v**2 = v1**2 * (-a0)
    # (Z**4 - v0)**2 + a0 * v1**2
    Fp8D = Fp.extension((z**4 - v0D)**2 + a0*v1D**2, names=('vD',)); (vD,) = Fp8D._first_ngens(1)
    E_Fp8D = EllipticCurve([Fp8D(a), Fp8D(0)])

    def map_Fq4D_Fp8D(X, aD=None):
        if aD is None:
            aD = vD
        return sum([xi.polynomial()((aD**4-v0D)/v1D) * aD**e for e,xi in enumerate(X.list())])

    def map_ED_E_Fq4D(Q2):
        return E_Fq4D(psi_sextic_d_twist(Q2, wD))

    # now find a 4-th M-twist
    xiM, aM = find_twist_curve_parameter_xi_ab(a, Fq, r, g2c, d=4, D_twist=False)
    print("xiM = {}".format(xiM))
    EM = EllipticCurve([aM, Fq(0)])
    Q = EM.random_element()
    orderM = g2c * r
    assert orderM * Q == EM(0)
    print("EM ok")
    Fq4M = Fq.extension(Z**4 - xiM, names=('wM',)); (wM,) = Fq4M._first_ngens(1)
    E_Fq4M = EllipticCurve([Fq4M(a), Fq4M(0)])
    # absolute extension
    # Z**4 == xiM
    # xiM has the form v0 + v1*v in Fq
    v0M, v1M = xiM.polynomial().list()
    assert xiM == v0M + v1M*v
    # (Z**4 - v0) == v1*v
    # (Z**4 - v0)**2 == v1**2 * v**2 = v1**2 * (-a0)
    # (Z**4 - v0)**2 + a0 * v1**2
    Fp8M = Fp.extension((z**4 - v0M)**2 + a0*v1M**2, names=('vM',)); (vM,) = Fp8M._first_ngens(1)
    E_Fp8M = EllipticCurve([Fp8M(a), Fp8M(0)])
    
    def map_Fq4M_Fp8M(X, aM=None):
        if aM is None:
            aM = vM
        return sum([xi.polynomial()((aM**4-v0M)/v1M) * aM**e for e,xi in enumerate(X.list())])

    def map_EM_E_Fq4M(Q2):
        return E_Fq4M(psi_sextic_m_twist(Q2, wM))

    a0, a1, a2, a3 = miller_scalars
    #final_exp_hard
    ee = (p**4+1)//r
    M = Matrix(ZZ, 4, 4, [p,0,0,0, ee,1,0,0, 0,ee,1,0, 0,0,ee,1])
    R = M.LLL()
    e0, e1, e2, e3 = (R[0]).list()
    # what would be a Tate pairing? There is an endomorphism on G1, of eigenvalue sqrt(-1) mod r.
    Fr = GF(r, proof=False)
    sqrt_1_mod_r = Fr(tr0-2) / Fr(y0)
    assert sqrt_1_mod_r**2 + 1 == 0
    sqrt_1_mod_r = ZZ(sqrt_1_mod_r)
    p2_mod_r = p**2 % r
    assert (p2_mod_r == sqrt_1_mod_r) or (p2_mod_r == r-sqrt_1_mod_r)
    M = Matrix(ZZ, 2, 2, [r, 0, -p2_mod_r, 1])
    R = M.LLL()
    lambda_1 = R[0][0]
    lambda_2 = R[0][1]
    assert ((lambda_1 + lambda_2*p2_mod_r) % r) == 0
    cost_pairing_cp8_cp12(r, a0, a1, a2, a3, lambda_1, lambda_2, e0, e1, e2, e3, k)
    t_1 = tr0 - 1
    print("tests with D-twist")

    test_miller_function_ate(E, E_Fq4D, ED, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_2naf(E, E_Fq4D, ED, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_csb(E, E_Fq4D, ED, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_cln_b0(E, E_Fq4D, ED, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_2naf_cln_b0(E, E_Fq4D, ED, r, c, c2, t_1, D_twist=True)

    test_sparse_sparse_mult_d4_twist(Fq4D)
    test_sparse_mult_d4_twist(Fq4D)

    test_double_line_ate_cln_b0(E, ED, Fq4D, D_twist=True)
    test_add_line_ate_cln_b0(E, ED, Fq4D, D_twist=True)

    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist, E, ED, Fq4D, r, c, c2, t_1, D_twist=True, verbose=False)
    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist_acc_line, E, ED, Fq4D, r, c, c2, t_1, D_twist=True, verbose=False)

    print("tests with M-twist")

    test_miller_function_ate(E, E_Fq4M, EM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_2naf(E, E_Fq4M, EM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_csb(E, E_Fq4M, EM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_cln_b0(E, E_Fq4M, EM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_2naf_cln_b0(E, E_Fq4M, EM, r, c, c2, t_1, D_twist=False)

    test_sparse_sparse_mult_m4_twist(Fq4M)
    test_sparse_mult_m4_twist(Fq4M)

    test_double_line_ate_cln_b0(E, EM, Fq4M, D_twist=False)
    test_add_line_ate_cln_b0(E, EM, Fq4M, D_twist=False)

    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist, E, EM, Fq4M, r, c, c2, t_1, D_twist=False, verbose=False)
    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist_acc_line, E, EM, Fq4M, r, c, c2, t_1, D_twist=False, verbose=False)

    test_miller_function_ate_cln_b0_4_multi_scalar(E, ED, E_Fq4D, E_Fp8D, p, r, tr, c, c2, miller_scalars, map_Fq4D_Fp8D, D_twist=True, function_name=miller_function_ate_cln_b0_4_parallel)
    test_miller_function_ate_cln_b0_4_multi_scalar(E, EM, E_Fq4M, E_Fp8M, p, r, tr, c, c2, miller_scalars, map_Fq4M_Fp8M, D_twist=False, function_name=miller_function_ate_cln_b0_4_parallel)
    test_miller_function_ate_cln_b0_4_multi_scalar(E, ED, E_Fq4D, E_Fp8D, p, r, tr, c, c2, miller_scalars, map_Fq4D_Fp8D, D_twist=True, function_name=miller_function_ate_cln_b0_4_parallel_alt)
    test_miller_function_ate_cln_b0_4_multi_scalar(E, EM, E_Fq4M, E_Fp8M, p, r, tr, c, c2, miller_scalars, map_Fq4M_Fp8M, D_twist=False, function_name=miller_function_ate_cln_b0_4_parallel_alt)
    test_miller_function_ate_cln_b0_4_multi_scalar(E, ED, E_Fq4D, E_Fp8D, p, r, tr, c, c2, miller_scalars, map_Fq4D_Fp8D, D_twist=True, function_name=miller_function_ate_cln_b0_4_multi_scalar)
    test_miller_function_ate_cln_b0_4_multi_scalar(E, EM, E_Fq4M, E_Fp8M, p, r, tr, c, c2, miller_scalars, map_Fq4M_Fp8M, D_twist=False, function_name=miller_function_ate_cln_b0_4_multi_scalar)

    print("l1 + l2*i = 0 mod r, {}, {} bits, l1 = {:#x}, l2 = {:#x}".format(lambda_1.nbits(), lambda_2.nbits(), lambda_1, lambda_2))
    # compute omega pod p such that omega^2 = -1
    omega = Fp(t)/Fp(y)
    P = c*E.random_element()
    while P == E(0):
        P = c*E.random_element()
    psiP = E([-P[0], P[1]*omega])
    lambdaP = p2_mod_r * P
    if psiP[0] == lambdaP[0] and psiP[1] == -lambdaP[1]:
        omega = -omega
    elif psiP[0] != lambdaP[0] or psiP[1] != lambdaP[1]:
        raise ValueError("Error eigenvalue of endomorphism on G1")

    miller_scalars = [lambda_1, lambda_2, omega]
    test_miller_function_tate_cln_b0_2_multi_scalar(E, ED, E_Fq4D, E_Fp8D, p, r, tr, c, c2, miller_scalars, map_Fq4D_Fp8D, D_twist=True, function_name=miller_function_tate_cln_b0_2_parallel)
    test_miller_function_tate_cln_b0_2_multi_scalar(E, EM, E_Fq4M, E_Fp8M, p, r, tr, c, c2, miller_scalars, map_Fq4M_Fp8M, D_twist=False, function_name=miller_function_tate_cln_b0_2_parallel)

    test_miller_function_tate_cln_b0_2_multi_scalar(E, ED, E_Fq4D, E_Fp8D, p, r, tr, c, c2, miller_scalars, map_Fq4D_Fp8D, D_twist=True, function_name=miller_function_tate_cln_b0_2_parallel_alt)
    test_miller_function_tate_cln_b0_2_multi_scalar(E, EM, E_Fq4M, E_Fp8M, p, r, tr, c, c2, miller_scalars, map_Fq4M_Fp8M, D_twist=False, function_name=miller_function_tate_cln_b0_2_parallel_alt)

    test_miller_function_tate_cln_b0_2_multi_scalar(E, ED, E_Fq4D, E_Fp8D, p, r, tr, c, c2, miller_scalars, map_Fq4D_Fp8D, D_twist=True, function_name=miller_function_tate_cln_b0_2_multi_scalar)
    test_miller_function_tate_cln_b0_2_multi_scalar(E, EM, E_Fq4M, E_Fp8M, p, r, tr, c, c2, miller_scalars, map_Fq4M_Fp8M, D_twist=False, function_name=miller_function_tate_cln_b0_2_multi_scalar)

if __name__ == "__main__":
    arithmetic(False)
    print("BLS12-377")
    for v in testvector_bls12_377_cp8d1_768:
        if v['pnbits'] <= 764:
            test_curve(v)
    print("##################################################")
