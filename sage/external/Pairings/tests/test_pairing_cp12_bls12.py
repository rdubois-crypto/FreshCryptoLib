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
from external.Pairings.pairing_cocks_pinch import *
from external.Pairings.tests.test_pairing import *

from external.Pairings.tests.testvector_bls12_377_cp12d3_768 import testvector_bls12_377_cp12d3_768

from external.Pairings.tests.cost_pairing import cost_pairing_cp8_cp12, cost_tate_pairing_cp12_disc3

def test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, E2, E_Fqd, E_Fpk, p, r, tr, c, c2, miller_scalars, map_Fqd_Fpk, D_twist=True, function_name=miller_function_ate_a0_twist6_aklgl_4_multi_scalar):
    # test miller_function_ate_a0_twist6_aklgl_4_multi_scalar

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
    Fq = E2.base_field()
    if Fq.degree() == 1:
        xi = -Fqd.polynomial().constant_coefficient()
    else:
        xi = -Fqd.modulus().constant_coefficient()
    # precompute xi^((p-1)/6) and powers 1, 2, 3, 4, 5
    p = E.base_field().characteristic()
    xi1p = xi**((p-1)//6)
    xi2p = xi1p**2
    xi3p = xi1p * xi2p
    #xi4p = xi2p**2
    #xi5p = xi2p*xi3p
    xi1p2 = xi1p.frobenius() * xi1p # power (p+1)
    xi2p2 = xi1p2**2
    xi3p2 = xi1p2 * xi2p2
    #xi4p2 = xi2p2**2
    #xi5p2 = xi2p2*xi3p2
    xi1p3 = xi1p.frobenius(2) * xi1p2 # power (p^2+p+1)
    xi2p3 = xi1p3**2
    xi3p3 = xi1p3 * xi2p3
    #xi4p3 = xi2p3**2
    #xi5p3 = xi2p3*xi3p3
    if D_twist:
        pQ = E2([Q[0].frobenius() * xi2p, Q[1].frobenius() * xi3p])
        p2Q = E2([Q[0].frobenius(2) * xi2p2, Q[1].frobenius(2) * xi3p2])
        p3Q = E2([Q[0].frobenius(3) * xi2p3, Q[1].frobenius(3) * xi3p3])
    else:
        pQ = E2([Q[0].frobenius() / xi2p, Q[1].frobenius() / xi3p])
        p2Q = E2([Q[0].frobenius(2) / xi2p2, Q[1].frobenius(2) / xi3p2])
        p3Q = E2([Q[0].frobenius(3) / xi2p3, Q[1].frobenius(3) / xi3p3])
        
    a0, a1, a2, a3 = miller_scalars

    f, scalarsQ = function_name(Q, P, E2.a6(), miller_scalars, Fqd, xi2p, xi3p, xi2p2, xi3p2, D_twist=D_twist)
    ok = f != Fpk(1)
    if not ok:
        print("test bilinearity {}: {} pairing result is 1".format(function_name.__name__, ok))
        return False
    #print("type(scalarsQ) = {}".format(type(scalarsQ)))
    check_ai_piQ = ((type(scalarsQ) is list) or (type(scalarsQ) is tuple)) and len(scalarsQ) == 4
    if check_ai_piQ:
        a0Q, a1Q, a2Q, a3Q = scalarsQ
        # aklgl uses homogeneous projective coordinates (x,y) = (X/Z,Y/Z)
        a0Q = (a0Q[0]/a0Q[2], a0Q[1]/a0Q[2])
        a1Q = (a1Q[0]/a1Q[2], a1Q[1]/a1Q[2])
        a2Q = (a2Q[0]/a2Q[2], a2Q[1]/a2Q[2])
        a3Q = (a3Q[0]/a3Q[2], a3Q[1]/a3Q[2])
        ok_a0Q = E2(a0Q) == a0*Q
        ok_a1Q = E2(a1Q) == a1*pQ
        ok_a2Q = E2(a2Q) == a2*p2Q
        ok_a3Q = E2(a3Q) == a3*p3Q
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
            if D_twist:
                pQb = E2([Qb[0].frobenius() * xi2p, Qb[1].frobenius() * xi3p])
                p2Qb = E2([Qb[0].frobenius(2) * xi2p2, Qb[1].frobenius(2) * xi3p2])
                p3Qb = E2([Qb[0].frobenius(3) * xi2p3, Qb[1].frobenius(3) * xi3p3])
            else:
                pQb = E2([Qb[0].frobenius() / xi2p, Qb[1].frobenius() / xi3p])
                p2Qb = E2([Qb[0].frobenius(2) / xi2p2, Qb[1].frobenius(2) / xi3p2])
                p3Qb = E2([Qb[0].frobenius(3) / xi2p3, Qb[1].frobenius(3) / xi3p3])
            fab, scalars_abQ = function_name(Qb, Pa, E2.a6(), miller_scalars, Fqd, xi2p, xi3p, xi2p2, xi3p2, D_twist=D_twist)
            gab = fab**exponent
            ok = ok and gab != Fpk(1) and gab == g**(aa*bb)
            if check_ai_piQ:
                a0Q, a1Q, a2Q, a3Q = scalars_abQ
                a0Q = (a0Q[0]/a0Q[2], a0Q[1]/a0Q[2])
                a1Q = (a1Q[0]/a1Q[2], a1Q[1]/a1Q[2])
                a2Q = (a2Q[0]/a2Q[2], a2Q[1]/a2Q[2])
                a3Q = (a3Q[0]/a3Q[2], a3Q[1]/a3Q[2])
                ok_a0Q = E2(a0Q) == a0*Qb
                ok_a1Q = E2(a1Q) == a1*pQb
                ok_a2Q = E2(a2Q) == a2*p2Qb
                ok_a3Q = E2(a3Q) == a3*p3Qb
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

def test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, E2, E_Fqd, E_Fpk, u, p, r, tr, c, c2, omega, map_Fqd_Fpk, k_bls=12, D_twist=True, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar):
    # test miller_function_tate_a0_twist6_aklgl_2_multi_scalar

    P = c*E.random_element()
    while P == E(0):
        P = c * E.random_element()
    Q = c2*E2.random_element()
    while Q == E2(0):
        Q = c2 * E2.random_element()
    Fqd = E_Fqd.base_field() # assume that the degree d extension is explicit
    Fpk = E_Fpk.base_field()
    w = Fqd.gen(0)
    if Fqd.base_ring().degree() == 1:
        xi = -Fqd.polynomial().constant_coefficient()
    else:
        xi = -Fqd.modulus().constant_coefficient()
    exponent = (Fqd.cardinality()-1) // r # (q^d-1)//r = (p^k-1)//r

    #a0, a1, omega = miller_scalars
    k = Fpk.degree()
    assert k % 3 == 0
    d = gcd(k, 6)
    e = k//d
    tre_mod_r = (tr-1)**e % r
    if r-tre_mod_r < tre_mod_r:
        tre_mod_r = tre_mod_r - r
    if k_bls == 12 and k == 12 and e == 2:
        print("(tr-1)^2 mod r == u^5 - 3*u^4 + 3*u^3 - u + 2: {}".format(tre_mod_r == u**5 - 3*u**4 + 3*u**3 - u + 2))
        print("(tr-1)^2 mod r == -u^5 + 3*u^4 - 3*u^3 + u - 1: {}".format(tre_mod_r == -u**5 + 3*u**4 - 3*u**3 + u - 1))
    elif k_bls == 24 and k == 12 and e == 2:
        print("(tr-1)^2 mod r == u^9-3*u^8+4*u^7-4*u^6+3*u^5-2*u^3+2*u^2-u+2: {}".format(tre_mod_r == u**9-3*u**8+4*u**7-4*u**6+3*u**5-2*u**3+2*u**2-u+2))
        print("(tr-1)^2 mod r == -u^9+3*u^8-4*u^7+4*u^6-3*u^5+2*u^3-2*u^2+u-1: {}".format(tre_mod_r == -u**9+3*u**8-4*u**7+4*u**6-3*u**5+2*u**3-2*u**2+u-1))
    tre_mod_r_mod_u = tre_mod_r % u
    print("tre_mod_r = {}".format(tre_mod_r))
    print("tre_mod_r_mod_u = {}".format(tre_mod_r_mod_u))
    if abs(tre_mod_r_mod_u + u) < abs(tre_mod_r_mod_u):
        tre_mod_r_mod_u = tre_mod_r_mod_u + u
    elif abs(tre_mod_r_mod_u - u) < abs(tre_mod_r_mod_u):
        tre_mod_r_mod_u = tre_mod_r_mod_u - u
    print("tre_mod_r_mod_u = {}".format(tre_mod_r_mod_u))
    if d == 6:
        assert ((tre_mod_r**2 - tre_mod_r + 1) % r) == 0
    else:
        assert ((tre_mod_r**2 + tre_mod_r + 1) % r) == 0
    if (k_bls == 12 or k_bls==24) and k == 12 and e == 2:
        print("(tr-1)^2 mod r mod u == 2: {}".format(tre_mod_r_mod_u == 2))
        print("(tr-1)^2 mod r mod u == -1: {}".format(tre_mod_r_mod_u == - 1))
    T_mod_r_mod_u = tre_mod_r_mod_u + 1
    print("k = {} d={} e={} tre_mod_r_mod_u = {} T_mod_r_mod_u = {}".format(k, d, e, tre_mod_r_mod_u, T_mod_r_mod_u))
    assert ((k % 6) == 0 and (T_mod_r_mod_u == 0 or T_mod_r_mod_u == 3)) or ((k % 6) == 3 and (T_mod_r_mod_u == 0 or T_mod_r_mod_u == -3))

    Fp = Fpk.base_ring()
    omega = Fp(omega)
    assert omega**2 + omega + 1 == 0
    psiP = E([omega * P[0], P[1]])
    # TODO check that psiP has the expected eigenvalue mod r!!!
    # if parallel, expected eigenvalue is p^4 mod r.
    # is multi-scalar, expected eigenvalue is (-p^2) mod r.
    if d == 6:
        if function_name is miller_function_tate_a0_twist6_aklgl_2_multi_scalar:
            lambda_mod_r = -tre_mod_r
        elif function_name is miller_function_tate_aklgl_a0_2_parallel:
            lambda_mod_r = tre_mod_r-1
        else: # default
            lambda_mod_r = -tre_mod_r
    else:
        lambda_mod_r = tre_mod_r
    check_psiP = lambda_mod_r*P
    # omega -> -omega-1?
    if check_psiP[0] == -psiP[0]-P[0] and check_psiP[1] == psiP[1]:
        omega = -omega - 1
        print("adjusting omega -> -omega-1")
    elif check_psiP[0] != psiP[0] or check_psiP[1] != psiP[1]:
        print("Problem psi(P) = [omega*x, y] has not eigenvalue (tr-1)^{} mod r".format((k//3)))
    if k_bls == 12 and d == 6:
        if T_mod_r_mod_u == 0:
            # case corresponding to t = 0 mod r mod u for BLS12-BW6 curves
            list_a1 = (-(u+1), u**3-u**2+1, "(-(u+1), u^3-u^2+1)")
            list_a2 = (u**3-u**2-u, u+1, "(u^3-u^2-u, u+1)")
        elif T_mod_r_mod_u == 3:
            # case corresponding to t = 3 mod r mod u for BLS12-BW6 curves
            list_a1 = (u+1, u**3-u**2-u, "(u+1, u^3-u^2-u)")
            list_a2 = (u**3-u**2+1, -(u+1), "(u^3-u^2+1, -(u+1))")
    elif k_bls == 12 and d == 3: # cubic twist instead of sextic twist, T is replaced by -T in a0 + a1*T = 0 mod r -> change one of ai to -ai
        if T_mod_r_mod_u == 0:
            list_a1 = (u+1, u**3-u**2+1, "(u+1, u^3-u^2+1)")
            list_a2 = (u**3-u**2-u, -(u+1), "(u^3-u^2-u, -(u+1))")
        elif T_mod_r_mod_u == -3:
            list_a1 = (u+1, u**3-u**2-u, "(u+1, u^3-u^2-u)")
            list_a2 = (u**3-u**2+1, -(u+1), "(u^3-u^2+1, -(u+1))")
    elif k_bls == 24 and d == 6:
        if T_mod_r_mod_u == 0:
            # case corresponding to t = 0 mod r mod u for BLS24-BW6 curves
            list_a1 = (-(u+1), u**5-u**4+1, "(-(u+1), u^5-u^4+1)")
            list_a2 = (u**5-u**4-u, u+1, "(u^5-u^4-u, u+1)")
        elif T_mod_r_mod_u == 3:
            # case corresponding to t = 3 mod r mod u for BLS24-BW6 curves
            list_a1 = (u+1, u**5-u**4-u, "(u+1, u^5-u^4-u)")
            list_a2 = (u**5-u**4+1, -(u+1), "(u^5-u^4+1, -(u+1))")
    elif k_bls == 24 and d == 3:
        if T_mod_r_mod_u == 0:
            list_a1 = (u+1, u**5-u**4+1, "(u+1, u^5-u^4+1)")
            list_a2 = (u**5-u**4-u, -(u+1), "(u^5-u^4-u, -(u+1))")
        elif T_mod_r_mod_u == -3:
            list_a1 = (-(u+1), u**5-u**4-u, "(-(u+1), u^5-u^4-u)")
            list_a2 = (u**5-u**4+1, -u+1, "(u^5-u^4+1, u+1)")

    all_ok = True
    for (a0, a1, str_a) in [list_a1, list_a2]:
        assert ((a0 + a1*tre_mod_r) % r) == 0
        ok = True
        if function_name is miller_function_tate_a0_twist6_aklgl_2_multi_scalar:
            a1 = -a1
            f, scalarsP = function_name(P, Q, E.a6(), a0, a1, omega, xi, Fqd, D_twist)
        elif function_name is miller_function_tate_aklgl_a0_2_parallel:
            f, scalarsP = function_name(P, Q, E.a6(), a0, a1, Fqd, map_Fqd_Fpk, D_twist)
            #f, scalarsP = function_name(P, Q, E.a6(), a0, a1, Fqd, map_Fqd_Fpk, D_twist, omega, xi)
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
            # should be zero?
            ok_ai_P = scalarsP[2] == 0
        g = f**exponent
        aa = 1
        bb = 1
        while ok and aa < 4:
            Pa = aa*P
            psiPa = E([Pa[0]*omega, Pa[1]])
            bb = 1
            while ok and bb < 4:
                Qb = bb*Q
                if function_name is miller_function_tate_a0_twist6_aklgl_2_multi_scalar:
                    fab, scalars_abP = function_name(Pa, Qb, E.a6(), a0, a1, omega, xi, Fqd, D_twist)
                elif function_name is miller_function_tate_aklgl_a0_2_parallel:
                    fab, scalars_abP = function_name(Pa, Qb, E.a6(), a0, a1, Fqd, map_Fqd_Fpk, D_twist)
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
        print("test bilinearity {} scalars {}: {} ({} tests)".format(function_name.__name__, str_a, ok, (aa-1)*(bb-1)))
        if check_ai_P:
            print("a0P == [a0]P, a1P == [a1](P): {}".format(ok_a0P and ok_a1P))
        else:
            print("a0*P + a1*psiP == 0: {}".format(ok_ai_P))
        all_ok = all_ok and ok
    return all_ok

def test_curve(test_vector_value, k_bls):
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
    b = v['b']
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
    # try with cubic twists
    print("gcd(y, (p+1-(3*y-t)//2)*(p+1-(-3*y-t)//2)) = {}".format( (gcd(y, (p+1-(3*y-t)//2)*(p+1-(-3*y-t)//2)).factor())))
    # try with sextic twists
    print("gcd(y, (p+1-(3*y+t)//2)*(p+1-(-3*y+t)//2)) = {}".format( (gcd(y, (p+1-(3*y+t)//2)*(p+1-(-3*y+t)//2)).factor())))
    print("--> isogenies of these degrees are defined over Fp.")
    
    curve_order = p+1-t
    assert curve_order == r*c
    Fp = GF(p, proof=False)
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    E = EllipticCurve([Fp(0), Fp(b)])
    b = Fp(b)
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
    # sextic twists are p^2 + 1 - (t2+3*y2)/2 = p^2 + p + 1 - t*(t+3*y)/2,
    #                   p^2 + 1 - (t2-3*y2)/2 = p^2 + p + 1 - t*(t-3*y)/2.
    assert (tr2+3*y2) % 2 == 0
    assert (tr2-3*y2) % 2 == 0
    if ((p2+1-(tr2+3*y2)//2) % r) == 0:
        twist_order = p2+1-(tr2+3*y2)//2
        g2c = twist_order // r
    elif ((p2+1-(tr2-3*y2)//2) % r) == 0:
        twist_order = p2+1-(tr2-3*y2)//2
        g2c = twist_order // r
    else:
        raise ValueError("Error no appropriate sextic twist")
    
    xiD, bD = find_twist_curve_parameter_xi_ab(b, Fq, r, g2c, d=6, D_twist=True)
    # where bD = b/xiD
    
    c2 = g2c

    print("xiD = {}".format(xiD))
    ED = EllipticCurve([Fq(0), bD])
    Q = ED.random_element()
    orderD = g2c * r
    assert orderD * Q == ED(0)
    print("ED ok")
    Fq6D = Fq.extension(Z**6 - xiD, names=('wD',)); (wD,) = Fq6D._first_ngens(1)
    E_Fq6D = EllipticCurve([Fq6D(0), Fq6D(b)])
    # absolute extension
    # Z**6 == xiD
    # xiD has the form v0 + v1*v in Fq
    v0D, v1D = xiD.polynomial().list()
    assert xiD == v0D + v1D*v
    # (Z**6 - v0) == v1*v
    # (Z**6 - v0)**2 == v1**2 * v**2 = v1**2 * (-a0)
    # (Z**6 - v0)**2 + a0 * v1**2
    Fp12D = Fp.extension((z**6 - v0D)**2 + a0*v1D**2, names=('vD',)); (vD,) = Fp12D._first_ngens(1)
    E_Fp12D = EllipticCurve([Fp12D(0), Fp12D(b)])

    def map_Fq6D_Fp12D(X, aD=None):
        if aD is None:
            aD = vD
        return sum([xi.polynomial()((aD**6-v0D)/v1D) * aD**e for e,xi in enumerate(X.list())])

    def map_ED_E_Fq6D(Q2):
        return E_Fq6D(psi_sextic_d_twist(Q2, wD))

    # now find a 6-th M-twist
    xiM, bM = find_twist_curve_parameter_xi_ab(b, Fq, r, g2c, d=6, D_twist=False)
    # where bM = b*xiM
    print("xiM = {}".format(xiM))
    EM = EllipticCurve([Fq(0), bM])
    Q = EM.random_element()
    orderM = g2c * r
    assert orderM * Q == EM(0)
    print("EM ok")
    Fq6M = Fq.extension(Z**6 - xiM, names=('wM',)); (wM,) = Fq6M._first_ngens(1)
    E_Fq6M = EllipticCurve([Fq6M(0), Fq6M(b)])
    # absolute extension
    # Z**6 == xiM
    # xiM has the form v0 + v1*v in Fq
    v0M, v1M = xiM.polynomial().list()
    assert xiM == v0M + v1M*v
    # (Z**6 - v0) == v1*v
    # (Z**6 - v0)**2 == v1**2 * v**2 = v1**2 * (-a0)
    # (Z**6 - v0)**2 + a0 * v1**2
    Fp12M = Fp.extension((z**6 - v0M)**2 + a0*v1M**2, names=('vM',)); (vM,) = Fp12M._first_ngens(1)
    E_Fp12M = EllipticCurve([Fp12M(0), Fp12M(b)])
    
    def map_Fq6M_Fp12M(X, aM=None):
        if aM is None:
            aM = vM
        return sum([xi.polynomial()((aM**6-v0M)/v1M) * aM**e for e,xi in enumerate(X.list())])

    def map_EM_E_Fq6M(Q2):
        return E_Fq6M(psi_sextic_m_twist(Q2, wM))

    a0, a1, a2, a3 = miller_scalars
    #final_exp_hard
    ee = (p**4-p**2+1)//r
    M = Matrix(ZZ, 4, 4, [p,0,0,0, ee,1,0,0, 0,ee,1,0, 0,0,ee,1])
    R = M.LLL()
    e0, e1, e2, e3 = (R[0]).list()
    # what would be a Tate pairing? There is an endomorphism on G1, of eigenvalue zeta_6 = (1+sqrt(-3))/2 mod r.
    # if M-twist, E': y^2 = x^3 + b*beta, this is (x,y) -> (x^(p^2)*w^2^(p^2)/w^2, y^(p^2)*w^3^(p^2)/w^3)
    # where w^(p^2) = w^(p+1)*(p-1) * w = (w^6)^((p+1)*(p-1)/6) * w = Trace(beta)^(p-1)/6 * w
    # finally, we have (x * Trace(beta)^(p-1)/3, -y) because Trace(beta)^(p-1)/2 = -1 (it cannot be 1).
    #                  (x*beta^((p+1)*(p-1)/3), -y), eigenvalue p^2 mod r = (t-1)^2 mod r
    # if D-twist, E': y^2 = x^3 + b/beta, this is (x,y) -> (x/beta^((p+1)*(p-1)/3), -y), eigenvalue p^2 mod r = (t-1)^2 mod r
    Fr = GF(r, proof=False)
    sqrt_3_mod_r = Fr(tr0-2) / Fr(y0)
    assert sqrt_3_mod_r**2 + 3 == 0
    eigenvalue3 = (-1+sqrt_3_mod_r)/2
    assert ((eigenvalue3**2 + eigenvalue3 + 1) % r) == 0
    eigenvalue6 = (1+sqrt_3_mod_r)/2
    assert ((eigenvalue6**2 - eigenvalue6 + 1) % r) == 0
    eigenvalue3_mod_r = ZZ(eigenvalue3)
    eigenvalue6_mod_r = ZZ(eigenvalue6)
    # for k == 12, p^4 - p^2 + 1 = 0 mod r => p^4 = p^2 - 1 mod r
    p2_mod_r = p**2 % r
    p4_mod_r = (p2_mod_r-1) % r
    assert (p4_mod_r == eigenvalue3_mod_r) or (p4_mod_r == r-eigenvalue3_mod_r-1)
    assert (p2_mod_r == eigenvalue6_mod_r) or (p2_mod_r == r-eigenvalue6_mod_r+1)
    if p4_mod_r == r-eigenvalue3_mod_r-1:
        eigenvalue3_mod_r = r-eigenvalue3_mod_r-1
    if p2_mod_r == r-eigenvalue6_mod_r+1:
        eigenvalue6_mod_r = r-eigenvalue6_mod_r+1
    t2 = p2_mod_r
    if r-t2 < t2:
        t2 = t2-r
    t2_mod_r_mod_u = t2 % u
    # and what if u < 0?
    if abs(t2_mod_r_mod_u-u) < abs(t2_mod_r_mod_u):
        t2_mod_r_mod_u = t2_mod_r_mod_u - u
    elif abs(t2_mod_r_mod_u+u) < abs(t2_mod_r_mod_u):
        t2_mod_r_mod_u = t2_mod_r_mod_u + u
    print("k_bls = {} t2_mod_r_mod_u = {}".format(k_bls, t2_mod_r_mod_u))
    if k_bls == 12 and t2_mod_r_mod_u == -1: # like t_mod_r_mod_u == 0 in BLS12-BW6 case
        print("(t-1)^2 mod r = (-u^5 + 3*u^4 - 3*u^3 + u) - 1: {}".format(t2 == -u**5 + 3*u**4 - 3*u**3 + u - 1))
        print("-(u+1) + (u^3-u^2+1)*(t-1)**2 = {} mod r".format((-(u+1) + (u**3-u**2+1)*t2) % r))
        print("u^3-u^2-u + (u+1)*(t-1)**2 = {} mod r".format((u**3-u**2-u + (u+1)*t2) % r))
        assert (-(u+1) + (u**3-u**2+1)*t2) % r == 0
        assert (u**3-u**2-u + (u+1)*t2) % r == 0
        lambda_1 = (u+1)
        lambda_2 = (u**3-u**2+1)
    elif k_bls == 12 and t2_mod_r_mod_u == 2: # like t_mod_r_mod_u == 3 in BLS12-BW6 case
        print("(t-1)^2 mod r = (u^5 - 3*u^4 + 3*u^3 - u + 3) - 1: {}".format(t2 == u**5 - 3*u**4 + 3*u**3 - u + 3 - 1))
        print("(u+1) + (u^3-u^2-u)*(t-1)**2 = {} mod r".format((u+1 + (u**3-u**2-u)*t2) % r))
        print("(u^3-u^2+1) - (u+1)*(t-1)**2 = {} mod r".format((u**3-u**2+1 - (u+1)*t2) % r))
        assert (u+1 + (u**3-u**2-u)*t2) % r == 0
        assert (u**3-u**2+1 - (u+1)*t2) % r == 0
        lambda_1 = -(u+1)
        lambda_2 = (u**3-u**2-u)
    elif k_bls == 24 and t2_mod_r_mod_u == -1: # like t_mod_r_mod_u == 0 in BLS24-BW6 case
        print("(t-1)^2 mod r = -u^9+3*u^8-4*u^7+4*u^6-3*u^5+2*u^3-2*u^2+u-1: {}".format(t2 == -u**9+3*u**8-4*u**7+4*u**6-3*u**5+2*u**3-2*u**2+u-1))
        print("-(u+1) + (u^5-u^4+1)*(t-1)**2 = {} mod r".format((-(u+1) + (u**5-u**4+1)*t2) % r))
        print("u^5-u^4-u + (u+1)*(t-1)**2 = {} mod r".format((u**5-u**4-u + (u+1)*t2) % r))
        assert (-(u+1) + (u**5-u**4+1)*t2) % r == 0
        assert (u**5-u**4-u + (u+1)*t2) % r == 0
        lambda_1 = (u+1)
        lambda_2 = (u**5-u**4+1)
        lambda_3 = (u**5-u**4-u)
        bits_2naf_l2 = bits_2naf(abs(lambda_2))
        bits_2naf_l3 = bits_2naf(abs(lambda_3))
        hw_bits_2naf_l2 = sum([1 for bi in bits_2naf_l2 if bi != 0])
        hw_bits_2naf_l3 = sum([1 for bi in bits_2naf_l3 if bi != 0])
        print("hw_bits_2naf(u^5-u^4+1) = {} hw_bits_2naf(u^5-u^4-u) = {}".format(hw_bits_2naf_l2, hw_bits_2naf_l3))
        if hw_bits_2naf_l2 > hw_bits_2naf_l3:
            lambda_2 = lambda_3
    elif k_bls == 24 and t2_mod_r_mod_u == 2: # like t_mod_r_mod_u == 3 in BLS24-BW6 case
        print("(tr-1)^2 mod r = u^9-3*u^8+4*u^7-4*u^6+3*u^5-2*u^3+2*u^2-u+2: {}".format(t2 == u**9-3*u**8+4*u**7-4*u**6+3*u**5-2*u**3+2*u**2-u+2))
        print("(u+1) + (u^5-u^4-u)*(t-1)**2 = {} mod r".format((u+1 + (u**5-u**4-u)*t2) % r))
        print("(u^5-u^4+1) - (u+1)*(t-1)**2 = {} mod r".format((u**5-u**4+1 - (u+1)*t2) % r))
        assert (u+1 + (u**5-u**4-u)*t2) % r == 0
        assert (u**5-u**4+1 - (u+1)*t2) % r == 0
        lambda_1 = -(u+1)
        lambda_2 = (u**5-u**4-u)
        lambda_3 = (u**5-u**4+1)
        bits_2naf_l2 = bits_2naf(abs(lambda_2))
        bits_2naf_l3 = bits_2naf(abs(lambda_3))
        hw_bits_2naf_l2 = sum([1 for bi in bits_2naf_l2 if bi != 0])
        hw_bits_2naf_l3 = sum([1 for bi in bits_2naf_l3 if bi != 0])
        print("hw_bits_2naf(u^5-u^4+1) = {} hw_bits_2naf(u^5-u^4-u) = {}".format(hw_bits_2naf_l3, hw_bits_2naf_l2))
        if hw_bits_2naf_l2 > hw_bits_2naf_l3:
            lambda_2 = lambda_3
    else:
        print("problem with t2_mod_r_mod_u = {}".format(t2_mod_r_mod_u))
        print("t2 = (t-1)^2 mod r = {}".format(t2))
        print("u^9-3*u^8+4*u^7-4*u^6+3*u^5-2*u^3+2*u^2-u+2 = {}".format(u**9-3*u**8+4*u**7-4*u**6+3*u**5-2*u**3+2*u**2-u+2))
        print("-u^9+3*u^8-4*u^7+4*u^6-3*u^5+2*u^3-2*u^2+u-1 = {}".format(-u**9+3*u**8-4*u**7+4*u**6-3*u**5+2*u**3-2*u**2+u-1))

    cost_pairing_cp8_cp12(r, a0, a1, a2, a3, lambda_1, lambda_2, e0, e1, e2, e3, k)
    if k_bls == 12:
        t_1_e_mod_r = p4_mod_r
    else:
        t_1_e_mod_r = p2_mod_r
    if r-t_1_e_mod_r < t_1_e_mod_r:
        t_1_e_mod_r = r-t_1_e_mod_r
    t_1_e_mod_r_mod_u = t_1_e_mod_r % u
    if abs(t_1_e_mod_r_mod_u-u) < t_1_e_mod_r_mod_u:
        t_1_e_mod_r_mod_u = t_1_e_mod_r_mod_u - u
    elif abs(t_1_e_mod_r_mod_u+u) < t_1_e_mod_r_mod_u:
        t_1_e_mod_r_mod_u = t_1_e_mod_r_mod_u + u
    cost_tate_pairing_cp12_disc3(r, u, k, k_bls, t_1_e_mod_r_mod_u)
    t_1 = tr0 - 1

    print("test sparse-dense and sparse-sparse multiplications")
    test_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_d6_twist(Fq6D)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_sparse_mult_d6_twist(Fq6D)

    print("tests with D-twist")
    test_miller_function_tate(E, E_Fq6D, ED, r, c, c2, D_twist=True)
    test_miller_function_tate_2naf(E, E_Fq6D, ED, r, c, c2, D_twist=True)

    assert ((t_1**4 - t_1**2 + 1) % r) == 0
    test_miller_function_ate(E, E_Fq6D, ED, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_2naf(E, E_Fq6D, ED, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_aklgl(E, ED, Fq6D, xiD, r, c, c2, t_1, D_twist=True)
    test_miller_function_ate_2naf_aklgl(E, ED, Fq6D, xiD, r, c, c2, t_1, D_twist=True)

    print("tests with M-twist")
    test_miller_function_tate(E, E_Fq6M, EM, r, c, c2, D_twist=False)
    test_miller_function_tate_2naf(E, E_Fq6M, EM, r, c, c2, D_twist=False)

    test_miller_function_ate(E, E_Fq6M, EM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_2naf(E, E_Fq6M, EM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_aklgl(E, EM, Fq6M, xiM, r, c, c2, t_1, D_twist=False)
    test_miller_function_ate_2naf_aklgl(E, EM, Fq6M, xiM, r, c, c2, t_1, D_twist=False)

    """
    test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, ED, E_Fq6D, E_Fp12D, p, r, tr, c, c2, miller_scalars, map_Fq6D_Fp12D, D_twist=True, function_name=miller_function_ate_a0_twist6_aklgl_4_parallel)
    test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, EM, E_Fq6M, E_Fp12M, p, r, tr, c, c2, miller_scalars, map_Fq6M_Fp12M, D_twist=False, function_name=miller_function_ate_a0_twist6_aklgl_4_parallel)
    test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, ED, E_Fq6D, E_Fp12D, p, r, tr, c, c2, miller_scalars, map_Fq6D_Fp12D, D_twist=True, function_name=miller_function_ate_a0_twist6_aklgl_4_parallel_alt)
    test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, EM, E_Fq6M, E_Fp12M, p, r, tr, c, c2, miller_scalars, map_Fq6M_Fp12M, D_twist=False, function_name=miller_function_ate_a0_twist6_aklgl_4_parallel_alt)
    """

    test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, ED, E_Fq6D, E_Fp12D, p, r, tr, c, c2, miller_scalars, map_Fq6D_Fp12D, D_twist=True, function_name=miller_function_ate_a0_twist6_aklgl_4_multi_scalar)
    test_miller_function_ate_a0_twist6_aklgl_4_multi_scalar(E, EM, E_Fq6M, E_Fp12M, p, r, tr, c, c2, miller_scalars, map_Fq6M_Fp12M, D_twist=False, function_name=miller_function_ate_a0_twist6_aklgl_4_multi_scalar)

    # compute omega mod p such that omega^2 + omega + 1 = 0
    omega = Fp(t)/Fp(y)
    assert omega**2 + 3 == 0
    omega = (-1 + omega)/2
    assert omega**2 + omega + 1 == 0
    P = c*E.random_element()
    while P == E(0):
        P = c*E.random_element()
    psiP = E([omega*P[0], P[1]])
    # be careful: the scalars for BW6 assume that a0 + a1*(tr_bw-1) = 0 mod r, and
    # either tr_w-1 = -x^5 ... - 1 or +x^5 ... + 2 for BW6-BLS12,
    # or              -x^9 ... -1, or +x^9 ... + 2 for BW6-BLS24.
    # then re-using the same (a0,a1) here needs an adjustment.
    # because psi has eigenvalue either T or -T-1 with T = -(tr_bw-1) the one for BW6 curves...
    # for CP12 curves, p = tr-1 mod r on one side, (p^4-p^2+1) = 0 mod r on the other side,
    # hence (tr_cp12 - 1)^4 <-> p^4 = p^2-1
    # (tr_cp12 - 1)^2 = p^2 corresponds to one of (tr_bw0-1), (tr_bw3-1) = -(tr_bw0-1)+1
    # but then psi has eigenvalue -(tr_bw0-1) or (tr_bw0-1)-1.
    assert ((p2_mod_r**2 - p2_mod_r + 1) % r) == 0
    assert ((p4_mod_r**2 + p4_mod_r + 1) % r) == 0
    lambda_p2_P = (-p2_mod_r) * P
    lambdap4_P = (p4_mod_r) * P
    if psiP[0] == lambda_p2_P[0] and psiP[1] == lambda_p2_P[1]:
        print("psi(P) = (omega*xP, yP) of eigenvalue (-p^2 mod r)")
        omega_p2 = omega
        omegap4 = -omega-1
        assert -P[0]-psiP[0] == lambdap4_P[0] and psiP[1] == lambdap4_P[1]
        print("psi'(P) = ((-omega-1)*xP, yP) of eigenvalue (p^4 mod r)")
    elif -P[0]-psiP[0] == lambda_p2_P[0] and psiP[1] == lambda_p2_P[1]:
        omega_p2 = -1-omega
        omegap4 = omega
        assert psiP[0] == lambdap4_P[0] and psiP[1] == lambdap4_P[1]
        print("omega->-omega-1, now:")
        print("psi(P) = (omega*xP, yP) of eigenvalue (-p^2 mod r)")
        print("psi'(P) = ((-omega-1)*xP, yP) of eigenvalue (p^4 mod r)")
    elif psiP[0] != lambdaP[0] or psiP[1] != lambdaP[1]:
        raise ValueError("Error eigenvalue of endomorphism on G1")

    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fp12D, u, p, r, tr, c, c2, omegap4, map_Fq6D_Fp12D, k_bls, D_twist=True, function_name=miller_function_tate_aklgl_a0_2_parallel)
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fp12M, u, p, r, tr, c, c2, omegap4, map_Fq6M_Fp12M, k_bls, D_twist=False, function_name=miller_function_tate_aklgl_a0_2_parallel)

    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fp12D, u, p, r, tr, c, c2, omega_p2, map_Fq6D_Fp12D, k_bls, D_twist=True, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar)
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fp12M, u, p, r, tr, c, c2, omega_p2, map_Fq6M_Fp12M, k_bls, D_twist=False, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar)


if __name__ == "__main__":
    arithmetic(False)
    print("BLS12-377")
    for v in test_vector_BLS12_377_CP12D3_768:
        if v['pnbits'] <= 764:
            test_curve(v, k_bls=12)
    print("##################################################")
