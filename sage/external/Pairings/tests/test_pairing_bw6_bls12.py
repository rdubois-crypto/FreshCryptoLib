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
from external.Pairings.pairing_bw6_bls12 import *
from external.Pairings.pairing_cocks_pinch import miller_function_tate_aklgl_a0_2_parallel, miller_function_tate_a0_twist6_aklgl_2_multi_scalar
from external.Pairings.tests.test_pairing import *
from external.Pairings.tests.test_pairing_bw6 import *
from external.Pairings.tests.test_pairing_bw6_761 import test_g2_frobenius_eigenvalue_bw6
from external.Pairings.tests.test_pairing_bw6_761 import cost_pairing_bw6_761
from external.Pairings.tests.test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761
from external.Pairings.tests.test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761_aklgl
from external.Pairings.tests.test_pairing_bw6_761 import test_miller_loop_opt_ate_bw6_761_all
from external.Pairings.tests.test_pairing_bw6_761 import test_miller_loop_opt_ate_bw6_761_aklgl_all
from external.Pairings.tests.test_pairing_bw6_761 import test_miller_loop_opt_ate_all_bw6_761_and_aklgl
from external.Pairings.tests.test_pairing_bw6_761 import test_formula_miller_loop_opt_ate_bw6_761_aklgl
from external.Pairings.tests.test_pairing_cp12_bls12 import test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar

from external.Pairings.tests.test_scalar_mult import test_glv_scalar_mult_g1

from external.Pairings.tests.testvector_bls12_377_bw6_768 import testvector_bls12_377_bw6_768
from external.Pairings.tests.testvector_bls12_379_bw6_768 import testvector_bls12_379_bw6_768
from external.Pairings.tests.testvector_bls12_381_bw6_768 import testvector_bls12_381_bw6_768

from external.Pairings.tests.cost_pairing import cost_final_exp_hard_bw6_bls12, cost_pairing_opt_tate_bw6_bls12, cost_pairing_bw6_bls12

def test_bw6_bls12_test_cofactor_formula(E, r, c, t, y, u, omega, ht, hy):
    """
    check 
    t mod r mod u = 0:
    cx = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t-2) + (u^2 - 3*u + 3)*(u^2 - u + 1) - hy
    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx
    l0 = (u^3-u^2+1)*(ht^2+3*hy^2) - 4*ht*(u - 1)^2 - 2*(ht-3*hy)
    l1 = (ht^2+3*hy^2)*(u + 1) - 2*(ht + 3*hy)*(u - 1)^2 -4*ht

    t mod r mod u = 3: the same but with (hy) instead of (-hy)
    cx = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t-2) + (u^2-3*u+3)*(u^2-u+1) + hy

    because q+1-t = r*c = ((t-2)^2 + 3*y^2)/4,
    sqrt_3 = (t-2)/y mod r*c
    eigenvalue = (-1 + (t-2)*inv_y)/2 % c;
    """
    eigenvalue_mod_c = None
    g, inv_y, inv_c = xgcd(y, c)
    assert g == y*inv_y + inv_c*c
    if abs(g) != 1:
        mod_c = c // abs(g)
    else:
        mod_c = c
    numer = (-1 + (t-2)*inv_y) % mod_c
    if (numer % 2) == 0:
        eigenvalue_mod_c = numer // 2
    elif (mod_c % 2) == 1:
        eigenvalue_mod_c = (numer + mod_c) // 2
    else:
        numer = (-1 - (t-2)*inv_y) % mod_c
        if (numer % 2) == 0:
            eigenvalue_mod_c = numer // 2
        elif (mod_c % 2) == 1:
            eigenvalue_mod_c = (numer + mod_c)//2
        else:
            print("error cannot compute the eigenvalue (-1+sqrt(-3))/2 mod c, gcd(y,c) = {}, val_2(y) = {}, val_2(c) = {}, c % 3 = {}, val_3(c) = {}, val_3(c-1) = ".format(g, valuation(y,2), valuation(c,2), c % 3, valuation(c,3), valuation(c-1,3)))
            print("c = {}\ny = {}\ngcd_cy = {}\n inv_y = {}\ninv_c = {}\n".format(c, y, g, inv_y, inv_c))
            print("inv_y*y + inv_c*c = {}".format(inv_y*y + inv_c*c))
            #raise ValueError("cannot compute the eigenvalue (-1+sqrt(-3))/2 mod c")
            eigenvalue_mod_c = None
    P = E.random_element()
    E0 = E(0)
    rP = r*P
    omeg_rP = bw6_phi(rP, omega)
    omeg_rP_ = bw6_phi(rP, -omega-1)
    omeg_P = bw6_phi(P, omega)
    omeg_P_ = bw6_phi(P, -omega-1)

    if eigenvalue_mod_c is not None:
        lamb_rP = eigenvalue_mod_c * rP
        lamb_rP_ = (-eigenvalue_mod_c-1) * rP
        print("testing eigenvalue of phi on cofactor subgroup")
        one_match_found = False
        if lamb_rP == omeg_rP:
            print("lamb_rP == bw6_phi(rP, omega): True")
            one_match_found = True
        if lamb_rP_ == omeg_rP_:
            print("lamb_rP_== bw6_phi(rP, -omega-1): True")
            one_match_found = True
        if lamb_rP == omeg_rP_:
            print("lamb_rP == bw6_phi(rP, -omega-1): True")
            one_match_found = True
        if lamb_rP_ == omeg_rP:
            print("lamb_rP_== bw6_phi(rP, omega): True")
            one_match_found = True
        if not one_match_found:
            print("No compatibility lamb_rP, bw6_phi(rP, omega) found. Result is False")

        if not (lamb_rP == omeg_rP and lamb_rP_ == omeg_rP_) and (lamb_rP == omeg_rP_ and lamb_rP_ == omeg_rP):
            eigenvalue_mod_c = -eigenvalue_mod_c-1

    t_mod_r = (t % r)
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r = t_mod_r - r
    t_mod_r_mod_u = t_mod_r % u
    if abs(t_mod_r_mod_u-u) < abs(t_mod_r_mod_u):
        t_mod_r_mod_u = t_mod_r_mod_u - u
    #l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx
    if t_mod_r_mod_u == 0:
        l0 = (u**3-u**2+1)*(ht**2+3*hy**2)//4 - ht*(u - 1)**2 - (ht-3*hy)//2
        l1 = (ht**2+3*hy**2)//4*(u + 1) - (ht + 3*hy)//2*(u - 1)**2 - ht
    elif t_mod_r_mod_u == 3:
        l0 = (u**3-u**2+1)*(ht**2+3*hy**2)//4 + ht + (ht+3*hy)//2*(u - 1)**2
        l1 = (ht**2+3*hy**2)//4*(u + 1) - (ht-3*hy)//2*(u**2-2*u+2) + ht
    else:
        print("Error in test_bw6_bls12_test_cofactor_formula: (t mod r) mod u = {} not 0 not 3".format(t_mod_r_mod_u))
        return
    Q1 = l0*rP + l1*omeg_rP
    Q1_ = l0*rP + l1*omeg_rP_
    Q2 = l0*P + l1*omeg_P
    Q2_ = l0*P + l1*omeg_P_
    print("l0*rP + l1*omeg_rP = {}".format(Q1))
    print("l0*rP + l1*omeg_rP_= {}".format(Q1_))
    print("l0*P + l1*omeg_P == 0: {}".format(Q2 == E0))
    print("l0*P + l1*omeg_P_== 0: {}".format(Q2_ == E0))
    print("r*(l0*P + l1*omeg_P) == 0: {}".format(r*Q2 == E0))
    print("r*(l0*P + l1*omeg_P_)== 0: {}".format(r*Q2_ == E0))
    
def test_bw6_bls12_g1_mult_by_cofactor(E, r, c, t, u, omega, ht, hy, verbose=True):
    t_mod_r = t % r
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r -= r
    t_mod_r_mod_u = t_mod_r % u
    if t_mod_r_mod_u == 0:
        function = bw6_bls12_g1_mult_by_cofactor_trace_0_mod_r_u
        function_alt = bw6_bls12_g1_mult_by_cofactor_trace_0_mod_r_u_alt
    else:
        function = bw6_bls12_g1_mult_by_cofactor_trace_3_mod_r_u
        function_alt = bw6_bls12_g1_mult_by_cofactor_trace_3_mod_r_u_alt
    E0 = E(0)
    ok = True
    ok_alt = True
    i = 0
    while (ok or ok_alt) and i < 10:
        P = E.random_element()
        rP = r*P
        assert rP != E0 and c*rP == E0 and c*P != E0
        if ok:
            cP = function(P, omega, u, ht, hy)
            crP = function(rP, omega, u, ht, hy)
            ok = cP != E0 and r*cP == E0 and crP == E0
        if ok_alt:
            cP_alt = function_alt(P, omega, u, ht, hy)
            crP_alt = function_alt(rP, omega, u, ht, hy)
            ok_alt = cP_alt != E0 and r*cP_alt == E0 and crP_alt == E0
        i = i+1
    if verbose:
        if ok:
            print("test {}: {}".format(function.__name__, ok))
            return ok
        if ok_alt:
            print("test {}: {}".format(function_alt.__name__, ok_alt))
            return ok_alt
    if ok:
        return ok, "test {}: {}".format(function.__name__, ok)
    if ok_alt:
        return ok_alt, "test {}: {}".format(function_alt.__name__, ok_alt)
    return False, "test {}: {}\n".format(function.__name__, ok) + "test {}: {}".format(function_alt.__name__, ok_alt)

def test_bw6_bls12_g2_mult_by_cofactor(Et, r, c2, t, u, omega, ht, hy,verbose=True):
    t_mod_r = (t % r)
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r = t_mod_r - r
    t_mod_r_mod_u = t_mod_r % u
    if abs(t_mod_r_mod_u-u) < abs(t_mod_r_mod_u):
        t_mod_r_mod_u = t_mod_r_mod_u - u
    if t_mod_r_mod_u == 0:
        function = bw6_bls12_g2_mult_by_cofactor_trace_0_mod_r_u
        function_alt = bw6_bls12_g2_mult_by_cofactor_trace_0_mod_r_u_alt
    else:
        function = bw6_bls12_g2_mult_by_cofactor_trace_3_mod_r_u
        function_alt = bw6_bls12_g2_mult_by_cofactor_trace_3_mod_r_u_alt
    E0 = Et(0)
    ok = True
    ok_alt = True
    i = 0
    while (ok or ok_alt) and i < 10:
        P = Et.random_element()
        rP = r*P
        assert rP != E0 and c2*rP == E0 and c2*P != E0
        if ok:
            c2P = function(P, omega, u, ht, hy)
            c2rP = function(rP, omega, u, ht, hy)
            ok = c2P != E0 and r*c2P == E0 and c2rP == E0
        if ok_alt:
            c2P_alt = function_alt(P, omega, u, ht, hy)
            c2rP_alt = function_alt(rP, omega, u, ht, hy)
            ok_alt = c2P_alt != E0 and r*c2P_alt == E0 and c2rP_alt == E0
        i = i+1
    if verbose:
        if ok:
            print("test {}: {}".format(function.__name__, ok))
            return ok
        if ok_alt:
            print("test {}: {}".format(function_alt.__name__, ok_alt))
            return ok_alt
    if ok:
        return ok, "test {}: {}".format(function.__name__, ok)
    if ok_alt:
        return ok_alt, "test {}: {}".format(function_alt.__name__, ok_alt)
    return False, "test {}: {}\n".format(function.__name__, ok) + "test {}: {}".format(function_alt.__name__, ok_alt)

def test_bw6_bls12_r_subgroup_membersip_testing(E, r, c, t, u, omega):
    t_mod_r = t % r
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r -= r
    t_mod_r_mod_u = t_mod_r % u
    if t_mod_r_mod_u == 0:
        function = bw6_bls12_g1_mult_by_3r_trace_0_mod_r_u
    else:
        function = bw6_bls12_g1_mult_by_3r_trace_3_mod_r_u
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        cP = c*P
        assert cP != E0 and r*cP == E0 and r*P != E0
        R = function(P, omega, u)
        cR = function(cP, omega, u)
        ok = c*R == E0 and cR == E0
        i = i+1
    print("test {}: {}".format(function.__name__, ok))
    return ok

def test_curve(test_vector_value):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    k = 6

    v = test_vector_value

    u0 = ZZ(v['u'])
    u = u0
    ht = v['ht']
    hy = v['hy']
    px_coeffs = v['px']
    px_denom = v['px_denom']
    rx_coeffs = v['rx']
    rx_denom = v['rx_denom']
    cx_coeffs = v['cx']
    cx_denom = v['cx_denom']
    yx_coeffs = v['yx']
    yx_denom = v['yx_denom']
    tx_coeffs = v['tx']
    tx_denom = v['tx_denom']
    g2cx_coeffs = v['g2cx']
    g2cx_denom = v['g2cx_denom']

    betax_coeffs = v['betax']
    betax_denom = v['betax_denom']
    lambx_coeffs = v['lambx']
    lambx_denom = v['lambx_denom']
    
    px = QQx(px_coeffs)/QQ(px_denom)
    qx = px
    rx = QQx(rx_coeffs)/QQ(rx_denom)
    cx = QQx(cx_coeffs)/QQ(cx_denom)
    yx = QQx(yx_coeffs)/QQ(yx_denom)
    tx = QQx(tx_coeffs)/QQ(tx_denom)
    tr0x = tx % rx
    y0x = yx % rx
    g2cx = QQx(g2cx_coeffs)/QQ(g2cx_denom)
    betax = QQx(betax_coeffs)/QQ(betax_denom)
    lambx = QQx(lambx_coeffs)/QQ(lambx_denom)
    ex = (px**2 - px + 1) // rx
    
    p = ZZ(px(u))
    q = p
    r = ZZ(rx(u))
    c = ZZ(cx(u))
    y = ZZ(yx(u))
    t = ZZ(tx(u))
    tr = t
    tr0 = ZZ(tr0x(u))
    y0 = ZZ(y0x(u))
    g2c = ZZ(g2cx(u))
    beta = (betax(u))
    lamb = ZZ(lambx(u))
    ee = ZZ(ex(u))

    b = v['b']
    t_mod_r = t % r
    if abs(t_mod_r - r) < abs(t_mod_r):
        t_mod_r -= r
    t_mod_r_mod_u = (t_mod_r) % abs(u)
    print("t_mod_r_mod_u = {0} = {0:#x}".format(t_mod_r_mod_u))
    tx_mod_rx_mod_x = ((tx % rx) % x)

    print("Curve E: y^2 = x^3 {:+d} /Fp of {} bits, (t mod r) mod u = {}, (tx mod rx) mod x = {}".format(b, p.nbits(), t_mod_r_mod_u, tx_mod_rx_mod_x))
    print("ht= {}".format(ht))
    print("hy= {}".format(hy))
    print("u = {:#x}".format(u))
    print("p = {:#x}".format(p))
    print("r = {:#x}".format(r))
    print("c = {:#x}".format(c))
    print("y = {:#x}".format(y))
    print("t = {:#x}".format(t))

    print("ht^2 + 3*hy^2 = {}".format(ht**2 + 3*hy**2))
    print("tr == tr0 + ht*r: {}".format(t == tr0 + ht*r))
    print("y == y0 + hy*r: {}".format(y == y0 + hy*r))
    print("(tr^2 + 3*y^2) mod 4 = {}".format((t**2 + 3*y**2) % 4))
    print("(tr0^2 + 3*y0^2) % 4 = {}".format((tr0**2 + 3*y0**2) % 4))
    print("tr0 % 3 = {}".format(tr0 % 3))
    if 3*y0 == tr0:
        print("tr0 % 3 = {} y0 == tr0/3".format(tr0 % 3, y0 == tr0/3))
    elif 3*y0 == -tr0:
        print("tr0 % 3 = {} y0 == -tr0/3".format(tr0 % 3, y0 == -tr0/3))
    else:
        print("y0 is none of +/- tr0/3")
    print("c % 3 = {} c2 % 3 = {}".format(c % 3, g2c % 3))
    assert t_mod_r == tr0
    gr = gcd(p+1-tr0, p**2-p+1)
    assert (gr % r) == 0
    cgr = gr//r
    if cgr == 1 or cgr == -1:
        print("gcd(p+1-(tr mod r), p^2-p+1) = r, Scott's trick available")
    else:
        print("gcd(p+1-(tr mod r), p^2-p+1) = {}*r, Scott's trick not available".format(cgr))

    print("c mod 8 = {}\nc2 mod 8 = {}".format(c % 8, g2c % 8))

    # find smallest integer a so that given x mod p, x^a is a permutation
    # will work for gcd(a, p-1) == 1
    # since 2 | (p-1) and 3 | (p-1) by construction, start at 5.
    a = 5
    while gcd(a, p-1) > 1:
        a = a + 1
    print("Smallest a such that x^a mod p is a permutation: a={}".format(a))
    print("Elligator2 available: {} (<=> the curve has a 2-torsion point <=> cofactor is even)".format((c % 2) == 0))
    print("small factors of y that are also factors of (p+1-t) or (p+1+t):")
    print("gcd(y, (p+1-t)*(p+1+t)) = {}".format( (gcd(y, (p+1-t)*(p+1+t)).factor())))
    # try with cubic twists
    print("gcd(y, (p+1-(3*y-t)//2)*(p+1-(-3*y-t)//2)) = {}".format( (gcd(y, (p+1-(3*y-t)//2)*(p+1-(-3*y-t)//2)).factor())))
    # try with sextic twists
    print("gcd(y, (p+1-(3*y+t)//2)*(p+1-(-3*y+t)//2)) = {}".format( (gcd(y, (p+1-(3*y+t)//2)*(p+1-(-3*y+t)//2)).factor())))
    print("--> isogenies of these degrees are defined over Fp.")

    curve_order = q+1-t
    c = (q+1-t) // r
    Fq = GF(q, proof=False)

    if beta in ZZ:
        beta = ZZ(beta)
    else:
        beta = ZZ(Fq(beta))
        if abs(beta-q) < beta:
            beta = beta-q

    Fqz = Fq['z']; (z,) = Fqz._first_ngens(1)
    E = EllipticCurve([Fq(0),Fq(b)])
    b = Fq(b)
    print(E)
    P = E.random_element()
    assert r*c*P == E(0)

    print("\ncost pairing:")
    cost_pairing_bw6_bls12(u0, t_mod_r_mod_u)
    cost_pairing_opt_tate_bw6_bls12(u)
    cost_final_exp_hard_bw6_bls12(u0, ht, hy, trace_0_mod_r_mod_u = t_mod_r_mod_u==0)
    if p.nbits() == 761:
        print("\ncost pairing bw6_761:")
        cost_pairing_bw6_761()

    # beta mod p
    omega = Fq(beta)
    test_bw6_phi(E, r, c, omega)
    test_bw6_phi(E, r, c, -omega-1)
    P = E.random_element()
    phiP = bw6_phi(P, omega)
    phiP_ = bw6_phi(P, -omega-1)
    t_1P = -(t-1)*P
    if phiP == t_1P:
        print("phi with omega has eigenvalue -(t-1)")
    elif phiP_ == t_1P:
        print("adjusting: omega changed to -omega-1 so that Phi has eigenvalue -(t-1)")
        omega = -omega-1
    else:
        # print("error phi(P) has not eigenvalue -(t-1)")
        P = c*P
        phiP = bw6_phi(P, omega)
        phiP_ = bw6_phi(P, -omega-1)
        t_1P = -(t-1)*P
        if phiP == t_1P:
            print("phi with omega has eigenvalue -(t-1) on r-torsion points only")
        elif phiP_ == t_1P:
            print("adjusting: omega changed to -omega-1 so that Phi has eigenvalue -(t-1) on r-torsion points only")
            omega = -omega-1
        else:
            print("error phi(P) either with omega or -omega-1 has not eigenvalue -(t-1) even on r-torsion points")

    xiD, bD = find_twist_curve_parameter_xi_ab(b, Fq, r, g2c, D_twist=True)
    print("xiD = {}".format(xiD))
    ED = EllipticCurve([Fq(0), bD])
    orderD = ED.order()
    assert orderD == g2c*r
    print("ED ok")
    Fq6D = Fq.extension(z**6 - xiD, names=('wD',),proof=False); (wD,) = Fq6D._first_ngens(1)
    E_Fq6D = EllipticCurve([Fq6D(0), Fq6D(b)])
    # this is already an absolute extension
    c2 = orderD // r

    def map_Fq6D_Fp6D(X, aD=None):
        return X

    def map_ED_E_Fq6D(Q2):
        return E_Fq6D(psi_sextic_d_twist(Q2, wD))

    # now find a 6-th M-twist
    xiM, bM = find_twist_curve_parameter_xi_ab(b, Fq, r, g2c, D_twist=False)
    print("xiM = {}".format(xiM))
    EM = EllipticCurve([Fq(0), bM])
    orderM = EM.order()
    assert orderM == r*g2c
    print("EM ok")
    Fq6M = Fq.extension(z**6 - xiM, names=('wM',),proof=False); (wM,) = Fq6M._first_ngens(1)
    E_Fq6M = EllipticCurve([Fq6M(0), Fq6M(b)])
    c2 = orderM // r
    
    def map_Fq6M_Fp6M(X, aM=None):
        return X

    def map_EM_E_Fq6M(Q2):
        return E_Fq6M(psi_sextic_m_twist(Q2, wM))


    test_bw6_bls12_r_subgroup_membersip_testing(E, r, c, t, u, omega)

    test_bw6_bls12_test_cofactor_formula(E, r, c, t, y, u, omega, ht, hy)

    ok_c_omega, str_test = test_bw6_bls12_g1_mult_by_cofactor(E, r, c, t, u, omega, ht, hy, verbose=False)
    if ok_c_omega:
        print("with phi(x,y) = (omega*x, y):")
        print(str_test)
    ok_c_omega_, str_test_ = test_bw6_bls12_g1_mult_by_cofactor(E, r, c, t, u, -omega-1, ht, hy, verbose=False)
    if ok_c_omega_:
        print("with phi(x,y) = ((-omega-1)*x, y):")
        print(str_test_)
    if not ok_c_omega and not ok_c_omega_:
        print("both with phi(x,y) = (omega*x, y) or ((-omega-1)*x, y):")
        print(str_test)
        print(str_test_)

    ok_c2_omega, str_test = test_bw6_bls12_g2_mult_by_cofactor(EM, r, c2, t, u, omega, ht, hy,verbose=False)
    if ok_c2_omega:
        print("with phi(x,y) = (omega*x, y) and M-twist:")
        print(str_test)
    ok_c2_omega_, str_test_ = test_bw6_bls12_g2_mult_by_cofactor(EM, r, c2, t, u, -omega-1, ht, hy,verbose=False)
    if ok_c2_omega_:
        print("with phi(x,y) = ((-omega-1)*x, y) and M-twist:")
        print(str_test_)
    if not ok_c2_omega and not ok_c2_omega_:
        print("both with phi(x,y) = (omega*x, y) or ((-omega-1)*x, y) and M-twist:")
        print(str_test)
        print(str_test_)

    ok_c2_omega, str_test = test_bw6_bls12_g2_mult_by_cofactor(ED, r, c2, t, u, omega, ht, hy,verbose=False)
    if ok_c2_omega:
        print("with phi(x,y) = (omega*x, y) and D-twist:")
        print(str_test)
    ok_c2_omega_, str_test_ = test_bw6_bls12_g2_mult_by_cofactor(ED, r, c2, t, u, -omega-1, ht, hy,verbose=False)
    if ok_c2_omega_:
        print("with phi(x,y) = ((-omega-1)*x, y) and D-twist:")
        print(str_test_)
    if not ok_c2_omega and not ok_c2_omega_:
        print("both with phi(x,y) = (omega*x, y) or ((-omega-1)*x, y) and D-twist:")
        print(str_test)
        print(str_test_)

    print("\ntest pairings with M-twist")

    test_g2_frobenius_eigenvalue_bw6(E_Fq6M, EM, r, c2, D_twist=False)

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

    if t_mod_r_mod_u == 3:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761_2naf,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_miller_loop_opt_ate_bw6_761_all(E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_alt,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_2naf,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)

    test_double_line_h_a0_twist6_aklgl(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_add_line_h_a0_twist6_aklgl_test(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl_with_z(E, EM, wM, D_twist=False)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_m6_twist(Fq6M)

    test_miller_function_ate_aklgl(E,EM,Fq6M,xiM,r,c,c2,t-1,D_twist=False)

    if t_mod_r_mod_u == 3:
        test_formula_miller_loop_opt_ate_bw6_761_aklgl(E,EM,Fq6M,r,c,c2,u0,D_twist=False)

        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_miller_loop_opt_ate_bw6_761_aklgl_all(E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_miller_loop_opt_ate_all_bw6_761_and_aklgl(E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)

    omega_ = -omega-1
    k_bls = 12

    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fq6M, u, p, r, tr, c, c2, omega_, map_Fq6M_Fp6M, k_bls, D_twist=False, function_name=miller_function_tate_aklgl_a0_2_parallel)
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fq6M, u, p, r, tr, c, c2, omega, map_Fq6M_Fp6M, k_bls, D_twist=False, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar)

    print("\ntest pairings with D-twist")

    test_g2_frobenius_eigenvalue_bw6(E_Fq6D, ED, r, c2, D_twist=True)

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

    if t_mod_r_mod_u == 3:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_761_2naf,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_miller_loop_opt_ate_bw6_761_all(E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_2naf,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)

    test_double_line_h_a0_twist6_aklgl(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_add_line_h_a0_twist6_aklgl_test(E, ED, wD, D_twist=True)
    test_add_line_h_a0_twist6_aklgl_with_z(E,ED,wD,D_twist=True)
    test_sparse_sparse_mult_d6_twist(Fq6D)
    test_sparse_mult_d6_twist(Fq6D)

    test_miller_function_ate_aklgl(E,ED,Fq6D,xiD,r,c,c2,t-1,D_twist=True)

    if t_mod_r_mod_u == 3:
        test_formula_miller_loop_opt_ate_bw6_761_aklgl(E,ED,Fq6D,r,c,c2,u0,D_twist=True)

        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_miller_loop_opt_ate_bw6_761_aklgl_all(E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_miller_loop_opt_ate_all_bw6_761_and_aklgl(E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)

    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fq6D, u, p, r, tr, c, c2, omega_, map_Fq6D_Fp6D, k_bls, D_twist=True, function_name=miller_function_tate_aklgl_a0_2_parallel)
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fq6D, u, p, r, tr, c, c2, omega, map_Fq6D_Fp6D, k_bls, D_twist=True, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar)

    print("test final exp hard part:")
    if t_mod_r_mod_u == 0:
        ee0 = ee*(u+1)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u, Fq6M, u0, ht, hy, r, t, expected_exp=ee0)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u, Fq6D, u0, ht, hy, r, t, expected_exp=ee0)
        ee0_alt = ee*(u**3-u**2-u)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_alt, Fq6M, u0, ht, hy, r, t, expected_exp=ee0_alt)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_alt, Fq6D, u0, ht, hy, r, t, expected_exp=ee0_alt)
    else:
        ee3 = ee*(u+1)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u, Fq6M, u0, ht, hy, r, t, expected_exp=ee3)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u, Fq6D, u0, ht, hy, r, t, expected_exp=ee3)
        ee3_alt = ee*(u**3-u**2+1)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_alt, Fq6M, u0, ht, hy, r, t, expected_exp=ee3_alt)
        test_final_exp_hard_bw6(final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_alt, Fq6D, u0, ht, hy, r, t, expected_exp=ee3_alt)

if __name__ == "__main__":
    arithmetic(False)
    print("BLS12-379")
    for i in [8,4,0,1]:
        v = testvector_bls12_379_bw6_768[i]
        #if v['pnbits'] == 764 and v['ht'] == -25 and v['hy'] == 3:
        test_curve(v)
    print("##################################################")
    print("BLS12-377")
    for i in [0,1,2,3]:
        v = testvector_bls12_377_bw6_768[i]
        test_curve(v)

    print("##################################################")
    print("BLS12-381")
    for i in [0]:
        v = testvector_bls12_381_bw6_768[i]
        test_curve(v)
