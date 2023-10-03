from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from pairing import *
from pairing_bw6_bn import *
from pairing_cocks_pinch import miller_function_tate_aklgl_a0_2_parallel, miller_function_tate_a0_twist6_aklgl_2_multi_scalar

from cost_pairing import cost_pairing_bw6_bn
from test_pairing import *
from test_pairing_bw6 import *
from test_pairing_bw6_761 import test_g2_frobenius_eigenvalue_bw6
from test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761
from test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761_aklgl
from test_pairing_cp12_bls12 import test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar

from testvector_bn_253_bw6_512 import testvector_bn_253_bw6_512
from testvector_bn_253b_bw6_512 import testvector_bn_253_bw6_512 as testvector_bn_253b_bw6_512
from testvector_bn_254_bw6_512 import testvector_bn_254_bw6_512
from testvector_bn_254e_bw6_516 import testvector_bn_254_bw6_516
from testvector_bn_382_bw6_767 import testvector_bn_382_bw6_767
from testvector_bn_446_bw6_896 import testvector_bn_446_bw6_896

# test subgroup membership testing and cofactor multiplication
prod_primes = prod(prime_range(10**7))

def test_cofactor_multiplication_formulas_g1(E, r, t, u, ht, hy, omega):
    """
    trace is 0 mod r mod u, eigenvalue corresponds to -(t0-1) mod r
    fast subgroup multiplication c*P = c0*P + c1*phi(P)
    c0 = (ht^2+3*hy^2)/4*(2*u)         + (ht+3*hy)/2*u       + hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) - (ht-3*hy)/2*(u+1)   - hy
    k0 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) - (ht+hy)  /2*(2*u+1) + hy*u
    k1 = (ht^2+3*hy^2)/4*(-2*u)        - (ht+3*hy)/2*u       - hy
    with the other eigenvalue that corresponds to (t0-2) mod r
    c0 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) - (ht-3*hy)/2*(u+1)   - hy
    c1 = (ht^2+3*hy^2)/4*(2*u)         + (ht+3*hy)/2*u       + hy
    k0 = (ht^2+3*hy^2)/4*(-2*u)        - (ht+3*hy)/2*u       - hy
    k1 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) - (ht+hy)  /2*(2*u+1) + hy*u
    trace is 3 mod r mod u, eigenvalue corresponds to -(t3-1) mod r
    c0 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) + (ht+3*hy)/2*u       + hy
    c1 = (ht^2+3*hy^2)/4*(-2*u)        - (ht-3*hy)/2*(u+1)   - hy
    k0 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) + (ht+hy)  /2*(2*u+1) - hy*u
    k1 = (ht^2+3*hy^2)/4*(2*u)         + (ht-3*hy)/2*(u+1)   + hy
    with the other eigenvalue that corresponds to (t3-2) mod r
    c0 = (ht^2+3*hy^2)/4*(-2*u)        - (ht-3*hy)/2* (u+1)  - hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) + (ht+3*hy)/2* u      + hy
    k0 = (ht^2+3*hy^2)/4*(2*u)         + (ht-3*hy)/2*(u+1)   + hy
    k1 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) + (ht+hy)  /2*(2*u+1) - hy*u
    """
    assert ((ht**2+3*hy**2) % 4) == 0
    assert ((ht+3*hy) % 2) == 0
    assert ((ht-3*hy) % 2) == 0

    t_mod_r = (t % r)
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r = t_mod_r - r
    t_mod_r_mod_u = t_mod_r % u
    if abs(t_mod_r_mod_u-u) < abs(t_mod_r_mod_u):
        t_mod_r_mod_u = t_mod_r_mod_u - u

    if t_mod_r_mod_u == 0:
        l1 = (ht**2+3*hy**2)//4*(2*u)          + (ht+3*hy)//2*u       + hy
        l2 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) - (ht-3*hy)//2*(u+1)   - hy
        k1 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) - (ht+hy)  //2*(2*u+1) + hy*u
        k2 = (ht**2+3*hy**2)//4*(-2*u)         - (ht+3*hy)//2*u       - hy
        m1 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) - (ht-3*hy)//2*(u+1)   - hy
        m2 = (ht**2+3*hy**2)//4*(2*u)          + (ht+3*hy)//2*u       + hy
        n1 = (ht**2+3*hy**2)//4*(-2*u)         - (ht+3*hy)//2*u       - hy
        n2 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) - (ht+hy)  //2*(2*u+1) + hy*u
    if t_mod_r_mod_u == 3:
        l1 = (ht**2+3*hy**2)//4*(-2*u)         - (ht-3*hy)//2 *(u+1)  - hy
        l2 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) + (ht+3*hy)//2 *u      + hy
        k1 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) + (ht+hy)  //2*(2*u+1) - hy*u
        k2 = (ht**2+3*hy**2)//4*(2*u)          + (ht-3*hy)//2*(u+1)   + hy
        m1 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) + (ht+3*hy)//2* u      + hy
        m2 = (ht**2+3*hy**2)//4*(-2*u)         - (ht-3*hy)//2* (u+1)  - hy
        n1 = (ht**2+3*hy**2)//4*(2*u)          + (ht-3*hy)//2*(u+1)   + hy
        n2 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) + (ht+hy)  //2*(2*u+1) - hy*u

    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        rP = r*P
        omeg_rP = bw6_phi(rP, omega)
        omeg_rP_ = bw6_phi(rP, -omega-1)
        omeg_P = bw6_phi(P, omega)
        omeg_P_ = bw6_phi(P, -omega-1)
        if t_mod_r_mod_u == 0:
            ok0 = r*(l1*P + l2*omeg_P) == E0
            oka = r*(k1*P + k2*omeg_P) == E0
            # other -omega-1
            ok1 = r*(m1*P + m2*omeg_P_) == E0
            okb = r*(n1*P + n2*omeg_P_) == E0
            ok = ok0 and oka and ok1 and okb
        if t_mod_r_mod_u == 3:
            ok0 = r*(l1*P + l2*omeg_P) == E0
            oka = r*(k1*P + k2*omeg_P) == E0
            # other eigenvalue
            ok1 = r*(m1*P + m2*omeg_P_) == E0
            okb = r*(n1*P + n2*omeg_P_) == E0
            ok = ok0 and oka and ok1 and okb
        i = i+1
    print("test_cofactor_multiplication_formulas_g1: {} ({} {} {} {})".format(ok, ok0, oka, ok1, okb))

def test_cofactor_multiplication_formulas_g2(E2, r, t2, u, ht, hy, omega):
    """
    trace is 0 mod r mod u, eigenvalue corresponds to -(t0-1) mod r
    fast subgroup multiplication c*P = c0*P + c1*phi(P)
    with the other eigenvalue that corresponds to (t0-2) mod r
    trace is 3 mod r mod u, eigenvalue corresponds to -(t3-1) mod r
    with the other eigenvalue that corresponds to (t3-2) mod r
    """
    assert ((ht**2+3*hy**2) % 4) == 0
    assert ((ht+3*hy) % 2) == 0
    assert ((ht-3*hy) % 2) == 0
    assert ((ht+  hy) % 2) == 0
    assert ((ht-  hy) % 2) == 0

    t_mod_r = (t2 % r)
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r = t_mod_r - r
    t_mod_r_mod_u = t_mod_r % u
    if abs(t_mod_r_mod_u-u) < abs(t_mod_r_mod_u):
        t_mod_r_mod_u = t_mod_r_mod_u - u

    if t_mod_r_mod_u == 0:
        l1 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) - (ht-3*hy)//2*u     + hy
        l2 = (ht**2+3*hy**2)//4*(-2*u)         + (ht+3*hy)//2*(u+1) - hy
        k1 = (ht**2+3*hy**2)//4*(2*u)          - (ht+3*hy)//2*(u+1) + hy
        k2 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) - (ht-  hy)//2       - ht*u
        m1 = (ht**2+3*hy**2)//4*(-2*u)         + (ht+3*hy)//2*(u+1) - hy
        m2 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) - (ht-3*hy)//2*u     + hy
        n1 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) - (ht-  hy)//2       - ht*u
        n2 = (ht**2+3*hy**2)//4*(2*u)          - (ht+3*hy)//2*(u+1) + hy
    if t_mod_r_mod_u == 3:
        l1 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) + (ht+3*hy)//2*(u+1) - hy
        l2 = (ht**2+3*hy**2)//4*(2*u)          - (ht-3*hy)//2*u     + hy
        k1 = (ht**2+3*hy**2)//4*(-2*u)         + (ht-3*hy)//2*u     - hy
        k2 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) + (ht-hy)//2         + ht*u
        m1 = (ht**2+3*hy**2)//4*(2*u)          - (ht-3*hy)//2*u     + hy
        m2 = (ht**2+3*hy**2)//4*(6*u**2+4*u+1) + (ht+3*hy)//2*(u+1) - hy
        n1 = (ht**2+3*hy**2)//4*(6*u**2+2*u+1) + (ht-  hy)//2       + ht*u
        n2 = (ht**2+3*hy**2)//4*(-2*u)         + (ht-3*hy)//2*u     - hy

    E0 = E2(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E2.random_element()
        rP = r*P
        omeg_rP = bw6_phi(rP, omega)
        omeg_rP_ = bw6_phi(rP, -omega-1)
        omeg_P = bw6_phi(P, omega)
        omeg_P_ = bw6_phi(P, -omega-1)
        if t_mod_r_mod_u == 0:
            ok0 = r*(l1*P + l2*omeg_P) == E0
            oka = r*(k1*P + k2*omeg_P) == E0
            # other -omega-1
            ok1 = r*(m1*P + m2*omeg_P_) == E0
            okb = r*(n1*P + n2*omeg_P_) == E0
            ok = ok0 and oka and ok1 and okb
        if t_mod_r_mod_u == 3:
            ok0 = r*(l1*P + l2*omeg_P) == E0
            oka = r*(k1*P + k2*omeg_P) == E0
            # other eigenvalue
            ok1 = r*(m1*P + m2*omeg_P_) == E0
            okb = r*(n1*P + n2*omeg_P_) == E0
            ok = ok0 and oka and ok1 and okb
        i = i+1
    print("test_cofactor_multiplication_formulas_g2: {} ({} {} {} {})".format(ok, ok0, oka, ok1, okb))

def test_curve(v):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    k = 6
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
    g2c = ZZ(g2cx(u))
    c2 = g2c
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

    if t_mod_r_mod_u != 0 and t_mod_r_mod_u != 3:
        raise ValueError("error (t % r) % u = {} is neither 0 nor 3".format(t_mod_r_mod_u))

    print("Curve E: y^2 = x^3 {:+d} /Fp of {} bits, (t mod r) mod u = {}, (tx mod rx) mod x = {}".format(b, p.nbits(), t_mod_r_mod_u, tx_mod_rx_mod_x))
    print("u = {:#x} {} bits".format(u, u.nbits()))
    print("ht= {:#x} hy= {:#x}".format(ht, hy))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("c = {:#x} {} bits".format(c, c.nbits()))
    print("y = {:#x} {} bits".format(y, y.nbits()))
    print("t = {:#x} {} bits".format(t, t.nbits()))

    cost_pairing_bw6_bn(u, ht, hy, tr_mod=t_mod_r_mod_u)

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
        print("phi with -omega-1 has eigenvalue -(t-1)")
        omega = -omega-1
    else:
        print("error phi(P) has not eigenvalue -(t-1)")
        P = c*P
        phiP = bw6_phi(P, omega)
        phiP_ = bw6_phi(P, -omega-1)
        t_1P = -(t-1)*P
        if phiP == t_1P:
            print("phi with omega has eigenvalue -(t-1) on r-torsion points only")
        elif phiP_ == t_1P:
            print("phi with -omega-1 has eigenvalue -(t-1) on r-torsion points only")
            omega = -omega-1
        else:
            print("error phi(P) has not eigenvalue -(t-1) even on r-torsion points")
    xiD, bD = find_twist_curve_parameter_xi_ab(b, Fq, r, c2, D_twist=True)
    print("xiD = {}".format(xiD))
    ED = EllipticCurve([Fq(0), bD])
    orderD = ED.order()
    assert orderD == c2*r
    print("ED ok")
    Fq6D = Fq.extension(z**6 - xiD, names=('wD',),proof=False); (wD,) = Fq6D._first_ngens(1)
    E_Fq6D = EllipticCurve([Fq6D(0), Fq6D(b)])
    # this is already an absolute extension

    def map_Fq6D_Fp6D(X, aD=None):
        return X

    def map_ED_E_Fq6D(Q2):
        return E_Fq6D(psi_sextic_d_twist(Q2, wD))

    # now find a 6-th M-twist
    xiM, bM = find_twist_curve_parameter_xi_ab(b, Fq, r, c2, D_twist=False)
    print("xiM = {}".format(xiM))
    EM = EllipticCurve([Fq(0), bM])
    orderM = EM.order()
    assert orderM == r*c2
    print("EM ok")
    Fq6M = Fq.extension(z**6 - xiM, names=('wM',),proof=False); (wM,) = Fq6M._first_ngens(1)
    E_Fq6M = EllipticCurve([Fq6M(0), Fq6M(b)])

    def map_Fq6M_Fp6M(X, aM=None):
        return X

    def map_EM_E_Fq6M(Q2):
        return E_Fq6M(psi_sextic_m_twist(Q2, wM))

    print("c2 = {} # {} bits".format(c2, c2.nbits()))
    print("g2cx = {}\n     = {}".format(g2cx, g2cx.factor()))
    m = c2
    c2_0 = ZZ(1)
    g = gcd(c2, prod_primes)
    while g > 1:
        m = m // g
        c2_0 = c2_0 * g
        g = gcd(m, g)
    print("c2 = ({}) * {}".format(c2_0.factor(), m))

    test_cofactor_multiplication_formulas_g1(E, r, t, u, ht, hy, omega)
    if t_mod_r_mod_u == 3:
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u, E, r, c, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u_, E, r, c, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u_alt, E, r, c, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u_alt_, E, r, c, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_r(bw6_bn_g1_mult_by_r_trace_3_mod_r_u, E, r, c, u, omega)
        test_bw6_g1_mult_by_r(bw6_bn_g1_mult_by_r_trace_3_mod_r_u_alt, E, r, c, u, -omega-1)
        test_bw6_g1_r_subgroup_membersip_testing(bw6_bn_g1_check_membership_trace_3_mod_r_u, E, r, c, u, omega)
        test_bw6_g1_r_subgroup_membersip_testing(bw6_bn_g1_check_membership_trace_3_mod_r_u_alt, E, r, c, u, -omega-1)
    if t_mod_r_mod_u == 0:
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u, E, r, c, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u_, E, r, c, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u_alt, E, r, c, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u_alt_, E, r, c, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_r(bw6_bn_g1_mult_by_r_trace_0_mod_r_u, E, r, c, u, omega)
        test_bw6_g1_mult_by_r(bw6_bn_g1_mult_by_r_trace_0_mod_r_u_alt, E, r, c, u, -omega-1)
        test_bw6_g1_r_subgroup_membersip_testing(bw6_bn_g1_check_membership_trace_0_mod_r_u, E, r, c, u, omega)
        test_bw6_g1_r_subgroup_membersip_testing(bw6_bn_g1_check_membership_trace_0_mod_r_u_alt, E, r, c, u, -omega-1)

    print("test G2 cofactor multiplication, M-twist")
    test_cofactor_multiplication_formulas_g2(EM, r, t, u, ht, hy, -omega-1)
    if t_mod_r_mod_u == 3:
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u, EM, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_, EM, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_alt, EM, r, c2, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_alt_, EM, r, c2, u, omega, ht, hy)
    if t_mod_r_mod_u == 0:
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u, EM, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_, EM, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_alt, EM, r, c2, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_alt_, EM, r, c2, u, omega, ht, hy)

    print("test G2 cofactor multiplication, D-twist")
    test_cofactor_multiplication_formulas_g2(ED, r, t, u, ht, hy, -omega-1)
    if t_mod_r_mod_u == 3:
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u, ED, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_, ED, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_alt, ED, r, c2, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_alt_, ED, r, c2, u, omega, ht, hy)
    if t_mod_r_mod_u == 0:
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u, ED, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_, ED, r, c2, u, -omega-1, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_alt, ED, r, c2, u, omega, ht, hy)
        test_bw6_g1_mult_by_cofactor(bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_alt_, ED, r, c2, u, omega, ht, hy)

    print("\ntest pairings with M-twist")

    if t_mod_r_mod_u == 3:
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_3_mod_r_mod_u, Fq6M, u, ht, hy, r, t, expected_exp=ee*(2*u))
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_3_mod_r_mod_u_alt, Fq6M, u, ht, hy, r, t, expected_exp=ee*(6*u**2+2*u+1))
    if t_mod_r_mod_u == 0:
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_0_mod_r_mod_u, Fq6M, u, ht, hy, r, t, expected_exp=ee*(2*u))
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_0_mod_r_mod_u_alt, Fq6M, u, ht, hy, r, t, expected_exp=ee*(6*u**2+4*u+1))

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
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_2naf,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
    if t_mod_r_mod_u == 0:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_2naf,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)

    test_double_line_h_a0_twist6_aklgl(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_add_line_h_a0_twist6_aklgl_test(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl_with_z(E, EM, wM, D_twist=False)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_m6_twist(Fq6M)

    test_miller_function_ate_aklgl(E,EM,Fq6M,xiM,r,c,c2,t-1,D_twist=False)

    if t_mod_r_mod_u == 3:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
    if t_mod_r_mod_u == 0:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)

    print("\ntest pairings with D-twist")

    if t_mod_r_mod_u == 3:
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_3_mod_r_mod_u, Fq6D, u, ht, hy, r, t, expected_exp=ee*(2*u))
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_3_mod_r_mod_u_alt, Fq6D, u, ht, hy, r, t, expected_exp=ee*(6*u**2+2*u+1))
    if t_mod_r_mod_u == 0:
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_0_mod_r_mod_u, Fq6D, u, ht, hy, r, t, expected_exp=ee*(2*u))
        test_final_exp_hard_bw6(final_exp_hard_bw6_bn_trace_0_mod_r_mod_u_alt, Fq6D, u, ht, hy, r, t, expected_exp=ee*(6*u**2+4*u+1))

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
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_2naf,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
    if t_mod_r_mod_u == 0:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_2naf,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)

    test_double_line_h_a0_twist6_aklgl(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_add_line_h_a0_twist6_aklgl_test(E, ED, wD, D_twist=True)
    test_add_line_h_a0_twist6_aklgl(E, ED, wD, D_twist=True)
    test_add_line_h_a0_twist6_aklgl_with_z(E, ED, wD, D_twist=True)
    test_sparse_sparse_mult_d6_twist(Fq6D)
    test_sparse_mult_d6_twist(Fq6D)

    test_miller_function_ate_aklgl(E,ED,Fq6D,xiD,r,c,c2,t-1,D_twist=True)

    if t_mod_r_mod_u == 3:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
    if t_mod_r_mod_u == 0:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)


if __name__ == "__main__":
    arithmetic(False)
    print("BN-253 BW6-509")
    for v in testvector_bn_253_bw6_512[:2]:
        test_curve(v)
    print("BN-253b BW6-509")
    for v in testvector_bn_253b_bw6_512[:2]:
        test_curve(v)
    print("BN-254 BW6-512 (Geppetto)")
    for v in testvector_bn_254_bw6_512[:2]:
        test_curve(v)
    print("BN-254 Ethereum BW6-516")
    for v in testvector_bn_254_bw6_516[:2]:
        test_curve(v)
    print("BN-382 BW6-767")
    for v in testvector_bn_382_bw6_767[:2]:
        test_curve(v)
    print("BN-446 Pluto BW6-896")
    for v in testvector_bn_446_bw6_896:
        test_curve(v)
