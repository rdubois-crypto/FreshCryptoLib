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
from pairing_bw6_bls24 import *
from pairing_cocks_pinch import miller_function_tate_aklgl_a0_2_parallel, miller_function_tate_a0_twist6_aklgl_2_multi_scalar

from cost_pairing import cost_final_exp_hard_bw6_bls24, cost_pairing_bw6_bls24, cost_pairing_opt_tate_bw6_bls24
from test_pairing import *
from test_pairing_bw6_761 import test_g2_frobenius_eigenvalue_bw6
from test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761
from test_pairing_bw6_761 import test_bilinear_miller_loop_opt_ate_bw6_761_aklgl
from test_pairing_cp12_bls12 import test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar

from testvector_bls24_315_bw6_640 import testvector_bls24_315_bw6_640
from testvector_bls24_315_bw6_672 import testvector_bls24_315_bw6_672

# test subgroup membership testing and cofactor multiplication
prod_primes = prod(prime_range(10**7))

def test_bw6_bls24_g1_mult_by_cofactor(E, r, c, t, u, omega, ht, hy,verbose=True):
    t_mod_r = (t % r)
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r = t_mod_r - r
    t_mod_r_mod_u = t_mod_r % u
    if abs(t_mod_r_mod_u-u) < abs(t_mod_r_mod_u):
        t_mod_r_mod_u = t_mod_r_mod_u - u
    if t_mod_r_mod_u == 0:
        function = bw6_bls24_g1_mult_by_cofactor_trace_0_mod_r_u
        function_alt = bw6_bls24_g1_mult_by_cofactor_trace_0_mod_r_u_alt
    else:
        function = bw6_bls24_g1_mult_by_cofactor_trace_3_mod_r_u
        function_alt = bw6_bls24_g1_mult_by_cofactor_trace_3_mod_r_u_alt
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

def test_bw6_bls24_g2_mult_by_cofactor(Et, r, c2, t, u, omega, ht, hy,verbose=True):
    t_mod_r = (t % r)
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r = t_mod_r - r
    t_mod_r_mod_u = t_mod_r % u
    if abs(t_mod_r_mod_u-u) < abs(t_mod_r_mod_u):
        t_mod_r_mod_u = t_mod_r_mod_u - u
    if t_mod_r_mod_u == 0:
        function = bw6_bls24_g2_mult_by_cofactor_trace_0_mod_r_u
        function_alt = bw6_bls24_g2_mult_by_cofactor_trace_0_mod_r_u_alt
    else:
        function = bw6_bls24_g2_mult_by_cofactor_trace_3_mod_r_u
        function_alt = bw6_bls24_g2_mult_by_cofactor_trace_3_mod_r_u_alt
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

def test_bw6_bls24_r_subgroup_membersip_testing(E, r, c, t, u, omega):
    t_mod_r = t % r
    if abs(t_mod_r-r) < abs(t_mod_r):
        t_mod_r -= r
    t_mod_r_mod_u = t_mod_r % u
    if t_mod_r_mod_u == 0:
        function = bw6_bls24_g1_mult_by_3r_trace_0_mod_r_u
    else:
        function = bw6_bls24_g1_mult_by_3r_trace_3_mod_r_u
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

def test_optimal_ate_formula_bw6_bls24(E_Fpk, E2, w, u, r, t, c2, D_twist=False):
    """
    for Q in the trace-0 subgroup of E(Fpk) (G2):
    if ((t % r) % u) == 0:
    Test (u+1)*(-Q) + (u^5-u^4+1)*pi(Q) = 0
    Test (u^5-u^4-u)*Q + (u+1)*pi(Q) = 0
    
    if ((t % r) % u) == 3:
    Test (u+1)*Q + (u^5-u^4-u)*pi(Q) = 0 
    Test (u^5-u^4+1)*Q + (u+1)*pi(-Q) = 0

    INPUT:
    - `E_Fpk`: EllipticCurve instance defined over an absolute extension of Fp
    - `E2`: EllipticCurve instance defined over an absolute extension of Fp of degree p^{k/d}
    - `map_Fqd_Fpk`: map from the relative extension Fqd to the isomorphic absolute extension Fpk
    - `w`: the generator of Fpd (for the twisting map)
    - `u`: integer, parameter seed
    - `r`: prime integer, E(Fp) has order r*c
    - `t`: integer, trace of E(Fp)
    - `c2`: integer, twist cofactor, E2 has order r*c2
    - `D_twist`: whether the twist is a D-twist or M-twist
    """
    p = E_Fpk.base_field().characteristic()
    E2_0 = E2(0)
    E_Fpk_0 = E_Fpk(0)
    t_mod_r = t % r
    if abs(t_mod_r - r) < abs(t_mod_r):
        t_mod_r -= r
    t_mod_r_mod_u = (t_mod_r) % abs(u)
    if t_mod_r_mod_u == 0:
        tr_u = True
        tr_3u = False
    elif t_mod_r_mod_u == 3:
        tr_u = False
        tr_3u = True
    else:
        raise ValueError("Error the trace of the curve bw6_bls24 should be ((t mod r) mod u) = 0 or 3 but this is {}, and u = {}".format(t_mod_r_mod_u, u))
    ok1 = True
    ok2 = True
    i = 0
    while ok1 and ok2 and i < 10:
        Q = c2*E2.random_element()
        while Q == E2_0 or r*Q != E2_0:
            Q = c2*E2.random_element()
        if D_twist:
            Q2 = psi_sextic_d_twist(Q, w)
        else:
            Q2 = psi_sextic_m_twist(Q, w)
        Q = E_Fpk((Q2[0], Q2[1]))
        piQ = E_Fpk(((Q[0]).frobenius(), (Q[1]).frobenius()))
        if tr_u:
            ok1 = (u+1)*(-Q) + (u**5-u**4+1)*piQ == E_Fpk_0
            ok2 = u*(u**4-u**3-1)*Q + (u+1)*piQ == E_Fpk_0
        else:
            ok1 = (u+1)*Q + u*(u**4-u**3-1)*piQ == E_Fpk_0
            ok2 = (u**5-u**4+1)*Q + (u+1)*(-piQ) == E_Fpk_0
        i = i+1
    if D_twist:
        str_tw = "D"
    else:
        str_tw = "M"
    if tr_u:
        print("test optimal ate formula {}_twist (u+1)*(-Q) + (u^5-u^4+1)*pi(Q) = 0: {}".format(str_tw, ok1))
        print("test optimal ate formula {}_twist u*(u^4-u^3-1)*Q + (u+1)*pi(Q)  = 0: {}".format(str_tw, ok2))
    else:
        print("test optimal ate formula {}_twist (u+1)*Q + u*(u^4-u^3-1)*pi(Q)  = 0: {}".format(str_tw, ok1))
        print("test optimal ate formula {}_twist (u^5-u^4+1)*Q + (u+1)*pi(-Q)   = 0: {}".format(str_tw, ok2))
    return ok1 and ok2

def miller_loop_bw6_bls24_naive_trace_u(Q,P,u,verbose=False):
    """
    non-optimized optimal ate Miller loop for BW6 (coming from BLS24)
    where the trace mod r is a multiple of u:
    tr mod r = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x

    returns f_{-u-1,Q}(P)*Frobenius(f_{u^5-u^4+1,Q}(P))
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    # negative Q
    uu = -(u+1)
    if uu < 0:
        Q_ = -Q
        uu = -uu
    else:
        Q_ = Q
    m_u1,Su1 = miller_function_ate(Q_, P, 0, uu)
    if verbose:
        Z1 = 1/Su1[2]
        Z2 = Z1**2
        S1 = (Su1[0]*Z2, Su1[1]*Z1*Z2,1,1)
        print("f_u1 = {}\nQu1 = {}\nQu1 = {}".format(m_u1, Su1, S1))
    v = u**5-u**4+1
    if v < 0:
        Q_ = -Q
        v = -v
    else:
        Q_ = Q
    m_v,Sv = miller_function_ate(Q_, P, 0, v)
    #return m_u1 * m_v**p;
    return m_u1 * m_v.frobenius()

def miller_loop_bw6_bls24_naive_trace_u_alt(Q,P,u,verbose=False):
    """
    non-optimized optimal ate Miller loop for BW6 (coming from BLS24)
    where the trace mod r is a multiple of u:
    tr mod r = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x

    returns f_{u^5-u^4-u,Q}(P)*Frobenius(f_{u+1,Q}(P))
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    # negative Q
    uu = u**5-u**4-u
    if uu < 0:
        Q_ = -Q
        uu = -uu
    else:
        Q_ = Q
    m_u1,Su1 = miller_function_ate(Q_, P, 0, uu)
    if verbose:
        Z1 = 1/Su1[2]
        Z2 = Z1**2
        S1 = (Su1[0]*Z2, Su1[1]*Z1*Z2,1,1)
        print("f_u1 = {}\nQu1 = {}\nQu1 = {}".format(m_u1, Su1, S1))
    v = u+1
    if v < 0:
        Q_ = -Q
        v = -v
    else:
        Q_ = Q
    m_v,Sv = miller_function_ate(Q_, P, 0, v)
    #return m_u1 * m_v**p;
    return m_u1 * m_v.frobenius()

def miller_loop_bw6_bls24_naive_trace_3u(Q,P,u,verbose=False):
    """
    non-optimized optimal ate Miller loop for BW6 (coming from BLS24)
    where the trace mod r is 3 mod u:
    tr mod r = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3
    returns f_{u+1,Q}(P)*Frobenius(f_{u^5-u^4-u,Q}(P))
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    # negative Q
    uu = u+1
    if uu < 0:
        Q_ = -Q
        uu = -uu
    else:
        Q_ = Q
    m_u1,Su1 = miller_function_ate(Q_, P, 0, uu)
    if verbose:
        Z1 = 1/Su1[2]
        Z2 = Z1**2
        S1 = (Su1[0]*Z2, Su1[1]*Z1*Z2,1,1)
        print("f_u1 = {}\nQu1 = {}\nQu1 = {}".format(m_u1, Su1, S1))
    v = u**5-u**4-u
    if v < 0:
        Q_ = -Q
        v = -v
    else:
        Q_ = Q
    m_v,Sv = miller_function_ate(Q_, P, 0, v)
    #return m_u1 * m_v**p;
    return m_u1 * m_v.frobenius()

def miller_loop_bw6_bls24_naive_trace_3u_alt(Q,P,u,verbose=False):
    """
    non-optimized optimal ate Miller loop for BW6 (coming from BLS24)
    where the trace mod r is 3 mod u:
    tr mod r = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3

    returns f_{u^5-u^4+1,Q}(P)*Frobenius(f_{-u-1,Q}(P))
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    # negative Q
    uu = u**5-u**4+1
    if uu < 0:
        Q_ = -Q
        uu = -uu
    else:
        Q_ = Q
    m_u1,Su1 = miller_function_ate(Q_, P, 0, uu)
    if verbose:
        Z1 = 1/Su1[2]
        Z2 = Z1**2
        S1 = (Su1[0]*Z2, Su1[1]*Z1*Z2,1,1)
        print("f_u1 = {}\nQu1 = {}\nQu1 = {}".format(m_u1, Su1, S1))
    v = -(u+1)
    if v < 0:
        Q_ = -Q
        v = -v
    else:
        Q_ = Q
    m_v,Sv = miller_function_ate(Q_, P, 0, v)
    #return m_u1 * m_v**p;
    return m_u1 * m_v.frobenius()

def test_final_exp_hard_bw6_bls24(Fpk, u, ht, hy, r, function):
    """
    Test final_exp_hard_bw6_bls24_trace_{0,3}_mod_r_mod_u{,_alt} functions
    function is one of
    final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u(m, u0, ht, hy)
    final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt(m, u0, ht, hy)
    final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u(m, u0, ht, hy)
    final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt(m, u0, ht, hy)
    """
    i = 0
    ok = True
    p = Fpk.characteristic()
    if function is final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u:
        expected_exp = (u+1)*(p**2-p+1)//r
    elif function is final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt:
        expected_exp = (u**5-u**4-u)*(p**2-p+1)//r
    elif function is final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u:
        expected_exp = (u+1)*(p**2-p+1)//r
    elif function is final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt:
        expected_exp = (u**5-u**4+1)*(p**2-p+1)//r
    else:
        expected_exp = 3*(u+1)*(p**2-p+1)//r
    while ok and i < 10:
        m = Fpk.random_element()
        m = final_exp_easy_k6(m)
        s = function(m, u, ht, hy)
        f = m**expected_exp
        ok1 = f == s
        ok2 = s**r == 1
        ok = ok1 and ok2
        i = i+1
    print("test {}: {}".format(function.__name__, ok))
    if not ok:
        if ok1:
            print("error {}: expected_exp ok but s**r != 1".format(function.__name__))
        elif ok2:
            print("error {}: s**r == 1 but not expected_exp".format(function.__name__))
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
    g2cx = QQx(g2cx_coeffs)/QQ(g2cx_denom)
    betax = QQx(betax_coeffs)/QQ(betax_denom)
    lambx = QQx(lambx_coeffs)/QQ(lambx_denom)
    
    p = ZZ(px(u))
    q = p
    r = ZZ(rx(u))
    c = ZZ(cx(u))
    y = ZZ(yx(u))
    t = ZZ(tx(u))
    tr = t
    tr0 = ZZ(tr0x(u))
    g2c = ZZ(g2cx(u))
    c2 = g2c
    beta = (betax(u))
    lamb = ZZ(lambx(u))

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
    print("u = {:#x}".format(u))
    print("ht= {:#x} hy= {:#x}".format(ht, hy))
    print("p = {:#x}".format(p))
    print("r = {:#x}".format(r))
    print("c = {:#x}".format(c))
    print("y = {:#x}".format(y))
    print("t = {:#x}".format(t))

    assert t_mod_r == tr0
    gr = gcd(p+1-tr0, p**2-p+1)
    assert (gr % r) == 0
    cgr = gr//r
    if cgr == 1 or cgr == -1:
        print("gcd(p+1-(tr mod r), p^2-p+1) = r, Scott's trick available")
    else:
        print("gcd(p+1-(tr mod r), p^2-p+1) = {}*r, Scott's trick not available".format(cgr))

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

    print("cost pairing")
    cost_final_exp_hard_bw6_bls24(u, ht, hy, t_mod_r_mod_u)
    cost_pairing_bw6_bls24(u, t_mod_r_mod_u)
    cost_pairing_opt_tate_bw6_bls24(u)

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

    test_bw6_bls24_r_subgroup_membersip_testing(E, r, c, t, u, omega)

    ok_c_omega, str_test = test_bw6_bls24_g1_mult_by_cofactor(E, r, c, t, u, omega, ht, hy, verbose=False)
    if ok_c_omega:
        print("with phi(x,y) = (omega*x, y):")
        print(str_test)
    ok_c_omega_, str_test_ = test_bw6_bls24_g1_mult_by_cofactor(E, r, c, t, u, -omega-1, ht, hy, verbose=False)
    if ok_c_omega_:
        print("with phi(x,y) = ((-omega-1)*x, y):")
        print(str_test_)
    if not ok_c_omega and not ok_c_omega_:
        print("both with phi(x,y) = (omega*x, y) or ((-omega-1)*x, y):")
        print(str_test)
        print(str_test_)

    xiD, bD = find_twist_curve_parameter_xi_ab(b, Fq, r, c2, D_twist=True)
    print("xiD = {}".format(xiD))
    ED = EllipticCurve([Fq(0), bD])
    orderD = ED.order()
    assert orderD == c2*r
    print("ED ok")
    Fq6D = Fq.extension(z**6 - xiD, names=('wD',),proof=False); (wD,) = Fq6D._first_ngens(1)
    E_Fq6D = EllipticCurve([Fq6D(0), Fq6D(b)])

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

    ok_c_omega, str_test = test_bw6_bls24_g2_mult_by_cofactor(EM, r, c2, t, u, omega, ht, hy, verbose=False)
    if ok_c_omega:
        print("with phi(x,y) = (omega*x, y) and M-twist:")
        print(str_test)
    ok_c_omega_, str_test_ = test_bw6_bls24_g2_mult_by_cofactor(EM, r, c2, t, u, -omega-1, ht, hy, verbose=False)
    if ok_c_omega_:
        print("with phi(x,y) = ((-omega-1)*x, y) and M-twist:")
        print(str_test_)
    if not ok_c_omega and not ok_c_omega_:
        print("both with phi(x,y) = (omega*x, y) or ((-omega-1)*x, y) and M-twist:")
        print(str_test)
        print(str_test_)
    ok_c_omega, str_test = test_bw6_bls24_g2_mult_by_cofactor(ED, r, c2, t, u, omega, ht, hy, verbose=False)
    if ok_c_omega:
        print("with phi(x,y) = (omega*x, y) and D-twist:")
        print(str_test)
    ok_c_omega_, str_test_ = test_bw6_bls24_g2_mult_by_cofactor(ED, r, c2, t, u, -omega-1, ht, hy, verbose=False)
    if ok_c_omega_:
        print("with phi(x,y) = ((-omega-1)*x, y) and D-twist:")
        print(str_test_)
    if not ok_c_omega and not ok_c_omega_:
        print("both with phi(x,y) = (omega*x, y) or ((-omega-1)*x, y) and D-twist:")
        print(str_test)
        print(str_test_)

    print("\ntest pairings with M-twist")

    test_g2_frobenius_eigenvalue_bw6(E_Fq6M, EM, r, c2, D_twist=False)
    
    test_optimal_ate_formula_bw6_bls24(E_Fq6M, EM, wM, u, r, t, c2, D_twist=False)

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
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_2naf,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_2naf,E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)

    test_double_line_h_a0_twist6_aklgl(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,EM,Fq6M,r,c,c2,D_twist=False)
    test_add_line_h_a0_twist6_aklgl_test(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl_with_z(E, EM, wM, D_twist=False)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_m6_twist(Fq6M)

    test_miller_function_ate_aklgl(E,EM,Fq6M,xiM,r,c,c2,t-1,D_twist=False)

    if t_mod_r_mod_u == 3:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_aklgl,E,EM,Fq6M,r,c,c2,u0,D_twist=False)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_aklgl_2naf,E,EM,Fq6M,r,c,c2,u0,D_twist=False)

    omega_ = -omega-1
    k_bls = 24
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fq6M, u, p, r, tr, c, c2, omega_, map_Fq6M_Fp6M, k_bls, D_twist=False, function_name=miller_function_tate_aklgl_a0_2_parallel)
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, EM, E_Fq6M, E_Fq6M, u, p, r, tr, c, c2, omega, map_Fq6M_Fp6M, k_bls, D_twist=False, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar)

    print("\ntest pairings with D-twist")

    test_g2_frobenius_eigenvalue_bw6(E_Fq6D, ED, r, c2, D_twist=True)
    
    test_optimal_ate_formula_bw6_bls24(E_Fq6D, ED, wD, u, r, t, c2, D_twist=True)

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
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_2naf,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_2naf,E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)

    test_double_line_h_a0_twist6_aklgl(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_double_line_h_a0_twist6_aklgl_no_div2(E,ED,Fq6D,r,c,c2,D_twist=True)
    test_add_line_h_a0_twist6_aklgl_test(E, ED, wD, D_twist=True)
    test_add_line_h_a0_twist6_aklgl(E, ED, wD, D_twist=True)
    test_add_line_h_a0_twist6_aklgl_with_z(E, ED, wD, D_twist=True)
    test_sparse_sparse_mult_d6_twist(Fq6D)
    test_sparse_mult_d6_twist(Fq6D)

    test_miller_function_ate_aklgl(E,ED,Fq6D,xiD,r,c,c2,t-1,D_twist=True)

    if t_mod_r_mod_u == 3:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
    else:
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_aklgl,E,ED,Fq6D,r,c,c2,u0,D_twist=True)
        test_bilinear_miller_loop_opt_ate_bw6_761_aklgl(miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_aklgl_2naf,E,ED,Fq6D,r,c,c2,u0,D_twist=True)

    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fq6D, u, p, r, tr, c, c2, omega_, map_Fq6D_Fp6D, k_bls, D_twist=True, function_name=miller_function_tate_aklgl_a0_2_parallel)
    test_miller_function_tate_a0_twist6_aklgl_2_multi_scalar(E, ED, E_Fq6D, E_Fq6D, u, p, r, tr, c, c2, omega, map_Fq6D_Fp6D, k_bls, D_twist=True, function_name=miller_function_tate_a0_twist6_aklgl_2_multi_scalar)

    print("test final exp hard part:")
    if t_mod_r_mod_u == 0:
        test_final_exp_hard_bw6_bls24(Fq6D, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u)
        test_final_exp_hard_bw6_bls24(Fq6M, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u)
        test_final_exp_hard_bw6_bls24(Fq6D, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt)
        test_final_exp_hard_bw6_bls24(Fq6M, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt)
    else:
        test_final_exp_hard_bw6_bls24(Fq6D, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u)
        test_final_exp_hard_bw6_bls24(Fq6M, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u)
        test_final_exp_hard_bw6_bls24(Fq6D, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt)
        test_final_exp_hard_bw6_bls24(Fq6M, u, ht, hy, r, final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt)

    print("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    
if __name__ == "__main__":
    arithmetic(False)
    #print("BLS24-315 BW6-672 small hy")
    #for v in testvector_bls24_315_bw6_672_672:
    #    test_curve(v)
    print("BLS24-315 BW6-633")
    v = testvector_bls24_315_bw6_640[5]
    test_curve(v)
    print("BLS24-315 BW6-672 hy = 0")
    v = testvector_bls24_315_bw6_672[83]
    test_curve(v)

    for v in testvector_bls24_315_bw6_672:
        if v['ht_hw'] <= 5 or v['ht_hw2naf'] <= 4:
            test_curve(v)
    print("BLS24-315 BW6-633 to BW6-640")
    for v in testvector_bls24_315_bw6_640:
        if v['pnbits'] <= 637:
            test_curve(v)
