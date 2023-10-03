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
from pairing_kss16 import *
from test_pairing import *
from cost_pairing import cost_pairing_kss16

def test_optimal_ate_formula(E_Fpk, E_Fqd, E2, map_Fqd_Fpk, w, u, r, c2, D_twist=False):
    """Test 2*Q + u*pi_3(Q) + pi_4(Q) = 0 for Q in the trace-0 subgroup of E

    INPUT:
    - `E_Fpk`: EllipticCurve instance defined over an absolute extension of Fp
    - `E2`: EllipticCurve instance defined over an absolute extension of Fp of degree p^{k/d}
    - `map_Fqd_Fpk`: map from the relative extension Fqd to the isomorphic absolute extension Fpk
    - `w`: the generator of Fpd (for the twisting map)
    - `u`: integer, parameter seed
    - `r`: prime integer, E(Fp) has order r*c
    - `c2`: integer, twist cofactor, E2 has order r*c2
    - `D_twist`: whether the twist is a D-twist or M-twist
    """
    ok = True
    i = 0
    p = E_Fpk.base_field().characteristic()
    p3 = p**3
    p4 = p3*p
    while ok and i < 10:
        Q = c2*E2.random_element()
        while Q == E2(0) or r*Q != E2(0):
            Q = c2*E2.random_element()
        if D_twist:
            Q2 = psi_sextic_d_twist(Q, w)
        else:
            Q2 = psi_sextic_m_twist(Q, w)
        Qd = E_Fqd((Q2[0], Q2[1]))
        pi3Qd = E_Fqd(((Qd[0])**p3, (Qd[1])**p3))
        pi4Qd = E_Fqd(((Qd[0])**p4, (Qd[1])**p4))
        ok = (2*Qd + u*pi3Qd + pi4Qd) == E_Fqd(0)

        Qk = E_Fpk((map_Fqd_Fpk(Q2[0]), map_Fqd_Fpk(Q2[1])))
        pi3Qk = E_Fpk(((Qk[0]).frobenius(3), (Qk[1]).frobenius(3)))
        pi4Qk = E_Fpk(((Qk[0]).frobenius(4), (Qk[1]).frobenius(4)))
        ok = ok and (2*Qk + u*pi3Qk + pi4Qk) == E_Fpk(0)
        i = i+1
    if D_twist:
        print("test optimal ate formula D_twist 2*Q + u*pi_3(Q) + pi_4(Q) = 0: {}".format(ok))
    else:
        print("test optimal ate formula M_twist 2*Q + u*pi_3(Q) + pi_4(Q) = 0: {}".format(ok))
    return ok
    
def test_final_exp_hard_kss16(Fpk, r, u, function_name=final_exp_hard_kss16, expected_exponent=None):
    """testing final_exp_hard_kss16(m, u)"""
    p = Fpk.characteristic()
    if expected_exponent is None:
        expected_exponent = ((p**8 + 1)//r) * 14*(u//5)**3
    else:
        assert expected_exponent % ((p**8 + 1)//r) == 0
    ok_exp = True
    ok_r = True
    ok_inv = True
    i = 0
    while (ok_r and ok_inv and ok_exp) and i < 10:
        f0 = Fpk.random_element()
        #f1 = f0.frobenius(8)
        #f = f1/f0 # f^(p^8-1)
        f = final_exp_easy_k16(f0)
        g = function_name(f, u)
        ok_r = g**r == Fpk(1)
        ok_exp = g == f**expected_exponent
        ok_inv = g.frobenius(8) == 1/g
        i = i+1
    print("test {}: f^r == 1: {}, f == m^expected_exp: {}, f^8 == 1/f: {}".format(function_name.__name__, ok_r, ok_exp, ok_inv))
    return ok_r and ok_exp and ok_inv

def test_curve(u0):
    assert (u0 % 70) in [25, 45]
    # we need u=+-25 mod 70 to ensure px, rx to be integers.
    u0 = ZZ(u0)
    cost_pairing_kss16(u0)    

    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    k = 16
    D = 1
    px = (x**10 + 2*x**9 + 5*x**8 + 48*x**6 + 152*x**5 + 240*x**4 + 625*x**2 + 2398*x + 3125)/980
    rx = (x**8 + 48*x**4 + 625)/61250 # 625 = 5^4, 61250 = 2*5^4*7^2
    tx = (2*x**5 + 41*x + 35)/35
    yx = (x**5 + 5*x**4 + 38*x + 120)/70 # Y such that T^2 - 4*P = -4*Y^2
    betax = (x**9-11*x**8-16*x**7-120*x**6-32*x**5-498*x**4-748*x**3-2740*x**2-3115*x-5651)/4018
    lambx = (x**4 + 24)/7 # sqrt(-1) mod R
    exponent = (px**8+1)//rx
    cx = 125 * (x**2 + 2*x + 5)/2 # C such that P+1-T = C*R

    # Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
    # https://eprint.iacr.org/2020/875
    Tx = tx-1
    Phi_k = cyclotomic_polynomial(k)
    h2x = Phi_k(Tx) // rx
    assert exponent == cx * (px+Tx) * (px**2 + Tx**2) * (px**4 + Tx**4) + h2x
    #print("Hayashida-Hayasaka-Teruya eprint 2020/875")
    #print("exponent_hard == cx * (px + tx-1) * (px^2 + (tx-1)^2) * (px^4 + (tx-1)^4) + Phi_16(tx-1)//rx")

    # for G2, compute #E(Fp4) then compute its 4-th twist
    #print("tx^2-4*px/yx^2 = {}".format((tx**2 - 4*px)/yx**2))
    D = 4
    assert tx**2 - 4*px == -D*yx**2
    px2 = px**2
    tx2 = tx**2 - 2*px
    yx2 = yx*tx
    assert tx2**2 - 4*px2 == -D*yx2**2
    assert px2 == (tx2**2 + D*yx2**2)/4
    px3 = px*px2
    px4 = px2**2
    tx3 = tx*tx2 - px*tx
    tx4 = tx*tx3 - px*tx2
    tx4 = (tx**4 - 4*px*tx**2 + 2*px**2)
    assert (px4 + 1 - tx4) == (px+1-tx)*(px+1+tx)*(px2+1+tx2)
    yx4 = yx*tx*(tx**2-2*px)
    assert px4 == (tx4**2 + D*yx4**2)/4
    # now the 4-th twist
    G2_order = px4+1-(-2*yx4)
    G2_order_= px4+1-( 2*yx4)
    #print("px4+1-(-yx4+tx4)/2 % rx == 0: {}".format((G2_order_ % rx) == 0))
    #print("px4+1-( yx4+tx4)/2 % rx == 0: {}".format((G2_order % rx) == 0))
    assert (G2_order % rx) == 0
    G2_cofactor = G2_order // rx # irreducible
    # print("G2_cofactor={}\n=({})/15059072".format(G2_cofactor,15059072*G2_cofactor))
    c2x = (x**32 + 8*x**31 + 44*x**30 + 152*x**29 + 550*x**28 + 2136*x**27 + 8780*x**26 + 28936*x**25 + 83108*x**24 + 236072*x**23 + 754020*x**22 + 2287480*x**21 + 5986066*x**20 + 14139064*x**19 + 35932740*x**18 + 97017000*x**17 + 237924870*x**16 + 498534968*x**15 + 1023955620*x**14 + 2353482920*x**13 + 5383092978*x**12 + 10357467880*x**11 + 17391227652*x**10 + 31819075896*x**9 + 65442538660*x**8 + 117077934360*x**7 + 162104974700*x**6 + 208762740168*x**5 + 338870825094*x**4 + 552745197960*x**3 + 632358687500*x**2 + 414961135000*x + 126854087873)/15059072

    # optimal ate pairing has Miller loop (f_{u,Q}(P) l_{[u]Q,\pi(Q)}(P))^{p^3} l_{Q,Q}(P)
    assert ((2 + x*px**3 + px**4) % rx) == 0
    # so 2*Q + u*pi_3(Q) + pi_4(Q) = 0
    
    # final exp hard
    # Loubna Ghammam
    # https://tel.archives-ouvertes.fr/tel-01469981v1
    # Utilisation des Couplages en Cryptographie asymétrique pour la micro-électronique.
    # PhD thesis, Universite de Rennes I, France
    # Eq (4.9) p. 114 Section 4.3.

    s = 14*(x/5)**3
    # for KSS16 curve we have u = 25,45 mod 70 <=> +/- 25 mod 70 --> 5 | u
    m_0 = 2*x**8 + 4*x**7 + 10*x**6 + 55*x**4 + 222*x**3 + 275*x**2
    m_1 = -4*x**7 - 8*x**6 - 20*x**5 - 75*x**3 - 374*x**2 - 375*x
    m_2 = -2*x**6 - 4*x**5 - 10*x**4 - 125*x**2 - 362*x - 625
    m_3 = -x**9 - 2*x**8 - 5*x**7 - 24*x**5 - 104*x**4 - 120*x**3 + 196
    m_4 = x**8 + 2*x**7 + 5*x**6 + 10*x**4 + 76*x**3 + 50*x**2
    m_5 = 3*x**7 + 6*x**6 + 15*x**5 + 100*x**3 + 368*x**2 + 500*x
    m_6 = -11*x**6 - 22*x**5 - 55*x**4 - 250*x**2 - 1116*x - 1250
    m_7 = 7*x**5 + 14*x**4 + 35*x**3 + 392
    assert m_0 + m_1*px + m_2*px**2 + m_3*px**3 + m_4*px**4 + m_5*px**5 + m_6*px**6 + m_7*px**7 == s*exponent

    B = (x+1)**2 + 4
    A = x**3*B + 56
    assert m_0 == 2*x**3*A + 55*x**2*B
    assert m_1 == -4*x**2*A - 75*x*B
    assert m_2 == -2*x*A - 125*B
    assert m_3 == -x**4*A -24*x**3*B + 196
    assert m_4 == x**3*A + 10*x**2*B
    assert m_5 == 3*x**2*A + 100*x*B
    assert m_6 == -11*x*A-250*B
    assert m_7 == 7*A

    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    c2 = ZZ(c2x(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    a, E = find_curve_parameter_a(Fp, r, c)
    #E = EllipticCurve([Fp(a), Fp(0)])
    print("KSS16-{} curve seed u = {:#x}".format(p.nbits(), u0.nbits()))
    print("u = {:#x} {} bits".format(u0, u0.nbits()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    assert (p % 4) == 1
    b = 1
    while not (z**4 - b).is_irreducible():
        b = b+1
    print("Fp4 = Fp[x]/(x^4-{})".format(b))
    Fp4 = Fp.extension(z**4 - b, names=('i',)); (i,) = Fp4._first_ngens(1)
    Fp4s = Fp4['s']; (s,) = Fp4s._first_ngens(1)

    xiD, atwD = find_twist_curve_parameter_xi_ab(a, Fp4, r, c2, d=4, D_twist=True)
    print("Fq4 = Fq[y]/(y^4-({}))".format(xiD))
    Fq4D = Fp4.extension(s**4 - xiD, names=('wD',)); (wD,) = Fq4D._first_ngens(1)
    ED = EllipticCurve([Fp4(atwD), Fp4(0)])
    E_Fq4D = EllipticCurve([Fq4D(a), Fq4D(0)])
    # s^4=xiD = xiD0 + i*xiD1 => (s^4-xiD0)^4 = i^4*xiD1^4 = b*xiD1^4
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
    Fp16D = Fp.extension((z**4 - i0D)**4 - b*i1D**4, names=('SD',)); (SD,) = Fp16D._first_ngens(1)
    E16D = EllipticCurve([Fp16D(a), Fp16D(0)])

    def map_Fq4D_Fp16D(X, aD=None):
        # evaluate elements of Fq4D = Fp[i]/(i^4-2)[s]/(s^4-i) at i=S^4 and s=S
        # i <-> wD^4 = SD^4 and wD <-> SD
        if aD is None:
            aD = SD
        return sum([xi.polynomial()((aD**4-i0D)/i1D) * aD**e for e,xi in enumerate(X.list())])

    def map_Fp4_Fp16D(X):
        # evaluate elements of Fq=Fp[i] at i=wD^4 = SD^4
        return X.polynomial()((SD**4-i0D)/i1D)

    xiM, atwM = find_twist_curve_parameter_xi_ab(a, Fp4, r, c2, d=4, D_twist=False) # should find 
    print("Fq4 = Fq[y]/(y^4-({}))".format(xiM))
    Fq4M = Fp4.extension(s**4 - xiM, names=('wM',)); (wM,) = Fq4M._first_ngens(1)
    EM = EllipticCurve([Fp4(atwM), Fp4(0)])
    E_Fq4M = EllipticCurve([Fq4M(a), Fq4M(0)])
    # s^4=xiM = xiM0 + i*xiM1 => (s^4-xiM0)^4 = i^4*xiM1^4 = b*xiM1^4
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
    Fp16M = Fp.extension((z**4 - i0M)**4 - b*i1M**4, names=('SM',)); (SM,) = Fp16M._first_ngens(1)
    E16M = EllipticCurve([Fp16M(a), Fp16M(0)])

    def map_Fq4M_Fp16M(X, aM=None):
        # evaluate elements of Fp16M = Fp[i]/(i^4-2)[s]/(s^4-(i+10)) at i=S^4-10 and s=S
        # i <-> s^4-10 = SM^4-10 and s <-> SM
        if aM is None:
            aM = SM
        return sum([xi.polynomial()((aM**4-i0M)/i1M) * aM**e for e,xi in enumerate(X.list())])

    def map_Fp4_Fp16M(X):
        # evaluate elements of Fq=Fp[i] at i=s^4-10 = S^4-10
        return X.polynomial()((SM**4-i0M)/i1M)

    print("test map_Fq4M_Fp16M")
    x0 = Fq4M.random_element()
    x1 = map_Fq4M_Fp16M(x0)
    
    print("test map_Fq4D_Fp16D")
    x0 = Fq4D.random_element()
    x1 = map_Fq4D_Fp16D(x0)


    print("test optimal ate pairing formula")
    test_optimal_ate_formula(E16M, E_Fq4M, EM, map_Fq4M_Fp16M, wM, u0, r, c2, D_twist=False)
    test_optimal_ate_formula(E16D, E_Fq4D, ED, map_Fq4D_Fp16D, wD, u0, r, c2, D_twist=True)

    print("test E (G1)")
    test_order(E,r*c)
    print("test ED (G2)")
    test_order(ED,r*c2)
    print("test EM (G2)")
    test_order(EM,r*c2)

    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(E16M,EM,Fq4M,map_Fq4M_Fp16M,r,c2,D_twist=False)
    test_g2_frobenius_eigenvalue_alt(E16M,EM,map_Fp4_Fp16M,r,c2,D_twist=False)
    print("test Frobenius map on G2 with D-twist")
    test_g2_frobenius_eigenvalue(E16D,ED,Fq4D,map_Fq4D_Fp16D,r,c2,D_twist=True)
    test_g2_frobenius_eigenvalue_alt(E16D,ED,map_Fp4_Fp16D,r,c2,D_twist=True)

    print("tests with D-twist")

    test_sparse_sparse_mult_d4_twist(Fq4D)
    test_sparse_mult_d4_twist(Fq4D)

    test_double_line_ate_cln_b0(E, ED, Fq4D, D_twist=True)
    test_add_line_ate_cln_b0(E, ED, Fq4D, D_twist=True)

    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist, E, ED, Fq4D, r, c, c2, t-1, D_twist=True, verbose=False)
    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist_acc_line, E, ED, Fq4D, r, c, c2, t-1, D_twist=True, verbose=False)

    test_double_line_j(E, ED, Fq4D, D_twist=True)
    test_double_line_affine_j(E, ED, Fq4D, D_twist=True)
    test_add_line_j(E, ED, Fq4D, D_twist=True)
    test_add_line_affine_j(E, ED, Fq4D, D_twist=True)
    test_double_line_cln_b0(E, ED, Fq4D, D_twist=True)
    test_add_line_cln_b0(E, ED, Fq4D, D_twist=True)
    test_add_line_cln_b0_with_z(E, ED, Fq4D, D_twist=True)

    test_miller_function_ate(E, E_Fq4D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_2naf(E, E_Fq4D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_csb(E, E_Fq4D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_cln_b0(E, E_Fq4D, ED, r, c, c2, t-1, D_twist=True)
    test_miller_function_ate_2naf_cln_b0(E, E_Fq4D, ED, r, c, c2, t-1, D_twist=True)

    test_bilinearity_miller_loop_ate_absolute_extension(E, ED, Fq4D, Fp16D, map_Fq4D_Fp16D, r, c, c2, u0, D_twist=True, function_name=miller_loop_opt_ate_kss16)
    test_bilinearity_miller_loop_ate_absolute_extension(E, ED, Fq4D, Fp16D, map_Fq4D_Fp16D, r, c, c2, u0, D_twist=True, function_name=miller_loop_opt_ate_kss16_v2)

    test_bilinearity_miller_loop_ate_absolute_extension(E, ED, Fq4D, Fp16D, map_Fq4D_Fp16D, r, c, c2, u0, D_twist=True, function_name=miller_loop_opt_ate_kss16_cln_b0)
    test_bilinearity_miller_loop_ate_absolute_extension(E, ED, Fq4D, Fp16D, map_Fq4D_Fp16D, r, c, c2, u0, D_twist=True, function_name=miller_loop_opt_ate_kss16_cln_b0_v2)

    print("tests with M-twist")

    test_sparse_sparse_mult_m4_twist(Fq4M)
    test_sparse_mult_m4_twist(Fq4M)

    test_double_line_ate_cln_b0(E, EM, Fq4M, D_twist=False)
    test_add_line_ate_cln_b0(E, EM, Fq4M, D_twist=False)

    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist, E, EM, Fq4M, r, c, c2, t-1, D_twist=False, verbose=False)
    test_miller_function_ate_cln_b0_quartic_twist(miller_function_ate_cln_b0_quartic_twist_acc_line, E, EM, Fq4M, r, c, c2, t-1, D_twist=False, verbose=False)

    test_double_line_j(E, EM, Fq4M, D_twist=False)
    test_double_line_affine_j(E, EM, Fq4M, D_twist=False)
    test_add_line_j(E, EM, Fq4M, D_twist=False)
    test_add_line_affine_j(E, EM, Fq4M, D_twist=False)
    test_double_line_cln_b0(E, EM, Fq4M, D_twist=False)
    test_add_line_cln_b0(E, EM, Fq4M, D_twist=False)
    test_add_line_cln_b0_with_z(E, EM, Fq4M, D_twist=False)

    test_miller_function_ate(E, E_Fq4M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_2naf(E, E_Fq4M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_csb(E, E_Fq4M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_cln_b0(E, E_Fq4M, EM, r, c, c2, t-1, D_twist=False)
    test_miller_function_ate_2naf_cln_b0(E, E_Fq4M, EM, r, c, c2, t-1, D_twist=False)

    test_bilinearity_miller_loop_ate_absolute_extension(E, EM, Fq4M, Fp16M, map_Fq4M_Fp16M, r, c, c2, u0, D_twist=False, function_name=miller_loop_opt_ate_kss16)
    test_bilinearity_miller_loop_ate_absolute_extension(E, EM, Fq4M, Fp16M, map_Fq4M_Fp16M, r, c, c2, u0, D_twist=False, function_name=miller_loop_opt_ate_kss16_v2)

    test_bilinearity_miller_loop_ate_absolute_extension(E, EM, Fq4M, Fp16M, map_Fq4M_Fp16M, r, c, c2, u0, D_twist=False, function_name=miller_loop_opt_ate_kss16_cln_b0)
    test_bilinearity_miller_loop_ate_absolute_extension(E, EM, Fq4M, Fp16M, map_Fq4M_Fp16M, r, c, c2, u0, D_twist=False, function_name=miller_loop_opt_ate_kss16_cln_b0_v2)

    print("Test Final Exp")
    test_final_exp_easy_k16(Fp16D)
    test_final_exp_easy_k16(Fp16M)
    expected_exponent = ((p**8 + 1)//r) * (-14)*(u0//5)**3 # -14/125*x^3
    test_final_exp_hard_kss16(Fp16M, r, u0, function_name=final_exp_hard_kss16, expected_exponent=expected_exponent)
    test_final_exp_hard_kss16(Fp16D, r, u0, function_name=final_exp_hard_kss16, expected_exponent=expected_exponent)
    test_final_exp_hard_kss16(Fp16M, r, u0, function_name=final_exp_hard_kss16_ghammam, expected_exponent=-expected_exponent)
    test_final_exp_hard_kss16(Fp16D, r, u0, function_name=final_exp_hard_kss16_ghammam, expected_exponent=-expected_exponent)
    test_final_exp_hard_kss16(Fp16M, r, u0, function_name=final_exp_hard_kss16_v2, expected_exponent=expected_exponent)
    test_final_exp_hard_kss16(Fp16D, r, u0, function_name=final_exp_hard_kss16_v2, expected_exponent=expected_exponent)

if __name__ == "__main__":
    arithmetic(False)
        
    u0_330 = -2**34+2**27-2**23+2**20-2**11+1
    u0_339 = 2**35-2**32-2**18+2**8+1
    u0_766 = 2**78-2**76-2**28+2**14+2**7+1
    
    args = sys.argv
    if len(args) == 2:
        if args[1] == "330":
            test_curve(u0_330)
        elif args[1] == "339":
            test_curve(u0_339)
        elif args[1] == "766":
            test_curve(u0_766)
    else:
        test_curve(u0_330)
        test_curve(u0_339)
