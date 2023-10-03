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
from pairing_kss18 import *
from test_pairing import *

from cost_pairing import cost_pairing_kss18

def test_miller_loop_opt_ate_kss18(E, E2, Fqd, Fpk, map_Fqd_Fpk, r, c, c2, u, D_twist=False):
    """Testing Miller loop of optimal ate pairing on KSS18 curve

    Testing the function miller_loop_opt_ate_kss18(Q2, P, u)
    where Q2 and P are both of order r and in E(Fpk) but in distinct subgroups
    To obtain valid Q2, first Q of order r is sampled from E2(Fq) then
    a map (D-twist or M-twist) is applied to Q to obtain Q2 in E_Fqd,
    and then finally in E(Fpk).

    INPUT:
    - `E`: elliptic curve over ground field GF(p) of order r*c
    - `E2`: quartic twist over GF(q) where q = p^{k/d} of order r*c2
    - `Fqd`: degree d extension over Fq, where q = p^{k/d}
    - `Fpk`: absolute extension over Fp, isomorphic to Fqd
    - `map_Fqd_Fpk`: the isomorphism map from Fqd to Fpk
    - `r`: prime integer, order of subgroup of E and E2
    - `c`: cofactor of E, so that # E(Fp) = r*c
    - `c2`: cofactor of E2, so that # E2(Fq) = r*c2, and q = p^{k/d}
    - `u`: the seed of parameters
    - `D_twist`: whether E2(Fq) is a D-twist or an M-twist of E(Fq)

    RETURN: true or False
    """
    ok = test_bilinearity_miller_loop_ate_absolute_extension(E, E2, Fqd, Fpk, map_Fqd_Fpk, r, c, c2, u, D_twist=D_twist, function_name=miller_loop_opt_ate_kss18, curve_a_is_0=True)
    print("test bilinearity of KSS18 optimal ate pairing: {}".format(ok))
    return ok
    
def test_optimal_ate_formula_kss18(E_Fpk, E_Fqd, E2, map_Fqd_Fpk, w, u, r, c2, D_twist=False):
    """Test u*Q + 3*pi(Q) - pi_4(Q) = 0 for Q in the trace-0 subgroup of E
    Test Q + u*pi^2(Q) + 2*pi^3(Q) = 0

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
    p2 = p**2
    p3 = p2*p
    p4 = p2**2
    while ok and i < 10:
        Q = c2*E2.random_element()
        while Q == E2(0) or r*Q != E2(0):
            Q = c2*E2.random_element()
        if D_twist:
            Q2 = psi_sextic_d_twist(Q, w)
        else:
            Q2 = psi_sextic_m_twist(Q, w)
        Qd = E_Fqd((Q2[0], Q2[1]))
        piQd = E_Fqd(((Qd[0])**p, (Qd[1])**p))
        pi2Qd = E_Fqd(((Qd[0])**p2, (Qd[1])**p2))
        pi3Qd = E_Fqd(((Qd[0])**p3, (Qd[1])**p3))
        pi4Qd = E_Fqd(((Qd[0])**p4, (Qd[1])**p4))
        ok1 = (u*Qd + 3*piQd - pi4Qd) == E_Fqd(0)
        ok2 = (Qd + u*pi2Qd + 2*pi3Qd) == E_Fqd(0)

        Qk = E_Fpk((map_Fqd_Fpk(Q2[0]), map_Fqd_Fpk(Q2[1])))
        piQk = E_Fpk(((Qk[0]).frobenius(), (Qk[1]).frobenius()))
        pi2Qk = E_Fpk(((Qk[0]).frobenius(2), (Qk[1]).frobenius(2)))
        pi3Qk = E_Fpk(((Qk[0]).frobenius(3), (Qk[1]).frobenius(3)))
        pi4Qk = E_Fpk(((Qk[0]).frobenius(4), (Qk[1]).frobenius(4)))
        ok3 = (u*Qk + 3*piQk - pi4Qk) == E_Fpk(0)
        ok4 = (Qk + u*pi2Qk + 2*pi3Qk) == E_Fpk(0)
        ok = ok1 and ok2 and ok3 and ok4
        i = i+1
    if D_twist:
        print("test optimal ate formula D_twist u*Q + 3*pi(Q) - pi_4(Q) = 0: {}".format(ok1 and ok2))
        print("test optimal ate formula D_twist Q + u*pi_2(Q) + 2*pi_3(Q) = 0: {}".format(ok3 and ok4))
    else:
        print("test optimal ate formula M_twist u*Q + 3*pi(Q) - pi_4(Q) = 0: {}".format(ok1 and ok2))
        print("test optimal ate formula M_twist Q + u*pi_2(Q) + 2*pi_3(Q) = 0: {}".format(ok3 and ok4))
    return ok

def test_final_exp_hard_kss18(Fpk, r, u, function_name=final_exp_hard_kss18, expected_exponent=None):
    # test final_exp_hard_kss18(m, u)
    # expected exponent is (p^6 - p^3 + 1)/r * (-3/49*x^2)
    p = Fpk.characteristic()
    if expected_exponent is None:
        expected_exponent = ((p**6 - p**3 + 1)//r) * (-3)*(u//7)**2
    else:
        assert expected_exponent % ((p**6 - p**3 + 1)//r) == 0
    ok_exp = True
    ok_r = True
    ok_inv = True
    i = 0
    while (ok_r and ok_inv and ok_exp) and i < 10:
        f0 = Fpk.random_element()
        #f1 = f0.frobenius(9)
        #f = f1/f0 # f^(p^9-1)
        f = final_exp_easy_k18(f0)
        g = function_name(f, u)
        ok_r = g**r == Fpk(1)
        ok_exp = g == f**expected_exponent
        ok_inv = g.frobenius(9) == 1/g
        i = i+1
    print("test {}: f^r == 1: {}, f == m^expected_exp: {}, f^(p^9) == 1/f: {}".format(function_name.__name__, ok_r, ok_exp, ok_inv))
    return ok_r and ok_exp and ok_inv

def test_curve(u0):
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)

    cost_pairing_kss18(u0)
    # KSS18 polynomials
    px = (x**8 + 5*x**7 + 7*x**6 + 37*x**5 + 188*x**4 + 259*x**3 + 343*x**2 + 1763*x + 2401)/21
    rx = (x**6 + 37*x**3 + 343)/343 # 343 = 7^3
    cx = (x**2 + 5*x + 7)*49/3
    tx = (x**4 + 16*x + 7)/7
    yx = (5*x**4 + 14*x**3 + 94*x + 259)/21 # Y such that T^2 - 4*P = -3*Y^2
    betax = (x**7 + 3*x**6 + 4*x**5 + 44*x**4 + 118*x**3 + 71*x**2 + 483*x + 1118)/24
    lambx = x**3 + 18
    D = 3 # discriminant (-D = -3)
    assert tx**2 - 4*px == -D*yx**2
    px2 = px**2
    tx2 = tx**2 - 2*px
    yx2 = yx*tx
    assert tx2**2 - 4*px2 == -D*yx2**2
    assert px2 == (tx2**2 + D*yx2**2)/4
    px3 = px*px2
    tx3 = tx*tx2 - px*tx
    assert (px3 + 1 - tx3) == (px+1-tx)*(px2-px+1 + tx*(px+1+tx))
    yx3 = yx * (tx**2-px)
    G2_order = (px3 + 1 - (tx3+D*yx3)/2)
    G2_order_= (px3 + 1 - (tx3-D*yx3)/2)
    print("G2_order % rx = 0: {}".format((G2_order % rx) == 0))
    print("G2_order_ % rx = 0: {}".format((G2_order_ % rx) == 0))
    c2x = G2_order_ // rx # cofactor for G2
    c2x = (x**18 + 15*x**17 + 96*x**16 + 409*x**15 + 1791*x**14 + 7929*x**13 + 27539*x**12 + 81660*x**11 + 256908*x**10 + 757927*x**9 + 1803684*x**8 + 4055484*x**7 + 9658007*x**6 + 19465362*x**5 + 30860595*x**4 + 50075833*x**3 + 82554234*x**2 + 88845918*x + 40301641)/27
    print("c2x = {}\n = {}".format(c2x, c2x.factor()))
    exponent_easy = (px**9-1)*(px**3+1) # (x^18-1)/(x^6-x^3+1)
    exponent_hard = (px**6-px**3+1)//rx
    print("exponent_easy = (px**9-1)*(px**3+1)")
    print("exponent_hard = 3*(px**6-px**3+1)//rx = {}".format(exponent_hard))

    # Daiki Hayashida, Kenichiro Hayasaka, and Tadanori Teruya
    # https://eprint.iacr.org/2020/875
    Tx = tx - 1
    k = 18
    Phi_k = cyclotomic_polynomial(k)
    h2x = Phi_k(Tx)//rx
    assert exponent_hard == cx * (px**2 + px*Tx + Tx**2) * (Tx**3 + px**3 - 1) + h2x
    print("Hayashida-Hayasaka-Teruya eprint 2020/875")
    print("exponent_hard == cx * (px^2 + px*(tx-1) + (tx-1)^2)*((tx-1)^3 + px^3 - 1) + Phi_18(tx-1)//rx")

    #https://eprint.iacr.org/2021/1309
    l0 = 147*x + 108*x**2 + 21*x**3 + 7*x**4 + 5*x**5 + x**6
    l1 = -686 -505*x - 98*x**2 - 35*x**3 - 25*x**4 - 5*x**5
    l2 = 6 - 133*x**2 - 98*x**3 - 19*x**4 - 7*x**5 - 5*x**6 - x**7
    l3 = 245*x + 181*x**2 + 35*x**3 + 14*x**4 + 10*x**5 + 2*x**6
    l4 = -343 - 254*x - 49*x**2 - 21*x**3 - 15*x**4 - 3*x**5
    l5 = 3 + 7*x**2 + 5*x**3 + x**4
    eex = l0 + l1*px + l2*px**2 + l3*px**3 + l4*px**4 + l5*px**5
    print("eex % exponent_hard = {}".format(eex % exponent_hard))
    print("eex // exponent_hard = {}".format(eex // exponent_hard))
    print("eex / (3/49 * x**2 * exponent_hard) = {}".format(eex // (3 * x**2 / 49 * exponent_hard)))
    print("eex == (3/49 * x**2 * exponent_hard): {}".format(eex == (3 * x**2 / 49 * exponent_hard)))
    assert eex == 3 * x**2 / 49 * exponent_hard

    l6 = x**2 + 5*x + 7
    assert l5 == x**2*l6 + 3
    assert l4 == -3*x*l5 - 49*l6
    assert l3 == 2*x**2*l5 + 35*x*l6
    assert l1 == 2*l4 + x*l5
    assert l0 == 2*l3 + x*l4
    assert l2 == -x*l0 + 2*l5

    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    c2 = ZZ(c2x(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    b, E = find_curve_parameter_b(Fp, r, c)
    #E = EllipticCurve([Fp(0), Fp(b)])
    print("KSS18-{} E: y^2 = x^3 {:+d} /Fp of {} bits".format(p.nbits(), b, p.nbits()))
    print("u = {:#x} {} bits".format(u0, u0.nbits()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("c = {:#x} {} bits".format(c, c.nbits()))
    print("y = {:#x} {} bits".format(y, y.nbits()))
    print("t = {:#x} {} bits".format(t, t.nbits()))
    print("p = {} mod 4".format(p % 4))
    print("p-1 has 2-valuation {}".format((p-1).valuation(2)))
    print("r-1 has 2-valuation {}".format((r-1).valuation(2)))
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    assert (p % 3) == 1
    a = -1
    while not (z**3 - a).is_irreducible():
        a = a+1
    print("Fp3 = Fp[x]/(x^3 - {})".format(a))
    Fp3 = Fp.extension(z**3 - a, names=('i',)); (i,) = Fp3._first_ngens(1)
    Fq = Fp3
    Fp3s = Fp3['s']; (s,) = Fp3s._first_ngens(1)
    xiM, btwM = find_twist_curve_parameter_xi_ab(b, Fp3, r, c2, d=6, D_twist=False)
    EM = EllipticCurve([Fp3(0), Fp3(btwM)])
    Fp18M = Fp3.extension(s**6 - xiM, names=('wM',)); (wM,) = Fp18M._first_ngens(1)
    #Fp18M = Fp3.extension(s**6 - (i+7), names=('wM',)); (wM,) = Fp18M._first_ngens(1)
    #xiM = i+7
    Fq6M = Fp18M
    E18M = EllipticCurve([Fp18M(0), Fp18M(b)])
    E_Fq3M = E18M

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
    # s^6 - xi = 0 <=> s^6 - i0 = i1*i <=> (s^6 - i0)^3 = i1^3*a
    Fp18M_A = Fp.extension((z**6 - i0M)**3 - i1M**3*a, names=('SM',)); (SM,) = Fp18M_A._first_ngens(1)
    E18M_A = EllipticCurve([Fp18M_A(0), Fp18M_A(b)])

    def map_Fp18M_Fp18M_A(x):
        # evaluate elements of Fp18M = Fp[i]/(i^3-3)[s]/(s^6-(i+7)) at i=S^6-7 and s=S
        # i <-> s^6-7 = SM^6-7 and s <-> SM
        #return sum([sum([yj*(SM**6-7)**j for j,yj in enumerate(xi.polynomial())]) * SM**e for e,xi in enumerate(x.list())])
        return sum([xi.polynomial()((SM**6-i0M)/i1M) * SM**e for e,xi in enumerate(x.list())])
    def map_Fq6M_Fp18M_A(x, aM):
        return sum([xi.polynomial()((aM**6-i0M)/i1M) * aM**e for e,xi in enumerate(x.list())])
    def map_Fp3_Fp18M_A(x):
        # evaluate elements of Fq=Fp[i] at i=s^6-7 = S^6-7
        return x.polynomial()((SM**6-i0M)/i1M)
    P = EM.random_element()
    print("EM order c2*r: {}".format(r*c2*P==0))

    xiD, btwD = find_twist_curve_parameter_xi_ab(b, Fp3, r, c2, d=6, D_twist=True)
    ED = EllipticCurve([Fp3(0), Fp3(btwD)])
    Fp18D = Fp3.extension(s**6 - xiD, names=('wD',)); (wD,) = Fp18D._first_ngens(1)
    Fq6D = Fp18D
    E18D = EllipticCurve([Fp18D(0), Fp18D(b)])
    E_Fq3D = E18D

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
    # s^6 - xiD = 0 <=> s^6 - i0D = i1D*i <=> (s^6 - i0D)^3 = i1D^3*a
    # resultant(s^6-xiD, z^3 - a)
    Fp18D_A = Fp.extension((z**6 - i0D)**3 - i1D**3*a, names=('SD',)); (SD,) = Fp18D_A._first_ngens(1)
    E18D_A = EllipticCurve([Fp18D_A(0), Fp18D_A(b)])

    P = ED.random_element()
    print("ED order c2*r:{}".format(r*c2*P==0)) # ok
    def map_Fp18D_Fp18D_A(x):
        # evaluate elements of Fp18D = Fp[i]/(i^3-3)[s]/(s^6-i) at i=S^6 and s=S
        # i <-> wD^6 = SD^6 and wD <-> SD
        #return sum([sum([yj*(SD**6)**j for j,yj in enumerate(xi.polynomial())]) * SD**i for i,xi in enumerate(x.list())])
        return sum([xi.polynomial()((SD**6-i0D)/i1D) * SD**e for e,xi in enumerate(x.list())])
    def map_Fq6D_Fp18D_A(x, aD):
        return sum([xi.polynomial()((aD**6-i0D)/i1D) * aD**e for e,xi in enumerate(x.list())])
    def map_Fp3_Fp18D_A(x):
        # evaluate elements of Fq=Fp[i] at i=wD^6 = SD^6
        return x.polynomial()((SD**6-i0D)/i1D)

    print("test E (G1)")
    test_order(E,r*c)
    print("test E' (G2) M-twist")
    test_order(EM,r*c2)
    print("test E' (G2) D-twist")
    test_order(ED,r*c2)
    
    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(E18M_A,EM,Fp18M,map_Fq6M_Fp18M_A,r,c2,D_twist=False)
    test_g2_frobenius_eigenvalue_alt(E18M_A,EM,map_Fp3_Fp18M_A,r,c2,D_twist=False)
    print("test Frobenius map on G2 with D-twist")
    test_g2_frobenius_eigenvalue(E18D_A,ED,Fp18D,map_Fq6D_Fp18D_A,r,c2,D_twist=True)
    test_g2_frobenius_eigenvalue_alt(E18D_A,ED,map_Fp3_Fp18D_A,r,c2,D_twist=True)

    print("test ate pairing")
    print("tests with M-twist")
    test_miller_function_ate_aklgl(E,EM,Fp18M,xiM,r,c,c2,t-1,D_twist=False,verbose=False)
    test_miller_function_ate_2naf_aklgl(E,EM,Fp18M,xiM,r,c,c2,t-1,D_twist=False,verbose=False)
    print("tests with D-twist")
    test_miller_function_ate_aklgl(E,ED,Fp18D,xiD,r,c,c2,t-1,D_twist=True,verbose=False)
    test_miller_function_ate_2naf_aklgl(E,ED,Fp18D,xiD,r,c,c2,t-1,D_twist=True,verbose=False)

    print("test optimal ate pairing formula")
    print("tests with M-twist")
    test_optimal_ate_formula_kss18(E18M_A, E_Fq3M, EM, map_Fp18M_Fp18M_A, wM, u0, r, c2, D_twist=False)
    print("tests with D-twist")
    test_optimal_ate_formula_kss18(E18D_A, E_Fq3D, ED, map_Fp18D_Fp18D_A, wD, u0, r, c2, D_twist=True)

    test_add_line_h_a0_twist6_aklgl_with_z(E, EM, wM, D_twist=False)
    test_add_line_h_a0_twist6_aklgl_with_z(E, ED, wD, D_twist=True)

    print("tests with M-twist")
    test_miller_loop_opt_ate_kss18(E, EM, Fp18M, Fp18M_A, map_Fp18M_Fp18M_A, r, c, c2, u0, D_twist=False)
    test_miller_loop_opt_ate_aklgl(miller_loop_opt_ate_aklgl_kss18, E, EM, Fp18M, xiM, map_Fp18M_Fp18M_A, r, c, c2, u0, D_twist=False)
    print("tests with D-twist")
    test_miller_loop_opt_ate_kss18(E, ED, Fp18D, Fp18D_A, map_Fp18D_Fp18D_A, r, c, c2, u0, D_twist=True)
    test_miller_loop_opt_ate_aklgl(miller_loop_opt_ate_aklgl_kss18, E, ED, Fp18D, xiD, map_Fp18D_Fp18D_A, r, c, c2, u0, D_twist=True)

    print("\nFinal exponentiation")
    ee = ((px**9-1)*(px**3+1)*(px**6-px**3+1)//rx)(u0)
    test_final_exp_easy_k18(Fp18D_A)
    test_final_exp_hard_kss18(Fp18D_A, r, u0)
    expected_exponent = ((p**6 - p**3 + 1)//r) * (-3)*(u0//7)**2
    test_final_exp_hard_kss18(Fp18D_A, r, u0, function_name=final_exp_hard_kss18_v0, expected_exponent=-expected_exponent)
    test_final_exp_hard_kss18(Fp18D_A, r, u0, function_name=final_exp_hard_kss18_v1, expected_exponent=-expected_exponent)

    test_final_exp_easy_k18(Fp18M_A)
    test_final_exp_hard_kss18(Fp18M_A, r, u0)
    test_final_exp_hard_kss18(Fp18M_A, r, u0, function_name=final_exp_hard_kss18_v0, expected_exponent=-expected_exponent)
    test_final_exp_hard_kss18(Fp18M_A, r, u0, function_name=final_exp_hard_kss18_v1, expected_exponent=-expected_exponent)

if __name__ == "__main__":
    # KSS18-348 seed
    arithmetic(False)
    #u0 = ZZ(2**44+2**22-2**9+2)
    #test_curve(u0)
    # seed with high 2-valuation of (r-1)
    #u0 = ZZ(0xc0c44000000)
    #test_curve(u0)
    # seed for 192-bit security level
    u0 = ZZ(0x12fffdfdfffffffffc000)
    test_curve(u0)
