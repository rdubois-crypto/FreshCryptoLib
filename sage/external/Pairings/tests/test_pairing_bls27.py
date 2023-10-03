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
from test_pairing import *
from cost_pairing import cost_pairing_bls27

def test_curve(u0, b=None, alphaM=None, alphaD=None):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    k = 27
    D = 3
    print("BLS27 u={:#x}".format(u0))
    cost_pairing_bls27(u0)

    # polynomials for BLS 27
    rx = cyclotomic_polynomial(27)
    tx = x+1
    yx = (tx-2)*(2*x**(k//3)+1)/3
    px = (tx**2 + 3*yx**2)/4
    lambx = QQx(x**(k//3))
    cx = (px+1-tx) // rx
    # because k is a power of 3:
    rx = rx/3
    cx = cx*3
    # c2x?
    tx9 = tx**9 - 9*px*tx**7 + 27*px**2*tx**5 - 30*px**3*tx**3 + 9*px**4*tx
    px9 = px**9
    yx9 = yx * (tx**2 - px) * (tx**6 - 6*px*tx**4 + 9*px**2*tx**2 - px**3)
    assert tx9**2-4*px9 == -3*yx9**2
    E2_orderx = px9 + 1 + (tx9+3*yx9)/2
    assert E2_orderx % rx == 0
    c2x = E2_orderx // rx

    print("test formula final exp hard")
    assert (px**18+px**9+1)//rx == 3+(x-1)**2*(px**2+px*x+x**2)*(px**6+px**3*x**3+x**6)*(px**9+x**9+1)
    assert (px**18+px**9+1)//rx == 3+(x-1)**2*(px**2+(px+x)*x)*(px**6+(px**3+x**3)*x**3)*(px**9+x**9+1)

    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    c2 = ZZ(c2x(u0))
    E2_order = ZZ(E2_orderx(u0))
    assert c2 * r == E2_order
    lambda_mod_r = ZZ(lambx(u0))
    Fp = GF(p, proof=False)
    if b is None:
        b, E = find_curve_parameter_b(Fp, r, c) #E/Fp: y^2 = x^3 + 15
    else:
        E = EllipticCurve(Fp, [0, b])

    print("BLS27-{} E: y^2 = x^3 {:+d} /Fp of {} bits".format(p.nbits(), b, p.nbits()))
    print("u = {:#x} {} bits".format(u0, u0.nbits()))
    print("p = {:#x} {} bits {} mod 3, {} mod {}, {} mod k={}".format(p, p.nbits(), p % 3, p % (k //3), (k // 3), p % k, k))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("c = {:#x} {} bits".format(c, c.nbits()))
    print("t = {:#x} {} bits".format(t, t.nbits()))
    print("c2 = {:#x} {} bits".format(c2, c2.nbits()))

    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    assert (p % 3) == 1
    e = k//3
    assert (p % e) == 1
    a = -1
    while not (z**e - a).is_irreducible():
        a = a+1

    print("Fp{} = Fp[x]/(x^{} - {}), p = {} mod {}".format(e, e, a, p%k, k))
    Fp9 = Fp.extension(z**e - a, names=('i',)); (i,) = Fp9._first_ngens(1)

    Fq = Fp9
    Fp9s = Fp9['s']; (s,) = Fp9s._first_ngens(1)

    print("find_twist_curve_parameter_xi_ab")
    xiM, btwM = find_twist_curve_parameter_xi_ab(b, Fp9, r, c2, d=3, D_twist=False)# Etw: y^2 = x^3 + b/xi^2 for cubic twist

    EM = EllipticCurve([Fp9(0), Fp9(btwM)])
    Fq3M = Fp9.extension(s**3 - xiM, names=('wM',)); (wM,) = Fq3M._first_ngens(1)
    Eq3M = EllipticCurve([Fq3M(0), Fq3M(b)])

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

    Fp27M = Fp.extension((z**3 - i0M)**9 - i1M**9*a, names=('SM',)); (SM,) = Fp27M._first_ngens(1)
    E27M = EllipticCurve([Fp27M(0), Fp27M(b)])

    def map_Fq3M_Fp27M(x):
        return sum([xi.polynomial()((SM**3-i0M)/i1M) * SM**e for e,xi in enumerate(x.list())])
    def map_Fq3M_Fp27M(x, aM = None):
        if aM is None:
            aM = SM
        return sum([xi.polynomial()((aM**3-i0M)/i1M) * aM**e for e,xi in enumerate(x.list())])
    def map_Fq_Fp27M(x):
        return x.polynomial()((SM**3-i0M)/i1M)

    print("test map_Fq3M_Fp27M")
    x0 = Fq3M.random_element()
    x1 = map_Fq3M_Fp27M(x0)

    print("test map_Fq_Fp27M")
    x0 = Fq.random_element()
    x1 = map_Fq_Fp27M(x0)

    P = EM.random_element()

    print("EM order c2*r: {}".format(r*c2*P==0))

    print("test E (G1)")
    test_order(E,r*c)
    print("test E' (G2) M-twist")
    test_order(EM,r*c2)

    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(E27M, EM, Fq3M, map_Fq3M_Fp27M, r, c2, D_twist=False)
    test_g2_frobenius_eigenvalue_alt(E27M, EM, map_Fq_Fp27M, r, c2, D_twist=False)

    print("tests lines with M-twist and cubic twist (denominators do not vanish easily)")
    test_double_line_cln_a0_cubic_twist(E, EM, Fq3M, D_twist=False)
    test_add_line_cln_a0_cubic_twist(E, EM, Fq3M, D_twist=False)
    test_add_line_cln_a0_cubic_twist_with_z(E, EM, Fq3M, D_twist=False)

    print("test ate pairing with M-twist and verticals because of cubic twist")
    test_bilinearity_miller_loop_ate_cln_a0_cubic_twist(E, Eq3M, EM, r, c, c2, r, D_twist=False, function_name=miller_function_tate_cln_a0_cubic_twist, Tate=True)
    test_bilinearity_miller_loop_ate_cln_a0_cubic_twist(E, Eq3M, EM, r, c, c2, r, D_twist=False, function_name=miller_function_tate_cln_a0_cubic_twist, Tate=False)
    test_bilinearity_miller_loop_ate_cln_a0_cubic_twist(E, Eq3M, EM, r, c, c2, u0, D_twist=False, function_name=miller_function_ate_cln_a0_cubic_twist)

    print("test ate pairing with M-twist and verticals, explicit cubic twist")
    test_bilinearity_miller_loop_ate_cln_a0_cubic_twist_explicit(E, Eq3M, EM, r, c, c2, u0, D_twist=False, function_name=miller_function_ate_cln_a0_cubic_twist)

    print("\nFinal exponentiation")
    test_final_exp_easy_k27(Fp27M)
    expected_exponent = (p**18 + p**9 + 1)//r
    test_final_exp_hard_bls27(Fp27M, u0, r, function_name=final_exp_hard_bls27, exponent_hard=expected_exponent)
    test_final_exp_hard_bls27(Fp27M, u0, r, function_name=final_exp_hard_bls27_zhang_lin, exponent_hard=expected_exponent)

    print("")

def test_curve_266():
    u = ZZ(-2**13-2**11-2**9+2**7+2**2+1)
    assert u == ZZ(-0x297b) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_276():
    u = ZZ(-2**14+2**11-2**9+2**6-2**4-2**2)
    assert u == ZZ(-0x39d4) # 1 mod 9
    b = 2
    test_curve(u, b)

def test_curve_290():
    u = ZZ(2**15-2**13-2**9+2**7+1)
    assert u == ZZ(0x5e81) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_422():
    u = ZZ(-2**21-2**18+2**7+2**5+2**2-1)
    assert u == ZZ(-0x23ff5d) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_423():
    u = ZZ(-2**21-2**18-2**16-2**14+2**7+1)
    assert u == ZZ(-0x253f7f) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_425a():
    u = ZZ(-2**21-2**19+2**12-2**9-2**7-2**4)
    assert u == ZZ(-0x27f290) # 1 mod 9
    b = -2
    test_curve(u, b)

def test_curve_425b():
    u = ZZ(-2**21-2**19+2**16+2**12+2**9-2**2-1)
    assert u == ZZ(-0x26ee05) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_425c():
    u = ZZ(2**21+2**19-2**15-2**12+2**7-2**4-2**2)
    assert u == ZZ(0x27706c) # 1 mod 9
    b = 2
    test_curve(u, b)

def test_curve_425d():
    u = ZZ(2**21+2**19-2**15+2**13+2**10-2**5+2**2)
    assert u == ZZ(0x27a3e4) # 1 mod 9
    b = 2
    test_curve(u, b)

def test_curve_425e():
    u = ZZ(2**21+2**19-2**14+2**12+2**10-2**7+2**4)
    assert u == ZZ(0x27d390)
    b = -2
    test_curve(u, b)

def test_curve_425f():
    u = ZZ(-2**21-2**19+2**13+2**9-2**5+2**3-2)
    assert u == ZZ(-0x27de1a)
    b = 2
    test_curve(u, b)

def test_curve_425g():
    u = ZZ(2**21+2**19+2**12-2**10+2**7-2**2-1)
    assert u == ZZ(0x280c7b) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_426():
    u = ZZ(-2**21-2**19-2**15+2**10+2**4+2**2+1)
    assert u == ZZ(-0x287beb) # 1 mod 9
    b = 1
    test_curve(u, b)

def test_curve_427():
    u = ZZ(2**21+2**19+2**17-2**14-2**8-2**5+2**3)
    assert u == ZZ(0x29bee8) # 1 mod 9
    b = -2
    alphaM = 2
    alphaD = 4
    test_curve(u, b, alphaM, alphaD)

if __name__ == "__main__":
    arithmetic(False)
    test_curve_426()
    test_curve_427()
    test_curve_425a()
    test_curve_423()
    test_curve_422()
    test_curve_266()
    test_curve_276()
    test_curve_290()
