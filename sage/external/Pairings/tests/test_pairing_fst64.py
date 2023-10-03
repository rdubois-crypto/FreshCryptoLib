from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.arith.misc import is_prime

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from pairing import *
from pairing_fst64 import *
from test_pairing import *
from cost_pairing import cost_pairing_fst64, cost_pairing_fst64_k20, cost_pairing_fst64_k28

from test_scalar_mult import test_glv_scalar_mult_g1
from test_pairing_bls12 import test_G2_endomorphism

def test_final_exp_hard_fst64(Fpk, k, u, r, function_name, final_exp_easy_k, expected_exponent=None):
    assert k == 20 or k == 28
    p = Fpk.characteristic()
    if expected_exponent is None:
        expected_exponent = cyclotomic_polynomial(k)(p)//r
    ok_exp = True
    ok_r = True
    ok_inv = True
    i = 0
    while (ok_r and ok_inv and ok_exp) and i < 10:
        f0 = Fpk.random_element()
        f = final_exp_easy_k(f0)
        g = function_name(f, u)
        ok_r = g**r == Fpk(1)
        ok_exp = g == f**expected_exponent
        ok_inv = g.frobenius((k//2)) == 1/g
        i += 1
    print("test {}: f^r == 1: {}, f == m^expected_exp: {}, f^(p^{}) == 1/f: {}".format(function_name.__name__, ok_r, ok_exp, k//2, ok_inv))
    return ok_r and ok_exp and ok_inv

def test_final_exp_hard_fst64_k20(Fpk, u, r, function_name=final_exp_hard_fst64_k20, expected_exponent=None):
    test_final_exp_hard_fst64(Fpk, 20, u, r, function_name, final_exp_easy_k20, expected_exponent)

def test_final_exp_hard_fst64_k28(Fpk, u, r, function_name=final_exp_hard_fst64_k28, expected_exponent=None):
    test_final_exp_hard_fst64(Fpk, 28, u, r, function_name, final_exp_easy_k28, expected_exponent)

def test_curve(k, u0, a=None):

    print("FST 6.4 k={} u0 = {:#x}".format(k, u0))
    u0 = ZZ(u0)
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # FST64 polynomials

    px = QQx(x**(k//2+2) - 2*x**(k//2+1) + x**(k//2) + x**2 + 2*x + 1)/4
    rx = QQx(cyclotomic_polynomial(k))
    tx = QQx(x+1)
    assert (px+1-tx) % rx == 0
    cx = (px+1-tx) // rx
    # c2x cofactor of G2
    d = 4
    c2x = get_g2_cofactor(k, d, px, rx, tx)
    yx = QQx(x**(k//4)*(x-1))
    assert 4*px == tx**2 + yx**2
    lambx = QQx(x**(k//4))

    cost_pairing_fst64(u0, k)

    # check formula for final exponentiation:
    assert ((x**k-1) // cyclotomic_polynomial(k)) % (x**(k//2)-1) == 0
    exponent_x = cyclotomic_polynomial(k)(px)//rx
    exponent = ZZ(exponent_x(u0))
    if is_prime(k // 4):
        assert ((x**k-1) // cyclotomic_polynomial(k)) == (x**(k//2)-1)*(x**2+1)
    if k == 20:
        # ((u-1)//2)^2 (u^2+1) (q + u) ((q^2 + u^2 - 1) (q^4 + u^4 - u^2 + 1) + (u^4 - u^2)) + 1
        l0 = QQx (((x-1)/2)**2)
        l1 = (x**2+1)
        l2 = (px + x)
        l3 = (px**2 + x**2 - 1)
        l5 = x**4 - x**2
        l4 = (px**4 + (l5 + 1))
        assert (l0 * l1 * l2 * (l3 * l4 + l5) + 1 == exponent_x)
        print("l0 * l1 * l2 * (l3 * l4 + l5) + 1 == exponent_x")
        exponent_easy = (px**(k//2)-1)*(px**2+1)
        exponent_hard = exponent_x
        print("exponent_easy = (px**{}-1)*(px**2+1)".format(k//2))
        print("exponent_hard = Phi_{}(px)//rx = ((x-1)/2)^2 * (x^2+1) * (p+x) * ((p^2+x^2-1) * (p^4+x^4-x^2+1) + x^4-x^2) + 1".format(k))
    print("Phi_{}(p)//r:".format(k))
    print("exponent = {}".format(exponent))
    print("exponent is prime: {}\n".format(exponent.is_prime()))

    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    c2 = ZZ(c2x(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    if a is None:
        a, E = find_curve_parameter_a(Fp, r, c)
    else:
        E = EllipticCurve(Fp, [a, 0])
    print("curve parameter a = {}".format(a))
    print("p = {} mod 4".format(p % 4))
    print("p = {} mod k".format(p % k))
    print("p = {} mod {}".format(p % (k//4), (k//4)))
    lambda_mod_r = ZZ(lambx(u0))
    #beta_mod_p = Fp(betax(u0))

    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    # extension of degree k//4
    e = k//4
    assert (p % k) == 1
    assert (p % 4) == 1
    assert (p % (k//4)) == 1
    # we should have p = 1 mod k//4 to be able to define the extension with a binomial x^(k//4)+a0, and because D=1, we have p=1 mod 4, consequently p = 1 mod k
    a0 = 2
    while not (z**e - a0).is_irreducible():
        a0 = a0+1
    print("Fp{} = Fp[x]/(x^{} - {})".format(e, e, a0))
    Fpe = Fp.extension(z**e - a0, names=('i',)); (i,) = Fpe._first_ngens(1)
    Fq = Fpe
    Fpes = Fpe['s']; (s,) = Fpes._first_ngens(1)
    xiM, atwM = find_twist_curve_parameter_xi_ab(a, Fpe, r, c2, d=4, D_twist=False)
    EM = EllipticCurve([Fpe(atwM), Fpe(0)])
    Fq4M = Fpe.extension(s**4 - xiM, names=('wM',)); (wM,) = Fq4M._first_ngens(1)
    Eq4M = EllipticCurve([Fq4M(a), Fq4M(0)])

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
    # s^4 - xi = 0 <=> s^4 - i0 = i1*i <=> (s^4 - i0)^e = i1^e*a0
    # resultant(s^4-xi, z^e - a0)

    FpkM = Fp.extension((z**4 - i0M)**e - i1M**e * a0, names=('SM',)); (SM,) = FpkM._first_ngens(1)
    EkM = EllipticCurve([FpkM(a), FpkM(0)])

    try:
        test_xiM = -Fq4M.modulus().constant_coefficient()
        print("xiM == -Fq4M.modulus().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fq4M.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiM = -Fq4M.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiM == -Fq4M.polynomial().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fq4M.polynomial().constant_coefficient() raised an error:\n{}".format(err))

    def map_Fq4M_FpkM(x, aM=None):
        if aM is None:
            # evaluate elements of Fq4M = Fp[i]/(i^e-a0)[s]/(s^4-xiM) at i=S^4-i0 and s=S
            return sum([xi.polynomial()((SM**4-i0M)/i1M) * SM**ei for ei,xi in enumerate(x.list())])
        else:
            return sum([xi.polynomial()((aM**4-i0M)/i1M) * aM**ei for ei,xi in enumerate(x.list())])

    def map_Fq_FpkM(x):
        # evaluate elements of Fq=Fp[i] at i0+i1*i=wM^4 <=> i = (wM^4 - i0)/i1
        return x.polynomial()((SM**4-i0M)/i1M)

    print("test E (G1)")
    test_order(E,r*c)
    print("test E' (G2) M-twist")
    test_order(EM,r*c2)

    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(EkM, EM, Fq4M, map_Fq4M_FpkM, r, c2, D_twist=False)
    test_g2_frobenius_eigenvalue_alt(EkM, EM, map_Fq_FpkM, r, c2, D_twist=False)

    print("Test endomorphism on G2")
    test_G2_endomorphism(EM, Eq4M, D_twist=False)

    #print("test GLV on G1")
    #test_glv_scalar_mult_g1(E, lambda_mod_r, beta_mod_p, r, c)

    print("test Miller M-twist")
    test_miller_function_tate(E, Eq4M, EM, r, c, c2, D_twist=False)
    test_miller_function_tate_2naf(E, Eq4M, EM, r, c, c2, D_twist=False)

    test_miller_function_ate(E, Eq4M, EM, r, c, c2, u0, D_twist=False)
    test_miller_function_ate_2naf(E, Eq4M, EM, r, c, c2, u0, D_twist=False)

    print("\nFinal exponentiation")
    #ee = ((px**(k//2)-1)*(px**2+1)*cyclotomic_polynomial(k)(px)//rx)(u0)
    if k == 20:
        test_final_exp_easy_k20(FpkM)
        test_final_exp_hard_fst64_k20(FpkM, u0, r, function_name=final_exp_hard_fst64_k20, expected_exponent=exponent)
        test_final_exp_hard_fst64_k20(FpkM, u0, r, function_name=final_exp_hard_fst64_k20_v1, expected_exponent=exponent)
    elif k == 28:
        test_final_exp_easy_k28(FpkM)
        test_final_exp_hard_fst64_k28(FpkM, u0, r, function_name=final_exp_hard_fst64_k28, expected_exponent=exponent)
        test_final_exp_hard_fst64_k28(FpkM, u0, r, function_name=final_exp_hard_fst64_k28_v1, expected_exponent=exponent)

    print("\ncurve done.\n")

def test_curve_fst64_k20_574():
    k = 20
    u0 = ZZ(-0xffffffc08001)
    a = 29
    test_curve(k, u0, a)

def test_curve_fst64_k20_636():
    k = 20
    u0 = ZZ(-0x21fffffffffbff)
    a = 1
    test_curve(k, u0, a)

def test_curve_fst64_k20_670():
    k = 20
    u0 = ZZ(-0xffefffffffffff)
    a = 1
    test_curve(k, u0, a)

def test_curve_fst64_k28_510():
    k = 28
    u0 = ZZ(2**32-2**25+2**22+2**15+1)
    assert u0 == ZZ(0xfe408001)
    a = 1
    # and p = 1 mod 7 and p=1 mod 28, possible to define GF(p^7) with a binomial
    test_curve(k, u0, a)

if __name__ == "__main__":
    arithmetic(False)
    test_curve_fst64_k20_574()
    test_curve_fst64_k20_636()
    test_curve_fst64_k20_670()
    test_curve_fst64_k28_510() # ok, p=1 mod 28
