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
from cost_pairing import cost_pairing_bls12

from test_scalar_mult import test_glv_scalar_mult_g1

def test_G2_endomorphism(E2, E_Fqd, D_twist=True):
    # the endomorphism on G2 is: twist o Frob o untwist, where
    # twist: E'(Fp2) -> E(Fqd), (x,y) -> (x * wD**2, y * wD**3)
    # Frob: E(Fpk) -> E(Fpk), (x, y) -> (x^p, y^p)
    # untwist: E(Fqd) -> E'(Fqd), (x,y) -> (x / wD**2, y / wD**3)
    Q2 = E2.random_element()
    Fqd = E_Fqd.base_field()
    omega = Fqd.gen(0)
    xi = -Fqd.modulus().constant_coefficient()
    d = Fqd.degree()
    assert xi == omega**d
    Fq = E2.base_field()
    Fp = Fq.base()
    p = Fp.characteristic()
    if D_twist:
        Q = E_Fqd([Q2[0]*omega**2, Q2[1]*omega**3])
    else:
        Q = E_Fqd([Q2[0]/omega**2, Q2[1]/omega**3])
    # because SageMath does not know about field representation isomorphisms...
    #Q_Fpk = E_Fpk([map_Fq6_Fp12_A(Q[0]), map_Fq6_Fp12_A(Q[1])])
    #Q_Fpk_frob = E_Fpk([(Q_Fpk[0]).frobenius(), (Q_Fpk[1]).frobenius()])
    #Q_Fqd_frob = E_Fqd([map_Fpk_Fqd(Q_Fpk_frob[0]), map_Fpk_Fqd(Q_Fpk_frob[1])])
    # .frobenius() is not available on a relative extension like Fqd
    Q_frob = E_Fqd([Q[0]**p, Q[1]**p])
    if D_twist:
        psiQ = ([Q_frob[0]/omega**2, Q_frob[1]/omega**3])
    else:
        psiQ = ([Q_frob[0]*omega**2, Q_frob[1]*omega**3])
    xi_1 = xi**((p-1)//d)
    xi_2 = xi**(2*(p-1)//d)
    xi_3 = xi**(3*(p-1)//d)
    if D_twist:
        assert psiQ[0] == Q2[0].frobenius() * xi_2
        assert psiQ[1] == Q2[1].frobenius() * xi_3
    else:
        assert psiQ[0] == Q2[0].frobenius() / xi_2
        assert psiQ[1] == Q2[1].frobenius() / xi_3
    # now check that psi o psi is -phi in E'(Fp2):
    # psi o psi is
    # (x^(p^2) * xi_2^p * xi_2, y^(p^2) * xi_3^p * xi_3)
    xi_p1 = (xi.frobenius()*xi)
    if Fq.degree() == 2 and d==6:
        assert xi_p1 in Fp
        w = xi_1.frobenius() * xi_1
        assert w**2 - w + 1 == 0
        beta = xi_2.frobenius() * xi_2
        assert beta**2 + beta + 1 == 0
        # check xi_3**(p+1) == -1
        assert xi_3.frobenius() * xi_3 == -1

    if D_twist:
        print("test G2 endormorphism with D-twist ok")
    else:
        print("test G2 endormorphism with M-twist ok")

def test_curve(u0):
    print("u0 = {:#x}".format(u0))
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # BLS12 polynomials
    px = (x**6 - 2*x**5 + 2*x**3 + x + 1)/3
    rx = x**4-x**2+1
    tx = x+1
    cx = (x-1)**2/3
    yx = (x-1)*(2*x**2 - 1)/3
    betax = x**5 - 3*x**4 + 3*x**3 - x + 1
    lambx = x**2 -1
    c2x = (x**8 - 4*x**7 + 5*x**6 - 4*x**4 + 6*x**3 - 4*x**2 - 4*x + 13)/9 # cofactor for G2

    cost_pairing_bls12(u0)
    
    # check formula for final exponentiation
    assert ((x**12-1) // cyclotomic_polynomial(12)) == (x**6-1)*(x**2+1)
    exponent_x = (px**4-px**2+1)//rx
    exponent = ZZ(exponent_x(u0))
    l3 = (x-1)**2
    l2 = l3*x
    l1 = l2*x-l3
    l0 = l1*x+3
    assert l0+px*(l1+px*(l2+px*l3)) == 3*exponent_x
    exponent_easy = (px**6-1)*(px**2+1)
    exponent_hard = 3*(px**4-px**2+1)//rx
    print("exponent_easy = (px**6-1)*(px**2+1)")
    print("exponent_hard = 3*(px**4-px**2+1)//rx = {}".format(exponent_hard))
    print("(p^4-p^2+1)//r:")
    print("exponent = {}".format(exponent))
    print("exponent is prime: {}\n".format(exponent.is_prime()))

    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    c2 = ZZ(c2x(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    b, E = find_curve_parameter_b(Fp, r, c)
    print("curve parameter b = {}".format(b))
    print("p = {} mod 4".format(p % 4))
    lambda_mod_r = ZZ(lambx(u0))
    beta_mod_p = Fp(betax(u0))
    # sextic twist over GF(p^2)
    #Fpz.<z> = Fp[]
    #Fp2.<i> = Fp.extension(z**2+1)
    #Fp2s.<s> = Fp2[]
    #Fp12.<j> = Fp2.extension(s**6 - (1+i))

    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    if (p % 4) == 3:
        Fp2 = Fp.extension(z**2 + 1, names=('i',)); (i,) = Fp2._first_ngens(1)
        a = -1
        print("Fp2 = Fp[x]/(x^2 + 1)")
    else:
        a = 2
        while not (z**2 - a).is_irreducible():
            a = a+1
        print("Fp2 = Fp[x]/(x^2 - {})".format(a))
        Fp2 = Fp.extension(z**2 - a, names=('i',)); (i,) = Fp2._first_ngens(1)
    Fq = Fp2
    Fp2s = Fp2['s']; (s,) = Fp2s._first_ngens(1)
    xiM, btwM = find_twist_curve_parameter_xi_ab(b, Fp2, r, c2, d=6, D_twist=False)
    EM = EllipticCurve([Fp2(0), Fp2(btwM)])
    Fp12M = Fp2.extension(s**6 - xiM, names=('wM',)); (wM,) = Fp12M._first_ngens(1)
    Fq6M = Fp12M
    E12M = EllipticCurve([Fp12M(0), Fp12M(b)])

    # note: s^6 = i+1 => (s^6-1)^2 = i^2 = -1 so s^12 - 2*s^6 + 2 = 0
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
    # s^6 - xi = 0 <=> s^6 - i0 = i1*i <=> (s^6 - i0)^2 = i1^2*a
    # resultant(s^6-xi, z^2 - a)
    Fp12M_A = Fp.extension((z**6 - i0M)**2 - i1M**2*a, names=('SM',)); (SM,) = Fp12M_A._first_ngens(1)
    E12M_A = EllipticCurve([Fp12M_A(0), Fp12M_A(b)])

    try:
        test_xiM = -Fp12M.modulus().constant_coefficient()
        print("xiM == -Fp12M.modulus().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fp12M.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiM = -Fp12M.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiM == -Fp12M.polynomial().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fp12M.polynomial().constant_coefficient() raised an error:\n{}".format(err))

    def map_Fp12M_Fp12M_A(x):
        # evaluate elements of Fp12M = Fp[i]/(i^2+1)[s]/(s^6-(i+1)) at i=S^6-1 and s=S
        # i <-> wM^6-1 = SM^6-1 and wM <-> SM
        #return sum([sum([yj*(SM**6-1)**j for j,yj in enumerate(xi.polynomial())]) * SM**i for i,xi in enumerate(x.list())])
        return sum([xi.polynomial()((SM**6-i0M)/i1M) * SM**e for e,xi in enumerate(x.list())])
    def map_Fq6M_Fp12M_A(x, aM):
        return sum([xi.polynomial()((aM**6-i0M)/i1M) * aM**e for e,xi in enumerate(x.list())])

    def map_Fp2_Fp12M_A(x):
        # evaluate elements of Fq=Fp[i] at i=wM^6-1 = SM^6-1
        return x.polynomial()((SM**6-i0M)/i1M)

    xiD, btwD = find_twist_curve_parameter_xi_ab(b, Fp2, r, c2, d=6, D_twist=True)
    ED = EllipticCurve([Fp2(0), Fp2(btwD)])
    Fp12D = Fp2.extension(s**6 - xiD, names=('wD',)); (wD,) = Fp12D._first_ngens(1)
    Fq6D = Fp12D
    E12D = EllipticCurve([Fp12D(0), Fp12D(b)]) # but un Sage, an elliptic curve over an extension field in two layers is not handled easily
    # note: s^6 = i+2 <=> s^6-2 = i => s^12 -4*s^6 + 4 = i^2 = -1 so s^12 - 4*s^6 + 5 = 0
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
    # s^6 - xiD = 0 <=> s^6 - i0D = i1D*i <=> (s^6 - i0D)^2 = i1D^2*a
    # resultant(s^6-xiD, z^2 - a)
    Fp12D_A = Fp.extension((z**6 - i0D)**2 - i1D**2*a, names=('SD',)); (SD,) = Fp12D_A._first_ngens(1)
    E12D_A = EllipticCurve([Fp12D_A(0), Fp12D_A(b)])

    try:
        test_xiD = -Fp12D.modulus().constant_coefficient()
        print("xiD == -Fp12D.modulus().constant_coefficient(): {}".format(xiD == test_xiD))
    except AttributeError as err:
        print("xiD = -Fp12D.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiD = -Fp12D.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiD == -Fp12D.polynomial().constant_coefficient(): {}".format(xiD == test_xiD))
    except AttributeError as err:
        print("xiD = -Fp12D.polynomial().constant_coefficient() raised an error:\n{}".format(err))

    def map_Fp12D_Fp12D_A(x):
        # evaluate elements of Fp12D = Fp[i]/(i^2+1)[s]/(s^6-(i+2)) at i=S^6-2 and s=S
        # i <-> s^6-2 = SD^6-2 and s <-> SD
        #return sum([sum([yj*(SD**6-2)**j for j,yj in enumerate(xi.polynomial())]) * SD**e for e,xi in enumerate(x.list())])
        return sum([xi.polynomial()((SD**6-i0D)/i1D) * SD**e for e,xi in enumerate(x.list())])
    def map_Fq6D_Fp12D_A(x, aD):
        return sum([xi.polynomial()((aD**6-i0D)/i1D) * aD**e for e,xi in enumerate(x.list())])
    def map_Fp2_Fp12D_A(x):
        # evaluate elements of Fq=Fp[i] at i=s^6-1 = S^6-1
        return x.polynomial()((SD**6-i0D)/i1D)

    print("test E (G1)")
    test_order(E,r*c)
    print("test E' (G2) M-twist")
    test_order(EM,r*c2)

    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(E12M_A,EM,Fp12M,map_Fq6M_Fp12M_A,r,c2,D_twist=False)
    test_g2_frobenius_eigenvalue_alt(E12M_A,EM,map_Fp2_Fp12M_A,r,c2,D_twist=False)
    print("test Frobenius map on G2 with D-twist")
    test_g2_frobenius_eigenvalue(E12D_A,ED,Fp12D,map_Fq6D_Fp12D_A,r,c2,D_twist=True)
    test_g2_frobenius_eigenvalue_alt(E12D_A,ED,map_Fp2_Fp12D_A,r,c2,D_twist=True)

    print("Test endomorphism on G2")
    test_G2_endomorphism(ED, E12D, D_twist=True)
    test_G2_endomorphism(EM, E12M, D_twist=False)

    print("test GLV on G1")
    test_glv_scalar_mult_g1(E, lambda_mod_r, beta_mod_p, r, c)
    
    print("test sparse-dense and sparse-sparse multiplications")
    test_sparse_mult_m6_twist(Fq6M)
    test_sparse_mult_d6_twist(Fq6D)
    test_sparse_sparse_mult_m6_twist(Fq6M)
    test_sparse_sparse_mult_d6_twist(Fq6D)

    print("test Miller M-twist")
    test_miller_function_tate(E, E12M, EM, r, c, c2, D_twist=False)
    test_miller_function_tate_2naf(E, E12M, EM, r, c, c2, D_twist=False)

    test_miller_function_ate(E,E12M,EM,r,c,c2,u0,D_twist=False)
    test_miller_function_ate_2naf(E,E12M,EM,r,c,c2,u0,D_twist=False)
    test_miller_function_ate_aklgl(E,EM,Fp12M,xiM,r,c,c2,u0,D_twist=False)
    test_miller_function_ate_2naf_aklgl(E,EM,Fp12M,xiM,r,c,c2,u0,D_twist=False)

    #P12 = E12(P)
    #Q12 = E12(Q2)
    #P12.ate_pairing(Q12, r, 12, t, q=p) # very very slow

    print("test E'' (G2) D-twist")
    test_order(ED,r*c2)
    print("test Miller D-twist")
    test_miller_function_tate(E, E12D, ED, r, c, c2, D_twist=True)
    test_miller_function_tate_2naf(E, E12D, ED, r, c, c2, D_twist=True)

    test_miller_function_ate(E,E12D,ED,r,c,c2,u0,D_twist=True)
    test_miller_function_ate_2naf(E,E12D,ED,r,c,c2,u0,D_twist=True)
    test_miller_function_ate_aklgl(E,ED,Fp12D,xiD,r,c,c2,u0,D_twist=True)
    test_miller_function_ate_2naf_aklgl(E,ED,Fp12D,xiD,r,c,c2,u0,D_twist=True)

    print("\nFinal exponentiation")
    ee = ((px**6-1)*(px**2+1)*3*(px**4-px**2+1)//rx)(u0)
    test_final_exp_easy_k12(Fp12D_A)
    test_final_exp_bls12(Fp12D_A,r,u0,expected_exponent=ee)
    test_final_exp_2naf_bls12(Fp12D_A,r,u0,expected_exponent=ee)

    test_final_exp_easy_k12(Fp12M_A)
    test_final_exp_bls12(Fp12M_A,r,u0,expected_exponent=ee)
    test_final_exp_2naf_bls12(Fp12M_A,r,u0,expected_exponent=ee)

    print("\npairing")
    test_ate_pairing_bls12_aklgl(E,EM,r,c,c2,u0,Fp12M,map_Fp12M_Fp12M_A,D_twist=False)
    test_ate_pairing_bls12_aklgl(E,ED,r,c,c2,u0,Fp12D,map_Fp12D_Fp12D_A,D_twist=True)

    test_ate_pairing_bls12_2naf_aklgl(E,EM,r,c,c2,u0,Fp12M,map_Fp12M_Fp12M_A,D_twist=False)
    test_ate_pairing_bls12_2naf_aklgl(E,ED,r,c,c2,u0,Fp12D,map_Fp12D_Fp12D_A,D_twist=True)

    print("")

def test_curve_bls12_377():
    # BLS12-377 seed
    u0 = ZZ(0x8508C00000000001)
    print("BLS12-377")
    r = (u0**4 - u0**2 + 1)
    c = (u0-1)**2//3
    t = u0 + 1
    p = c*r + t - 1
    assert p-1 == (2**46 * 3 * 7 * 13 * 499) * 53 * 409 * 2557 * 6633514200929891813 * 73387170334035996766247648424745786170238574695861388454532790956181
    print("p-1 = (u-1)*(u^5-u^4-u^3+u^2+u+2)/3\n    = (2**46 * 3 * 7 * 13 * 499) * (53 * 409 * 2557 * 6633514200929891813 * 73387170334035996766247648424745786170238574695861388454532790956181)")
    # rx-1 == (x - 1) * (x + 1) * x^2
    assert r-1 == (2**46 * 3 * 7 * 13 * 499) * (2 * 5 * 958612291309063373) * 9586122913090633729**2
    print("r-1 = (u-1)*(u+1)*u^2\n    = (2**46 * 3 * 7 * 13 * 499) * (2 * 5 * 958612291309063373) * 9586122913090633729**2")
    # co-factor of GT
    cofactorGT = (p**4 - p**2 + 1)//r
    cofactorGT == 73 * 733 * 9907133392668042525879507032505736467897843055444091915449145127341979335245297039623613430249313493610313967651813052567591541348837126946242723119127718146453611530031520691590326560874282582459918946440322689054712022836308086271423008860459900177148220797027656526576202902477337430474960665692030443666345467155806114621584830666273889800768834356323344510407388976077
    test_curve(u0)

def test_curve_bls12_379():
    # BLS12-379
    u0 = ZZ(0x9b04000000000001)
    test_curve(u0)

def test_curve_bls12_381():
    # BLS12-381 seed from Zcash https://electriccoin.co/blog/new-snark-curve/
    u0 = ZZ(-(2**63+2**62+2**60+2**57+2**48+2**16))
    test_curve(u0)

def test_curve_bls12_383():
    # BLS12-383 seed from https://eprint.iacr.org/2019/077 Mike Scott
    u0 = ZZ(0x10008000001001200)
    test_curve(u0)
    
if __name__ == "__main__":
    arithmetic(False)
    test_curve_bls12_377()
    test_curve_bls12_379()
    test_curve_bls12_381()
    test_curve_bls12_383()
