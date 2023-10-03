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
from cost_pairing import cost_pairing_bls24

from test_scalar_mult import test_glv_scalar_mult_g1

def test_curve(u0):
    print("u0 = {:#x}".format(u0))
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # the problem with BLS24 is the cost of the final exponentiation, and the size of G2.
    # BLS24 polynomials
    px = (x-1)**2*(x**8 - x**4 + 1)/3 + x
    rx = x**8-x**4+1
    tx = x+1
    cx = (x-1)**2/3
    yx = (x-1)*(2*x**4 - 1)/3
    betax = x**9 - 3*x**8 + 4*x**7 - 4*x**6 + 3*x**5 - 2*x**3 + 2*x**2 - x + 1
    lambx = x**4 -1
    c2x = (x**32 - 8*x**31 + 28*x**30 - 56*x**29 + 67*x**28 - 32*x**27 - 56*x**26 + 160*x**25 - 203*x**24 + 132*x**23 + 12*x**22 - 132*x**21 + 170*x**20 - 124*x**19 + 44*x**18 - 4*x**17 + 2*x**16 + 20*x**15 - 46*x**14 + 20*x**13 + 5*x**12 + 24*x**11 - 42*x**10 + 48*x**9 - 101*x**8 + 100*x**7 + 70*x**6 - 128*x**5 + 70*x**4 - 56*x**3 - 44*x**2 + 40*x + 100)/81 # cofactor for G2
    D = 3 # discriminant (-D = -3)
    k = 24

    cost_pairing_bls24(u0)

    exponent_easy = (px**12-1)*(px**4+1)
    exponent_hard = 3*(px**8-px**4+1)//rx
    print("exponent_easy = (px**12-1)*(px**4+1)")
    print("exponent_hard = 3*(px**8-px**4+1)//rx = {}".format(exponent_hard))
    
    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    c = ZZ(cx(u0))
    c2 = ZZ(c2x(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    b, E = find_curve_parameter_b(Fp, r, c)
    #E = EllipticCurve(Fp, [0, b])
    print("BLS{}-{} E: y^2 = x^3 {:+d} /Fp of {} bits".format(k, p.nbits(), b, p.nbits()))
    print("u = {:#x}".format(u0))
    print("p = {:#x} # {} bits, log_2 p^{} = {:.2f}, {} bits".format(p, p.nbits(), k, float((k*log(p)/log(2)).n()), ceil((k*log(p)/log(2)).n())))
    print("r = {:#x} # {} bits".format(r, r.nbits()))
    print("c = {:#x} # {} bits".format(c, c.nbits()))
    print("y = {:#x}".format(y))
    print("t = {:#x}".format(t))

    print("p = {} mod 4".format(p % 4))
    print("p-1 has 2-valuation {}".format((p-1).valuation(2)))
    print("r-1 has 2-valuation {}".format((r-1).valuation(2)))
    lambda_mod_r = ZZ(lambx(u0))
    beta_mod_p = Fp(betax(u0))

    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)
    # define Fp2 then Fp4
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
    Fp2w = Fp2['w']; (w,) = Fp2w._first_ngens(1)
    a2 = i
    while not (w**2 - a2).is_irreducible():
        a2 = a2+1
    print("Fp4 = Fp2[w]/(w^2 - {})".format(a2))
    Fp22 = Fp2.extension(w**2-a2 , names=('j',)); (j,) = Fp22._first_ngens(1)
    try:
        coeffs_a2 = a2.polynomial().list()
    except AttributeError as err:
        coeffs_a2 = a2.list()
    a20 = coeffs_a2[0]
    a21 = coeffs_a2[1]
    a20m = ZZ(a20)
    if abs(a20m - p) < abs(a20m):
        a20m = a20m - p
    a21m = ZZ(a21)
    if abs(a21m - p) < abs(a21m):
        a21m = a21m - p
    if a20m == 0:
        str_a2 = ""
    else:
        str_a2 = "{}".format(a20m)
    if a21m == 1:
        if len(str_a2) == 0:
            str_a2 = "i"
        else:
            str_a2 += "+i"
    elif a21m == -1:
        str_a2 += "-i"
    elif a21m != 0:
        if len(str_a2) == 0:
            str_a2 = "{}*i".format(a21m)
        else:
            str_a2 += "{:+}*i".format(a21m)
    print("F_p22 = Fp2[w]/(w^2 + {})".format(str_a2))
    # w^2 = a21*i+a20 <=> w^2-a20 = a21*i <=> (w^2-a20)^2 = a21*a <=> w^4 -2*a20*w^2 + a20^2 - a21*a = 0
    Fp4 = Fp.extension(z**4-2*a20*z**2+a20**2-a21*a, names=('ii',)); (ii,) = Fp4._first_ngens(1)
    print("Fp4 = Fp[z]/(z^4{:+d}*z^2{:+d})".format(int(-2*a20m), int(a20**2-a21*a)))
    Fq = Fp4
    Fp4s = Fp4['s']; (s,) = Fp4s._first_ngens(1)
    # find irreducible polynomial
    xiM, btwM = find_twist_curve_parameter_xi_ab(b, Fp4, r, c2, d=6, D_twist=False)
    EM = EllipticCurve(Fp4, [0, btwM])
    Fq6M = Fp4.extension(s**6 - xiM, names=('wM',)); (wM,) = Fq6M._first_ngens(1)
    E_Fq6M = E.base_extend(Fq6M)
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
    # w^6 = ii+4 <=> w^6-4 = ii <=> (w^6-4)^4 - 2*(w^6-4)^2 + 2 = 0
    # w^6 = i1M*ii + i0M <=> w^6-i0M = i1M*ii => (w^6-i0M)^2 = i1M^2*(a21*i+a20)
    #                                         => ((w^6-i0M)^2 - i1M^2*a20)^2 = i1M^4*a21^2 * a
    assert ((wM**6-i0M)**2 - a20*i1M**2)**2  == a*a21**2 * i1M**4
    polyM = ((z**6-i0M)**2 - a20*i1M**2)**2 - a*a21**2 * i1M**4
    Fp24M = Fp.extension(polyM, names=('SM',)); (SM,) = Fp24M._first_ngens(1)
    E_Fp24M = E.base_extend(Fp24M)
    print("Fq6M = Fq[s]/(s^6{:+d}{:+d}*ii)".format(int(-i0m), int(-i1m)))
    print("Fp24M = Fp[z]/({})".format(polyM))

    try:
        test_xiM = -Fq6M.modulus().constant_coefficient()
        print("xiM == -Fq6M.modulus().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fq6M.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiM = -Fq6M.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiM == -Fq6M.polynomial().constant_coefficient(): {}".format(xiM == test_xiM))
    except AttributeError as err:
        print("xiM = -Fq6M.polynomial().constant_coefficient() raised an error:\n{}".format(err))

    def map_Fq6M_Fp24M(x):
        return sum([xi.polynomial()((SM**6-i0M)/i1M) * SM**e for e,xi in enumerate(x.list())])
    def map_Fq6M_Fp24M(x, aM):
        return sum([xi.polynomial()((aM**6-i0M)/i1M) * aM**e for e,xi in enumerate(x.list())])
    def map_Fp4_Fp24M(x):
        # evaluate elements of Fq=Fp[ii] at ii=(s^6-i0M)/i1M
        return x.polynomial()((SM**6-i0M)/i1M)

    xiD, btwD = find_twist_curve_parameter_xi_ab(b, Fp4, r, c2, d=6, D_twist=True)
    ED = EllipticCurve(Fp4, [0, btwD])
    Fq6D = Fp4.extension(s**6 - xiD, names=('wD',)); (wD,) = Fq6D._first_ngens(1)
    E_Fq6D = E.base_extend(Fq6D)
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
    assert ((wD**6-i0D)**2 - a20*i1D**2)**2  == a*a21**2 * i1D**4
    polyD = ((z**6-i0D)**2 - a20*i1D**2)**2 - a*a21**2 * i1D**4
    Fp24D = Fp.extension(polyD, names=('SD',)); (SD,) = Fp24D._first_ngens(1)
    E_Fp24D = E.base_extend(Fp24D)
    print("Fq6D = Fq[s]/(s^6{:+d}{:+d}*ii)".format(int(-i0d), int(-i1d)))
    print("Fp24D = Fp[z]/({})".format(polyD))

    try:
        test_xiD = -Fq6D.modulus().constant_coefficient()
        print("xiD == -Fq6D.modulus().constant_coefficient(): {}".format(xiD == test_xiD))
    except AttributeError as err:
        print("xiD = -Fq6D.modulus().constant_coefficient() raised an error:\n{}".format(err))
    try:
        test_xiD = -Fq6D.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        print("xiD == -Fq6D.polynomial().constant_coefficient(): {}".format(xiD == test_xiD))
    except AttributeError as err:
        print("xiD = -Fq6D.polynomial().constant_coefficient() raised an error:\n{}".format(err))

    def map_Fq6D_Fp24D(x):
        # evaluate elements of Fq6D = Fp[ii]/(ii^4-2*a20*ii^2+a20^2-a21*a)[s]/(s^6-(i1D*ii+i0D)) at ii=(S^6-i0D)/i1D and s=S
        # ii <-> (wD^6-i0D)/i1D = (SD^6-i0D)/i1D and wD <-> SD
        return sum([xi.polynomial()((SD**6-i0D)/i1D) * SD**e for e,xi in enumerate(x.list())])
    def map_Fq6D_Fp24D(x, aD):
        return sum([xi.polynomial()((aD**6-i0D)/i1D) * aD**e for e,xi in enumerate(x.list())])
    def map_Fp4_Fp24D(x):
        return x.polynomial()((SD**6-i0D)/i1D)

    print("test E (G1)")
    test_order(E,r*c)
    print("test E' (G2) M-twist")
    test_order(EM,r*c2)
    print("test E' (G2) D-twist")
    test_order(ED,r*c2)
    
    print("test Frobenius map on G2 with M-twist")
    test_g2_frobenius_eigenvalue(E_Fp24M,EM,Fq6M,map_Fq6M_Fp24M,r,c2,D_twist=False)
    test_g2_frobenius_eigenvalue_alt(E_Fp24M,EM,map_Fp4_Fp24M,r,c2,D_twist=False)
    print("test Frobenius map on G2 with D-twist")
    test_g2_frobenius_eigenvalue(E_Fp24D,ED,Fq6D,map_Fq6D_Fp24D,r,c2,D_twist=True)
    test_g2_frobenius_eigenvalue_alt(E_Fp24D,ED,map_Fp4_Fp24D,r,c2,D_twist=True)
    
    print("test GLV on G1")
    test_glv_scalar_mult_g1(E, lambda_mod_r, beta_mod_p, r, c)

    print("test Miller M-twist")
    test_miller_function_tate(E, E_Fq6M, EM, r, c, c2, D_twist=False)
    test_miller_function_tate_2naf(E, E_Fq6M, EM, r, c, c2, D_twist=False)

    test_miller_function_ate(E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
    test_miller_function_ate_2naf(E,E_Fq6M,EM,r,c,c2,u0,D_twist=False)
    test_miller_function_ate_aklgl(E,EM,Fq6M,xiM,r,c,c2,u0,D_twist=False,verbose=False)
    test_miller_function_ate_2naf_aklgl(E,EM,Fq6M,xiM,r,c,c2,u0,D_twist=False,verbose=False)

    print("test Miller D-twist")
    test_miller_function_tate(E, E_Fq6D, ED, r, c, c2, D_twist=True)
    test_miller_function_tate_2naf(E, E_Fq6D, ED, r, c, c2, D_twist=True)

    test_miller_function_ate(E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
    test_miller_function_ate_2naf(E,E_Fq6D,ED,r,c,c2,u0,D_twist=True)
    test_miller_function_ate_aklgl(E,ED,Fq6D,xiD,r,c,c2,u0,D_twist=True,verbose=False)
    test_miller_function_ate_2naf_aklgl(E,ED,Fq6D,xiD,r,c,c2,u0,D_twist=True,verbose=False)
    
    print("\nFinal exponentiation")
    ee = ((px**12-1)*(px**4+1)*3*(px**8-px**4+1)//rx)(u0)
    test_final_exp_easy_k24(Fp24D)
    test_final_exp_bls24(Fp24D,r,u0,expected_exponent=ee)
    test_final_exp_2naf_bls24(Fp24D,r,u0,expected_exponent=ee)

    test_final_exp_easy_k24(Fp24M)
    test_final_exp_bls24(Fp24M,r,u0,expected_exponent=ee)
    test_final_exp_2naf_bls24(Fp24M,r,u0,expected_exponent=ee)
    
    print("\npairing")
    test_ate_pairing_bls24_aklgl(E,EM,r,c,c2,u0,Fq6M,map_Fq6M_Fp24M,D_twist=False)
    test_ate_pairing_bls24_aklgl(E,ED,r,c,c2,u0,Fq6D,map_Fq6D_Fp24D,D_twist=True)

    test_ate_pairing_bls24_2naf_aklgl(E,EM,r,c,c2,u0,Fq6M,map_Fq6M_Fp24M,D_twist=False)
    test_ate_pairing_bls24_2naf_aklgl(E,ED,r,c,c2,u0,Fq6D,map_Fq6D_Fp24D,D_twist=True)

if __name__ == "__main__":
    arithmetic(False)
    u_315 = ZZ(-0xbfcfffff) # https://eprint.iacr.org/2021/1359 Table 11
    u_317a = ZZ(0xe19c0001) # https://eprint.iacr.org/2021/1359 Table 11
    u_317b = ZZ(0xd9018000) # https://eprint.iacr.org/2022/586 Table 8
    #test_curve(u_315)
    #test_curve(u_317a)
    test_curve(u_317b)

    # BLS24-509 seed
    u1 = ZZ(-2**51 - 2**28 + 2**11 - 1)
    test_curve(u1)
    # BLS24-318 seed
    u0 = ZZ(-0xeffff000)
    # and E: y^2 = x^3 + 5 (a,b) = (0,5)
    assert u0 == ZZ(-2**32+2**28+2**12)
    test_curve(u0)
