from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

def check_distorsion_map_bls(k, u0, b):
    # BLS polynomials
    assert (k % 3) == 0
    rx = cyclotomic_polynomial(k)
    tx = x+1
    if (k % 6) == 0:
        yx = (x-1)*(2*x**(k//6) - 1)/3
        # x**(k//6) is a primitive 6-th root of unity
        lambx = x**(k//6) - 1 # -x**(k//6) and x**(k//6) - 1 are primitive 3rd roots of unity
    else:
        yx = (x-1)*(2*x**(k//3) + 1)/3
        lambx = x**(k//3)
    px = (tx**2 + 3*yx**2)/4
    assert (px+1-tx) % rx == 0
    cx = (px+1-tx) // rx
    assert cx % (x-1)**2/3 == 0
    # yx <-> 1/sqrt(-3) mod rx 
    # 1/yx <-> sqrt(-3) mod rx
    # compute betax
    g, inv_yx, inv_px = xgcd(yx, px)
    # g might be an integer.
    assert inv_yx*yx + inv_px*px == g
    assert g in QQ
    inv_yx = inv_yx / QQ(g)
    betax = (-1 + tx * inv_yx)/2 % px
    assert (betax**2 + betax + 1) % px == 0
    
    p = ZZ(px(u0))
    r = ZZ(rx(u0))
    if r % 3 == 0:
        r = r // 3
        rx = rx/3
        cx = cx*3
    c = ZZ(cx(u0))
    t = ZZ(tx(u0))
    y = ZZ(yx(u0))
    Fp = GF(p, proof=False)
    E = EllipticCurve([Fp(0), Fp(b)])
    lambda_mod_r = ZZ(lambx(u0))
    beta_mod_p = Fp(betax(u0))

    print(E)
    print("BLS-{} defined over {}-bit prime p".format(k, p.nbits()))
    
    #Kp.<pi> = NumberField(x**2 -t*x+p) # (t - sqrt(-3)*y)/2
    #Kp = NumberField(x**2-t*x+p, names=('pi',)); (pi,) = Kp._first_ngens(1)
    #K.<w> = NumberField(x**2+x+1)  # (-1 + sqrt(-3))/2
    #K = NumberField(x**2+x+1, names=('w',)); (w,) = K._first_ngens(1)
    #O_K = K.maximal_order()
    # On the existence of distortion maps on ordinary elliptic curves
    # Denis Charles

    def psi(P):
        return P.curve()(P[0]*beta_mod_p, P[1])
    
    c0 = (u0-1)//3
    print("c0 = (u-1)/3 = {} = {}".format(c0, c0.factor()))
    for i in range(10):
        P = 3*r*E.random_element()
        while P == E(0) or not c0*P == 0:
            P = 3*r*E.random_element()
        psiP = psi(P)
        assert c0*psiP == E(0)
        # compute a Weil or Tate pairing?
        e_weil = P.weil_pairing(psiP, c0)
        assert e_weil**c0 == 1
        assert e_weil != 1
        assert P.weil_pairing(randint(1, abs(c0))*P, c0) == 1
        assert psiP.weil_pairing(randint(1, abs(c0))*psiP, c0) == 1
        # check MOV attack: computing a DL on the curve and in Fp is equivalent thanks to the Weil pairing
        a = randint(2, abs(c0))
        assert P.weil_pairing(psi(a*P), c0) == (P.weil_pairing(psiP, c0))**a
        assert psiP.weil_pairing(a*P, c0) == (psiP.weil_pairing(P, c0))**a

        # looks like it does not work with Tate pairing
        #e_tate = P.tate_pairing(psiP, c0, 1)
        #assert e_tate**c0 == 1
        #assert e_tate != 1
        #assert P.tate_pairing(randint(1, abs(c0))*P, c0, 1) == 1
        #assert psiP.tate_pairing(randint(1, abs(c0))*psiP, c0, 1) == 1

    for (ci, ei) in c0.factor():
        if ci > 20:
            continue
        print("subgroup of order {}".format(ci))
        P = (c0//ci)*3*r*E.random_element()
        while P == E(0):
            P = (c0//ci)*3*r*E.random_element()
        assert ci*P == E(0)
        psiP = psi(P)
        assert ci*psiP == E(0)
        P0 = P-psiP
        psiP0 = psi(P0)
        for i in range(1, ci+1):
            print("{}*P0 = {}".format(i, i*P0))
        print()
        for i in range(1, ci+1):
            print("{}*psiP0 = {}".format(i, i*psiP0))
        for i in range(1, ci+1):
            print("psi^2({}P0)+psi({}P)+P = {}".format(i,i, psi(psiP0)+psiP0+P0))
        


if __name__ == "__main__":
    arithmetic(False)
    #preparse("QQx.<x> = QQ[]")
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    # BLS12-381 seed
    u0 = ZZ(-(2**63+2**62+2**60+2**57+2**48+2**16))
    check_distorsion_map_bls(12, u0, 4)
    # BLS12-377
    check_distorsion_map_bls(12, ZZ(0x8508C00000000001), 1)
    # BLS-24
    check_distorsion_map_bls(24, ZZ(-0xeffff000), 5)
