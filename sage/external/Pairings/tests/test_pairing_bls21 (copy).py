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
from cost_pairing import cost_pairing_bls21


def test_curve(u0, b):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    k = 21
    D = 3
    rx = x**12 - x**11 + x**9 - x**8 + x**6 - x**4 + x**3 - x + 1
    tx = x + 1
    yx = (x - 1)/3 * (2*x**7 + 1)
    cx = (x - 1)**2/3 * (x**2 + x + 1)
    #px = (x^16 - 2*x^15 + x^14 + x^9 - 2*x^8 + x^7 + x^2 + x + 1)/3
    px = cx*rx + tx - 1
    assert px == (tx**2 + 3*yx**2)/4
    px7 = px**7
    tx7 = tx**7 - 7*px*tx**5 + 14*px**2*tx**3 - 7*px**3*tx
    yx7 = yx * (tx**6 - 5*px*tx**4 + 6*px**2*tx**2 - px**3)
    assert tx7**2 - 4*px7 == -D*yx7**2
    # now the 3-rd twist
    E2_order = px7+1+( 3*yx7+tx7)/2
    assert (E2_order % rx) == 0
    g2cx = E2_order // rx # 3 factors of degrees 2, 14, 84

    print("BLS21 u={:#x}".format(u0))
    cost_pairing_bls21(u0)
    # TODO:
    # build the finite field extension
    # define the curve and its twists
    # test the Miller loop, final exponentiation easy, final exponentiation hard
    print("")

def test_bls21_511a():
    u0 = ZZ(2**32 - 2**20 + 2**7 -2**4)
    assert u0 == ZZ(0xfff00070)
    b = -2
    test_curve(u0, b)

def test_bls21_511b():
    u0 = ZZ(-2**32+2**24-2**22+2**3)
    assert u0 == ZZ(-0xff3ffff8)
    b = -2
    test_curve(u0, b)

def test_bls21_511c():
    u0 = ZZ(-2**32+2**25+2**6+2)
    assert u0 == ZZ(-0xfdffffbe)
    b = 16
    test_curve(u0, b)

if __name__ == "__main__":
    arithmetic(False)
    test_bls21_511a()
    test_bls21_511b()
    test_bls21_511c()
