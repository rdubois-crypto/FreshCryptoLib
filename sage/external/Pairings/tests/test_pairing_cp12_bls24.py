from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.arith.misc import XGCD, xgcd
from sage.arith.misc import valuation

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from pairing import *
from pairing_cocks_pinch import *
from test_pairing import *
from testvector_bls24_315_cp12d3_640 import testvector_bls24_315_cp12d3_640

from cost_pairing import cost_pairing_cp8_cp12, cost_tate_pairing_cp12_disc3
from test_pairing_cp12_bls12 import test_curve

if __name__ == "__main__":
    arithmetic(False)
    print("BLS24-315")
    for v in testvector_bls24_315_cp12d3_640:
        if v['pnbits'] <= 630:
            test_curve(v, k_bls = 24)
    print("##################################################")
