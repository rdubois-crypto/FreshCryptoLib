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
from pairing_bw13_bw19 import *
from test_pairing import *
from test_scalar_mult import test_glv_scalar_mult_g1

def trace_point(E_Fpk, Q, k):
    pQ = E_Fpk((Q[0].frobenius(), Q[1].frobenius()))
    p2Q = E_Fpk((pQ[0].frobenius(), pQ[1].frobenius()))
    traceQ = Q + pQ + p2Q
    piQ = p2Q
    for i in range(3, k):
        piQ = E_Fpk((piQ[0].frobenius(), piQ[1].frobenius()))
        traceQ += piQ
    return traceQ

def trace0(E_Fpk, Q, k):
    return k*Q - trace_point(E_Fpk, Q, k)

def test_cyclotomic_inverse(function_name, Fpk):
    i = 0
    ok = True
    while ok and i < 10:
        m = Fpk.random_element()
        m = m.frobenius() / m
        inv_m = function_name(m)
        ok = inv_m * m == 1
        i = i+1
    print("test {}: {}".format(function_name.__name__, ok))

def test_final_exp_hard(final_exp_hard_function, Fqk, u, r, expected_exp):
    """
    test the functions
    final_exp_hard_k13
    final_exp_hard_k19
    whose exponent is for k=13:
    (p^9-u^3*p^6+u^6*p^3-u^9)*u*(u + p)*((u^2-u+1)*(p*u^13 - p + u) - 3*p*u) + 3 + (u^2-u+1) * (p^3+1)*(p^6+1)*(p^2+p+1);
    """
    ok = True
    ok_r = True
    i = 0
    while ok and ok_r and i < 10:
        f = Fqk.random_element()
        g = final_exp_easy_k_prime(f)
        h = final_exp_hard_function(g, u)
        ok = g**expected_exp == h
        ok_r = h**r == 1
        i += 1
    print("test {}: {} g^r = 1: {}".format(final_exp_hard_function.__name__,ok, ok_r))
    return ok

def test_formula_ate_pairing(E_Fpk, r, ck, u, k):
    """
    the formula is u^2 + u*p + p^2 = 0 mod r for BW13 and BW19
    Test that [u^2]Q - [u]pi(Q) + pi^2(Q) = 0 on G2 = E(Fpk)[r] and such that Trace(Q) = 0
    """
    E0 = E_Fpk(0)
    p = E_Fpk.base_field().characteristic()
    p_mod_r = p % r
    if (r-p_mod_r) < p_mod_r:
        p_mod_r = p_mod_r - r
    ok = True
    i = 0
    while (i < 10) and ok:
        Q = E_Fpk.random_element()
        Q = ck*Q
        okr = r*Q == E0
        while Q == E0 or not okr:
            Q = E_Fpk.random_element()
            Q = ck*Q
            okr = r*Q == E0
        # make Q to have trace 0
        pQ = E_Fpk((Q[0].frobenius(), Q[1].frobenius()))
        p2Q = E_Fpk((pQ[0].frobenius(), pQ[1].frobenius()))
        traceQ = Q + pQ + p2Q
        piQ = p2Q
        for i in range(3, k):
            piQ = E_Fpk((piQ[0].frobenius(), piQ[1].frobenius()))
            traceQ += piQ
        Q_ = k*Q - traceQ
        print("Q_ != O: {}, r*Q_ == O: {}".format(Q_ != E0, r*Q_ == E0))
        trace_Q_= trace_point(E_Fpk, Q_, k)
        print("trace(Q_) = 0: {}".format(trace_Q_ == E0))
        pQ_ = E_Fpk((Q_[0].frobenius(), Q_[1].frobenius()))
        p2Q_ = E_Fpk((pQ_[0].frobenius(), pQ_[1].frobenius()))
        okq = (p_mod_r)*Q_ == pQ_
        ok = u**2*Q_ - u*pQ_ + p2Q_ == E0
        i = i+1
    print("test formula ate pairing: {} (u^2 - u*p + p^2) pQ = pi(Q): {} r*Q = O: {}".format(ok, okq, okr))
    return ok

def test_formulas_tate_pairing(E, r, c, u, omega, k):
    """
    QQ := Rationals();
    QQx<x> := PolynomialRing(QQ);
    for k in [13, 19] do
        rx := QQx ! CyclotomicPolynomial(6*k);
        tx := QQx ! (-x^(k+1) + x + 1);
        px := QQx ! ((x+1)^2 * (x^(2*k) - x^k + 1)/3 - x^(2*k+1));
        lambx := QQx ! (x^k - 1);
        M := Matrix(QQx, 2, 2, [rx, 0, -lambx, 1]);
        R := LLL(M);
        printf "%o\n", R;
        M := Matrix(QQx, 3, 3, [rx,0,0, -(px mod rx),1,0, -(px^2 mod rx),0,1]);
        R := LLL(M);
        printf "%o\n", R;
    end for;

    [   -x^11-x^10+x^8+x^7-x^5-x^4+x^2+x   x^12-x^10-x^9+x^7+x^6-x^4-x^3+x+1]
    [x^12+x^11-x^9-x^8+x^6+x^5-x^3-x^2+1     x^11+x^10-x^8-x^7+x^5+x^4-x^2-x]
    [x^2   -x   1]

    [   -x^17-x^16+x^14+x^13-x^11-x^10+x^8+x^7-x^5-x^4+x^2+x   x^18-x^16-x^15+x^13+x^12-x^10-x^9+x^7+x^6-x^4-x^3+x+1]
    [x^18+x^17-x^15-x^14+x^12+x^11-x^9-x^8+x^6+x^5-x^3-x^2+1     x^17+x^16-x^14-x^13+x^11+x^10-x^8-x^7+x^5+x^4-x^2-x]
    [x^2   -x   1]
    """
    if k == 13:
        a0 = -u**11-u**10+u**8+u**7-u**5-u**4+u**2+u
        a1 = u**12-u**10-u**9+u**7+u**6-u**4-u**3+u+1
        b0 = u**12+u**11-u**9-u**8+u**6+u**5-u**3-u**2+1
        b1 = u**11+u**10-u**8-u**7+u**5+u**4-u**2-u
    elif k == 19:
        a0 = -u**17-u**16+u**14+u**13-u**11-u**10+u**8+u**7-u**5-u**4+u**2+u
        a1 = u**18-u**16-u**15+u**13+u**12-u**10-u**9+u**7+u**6-u**4-u**3+u+1
        b0 = u**18+u**17-u**15-u**14+u**12+u**11-u**9-u**8+u**6+u**5-u**3-u**2+1
        b1 = u**17+u**16-u**14-u**13+u**11+u**10-u**8-u**7+u**5+u**4-u**2-u
    E0 = E(0)
    ok = True
    i = 0
    while (i < 10) and ok:
        P = E.random_element()
        P = c*P
        while P == E0 or not r*P == E0:
            P = E.random_element()
            P = c*P
        oka = a0*P + a1*bw6_phi(P, omega) == E0
        okb = b0*P + b1*bw6_phi(P, omega) == E0
        ok = oka and okb
        i = i+1
    print("test formulas Tate pairing: {} ({} {})".format(ok, oka, okb))
    return ok

def test_vertical_j_denom(E, E_Fpk):
    E0 = E(0)
    Ek0 = E_Fpk(0)
    Fpk = E_Fpk.base_field()
    list_a = [1, 2, 3, 11]
    # 1. test with S in E/Fp
    ok = True
    i = 0
    while ok and i < 10:
        S = E.random_element()
        for a in list_a:
            Sa = (S[0]*a**2, S[1]*a**3, a, a**2)
            vn, vd = vertical_j_denom(Sa, S)
            # it should be 0
            ok = ok and vn == 0
        i = i+1
    print("test vertical_j_denom(S, P) on E(Fp): {}".format(ok))
    ok = True
    i = 0
    while ok and i < 10:
        S = E_Fpk.random_element()
        for a in list_a:
            Sa = (S[0]*a**2, S[1]*a**3, a, a**2)
            vn, vd = vertical_j_denom(Sa, S)
            # it should be 0
            ok = ok and vn == 0
        i = i+1
    print("test vertical_j_denom(S, P) on E(Fpk): {}".format(ok))

def test_bilinear_miller_loop_with_frobenius(miller_loop_function, E, E_Fpk, r, T, c, ck, u, k):
    """Test the bilinearity of the arg function

    INPUT:
    -`miller_loop_function`: a function
    -`E`: elliptic curve over prime field Fp of order c*r (subgroup of order r is G1)
    -`E_Fpk`: E(Fpk) (just a change in the coordinate domain)
    -`r`: prime integer
    -`c`: curve cofactor (for G1)
    -`ck`: twist cofactor (for G2)
    -`u`: seed
    -`Tate`: input points are (P, Q) and loop length is r, instead of (Q, P) and T (=tr-1 for example) or u

    Functions:
    miller_function_tate_with_denom(P, Q, a, r)
    miller_loop_tate_with_denom(P, Q, a, r)
    miller_function_ate_with_denom(Q, P, a, T)
    miller_loop_ate_with_denom(Q, P, a, T)
    miller_loop_opt_ate_bw13_bw19_seq(Q, P, u, k)
    miller_loop_opt_ate_bw13_bw19(Q, P, u, k)
    """
    Tate = miller_loop_function is miller_loop_tate_with_denom or miller_loop_function is miller_function_tate_with_denom
    E0 = E(0)
    Ek0 = E_Fpk(0)
    Fpk = E_Fpk.base_field()
    P = c * E.random_element()
    while P == E0 or r*P != E0:
        P = c * E.random_element()
    Q2 = ck*E_Fpk.random_element()
    while Q2 == Ek0 or r*Q2 != Ek0:
        Q2 = ck * E_Fpk.random_element()
    # mow apply the trace map formula
    Q = trace0(E_Fpk, Q2, k)
    if miller_loop_function is miller_loop_tate_with_denom:
        f = miller_loop_function(P, Q, 0, r)
    elif miller_loop_function is miller_loop_ate_with_denom:
        f = miller_loop_function(Q, P, 0, T)
    elif miller_loop_function is miller_loop_opt_ate_bw13_bw19 or miller_loop_function is miller_loop_opt_ate_bw13_bw19_seq:
        f = miller_loop_function(Q, P, u, k)
    elif miller_loop_function is miller_function_tate_with_denom:
        m, md, S = miller_function_tate_with_denom(P, Q, 0, r, m0=1)
        f = m/md
        print("S = [r]P == O: {}".format(S[2] == 0))
    elif miller_loop_function is miller_function_ate_with_denom:
        m, md, S = miller_function_ate_with_denom(Q, P, 0, T, m0=1)
        f = m/md
        if S[2] != 0 and S[3] != 0:
            SS = E_Fpk((S[0]/S[3], S[1]/(S[2]*S[3])))
            print("S == [T]Q: {}".format(SS == T*Q))
        else:
            print("S = O, [T]Q = O: {}".format(T*Q == E_Fpk(0)))

    assert (Fpk.cardinality()-1) % r == 0
    exponent = (Fpk.cardinality()-1) // r
    assert exponent % r != 0
    g = f**exponent
    assert g != 1 and g**r == 1
    ok = True
    okS = True
    ok1 = True
    bb = 1
    T_Q = T*Q
    while ok and ok1 and okS and bb < 4:
        aa = 1
        while ok and okS and ok1 and aa < 4:
            if miller_loop_function is miller_loop_tate_with_denom:
                fij = miller_loop_function(aa*P, bb*Q, 0, r)
            elif miller_loop_function is miller_loop_ate_with_denom:
                fij = miller_loop_function(aa*Q, bb*P, 0, T)
            elif miller_loop_function is miller_loop_opt_ate_bw13_bw19 or miller_loop_function is miller_loop_opt_ate_bw13_bw19_seq:
                fij = miller_loop_function(aa*Q, bb*P, u, k)
            elif miller_loop_function is miller_function_tate_with_denom:
                m, md, S = miller_function_tate_with_denom(aa*P, bb*Q, 0, r, m0=1)
                okS = S[2] == 0
                fij = m/md
            elif miller_loop_function is miller_function_ate_with_denom:
                m, md, S = miller_function_ate_with_denom(aa*Q, bb*P, 0, T, m0=1)
                if S[2] != 0 and S[3] != 0:
                    SS = E_Fpk((S[0]/S[3], S[1]/(S[2]*S[3])))
                    okS = SS == aa*T_Q
                else:
                    okS = T_Q[2] == 0
                fij = m/md
            gij = fij**exponent
            ok1 = gij != 1 and gij**r == 1
            if not okS:
                if Tate:
                    print("error S is not O (Tate)")
                else:
                    print("error S is not [T]Q (ate)")
            if not ok1:
                print("error gij == 1: {}, gij^r == 1: {}".format(gij != 1, gij**r == 1))
            ok = gij == g**(aa*bb)
            if not ok:
                print("error gij ! g^(a*b), a={}, b={}".format(aa,bb))
                #print("gij = {}\n1/gij = {}".format(gij, 1/gij))
                #print("g = {}\ng^(a*b) = {}".format(g, g**(aa*bb)))
            aa += 1
        bb += 1
    if okS:
        if Tate:
            print("S is ok: O (Tate)")
        else:
            print("S is ok: [T]Q (ate)")
    print("test bilinearity {}: {}".format(miller_loop_function.__name__, ok))
    return ok

def test_curve(u, k):
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    rx = QQx(cyclotomic_polynomial(6*k))
    tx = QQx(-x**(k+1) + x + 1)
    px = QQx((x+1)**2 * (x**(2*k) - x**k + 1)/3 - x**(2*k+1))
    yx = QQx((x**(k+1) - 2*x**k + x + 1)/3)
    cx = QQx((x**2 - x + 1)**2/3)
    assert px == (tx**2 + 3*yx**2)/4
    assert (px + 1 - tx) % rx == 0
    assert cx * rx == px + 1 - tx
    # find an alternative to (tx-1) for the ate pairing
    Tx = tx-1
    for i in [j for j in range(2, k) if gcd(j, k) == 1]:
        Txi = ((tx-1)**i) % rx
        if Txi.degree() < Tx.degree():
            Tx = Txi
            Tx_exponent = i
    lambx = QQx(x**k - 1)
    lambx_ = -lambx - 1
    assert (lambx**2 + lambx + 1) % rx == 0
    assert (lambx_**2 + lambx_ + 1) % rx == 0
    #betax = QQx((x**2+x+1)*(x**6+x**3+1) * (x**18+x**9+1) + 2*x**26 + 3*(x**14 * (x+1)*(x**8+x**4+1)))
    if k == 13:
        betax = QQx([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 1])
    elif k == 19:
        betax = QQx([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 2, -1, 1])

    betax_ = -betax - 1
    assert (betax**2 + betax + 1) % px == 0
    assert (betax_**2 + betax_ + 1) % px == 0

    # the trace and the cofactor for G2
    t2x = tx**2 - 2*px
    y2x = tx*yx
    p2x = px**2
    tk_1 = t2x ; tk_2 = tx
    for i in range(3, k+1):
        tk = tx*tk_1 - px*tk_2
        tk_2 = tk_1
        tk_1 = tk
    tkx = tk
    pkx = px**k
    assert (pkx+1-tkx) % (cx * rx**2) == 0
    ckx = (pkx+1-tkx)//rx**2

    p = ZZ(px(u))
    r = ZZ(rx(u))
    c = ZZ(cx(u))
    t = ZZ(tx(u))
    T = ZZ(Tx(u))
    y = ZZ(yx(u))
    ck = ZZ(ckx(u))
    exponent_hard = cyclotomic_polynomial(k)(px) // rx
    ex = exponent_hard

    Fp = GF(p, proof=False)
    b, E = find_curve_parameter_b(Fp, r, c)

    print("BW{}-{} curve u={:#x}".format(k, p.nbits(), u))
    print("curve parameter b = {}".format(b))
    print("p = {} mod 4".format(p % 4))
    print("p = {} mod 3".format(p % 3))
    print("p = {} mod k={}".format(p % k, k))
    print("r is prime: {}".format(r.is_prime()))
    print("r-1 has 2-adicity {}".format((r-1).valuation(2)))
    print("p-1 has 2-adicity {}".format((p-1).valuation(2)))
    print("exponent_easy = (px-1)")
    print("exponent_hard = Phi_{}(px)/rx has degree {}".format(k, ex.degree()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("tr= {:#x} {} bits".format(t, t.nbits()))
    print("T = {:#x} {} bits, = (t-1)^{} mod r, = {}".format(T, T.nbits(), Tx_exponent, Tx))
    print("c = {:#x} {} bits".format(c, c.nbits()))
    print("y = {:#x} {} bits".format(y, y.nbits()))
    print("u^2-u*p+p^2 = {:#x}".format(u**2-u*p+p**2))
    print("u^2-u*p+p^2 = {} mod r".format((u**2-u*p+p**2) % r))
    print("(u^2-u*p+p^2)//r = {:#x}".format((u**2-u*p+p**2) // r))

    lambda_mod_r = ZZ(lambx(u))
    beta_mod_p = Fp(betax(u))
    omega = beta_mod_p
    omega_ = -omega-1
    Fpz = Fp['z']; (z,) = Fpz._first_ngens(1)

    if (p % k) == 1:
        # find an extension with a binomial polynomial
        a = -1
    else:
        a = z-1
    while not (z**k - a).is_irreducible():
        a = a+1
    Fpk = Fp.extension(z**k - a, names=('w',)); (w,) = Fpk._first_ngens(1)
    print("Fp{} = Fp[x]/(x^{} - {})".format(k, k, a))
    E_Fpk = EllipticCurve([Fpk(0), Fpk(b)])

    print("test E (G1)")
    test_order(E, r*c)
    print("test E (G2)")
    test_order(E_Fpk,r**2 * ck)

    if k == 13:
        test_cyclotomic_inverse(cyclotomic_inverse_k13, Fpk)
    elif k == 19:
        test_cyclotomic_inverse(cyclotomic_inverse_k19, Fpk)

    ee = ex(u)
    if k == 13:
        cex = 3*(x**12 + x**11 - x**9 - x**8 + x**6 + x**5 - x**3 - x**2 + 1)
        ce = cex(u)
        expected_exp = ce * ee
        test_final_exp_hard(final_exp_hard_k13, Fpk, u, r, expected_exp)
    elif k == 19:
        cex = 3*(x**18+x**17-x**15-x**14+x**12+x**11-x**9-x**8+x**6+x**5-x**3-x**2+1)
        ce = cex(u)
        expected_exp = ce * ee
        test_final_exp_hard(final_exp_hard_k19, Fpk, u, r, expected_exp)

    P = c*E.random_element()
    psiP = bw6_phi(P, omega)
    lambdaP = lambda_mod_r * P
    if not lambdaP == psiP:
        omega, omega_ = omega_, omega
        print("swapping omega, omega_ wrt lambda, lambda_")
        psiP = bw6_phi(P, omega)
        if not lambdaP == psiP:
            print("problem of omega and lambda")
            return

    test_double_line_j_no_twist(E, str_E=" on E(Fp)")
    test_double_line_j_no_twist(E_Fpk, str_E=" on E(Fpk)")
    test_add_line_j_no_twist(E, str_E=" on E(Fp)")
    test_add_line_j_no_twist(E_Fpk, str_E=" on E(Fpk)")
    test_add_line_j_with_z_no_twist(E, str_E=" on E(Fp)")
    test_add_line_j_with_z_no_twist(E_Fpk, str_E=" on E(Fpk)")
    test_vertical_j_denom(E, E_Fpk)

    test_formulas_tate_pairing(E, r, c, u, omega, k)
    test_formula_ate_pairing(E_Fpk, r, ck, u, k)

    print("Tate pairing f_{r, P}(Q)")
    test_bilinear_miller_loop_with_frobenius(miller_function_tate_with_denom, E, E_Fpk, r, t-1, c, ck, u, k)
    print("ate pairing with T = (tr-1)^{} mod r:".format(Tx_exponent))
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, T, c, ck, u, k)
    print("ate pairing with T = tr-1")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, t-1, c, ck, u, k)
    print("ate pairing with f_{u^2-u*p+p^2,Q}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, u**2-u*p+p**2, c, ck, u, k)
    print("ate pairing with f_{2*r,Q}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, 2*r, c, ck, u, k)
    print("ate pairing with f_{-2*r,Q}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, -2*r, c, ck, u, k)
    print("ate pairing with f_{-3*r,Q}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, -3*r, c, ck, u, k)

    print("ate pairing with f_{p,Q}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, p, c, ck, u, k)
    print("ate pairing with f_{p^2,Q}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_function_ate_with_denom, E, E_Fpk, r, p**2, c, ck, u, k)

    print("optimal ate pairing with f_{u^2,Q}(P) * f_{-u, pi(Q)}(P) * l_{[-u]pi(Q), pi^2(Q)}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_ate_bw13_bw19_seq, E, E_Fpk, r, t, c, ck, u, k)
    print("optimal ate pairing with f_{u^2,Q}(P) * f_{-u, pi(Q)}(P) * l_{[-u]pi(Q), pi^2(Q)}(P)")
    test_bilinear_miller_loop_with_frobenius(miller_loop_opt_ate_bw13_bw19, E, E_Fpk, r, t, c, ck, u, k)
    print("\n")

    
if __name__ == "__main__":
    u = ZZ(-2224)
    test_curve(u, k=13)
    u = ZZ(0x1ec5)
    test_curve(u, k=13)
    #u = ZZ(-0x1c21)
    #test_curve(u, k=13)
    #u = ZZ(0x2c90)
    #test_curve(u, k=13)

    u = ZZ(-145)
    test_curve(u, k=19)
    u = ZZ(0x4fa)
    test_curve(u, k=19)

"""
QQty.<t,y> = QQ[]
p = (t^2 + 3*y^2)/4 ; t2 = t^2 - 2*p ; y2 = t*y ; tk_1 = t2 ; tk_2 = t ; p2 = p**2 ; pk_1 = p2
tk_list = [t, t2]
for k in range(3, 14):
    tk = t*tk_1 - p*tk_2
    pk = pk_1*p
    tk_list.append(tk)
    if k in [13,19]:
        print("k = {}".format(k))
        print("tk = {}".format(tk))
        print("tk = {}".format(tk.factor()))
        print("tk^2 - 4*pk = {}".format((tk^2 - 4*pk).factor()))
        print("(p^k+1-tk) = {}".format( (pk+1-tk) ))
    tk_2 = tk_1
    tk_1 = tk
    pk_1 = pk

t13 = tk_list(12)
t19 = tk_list(18)
t13x = t13([tx, yx])
t19x = t19([tx, yx])

"""
