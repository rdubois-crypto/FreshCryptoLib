"""
Build Brezing-Weng and Cocks-Pinch curves on top of BN, BLS12, BLS24 curves
Seeds of inner BN curves obtained from Geppetto paper eprint 2014/976, ethereum parameters, and Pluto/Eris curves
"""
testvector_bn = [
    {'u':-(2**62+2**55+1)},    # Geppetto, eprint 2014/976
    {'u':4965661367192848881}, # ethereum, r-1 has 2-adicity 28 0x44e992b44a6909f1 but no BW6-512, only BW6-515
    {'u':0x36ab400000000000, 'b':11, 'pnbits':253,'rnbits':253, 'label':"u=+2^62-2^59-2^56-2^54-2^52-2^50-2^48+2^46 Hw2NAF=8"},#twist-secure, BW6-509, BW6-510, BW6-512
    {'u':0x381edc0000000000, 'b':15, 'pnbits':253,'rnbits':253, 'label':"u=+2^62-2^59+2^53-2^48-2^45-2^42     Hw2NAF=6"},#twist-secure, BW6-512
    {'u':0x3e59c4000000000000000000, 'b':10, 'pnbits':382,'rnbits':382, 'label':"u=+2^94-2^89+2^87-2^85-2^83+2^81-2^78+2^74 Hw2NAF=8"},#twist-secure, BW6-767
    {'u':0x467eae000000000000000000, 'b': 7, 'pnbits':382,'rnbits':382, 'label':"u=+2^94+2^91-2^89+2^87-2^80-2^78-2^76-2^73 Hw2NAF=8"},#twist-secure, BW6-768
    {'u':-0x4000000000001000008780000000}, # Pluto curve https://github.com/daira/pluto-eris BN-446
]

# with 2^n | (u-1) and n is maximal
testvector_bls24 = [
    {'u': 2**31-2**29+2**22-2**20+1            ,'b':1,'pnbits':305, 'rnbits':245, 'cost_S':158, 'deg_h_S':24, 'label':" 2^31-2^29+2^22-2^20+1"},
    {'u':-2**31-2**28-2**26-2**24-2**20+1      ,'b':1,'pnbits':311, 'rnbits':250, 'cost_S':159, 'deg_h_S':24, 'label':"-2^31-2^28-2^26-2^24-2^20+1"},
    {'u': 2**31+2**29-2**23+2**21-2**18+1      ,'b':1,'pnbits':312, 'rnbits':251, 'cost_S':159, 'deg_h_S':24, 'label':"2^31+2^29-2^23+2^21-2^18+1"},
    {'u':-2**32+2**30+2**22-2**20+1            ,'b':1,'pnbits':315, 'rnbits':253, 'cost_S':160, 'deg_h_S':24, 'label':"-2^32+2^30+2^22-2^20+1"},
    {'u':-2**32+2**30-2**27-2**24-2**20+2**18+1,'b':1,'pnbits':315, 'rnbits':254, 'cost_S':160, 'deg_h_S':24, 'label':"-2^32+2^30-2^27-2^24-2^20+2^18+1"},
    {'u': 2**32-2**29+2**25-2**23+2**21-2**18+1,'b':1,'pnbits':317, 'rnbits':255, 'cost_S':160, 'deg_h_S':24, 'label':"2^32-2^29+2^25-2^23+2^21-2^18+1"},
    {'u':-2**32-2**26-2**23-2**19+1            ,'b':1,'pnbits':319, 'rnbits':257, 'cost_S':161, 'deg_h_S':24, 'label':"-2^32-2^26-2^23-2^19+1"},
]

testvector_bls12 = [
    {'u':   2**63+2**58+2**56+2**51+2**47+2**46+1, 'b': 1,'pnbits': 377, 'rnbits':253, 'cost_S':126, 'deg_h_S':6, 'label':"Zexe eprint 2018/962 Table 16"},
    {'u':  -(2**63+2**62+2**60+2**57+2**48+2**16), 'b': 4,'pnbits': 381, 'rnbits':255, 'cost_S':126, 'deg_h_S':6, 'label':"Zcash https://blog.z.cash/new-snark-curve/"},
    # other curves with high 2-adicity
    {'u':-2**63+2**54+2**51-2**44+1,               'b': 1, 'pnbits':377, 'rnbits':252, 'cost_S':126, 'deg_h_S':6, 'label':"-2^63+2^54+2^51-2^44+1"},
    {'u': 2**63+2**61-2**58-2**56+2**50+1,         'b': 1, 'pnbits':379, 'rnbits':254, 'cost_S':126, 'deg_h_S':6, 'label':"2^63+2^61-2^58-2^56+2^50+1"},
    {'u':-2**64+2**50+2**46-2**42+1,               'b': 1, 'pnbits':383, 'rnbits':256, 'cost_S':126, 'deg_h_S':6, 'label':"-2^64+2^50+2^46-2^42+1"},
    {'u':-2**64+2**51+2**46-2**41+1,               'b': 1, 'pnbits':383, 'rnbits':256, 'cost_S':126, 'deg_h_S':6, 'label':"-2^64+2^51+2^46-2^41+1"},
    {'u':-2**64+2**54-2**50+2**46+1,               'b': 1, 'pnbits':383, 'rnbits':256, 'cost_S':126, 'deg_h_S':6, 'label':"-2^64+2^54-2^50+2^46+1"},
    {'u': 2**64+2**59-2**57-2**55+2**53+2**51+1,   'b': 1, 'pnbits':383, 'rnbits':257, 'cost_S':126, 'deg_h_S':6, 'label':"2^64+2^59-2^57-2^55+2^53+2^51+1"},
]

QQx.<x> = QQ[]
prod_primes = prod(prime_range(10^7))

def compute_beta_lambda(px, rx, tx, yx, D):
    """
    Compute the eigenvalues beta mod px and lambda mod rx such that
    if D!=3 mod 4: beta^2 + D = 0 mod px, lambda^2 + D = 0 mod rx,
    if D=3 mod 4:  beta^2 + beta + (1+D)/4 = 0 mod px, i.e.,
    beta = (-1+sqrt(-D))/2 mod px, and
    lambda^2 + lambda + (1+D)/4 = 0 mod rx, i.e.,
    lambda = (-1+sqrt(-D))/2 mod rx.

    INPUT: polynomial parameters of a pairing-friendly elliptic curve.

    from sage/tnfs/curve/pairing_friendly_curve.py
    at https://gitlab.inria.fr/tnfs-alpha/alpha
    """
    QQx = QQ['x']; #(x,) = QQx._first_ngens(1)
    if D < 0:
        D = -D
    g, u, v = px.xgcd(yx) # u*px + v*yx == g
    inv_yx_modp = QQx(v)/QQ(g)
    g, u, v = rx.xgcd(yx)
    inv_yx_modr = QQx(v)/QQ(g)
    if (D % 4) != 3:
        betax = (tx * inv_yx_modp) % px
        if betax.leading_coefficient() < 0:
            betax = -betax
        assert ((betax**2 + D) % px) == 0
        lambx = ((tx-2)*inv_yx_modr) % rx
        if lambx.leading_coefficient() < 0:
            lambx = -lambx
        assert ((lambx**2 + D) % rx) == 0
    else: #  (-1+sqrt(-D))/2, minimal poly is x^2 + x + (1+D)//4
        betax = ((-1+(tx * inv_yx_modp))/2) % px
        if betax.leading_coefficient() < 0:
            betax = -betax-1
        assert ((betax**2 + betax + (1+D)//4) % px) == 0
        lambx = ((-1+(tx-2)*inv_yx_modr)/2) % rx
        if lambx.leading_coefficient() < 0:
            lambx = -lambx-1
        assert ((lambx**2 + lambx + (1+D)//4) % rx) == 0
    return betax, lambx

def poly_cofactor_twist_g1_g2(k: int, px, rx, tx, cx, yx, D):
    """
    Computes the curve co-factors for the twists

    from sage/tnfs/curve/pairing_friendly_curve.py
    at https://gitlab.inria.fr/tnfs-alpha/alpha
    """
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    twx = px+1+tx
    # does it factor?
    if (D == 3) and (k % 6) == 0:
        # sextic twist
        d = 6
    elif (D == 3) and (k % 3) == 0:
        # cubic twist
        d = 3
    elif ((D == 1) or (D==4)) and (k % 4) == 0:
        # quartic twist
        d = 4
    elif (k % 2) == 0:
        # quadratic twist
        d = 2
    else:
        # no twist
        d = 1
    k1 = k // d
    # compute the curve order up to k1 before computing the d-twist order
    t1 = tx
    # (p+1-t)*(p+1+t) = p^2 + 1 + 2*p - t^2 = p^2 + 1 - (t^2 - 2*p)
    t2 = tx**2 - 2*px
    # tk = t1*t_{k-1} -p*t_{k-2}
    i = 3
    t_im1 = t2 # t_{i-1}
    t_im2 = t1 # t_{i-2}
    while i <= k1:
        t_i = t1*t_im1 - px*t_im2
        t_im2 = t_im1
        t_im1 = t_i
        i = i+1
    if k1 == 1:
        tx_k1 = t1
    elif k1 == 2:
        tx_k1 = t2
    else:
        tx_k1 = t_i
    px_k1 = px**k1
    if d==3 or d==6 or d==4:
        yx_k1_square = (tx_k1**2 - 4*px_k1)/(-D)
        lc = yx_k1_square.leading_coefficient()
        assert lc.is_square()
        yx_k1_square_monic = yx_k1_square / lc
        yx_k1_factors = yx_k1_square_monic.factor()
        yx_k1 = lc.sqrt()
        for fa, ee in yx_k1_factors:
            assert (ee % 2) == 0
            yx_k1 = yx_k1 * fa**(ee//2)
        assert yx_k1**2 == yx_k1_square
    else:
        yx_k1 = 1
    if d==3 or d==6:
        if d==6:
            E2_order = px_k1+1-(-3*yx_k1+tx_k1)/2
            E2_order_= px_k1+1-( 3*yx_k1+tx_k1)/2
            g2twx = px_k1+1+(-3*yx_k1+tx_k1)/2
            g2twx_= px_k1+1+( 3*yx_k1+tx_k1)/2
        elif d==3:
            E2_order = px_k1+1-(-3*yx_k1-tx_k1)/2
            E2_order_= px_k1+1-( 3*yx_k1-tx_k1)/2
            g2twx = px_k1+1+(-3*yx_k1-tx_k1)/2
            g2twx_= px_k1+1+( 3*yx_k1-tx_k1)/2
        if (E2_order % rx) != 0 and (E2_order_ % rx) == 0:
            E2_order = E2_order_
            g2twx = g2twx_
    elif d==4:
        if D==1:
            E2_order = px_k1 + 1 + yx_k1
            g2twx = px_k1 + 1 - yx_k1 # quadratic twist of G2
        elif D==4:
            E2_order = px_k1 + 1 + 2*yx_k1
            g2twx = px_k1 + 1 - 2*yx_k1 # quadratic twist of G2
        if (E2_order % rx) != 0 and (g2twx % rx) == 0:
            E2_order, g2twx = g2twx, E2_order
    elif d == 2:
        E2_order = px_k1 + 1 + tx_k1
        g2twx = px_k1 + 1 - tx_k1 # quadratic twist of G2
    else: # d==1
        assert d==1
        E2_order = px_k1 + 1 - tx_k1
        g2twx = px_k1 + 1 + tx_k1 # quadratic twist of G2

    assert (E2_order % rx) == 0
    g2cx = E2_order // rx # irreducible
    # do cx, twx, g2cx, g2twx factor?
    polys_cofact_twists = [cx, twx, g2cx, g2twx]
    label_factors = ["cx", "twx", "g2cx", "g2twx"]
    small_cofactors = [1, 1, 1, 1]
    return twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors

def compute_curve_coefficient_D1D3(p, t, r, D):
    if D != 4 and D != 1 and D != 3:
        return 0
    rc = p+1-t
    Fp = GF(p, proof=False)
    ab = 1
    ok = False
    while not ok:
        if D == 3:
            E = EllipticCurve([Fp(0), Fp(ab)])
        elif D == 1 or D ==4:
            E = EllipticCurve([Fp(ab), Fp(0)])
        ok=True
        i = 0
        while ok and i < 4:
            P = E.random_element()
            ok = rc*P == E(0)
            i = i+1
        if not ok:
            if ab > 0:
                ab = -ab
            else:
                ab = -ab + 1
    return ab

def get_int_coeffs_denom_from_poly(px):
    c = px.coefficients(sparse=False)
    #c = px.list()
    d = lcm([ci.denom() for ci in c])
    int_coeffs = [int(ci*d) for ci in c]
    return int_coeffs, d

def write_test_vector(outfile, magma_outfile, inner_n, k, u, D, p, t, r, px, rx, cx, yx, tx, ht, hy, inner_curve_label="bls", with_prec_comma:bool=True):
    # compute cofactors
    twx, g2cx, g2twx, polys_cofact_twists, label_factors, small_cofactors = poly_cofactor_twist_g1_g2(k, px, rx, tx, cx, yx, D)
    # compute betax and lambx
    betax, lambx = compute_beta_lambda(px, rx, tx, yx, D)
    
    # also compute a curve coefficient
    # write:
    #    {'u': u, 'b': b, 'px': [...], ...}
    ab = compute_curve_coefficient_D1D3(p, t, r, D)
    if with_prec_comma:
        s = ",\n"
        m = ",\n"
    else:
        s = ""
        m = ""
    if hy == 0:# write ht in hexa if hy is 0
        s = s+"    {{'u':{:#x}, 'D':{}, 'ht':{:#x}, 'hy':{}, ".format(u, D, ht, hy)
    else:
        s = s+"    {{'u':{:#x}, 'D':{}, 'ht':{}, 'hy':{}, ".format(u, D, ht, hy)
    if D == 3:
        s = s+"'b':{}, ".format(ab)
    elif D==1 or D==4:
        s = s+"'a':{}, ".format(ab)
    s = s + "'pnbits':{}, 'rnbits':{}, ".format(p.nbits(), r.nbits())

    m = m + "    <{},{},{},{},{},{},{},".format(u, D, ht, hy, ab, p.nbits(), r.nbits())
    # compute denominator of polynomial coefficients, and write 1) poly with int coeffs, 2) denom.
    for poly, str_poly in [(px, "px"), (rx, "rx"), (cx, "cx"), (yx, "yx"), (tx, "tx"), (betax, "betax"), (lambx, "lambx"), (g2cx, "g2cx")]:
        c, d = get_int_coeffs_denom_from_poly(poly)
        s = s + "'{0}':[{1}], '{0}_denom':{2}, ".format(str_poly, ",".join([str(ci) for ci in c]), d)
        m +="[{}],{},".format(",".join([str(ci) for ci in c]), d)
    s = s + "'label':\"{}{}_BW{}D{}_{}\"}}".format(inner_curve_label, inner_n, k, D, p.nbits())
    m += "\"{}{}_BW{}D{}_{}\">".format(inner_curve_label, inner_n, k, D, p.nbits())
    outfile.write(s)
    magma_outfile.write(m)

def write_test_vector_cpk(outfile, magma_outfile, inner_n, k, u, D, p, t0, y0, r, rx, tr, y, ht, sign_ht, hy, sign_y0, miller_scalars, miller_scalars_tate, inner_curve_label="bls", with_prec_comma:bool=True):
    # compute cofactors: not implemented for Cocks-Pinch curves
    # compute betax and lambx: not implemented for Cocks-Pinch curves
    
    # also compute a curve coefficient
    # write:
    #    {'u': u, 'b': b, 'rx': [...], ...}
    ab = compute_curve_coefficient_D1D3(p, tr, r, D)
    if with_prec_comma:
        s = ",\n"
        m = ",\n"
    else:
        s = ""
        m = ""
    s = s+"    {{'k':{}, 'u':{:#x}, 'D':{}, 'ht':{}, 'hy':{}, 'sign_ht':{}1, 'sign_y0':{}1, 'tr0':{}, 'y0':{}, 'tr':{}, 'y':{}, 'p':{}, 'r':{}, 'miller':{}, ".format(k, u, D, ht, hy, sign_ht, sign_y0, t0, y0, tr, y, p, r, miller_scalars)
    if miller_scalars_tate is not None:
        s = s + "'miller_tate':{}, ".format(miller_scalars_tate)
    if D == 3:
        s = s+"'b':{}, ".format(ab)
    elif D==1 or D==4:
        s = s+"'a':{}, ".format(ab)
    s = s + "'pnbits':{}, 'rnbits':{}, ".format(p.nbits(), r.nbits())

    m = m + "    <{},{},{},{},{},{}1,{}1,{},{},{},{},{},{},{},{},{},".format(k, u, D, ht, hy, sign_ht, sign_y0, t0, y0, tr, y, p, r, ab, p.nbits(), r.nbits())
    # compute denominator of polynomial coefficients, and write 1) poly with int coeffs, 2) denom.
    for poly, str_poly in [(rx, "rx")]:
        c, d = get_int_coeffs_denom_from_poly(poly)
        s = s + "'{0}':[{1}], '{0}_denom':{2}, ".format(str_poly, ",".join([str(ci) for ci in c]), d)
        m +="[{}],{},".format(",".join([str(ci) for ci in c]), d)
    s = s + "'label':\"{}{}_CP{}D{}_{}\"}}".format(inner_curve_label, inner_n, k, D, p.nbits())
    m += "\"{}{}_CP{}D{}_{}\">".format(inner_curve_label, inner_n, k, D, p.nbits())
    outfile.write(s)
    magma_outfile.write(m)
    
def compute_min_max_ht_hy(R_BW6, Tr0, Y0, max_pbits, min_pbits=None):
    # to find the max values for ht hy:
    QQhthy.<h_t,h_y> = QQ[]
    Q0_BW6_ht_hy = ((Tr0(u) + h_t*R_BW6(u))^2 + 3*(Y0(u) + h_y*R_BW6(u))^2)/4
    Q0_BW6_ht = Q0_BW6_ht_hy([x, 0])
    Q0_BW6_hy = Q0_BW6_ht_hy([0, x])
    RRa = RealField(ceil(max_pbits*log(2) + 2))
    roots_Q0_ht = (Q0_BW6_ht - (2**max_pbits-1)).roots(RRa)
    roots_Q0_hy = (Q0_BW6_hy - (2**max_pbits-1)).roots(RRa)
    ht_pos = [floor(ri) for (ri, ei) in roots_Q0_ht if ri > 0]
    hy_pos = [floor(ri) for (ri, ei) in roots_Q0_hy if ri > 0]
    ht_neg = [ceil(ri) for (ri, ei) in roots_Q0_ht if ri < 0]
    hy_neg = [ceil(ri) for (ri, ei) in roots_Q0_hy if ri < 0]
    max_ht_pos = max(ht_pos); min_ht_neg = min(ht_neg); max_hy_pos = max(hy_pos); min_hy_neg = min(hy_neg)
    if min_pbits is None:
        return max_ht_pos, min_ht_neg, max_hy_pos, min_hy_neg
    # if we want e.g. BW6-672 with BLS24-315: we need to increase a lot the sizes of one of ht, hy for p to reach the target size
    # compute max values of ht, hy that reach 2^(min_pbits-1)
    roots_Q0_ht = (Q0_BW6_ht - (2**(max_pbits-1))).roots(RRa)
    roots_Q0_hy = (Q0_BW6_hy - (2**(max_pbits-1))).roots(RRa)
    ht_pos_ = [ceil(ri) for (ri, ei) in roots_Q0_ht if ri > 0]
    hy_pos_ = [ceil(ri) for (ri, ei) in roots_Q0_hy if ri > 0]
    ht_neg_ = [floor(ri) for (ri, ei) in roots_Q0_ht if ri < 0]
    hy_neg_ = [floor(ri) for (ri, ei) in roots_Q0_hy if ri < 0]
    min_ht_pos = min(ht_pos_); max_ht_neg = max(ht_neg_); min_hy_pos = min(hy_pos_); min_hy_neg = min(hy_neg)
    return max_ht_pos, min_ht_neg, max_hy_pos, min_hy_neg, min_ht_pos, max_ht_neg, min_hy_pos, min_hy_neg

def get_small_factors(n, prod_pr):
    gcd_n = gcd(n, prod_pr)
    cn = 1
    N = n
    while gcd_n > 1:
        cn *= gcd_n
        N = N // gcd_n
        gcd_n = gcd(N, gcd_n)
    return cn, N

def compute_params_BN_BLSn_BWk(n, k, D, u, inner:str="bls", max_pbits:int=None, min_pbits:int=None, min_ht:int=0, max_ht:int=None, min_hy:int=0, max_hy:int=None, outfilename:str=None, max_curves:int=None, smallest_hy:int=False, hy_0:int=False):
    print("inner {} n = {} k = {} D = {} u = {:#x}, max_pbits = {}, min_pbits = {}, min_ht = {}, max_ht = {}, min_hy = {}, max_hy = {}, max_curves = {}, smallest_hy = {}, hy_0 = {}".format(inner, n, k, D, u, max_pbits, min_pbits, min_ht, max_ht, min_hy, max_hy, max_curves, smallest_hy, hy_0))
    write = outfilename is not None
    if n == 24 and inner == "bls":
        R_inner = cyclotomic_polynomial(24)
        P_inner = (x-1)**2*R_inner/3 + x # (x^10 - 2*x^9 + x^8 - x^6 + 2*x^5 - x^4 + x^2 + x + 1)/3
        Y_inner = (x-1)*(2*x**4 - 1)/3
        C_inner = (x-1)**2/3
        T_inner = x+1
    elif n == 12 and inner == "bls":
        R_inner = cyclotomic_polynomial(12)
        P_inner = (x-1)**2*R_inner/3 + x # (x^6 - 2*x^5 + 2*x^3 + x + 1)/3
        Y_inner = (x-1)*(2*x**2 - 1)/3
        C_inner = (x-1)**2/3
        T_inner = x+1
    elif n == 12 and inner == "bn":
        R_inner = 36*x**4 + 36*x**3 + 18*x**2 + 6*x + 1
        T_inner = 6*x**2+1
        P_inner = 36*x**4 + 36*x**3 + 24*x**2 + 6*x + 1
        Y_inner = 6*x**2 + 4*x + 1
        C_inner = QQx(1)

    R_BW6 = P_inner
    u = ZZ(u)
    p = ZZ(P_inner(u))
    r = ZZ(R_inner(u))
    y = ZZ(Y_inner(u))
    c = ZZ(C_inner(u))
    t = ZZ(T_inner(u))
    if write:
        outfile = open(outfilename + ".py", "w")
        magma_outfile = open(outfilename + ".mag", "w")
        if min_pbits is not None:
            tag_min_pbits = "_{}".format(min_pbits)
        else:
            tag_min_pbits = ""
        if max_pbits is not None:
            tag_max_pbits = "_{}".format(max_pbits)
        else:
            tag_max_pbits = ""
        outfile.write("testvector_{}{}_{}_bw{}{}{} = [\n".format(inner, n, p.nbits(), k, tag_min_pbits, tag_max_pbits))
        magma_outfile.write("testvector_{}{}_{}_bw{}{}{} := [\n".format(inner, n, p.nbits(), k, tag_min_pbits, tag_max_pbits))
        write_first_item = True # because the last item should not have a comma after it, but we don't know in advance if an item is the last item
            
    print("u={:#x} {} bits".format(u, u.nbits()))
    print("p={:#x} {} bits".format(p, p.nbits()))
    print("r={:#x} {} bits".format(r, r.nbits()))
    print("y={:#x} {} bits p = (t^2+3*y^2)/4".format(y, y.nbits()))
    print("c={:#x} {} bits p+1-t = r*c".format(c, c.nbits()))
    cr1, R1 = get_small_factors(r-1, prod_primes)
    cp1, P1 = get_small_factors(p-1, prod_primes)
    print("r-1 = {} * R1 # R1 {} bits".format(cr1.factor(), R1.nbits()))
    print("p-1 = {} * P1 # P1 {} bits".format(cp1.factor(), P1.nbits()))
    # for simplified SWU (Wahby--Boneh hash-on-curve: needs an isogeny)
    cy, Y = get_small_factors(y, prod_primes)
    print("y = {} * Y # Y {} bits".format(cy.factor(), Y.nbits()))
    cc, C = get_small_factors(c, prod_primes)
    print("c = {} * C # C {} bits".format(cc.factor(), C.nbits()))
    print("gcd(y, c) = {}".format((gcd(cc, cy)).factor()))
    tw = (p + 1 + t)
    ctw, W = get_small_factors(tw, prod_primes)
    print("tw = {} * W # twist order".format(ctw.factor()))
    print("gcd(y, tw) = {}".format((gcd(ctw, cy)).factor()))

    # now find a k-th root of unity modulo R_BW6
    # the following assumes that k=6
    K.<a> = NumberField(R_BW6)
    Phi_k = cyclotomic_polynomial(k)
    rr = Phi_k.roots(K)
    if len(rr) == 0:
        return
    # yes there are k-th roots of unity mod R_CPk
    zeta_k = [z_i.polynomial() for (z_i, e_i) in rr]
    if len(zeta_k) == 2:
        Tr0 = zeta_k[0]+1
        Tr1 = zeta_k[1]+1
    print("traces are {} and {}".format(Tr0, Tr1))
    if (k == 6 or k == 3) and Tr1.constant_coefficient() == 3:
        Tr0, Tr1 = Tr1, Tr0
    if k == 4 and Tr0.leading_coefficient() < 0:
        Tr0, Tr1 = Tr1, Tr0
        
    #zeta_6 = ((1 + sqrt(K(-3)))/2)
    #Tr0 = zeta_6.polynomial() + 1
    #Tr1 = (zeta_6^5).polynomial() + 1
    print("zeta_k = {}".format(zeta_k))
    
    print("R_BW{}   = {}".format(R_BW6, k))
    print("Tr0_BW{} = {}".format(Tr0, k))
    print("Tr1_BW{} = {}".format(Tr1, k))
    inv_sqrt_D = (K(-1/D)).sqrt().polynomial()
    assert (D*(inv_sqrt_D^2 % R_BW6) == -1)
    Y0 = ((Tr0-2)*inv_sqrt_D) % R_BW6
    Y1 = ((Tr1-2)*inv_sqrt_D) % R_BW6
    print("1/sqrt(-{}) = ({})/{}".format(D, D*inv_sqrt_D, D))
    print("Y0_BW{} = ({})/{}".format(k, Y0*D, D))
    print("Y1_BW{} = ({})/{}".format(k, Y1*D, D))
    # convention: positive leading term
    if Y0.leading_coefficient() < 0:
        Y0 = -Y0
    if Y1.leading_coefficient() < 0:
        Y1 = -Y1
        print("Y0 = ({})/{}".format(Y0*D, D))
        print("Y1 = ({})/{}".format(Y1*D, D))

    if min_pbits is None:
        max_ht_pos, max_ht_neg, max_hy_pos, max_hy_neg = compute_min_max_ht_hy(R_BW6, Tr0, Y0, max_pbits)
        min_ht_pos = 0; min_ht_neg = 0; min_hy_pos = 0; min_hy_neg = 0 # start at 0 by default
    else:
        max_ht_pos, max_ht_neg, max_hy_pos, max_hy_neg, min_ht_pos, min_ht_neg, min_hy_pos, min_hy_neg = compute_min_max_ht_hy(R_BW6, Tr0, Y0, max_pbits, min_pbits)
    if max_ht is not None:
        max_ht_pos = min(abs(max_ht), max_ht_pos)
        max_ht_neg = max(-abs(max_ht), max_ht_neg)
    if max_hy is not None:
        max_hy_pos = min(abs(max_hy), max_hy_pos)
        max_hy_neg = max(-abs(max_hy), max_hy_neg)
    min_Q_bits = max_pbits
    print("ht in [{}, {}] hy in [{}, {}]".format(max_ht_neg, max_ht_pos, max_hy_neg, max_hy_pos))
    if min_pbits is None:
        ht = min_ht # 0 by default
    else:
        ht = min(min_ht_pos, abs(min_ht_neg))
        min_ht = abs(min_ht)
        if min_ht >= ht and min_ht <= max_ht_pos:
            ht = min_ht
    # 3 does NOT divide ht, and 2 divides (ht+hy) and (ht - hy)
    ht = ht - (ht % 3) + 1
    print("starting at ht = {}".format(ht))
    curves = 0
    while (max_curves is None or curves < max_curves) and ((ht >= 0 and ht <= max_ht_pos) or (ht < 0 and ht >= max_ht_neg)):
        hy = 0
        if (ht % 2) == 1:
            hy = 1
        if max_ht_pos * max_hy_pos > 10**5:
            hy_max_here = min(max(ht+12, ceil(ht*1.1)), max_hy_pos) # why????
        else:
            hy_max_here = max_hy_pos
        found_hy = False
        while (not hy_0 or hy == 0) and (not smallest_hy or not found_hy) and (max_curves is None or curves < max_curves) and (hy >= 0 and hy <= hy_max_here):
            Q0_BW6 = ((Tr0 + ht*R_BW6)^2 + D*(Y0 + hy*R_BW6)^2)/4
            Q0_BW6_ = ((Tr0 + ht*R_BW6)^2 + D*(-Y0 + hy*R_BW6)^2)/4
            Q1_BW6 = ((Tr1 + ht*R_BW6)^2 + D*(Y1 + hy*R_BW6)^2)/4
            Q1_BW6_ = ((Tr1 + ht*R_BW6)^2 + D*(-Y1 + hy*R_BW6)^2)/4
            if hy == 0:
                L = [(Q0_BW6, Tr0, Y0, 1), (Q1_BW6, Tr1, Y1, 1)]
            else:
                L = [(Q0_BW6, Tr0, Y0, 1), (Q0_BW6_, Tr0, -Y0, -1), (Q1_BW6, Tr1, Y1, 1), (Q1_BW6_, Tr1, -Y1, -1)]
            for (Qa, Tra, Ya, sign_hy) in L:
                if Qa.is_irreducible():
                    q0_BW6 = Qa(u)
                    if q0_BW6 in ZZ:
                        q0_BW6 = ZZ(q0_BW6)
                        if q0_BW6.nbits() <= max_pbits and gcd(q0_BW6, prod_primes) == 1:
                            if q0_BW6.is_prime():
                                curves += 1
                                found_hy = True
                                Ya_hy = Ya + hy*R_BW6
                                Tra_ht = Tra + ht*R_BW6
                                C_BW6 = (Qa + 1 - Tra_ht) // R_BW6
                                min_Q_bits = min(min_Q_bits, q0_BW6.nbits())
                                print("curve #{}".format(curves))
                                if hy_0:
                                    print("ht={}={:#x},hy={} {} bits q={:#x}".format(ht,ht,sign_hy*hy,q0_BW6.nbits(),q0_BW6))
                                else:
                                    print("ht={:3d},hy={:3d} {} bits q={:#x}".format(ht,sign_hy*hy,q0_BW6.nbits(),q0_BW6))
                                print("Q = ({})/9".format(Qa*9))
                                print("T = ({})".format(Tra_ht))
                                print("Y = ({})/3".format(Ya_hy*3))
                                print("C = ({})/9".format(C_BW6*9))
                                # print cofactors, small factors of y, of cofactors, of twist order.
                                y0 = ZZ(Ya_hy(u))
                                tr0 = ZZ(Tra_ht(u))
                                gcd_y = gcd(y0, prod_primes)
                                cy, Y = get_small_factors(y0, prod_primes)
                                print("    y = {} * Y # Y {} bits".format(cy.factor(), Y.nbits()))
                                # and also the curve cofactor, or the twist order, should have a small cofactor, and finally there should have a common small factor so that an isogeny is defined.
                                curve_order = (q0_BW6 + 1 - tr0)
                                assert (curve_order % p) == 0
                                curve_cofactor = curve_order // p
                                cc, C = get_small_factors(curve_cofactor, prod_primes)
                                twist_order = (q0_BW6 + 1 + tr0)
                                ctw, W = get_small_factors(twist_order, prod_primes)
                                print("    c = {} * C # curve cofactor".format(cc.factor()))
                                print("    gcd(y, c) = {}".format((gcd(cc, cy)).factor()))
                                print("    tw = {} * W # twist order".format(ctw.factor()))
                                print("    gcd(y, tw) = {}".format((gcd(ctw, cy)).factor()))
                                if write:
                                    assert C_BW6 * R_BW6 == (Qa + 1 - Tra_ht)
                                    assert Tra_ht**2 + D*Ya_hy**2 == 4*Qa
                                    if write_first_item:
                                        write_test_vector(outfile, magma_outfile, n, k, u, D, q0_BW6, tr0, p, Qa, R_BW6, C_BW6, Ya_hy, Tra_ht, ht, sign_hy*hy, inner_curve_label=inner, with_prec_comma=False)
                                        write_first_item = False
                                    else:
                                        write_test_vector(outfile, magma_outfile, n, k, u, D, q0_BW6, tr0, p, Qa, R_BW6, C_BW6, Ya_hy, Tra_ht, ht, sign_hy*hy, inner_curve_label=inner, with_prec_comma=True)
            # the sign is applied to y0 and y1, and hy > 0
            hy = hy+2 # because we need 2 | (ht + hy)
        # ht can be positive or negative
        if ht <= 0:
            ht = -ht+1
            if (ht % 3) == 0:
                ht = ht + 1
        else:
            ht = -ht # we should not have ht = 0 mod 3
            if (ht % 3) == 0:
                ht = ht - 1
    if write:
        outfile.write("\n]\n")
        outfile.close()
        magma_outfile.write("];\n")
        magma_outfile.close()

def compute_params_BLSn_CPk(n, k, D, u, max_pbits:int=None, min_ht:int=0, max_ht:int=None, min_hy:int=0, max_hy:int=None, outfilename:str=None):
    write = outfilename is not None
    if n == 24:
        R_BLS = cyclotomic_polynomial(24)
        P_BLS = (x-1)**2*R_BLS/3 + x # (x^10 - 2*x^9 + x^8 - x^6 + 2*x^5 - x^4 + x^2 + x + 1)/3
        Y_BLS = (x-1)*(2*x**4 - 1)/3
        C_BLS = (x-1)**2/3
        T_BLS = x+1
    elif n == 12:
        R_BLS = cyclotomic_polynomial(12)
        P_BLS = (x-1)**2*R_BLS/3 + x # (x^6 - 2*x^5 + 2*x^3 + x + 1)/3
        Y_BLS = (x-1)*(2*x**2 - 1)/3
        C_BLS = (x-1)**2/3
        T_BLS = x+1
    R_CPk = P_BLS
    u = ZZ(u)
    p = ZZ(P_BLS(u))
    r = ZZ(R_BLS(u))
    y = ZZ(Y_BLS(u))
    c = ZZ(C_BLS(u))
    t = ZZ(T_BLS(u))
    if write:
        outfile = open(outfilename + ".py", "w")
        magma_outfile = open(outfilename + ".mag", "w")
        if max_pbits is not None:
            outfile.write("testvector_bls{}_{}_cp{}d{}_{} = [\n".format(n, p.nbits(), k, D, max_pbits))
            magma_outfile.write("testvector_bls{}_{}_cp{}d{}_{} := [\n".format(n, p.nbits(), k, D, max_pbits))
        else:
            outfile.write("test_vector_bls{}_{}_cp{}d{} = [\n".format(n, p.nbits(), k, D))
            magma_outfile.write("test_vector_bls{}_{}_cp{}d{} := [\n".format(n, p.nbits(), k, D))
        write_first_item = True # because the last item should not have a comma after it, but we don't know in advance if an item is the last item
            
    print("u={:#x} {} bits".format(u, u.nbits()))
    print("p={:#x} {} bits".format(p, p.nbits()))
    print("r={:#x} {} bits".format(r, r.nbits()))
    print("y={:#x} {} bits p = (t^2+3*y^2)/4".format(y, y.nbits()))
    print("c={:#x} {} bits p+1-t = r*c".format(c, c.nbits()))
    cr1, R1 = get_small_factors(r-1, prod_primes)
    cp1, P1 = get_small_factors(p-1, prod_primes)
    print("r-1 = {} * R1 # R1 {} bits".format(cr1.factor(), R1.nbits()))
    print("p-1 = {} * P1 # P1 {} bits".format(cp1.factor(), P1.nbits()))
    # for simplified SWU (Wahby--Boneh hash-on-curve: needs an isogeny)
    cy, Y = get_small_factors(y, prod_primes)
    print("y = {} * Y # Y {} bits".format(cy.factor(), Y.nbits()))
    cc, C = get_small_factors(c, prod_primes)
    print("c = {} * C # C {} bits".format(cc.factor(), C.nbits()))
    print("gcd(y, c) = {}".format((gcd(cc, cy)).factor()))
    tw = (p + 1 + t)
    ctw, W = get_small_factors(tw, prod_primes)
    print("tw = {} * W # twist order".format(ctw.factor()))
    print("gcd(y, tw) = {}".format((gcd(ctw, cy)).factor()))

    # now find a k-th root of unity modulo R_CPk if possible,
    # otherwise modulo r_cpk (as a prime number)
    r_cpk = p
    Fp = GF(p)
    K.<a> = NumberField(R_CPk)
    Phi_k = cyclotomic_polynomial(k)
    rr = Phi_k.roots(K)
    alg_zeta_k = len(rr) > 0
    mod_zeta_k = False
    if len(rr) > 0: # yes there are k-th roots of unity mod R_CPk
        zeta_k = [z_i.polynomial() for (z_i, e_i) in rr]
        print("algebraic {}-roots of unity found mod R_CP{}".format(k,k))
    elif ((p-1) % k) == 0:
        rr = Phi_k.roots(Fp)
        zeta_k = [ZZ(z_i) for (z_i, e_i) in rr]
        for i in range(len(zeta_k)):
            if abs(zeta_k[i]-p) < abs(zeta_k[i]):
                zeta_k[i] = zeta_k[i]-p
            assert (Phi_k(zeta_k[i]) % p) == 0
        mod_zeta_k = True
        print("{} {}-roots of unity found mod r_cp{}".format(len(zeta_k), k,k))
    if not alg_zeta_k and not mod_zeta_k:
        return
    
    # now find inv_sqrt
    alg_sqrt_D = False
    mod_sqrt_D = False
    if (K(-D)).is_square():
        alg_sqrt_D = True
        inv_sqrt_D = (K(-1/D)).sqrt().polynomial()
    elif Fp(-D).is_square():
        inv_sqrt_D = ZZ((Fp(-1/D)).sqrt())
        mod_sqrt_D = True
        if p-inv_sqrt_D < inv_sqrt_D:
            inv_sqrt_D = p-inv_sqrt_D
    else: # no sqrt(-D) mod p, impossible to continue
        print("no sqrt(-{}) mod p".format(D))
        return

    min_q_bits = max_pbits

    # now, Cocks-Pinch
    # 1. compute the set of zeta_k_i values and tr_i = zeta_k_i + 1
    # 2. compute the set of y_i
    # 3. iterate over ht, hy
    #    3.1 compute tr_i_ht
    #    3.2 compute y_i_hy
    list_tr_y = [None]*len(zeta_k)
    i = 0
    for zeta_k_i in zeta_k:
        if alg_zeta_k: # it means zeta_k_i is a polynomial
            tr_i = zeta_k_i(u)+1
        else: # zeta_k_i is an interger (lifted mod p)
            tr_i = zeta_k_i + 1
        # normalize
        if abs(tr_i+p) < abs(tr_i):
            tr_i = tr_i+p
        elif abs(tr_i-p) < abs(tr_i):
            tr_i = tr_i-p
        if alg_zeta_k and alg_sqrt_D:
            y_i = ((zeta_k_i-1)*inv_sqrt_D) % R_CPk # (t-2)/sqrt(-D)
            y_i = ZZ(y_i(u))
        elif alg_sqrt_D:
            y_i = ((tr_i-2)*ZZ(inv_sqrt_D(u))) % p # (t-2)/sqrt(-D)
        else:
            y_i = ((tr_i-2)*inv_sqrt_D) % p # (t-2)/sqrt(-D)
        # normalize
        if abs(y_i+p) < abs(y_i):
            y_i = y_i+p
        elif abs(y_i-p) < abs(y_i):
            y_i = y_i-p
        # estimate max_ht, max_hy
        # (tr_i + ht*p)^2 <= q <= 2^max_pbits
        # tr_i + ht*p <= 2^(max_pbits/2)
        # ht*p <= (2^(max_pbits/2) - tr_i)
        # ht <= (2^(max_pbits/2) - tr_i)/p
        max_ht = ceil(((sqrt(RR(0x1 << (max_pbits)))) - tr_i)/p)
        # (-y_i + D*hy*p)^2 <= q <= 2^max_pbits
        # -y_i + D*hy*p <= 2^(max_pbits/2)
        # hy*D*p <= (2^(max_pbits/2) + y_i)
        # hy <= (2^(max_pbits/2) + y_i)/(p*D)
        max_hy = ceil((RR(sqrt(0x1 << (max_pbits))) + abs(y_i))/(p*D))
        list_tr_y[i] = (tr_i, y_i, max_ht, max_hy)
        i += 1
    print("{} combinations of tr_0 and y_0".format(len(list_tr_y)))
    max_ht = max([li[2] for li in list_tr_y])
    max_hy = max([li[3] for li in list_tr_y])
    print("max_ht = {}".format([li[2] for li in list_tr_y]))
    print("max_hy = {}".format([li[3] for li in list_tr_y]))
    ht = min_ht
    while ht <= max_ht:
        hy = min_hy
        while hy <= max_hy:
            i = 0
            for tr_i, y_i, _, _ in list_tr_y:
                if ht == 0:
                    list_tr_i_ht_sign_ht = [(tr_i,'+')]
                else:
                    list_tr_i_ht_sign_ht = [(tr_i+ht*p,'+'), (tr_i-ht*p, '-')]
                for tr_i_ht, sign_ht in list_tr_i_ht_sign_ht:
                    if hy == 0:
                        list_y_i_hy_sign_y0 = [(y_i, ' ')]
                    else:
                        list_y_i_hy_sign_y0 = [(y_i+hy*p, ' '), (-y_i+hy*p, '-')]
                    for (y_i_hy, sign_y0) in list_y_i_hy_sign_y0:
                        p_cpk_4 = (tr_i_ht**2 + D*y_i_hy**2)
                        if (p_cpk_4 % 4) != 0:
                            continue
                        p_cpk = p_cpk_4 // 4
                        # now test for primality
                        if p_cpk.nbits() <= max_pbits and gcd(p_cpk, prod_primes) == 1:
                            if p_cpk.is_prime():
                                min_q_bits = min(min_q_bits, p_cpk.nbits())
                                print("ht={:3d}, tr=tr0{}ht*p, hy={:3d}, y={}y0+hy*p, {} bits q={:#x}".format(ht,sign_ht,hy,sign_y0,p_cpk.nbits(),p_cpk))
                                curve_order = p_cpk + 1 - tr_i_ht
                                assert (curve_order % p) == 0
                                curve_cofactor = curve_order // p
                                cc, C = get_small_factors(curve_cofactor, prod_primes)
                                twist_order = (p_cpk + 1 + tr_i_ht)
                                ctw, W = get_small_factors(twist_order, prod_primes)
                                print("    c = {} * C # curve cofactor".format(cc.factor()))
                                print("    gcd(y, c) = {}".format((gcd(cc, cy)).factor()))
                                print("    tw = {} * W # twist order".format(ctw.factor()))
                                print("    gcd(y, tw) = {}".format((gcd(ctw, cy)).factor()))
                                # compute shortest vectors for optimal pairing
                                if k == 8 or k == 12:
                                    assert ((p_cpk + 1 - tr_i) % p) == 0
                                    if k == 8:
                                        assert (((tr_i-1)**4 + 1) % p) == 0
                                    elif k == 12:
                                        assert (((tr_i-1)**4 - (tr_i-1)**2 + 1) % p) == 0
                                    M = Matrix(ZZ, 4, 4, [p, 0,0,0,  -(tr_i-1),1,0,0,  0,-(tr_i-1),1,0, 0,0,-(tr_i-1),1])
                                    R = M.LLL()
                                    miller_smallest_i = -1
                                    miller_smallest_bits = ceil(p.nbits()/4)
                                    expected_miller_smallest_bits = miller_smallest_bits
                                    for i_ in range(4):
                                        miller_bits = max([abs(R[i_][j]).nbits() for j in range(4)])
                                        assert (sum([R[i_][j]*p_cpk**j for j in range(4)]) % p) == 0
                                        if miller_bits < miller_smallest_bits:
                                            miller_smallest_bits = miller_bits
                                            miller_smallest_i = i_
                                    print("optimal ate pairing computation: expected miller scalar length is {}, obtained {}".format(expected_miller_smallest_bits, miller_smallest_bits))
                                    miller_scalars_ate = [R[miller_smallest_i][j] for j in range(4)]
                                else:
                                    miller_scalars_ate = None
                                # compute shortest scalars for optimal Tate pairing
                                if k == 8: # eigenvalue is q^2 mod r
                                    M = Matrix(ZZ, 2, 2, [p, 0, (-(tr_i-1)**2) % p,1])
                                    R = M.LLL()
                                    miller_smallest_i = -1
                                    miller_smallest_bits = ceil(p.nbits()/2)
                                    expected_miller_smallest_bits = miller_smallest_bits
                                    p_cpk2 = p_cpk**2
                                    for i_ in range(2):
                                        miller_bits = max([abs(R[i_][j]).nbits() for j in range(2)])
                                        assert (sum([R[i_][j]*p_cpk2**j for j in range(2)]) % p) == 0
                                        if miller_bits < miller_smallest_bits:
                                            miller_smallest_bits = miller_bits
                                            miller_smallest_i = i_
                                    print("optimal Tate pairing computation: expected miller scalar length is {}, obtained {}".format(expected_miller_smallest_bits, miller_smallest_bits))
                                    miller_scalars_tate = [R[miller_smallest_i][j] for j in range(2)]
                                else:
                                    miller_scalars_tate = None
                                if write:
                                    if write_first_item:
                                        write_test_vector_cpk(outfile, magma_outfile, n, k, u, D, p_cpk, tr_i, y_i, p, P_BLS, tr_i_ht, y_i_hy, ht, sign_ht, hy, sign_y0, miller_scalars_ate, miller_scalars_tate, with_prec_comma=False)
                                        write_first_item = False
                                    else:
                                        write_test_vector_cpk(outfile, magma_outfile, n, k, u, D, p_cpk, tr_i, y_i, p, P_BLS, tr_i_ht, y_i_hy, ht, sign_ht, hy, sign_y0, miller_scalars_ate, miller_scalars_tate, with_prec_comma=True)
                                
                i += 1
            hy += 1
        ht += 1
    if write:
        outfile.write("\n]\n")
        outfile.close()
        magma_outfile.write("];\n")
        magma_outfile.close()


if __name__ == "__main__":

    inner = "bls"
    n = None
    list_n = [12, 24]
    k = 6
    D = 3
    max_pbits = None
    min_pbits = None
    max_curves = None
    smallest_hy = False
    min_ht = 0
    hy_0 = False
    i = None
    outfile = None
    args=sys.argv
    if len(args) <= 1:
        n=24
        print("BLS24 by default, switch to BLS12 with --bls12, switch to BN with --bn")
    options = ["--bn", "--bls12", "--bls24", "-n", "-k", "-D", "-i", "--pbits", "--minht", "-o", "--minpbits", "--maxcurves", "--smallest-hy", "--hy0"]
    if len(args) == 1:
        print("examples:")
        print("BW6-BN:")
        print("sage find_all_2chain_BLSn_BWk.sage --bn -i 0 -k 6 --pbits 512 -o testvector_bn_254_bw6_512")
        print("sage find_all_2chain_BLSn_BWk.sage --bn -i 1 -k 6 --pbits 768 -o testvector_bn_383_bw6_768")
        print("BW6-BLS12:")
        print("sage find_all_2chain_BLSn_BWk.sage -i 0 -n 12 -k 6 --pbits 768 -o testvector_bls12_377_bw6_768")
        print("sage find_all_2chain_BLSn_BWk.sage -i 3 -n 12 -k 6 --pbits 768 -o testvector_bls12_379_bw6_768")
        print("BW6-BLS24:")
        print("sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 6 --pbits 640 -o testvector_bls24_315_bw6_640")
        print("sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 6 --pbits 672 --minpbits 672 --hy0 --minht 5006110 -o testvector_bls24_315_bw6_672_hy0")
        print("CPk BLS12:")
        print("sage find_all_2chain_BLSn_BWk.sage -i 0 -n 12 -k 12 -D 3 --pbits 768 -o testvector_bls12_377_cp12D3_768")
        print("sage find_all_2chain_BLSn_BWk.sage -i 0 -n 12 -k 8 -D 1 --pbits 768 -o testvector_bls12_377_cp8D1_768")
        print("CPk BLS24:")
        print("sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 12 -D 3 --pbits 640 -o testvector_bls24_315_cp12d3_640")
        print("sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 8 -D 1 --pbits 640 -o testvector_bls24_315_cp8d1_640")

    j = 1
    while j < len(args):
        if args[j] == "--bls12":
            n=12
            inner = "bls"
            j = j+1
            continue
        if args[j] == "--bls24":
            n=24
            inner = "bls"
            j=j+1
            continue
        if args[j] == "--bn":
            n=12
            inner = "bn"
            j = j+1
            continue
        if args[j] == "--smallest-hy":
            smallest_hy=True
            j=j+1
            continue
        if args[j] == "--hy0":
            hy_0=True
            j=j+1
            continue
        if args[j] == "-n" and j+1 < len(args):
            n=int(args[j+1])
            j=j+2
            continue
        if args[j] == "-k" and j+1 < len(args):
            k=int(args[j+1])
            j=j+2
            continue
        if args[j] == "-D" and j+1 < len(args):
            D=int(args[j+1])
            j=j+2
            continue
        if args[j] == "-i" and j+1 < len(args):
            i=int(args[j+1])
            j=j+2
            continue
        if args[j] == "--pbits" and j+1 < len(args):
            max_pbits=int(args[j+1])
            j=j+2
            continue
        if args[j] == "--minpbits" and j+1 < len(args):
            min_pbits=int(args[j+1])
            j=j+2
            continue
        if args[j] == "--minht" and j+1 < len(args):
            min_ht=int(args[j+1])
            j=j+2
            continue
        if args[j] == "--maxcurves" and j+1 < len(args):
            max_curves=int(args[j+1])
            j=j+2
            continue
        if args[j] == "-o" and j+1 < len(args):
            outfile=args[j+1]
            j=j+2
            continue
        j = j+1
    if k==6 and D==3:
        BW = True
        CP = False
        print("{}{} with Brezing-Weng k={}".format(inner, n, k))
    else:
        BW = False
        CP = True
        print("{}{} with Cocks-Pinch k={} D={}".format(inner, n, k, D))

    if inner == "bn":
        test_vector = testvector_bn
        n = 12
    elif inner == "bls" and n is not None and n==12:
        test_vector = testvector_bls12
    else:
        test_vector = testvector_bls24
        n = 24
    if i is None:
        i_start = 0
        i_stop = len(test_vector)
    else:
        i_start = i
        i_stop = i+1

    if max_pbits is None:
        if n == 12:
            max_pbits = 768
        else:
            max_pbits = 640
    for i in range(i_start, i_stop):
        vec = test_vector[i]
        u = vec['u']
        if BW:
            compute_params_BN_BLSn_BWk(n, k, D, u, inner, max_pbits, min_pbits=min_pbits, min_ht=min_ht,  outfilename=outfile, max_curves=max_curves, smallest_hy=smallest_hy, hy_0=hy_0)
        else:
            compute_params_BN_BLSn_CPk(n, k, D, u, inner, max_pbits, outfilename=outfile)


"""
BW6-BN
sage find_all_2chain_BLSn_BWk.sage -i 0 --bn -k 6 --pbits 512 -o testvector_bn_254_bw6_512 > testvector_bn_254_bw6_512.out
sage find_all_2chain_BLSn_BWk.sage -i 1 --bn -k 6 --pbits 516 -o testvector_bn_254e_bw6_516 > testvector_bn_254e_bw6_516.out
sage find_all_2chain_BLSn_BWk.sage -i 6 --bn -k 6 --pbits 896 -o testvector_bn_446_bw6_896 > testvector_bn_446_bw6_896.out

sage find_all_2chain_BLSn_BWk.sage -i 2 --bn -k 6 --pbits 512 -o testvector_bn_253_bw6_512 > testvector_bn_253_bw6_512.out
sage find_all_2chain_BLSn_BWk.sage -i 3 --bn -k 6 --pbits 512 -o testvector_bn_253b_bw6_512 > testvector_bn_253b_bw6_512.out
sage find_all_2chain_BLSn_BWk.sage -i 4 --bn -k 6 --pbits 768 -o testvector_bn_382_bw6_768 > testvector_bn_382_bw6_768.out
sage find_all_2chain_BLSn_BWk.sage -i 5 --bn -k 6 --pbits 770 -o testvector_bn_382b_bw6_770 > testvector_bn_382b_bw6_768.out

BW6-BLS12
sage find_all_2chain_BLSn_BWk.sage -i 0 -n 12 -k 6 --pbits 768 -o testvector_bls12_377_bw6_768 > testvector_bls12_377_bw6_768.out
sage find_all_2chain_BLSn_BWk.sage -i 3 -n 12 -k 6 --pbits 768 -o testvector_bls12_379_bw6_768 > testvector_bls12_379_bw6_768.out

sage find_all_2chain_BLSn_BWk.sage -i 1 -n 12 -k 6 --pbits 768 -o testvector_bls12_381_bw6_768 > testvector_bls12_381_bw6_768.out

BW6-BLS24
sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 6 --pbits 640 -o testvector_bls24_315_bw6_640 > testvector_bls24_315_bw6_640.out
sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 6 --pbits 672 --minpbits 672 --hy0 --minht 5006110 -o testvector_bls24_315_bw6_672_hy0 > testvector_bls24_315_bw6_672_hy0.out

CPk BLS12
sage find_all_2chain_BLSn_BWk.sage -i 0 -n 12 -k 12 -D 3 --pbits 768 -o testvector_bls12_377_cp12d3_768 > testvector_bls12_377_cp12d3_768.out
sage find_all_2chain_BLSn_BWk.sage -i 0 -n 12 -k 8 -D 1 --pbits 768 -o testvector_bls12_377_cp8d1_768 > testvector_bls12_377_cp8d1_768.out

sage find_all_2chain_BLSn_BWk.sage -i 3  -n 12 -k 12 -D 3 --pbits 768 -o testvector_bls12_379_cp12d3_768 > testvector_bls12_379_cp12d3_768.out
sage find_all_2chain_BLSn_BWk.sage -i 3  -n 12 -k 8 -D 1 --pbits 768 -o testvector_bls12_379_cp8d1_768 > testvector_bls12_379_cp8d1_768.out

CPk BLS24
sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 12 -D 3 --pbits 640 -o testvector_bls24_315_cp12d3_640 > testvector_bls24_315_cp12d3_640.out
sage find_all_2chain_BLSn_BWk.sage -i 3 -n 24 -k 8 -D 1 --pbits 640 -o testvector_bls24_315_cp8d1_640 > testvector_bls24_315_cp8d1_640.out

BLS12-377:
(0,ht= 13,hy=  9) 761 bits q=0x122e824fb83ce0ad187c94004faff3eb926186a81d14688528275ef8087be41707ba638e584e91903cebaff25b423048689c8ed12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf0a437e6913e6870000082f49d00000000008b
Q = (103*x^12 - 379*x^11 + 250*x^10 + 691*x^9 - 911*x^8 - 79*x^7 + 623*x^6 - 640*x^5 + 274*x^4 + 763*x^3 + 73*x^2 + 254*x + 229)/9

(0,ht= 20,hy=  2) 761 bits q=0x122e824fb83ce0ad187c94004faff3eb926186a81d14688528275ef8087be41707ba638e584e91903cebaff25b423048689c8ed12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf0a437e6913e6870000082f49d00000000008b
Q = (103*x^12 - 379*x^11 + 250*x^10 + 691*x^9 - 911*x^8 - 79*x^7 + 623*x^6 - 640*x^5 + 274*x^4 + 763*x^3 + 73*x^2 + 254*x + 229)/9

BLS12-381: only 2 solutions (obtained each with two combinations of ht, hy)
(1,ht= -4,hy= -6) 767 bits q=0x51e2bcf25fa8992238259ea59a063294c36dc4098befce4230f8d18f41e3fc19665e4360b872007d3dd5a1b865cbe8dadc2ce0c034926d18fe0ef8c1c63df7d97cbc118805598e5c31732000974254c83a38b08e7179beb96896aaaec71538e7
(1,ht=  5,hy=  9) 768 bits q=0xb0fa901c5b221128c4c3b0f085b5b5b1ed9d5874077f657dd6c9c1fffadb1e7a2a95a82321132a2d44ec3ba76fdf824db287dcc1e911bacc2aeff4bcee70d9ba88b5d1c34f1d79820c42f9f9eeca76288f578cabaf57aac7aa7caab1c70e38eb
(1,ht=  7,hy=  5) 767 bits q=0x51e2bcf25fa8992238259ea59a063294c36dc4098befce4230f8d18f41e3fc19665e4360b872007d3dd5a1b865cbe8dadc2ce0c034926d18fe0ef8c1c63df7d97cbc118805598e5c31732000974254c83a38b08e7179beb96896aaaec71538e7
(1,ht=-11,hy= -7) 768 bits q=0xb0fa901c5b221128c4c3b0f085b5b5b1ed9d5874077f657dd6c9c1fffadb1e7a2a95a82321132a2d44ec3ba76fdf824db287dcc1e911bacc2aeff4bcee70d9ba88b5d1c34f1d79820c42f9f9eeca76288f578cabaf57aac7aa7caab1c70e38eb


"""
