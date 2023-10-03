from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from external.Pairings.pairing import bits_2naf

def karabina_cost_exp(u0, mk, mkd, skd, ikd):
    """ Karabina cost of compressed squaring
    https://www.ams.org/journals/mcom/2013-82-281/S0025-5718-2012-02625-1/
    mk = multiplication in Fpk
    mkd = multiplication in F_{p^{k/d}}
    skd = squaring in F_{p^{k/d}}
    ikd = inversion in F_{p^{k/d}}
    """
    bits_u0 = (abs(u0)).digits(2)
    len_u = len(bits_u0)
    hw_u = sum(bits_u0)
    bits_2naf_u0 = bits_2naf(abs(u0))
    len2naf_u = len(bits_2naf_u0)
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])
    cost_exp_to_u = 4*(len_u-1)*mkd + (6*(hw_u - 1) -3)*mkd + (hw_u - 1)*mk + 3*(hw_u - 1)*skd + ikd
    cost_exp_to_u_2naf = 4*(len2naf_u-1)*mkd + (6*(hw2naf_u - 1) -3)*mkd + (hw2naf_u - 1)*mk + 3*(hw2naf_u - 1)*skd + ikd
    return min(cost_exp_to_u, cost_exp_to_u_2naf)

def cost_final_exp_hard_bw6_bls12(u0, ht, hy, trace_0_mod_r_mod_u=True):
    """
    trace_0_mod_r_mod_u is True:
    final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_alt(m, u, ht, hy):
    cost e((u-1)/3) + 5e(u-1) +3e(u+1) + e(d1) + e(d2) + 14M + 2S + frb + 4cj
    final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u(m, u, ht, hy):
    cost e((u-1)/3) + 5e(u-1) +3e(u+1) + e(d1) + e(d2) + 15M + 2S + frb + 4cj

    trace_0_mod_r_mod_u is False: one less conjugate (which is actully free)
    final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u(m, u, ht, hy):
    cost e((u-1)/3) + 5e(u-1) + 3e(u+1) + e(d1) + e(d2) + 14M + 2S + frb + 3cj
    final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_alt(m, u, ht, hy):
    cost e((u-1)/3) + 5e(u-1) + 3e(u+1) + e(d1) + e(d2) + 15M + 2S + frb + 4cj

    5*exp(u-1) + exp((u-1)//3) + 3*exp(u+1) + exp(d2) + exp(d1) + 14 M + 2 S + frob + 3(or 4) conj
    """
    c1 = ZZ(abs(ht**2 + 3*hy**2)//4)
    if trace_0_mod_r_mod_u:
        c2 = ZZ(abs(ht - hy)//2)
    else:
        c2 = ZZ(abs(ht + hy)//2)
    u0 = ZZ(abs(u0))
    m = 1
    s = 1
    i = 25*m
    i6 = 34*m + i
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    s6_cyclo = 3*s2 # sk_cyclo = 3*s_{k/3} for 6 | k
    f3 = 0 # Frobenius powers x^(q^3) in Fp6 are only negation
    f = 4*m # assume a frobenius power costs 4 m because the 6-th roots of unity are of the form (1, w, w^2=-w, w^3=-1, w^4, w^5)
    bits_c1 = ZZ(abs(c1)).bits()
    hw_c1 = sum([1 for bi in bits_c1 if bi != 0])
    bits_c2 = ZZ(abs(c2)).bits()
    hw_c2 = sum([1 for bi in bits_c2 if bi != 0])
    bits_2naf_c1 = bits_2naf(c1)
    hw2naf_c1 = sum([1 for bi in bits_2naf_c1 if bi != 0])
    bits_2naf_c2 = bits_2naf(c2)
    hw2naf_c2 = sum([1 for bi in bits_2naf_c2 if bi != 0])

    cost_exp_c1 = (len(bits_c1)-1) * s6_cyclo + (hw_c1-1) * m6
    cost_exp_c2 = (len(bits_c2)-1) * s6_cyclo + (hw_c2-1) * m6

    cost_exp_c1_2naf = (len(bits_2naf_c1)-1) * s6_cyclo + (hw2naf_c1-1) * m6
    cost_exp_c2_2naf = (len(bits_2naf_c2)-1) * s6_cyclo + (hw2naf_c2-1) * m6

    cost_easy_part = f3 + f + 2*m6 + i6 # f3 + i6 + 2*m6 + f
    print("cost final exp easy: {}m".format(cost_easy_part))
    list_ui = [u0+1, u0-1, (u0-1)//3]
    idx = {"u+1": 0, "u-1":1, "(u-1)/3":2}
    bits_ui = [ZZ(abs(ui)).bits() for ui in list_ui]
    hw_ui = [sum((abs(ui)).digits(2)) for ui in list_ui]
    bits_2naf_ui = [bits_2naf(abs(ui)) for ui in list_ui]
    hw2naf_ui = [sum([1 for bi in bui if bi != 0]) for bui in bits_2naf_ui]
    cost_exp_ui = [(len(bits_ui[j])-1) * s6_cyclo + (hw_ui[j]-1) * m6 for j in range(len(list_ui))]
    cost_exp_ui_2naf = [(len(bits_2naf_ui[j])-1) * s6_cyclo + (hw2naf_ui[j]-1) * m6 for j in range(len(list_ui))]

    if trace_0_mod_r_mod_u:
        # 5*exp(u-1) + exp((u-1)//3) + 3*exp(u+1) + exp(d2) + exp(d1) + 14 M + 2 S + 1 frob + 4 conj
        cost_hard_a      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 14*m6 + 2*s6 + f
        cost_hard_2naf_a = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 14*m6 + 2*s6 + f
        # 5*exp(u-1) + exp((u-1)//3) + 3*exp(u+1) + exp(d2) + exp(d1) + 15 M + 2 S + 1 frob + 4 conj
        cost_hard_b      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 15*m6 + 2*s6 + f
        cost_hard_2naf_b = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 15*m6 + 2*s6 + f
    else:
        # 5*exp(u-1) + exp((u-1)//3) + 3*exp(u+1) + exp(d2) + exp(d1) + 14 M + 2 S + 1 frob + 4 conj
        cost_hard_a      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 14*m6 + 2*s6 + f
        cost_hard_2naf_a = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 14*m6 + 2*s6 + f
        # 5*exp(u-1) + exp((u-1)//3) + 3*exp(u+1) + exp(d2) + exp(d1) + 15 M + 2 S + 1 frob + 3 conj
        cost_hard_b      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 15*m6 + 2*s6 + f
        cost_hard_2naf_b = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 15*m6 + 2*s6 + f

    print("cost final exp hard: {}m (a) (with cyclotomic squaring)".format(cost_hard_a))
    print("cost final exp hard: {}m (b) (with cyclotomic squaring)".format(cost_hard_b))
    print("cost final exp:      {}m".format(cost_easy_part + min(cost_hard_a, cost_hard_b)))
    print("cost final exp hard: {}m (a, 2-naf) (with cyclotomic squaring)".format(cost_hard_2naf_a))
    print("cost final exp hard: {}m (b, 2-naf) (with cyclotomic squaring)".format(cost_hard_2naf_b))
    print("cost final exp 2naf: {}m".format(cost_easy_part + min(cost_hard_2naf_a, cost_hard_2naf_b)))

def cost_final_exp_hard_bw6_bls24(u0, ht, hy, trace_0_mod_r_mod_u=True):
    """
    trace_0_mod_r_mod_u is True:
    final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u
    cost e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 15M + 2S + frb + 4cj
    final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt
    cost e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 14M + 2S + frb + 4cj

    final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u
    cost e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 15M + 2S + frb + 2cj
    final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt
    cost e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 16M + 2S + frb + 3cj
    """
    c1 = ZZ(abs(ht**2 + 3*hy**2)//4)
    if trace_0_mod_r_mod_u:
        c2 = ZZ(abs(ht - hy)//2)
    else:
        c2 = ZZ(abs(ht + hy)//2)
    u0 = ZZ(abs(u0))
    m = 1
    s = 1
    i = 25*m
    i6 = 34*m + i
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    s6_cyclo = 3*s2 # sk_cyclo = 3*s_{k/3} for 6 | k
    f3 = 0 # Frobenius powers x^(q^3) in Fp6 are only negation
    f = 4*m # assume a frobenius power costs 4 m because the 6-th roots of unity are of the form (1, w, w^2=-w, w^3=-1, w^4, w^5)
    bits_c1 = ZZ(abs(c1)).bits()
    hw_c1 = sum([1 for bi in bits_c1 if bi != 0])
    bits_c2 = ZZ(abs(c2)).bits()
    hw_c2 = sum([1 for bi in bits_c2 if bi != 0])
    bits_2naf_c1 = bits_2naf(c1)
    hw2naf_c1 = sum([1 for bi in bits_2naf_c1 if bi != 0])
    bits_2naf_c2 = bits_2naf(c2)
    hw2naf_c2 = sum([1 for bi in bits_2naf_c2 if bi != 0])

    cost_exp_c1 = (len(bits_c1)-1) * s6_cyclo + (hw_c1-1) * m6
    cost_exp_c2 = (len(bits_c2)-1) * s6_cyclo + (hw_c2-1) * m6

    cost_exp_c1_2naf = (len(bits_2naf_c1)-1) * s6_cyclo + (hw2naf_c1-1) * m6
    cost_exp_c2_2naf = (len(bits_2naf_c2)-1) * s6_cyclo + (hw2naf_c2-1) * m6

    cost_easy_part = f3 + f + 2*m6 + i6 # f3 + i6 + 2*m6 + f
    print("cost final exp easy: {}m".format(cost_easy_part))

    list_ui = [u0+1, u0-1, (u0-1)//3, u0, u0**2+1]
    idx = {"u+1": 0, "u-1":1, "(u-1)/3":2, "u":3, "u^2+1":4}
    bits_ui = [ZZ(abs(ui)).bits() for ui in list_ui]
    hw_ui = [sum((abs(ui)).digits(2)) for ui in list_ui]
    bits_2naf_ui = [bits_2naf(abs(ui)) for ui in list_ui]
    hw2naf_ui = [sum([1 for bi in bui if bi != 0]) for bui in bits_2naf_ui]
    cost_exp_ui = [(len(bits_ui[j])-1) * s6_cyclo + (hw_ui[j]-1) * m6 for j in range(len(list_ui))]
    cost_exp_ui_2naf = [(len(bits_2naf_ui[j])-1) * s6_cyclo + (hw2naf_ui[j]-1) * m6 for j in range(len(list_ui))]

    if trace_0_mod_r_mod_u:
        # e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 14M + 2S + frb + 4cj
        cost_hard_a      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 14*m6 + 2*s6 + f
        cost_hard_2naf_a = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 14*m6 + 2*s6 + f
        cost_hard_a      += 3*min(cost_exp_ui[idx["u^2+1"]], 2*cost_exp_ui[idx["u"]]+m6)
        cost_hard_2naf_a += 3*min(cost_exp_ui_2naf[idx["u^2+1"]], 2*cost_exp_ui_2naf[idx["u"]]+m6)
        # e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 15M + 2S + frb + 4cj
        cost_hard_b      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 15*m6 + 2*s6 + f
        cost_hard_2naf_b = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 15*m6 + 2*s6 + f
        cost_hard_b      += 3*min(cost_exp_ui[idx["u^2+1"]], 2*cost_exp_ui[idx["u"]]+m6)
        cost_hard_2naf_b += 3*min(cost_exp_ui_2naf[idx["u^2+1"]], 2*cost_exp_ui_2naf[idx["u"]]+m6)
    else:
        # e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 15M + 2S + frb + 2cj
        cost_hard_a      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 15*m6 + 2*s6 + f
        cost_hard_2naf_a = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 15*m6 + 2*s6 + f
        cost_hard_a      += 3*min(cost_exp_ui[idx["u^2+1"]], 2*cost_exp_ui[idx["u"]]+m6)
        cost_hard_2naf_a += 3*min(cost_exp_ui_2naf[idx["u^2+1"]], 2*cost_exp_ui_2naf[idx["u"]]+m6)
        # e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 16M + 2S + frb + 3cj
        cost_hard_b      = 5*cost_exp_ui[idx["u-1"]]      + cost_exp_ui[idx["(u-1)/3"]]      + 3*cost_exp_ui[idx["u+1"]]      + cost_exp_c1      + cost_exp_c2      + 16*m6 + 2*s6 + f
        cost_hard_2naf_b = 5*cost_exp_ui_2naf[idx["u-1"]] + cost_exp_ui_2naf[idx["(u-1)/3"]] + 3*cost_exp_ui_2naf[idx["u+1"]] + cost_exp_c1_2naf + cost_exp_c2_2naf + 16*m6 + 2*s6 + f
        cost_hard_b      += 3*min(cost_exp_ui[idx["u^2+1"]], 2*cost_exp_ui[idx["u"]]+m6)
        cost_hard_2naf_b += 3*min(cost_exp_ui_2naf[idx["u^2+1"]], 2*cost_exp_ui_2naf[idx["u"]]+m6)

    print("min(e(u^2+1), 2*e(u)+m6) = min({}, {}) = {}".format(cost_exp_ui[idx["u^2+1"]], 2*cost_exp_ui[idx["u"]]+m6, min(cost_exp_ui[idx["u^2+1"]], 2*cost_exp_ui[idx["u"]]+m6)))
    print("min(e(u^2+1), 2*e(u)+m6) = min({}, {}) = {} (2-naf)".format(cost_exp_ui_2naf[idx["u^2+1"]], 2*cost_exp_ui_2naf[idx["u"]]+m6, min(cost_exp_ui_2naf[idx["u^2+1"]], 2*cost_exp_ui_2naf[idx["u"]]+m6)))
    print("cost final exp hard: {}m (a) (with cyclotomic squaring)".format(cost_hard_a))
    print("cost final exp hard: {}m (b) (with cyclotomic squaring)".format(cost_hard_b))
    print("cost final exp:      {}m".format(cost_easy_part + min(cost_hard_a, cost_hard_b)))
    print("cost final exp hard: {}m (a, 2-naf) (with cyclotomic squaring)".format(cost_hard_2naf_a))
    print("cost final exp hard: {}m (b, 2-naf) (with cyclotomic squaring)".format(cost_hard_2naf_b))
    print("cost final exp 2naf: {}m".format(cost_easy_part + min(cost_hard_2naf_a, cost_hard_2naf_b)))

def cost_pairing_bn(u):
    u = ZZ(abs(u))
    m = 1
    s = 1
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    m12 = 3*m6
    s12 = 2*m6
    f2 = 8*m
    f12 = 10*m
    inv = 25*m
    i12 = 97*m+inv # 94*m + i in Guillevic Masson Thome DCC'20 is not compatible with the tower
    i2 = 2*m + 2*s + inv
    k = 12
    e = 2
    d = 6
    assert e*d == k
    # AKLGL
    mk=m12; sk=s12; ik=i12
    me=m2 ; se=s2; sk=s12 # e = k/d
    double_line_ate = 3*me+6*se+(k//3)*m
    add_line_ate    = 11*me+2*se+(k//3)*m
    # Costello Lange Naehrig
    double_line_ate_cln = 2*me+7*se+(k//3)*m
    add_line_ate_cln    = 10*me+2*se+(k//3)*m
    sparse_mult        = 13*me
    sparse_sparse      = 6*me

    v = 6*u+2
    hw_v = sum((abs(v)).digits(2))
    bits_2naf_v = bits_2naf(abs(v))
    len2naf_v = len(bits_2naf_v)
    hw2naf_v = sum([1 for bi in bits_2naf_v if bi != 0])

    cost_ate = (v.nbits()-1)*(double_line_ate+sk) - sk + (v.nbits()-2 -(hw_v-1))*sparse_mult + (hw_v-1)*(add_line_ate+sparse_sparse+mk)
    cost_ate_2naf = (len(bits_2naf_v)-1)*(double_line_ate+sk) - sk + (len(bits_2naf_v)-2 -(hw2naf_v-1))*sparse_mult + (hw2naf_v-1)*(add_line_ate+sparse_sparse+mk)
    cost_ate_cln = (v.nbits()-1)*(double_line_ate_cln+sk) - sk + (v.nbits()-2 - (hw_v-1))*sparse_mult + (hw_v-1)*(add_line_ate_cln+sparse_sparse+mk)
    cost_ate_cln_2naf = (len(bits_2naf_v)-1)*(double_line_ate_cln+sk) - sk + (len(bits_2naf_v)-2 -(hw2naf_v-1))*sparse_mult + (hw2naf_v-1)*(add_line_ate_cln+sparse_sparse+mk)
    fp_Fp2 = 0 # this is a conjugation
    pi_Q = 2*me
    pi2_Q = 2*(2*m)
    additional_terms = pi_Q + pi2_Q + 2*(add_line_ate) + sparse_sparse + mk
    cost_ate += additional_terms
    cost_ate_2naf += additional_terms
    cost_ate_cln += additional_terms
    cost_ate_cln_2naf += additional_terms

    print("m12 = {}m s12 = {}m".format(m12,s12))
    print("m2 = {}m s2 = {}m".format(m2,s2))

    print("cost ate Miller       = {}m".format(cost_ate))
    print("cost ate Miller 2-naf = {}m".format(cost_ate_2naf))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk)".format(v.nbits(),hw_v,e))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk) (2-naf)".format(len(bits_2naf_v),hw2naf_v,e))
    print("cost ate Miller       = {}m (Costello Lange Naehrig)".format(cost_ate_cln))
    print("cost ate Miller 2-naf = {}m (Costello Lange Naehrig)".format(cost_ate_cln_2naf))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk)".format(v.nbits(),hw_v,e))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk) (2-naf)".format(len(bits_2naf_v),hw2naf_v,e))

    min_cost_miller = min(cost_ate, cost_ate_2naf, cost_ate_cln, cost_ate_cln_2naf)

    # final exponentiation
    # (p^12-1)/r = (p^12-1)/Phi_12(u) * Phi_12(u)/r = (p^6-1)*(p^2+1)*(p^4-p^2+1)/r
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    assert ((x**12-1) // cyclotomic_polynomial(12)) == (x**6-1)*(x**2+1)
    # now cost of exponentiation
    # easy part: exponent = (px**6-1)*(px**2+1)
    # px**6 is a conjugation: costs only 6 negations
    # one inversion in Fp12, one mult
    #f2_12 = 8*m # from cost_pairing.py DCC'2020
    # Simon Masson PhD thesis section 4.4.2
    f2_12 = 5*m2
    easy_part = i12 + m12 + f2_12 + m12
    # hard part: exponent = (px**4-px**2+1)//rx
    # assume that cyclotomic squarings are implemented (PKC'10 Granger Scott)
    # cost_exp_to_x = 4*(len2naf_u-1)*m2 + (6*(hw2naf_u - 1) -3)*m2 + (hw2naf_u - 1)*m12 + 3*(hw2naf_u - 1)*s2 + i2
    # Formulas from Fuentes-Castaneda et al, SAC'2011
    # cost: 3 exp(u) + 3 S + 10 M + 3 frob
    cost_exp_to_x = karabina_cost_exp(u, m12, m2, s2, i2)
    s12cyclo = 18*m
    hard_part = 3*cost_exp_to_x + 10*m12 + 3*f12 + 3*s12cyclo
    print("cost final exp easy: {}m".format(easy_part))
    print("cost final exp hard: {}m (with Karabina technique)".format(hard_part))
    print("cost final exp:      {}m".format(easy_part + hard_part))
    print("\ncost pairing (total) {}m".format(easy_part + hard_part + min_cost_miller))

def cost_pairing_bls12(u0):
    # cost
    # a=0
    # AKLGL'11
    # there are three formulas in https://www.lirmm.fr/arith18/papers/Chung-Squaring.pdf for squaring in cubic extensions
    # S3 = m + 4*s; S3 = 2*m+3*s; S3 = 3*m+2*s
    u0 = ZZ(abs(u0))
    m = 1
    s = 1
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    m12 = 3*m6
    s12 = 2*m6
    f2 = 8*m
    f12 = 10*m
    inv = 25*m
    i12 = 97*m+inv # 94*m + i in Guillevic Masson Thome DCC'20 is not compatible with the tower
    i2 = 2*m + 2*s + inv
    k = 12
    e = 2
    d = 6
    assert e*d == k
    # AKLGL
    mk=m12; sk=s12; ik=i12
    me=m2 ; se=s2; sk=s12 # e = k/d
    double_line_ate = 3*me+6*se+(k//3)*m
    add_line_ate    = 11*me+2*se+(k//3)*m
    # Costello Lange Naehrig
    double_line_ate_cln = 2*me+7*se+(k//3)*m
    add_line_ate_cln    = 10*me+2*se+(k//3)*m
    sparse_dense        = 13*me
    sparse_sparse       = 6*me
    update1        = 13*me+sk
    update2        = 13*me
    
    hw_u = sum((abs(u0)).digits(2))
    bits_2naf_u0 = bits_2naf(abs(u0))
    len2naf_u = len(bits_2naf_u0)
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])

    cost_ate = (u0.nbits()-1)*(double_line_ate+sk) - sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate+sparse_sparse+mk)
    cost_ate_2naf = (len(bits_2naf_u0)-1)*(double_line_ate+sk) - sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate+sparse_sparse+mk)
    cost_ate_cln = (u0.nbits()-1)*(double_line_ate_cln+sk) - sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate_cln+sparse_sparse+mk)
    cost_ate_cln_2naf = (len(bits_2naf_u0)-1)*(double_line_ate_cln+sk) - sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate_cln+sparse_sparse+mk)
    print("m12 = {}m s12 = {}m".format(m12,s12))
    print("m2 = {}m s2 = {}m".format(m2,s2))

    print("cost ate Miller       = {}m".format(cost_ate))
    print("cost ate Miller 2-naf = {}m".format(cost_ate_2naf))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk)".format(u0.nbits(),hw_u,e))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk) (2-naf)".format(len(bits_2naf_u0),hw2naf_u,e))
    print("cost ate Miller       = {}m (Costello Lange Naehrig)".format(cost_ate_cln))
    print("cost ate Miller 2-naf = {}m (Costello Lange Naehrig)".format(cost_ate_cln_2naf))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk)".format(u0.nbits(),hw_u,e))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m + sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m + 6*m{2}+mk) (2-naf)".format(len(bits_2naf_u0),hw2naf_u,e))

    min_cost_miller = min(cost_ate, cost_ate_2naf, cost_ate_cln, cost_ate_cln_2naf)
    
    # final exponentiation
    # (p^12-1)/r = (p^12-1)/Phi_12(u) * Phi_12(u)/r = (p^6-1)*(p^2+1)*(p^4-p^2+1)/r
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    assert ((x**12-1) // cyclotomic_polynomial(12)) == (x**6-1)*(x**2+1)
    px = (x-1)**2*(x**4 - x**2 + 1)/3 + x
    rx = x**4-x**2+1
    l3 = (x-1)**2
    l2 = l3*x
    l1 = l2*x-l3
    l0 = l1*x + 3
    exponent = (px**4-px**2+1)//rx
    assert l0+px*(l1+px*(l2+px*l3)) == 3*exponent
    # now cost of exponentiation
    # easy part: exponent = (px**6-1)*(px**2+1)
    # px**6 is a conjugation: costs only 6 negations
    # one inversion in Fp12, one mult
    #f2_12 = 8*m # from cost_pairing.py DCC'2020
    # Simon Masson PhD thesis section 4.4.2
    f2_12 = 5*m2
    easy_part = i12 + m12 + f2_12 + m12
    # hard part: exponent = (px**4-px**2+1)//rx
    # assume that cyclotomic squarings are implemented (PKC'10 Granger Scott)
    # cost_exp_to_x = 4*(len2naf_u-1)*m2 + (6*(hw2naf_u - 1) -3)*m2 + (hw2naf_u - 1)*m12 + 3*(hw2naf_u - 1)*s2 + i2
    cost_exp_to_x = karabina_cost_exp(u0, m12, m2, s2, i2)
    s12cyclo = 18*m
    hard_part = 5*cost_exp_to_x + 10*m12 + 3*f12 + 2*s12cyclo
    print("cost final exp easy: {}m".format(easy_part))
    print("cost final exp hard: {}m (with Karabina technique)".format(hard_part))
    print("cost final exp:      {}m".format(easy_part + hard_part))
    print("\ncost pairing (total) {}m".format(easy_part + hard_part + min_cost_miller))

def cost_pairing_bls12_377():
    u0 = ZZ(0x8508C00000000001)
    cost_pairing_bls12(u0)
def cost_pairing_bls12_379():
    u0 = ZZ(0x9b04000000000001)
    cost_pairing_bls12(u0)
def cost_pairing_bls12_381():
    u0 = ZZ(-(2**63+2**62+2**60+2**57+2**48+2**16))
    cost_pairing_bls12(u0)
def cost_pairing_bls12_383():
    u0 = ZZ(0x105a8000000000001)
    cost_pairing_bls12(u0)

def cost_pairing_bls12_440():
    u0 = ZZ(-(2**73+2**72+2**50+2**24))
    cost_pairing_bls12(u0)

def cost_pairing_bls12_442():
    u0 = ZZ(-(2**12-2**48-2**71+2**74))
    cost_pairing_bls12(u0)

def cost_pairing_bls12_446():
    u0 = ZZ(-(2**74+2**73+2**63+2**57+2**50+2**17+1))
    cost_pairing_bls12(u0)

def cost_pairing_bls24(u0):
    # cost
    # a=0
    # AKLGL'11
    # there are three formulas in https://www.lirmm.fr/arith18/papers/Chung-Squaring.pdf for squaring in cubic extensions
    # S3 = m + 4*s; S3 = 2*m+3*s; S3 = 3*m+2*s
    u0 = ZZ(abs(u0))
    m = 1
    s = 1
    m2 = 3*m
    s2 = 2*m
    m4 = 3*m2
    s4 = 2*m2
    m12 = 6*m4 # multiplication in cubic extension
    #s12 = 1*m4+4*s4 # squaring in cubic extension
    s12 = 2*m4+3*s4 # squaring in cubic extension
    m24 = 3*m12
    s24 = 2*m12
    s24cyclo = 18*m2 # = 6*m4
    f24 = 22*m # 24 - 2 TODO
    inv = 25*m
    i2 = 4*m + inv
    i4 = 2*m2 + 2*s2 + i2
    i12 = 9*m4 + 3*s4 + i4
    i24 = 2*m12 + 2*s12 + i12
    print("m24 = {}m s24 = {}m, f24 = {}m, s24cyclo = {}m, i24={}m+i".format(m24, s24, f24, s24cyclo, i24))
    print("i4 = {}m".format(i4))
    k = 24
    d = 6
    e = 4
    assert e*d == k
    # AKLGL
    me=m4 ; se=s4; mk=m24; sk=s24 # e = k/d
    double_line_ate = 3*me+6*se+(k//3)*m
    add_line_ate    = 11*me+2*se+(k//3)*m
    # Costello Lange Naehrig
    # double_line_ate_cln = 3*me+5*se+(k//3)*m    ???
    double_line_ate_cln = 2*me+7*se+(k//3)*m
    add_line_ate_cln    = 10*me+2*se+(k//3)*m
    sparse_dense        = 13*me
    sparse_sparse       = 6*me
    update1        = 13*me+sk
    update2        = 13*me
    
    hw_u = sum((abs(u0)).digits(2))
    bits_2naf_u0 = bits_2naf(abs(u0))
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])

    cost_ate = (u0.nbits()-1)*(double_line_ate+sk) -sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate+sparse_sparse+mk)
    cost_ate_2naf = (len(bits_2naf_u0)-1)*(double_line_ate+sk) -sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate+sparse_sparse+mk)
    cost_ate_cln = (u0.nbits()-1)*(double_line_ate_cln+sk) -sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate_cln+sparse_sparse+mk)
    cost_ate_cln_2naf = (len(bits_2naf_u0)-1)*(double_line_ate_cln+sk) -sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate_cln+sparse_sparse+mk)
    print("m24 = {}m s24 = {}m".format(m24,s24))
    print("m4 = {}m s4 = {}m".format(m4,s4))
    print("cost ate Miller       = {}m".format(cost_ate))
    print("cost ate Miller 2-naf = {}m".format(cost_ate_2naf))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(u0.nbits(),hw_u,e))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(len(bits_2naf_u0),hw2naf_u,e))
    print("cost ate Miller       = {}m (Costello Lange Naehrig)".format(cost_ate_cln))
    print("cost ate Miller 2-naf = {}m (Costello Lange Naehrig)".format(cost_ate_cln_2naf))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(u0.nbits(),hw_u,e))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(len(bits_2naf_u0),hw2naf_u,e))

    min_cost_miller = min(cost_ate, cost_ate_2naf, cost_ate_cln, cost_ate_cln_2naf)
    # final exponentiation
    # (p^24-1)/r = (p^24-1)/Phi_24(u) * Phi_24(u)/r = (p^12-1)*(p^4+1)*(p^8-p^4+1)/r
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    assert ((x**24-1) // cyclotomic_polynomial(24)) == (x**12-1)*(x**4+1)
    px = (x-1)**2*(x**8 - x**4 + 1)/3 + x
    rx = x**8-x**4+1
    l7 = (x-1)**2 # x**2-2*x+1 # this is (x-1)**2
    l6 = l7*x
    l5 = l6*x
    l4 = l5*x
    l3 = l4*x-l7
    l2 = l3*x
    l1 = l2*x
    l0 = l1*x + 3
    exponent = (px**8-px**4+1)//rx
    assert l0+px*(l1+px*(l2+px*(l3+px*(l4+px*(l5+px*(l6+px*l7)))))) == 3*exponent
    assert (l7*(x*x*x*x-1)*x*x*x+3)+px*(l7*(x*x*x*x-1)*x*x+px*(l7*(x*x*x*x-1)*x+px*(l7*(x*x*x*x-1)+px*(l7*x*x*x+px*(l7*x*x+px*(l7*x+px*l7)))))) == 3*exponent    
    # cost from GF08:
    # 8 exponentiations by u, 1 exponentiation by u/2, 1 squaring, 10 multiplications, 7 Frobenius in Fp24.
    # easy part: exponent = (px**12-1)*(px**4+1)
    # f^p12 is free (conjugation)
    # 1 Frobenius, 1 inversion, 2 multiplications
    f4_24 = (24//4 - 2) * 4 # 16m
    easy_part = f4_24 + i24 + 2*m24
    # hard part
    # exponent = (x-1)
    cost_exp_to_x_1 = karabina_cost_exp(u0-1, m24, m4, s4, i4)
    cost_exp_to_x = karabina_cost_exp(u0, m24, m4, s4, i4)
    hard_part = (2*cost_exp_to_x_1) + (cost_exp_to_x + m24 + f24)*7 + m24 + (s24cyclo + 2*m24)
    # total 18732 M + 10 I with Karabina compressed squarings
    # 23400 M + I with cyclotomic squaring of Granger-Scott.
    print("cost final exp easy: {}m".format(easy_part))
    print("cost final exp hard: {}m (with Karabina technique)".format(hard_part))
    print("cost final exp:      {}m".format(easy_part + hard_part))
    print("\ncost pairing (total) {}m".format(easy_part + hard_part + min_cost_miller))

def cost_pairing_bls24_318():
    u0 = ZZ(-2**32 + 2**28 + 2**12)
    cost_pairing_bls24(u0)
def cost_pairing_bls24_318b():
    u0 = ZZ(-2**32 + 2**28 - 2**23 + 2**21 + 2**18 + 2**12 - 1)
    cost_pairing_bls24(u0)
def cost_pairing_bls24_315():
    u0 = ZZ(-2**32 + 2**30 + 2**22 - 2**20 + 1)
    cost_pairing_bls24(u0)

def cost_pairing_kss16(u0):
    # cost
    # b=0
    u0 = ZZ(u0)
    m = 1
    s = 1
    m2 = 3*m
    s2 = 2*m
    m4 = 3*m2   # = 9m
    s4 = 2*m2   # = 6m
    m8 = 3*m4   # = 27m
    s8 = 2*m4   # = 18m
    m16 = 3*m8  # = 81m
    s16 = 2*m8  # = 54m
    # s16_cyclo = 4*m4 # m4=9*m
    # Granger Scott, Faster Squaring in the Cyclotomic Subgroup of Sixth Degree Extensions,
    # PKC 2010, pages 209-223.
    # http://www.iacr.org/archive/pkc2010/60560212/60560212.pdf
    s16_cyclo = 2*s8 # = 4*m4 indeed
    f16 = 15*m # Ghammam PhD thesis Table 4.12 p. 113
    c16_cyclo = s16_cyclo + m16 #8*m4 # cube, to check: s16_cyclo + m16 = 4*m4 + 9*m4, so where is 8*m4 it from?
    # INDOCRYPT 2012, Analysis of Optimum Pairing Products at High Security Levels
    # Xusheng Zhang and Dongdai Lin, but no formula provided
    inv = 25*m
    i2 = 2*m + 2*s + inv
    i4 = 2*m2 + 2*s2 + i2
    i8 = 2*m4 + 2*s4 + i4
    i16 = 2*m8 + 2*s8 + i8
    assert i16 == 134*m+inv # cf Guillevic Masson Thome DCC'20
    k = 16
    e = 4
    d = 4
    assert e*d == k
    me=m4; se=s4; mk=m16; sk=s16 # e = k/d
    # Costello Lange Naehrig, with d=4, b=0
    double_line_ate_cln = 2*me+8*se+(k//2)*m
    add_line_ate_cln    = 9*me+5*se+(k//2)*m
    sparse_dense  = 8*me
    sparse_sparse = 6*me
    update1       = 8*me+sk
    update2       = 8*me    # in other words, a sparse multiplication costs 8 m4 (8, not 7)

    hw_u = sum((abs(u0)).digits(2))
    bits_2naf_u0 = bits_2naf(abs(u0))
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])
    # cost_ate = (u0.nbits()-1)*(double_line_ate +sk) - sk + (u0.nbits()-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate+sparse_sparse+mk)
    # Frobenius(Q) : 2 sparse Frobenius in Fp16 (actually cheaper than 2*f16)
    additional_terms = 2*f16 + add_line_ate_cln + double_line_ate_cln + 2*m16 + f16
    cost_opt_ate_cln = (u0.nbits()-1)*(double_line_ate_cln +sk) -sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate_cln+sparse_sparse+mk) + additional_terms
    cost_opt_ate_cln_2naf = (len(bits_2naf_u0)-1)*(double_line_ate_cln +sk) -sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate_cln+sparse_sparse+mk) + additional_terms
    min_cost_miller = min(cost_opt_ate_cln, cost_opt_ate_cln_2naf)
    print("m16 = {}m s16 = {}m".format(m16,s16))
    print("m4 = {}m s4 = {}m".format(m4,s4))
    #print("cost ate Miller = {}m".format(cost_ate))
    print("cost opt ate Miller = {}m (Costello Lange Naehrig)".format(cost_opt_ate_cln))
    print("({0}-1)*(2*m{2}+8*s{2}+(k//2)*m +sk) -sk + ({0}-2-({1}-1))*(8*m{2}) + ({1}-1)*(9*m{2}+5*s{2}+(k//2)*m +6*m{2}+mk)".format(u0.nbits(),hw_u,e))
    print("cost opt ate Miller = {}m (2-naf, Costello Lange Naehrig)".format(cost_opt_ate_cln_2naf))
    print("({0}-1)*(2*m{2}+8*s{2}+(k//2)*m +sk) -sk + ({0}-2-({1}-1))*(8*m{2}) + ({1}-1)*(9*m{2}+5*s{2}+(k//2)*m +6*m{2}+mk)".format(len(bits_2naf_u0),hw2naf_u,e))

    # final exponentiation
    # (p^16-1)/r = (p^16-1)/Phi_16(u) * Phi_16(u)/r = (p^8-1)*(p^8+1)/r
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    assert ((x**16-1) // cyclotomic_polynomial(16)) == (x**8-1)

    px = (x**10 + 2*x**9 + 5*x**8 + 48*x**6 + 152*x**5 + 240*x**4 + 625*x**2 + 2398*x + 3125)/980
    rx = (x**8 + 48*x**4 + 625)/61250 # 625 = 5^4, 61250 = 2*5^4*7^2
    tx = (2*x**5 + 41*x + 35)/35
    cx = 125 * (x**2 + 2*x + 5)/2 # C such that P+1-T = C*R
    yx = (x**5 + 5*x**4 + 38*x + 120)/70 # Y such that T^2 - 4*P = -4*Y^2
    betax = (x**9-11*x**8-16*x**7-120*x**6-32*x**5-498*x**4-748*x**3-2740*x**2-3115*x-5651)/4018
    lambx = (x**4 + 24)/7 # sqrt(-1) mod R
    k = 16
    D = 1
    exponent = (px**8+1)//rx
    # easy part: px**8 - 1
    easy_part = i16 + m16
    # hard part: (px**8+1) // rx
    # exponentiation to the power u
    exp_u = (u0.nbits()-1)*s16_cyclo + (hw_u-1)*m16
    # exponentiation to the power u+1
    hw_u1 = sum((abs(u0+1)).digits(2))
    exp_u1 = ((u0+1).nbits()-1)*s16_cyclo + (hw_u1-1)*m16
    hard_part = 7*exp_u + 2*exp_u1 + 34*s16_cyclo + 32*m16 + 7*f16 + 3*c16_cyclo

    exp_u_2naf = (len(bits_2naf_u0)-1)*s16_cyclo + (hw2naf_u-1)*m16
    bits_2naf_u01 = bits_2naf(abs(u0+1))
    hw2naf_u1 = sum([1 for bi in bits_2naf_u01 if bi != 0])
    exp_u1_2naf = (len(bits_2naf_u01)-1)*s16_cyclo + (hw2naf_u1-1)*m16
    hard_part_2naf = 7*exp_u_2naf + 2*exp_u1_2naf + 34*s16_cyclo + 32*m16 + 7*f16 + 3*c16_cyclo

    print("cost final exp easy:  {}m".format(easy_part))
    print("cost final exp hard:  {}m (with cyclotomic squarings)".format(hard_part))
    print("cost final exp hard:  {}m (2-naf with cyclotomic squarings)".format(hard_part_2naf))
    print("cost final exp:       {}m".format(easy_part + hard_part))
    print("cost final exp 2-naf: {}m".format(easy_part + hard_part_2naf))
    print("\ncost pairing (total)  {}m".format(easy_part + min(hard_part, hard_part_2naf) + min_cost_miller))

    print("new formula")
    # new formula     # total cost 33 S + 31 M + 2 exp(u+1) + 7 exp(u) + 14*f16 + 8 inv_frob_8
    hard_part = 7*exp_u + 2*exp_u1 + 33*s16_cyclo + 31*m16 + 14*f16
    hard_part_2naf = 7*exp_u_2naf + 2*exp_u1_2naf + 33*s16_cyclo + 31*m16 + 14*f16
    print("cost final exp hard:  {}m (with cyclotomic squarings)".format(hard_part))
    print("cost final exp hard:  {}m (2-naf with cyclotomic squarings)".format(hard_part_2naf))
    print("cost final exp:       {}m".format(easy_part + hard_part))
    print("cost final exp 2-naf: {}m".format(easy_part + hard_part_2naf))
    print("\ncost pairing (total)  {}m".format(easy_part + min(hard_part, hard_part_2naf) + min_cost_miller))

def cost_pairing_kss16_330():
    u0 = ZZ(-2**34 + 2**27 - 2**23 + 2**20 - 2**11 + 1)
    cost_pairing_kss16(u0)

def cost_pairing_kss16_330b():
    u0 = ZZ(2**34 - 2**30 + 2**26 + 2**23 + 2**14 - 2**5 + 1)
    cost_pairing_kss16(u0)

def cost_pairing_kss16_339():
    u0 = ZZ(2**35 - 2**32 - 2**18 + 2**8 + 1)
    cost_pairing_kss16(u0)

def cost_pairing_kss18(u0):
    # extension: Fp - Fp3 - Fp9 - Fp18
    # a=0
    # AKLGL'11
    # there are three formulas in https://www.lirmm.fr/arith18/papers/Chung-Squaring.pdf for squaring in cubic extensions
    # S3 = m + 4*s; S3 = 2*m+3*s; S3 = 3*m+2*s
    u0 = ZZ(abs(u0))
    m = 1
    s = 1
    m3 = 6*m
    s3 = 2*m+3*s
    m9 = 6*m3
    s9 = 2*m3+3*s3
    m18 = 3*m9
    s18 = 2*m9
    s18_cyclo = 6*m3
    f18 = 16*m # to check
    inv = 25*m
    i3 = 9*m + 3*s + inv
    i9 = 9*m3 + 3*s3 + i3
    i6 = 2*m3 + 2*s3 + i3
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2
    s6 = 2*m2+3*s2

    i18 = min(2*m9 + 2*s9 + i9, 9*m6 + 3*s6 + i6) # 213+inv, 232 + inv
    i18_cyclo = 0 # Frobenius power p^9 is only conjugation
    k = 18

    # AKLGL
    # k=18, d=6, e=k/d = 3
    k = 18
    d = 6
    k_d = k//d
    me=m3 ; se=s3; mk=m18; sk=s18 # e = k/d
    double_line_ate = 3*me+6*se+(k//3)*m
    add_line_ate    = 11*me+2*se+(k//3)*m
    add_line_ate_with_z = 11*me+2*se+(k//3)*m + 5*me # +5*me seems the upper bound
    # Costello Lange Naehrig
    # double_line_ate_cln = 3*me+5*se+(k//3)*m    ???
    double_line_ate_cln = 2*me+7*se+(k//3)*m
    add_line_ate_cln    = 10*me+2*se+(k//3)*m
    add_line_ate_with_z_cln = 14*me+2*se+(k//3)*m
    sparse_dense  = 13*me
    sparse_sparse = 6*me
    update1       = 13*me+sk
    update2       = 13*me

    hw_u = sum((abs(u0)).digits(2))
    bits_2naf_u0 = bits_2naf(abs(u0))
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])

    additional_terms = double_line_ate + add_line_ate + add_line_ate_with_z + 5*f18 + 3*m18
    additional_terms_cln = double_line_ate_cln + add_line_ate_cln + add_line_ate_with_z_cln + 5*f18 + 3*m18

    cost_ate = (u0.nbits()-1)*(double_line_ate +sk) -sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate+sparse_sparse+mk) + additional_terms
    cost_ate_2naf = (len(bits_2naf_u0)-1)*(double_line_ate +sk) -sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate+sparse_sparse+mk) + additional_terms
    cost_ate_cln = (u0.nbits()-1)*(double_line_ate_cln +sk) -sk + (u0.nbits()-2-(hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate_cln+sparse_sparse+mk) + additional_terms_cln
    cost_ate_cln_2naf = (len(bits_2naf_u0)-1)*(double_line_ate_cln +sk) -sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate_cln+sparse_sparse+mk) + additional_terms_cln
    print("m18 = {}m s18 = {}m".format(m18,s18))
    print("m9 = {}m s9 = {}m".format(m9,s9))
    print("m3 = {}m s3 = {}m".format(m3,s3))
    print("cost ate Miller       = {}m".format(cost_ate))
    print("cost ate Miller 2-naf = {}m".format(cost_ate_2naf))

    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(u0.nbits(), hw_u, k_d))
    print("({0}-1)*(3*m{2}+6*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(11*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(len(bits_2naf_u0), hw2naf_u, k_d))
    print("cost ate Miller       = {}m (Costello Lange Naehrig)".format(cost_ate_cln))
    print("cost ate Miller 2-naf = {}m (Costello Lange Naehrig)".format(cost_ate_cln_2naf))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(u0.nbits(), hw_u, k_d))
    print("({0}-1)*(2*m{2}+7*s{2}+(k//3)*m +sk) -sk + ({0}-2-({1}-1))*(13*m{2}) + ({1}-1)*(10*m{2}+2*s{2}+(k//3)*m +6*m{2}+mk)".format(len(bits_2naf_u0), hw2naf_u, k_d))

    min_cost_miller = min(cost_ate, cost_ate_2naf, cost_ate_cln, cost_ate_cln_2naf)

    # final exponentiation
    # (p^18-1)/r = (p^18-1)/Phi_18(u) * Phi_18(u)/r = (p^9-1)*(p^3 + 1) * (p^6-p^3+1)/r
    QQx = QQ['x']; (x,) = QQx._first_ngens(1)
    assert ((x**18-1) // cyclotomic_polynomial(18)) == (x**9-1)*(x**3+1)
    px = (x**8 + 5*x**7 + 7*x**6 + 37*x**5 + 188*x**4 + 259*x**3 + 343*x**2 + 1763*x + 2401)/21
    rx = (x**6 + 37*x**3 + 343)/343 # 343 = 7^3
    # Algorithm 7.6 in Guide to Pairing-based cryptography: 6 S + 53 M + 29 F + 7 exp(u) + 8 F9
    # Variant a la https://eprint.iacr.org/2020/875:       19 S + 26 M + 10 F + 7 exp(u) + 6 F9
    # New variant https://eprint.iacr.org/2021/1309:       11 S + 24 M +  5 F + 7 exp(u) + 7 F9

    # easy part: (px**9 - 1)*(p**3+1)
    easy_part = 2*f18 + i18 + 2*m18
    # exponentiation to the power u
    exp_u = (u0.nbits()-1)*s18_cyclo + (hw_u-1)*m18
    exp_u_2naf = (len(bits_2naf_u0)-1)*s18_cyclo + (hw2naf_u-1)*m18

    hard_part = 7*exp_u + 6*s18_cyclo + 53*m18 + 29*f18 + 8 * i18_cyclo
    hard_part_2naf = 7*exp_u_2naf + 6*s18_cyclo + 53*m18 + 29*f18 + 8 * i18_cyclo

    new_hard_part = 7*exp_u + 19*s18_cyclo + 26*m18 + 10*f18 + 6 * i18_cyclo
    new_hard_part_2naf = 7*exp_u_2naf + 19*s18_cyclo + 26*m18 + 10*f18 + 6 * i18_cyclo
    #Cost: 7 exp(u) + 19 S + 26 M + 10 Frobenius powers + 6 Inversions with Frobenius power p^9

    hard_part_2021_1309 = 7*exp_u + 11*s18_cyclo + 24*m18 + 5*f18 + 7 * i18_cyclo
    hard_part_2naf_2021_1309 = 7*exp_u_2naf + 11*s18_cyclo + 24*m18 + 5*f18 + 7 * i18_cyclo

    print("cost final exp easy:  {}m".format(easy_part))
    print("cost final exp hard:  {}m (with cyclotomic squarings)".format(hard_part))
    print("cost final exp hard:  {}m (2-naf with cyclotomic squarings)".format(hard_part_2naf))
    print("cost final exp:       {}m".format(easy_part + hard_part))
    print("cost final exp 2-naf: {}m".format(easy_part + hard_part_2naf))
    print("\ncost pairing (total)  {}m".format(easy_part + min(hard_part, hard_part_2naf) + min_cost_miller))

    print("\nnew formula a la ePrint 2020/875:")
    print("cost final exp hard:  {}m (with cyclotomic squarings)".format(new_hard_part))
    print("cost final exp hard:  {}m (2-naf with cyclotomic squarings)".format(new_hard_part_2naf))
    print("cost final exp:       {}m".format(easy_part + new_hard_part))
    print("cost final exp 2-naf: {}m".format(easy_part + new_hard_part_2naf))
    print("\ncost pairing (total)  {}m".format(easy_part + min(new_hard_part, new_hard_part_2naf) + min_cost_miller))

    print("\nnew formula eprint 2021/1309:")
    print("cost final exp hard:  {}m (with cyclotomic squarings)".format(hard_part_2021_1309))
    print("cost final exp hard:  {}m (2-naf with cyclotomic squarings)".format(hard_part_2naf_2021_1309))
    print("cost final exp:       {}m".format(easy_part + hard_part_2021_1309))
    print("cost final exp 2-naf: {}m".format(easy_part + hard_part_2naf_2021_1309))
    print("\ncost pairing (total)  {}m".format(easy_part + min(hard_part_2021_1309, hard_part_2naf_2021_1309) + min_cost_miller))

def cost_pairing_fst64(u0, k):
    """
    The arithmetic cost over GF(p^5) is assumed with the formulas of
    Five, Six, and Seven-Term Karatsuba-Like Formulae, Peter L. Montgomery,
    IEEE Transactions on Computers, vol. 54, no. 3, March 2005
    M5 = 13*m,
    M7 = 22*m, and because the formulas are symmetric in the variables a_i, b_i,
    S5 = 13*s
    S7 = 22*s
    """
    assert (k % 8) == 4
    assert k in [20, 28]
    # cost
    # b=0
    #k = 20
    d = 4
    e = k // d
    u0 = ZZ(u0)
    m = 1
    s = 1
    inv = 25*m
    #
    # Granger Scott, Faster Squaring in the Cyclotomic Subgroup of Sixth Degree Extensions,
    # PKC 2010, pages 209-223.
    # http://www.iacr.org/archive/pkc2010/60560212/60560212.pdf
    if k == 20:
        m5 = 13*m # Montgomery formula
        s5 = 13*s # Montgomery formula
        f5 = 4*m
        m10 = 3*m5 # Karastuba in quadratic extensions
        s10 = 2*m5 # squaring in quadratic extensions
        m20 = 3*m10 # Karastuba in quadratic extensions
        s20 = 2*m10 # squaring in quadratic extensions
        s20_cyclo = 2*s10 # = 4*m5
        f20 = (k-2)*m
        i5 = 48*m + inv # Masson PhD thesis page ix
        i10 = 2*m5 + 2*s5 + i5
        i20 = 2*m10 + 2*s10 + i10
        me=m5; se=s5; mk = m20; sk=s20; sk_cyclo = s20_cyclo; fk = f20; ik = i20 # e = k/d with d=4
    elif k == 28:
        m7 = 22*m
        s7 = 22*s
        f7 = 6*m
        m14 = 3*m7
        s14 = 2*m7
        m28 = 3*m14
        s28 = 2*m14
        s28_cyclo = 2*s14
        f28 = (k-2)*m
        i7 = 104*m + inv # Masson PhD thesis page ix
        i14 = 2*m7 + 2*s7 + i7
        i28 = 2*m14 + 2*s14 + i14
        me=m7; se=s7; mk = m28; sk=s28; sk_cyclo = s28_cyclo; fk = f28; ik = i28 # e = k/d with d=4
    # Costello Lange Naehrig, with d=4, b=0
    # Table 6 and 7 page 20 of Eurocrypt'22, El Housni Guillevic, eprint.iacr.org/2021/1359
    double_line_ate_cln = 2*me+8*se+(k//2)*m
    add_line_ate_cln    = 9*me+5*se+(k//2)*m
    sparse_dense = 8*me
    sparse_sparse = 6*me

    hw_u = sum((abs(u0)).digits(2))
    bits_2naf_u0 = bits_2naf(abs(u0))
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])

    cost_opt_ate_cln = (u0.nbits()-1)*(double_line_ate_cln + sk) - sk + (u0.nbits()-2 - (hw_u-1))*sparse_dense + (hw_u-1)*(add_line_ate_cln + sparse_sparse + mk)
    cost_opt_ate_cln_2naf = (len(bits_2naf_u0)-1)*(double_line_ate_cln + sk) - sk + (len(bits_2naf_u0)-2-(hw2naf_u-1))*sparse_dense + (hw2naf_u-1)*(add_line_ate_cln + sparse_sparse + mk)
    min_cost_miller = min(cost_opt_ate_cln, cost_opt_ate_cln_2naf)

    print("m{} = {}m s{} = {}m".format(k, mk, k, sk))
    print("m{} = {}m s{} = {}m".format(e, me, e, se))
    #print("cost ate Miller = {}m".format(cost_ate))
    print("cost opt ate Miller = {}m (Costello Lange Naehrig)".format(cost_opt_ate_cln))
    print("({0}-1)*(2*m{2}+8*s{2}+(k//2)*m + sk) -sk + ({0}-2-({1}-1))*(8*m{2}) + ({1}-1)*(9*m{2}+5*s{2}+(k//2)*m + 6*m{2}+mk)".format(u0.nbits(), hw_u, e))
    print("cost opt ate Miller = {}m (2-naf, Costello Lange Naehrig)".format(cost_opt_ate_cln_2naf))
    print("({0}-1)*(2*m{2}+8*s{2}+(k//2)*m + sk) -sk + ({0}-2-({1}-1))*(8*m{2}) + ({1}-1)*(9*m{2}+5*s{2}+(k//2)*m + 6*m{2}+mk)".format(len(bits_2naf_u0), hw2naf_u, e))

    # final exponentiation
    # k=20 (p^20-1)/r = (p^20-1)/Phi_20(u) * Phi_20(u)/r = (p^10-1)*(p^2+1)*(p^8-p^6+p^4-p^2+1)/r
    # k=28 (p^28-1)/r = (p^28-1)/Phi_28(u) * Phi_28(u)/r = (p^14-1)*(p^2+1)*(p^12-p^10+p^8-p^6+p^4-p^2+1)/r
    # easy part: m^((p^(k//2-1))*(p^2+1)) where m^(k//2) is a conjugation and costs nothing
    easy_part = ik + mk + fk + mk
    # exponentiation to the power u
    exp_u = (u0.nbits()-1)*sk_cyclo + (hw_u-1)*mk
    # exponentiation to the power (u-1)//2
    u1 = abs(u0-1)//2
    hw_u1 = sum(u1.digits(2))
    exp_u1 = (u1.nbits()-1)*sk_cyclo + (hw_u1-1)*mk
    # 2-naf
    exp_u_2naf = (len(bits_2naf_u0)-1)*sk_cyclo + (hw2naf_u-1)*mk
    bits_2naf_u1 = bits_2naf(u1)
    hw2naf_u1 = sum([1 for bi in bits_2naf_u1 if bi != 0])
    exp_u1_2naf = (len(bits_2naf_u1)-1)*sk_cyclo + (hw2naf_u1-1)*mk

    if k == 20:
        #total cost 2*exp((u-1)//2) + 9*exp(|u|) + 8 M + frob + 2*frob(2) + frob(4) + 2 conj + 1 conj if u < 0
        hard_part = 2*exp_u1 + 9*exp_u + 8*mk + fk + 2*fk + fk
        hard_part_2naf = 2*exp_u1_2naf + 9*exp_u_2naf + 8*mk + fk + 2*fk + fk
    elif k == 28:
        #total cost 2*exp((u-1)//2) + 12*exp(|u|) + 13 M + frob + 5*frob(2) + conj (+ 1 conj if u < 0)
        hard_part = 2*exp_u1 + 13*exp_u + 12*mk + fk + 5*fk
        hard_part_2naf = 2*exp_u1_2naf + 13*exp_u_2naf + 12*mk + fk + 5*fk

    print("cost final exp easy:  {}m".format(easy_part))
    print("cost final exp hard:  {}m (with cyclotomic squarings)".format(hard_part))
    print("cost final exp hard:  {}m (2-naf with cyclotomic squarings)".format(hard_part_2naf))
    print("cost final exp:       {}m".format(easy_part + hard_part))
    print("cost final exp 2-naf: {}m".format(easy_part + hard_part_2naf))
    print("\ncost pairing (total)  {}m".format(easy_part + min(hard_part, hard_part_2naf) + min_cost_miller))

def cost_pairing_fst64_k20(u0):
    cost_pairing_fst64(u0, 20)

def cost_pairing_fst64_k28(u0):
    cost_pairing_fst64(u0, 28)

def cost_miller_loop_bls_k_odd(u, k, m, me, se, mk, sk, ik_cyclo):
    """
    The Miller loop for optimal ate pairing is f_{u,Q}(P)
    but there is no denominator elimination, hence
    the lines and tangents simplify less
    Formulas from Costello Lange Naehrig PKC 2010, Section 6
    Dbl+line evaluation: 6*me + 7*se + k*m
    mxAdd+line evaluation: 13*me + 3*se + k*m
    Add+line evaluation: 16*me + 3*se + k*m
    """
    d = 3
    assert (k % d) == 0
    e = k//d
    u0 = ZZ(u)
    double_line_ate_cln = 6*me + 7*se + k*m
    mx_add_line_ate_cln = 13*me + 3*se + k*m

    hw_u = sum((abs(u0)).digits(2))
    bits_2naf_u0 = bits_2naf(abs(u0))
    hw2naf_u = sum([1 for bi in bits_2naf_u0 if bi != 0])
    hw_neg_u = sum([1 for bi in bits_2naf_u0 if bi == -1])

    cost_ate_cln = (u0.nbits()-1)*(double_line_ate_cln+sk+mk) -sk -mk + (hw_u-1)*(mx_add_line_ate_cln+mk)
    cost_ate_cln_2naf = (len(bits_2naf_u0)-1)*(double_line_ate_cln+sk+mk) -sk -mk + hw_neg_u*ik_cyclo + (hw2naf_u-1)*(mx_add_line_ate_cln+mk)
    # I think there is another extra cost when the bit is -1
    print("m{} = {}m s{} = {}m".format(k, mk, k, sk))
    print("m{} = {}m s{} = {}m".format(e, me, e, se))
    print("cost ate Miller       = {}m (Costello Lange Naehrig)".format(cost_ate_cln))
    print("cost ate Miller 2-naf = {}m (Costello Lange Naehrig)".format(cost_ate_cln_2naf))
    print("({0}-1)*(6*m{2}+7*s{2}+k*m +sk+mk) -sk-mk + ({1}-1)*(13*m{2}+3*s{2}+k*m +mk)".format(u0.nbits(), hw_u, e))
    print("({0}-1)*(6*m{2}+7*s{2}+k*m +sk+mk) -sk-mk + ({1}-1)*(13*m{2}+3*s{2}+k*m +mk) + ({3})*({4})".format(len(bits_2naf_u0), hw2naf_u, e, hw_neg_u, ik_cyclo))
    min_cost_miller = min(cost_ate_cln, cost_ate_cln_2naf)
    return min_cost_miller

def cost_pairing_bls21(u):
    """
    The k=21 extension is built as a degree 3 extension on top of a degree 7 extension
    The arithmetic cost over GF(p^7) is assumed with the formulas of
    Five, Six, and Seven-Term Karatsuba-Like Formulae, Peter L. Montgomery,
    IEEE Transactions on Computers, vol. 54, no. 3, March 2005
    M7 = 22*m and because the formulas are symmetric in the variables a_i, b_i,
    S7 = 22*s
    """
    # cost
    # a=0 E: y^2 = x^3 + b
    k = 21
    d = 3
    e = k // d
    u0 = ZZ(u)
    m = 1
    s = 1
    inv = 25*m
    # no cyclotomic squaring
    m7 = 22*m
    s7 = 22*s
    f7 = 6*m
    m21 = 6*m7
    s21 = 2*m7 + 3*s7
    s21_cyclo = s21 # no better formula known
    f21 = (k-1)*m
    # f^q, f^(q^2) where q=p^7, in GF(p^21)
    # there is a better way than (k-1)*m
    # f = f0 + f1*alpha + f2*alpha^2 where fi in GF(p^7)
    # f^q = f0^q + f1^q*alpha^q + f2^q*alpha^2^q but fi^q = fi as fi in GF(p^7) = GF(q)
    # f^q = f0 + f1*alpha^q + f2*alpha^2^q
    # alpha^3 = residue in GF(p^7)
    # (alpha^3)^((q-1)/3) * alpha = residue^((q-1)/3) * alpha
    # = residue^((p-1)/3*(1+p^2+p^3+p^4+p^5+p^6)) * alpha
    # = (Norm_{GF(p^7)/GF(p)}(residue))^((p-1)/3) * alpha
    # = omega * alpha where omega is one of the two primitive 3-rd roots of unity in GF(p)
    #                       as p = 1 mod 3
    f21_14 = 14*m
    f21_7 = 14*m
    i7 = 104*m + inv # Masson PhD thesis page ix
    i21 = 3*m7 + 3*s7 + i7
    i21_cyclo = f21_14 + f21_7 + m21 # inversion in cyclotomic subgroup
    me=m7; se=s7; mk = m21; sk=s21; sk_cyclo = s21_cyclo; fk = f21; ik = i21; ik_cyclo = i21_cyclo # e = k/d with d=3
    cost_miller_loop_bls_k_odd(u0, k, m, me, se, mk, sk, ik_cyclo)
    # final exp:
    # total cost exp((u-1)/3) + exp(u^3-1) + exp(u-1) + 10 exp(u) + 18 M + 11 f + inv_cyclo_k21

def cost_pairing_bls27(u):
    """
    The k=27 extension is built as three degree 3 extensions on top of GF(p)
    """
    # cost
    # a=0 E: y^2 = x^3 + b
    k = 27
    d = 3
    e = k // d
    u0 = ZZ(u)
    m = 1
    s = 1
    inv = 25*m
    m3 = 6*m
    s3 = 2*m+3*s
    m9 = 6*m3
    s9 = 2*m3+3*s3
    m27 = 6*m9
    s27 = 2*m9+3*s9
    s27_cyclo = s27 # no formula
    i3 = 9*m + 3*s + inv
    i9 = 9*m3 + 3*s3 + i3
    i27 = 9*m9 + 3*s9 + i9
    f27 = (k-1)*m
    f27_9 = 2*(k//3)*m
    f27_18 = 2*(k//3)*m
    i27_cyclo = f27_18 + f27_9 + m27 # inversion in cyclotomic subgroup
    me=m9; se=s9; mk=m27; sk=s27; sk_cyclo=s27_cyclo; fk=f27; ik=i27; ik_cyclo=i27_cyclo
    cost_miller_loop_bls_k_odd(u0, k, m, me, se, mk, sk, ik_cyclo)

def cost_pairing_bw6_bls(u, k_bls, t_mod_r_mod_u=0):
    # cost
    # a=0
    # AKLGL'11
    m = 1
    s = 1
    m3 = 6*m
    s3 = 5*m
    m6 = 18*m
    s6 = 12*m
    s2 = 2*m
    s6_cyclo = 3*s2
    f6 = 4*m
    f3 = 0 # frobenius to the power p^3 costs 3 negations
    inv = 25*m
    i6 = 34*m+inv
    k = 6
    DoubleLine_ate = 3*m+6*s+(k//3)*m
    AddLine_ate    = 11*m+2*s+(k//3)*m
    Update1        = 13*m+s6
    Update2        = 13*m
    sparse_mult    = 13*m
    sparse_sparse_mult = 6*m

    u = ZZ(abs(u))
    # Miller optimal ate pairing formula is
    # tr = 0 mod r mod u:
    #f_{x+1, -Q}(P) f_{x^5-x^4+1, pi(Q)}(P)
    # tr = 3 mod r mod u:
    #f_{x+1, Q}(P) f_{x^5-x^4-x, pi(Q)}(P)
    # note that sparse_sparse_mult + m6 = 6*m + 18*m = 24*m,
    # while 2*(sparse_mult) = 2*13*m = 26*m
    # 1. f_u:
    #cost_ate1 = (len(bits_2naf_u)-1)*DoubleLine_ate + (len(bits_2naf_u)-2)*s6 + (HW2naf_u-1)*(AddLine_ate+sparse_sparse_mult+m6) + (len(bits_2naf_u)-2 - (HW2naf_u-1))*sparse_mult
    # 2. (f_u)_{x^4-x^3-1}
    if k_bls == 24:
        if t_mod_r_mod_u==0:
            print("formula with (u^5-u^4-u = u*(u^4-u^3-1) = (u+1)*(u^4-2u^3+2u^2-2u+1)-1, u+1)")
            print("formula with (-(u+1), u^5-u^4+1 = u*(u^4-u^3)+1 = (u+1)*(u^4-2u^3+2u^2-2u+2)-1)")
        else:
            print("formula with (u+1, u^5-u^4-u = u*(u^4-u^3-1) = (u+1)*(u^4-2u^3+2u^2-2u+1)-1")
            print("formula with (u^5-u^4+1 = u*(u^4-u^3)+1 = (u+1)*(u^4-2u^3+2u^2-2u+2)-1, -u-1)")
        u1 = u
        u1_str = "u"
        v1 = u**4-u**3-1
        v1_str = "u^4-u^3-1"
        u2 = u+1
        u2_str = "u+1"
        v2 = u**4-2*u**3+2*u**2-2*u+1
        v2_str = "u^4-2u^3+2u^2-2u+1"
        u3 = u
        u3_str = "u"
        v3 = u**4-u**3
        v3_str = "u^4-u^3"
        u4 = u+1
        u4_str = "u+1"
        v4 = u**4-2*u**3+2*u**2-2*u+2
        v4_str = "u^4-2u^3+2u^2-2u+2"
    else:
        if t_mod_r_mod_u==0:
            print("formula with (u^3-u^2-u = u*(u^2-u-1) = (u+1)*(u^2-2u+1)-1, u+1)")
            print("alternative formula with (-(u+1), u^3-u^2+1 = u*(u^2-u) + 1 = (u+1)*(u^2-2u+2)-1)")
        else:
            print("formula with (u+1, u^3-u^2-u = u*(u^2-u-1) = (u+1)*(u^2-2u+1)-1")
            print("alternative formula with (u^3-u^2+1 = u*(u^2-u) + 1 = (u+1)*(u^2-2u+2)-1, -(u+1))")
        u1 = u
        u1_str = "u"
        v1 = u**2-u-1
        v1_str = "u^2-u-1"
        u2 = u+1
        u2_str = "u+1"
        v2 = u**2-2*u+1
        v2_str = "u^2-2u+1"
        u3 = u
        u3_str = "u"
        v3 = u**2-u
        v3_str = "u^2-u"
        u4 = u+1
        u4_str = "u+1"
        v4 = u**2-2*u+2
        v4_str = "u^2-2u+2"
    for (uu, vv, uu_str, vv_str) in [(u1, v1, u1_str, v1_str), (u2, v2, u2_str, v2_str), (u3, v3, u3_str, v3_str), (u4, v4, u4_str, v4_str)]:
        bits_2naf_uu = bits_2naf(uu)
        HW2naf_uu = sum([1 for bi in bits_2naf_uu if bi != 0])
        cost_ate1 = (len(bits_2naf_uu)-1)*DoubleLine_ate + (len(bits_2naf_uu)-2)*s6 + (HW2naf_uu-1)*(AddLine_ate+sparse_sparse_mult+m6) + (len(bits_2naf_uu)-2 - (HW2naf_uu-1))*sparse_mult
        bits_2naf_v = bits_2naf(vv)
        HW2naf_v = sum([1 for bi in bits_2naf_v if bi != 0])
        cost_ate2 = (len(bits_2naf_v)-1)*(DoubleLine_ate + s6) + (HW2naf_v-1)*(AddLine_ate+sparse_sparse_mult+m6+m6) + (len(bits_2naf_v)-1 - (HW2naf_v-1))*sparse_mult
        # no inversion required anymore
        cost_ate = cost_ate1 + inv + 2*m + AddLine_ate + Update2 + cost_ate2 + f6 + m6

        print("cost ate Miller 2naf = {}m".format(cost_ate))
        print("{0:20s}, len(bits_2naf({0:3})) = {1:3}, Hw2naf({0:3}) = {2:2}, (Hw2naf({0:3})-1) * m6 = {3}".format(uu_str, len(bits_2naf_uu), HW2naf_uu, (HW2naf_uu-1)*m6))
        print("v={:18s}, len(bits_2naf(v  )) = {:3}, Hw2naf(v  ) = {:2}, (Hw2naf(v  )-1) * m6 = {}".format(vv_str, len(bits_2naf_v), HW2naf_v, (HW2naf_v-1)*m6))

    # alternative formula:
    # t_bw0
    #-u-1 + (u^3-u^2+1)*q -> fx = f_{u+1}, (fx)_{u^2-2*u-2} - 1 or fx = f_u, (fx)_{u^2-u} + 1
    # t_bw3
    # u^3-u^2+1 - (u+1)*q
    print("add step costs extra AddLine_ate+sparse_sparse_mult+m6+m6-sparse_mult = {}m".format(AddLine_ate+sparse_sparse_mult+m6+m6-sparse_mult))
    print("add step costs extra AddLine_ate+sparse_sparse_mult+m6-sparse_mult = {}m".format(AddLine_ate+sparse_sparse_mult+m6-sparse_mult))

def cost_pairing_bw6_bls12(u, t_mod_r_mod_u=0):
    return cost_pairing_bw6_bls(u, 12, t_mod_r_mod_u=t_mod_r_mod_u)
def cost_pairing_bw6_bls24(u, t_mod_r_mod_u=0):
    return cost_pairing_bw6_bls(u, 24, t_mod_r_mod_u=t_mod_r_mod_u)

def cost_miller_optimal_twisted_ate_bw6(b0, b1, sk, mk, m, double_line, add_line, add_line_affine, sparse_sparse_mult, sparse_mult):
    """ algorithm 2 to compute optimal twisted ate on BW6 curves
    """
    list_b = [ZZ(b0), ZZ(b1)]
    # optimal Tate pairing with scalars b0, b1
    abs_b = [abs(bi) for bi in list_b]
    abs_b.sort(reverse = True)
    digits_b = [bi.digits(2) for bi in abs_b]
    len_bi = [len(di) for di in digits_b]
    ll0, ll1 = len_bi
    max_len = max(len_bi)
    digits_b = [digits_b[i] + [0]*(max_len - len(digits_b[i])) for i in range(2)] # make all the same length
    hw_b_idx_i = [sum([digits_b[j][i] for j in range(2)]) for i in range(max_len)]
    hw_b = [sum(digits_bi) for digits_bi in digits_b]

    bits_2naf_b = [bits_2naf(abs(bi)) for bi in list_b]
    len_bi_2naf = [len(bits_2naf_b[i]) for i in range(2)]
    max_len_2naf = max(len_bi_2naf)
    bits_2naf_b = [bits_2naf_b[i] + [0]*(max_len - len(bits_2naf_b[i])) for i in range(2)] # make all the same length

    hw2naf_b = [sum([1 for bi in bits_2naf_bi if bi != 0]) for bits_2naf_bi in bits_2naf_b]
    hw_b_idx_i_2naf = [sum([abs(bits_2naf_b[j][i]) for j in range(2)]) for i in range(max_len)]

    cost_miller_f = 0
    for i in range(max_len-2, -1, -1):
        cost_miller_f += sk + double_line
        if hw_b_idx_i[i] > 0:
            cost_miller_f += add_line + sparse_sparse_mult + mk
            if hw_b_idx_i[i] == 2:
                cost_miller_f += sparse_mult
        else:
            cost_miller_f += sparse_mult

    cost_miller_f_2naf = 0
    additional_terms = m + 2*add_line_affine # compute psi(P), P+psi(P), P-psi(P) and the line through P, psi(P) evaluated at Q, the line through P, P-psi(P) evaluated at Q.
    for i in range(max_len_2naf-2, -1, -1):
        cost_miller_f_2naf += sk + double_line
        if hw_b_idx_i_2naf[i] > 0:
            cost_miller_f_2naf += add_line + sparse_sparse_mult + mk
            if hw_b_idx_i_2naf[i] == 2:
                cost_miller_f_2naf += sparse_mult
        else:
            cost_miller_f_2naf += sparse_mult
    return cost_miller_f, cost_miller_f_2naf

def cost_pairing_opt_tate_bw6_bls(u, bls_k, trace_mod_r_mod_u=0, k=6):
    m = 1
    s = 1
    inv = 25*m
    if k == 6 or k == 12:
        i2 = 2*m +2*s + inv
        i6 = 34*m + inv
        m2 = 3*m
        s2 = 2*m
        m6 = 6*m2 # multiplication in cubic extension
        #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
        s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
        s6_cyclo = 3*s2 # sk_cyclo = 3*s_{k/3} for 6 | k
        f6 = 4*m
        if k == 6:
            # CLN PKC 2010 section 5
            sk = s6; mk = m6; sk_cyclo = s6_cyclo; ik = i6; fk = f6
            me = m; se = s; ie = inv
        elif k == 12:
            m12 = 3*m6
            s12 = 2*m6
            f12 = 10*m
            f12_p2 = 8*m # p^2 power in Fp12
            s12_cyclo = 18*m # 3*s4
            i12 = 97*m+inv # 94*m + i in Guillevic Masson Thome DCC'20 is not compatible with the tower
            sk = s12; mk = m12; sk_cyclo = s12_cyclo; ik = i12; fk = f12
            me = m; se = s; ie = inv
        # for add_line_affine for k = 0 mod 6, take INDOCRYPT 2012 paper by Zhang and Lin
        add_line_affine = inv + 3*m + 2*s + (k//3)*m
        double_line_affine = inv + 3*m + 2*s + (k//3)*m
        add_line_tate_cln = 10*m + 2*s + (k//3)*m
        double_line_tate_cln = 2*m + 7*s + (k//3)*m
        add_line_tate_aklgl = 11*m + 2*s + (k//3)*m
        double_line_tate_aklgl = 3*m + 6*s + (k//3)*m
        if k == 6:
            sparse_sparse_line_mult = 6*m # aklgl
            sparse_line_mult = 13*m # aklgl
        else:
            sparse_sparse_line_mult = 5*me + m # aklgl
            sparse_line_mult = 10*me + 3*(k//6)*m # aklgl
        update1 = sparse_line_mult
        additional_terms = m + add_line_affine # compute psi(P), P+psi(P) and the line through P, psi(P) evaluated at Q.
        add_line_tate = add_line_tate_aklgl
        double_line_tate = double_line_tate_aklgl

    # it is like r-order subgroup membership testing
    # with BL12:
    # there is one endomorphism on G1 of characteristic polynomial x^2+x+1 and eigenvalue beta/lambda mod r / mod p
    # 3*r*P = (x+1)*phi(P) + (-x^3 + x^2 + x)*P
    # Miller loop of optimal Tate has scalars (b0,b1) = (x+1, -x^3 + x^2 + x)
    if bls_k == 12:
        if trace_mod_r_mod_u == 0:
            # bw0
            b0 = -(u+1)
            b1 = (u**3-u**2+1)
            b_str = "(-u-1,u^3-u^2+1) = (-u-1, (u+1)(u^2-2u+2)-1" # x^3-x^2+1 == (x+1)*(x^2-2*x+2) - 1
            bb0 = u**3 - u**2 - u
            bb1 = u+1
            bb_str = "(u^3-u^2-u,u+1 = ((u+1)(u^2-2u+1)-1, u+1)"
        else:
             # bw3
            b0 = u**3-u**2+1
            b1 = -(u+1)
            b_str = "(u^3-u^2+1,-u-1) = ((u+1)(u^2-2u+2)-1, -u-1)"
            bb0 = u+1
            bb1 = u**3 - u**2 - u
            bb_str = "(u+1,u^3-u^2-u) = (u+1, (u+1)(u^2-2u+1)-1)"
    elif bls_k == 24:
        if trace_mod_r_mod_u == 0:
            # bw0
            b0 = -(u+1)
            b1 = u**5 - u**4 + 1
            b_str = "(-u-1,u^5-u^4+1)"
            bb0 = u**5 - u**4 - u
            bb1 = u+1
            bb_str = "(u^5-u^4-u,u+1)"
        else:
            # bw3
            b0 = u**5-u**4+1
            b1 = -(u+1)
            b_str = "(u^5-u^4+1,-u-1)"
            bb0 = u+1
            bb1 = u**5-u**4-u
            bb_str = "(u+1,u^5-u^4-u)"
    list_scalars = [(b0,b1,b_str), (bb0, bb1, bb_str)]
    print("optimal Tate pairing with one endomorphism")
    for (b0,b1,b_str) in list_scalars:
        print("with scalars {}".format(b_str))
        cost, cost_2naf = cost_miller_optimal_twisted_ate_bw6(b0, b1, sk, mk, m, double_line_tate, add_line_tate, add_line_affine, sparse_sparse_line_mult, sparse_line_mult)
        print("cost Miller with 2-multi-scalar:       {} m".format(cost))
        print("cost Miller with 2-multi-scalar, 2naf: {} m".format(cost_2naf))

def cost_pairing_opt_tate_bw6_bls12(u):
    cost_pairing_opt_tate_bw6_bls(u, 12)

def cost_pairing_opt_tate_bw6_bls24(u):
    cost_pairing_opt_tate_bw6_bls(u, 24)

def cost_pairing_cp8_cp12(r, a0, a1, a2, a3, b0, b1, e0, e1, e2, e3, k=8):
    """cost of Cocks-Pinch curves of embedding degree 8, 12

    INTPUT:
    - `r`: subgroup order, prime
    - `ai`: scalars for the optimal ate pairing computation
    - `bi`: scalars for the optimal Tate pairing computation
    - `ei`: scalars for the hard part of final exp, decomposed in basis 1, p, p^2, p^3
    """
    assert k == 8 or k == 12
    m = 1
    s = 1
    m2 = 3*m
    s2 = 2*m
    inv = 25*m
    i2 = 2*m + 2*s + inv
    i4 = 2*m2 + 2*s2 + i2
    i6 = 9*m2 + 3*s2 + i2 # because one starts with a quadratic extension for Fp12 anyway
    if k == 8:
        # curves with b=0
        m4 = 3*m2   # = 9m
        s4 = 2*m2   # = 6m
        m8 = 3*m4   # = 27m
        s8 = 2*m4   # = 18m
        s8_cyclo = 2*s4 # = 4*m2
        f8 = 7*m
        f2 = 4*m
        i8 = 2*m4 + 2*s4 + i4 # 44*m + inv
        assert i8 == 44*m + inv
        me=m2; se=s2; ie = i2; sk=s8; mk = m8 # e = k/d = 8/4 = 2
        fk = f8
        ik = i8
        sk_cyclo = s8_cyclo
        # Costello Lange Naehrig, with d=4, b=0
        double_line_ate_cln = 2*me+8*se+(k//2)*m
        add_line_ate_cln    = 9*me+5*se+(k//2)*m
        add_line_ate_cln_with_z = 15*me+7*se+(k//2)*m
        add_line_cln_affine = 2*i2 + s2 + 4*m2
        add_line_affine = add_line_cln_affine
        sparse_sparse = 6*me
        sparse_mult = 8*me
        update1        = 8*me+sk
        update2        = 8*me    # in other words, a sparse multiplication costs 8 m4 (8, not 7)
    else:
        m6 = 6*m2 # multiplication in cubic extension
        #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
        s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
        m12 = 3*m6
        s12 = 2*m6
        f2 = 8*m
        f12 = 10*m
        s12_cyclo = 18*m # 3*s4
        i12 = 97*m+inv # 94*m + i in Guillevic Masson Thome DCC'20 is not compatible with the tower
        #i12 = min(9*m4 + 3*s4 + i4, 2*m6 + 2*s6 + i6)
        me=m2 ; se=s2; ie = i2; sk=s12; mk = m12 # e = k/d
        fk = f12
        ik = i12
        sk_cyclo = s12_cyclo
        # aklgl
        double_line_ate = 3*me+6*se+(k//3)*m
        add_line_ate    = 11*me+2*se+(k//3)*m
        add_line_ate_aklgl_with_z = 16*me + 2*se + 2*(k//6)*m
        # Costello Lange Naehrig
        double_line_ate_cln = 2*me+7*se+(k//3)*m
        add_line_ate_cln    = 10*me+2*se+(k//3)*m
        add_line_ate_cln_with_z = 15*me+7*se+(k//3)*m # TODO, this is the formula for CP8
        # for add_line_affine for k = 0 mod 6, take INDOCRYPT 2012 paper by Zhang and Lin
        add_line_affine = ie + 3*me + 2*se + (k//3)*m
        sparse_mult = 13*me
        sparse_sparse = 6*me
        sparse_mult_tate = 10*me + 3*(k//6)*m
        sparse_sparse_tate = 5*me + m
        update1        = sparse_mult+sk

    frobenius_pQ = 2*m2 # upper bound, but should be cheaper
    frobenius_p2Q = 2*m # upper bound, but should be cheaper
    frobenius_p3Q = 2*m2 # upper bound, but should be cheaper
    # Frobenius(Q) : sparse Frobenius in Fpk (actually cheaper than 2*fk)

    list_a = [ZZ(a0), ZZ(a1), ZZ(a2), ZZ(a3)]
    list_b = [ZZ(b0), ZZ(b1)]
    list_e = [ZZ(e0), ZZ(e1), ZZ(e2), ZZ(e3)]

    # 1. optimal ate pairing
    # for the multi-scalar mult: counts the number of non-zero bits for the a_i at each index i
    abs_a = [abs(ai) for ai in list_a]
    abs_a.sort(reverse = True)
    digits_a = [ai.digits(2) for ai in abs_a]
    len_ai = [len(di) for di in digits_a]
    ll0, ll1, ll2, ll3 = len_ai
    max_len = max(len_ai)
    digits_a = [digits_a[i] + [0]*(max_len - len(digits_a[i])) for i in range(4)] # make all the same length
    hw_a_idx_i = [sum([digits_a[j][i] for j in range(4)]) for i in range(max_len)]
    hw_a = [sum(digits_ai) for digits_ai in digits_a]

    bits_2naf_a = [bits_2naf(abs(ai)) for ai in list_a]
    len_ai_2naf = [len(bits_2naf_a[i]) for i in range(4)]
    max_len_2naf = max(len_ai_2naf)
    bits_2naf_a = [bits_2naf_a[i] + [0]*(max_len - len(bits_2naf_a[i])) for i in range(4)] # make all the same length

    hw2naf_a = [sum([1 for bi in bits_2naf_ai if bi != 0]) for bits_2naf_ai in bits_2naf_a]
    hw_a_idx_i_2naf = [sum([abs(bits_2naf_a[j][i]) for j in range(4)]) for i in range(max_len)]

    # sparse-sparse-mult CP8 :
    # tangents are of the form l0 + l1*w + l3*w^3 where w^4 = v
    # lines are of the form    l0 + l1*w + l3*w^3 where w^4 = v
    # alternatively, both are of the form  l0 + l2*w^2 + l3*w^3 where w^4 = v
    # in both cases, the sparse-sparse costs 6 m2.
    # then the total cost is: 6*m2 + m8 = 18+27 m = 45m (15*m2)
    # to be compared to 2*(8m2) = 16m2 = 48m
    # now estimate the cost of the Miller loop
    # with multi-scalar mult
    # new version with shared point doubling ans shared squaring in Fp8
    pi_pQ = frobenius_pQ
    pi_p2Q = frobenius_p2Q
    pi_p3Q = frobenius_p3Q
    precomputations = pi_pQ + pi_p2Q + pi_p3Q + 11 * add_line_affine + 4 * sparse_sparse + update1
    cost_miller_f = 1
    cost_miller_f_2naf = 1
    for i in range(max_len-2, -1, -1):
        cost_miller_f += sk + double_line_ate_cln
        cost_miller_f_2naf += sk + double_line_ate_cln
        j = hw_a_idx_i[i] # number of non-zero bits at index i
        j2 = hw_a_idx_i_2naf[i] # number of non-zero signed bits at index i
        if j == 0:
            cost_miller_f += sparse_mult
        else:
            cost_miller_f += add_line_ate_cln + sparse_sparse + mk
            # if j == 1 there is no precomputed line to accumulate
            # but if j > 1 there are.
            if j == 2:
                cost_miller_f += sparse_mult # the precomputed line has a sparse form
            elif j > 2:
                cost_miller_f += mk # the precomputed line does not have a sparse form
        if j2 == 0:
            cost_miller_f_2naf += sparse_mult
        else:
            cost_miller_f_2naf += add_line_ate_cln + sparse_sparse + mk
            if j2 == 2:
                cost_miller_f_2naf += sparse_mult
            elif j2 > 2:
                cost_miller_f_2naf += mk
    print("cost Miller ate with 4-scalar-mult     : {} m, precomp. {} m, total {} m".format(cost_miller_f, precomputations, cost_miller_f+precomputations))
    print("cost Miller ate with 4-scalar-mult 2naf: {} m, precomp. {} m, total {} m (not sure for precomp)".format(cost_miller_f_2naf, precomputations, cost_miller_f_2naf+precomputations))

    # now if each of the 4 Miller functions is done in parallel:
    cost_Miller_f = [(len_ai[i]-1)*(sk + double_line_ate_cln) + (hw_a[i] - 1) * (add_line_ate_cln + sparse_sparse + mk) + (len_ai[i]-1 - (hw_a[i] - 1))*sparse_mult for i in range(4)]
    print("cost Miller each function:      {} m".format(cost_Miller_f))
    cost_Miller_f_2naf = [(len_ai[i]-1)*(sk + double_line_ate_cln) + (hw2naf_a[i] - 1) * (add_line_ate_cln + sparse_sparse + mk) + (len_ai[i]-1 - (hw2naf_a[i] - 1))*sparse_mult for i in range(4)]
    print("cost Miller each function 2naf: {} m".format(cost_Miller_f_2naf))
    additional_terms = max([frobenius_pQ, frobenius_p2Q, frobenius_p3Q]) + 2*add_line_ate_cln_with_z + sparse_sparse + mk
    print("additional terms (lines, Frobenius in parallel (pQ, p2Q, p3Q)): {} m".format(additional_terms))

    # optimal Tate pairing with scalars b0, b1
    abs_b = [abs(bi) for bi in list_b]
    abs_b.sort(reverse = True)
    digits_b = [bi.digits(2) for bi in abs_b]
    len_bi = [len(di) for di in digits_b]
    ll0, ll1 = len_bi
    max_len = max(len_bi)
    digits_b = [digits_b[i] + [0]*(max_len - len(digits_b[i])) for i in range(2)] # make all the same length
    hw_b_idx_i = [sum([digits_b[j][i] for j in range(2)]) for i in range(max_len)]
    hw_b = [sum(digits_bi) for digits_bi in digits_b]

    bits_2naf_b = [bits_2naf(abs(bi)) for bi in list_b]
    len_bi_2naf = [len(bits_2naf_b[i]) for i in range(2)]
    max_len_2naf = max(len_bi_2naf)
    bits_2naf_b = [bits_2naf_b[i] + [0]*(max_len - len(bits_2naf_b[i])) for i in range(2)] # make all the same length

    hw2naf_b = [sum([1 for bi in bits_2naf_bi if bi != 0]) for bits_2naf_bi in bits_2naf_b]
    hw_b_idx_i_2naf = [sum([abs(bits_2naf_b[j][i]) for j in range(2)]) for i in range(max_len)]

    print("optimal Tate pairing with one endomorphism")
    if k == 8:
        add_line_tate_cln_affine = inv + i2 + s + (k//4+2)*m + m2
        add_line_tate_cln = 9*m + 5*s + (k//2)*m
        double_line_tate_cln = 2*m + 8*s + (k//2)*m
        add_line_tate_affine = add_line_tate_cln_affine
        add_line_tate = add_line_tate_cln
        double_line_tate = double_line_tate_cln
        sparse_mult_tate = sparse_mult
        sparse_sparse_tate = sparse_sparse
    else:
        add_line_tate_zl_affine = inv + 3*m + 2*s + 2*(k//6)*m # Zhang Lin INDOCRYPT 2012
        add_line_tate_aklgl = 11*m + 2*s + 2*(k//6)*m
        double_line_tate_aklgl = 3*m + 6*s + 2*(k//6)*m
        add_line_tate_affine = add_line_tate_zl_affine
        add_line_tate = add_line_tate_aklgl
        double_line_tate = double_line_tate_aklgl
        sparse_mult_tate = 10*me + 3*(k//6)*m
        sparse_sparse_tate = 5*me + m
        print("k={} sk = {}m mk = {}m double_line_tate = {}m add_line_tate = {}m sparse_sparse_tate = {}m sparse_tate = {}m".format(k, sk, mk, double_line_tate, add_line_tate, sparse_sparse_tate, sparse_mult_tate))

    additional_terms = m + add_line_tate_affine # compute psi(P), P+psi(P) and the line through P, psi(P) evaluated at Q.

    print("sk = {}m mk = {}m double_line_tate = {}m add_line_tate = {}m sparse_sparse_tate = {}m sparse_tate = {}m".format(sk, mk, double_line_tate, add_line_tate, sparse_sparse_tate, sparse_mult_tate))
    cost_Miller_f = 0
    for i in range(max_len-2, -1, -1):
        cost_Miller_f += sk + double_line_tate
        if hw_b_idx_i[i] > 0:
            cost_Miller_f += add_line_tate + sparse_sparse_tate + mk
            if hw_b_idx_i[i] == 2:
                cost_Miller_f += sparse_mult_tate
        else:
            cost_Miller_f += sparse_mult_tate
    print("cost Miller with 2-multi-scalar: {} m".format(cost_Miller_f))
  
    cost_Miller_f_2naf = 0
    additional_terms = m + 2*add_line_tate_affine # compute psi(P), P+psi(P), P-psi(P) and the line through P, psi(P) evaluated at Q, the line through P, P-psi(P) evaluated at Q.
    for i in range(max_len_2naf-2, -1, -1):
        cost_Miller_f_2naf += sk + double_line_tate
        if hw_b_idx_i_2naf[i] > 0:
            cost_Miller_f_2naf += add_line_tate + sparse_sparse_tate + mk
            if hw_b_idx_i_2naf[i] == 2: # multiply by a precomputed sparse line
                cost_Miller_f_2naf += sparse_mult_tate
        else:
            cost_Miller_f_2naf += sparse_mult_tate
    print("cost Miller with 2-multi-scalar, 2naf: {} m".format(cost_Miller_f_2naf))

    # now if each of the 2 Miller functions is done in parallel:
    # if the i-th bit is 0: sk + double_line + sparse_mult
    # if the i-th bit is 1: sk + double_line + add_line + sparse_sparse + mk
    cost_Miller_f = [(len_bi[i]-1)*(sk + double_line_tate) + (hw_b[i] - 1) * (add_line_tate + sparse_sparse_tate + mk) + (len_bi[i]-1 - (hw_b[i] - 1))*(sparse_mult_tate) for i in range(2)]
    print("cost Miller each function:      {} m".format(cost_Miller_f))
    cost_Miller_f_2naf = [(len_bi[i]-1)*(sk + double_line_tate) + (hw2naf_b[i] - 1) * (add_line_tate + sparse_sparse_tate + mk) + (len_bi[i]-1 - (hw2naf_b[i] - 1))*(sparse_mult_tate) for i in range(2)]
    print("cost Miller each function 2naf: {} m".format(cost_Miller_f_2naf))
    additional_terms = m + mk # compute psi(P) = (-x, y*omega) and multiply the two Miller functions at the end
    print("additional terms (psi(P) = (-x, omega*y) and 1 m{}: {} m".format(k, additional_terms))

    print("cost final exp:")
    # in parallel:
    abs_e = [abs(ei) for ei in list_e]
    abs_e.sort(reverse = True)
    digits_e = [ei.digits(2) for ei in abs_e]
    len_ei = [len(di) for di in digits_e]
    ll0, ll1, ll2, ll3 = len_ei
    max_len = max(len_ei)
    digits_e = [digits_e[i] + [0]*(max_len - len(digits_e[i])) for i in range(4)] # make all the same length
    hw_e_idx_i = [sum([digits_e[j][i] for j in range(4)]) for i in range(max_len)]

    hw_e = [sum(digits_ei) for digits_ei in digits_e]
    bits_2naf_e = [bits_2naf(abs(ei)) for ei in list_e]
    len_ei_2naf = [len(bits_2naf_e[i]) for i in range(4)]
    max_len_2naf = max(len_ei_2naf)
    bits_2naf_e = [bits_2naf_e[i] + [0]*(max_len_2naf - len(bits_2naf_e[i])) for i in range(4)]
    hw2naf_e = [sum([1 for bi in bits_2naf_ei if bi != 0]) for bits_2naf_ei in bits_2naf_e]
    cost_final_exp_hard = [(len_ei[i]-1)*sk_cyclo + (hw_e[i]-1) * mk for i in range(4)]
    cost_final_exp_hard[1] += fk
    cost_final_exp_hard[2] += fk
    cost_final_exp_hard[3] += fk

    cost_final_exp_hard_2naf = [(len_ei[i]-1)*sk_cyclo + (hw2naf_e[i]-1) * mk for i in range(4)]
    cost_final_exp_hard_2naf[1] += fk
    cost_final_exp_hard_2naf[2] += fk
    cost_final_exp_hard_2naf[3] += fk

    #print("ei = [{:#x}, {:#x}, {:#x}, {:#x}]".format(e0,e1,e2,e3))
    print("ei of {} bits, HW {}, 2-naf {}".format(len_ei, hw_e, hw2naf_e))
    print("cost final exp hard with 4 exp in parallel:       {} m".format(cost_final_exp_hard))
    print("cost final exp hard with 4 exp in parallel, 2naf: {} m".format(cost_final_exp_hard_2naf))
    if k == 8:
        easy_part = ik + mk # f^(p^4) - 1
    else:
        easy_part = ik + mk + f2 + mk # (p^6-1)*(p^2+1)
    print("cost final exp easy:   {}".format(easy_part))
    print("total final exp in parallel:       {} m".format(easy_part + max(cost_final_exp_hard) + 2*mk)) # one can do 2 mults in parallel, then 1 mult, then the 3 mults are done as 2 mults
    print("total final exp in parallel, 2naf: {} m".format(easy_part + max(cost_final_exp_hard_2naf) + 2*mk)) # one can do 2 mults in parallel, then 1 mult, then the 3 mults are done as 2 mults

    cost_final_exp_hard_multi_exp = 0
    for i in range(len_ei[0]-2, len_ei[1]-2, -1):
        cost_final_exp_hard_multi_exp += sk_cyclo
        if digits_e[0][i] != 0:
            cost_final_exp_hard_multi_exp += mk
    for i in range(len_ei[1]-2, len_ei[2]-2, -1):
        cost_final_exp_hard_multi_exp += sk_cyclo
        if digits_e[0][i] != 0 or digits_e[1][i] != 0:
            cost_final_exp_hard_multi_exp += mk
    for i in range(len_ei[2]-2, len_ei[3]-2, -1):
        cost_final_exp_hard_multi_exp += sk_cyclo
        if digits_e[0][i] != 0 or digits_e[1][i] != 0 or digits_e[2][i] != 0:
            cost_final_exp_hard_multi_exp += mk
    for i in range(len_ei[3]-2, -1, -1):
        cost_final_exp_hard_multi_exp += sk_cyclo
        if digits_e[0][i] != 0 or digits_e[1][i] != 0 or digits_e[2][i] != 0 or digits_e[3][i] != 0:
            cost_final_exp_hard_multi_exp += mk
    print("cost final exp hard with multi-exp technique      : {} m".format(cost_final_exp_hard_multi_exp))
    print("estimated cost of precomputations: 3 f{0} + (15-4) m{0} = 3 f{0} + 11 m{0} = {1} m, or maybe replace some multiplications by Frobenius powers".format(k, 3*fk + 11*mk))
    cost_final_exp_hard_multi_exp_2naf = 0
    for i in range(len_ei_2naf[0]-2, len_ei_2naf[1]-2, -1):
        cost_final_exp_hard_multi_exp_2naf += sk_cyclo
        if bits_2naf_e[0][i] != 0:
            cost_final_exp_hard_multi_exp_2naf += mk
    for i in range(len_ei_2naf[1]-2, len_ei_2naf[2]-2, -1):
        cost_final_exp_hard_multi_exp_2naf += sk_cyclo
        if bits_2naf_e[0][i] != 0 or bits_2naf_e[1][i] != 0:
            cost_final_exp_hard_multi_exp_2naf += mk
    for i in range(len_ei_2naf[2]-2, len_ei_2naf[3]-2, -1):
        cost_final_exp_hard_multi_exp_2naf += sk_cyclo
        if bits_2naf_e[0][i] != 0 or bits_2naf_e[1][i] != 0 or bits_2naf_e[2][i] != 0:
            cost_final_exp_hard_multi_exp_2naf += mk
    for i in range(len_ei_2naf[3]-2, -1, -1):
        cost_final_exp_hard_multi_exp_2naf += sk_cyclo
        if bits_2naf_e[0][i] != 0 or bits_2naf_e[1][i] != 0 or bits_2naf_e[2][i] != 0 or bits_2naf_e[3][i] != 0:
            cost_final_exp_hard_multi_exp_2naf += mk
    print("cost final exp hard with multi-exp technique, 2naf: {} m".format(cost_final_exp_hard_multi_exp_2naf))

def cost_tate_pairing_cp12_disc3(r, u, k, k_bls, t_1_e_mod_r_mod_u):
    """cost of Cocks-Pinch curves of embedding degree 12 and discriminant 3

    INTPUT:
    - `r`: subgroup order, prime, r(x) = (x-1)^2/3*(x^4-x^2+1) + x (=p_bls)
    - `u`: BLS seed
    - `k_cp`: Cocks-Pinch embedding degree 12
    - `k_bls`: BLS embedding degree (12 or 24)

    if k == 12: (t-1)^(12/6) = (t-1)^2 is a primitive 6-th root of unity mod r, and has a polynomial form in u
    """
    assert k == 12
    assert k_bls == 12 or k_bls == 24
    m = 1
    s = 1
    inv = 25*m
    m2 = 3*m
    s2 = 2*m
    i2 = 4*m + inv
    m6 = 6*m2 # multiplication in cubic extension
    #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    m12 = 3*m6
    s12 = 2*m6
    f12_p2 = 8*m # p^2 power in Fp12
    f12 = 10*m
    s12_cyclo = 18*m # 3*s4
    i12 = 97*m+inv # 94*m + i in Guillevic Masson Thome DCC'20 is not compatible with the tower
    me=m2 ; se=s2; ie=i2; sk=s12; mk = m12 # e = k/d
    fk = f12
    ik = i12
    sk_cyclo = s12_cyclo
    double_line_ate = 3*me+6*se+(k//3)*m
    add_line_ate    = 11*me+2*se+(k//3)*m
    # Costello Lange Naehrig
    double_line_ate_cln = 2*me+7*se+(k//3)*m
    add_line_ate_cln    = 10*me+2*se+(k//3)*m
    add_line_ate_cln_with_z = 15*me+7*se+(k//3)*m # TODO, this is the formula for CP8
    # for add_line_affine for k = 0 mod 6, take INDOCRYPT 2012 paper by Zhang and Lin
    add_line_affine = ie + 3*me + 2*se + (k//3)*m
    double_line_affine = ie + 3*me + 2*se + (k//3)*m
    # AKLGL
    double_line_ate_aklgl = 3*me + 6*se + 2*(k//6)*m
    add_line_ate_aklgl = 11*me + 2*se + 2*(k//6)*m
    add_line_ate_aklgl_with_z = 16*me + 2*se + 2*(k//6)*m

    double_line_tate_aklgl = 3*m + 6*s + 2*(k//6)*m
    add_line_tate_aklgl = 11*m + 2*s + 2*(k//6)*m
    add_line_tate_aklgl_with_z = 16*m + 2*s + 2*(k//6)*m

    add_line_tate = add_line_tate_aklgl
    double_line_tate = double_line_tate_aklgl

    add_line_tate_affine = inv + 3*m + 2*s + 2*(k//6)*m
    sparse_line_mult = 10*me + 3*(k//6)*m
    sparse_sparse = 5*me + m
    print("k={} sk = {}m mk = {}m double_line_tate = {}m add_line_tate = {}m sparse_sparse_tate = {}m sparse_tate = {}m".format(k, sk, mk, double_line_tate, add_line_tate, sparse_sparse, sparse_line_mult))

    # t_1_e_mod_r_mod_u
    if k_bls == 12 and k == 12 and t_1_e_mod_r_mod_u == -1:
        # case corresponding to t = 0 mod r mod u for BLS12-BW6 curves
        list_a1 = (-(u+1), u**3-u**2+1, "(-(u+1), u^3-u^2+1)")
        list_a2 = (u**3-u**2-u, u+1, "(u^3-u^2-u, u+1)")
    elif k_bls == 12 and k == 12 and t_1_e_mod_r_mod_u == 2:
        # case corresponding to t = 3 mod r mod u for BLS12-BW6 curves
        list_a1 = (u+1, u**3-u**2-u, "(u+1, u^3-u^2-u)")
        list_a2 = (u**3-u**2+1, -(u+1), "(u^3-u^2+1, -(u+1))")
    # for k=9, (tr-1)^3 is a primitive cube root of unity and we should have
    # almost the same scalars, (tr-1)^3 plays the role of -(t-1) for BW6 curves (?)
    elif k_bls == 12 and k == 9 and t_1_e_mod_r_mod_u == 1:
        # case corresponding to t = 0 mod r mod u for BLS12-BW6 curves
        list_a1 = (u+1, u**3-u**2+1, "(u+1, u^3-u^2+1)")
        list_a2 = (u**3-u**2-u, -(u+1), "(u^3-u^2-u, -(u+1))")
    elif k_bls == 12 and k == 9 and t_1_e_mod_r_mod_u == -2:
        # case corresponding to t = 3 mod r mod u for BLS12-BW6 curves
        list_a1 = (-(u+1), u**3-u**2-u, "(-(u+1), u^3-u^2-u)")
        list_a2 = (u**3-u**2+1, u+1, "(u^3-u^2+1, u+1)")
    elif k_bls == 24 and k == 12 and abs(t_1_e_mod_r_mod_u) == 1: # +/- 1, not clear which one?
        # case corresponding to t = 0 mod r mod u for BLS24-BW6 curves
        list_a1 = (-(u+1), u**5-u**4+1, "(-(u+1), u^5-u^4+1)")
        list_a2 = (u**5-u**4-u, u+1, "(u^5-u^4-u, u+1)")
    elif k_bls == 24 and k == 12 and abs(t_1_e_mod_r_mod_u) == 2: # +/- 1, not clear which one?
        # case corresponding to t = 3 mod r mod u for BLS24-BW6 curves
        list_a1 = (u+1, u**5-u**4-u, "(u+1, u^5-u^4-u)")
        list_a2 = (u**5-u**4+1, -(u+1), "(u^5-u^4+1, -(u+1))")
    # for k=9, (tr-1)^3 is a primitive cube root of unity and we should have
    # almost the same scalars, (tr-1)^3 plays the role of -(t-1) for BW6 curves (?)
    elif k_bls == 24 and k == 9 and t_1_e_mod_r_mod_u == 1:
        # case corresponding to t = 0 mod r mod u for BLS24-BW6 curves
        list_a1 = (u+1, u**5-u**4+1, "(u+1, u^5-u^4+1)")
        list_a2 = (u**5-u**4-u, -(u+1), "(u^5-u^4-u, -(u+1))")
    elif k_bls == 24 and k == 9 and t_1_e_mod_r_mod_u == -2:
        # case corresponding to t = 3 mod r mod u for BLS24-BW6 curves
        list_a1 = (-(u+1), u**5-u**4-u, "(-(u+1), u^5-u^4-u)")
        list_a2 = (u**5-u**4+1, u+1, "(u^5-u^4+1, u+1)")
    else:
        print("k = {} k_bls = {} t_1_e_mod_r_mod_u = {} no match found".format(k, k_bls, t_1_e_mod_r_mod_u))
    for (b0,b1,b_str) in [list_a1, list_a2]:
        list_b = [ZZ(b0), ZZ(b1)]
        print("with scalars {}".format(b_str))
        # optimal Tate pairing with scalars b0, b1
        abs_b = [abs(bi) for bi in list_b]
        abs_b.sort(reverse = True)
        digits_b = [bi.digits(2) for bi in abs_b]
        len_bi = [len(di) for di in digits_b]
        ll0, ll1 = len_bi
        max_len = max(len_bi)
        digits_b = [digits_b[i] + [0]*(max_len - len(digits_b[i])) for i in range(2)] # make all the same length
        hw_b_idx_i = [sum([digits_b[j][i] for j in range(2)]) for i in range(max_len)]
        hw_b = [sum(digits_bi) for digits_bi in digits_b]

        bits_2naf_b = [bits_2naf(abs(bi)) for bi in list_b]
        len_bi_2naf = [len(bits_2naf_b[i]) for i in range(2)]
        max_len_2naf = max(len_bi_2naf)
        bits_2naf_b = [bits_2naf_b[i] + [0]*(max_len - len(bits_2naf_b[i])) for i in range(2)] # make all the same length

        hw2naf_b = [sum([1 for bi in bits_2naf_bi if bi != 0]) for bits_2naf_bi in bits_2naf_b]
        hw_b_idx_i_2naf = [sum([abs(bits_2naf_b[j][i]) for j in range(2)]) for i in range(max_len)]

        print("optimal Tate pairing with one endomorphism")
        print("k={} sk = {}m mk = {}m double_line_tate = {}m add_line_tate = {}m sparse_sparse_tate = {}m sparse_tate = {}m".format(k, sk, mk, double_line_tate, add_line_tate, sparse_sparse, sparse_line_mult))
        cost_Miller_f = 0
        additional_terms = m + add_line_tate_affine # compute psi(P) (m), P+psi(P) and the line through P, psi(P) evaluated at Q
        for i in range(max_len-2, -1, -1):
            cost_Miller_f += sk + double_line_tate
            if hw_b_idx_i[i] > 0:
                cost_Miller_f += add_line_tate + sparse_sparse + mk
                if hw_b_idx_i[i] == 2:
                    cost_Miller_f += sparse_line_mult
            else:
                cost_Miller_f += sparse_line_mult
        print("cost Miller with 2-multi-scalar: {} m".format(cost_Miller_f))
        print("cost Miller with 2-multi-scalar: {} m + {} m = {}".format(cost_Miller_f, additional_terms, cost_Miller_f+additional_terms))

        cost_Miller_f_2naf = 0
        additional_terms = m + 2*add_line_tate_affine # compute psi(P) (m), P+psi(P), P-psi(P) and the line through P, psi(P) evaluated at Q, the line through P, P-psi(P) evaluated at Q.
        for i in range(max_len_2naf-2, -1, -1):
            cost_Miller_f_2naf += sk + double_line_tate
            if hw_b_idx_i_2naf[i] > 0:
                cost_Miller_f_2naf += add_line_tate + sparse_sparse + mk
                if hw_b_idx_i_2naf[i] == 2:
                    cost_Miller_f_2naf += sparse_line_mult
            else:
                cost_Miller_f_2naf += sparse_line_mult
        print("cost Miller with 2-multi-scalar, 2naf: {} m + {} m = {}".format(cost_Miller_f_2naf, additional_terms, cost_Miller_f_2naf+additional_terms))

def cost_miller_function_ate_bw6_twist6(u, m6, s6, m, s, with_m0=False, naf=False):
    k = 6
    add_line_ate_aklgl = 11*m + 2*s + (k//3)*m
    double_line_ate_aklgl = 3*m + 6*s + (k//3)*m

    sparse_sparse_mult = 6*m # aklgl
    sparse_mult = 13*m # aklgl

    u = ZZ(abs(u))
    if naf:
        bits_u = bits_2naf(u)
        hw_u = sum([1 for bi in bits_u if bi != 0])
    else:
        bits_u = u.digits(2)
        hw_u = sum(bits_u)
    if not with_m0:
        cost_miller = (len(bits_u)-1)*(double_line_ate_aklgl + s6) - s6 + (hw_u-1)*(add_line_ate_aklgl+sparse_sparse_mult+m6) + (len(bits_u)-2 - (hw_u-1))*sparse_mult
    else:
        cost_miller = (len(bits_u)-1)*(double_line_ate_aklgl + s6) + (hw_u-1)*(add_line_ate_aklgl+sparse_sparse_mult+2*m6) + (len(bits_u)-1 - (hw_u-1))*sparse_mult

    return cost_miller

def cost_miller_loop_bw6_bn(u):
    """
    cost pairing

    Miller loop optimal ate:
    t = 0: f_{2u,Q}(P)*frobenius(f_{2u(3u+1)+1,Q}(P))
    t = 3: frobenius(f_{2u,Q}(P))*f_{2u(3u+1)+1,Q}(P)
    same cost in both cases t=0 and t=3.

    with algo 2: same formulas but different miller loop

    final exponentiation
    t = 0:
    t = 3:

    """
    m = 1
    s = 1
    inv = 25*m
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    f6 = 4*m
    k = 6

    add_line_ate_aklgl = 11*m + 2*s + (k//3)*m
    sparse_line_mult = 13*m # aklgl

    cost1 = cost_miller_function_ate_bw6_twist6(2*u, m6, s6, m, s, with_m0=False, naf=False)
    cost2 = cost_miller_function_ate_bw6_twist6(3*u+1, m6, s6, m, s, with_m0=True, naf=False)

    additional_terms = inv + 2*m + add_line_ate_aklgl + sparse_line_mult + f6 + m6
    cost = cost1 + cost2 + additional_terms

    cost1_2naf = cost_miller_function_ate_bw6_twist6(2*u, m6, s6, m, s, with_m0=False, naf=True)
    cost2_2naf = cost_miller_function_ate_bw6_twist6(3*u+1, m6, s6, m, s, with_m0=True, naf=True)

    cost_2naf = cost1_2naf + cost2_2naf + additional_terms

    print("cost Miller loop      (2u, (3u+1)*2u+1): miller_1 + miller_2 + i + 2m + add_line + sparse_mult + f6 + m6 = {} + {} + {} = {}m".format(cost1, cost2, additional_terms, cost))
    print("cost Miller loop 2naf (2u, (3u+1)*2u+1): miller_1 + miller_2 + i + 2m + add_line + sparse_mult + f6 + m6 = {} + {} + {} = {}m".format(cost1_2naf, cost2_2naf, additional_terms, cost_2naf))
    return cost, cost_2naf

def cost_miller_optimal_twisted_ate_bw6_bn(u):
    b0 = 2*u
    b1 = (3*u+1)*b0 + 1
    k = 6
    m = 1
    s = 1
    inv = 25*m
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    double_line_tate_aklgl = 3*m + 6*s + (k//3)*m
    add_line_tate_aklgl = 11*m + 2*s + (k//3)*m
    add_line_affine = inv + 3*m + 2*s + (k//3)*m
    sparse_sparse_mult = 6*m
    sparse_mult = 13*m
    cost, cost_2naf = cost_miller_optimal_twisted_ate_bw6(b0, b1, s6, m6, m, double_line_tate_aklgl, add_line_tate_aklgl, add_line_affine, sparse_sparse_mult, sparse_mult)

    print("optimal Tate pairing with one endomorphism")
    print("with scalars (2*u, (3*u+1)*(2*u) + 1")
    print("cost Miller with 2-multi-scalar:       {} m".format(cost))
    print("cost Miller with 2-multi-scalar, 2naf: {} m".format(cost_2naf))
    return cost, cost_2naf

def cost_final_exp_easy_k6():
    m = 1
    s = 1
    inv = 25*m
    i6 = 34*m + inv
    m2 = 3*m
    m6 = 6*m2 # multiplication in cubic extension
    f6 = 4*m
    cost = i6 + 2*m6 + f6
    print("final exp easy k=6: {}m".format(cost))
    return cost

def cost_exp_u_bw6(u, m6, s6_cyclo, naf=False):
    u = ZZ(abs(u))
    if naf:
        bits_u = bits_2naf(u)
        hw_u = sum([1 for bi in bits_u if bi != 0])
    else:
        bits_u = u.digits(2)
        hw_u = sum(bits_u)
    cost = max(0, len(bits_u)-1)*s6_cyclo + max(0, hw_u-1) * m6
    return cost

def cost_final_exp_hard_bw6_bn(u, ht, hy, tr_mod=0):
    m = 1
    s = 1
    inv = 25*m

    i2 = 2*m +2*s + inv
    i6 = 34*m + inv
    m2 = 3*m
    s2 = 2*m
    m6 = 6*m2 # multiplication in cubic extension
    #s6 = m2+4*s2 # =3*m+4*2*m = 11m squaring in cubic extension
    s6 = 2*m2+3*s2 # =2*3*m+3*2*m = 12m squaring in cubic extension
    s6_cyclo = 3*s2 # sk_cyclo = 3*s_{k/3} for 6 | k
    f6 = 4*m

    d2 = (ht**2+3*hy**2)//4
    exp_d2 = cost_exp_u_bw6(d2, m6, s6_cyclo, naf=False)
    exp_d2_2naf = cost_exp_u_bw6(d2, m6, s6_cyclo, naf=True)
    exp_u = cost_exp_u_bw6(u, m6, s6_cyclo, naf=False)
    exp_u_2naf = cost_exp_u_bw6(u, m6, s6_cyclo, naf=True)

    if tr_mod == 0:
        # (2*u) * Phi_k(p)/r
        d1 = (ht-hy)//2
        print("ht={} hy={} (ht^2+3*hy^2)/4={} (ht-hy)/2={}".format(ht, hy, d2, d1))
        # 6 exp(u) + exp(d2) + exp(d1) + 6 S + 15 M + 2 Frobenius + 3 conjugate
        other_terms = 6*s6_cyclo + 15*m6 + 2*f6
        #
        # (6*u^2+4*u+1) * Phi_k(p)/r
        # 6 exp(u) + exp(d2) + exp(d1) + 7 S + 15 M + 2 Frobenius + 4 conjugate
        # cost2 = 6*exp_u + exp_d1 + exp_d2 + 7*s6_cyclo + 15*m + 2*f6
    else:
        # 2*u * Phi_k(p)/r
        # 6 exp(u) + exp(d2) + exp(d1) + 7 S + 16 M + 2 Frobenius + 2 conjugate
        d1 = (ht+hy)//2
        print("ht={} hy={} (ht^2+3*hy^2)/4={} (ht+hy)/2={}".format(ht, hy, d2, d1))
        # (6*u^2+2*u+1) * Phi_k(p)/r
        # 6 exp(u) + exp(d2) + exp(d1) + 6 S + 16 M + 2 Frobenius
        other_terms = 6*s6_cyclo + 16*m6 + 2*f6

    if d1 == 0:
        other_terms -= m6
    if d2 == 0:
        other_terms -= m6
    exp_d1 = cost_exp_u_bw6(d1, m6, s6_cyclo, naf=False)
    exp_d1_2naf = cost_exp_u_bw6(d1, m6, s6_cyclo, naf=True)
    cost1 = 6*exp_u + exp_d1 + exp_d2 + other_terms
    cost1_2naf = 6*exp_u_2naf + exp_d1_2naf + exp_d2_2naf + other_terms
    print("final exp hard:       {} (6 exp u) + {} (exp (ht^2+3*hy^2)/4) + {} (exp (ht-hy)/2) + {} = {} m".format(6*exp_u, exp_d2, exp_d1, other_terms, cost1))
    print("final exp hard, 2naf: {} (6 exp u) + {} (exp (ht^2+3*hy^2)/4) + {} (exp (ht-hy)/2) + {} = {} m".format(6*exp_u_2naf, exp_d2_2naf, exp_d1_2naf, other_terms, cost1_2naf))
    return cost1, cost1_2naf

def cost_pairing_bw6_bn(u, ht, hy, tr_mod=0):
    cost, cost_2naf = cost_miller_loop_bw6_bn(u)
    cost_tate, cost_tate_2naf = cost_miller_optimal_twisted_ate_bw6_bn(u)
    cost_easy = cost_final_exp_easy_k6()
    cost_hard, cost_hard_2naf = cost_final_exp_hard_bw6_bn(u, ht, hy, tr_mod=tr_mod)
    cost_miller = min(cost, cost_2naf)
    cost_exp_hard = min(cost_hard, cost_hard_2naf)
    print("total: {} + {} + {} = {}".format(cost_miller, cost_easy, cost_exp_hard, cost_miller+cost_easy+cost_exp_hard))
