from pairing import *

def miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u(Q, P, u):
    """
    Return f_{u*(u^4-u^3-1),Q}(P)*frobenius(f_{u+1,Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 0
    trace = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x + r(x) * ht
    v = (u^4-u^3-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate((Q[0], -Q[1]), P, 0, -u)
    else:
        m_u, uQ = miller_function_ate(Q, P, 0, u)
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, u1Q = add_line_affine_j((uQ[0], uQ[1]), (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**4 - u**3 - 1
    m_uv, uvQ = miller_function_ate(uQ, P, 0, v, m0=m_u)
    return m_uv * m_u1.frobenius()

def miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_2naf(Q, P, u):
    """
    Return f_{u*(u^4-u^3-1),Q}(P)*frobenius(f_{u+1,Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 0
    trace = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x + r(x) * ht
    v = (u^4-u^3-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate_2naf((Q[0], -Q[1]), P, 0, -u)
    else:
        m_u, uQ = miller_function_ate_2naf(Q, P, 0, u)
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, u1Q = add_line_affine_j((uQ[0], uQ[1]), (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**4 - u**3 - 1 # v > 0 as long as u >= 2
    # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf(uQ, P, 0, v, m0=m_u, m0_inv=inv_m_u)
    return m_uv * m_u1.frobenius()

def miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u(Q,P,u):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{u*(u^4-u^3-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 3
    trace = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3 + r(x) * ht
    v = (u^4-u^3-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate((Q[0], -Q[1]), P, 0, -u)
    else:
        m_u, uQ = miller_function_ate(Q, P, 0, u)
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, u1Q = add_line_affine_j((uQ[0], uQ[1]), (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**4 - u**3 - 1
    m_uv, uvQ = miller_function_ate(uQ, P, 0, v, m0=m_u)
    return m_uv.frobenius() * m_u1

def miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_2naf(Q,P,u,verbose=False):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{u*(u^4-u^3-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 3
    trace = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3 + r(x) * ht
    v = (u^4-u^3-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate_2naf((Q[0], -Q[1]), P, 0, -u)
    else:
        m_u, uQ = miller_function_ate_2naf(Q, P, 0, u)
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, Qu1 = add_line_affine_j((uQ[0], uQ[1]), (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**4 - u**3 - 1 # v is positive as long as u >= 2
    # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf(uQ, P, 0, v, m0=m_u, m0_inv=inv_m_u)
    return m_uv.frobenius() * m_u1

def miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_aklgl(Q,P,b_t,u,Fq6,D_twist=False,xi=None):
    """Return frobenius(f_{u+1,Q}(P)) * f_{u*(u^4-u^3-1),Q}(P)

    INPUT:
    -`Q`: point on the twist of degree 6 defined over Fp
    -`P`: point on E(Fp)
    -`b_t`: coefficient of the twist equation in y^2 = x^3 + b_t
    -`u`: seed
    -`Fp6`: degree 6 extension over a prime field
    -`D_twist`: whether the twisted curve is a D-twist of an M-twist
    -`xi`: the residue in Fp so that Fp6 = Fp[x]/(x^6-xi)

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 0
    trace = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x + r(x) * ht
    v = (u^4-u^3-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique,
    and AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    P_ = (P[0], P[1])
    if u < 0:
        Q_ = (Q[0], -Q[1], 1)
        m_u, uQ = miller_function_ate_aklgl(Q_, P_, b_t, -u, Fq6, D_twist=D_twist)
    else:
        Q_ = (Q[0], Q[1], 1)
        m_u, uQ = miller_function_ate_aklgl(Q_, P_, b_t, u, Fq6, D_twist=D_twist)
    Z1 = 1/uQ[2]
    uQ = (uQ[0]*Z1, uQ[1]*Z1, 1) # homogeneous projective coordinates
    l, u1Q = add_line_h_a0_twist6_aklgl(uQ, (Q[0], Q[1]), P_, D_twist=D_twist)
    if D_twist:
        m_u1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_u, xi, Fq6)
    else:
        m_u1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_u, xi, Fq6)
    v = u**4-u**3-1
    m_uv, uvQ = miller_function_ate_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, D_twist=D_twist)
    return m_u1.frobenius() * m_uv

def miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_aklgl(Q,P,b_t,u,Fq6,D_twist=False,xi=None):
    """Return f_{u+1,Q}(P) * frobenius(f_{u*(u^4-u^3-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 3
    trace = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3 + r(x) * ht
    v = (u^4-u^3-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    and AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    P_ = (P[0], P[1])
    if u < 0:
        Q_ = (Q[0], -Q[1], 1)
        m_u, uQ = miller_function_ate_aklgl(Q_, P_, b_t, -u, Fq6, D_twist=D_twist)
    else:
        Q_ = (Q[0], Q[1], 1)
        m_u, uQ = miller_function_ate_aklgl(Q_, P_, b_t, u, Fq6, D_twist=D_twist)
    Z1 = 1/uQ[2]
    uQ = (uQ[0]*Z1, uQ[1]*Z1, 1) # homogeneous projective coordinates
    l, u1Q = add_line_h_a0_twist6_aklgl(uQ, (Q[0], Q[1]), P_, D_twist=D_twist)
    if D_twist:
        m_u1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_u, xi, Fq6)
    else:
        m_u1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_u, xi, Fq6)
    v = u**4-u**3-1
    m_uv, uvQ = miller_function_ate_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, D_twist=D_twist)
    return m_u1 * m_uv.frobenius()

def miller_loop_opt_ate_bw6_bls24_trace_0_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fq6, D_twist=False, xi=None):
    """Return frobenius(f_{u+1,Q}(P)) * f_{u*(u^4-u^3-1),Q}(P)

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 0
    trace = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x + r(x) * ht
    v = (u^4-u^3-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique,
    AKLGL formulas,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    P_ = (P[0], P[1])
    if u < 0:
        Q_ = (Q[0], -Q[1], 1)
        m_u, uQ = miller_function_ate_2naf_aklgl(Q_, P_, b_t, -u, Fq6, D_twist=D_twist)
    else:
        Q_ = (Q[0], Q[1], 1)
        m_u, uQ = miller_function_ate_2naf_aklgl(Q_, P_, b_t, u, Fq6, D_twist=D_twist)
    Z1 = 1/uQ[2]
    uQ = (uQ[0]*Z1, uQ[1]*Z1, 1) # homogeneous projective coordinates
    l, Qu1 = add_line_h_a0_twist6_aklgl(uQ, (Q[0], Q[1]), P_, D_twist=D_twist)
    if D_twist:
        m_u1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_u, xi, Fq6)
    else:
        m_u1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_u, xi, Fq6)
    v = u**4-u**3-1 # v > 0 as long as u >= 2
    # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, m0_inv=inv_m_u, D_twist=D_twist)
    f = m_u1.frobenius() * m_uv
    return f

def miller_loop_opt_ate_bw6_bls24_trace_3_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fq6, D_twist=False, xi=None):
    """Return f_{u+1,Q}(P) * frobenius(f_{u*(u^4-u^3-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS24
    with trace mod r mod u = 3
    trace = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3 + r(x) * ht
    v = (u^4-u^3-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    AKLGL formulas,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^4-u^3-1 > 0 for any u != 0,1 in ZZ.
    """
    if xi is None:
        #xi = -Fq6.polynomial().constant_coefficient() # works only for absolute extensions on prime fields
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    P_ = (P[0], P[1])
    if u < 0:
        Q_ = (Q[0], -Q[1], 1)
        m_u, uQ = miller_function_ate_2naf_aklgl(Q_, P_, b_t, -u, Fq6, D_twist=D_twist)
    else:
        Q_ = (Q[0], Q[1], 1)
        m_u, uQ = miller_function_ate_2naf_aklgl(Q_, P_, b_t, u, Fq6, D_twist=D_twist)
    Z1 = 1/uQ[2]
    uQ = (uQ[0]*Z1, uQ[1]*Z1, 1) # homogeneous projective coordinates
    l, u1Q = add_line_h_a0_twist6_aklgl(uQ, (Q[0], Q[1]), P_, D_twist=D_twist)
    if D_twist:
        m_u1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_u, xi, Fq6)
    else:
        m_u1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_u, xi, Fq6)
    v = u**4-u**3-1 # v > 0 as long as u >= 2
    # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, m0_inv=inv_m_u, D_twist=D_twist)
    f = m_u1 * m_uv.frobenius()
    return f

####### membership testing, co-factor multiplication
# Formulas from Fuentes-Castaneda, Knapp, Rodriguez-Henriquez
# Faster hashing to G2.
# In: Miri, Vaudenay (eds.) SAC 2011. LNCS, vol. 7118, pp. 412-430.
# Springer, Heidelberg (Aug 2012). https://doi.org/10.1007/978-3-642-28496-0_25

def bw6_bls24_g1_mult_by_3r_trace_0_mod_r_u(P, omega, u0):
    """Return 3*r*P assuming phi(P) = (omega*x, y) = (-(t-1) mod r)*P

    Multiply P by 3*r in order to be able to check if 3*r*P == 0 in E(Fp)
    r = (u-1)^2/3*(u^4-u^2+1) + u
    assume that (trace mod r) = 0 mod u

    with t   = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x
    and  3*r = x^10 - 2*x^9 + x^8 - x^6 + 2*x^5 - x^4 + x^2 + x + 1
    then 3*r = (-x-1)*(t-1) - x^5 + x^4 + x
         3*r = (x+1)*(-t+1) - x^5 + x^4 + x
    and  3*r*P = (x+1)*phi(P) + (-x^5 + x^4 + x)*P
    """
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    x4P = u0*x3P
    x5P = u0*x4P
    Q = xP + P
    phi_Q = bw6_phi(Q, omega)
    R = -x5P + x4P + xP + phi_Q
    return R

def bw6_bls24_g1_check_membership_trace_u(P, omega, u0):
    R = bw6_bls24_g1_mult_by_3r_trace_0_mod_r_u(P, omega, u0)
    return R == P.curve()(0)

def bw6_bls24_g1_mult_by_3r_trace_3_mod_r_u(P, omega, u0):
    """Return 3*r*P assuming phi(P) = (omega*x, y) = (-(t-1) mod r)*P

    Multiply P by 3*r in order to be able to check if 3*r*P == 0 in E(Fp)
    r = (u-1)^2/3*(u^4-u^2+1) + u
    assume that (trace mod r) = 3 mod u

    with t   = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3
    and  3*r = x^10 - 2*x^9 + x^8 - x^6 + 2*x^5 - x^4 + x^2 + x + 1
    then 3*r = -x^5 + x^4 - 1 + (x+1)*(t-1)
         3*r = -x^5 + x^4 - 1 - (x+1)*(-t+1)
    and 3*r*P = (-x^5 + x^4 - 1)*P - (x+1)*phi(P)
    """
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    x4P = u0*x3P
    x5P = u0*x4P
    Q = xP + P
    phi_Q = bw6_phi(Q, omega)
    R = -x5P + x4P - P - phi_Q
    return R

def bw6_bls24_g1_check_membership_trace_3u(P, omega, u0):
    R = bw6_bls24_g1_mult_by_3r_trace_3_mod_r_u(P, omega, u0)
    return R == P.curve()(0)

def bw6_bls24_g1_mult_by_cofactor_trace_0_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t'-2) + ((u-1)^2*(u^2+1)*((u-1)^2*(u^2+1)+1) + 1 - hy
    where t' = t mod r = -u^9 + 3*u^8 - 4*u^7 + 4*u^6 - 3*u^5 + 2*u^3 - 2*u^2 + u

    Assume that eigenvalue of endomorphism phi(x,y) = (omega*x,y) is exactly -(t-1) mod r.
    Phi with (-omega-1) of eigenvalue (t-2) gives the function
    bw6_bls24_g1_mult_by_cofactor_trace_0_mod_r_u_alt

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)*(-u-1) + 2*(ht+3*hy)*(u-1)^2*(u^2+1) + 4*ht
    l1 = (ht^2+3*hy^2)*u*(u^4-u^3-1) - 2*(ht-3*hy)*(u-1)^2*(u^2+1) + 2*(ht+3*hy)

    alternative formulas:
    l0 = (u^5-u^4+1)*(ht^2+3*hy^2) - 4*ht*(u-1)^2*(u^2+1) - 2*(ht-3*hy)
    l1 = (ht^2+3*hy^2)*(u+1) - 2*(ht+3*hy)*(u-1)^2*(u^2+1) -4*ht

    (u-1)^2*(u^2+1) = u^4 - 2*u^3 + 2*u^2 - 2*u + 1
    (u^5-u^4+1) == (u-1)^2*(u^2+1)*(u+1) + u
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    vP = (u0**2+1)*(u0-1)*(uP-P) # (u-1)^2*(u^2+1)*P
    wP = (u0+1)*vP + uP
    L0 = d1 * wP - ht*vP - d3*P
    L1 = d1*(uP+P) - d2*vP - ht*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls24_g1_mult_by_cofactor_trace_0_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t'-2) + ((u-1)^2*(u^2+1)*((u-1)^2*(u^2+1)+1) + 1 - hy
    where t' = t mod r = -u^9 + 3*u^8 - 4*u^7 + 4*u^6 - 3*u^5 + 2*u^3 - 2*u^2 + u

    Assume that eigenvalue of endomorphism phi(x,y) = (omega'*x,y) is exactly (t-2) mod r.
    Phi with the other choice of omega of eigenvalue -(t-1) gives the above function
    bw6_bls24_g1_mult_by_cofactor_trace_0_mod_r_u

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)/4*(u+1) -(ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) - ht
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4+1) -   ht*(u^4-2*u^3+2*u^2-2*u+1) - (ht-3*hy)/2
    u^4-2*u^3+2*u^2-2*u+1 == (u-1)^2*(u^2+1)
    u^5-u^4+1 == (u+1)*(u^4-2*u^3+2*u^2-2*u+1) + u
    equivalently
    l0 = (ht^2+3*hy^2)/4*(u^5-u^4-u) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) + ht
    l1 = -(ht^2+3*hy^2)/4*(u+1) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    vP = (u0**2+1)*(u0-1)*(uP-P) # (u-1)^2*(u^2+1)*P
    wP = (u0+1)*vP + uP
    L0 = d1*(uP+P) - d2*vP - ht*P
    L1 = d1 * wP - ht*vP - d3*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls24_g1_mult_by_cofactor_trace_3_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    cx = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t'-2) + ((u-1)^2*(u^2+1)*((u-1)^2*(u^2+1)+1) + 1 + hy
    where t' = t mod r = u^9 - 3*u^8 + 4*u^7 - 4*u^6 + 3*u^5 - 2*u^3 + 2*u^2 - u + 3

    Assume that eigenvalue of endomorphism phi(x,y) = (omega*x,y) is exactly -(t-1) mod r.
    Phi with (-omega-1) of eigenvalue (t-2) gives the function
    bw6_bls24_g1_mult_by_cofactor_trace_3_mod_r_u_alt

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) - (ht^2+3*hy^2)/4*(u+1) - ht
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4-u) + ht*(u^4-2*u^3+2*u^2-2*u+1) + (ht-3*hy)/2

    (u^4-2*u^3+2*u^2-2*u+1) == (u-1)^2*(u^2+1)
    (u^5-u^4-u) == (u-1)^2*(u^2+1) * (u+1) - 1

    equivalent:
    l0 = (ht^2+3*hy^2)/4*(u^5-u^4+1) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) + ht
    l1 = (ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    vP = (u0**2+1)*(u0-1)*(uP-P) # (u-1)^2*(u^2+1)*P
    wP = (u0+1)*vP + uP          # (u^5-u^4+1)*P
    htP = ht*P
    L0 = d1*wP + d2*vP + htP
    L1 = d1*(uP+P) - d3*(vP+P) + htP
    return L0 + bw6_phi(L1, omega)

def bw6_bls24_g1_mult_by_cofactor_trace_3_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    cx = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t'-2) + ((u-1)^2*(u^2+1)*((u-1)^2*(u^2+1)+1) + 1 + hy
    where t' = t mod r = u^9 - 3*u^8 + 4*u^7 - 4*u^6 + 3*u^5 - 2*u^3 + 2*u^2 - u + 3

    Assume that eigenvalue of endomorphism phi(x,y) = (omega'*x,y) is exactly (t-2) mod r.
    Phi with the other choice of omega of eigenvalue -(t-1) gives the above function
    bw6_bls24_g1_mult_by_cofactor_trace_3_mod_r_u

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) + ht
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4+1) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) + ht
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^5-u^4-u) + ht*(u^4-2*u^3+2*u^2-2*u+1) + (ht-3*hy)/2
    l1 = -(ht^2+3*hy^2)/4*(u+1) + (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) - ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    vP = (u0**2+1)*(u0-1)*(uP-P) # (u-1)^2*(u^2+1)*P
    wP = (u0+1)*vP + uP          # (u^5-u^4+1)*P
    L0 = d1*(uP+P) - d3*vP + d2*P
    L1 = d1*wP + d2*vP + ht*P
    return L0 + bw6_phi(L1, omega)


def bw6_bls24_g2_mult_by_cofactor_trace_0_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    Assume that eigenvalue of endomorphism phi(x,y) = (omega*x,y) is exactly -(t'-1) mod r.
    The other pair is Phi with (-omega-1) of eigenvalue (t'-2) but gives other formulas.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)/4*(u+1) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) -ht
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4+1) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) - ht

    l0 = (ht^2+3*hy^2)/4*(u^5-u^4-u) - ht*(u^4-2*u^3+2*u^2-2*u+2) + (ht-3*hy)/2
    l1 = -(ht^2+3*hy^2)/4*(u+1) - (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht-3*hy)//2
    d3 = (ht+3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    u4P = u0*u3P
    u5P = u0*u4P
    vP = u4P+P -2*(u3P-u2P+uP)
    wP = u5P-u4P+P
    xP = vP + P
    htP = ht*P
    L0 = d1*(uP+P) + d3*xP - htP
    L1 = d1*wP - d2*vP -ht*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls24_g2_mult_by_cofactor_trace_0_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    Assume that eigenvalue of endomorphism phi(x,y) = (omega'*x,y) is exactly (t'-2) mod r.
    The other pair is Phi with omega = (-omega'-1) of eigenvalue -(t'-1) but gives other formulas.

    l0 + l1*eigenvalue_phi_mod_cx_alt = 0 mod cx

    l0 = -(ht^2+3*hy^2)/4*(u+1) - (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) + ht;
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4-u) - ht*(u^4-2*u^3+2*u^2-2*u+1) - (ht+3*hy)/4;
    alt
    l0 = (ht^2+3*hy^2)/4*(u^5-u^4+1) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) - ht;
    l1 = (ht^2+3*hy^2)/4*(u+1) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) - ht;

    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht-3*hy)//2
    d3 = (ht+3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    u4P = u0*u3P
    u5P = u0*u4P
    vP = u4P+P -2*(u3P-u2P+uP)
    wP = u5P-u4P-uP
    xP = vP + P
    htP = ht*P
    L0 = -d1 * (uP+P) -d3*xP + htP
    L1 = d1*wP -ht*vP - d3*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls24_g2_mult_by_cofactor_trace_3_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    Assume that eigenvalue of endomorphism phi(x,y) = (omega*x,y) is exactly -(t'-1) mod r.
    The other pair is Phi with (-omega-1) of eigenvalue (t'-2) but gives other formulas.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = -(ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) - ht
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4-u) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) - ht

    alt:
    l0 = (ht^2+3*hy^2)/4*(u^5-u^4+1) + ht*(u^4-2*u^3+2*u^2-2*u+1) + (ht+3*hy)/2
    l1 = (ht^2+3*hy^2)/4*(u+1) + (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    u4P = u0*u3P
    u5P = u0*u4P
    vP = u4P+P -2*(u3P-u2P+uP)
    wP = u5P-u4P-uP
    xP = vP + P
    htP = ht*P
    L0 = -d1 * (uP+P) -d3*vP - htP
    L1 = d1*wP + d2*xP -htP
    return L0 + bw6_phi(L1, omega)

def bw6_bls24_g2_mult_by_cofactor_trace_3_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    Assume that eigenvalue of endomorphism phi(x,y) = (omega'*x,y) is exactly (t'-2) mod r.
    The other pair is Phi with omega = (-omega'-1) of eigenvalue -(t'-1) but gives other formulas.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)/4*(u+1) + (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) + ht
    l1 = (ht^2+3*hy^2)/4*(u^5-u^4+1) + ht*(u^4-2*u^3+2*u^2-2*u+1) + (ht + 3*hy)/2

    l0 = (ht^2+3*hy^2)/4*(u^5-u^4-u) + (ht+3*hy)/2*(u^4-2*u^3+2*u^2-2*u+2) - ht
    l1 = -(ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^4-2*u^3+2*u^2-2*u+1) - ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    u4P = u0*u3P
    u5P = u0*u4P
    vP = u4P+P -2*(u3P-u2P+uP)
    wP = u5P-u4P+P
    htP = ht*P
    L0 = d1*(uP+P) + d3*vP + htP
    L1 = d1*wP + ht*vP +d2*P
    return L0 + bw6_phi(L1, omega)

def final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_subfunction(B, u, ht, hy):
    """
    return B^((ht^2+3*hy^2)/4*r + 3*((ht+hy)/2*t0/3 + (d-1)/3) + 1)

    INPUT:
    - `B`: a finite field element, with B.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t0
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t0^2-4*q = -D*y^2

    RETURN: B^e

    parameters are
    r = (x-1)^2/3*(x^8-x^4+1) + x
    t0 = -x^9 + 3*x^8 - 4*x^7 + 4*x^6 - 3*x^5 + 2*x^3 - 2*x^2 + x
    d = x^8-4*x^7+8*x^6-12*x^5+15*x^4-14*x^3+10*x^2-6*x+3
    (d-1)/3 == (x - 1)^2/3 * (x^2 + 1) * ((x - 1)^2 * (x^2 + 1) + 1)
    t0/3 == (-x-1)*(d-1)/3 + (x - 1)^2/3 * (x^2 + 1) - (x-1)/3
    r == -(x+1)*(t0/3 + (x - 1)^2/3 * (x^2 + 1)) + (x-1)/3 + 1

    a sequence to obtain the exponents r, t0/3 and (d-1)/3 is
    C = (x-1)/3
    D = C * (x-1) * (x^2+1)
    E = D * ((x - 1)^2 * (x^2 + 1) + 1)
    E == (d-1)/3
    F = -(x+1)*E + D - C
    F == t0/3
    G = -(x+1)*(F + D) + C + 1
    G == r

    cost
    e((u-1)/3) + 3e(u-1) + 2min(2e(u)+M, e(u^2+1)) + 2e(u+1) + e(h1) + e(h2) + 10 M + S + 2cj
    note that (A**u)**u * A == A**(u**2+1),
    if (u^2+1) has a lower Hamming weight than 2*HW(u)+1, do A**(u**2+1) instead of (A**u)**u * A
    """
    C = B**((u-1)//3)
    D = C**(u-1)
    D = D**(u**2+1)                    # or (D**u)**u * D
    E = (D**(u-1))**(u-1)
    E = E**(u**2+1) * D                # or ((E**u)**u * E)
                                       # E = B^((d-1)/3)
    F = (E**(u+1) * C).conjugate() * D # B^(t0/3)
    G = (F * D).conjugate()
    H = G**(u+1) * C * B    # B^r

    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+hy)//2

    I = F**d1 * E
    I = I**2 * I * B * H**d2

    return I

def final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u(m, u, ht, hy):
    """
    Exponentiation to (u+1)*Phi_6(q(u))/r(u)

    INPUT:
    - `m`: a finite field element, with m.frobenius() and m.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t0
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t0^2-4*q = -D*y^2

    RETURN: m^((u+1)*Phi_6(q)/r)

    similar to final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u(m, u, ht, hy)
    see also   final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt

    parameters are
    QQh.<ht,hy,q> = QQ[]
    QQu.<u> = QQh[]
    cc0 = u^8 - 4*u^7 + 8*u^6 - 12*u^5 + 15*u^4 - 14*u^3 + 10*u^2 - 6*u + 3
    tr0 = -u^9 + 3*u^8 - 4*u^7 + 4*u^6 - 3*u^5 + 2*u^3 - 2*u^2 + u
    r   = (u^10 - 2*u^9 + u^8 - u^6 + 2*u^5 - u^4 + u^2 + u + 1)/3
    e0 = -(ht^2+3*hy^2)/4*r*(u^5-u^4+1) - (ht-hy)/2*tr0*(u^5-u^4+1) - cc0*(u^5-u^4+1) + 3*(u^4-2*u^3+2*u^2-2*u+2)
    e1 = (ht^2+3*hy^2)/4*r*(u+1) + (ht-hy)/2*tr0*(u+1) + cc0*(u+1) - 3
    e0 + e1*q == (-(u^5-u^4+1) + (u+1)*q)*((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0) - 3*(-(u^4-2*u^3+2*u^2-2*u+2) + q)
    #with
    (u^5-u^4+1) == (u^4-2*u^3+2*u^2-2*u+2)*(u+1) - 1
    (u^5-u^4+1) == ((u - 1)^2 * (u^2 + 1) + 1)*(u+1) - 1
    (-(u^5-u^4+1) + (u+1)*q) == (-((u-1)^2*(u^2+1)+1) + q)*(u+1) + 1
    e0 + e1*q == (1 + (u+1)*(-((u-1)^2*(u^2+1)+1) + q))*((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0) - 3*(-((u-1)^2*(u^2+1)+1) + q)
    c0 = (u-1)^2*(u^2+1)+1
    e0 + e1*q == (1 + (u+1)*(-c0 + q))*((ht^2+3*hy^2)/4*r + 3*((ht-hy)/2*tr0/3 + (cc0-1)/3) + 1) - 3*(-c0 + q)

    cost subfunction e((u-1)/3) + 3e(u-1) + 2min(2e(u)+M, e(u^2+1)) + 2e(u+1) + e(h1) + e(h2) + 10M +  S       + 2cj
    cost here                     2e(u-1) +  min(2e(u)+M, e(u^2+1)) +  e(u+1)                 +  5M +  S + frb + 2cj
    total cost       e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 15M + 2S + frb + 4cj
    """
    A = m**(u-1)
    A = A**(u-1)
    A = A**(u**2+1)
    #A = (A**u)**u * A                 # if 2e(u)+M < e(u^2+1)
    A = A * m
    A = A.conjugate() * m.frobenius()
    B = A**(u+1) * m                   # m**(1 + (u+1)*(-c0 + q))
    A = (A**2 * A).conjugate()         # m**(-3*(-c0 + q))

    C = final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_subfunction(B, u, ht, hy)

    return A * C

def final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_alt(m, u, ht, hy):
    """
    Exponentiation to (u^5-u^4-u)*Phi_6(q(u))/r(u)

    INPUT:
    - `m`: a finite field element, with m.frobenius() and m.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t0
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t0^2-4*q = -D*y^2

    RETURN: m^((u^5-u^4-u)*Phi_6(q)/r)

    see also final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u

    parameters are
    QQh.<ht,hy,q> = QQ[]
    QQu.<u> = QQh[]
    cc0 = u^8 - 4*u^7 + 8*u^6 - 12*u^5 + 15*u^4 - 14*u^3 + 10*u^2 - 6*u + 3
    tr0 = -u^9 + 3*u^8 - 4*u^7 + 4*u^6 - 3*u^5 + 2*u^3 - 2*u^2 + u
    r   = (u^10 - 2*u^9 + u^8 - u^6 + 2*u^5 - u^4 + u^2 + u + 1)/3
    e0 = (ht^2+3*hy^2)/4*r*(u+1) + (ht-hy)/2*tr0*(u+1) + cc0*(u+1) - 3
    e1 = (ht^2+3*hy^2)/4*r*(u^5-u^4-u) + (ht-hy)/2*tr0*(u^5-u^4-u) + cc0*(u^5-u^4-u) - 3*(u^4-2*u^3+2*u^2-2*u+1)
    e0 + e1*q == ((u+1) + (u^5-u^4-u)*q)*((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0) - 3*(1 + (u^4-2*u^3+2*u^2-2*u+1)*q)
    #with
    (u^5-u^4-u) == (u^4-2*u^3+2*u^2-2*u+1)*(u+1) - 1
    ((u+1) + (u^5-u^4-u)*q) == (1 + (u^4-2*u^3+2*u^2-2*u+1)*q)*(u+1) - q
    u^4-2*u^3+2*u^2-2*u+1 == (u - 1)^2 * (u^2 + 1)

    c0 = (u-1)^2*(u^2+1)
    e0 + e1*q == ((u+1)*(1 + c0*q) - q)*((ht^2+3*hy^2)/4*r + 3*((ht-hy)/2*tr0/3 + (cc0-1)/3) + 1) - 3*(1 + c0*q)

    cost subfunction e((u-1)/3) + 3e(u-1) + 2min(2e(u)+M, e(u^2+1)) + 2e(u+1) + e(h1) + e(h2) + 10M +  S       + 2cj
    cost here                     2e(u-1) +  min(2e(u)+M, e(u^2+1)) +  e(u+1)                 +  4M +  S + frb + 2cj
    total cost       e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 14M + 2S + frb + 4cj
    """
    mq = m.frobenius()
    A = mq**(u-1)
    A = A**(u-1)
    A = A**(u**2+1)
    #A = (A**u)**u * A            # if 2e(u)+M < e(u^2+1)
    A = m * A
    B = A**(u+1) * mq.conjugate() # m**((u+1)*(1 + c0*q) - q)
    A = (A**2 * A).conjugate()    # m**(- 3*(1+ c0*q))

    C = final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u_subfunction(B, u, ht, hy)

    return A * C

def final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_subfunction(B, u, ht, hy):
    """
    return B^((ht^2+3*hy^2)/4*r + 3*((ht-hy)/2*t3/3 + (d-1)/3) + 1)

    INPUT:
    - `B`: a finite field element, with B.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t3
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t3^2-4*q = -D*y^2

    RETURN: B^e

    parameters are
    r = (x-1)^2/3*(x^8-x^4+1) + x
    t3 = x^9 - 3*x^8 + 4*x^7 - 4*x^6 + 3*x^5 - 2*x^3 + 2*x^2 - x + 3
    d = x^8-4*x^7+8*x^6-12*x^5+15*x^4-14*x^3+10*x^2-6*x+3
    (d-1)/3 == (x - 1)^2/3 * (x^2 + 1) * ((x - 1)^2 * (x^2 + 1) + 1)
    t3/3 == (x+1)*(d-1)/3 - (x - 1)^2/3 * (x^2 + 1) + (x-1)/3 + 1
    r == (x+1)*(t3/3 - 1 - (x - 1)^2/3 * (x^2 + 1)) + (x-1)/3 + 1
    r == (x+1)*(t3/3 - (x-1)/3 - (x - 1)^2/3 * (x^2 + 1)) + (x-1)^2/3 - 1

    a sequence to obtain the exponents r, t3/3 and (d-1)/3 is
    C = (x-1)/3
    D = C * (x-1) * (x^2+1)
    E = D * ((x - 1)^2 * (x^2 + 1) + 1)
    E == (d-1)/3
    D = -D
    F = (x+1)*E + D + C
    G = F + 1
    G == t3/3
    J = (x+1)*(F + D) + C + 1
    J == r

    cost
    e((u-1)/3) + 3e(u-1) + 2min(2e(u)+M, e(u^2+1)) + 2e(u+1) + e(h1) + e(h2) + 11 M + S + cj
    note that (A**u)**u * A == A**(u**2+1),
    if (u^2+1) has a lower Hamming weight than 2*HW(u)+1, do A**(u**2+1) instead of (A**u)**u * A
    """
    C = B**((u-1)//3)
    D = C**(u-1)
    D = D**(u**2+1)               # or (D**u)**u * D
    E = (D**(u-1))**(u-1)
    E = E**(u**2+1) * D           # or ((E**u)**u * E) * D
                                  # E = B^((d-1)/3)
    D = D.conjugate()
    F = E**(u+1) * D * C
    G = F * B                     # B^(t3/3)
    H = (F * D)**(u+1) * C * B    # B^r

    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+hy)//2

    I = G**d1 * E
    I = I**2 * I * B * H**d2

    return I

def final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u(m, u, ht, hy):
    """
    exponentiation to (u+1)*Phi_6(q(u))/r(u)

    INPUT:
    - `m`: a finite field element, with m.frobenius() and m.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t3
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t3^2-4*q = -D*y^2

    RETURN: m^((u+1)*Phi_6(q)/r)

    similar to final_exp_hard_bw6_bls24_trace_0_mod_r_mod_u(m, u, ht, hy)
    see also   final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt

    parameters are
    a = (u-1)^2*(u^2+1)
    ((u+1)*(a+q)-1)*(c3+ht) + 3*(a+q)

    cost subfunction e((u-1)/3) + 3e(u-1) + 2min(2e(u)+M, e(u^2+1)) + 2e(u+1) + e(h1) + e(h2) + 11M +  S       +  cj
    cost here                     2e(u-1) +  min(2e(u)+M, e(u^2+1)) +  e(u+1)                 +  4M +  S + frb +  cj
    total cost       e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 15M + 2S + frb + 2cj
    """
    A = m**(u-1)
    A = A**(u-1)
    A = A**(u**2+1)
    #A = (A**u)**u * A            # if 2e(u)+M < e(u^2+1)
    A = A * m.frobenius()
    B = A**(u+1) * m.conjugate()
    A = A**2 * A

    C = final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_subfunction(B, u, ht, hy)

    return C * A

def final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_alt(m, u, ht, hy):
    """
    Exponentiation to (u^5-u^4+1)*Phi_6(q(u))/r(u)

    INPUT:
    - `m`: a finite field element, with m.frobenius() and m.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t3
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t3^2-4*q = -D*y^2

    RETURN: m^((u^5-u^4+1)*Phi_6(q)/r)

    similar to final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u(m, u, ht, hy)

    parameters are
    QQh.<ht,hy,q> = QQ[]
    QQu.<u> = QQh[]
    cc3 = u^8 - 4*u^7 + 8*u^6 - 12*u^5 + 15*u^4 - 14*u^3 + 10*u^2 - 6*u + 3
    tr3 = u^9 - 3*u^8 + 4*u^7 - 4*u^6 + 3*u^5 - 2*u^3 + 2*u^2 - u + 3
    r   = (u^10 - 2*u^9 + u^8 - u^6 + 2*u^5 - u^4 + u^2 + u + 1)/3
    cc3_1_3 = (cc3-1)/3
    tr3_3 = tr3/3

    e0 = (-u - 1)*((ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3) - 3
    e1 = (u^5 - u^4 + 1)*((ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3) + 3*(u^4 - 2*u^3 + 2*u^2 - 2*u + 2)

    e0 + e1*q == ((ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3)*((-u - 1)+(u^5 - u^4 + 1)*q) + 3*(-1 + (u^4 - 2*u^3 + 2*u^2 - 2*u + 2)*q)
    #where
    (u^5 - u^4 + 1) == (u^4-2*u^3+2*u^2-2*u+2)*(u+1) - 1

    e0 + e1*q == ((u+1)*(-1 + c0*q) - q)*((ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3) + 3*(-1 + c0*q)

    cost subfunction e((u-1)/3) + 3e(u-1) + 2min(2e(u)+M, e(u^2+1)) + 2e(u+1) + e(h1) + e(h2) + 11M +  S       +  cj
    cost here                     2e(u-1) +  min(2e(u)+M, e(u^2+1)) +  e(u+1)                 +  5M +  S + frb + 2cj
    total cost       e((u-1)/3) + 5e(u-1) + 3min(2e(u)+M, e(u^2+1)) + 3e(u+1) + e(h1) + e(h2) + 16M + 2S + frb + 3cj
    """
    mq = m.frobenius()
    A = mq**(u-1)
    A = A**(u-1)
    A = A**(u**2+1)
    #A = (A**u)**u * A            # if 2e(u)+M < e(u^2+1)
    A = A * mq
    A = m.conjugate() * A         # A = m**(q*((u-1)^2*(u^2+1) + 1) - 1)
    B = A**(u+1) * mq.conjugate() # m**(-(u+1)*(1 - c0*q) - q)
    A = (A**2 * A) # m**(3*(-1+ c0*q))

    C = final_exp_hard_bw6_bls24_trace_3_mod_r_mod_u_subfunction(B, u, ht, hy)

    return A * C
