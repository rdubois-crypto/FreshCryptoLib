from external.Pairings.pairing import *

def miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u(Q, P, u):
    """
    Return f_{u*(u^2-u-1),Q}(P)*frobenius(f_{u+1,Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 0
    trace = -x^5 + 3*x^4 - 3*x^3 + x + r(x) * ht
    v = (u^2-u-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2 - u - 1 # v is positive as long as u >= 2
    m_uv, uvQ = miller_function_ate(uQ, P, 0, v, m0=m_u)
    return m_uv * m_u1.frobenius()

def miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_2naf(Q, P, u):
    """
    Return f_{u*(u^2-u-1),Q}(P)*frobenius(f_{u+1,Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 0
    trace = -x^5 + 3*x^4 - 3*x^3 + x + r(x) * ht
    v = (u^2-u-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2 - u - 1 # v > 0 as long as u >= 2
        # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf(uQ, P, 0, v, m0=m_u, m0_inv=inv_m_u)
    return m_uv * m_u1.frobenius()

def miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u(Q,P,u):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 3
    trace = x^5 - 3*x^4 + 3*x^3 - x + 3 + r(x) * ht
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2 - u - 1
    m_uv, uvQ = miller_function_ate(uQ, P, 0, v, m0=m_u)
    return m_uv.frobenius() * m_u1

def miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_alt(Q,P,u):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{(u+1)*(u^2-2*u+1),[u+1]Q}(P) * l_{(u+1)*(u^2-2*u+1)Q,-Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 3
    trace = x^5 - 3*x^4 + 3*x^3 - x + 3 + r(x) * ht
    v = (u^2-2*u+1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u+1,Q}(P)*f_{v,[u+1]Q}(P) * l_{(u+1)*vQ,-Q}(P))
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-2u+1 = (u-1)^2 > 0 for any u != 0 in ZZ.
    """
    u1 = u+1
    if u1 < 0:
        m_u1, u1Q = miller_function_ate((Q[0], -Q[1]), P, 0, -u1)
    else:
        m_u1, u1Q = miller_function_ate(Q, P, 0, u1)
    Z1 = 1/u1Q[2]
    Z2 = Z1**2
    u1Q = (u1Q[0]*Z2, u1Q[1]*Z1*Z2, 1, 1)
    v = u**2 - 2*u + 1
    m_u1v, u1vQ = miller_function_ate(u1Q, P, 0, v, m0=m_u1)
    l, u1v1Q = add_line_j(u1vQ, (Q[0], -Q[1]), (P[0], P[1]))
    m_u1v1 = m_u1v * l
    return m_u1v1.frobenius() * m_u1

def miller_loop_opt_ate_bw6_761(Q,P,u):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6-761
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    return miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u(Q,P,u)

def miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_2naf(Q,P,u):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 3
    trace = x^5 - 3*x^4 + 3*x^3 - x + 3 + r(x) * ht
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
    """
    if u < 0:
        m_u, uQ = miller_function_ate_2naf((Q[0], -Q[1]), P, 0, -u)
    else:
        m_u, uQ = miller_function_ate_2naf(Q, P, 0, u)
    # uQ in affine coordinates is required anyway for the second miller loop
    Z1 = 1/uQ[2]
    Z2 = Z1**2
    uQ = (uQ[0]*Z2, uQ[1]*Z1*Z2, 1, 1)
    l, u1Q = add_line_affine_j((uQ[0], uQ[1]), (Q[0], Q[1]), (P[0], P[1]))
    m_u1 = m_u * l
    v = u**2 - u - 1 # v > 0 as long as u >= 2
    # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf(uQ, P, 0, v, m0=m_u, m0_inv=inv_m_u)
    return m_uv.frobenius() * m_u1

def miller_loop_opt_ate_bw6_761_2naf(Q, P, u):
    """
    Return f_{u+1,Q}(P)*frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6-761
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    return miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_2naf(Q, P, u)

def miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_aklgl(Q,P,b_t,u,Fq6,D_twist=False,xi=None):
    """Return frobenius(f_{u+1,Q}(P)) * f_{u*(u^2-u-1),Q}(P)

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 0
    trace = -x^5 + 3*x^4 - 3*x^3 + x + r(x) * ht
    v = (u^2-u-1)
    as frobenius(f_{u+1,Q}(P)) * f^v_{u,Q}(P)*f_{v,[u]Q}(P)
    with a multi-exponentiation-like technique,
    and AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2-u-1
    m_uv, uvQ = miller_function_ate_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, D_twist=D_twist)
    return m_u1.frobenius() * m_uv

def miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl(Q,P,b_t,u,Fq6,D_twist=False,xi=None):
    """Return f_{u+1,Q}(P) * frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12 (incl. BW6_761)
    with trace mod r mod u = 3
    trace = x^5 - 3*x^4 + 3*x^3 - x + 3 + r(x) * ht
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    and AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2-u-1
    m_uv, uvQ = miller_function_ate_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, D_twist=D_twist)
    return m_u1 * m_uv.frobenius()

def miller_loop_opt_ate_bw6_761_aklgl(Q, P, b_t, u, Fq6, D_twist=False, xi=None):
    """Return f_{u+1,Q}(P) * frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6-761
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique
    and AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    return miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl(Q,P,b_t,u,Fq6,D_twist=D_twist,xi=xi)

def miller_loop_opt_ate_bw6_bls12_trace_0_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fq6, D_twist=False, xi=None):
    """Return frobenius(f_{u+1,Q}(P)) * f_{u*(u^2-u-1),Q}(P)

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 0
    trace = -x^5 + 3*x^4 - 3*x^3 + x + r(x) * ht
    v = (u^2-u-1)
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

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2-u-1 # v > 0 as long as u >= 2
    # compute the inverse of m_u with a conjugate m_u^(p^3)
    try:
        inv_m_u = m_u.conjugate()
    except AttributeError as err:
        print(err)
        inv_m_u = frobenius_fp6_p3(m_u)
    m_uv, uvQ = miller_function_ate_2naf_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, m0_inv=inv_m_u, D_twist=D_twist)
    f = m_u1.frobenius() * m_uv
    return f

def miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fq6, D_twist=False, xi=None):
    """Return f_{u+1,Q}(P) * frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6_BLS12
    with trace mod r mod u = 3
    trace = x^5 - 3*x^4 + 3*x^3 - x + 3 + r(x) * ht
    v = (u^2-u-1)
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

    Note that v = u^2-u-1 > 0 for any u != 0,1 in ZZ.
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
    v = u**2-u-1
    m_uv, uvQ = miller_function_ate_2naf_aklgl(uQ, P_, b_t, v, Fq6, m0=m_u, D_twist=D_twist)
    f = m_u1 * m_uv.frobenius()
    return f

def miller_loop_opt_ate_bw6_761_aklgl_2naf(Q, P, b_t, u, Fq6, D_twist=False, xi=None):
    """Return f_{u+1,Q}(P) * frobenius(f_{u*(u^2-u-1),Q}(P))

    Optimized optimal ate Miller loop for BW6-761
    v = (u^2-u-1)
    as f_{u+1,Q}(P) * frobenius(f^v_{u,Q}(P)*f_{v,[u]Q}(P))
    with a multi-exponentiation-like technique,
    AKLGL formulas,
    and 2-NAF representation of scalars
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients
    """
    return miller_loop_opt_ate_bw6_bls12_trace_3_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fq6, D_twist=D_twist, xi=xi)

def final_exp_easy_bw6_761(m):
    return final_exp_easy_k6(m)

def final_exp_hard_bw6_761(m, u0):
    """
    m^(3*(u^3-u^2+1)*(p^2-p+1)/r)
    """
    f1 = m**u0
    f2 = f1**u0
    f3 = f2**u0
    f4 = f3**u0
    f5 = f4**u0
    f6 = f5**u0
    f7 = f6**u0
    fR0 = (m**220*f1**263*f2**73*f3**314*f4**197*f7**103).frobenius(3) * f5**269*f6**70
    fp = m.frobenius()
    f1p = f1.frobenius()
    f2p = f2.frobenius()
    f3p = f3.frobenius()
    f4p = f4.frobenius()
    f5p = f5.frobenius()
    f6p = f6.frobenius()
    f7p = f7.frobenius()
    f8p = f7p**u0
    f9p = f8p**u0
    fR1 = fp**229*f1p**34*f3p**452*f6p**492*f7p**77*f9p**103 * (f2p**181*f4p**65*f5p**445*f8p**276).frobenius(3)
    ff = fR0*fR1
    return ff

def final_exp_hard_bw6_761_multi_2naf(m,u0):
    """ returns m^exponent where exponent is a multiple of (q^2-q+1)/r"""
    f0 = m
    f1 = m**u0
    f2 = f1**u0
    f3 = f2**u0
    f4 = f3**u0
    f5 = f4**u0
    f6 = f5**u0
    f7 = f6**u0
    f0p = m.frobenius()
    f1p = f1.frobenius()
    f2p = f2.frobenius()
    f3p = f3.frobenius()
    f4p = f4.frobenius()
    f5p = f5.frobenius()
    f6p = f6.frobenius()
    f7p = f7.frobenius()
    f8p = f7p**u0
    f9p = f8p**u0

    f = f3p*f6p*(f5p).frobenius(3)                         # 2M
    f = f**2
    f4f2p = f4*f2p                                         # 1M
    f *= f5*f0p*(f0*f1*f3*f4f2p*f8p).frobenius(3)          # 7M
    f = f**2
    f *= f9p*(f7).frobenius(3)                             # 2M
    f = f**2
    f2f4p = f2*f4p                                         # 1M
    f4f2pf5p = f4f2p*f5p                                   # 1M
    f *= f4f2pf5p*f6*f7p*(f2f4p*f3*f3p).frobenius(3)       # 6M
    f = f**2
    f *= f0*f7*f1p*(f0p*f9p).frobenius(3)                  # 5M
    f = f**2
    f6pf8p = f6p*f8p                                       # 1M
    f5f7p = f5*f7p                                         # 1M
    f *= f5f7p*f2p*(f6pf8p).frobenius(3)                   # 3M
    f = f**2
    f3f6 = f3*f6                                           # 1M
    f1f7 = f1*f7                                           # 1M
    f *= f3f6*f9p*(f1f7*f2).frobenius(3)                   # 4M
    f = f**2
    f *= f0*f0p*f3p*f5p*(f4f2p*f5f7p*f6pf8p).frobenius(3)  # 7M
    f = f**2
    f *= f1p*(f3f6).frobenius(3)                           # 2M
    f = f**2
    f *= f1f7*f5f7p*f0p*(f2f4p*f4f2pf5p*f9p).frobenius(3)  # 6M
    # 51 M
    #  9 S
    # 9 exponentiations to u0
    # 8 frobenius (cost 4 m)
    # 10 frobenius(3) (cost 3 negations)
    return f

def final_exp_bw6_761(m, u0):
    f = final_exp_easy_bw6_761(m)
    #f = final_exp_hard_bw6_761(f, u0)
    f = final_exp_hard_bw6_761_multi_2naf(f, u0)
    return f

def ate_pairing_bw6_761(Q, P, u0, tr):
    T = tr-1
    m,S1 = miller_function_ate(Q, P, 0, T)
    f = final_exp_bw6_761(m, u0)
    return f

def ate_pairing_bw6_761_csb(Q, P, u0, tr):
    T = tr-1
    m,S1 = miller_function_ate_csb(Q, P, 0, T)
    f = final_exp_bw6_761(m, u0)
    return f

def ate_pairing_bw6_761_2naf(Q, P, u0, tr):
    T = tr-1
    m,S1 = miller_function_ate_2naf(Q, P, 0, T)
    f = final_exp_bw6_761(m, u0)
    return f

def tate_pairing_bw6_761(P, Q, u0, r):
    m,S1 = miller_function_tate(P, Q, 0, r)
    f = final_exp_bw6_761(m, u0)
    return f

def tate_pairing_bw6_761_csb(P, Q, u0, r):
    m,S1 = miller_function_tate_csb(P, Q, 0, r)
    f = final_exp_bw6_761(m, u0)
    return f

def tate_pairing_bw6_761_2naf(P, Q, u0, r):
    m,S1 = miller_function_tate_2naf(P, Q, 0, r)
    f = final_exp_bw6_761(m, u0)
    return f

def optimal_ate_pairing_bw6_761(Q, P, u0):
    m = miller_loop_opt_ate_bw6_761(Q, P, u0)
    f = final_exp_bw6_761(m, u0)
    return f

####### membership testing, co-factor multiplication
# Formulas from Fuentes-Castaneda, Knapp, Rodriguez-Henriquez
# Faster hashing to G2.
# In: Miri, Vaudenay (eds.) SAC 2011. LNCS, vol. 7118, pp. 412-430.
# Springer, Heidelberg (Aug 2012). https://doi.org/10.1007/978-3-642-28496-0_25

def bw6_bls12_g1_mult_by_3r_trace_0_mod_r_u(P, omega, u0):
    """Return 3*r*P assuming phi(P) = (omega*x, y) = (-(t-1) mod r)*P

    Multiply P by 3*r in order to be able to check if 3*r*P == 0 in E(Fp)
    r = (u-1)^2/3*(u^2-u+1) + u
    assume that (trace mod r) = 0 mod u

    with  t   = -x^5 + 3*x^4 - 3*x^3 + x
    and  -t+1 = x^5 - 3*x^4 + 3*x^3 - x + 1
    and   3*r = x^6 - 2*x^5 + 2*x^3 + x + 1
    then -3*r = (-x-1)*(-t+1) + (x^3 - x^2 - x)
          3*r = (x+1)*(-t+1) - x^3 + x^2 + x
    and  3*r*P = (x+1)*phi(P) + (-x^3 + x^2 + x)*P

    alternatively if phi has eigenvalue (t-2) mod r:
    t-2 = -x^5 + 3*x^4 - 3*x^3 + x - 2
    3*r = (-x-1)*(t-2) - x^3 + x^2 - 1

    (u+1)*P + (u^3-u^2+1)*phi(P)  = r*(3*u^2-6*u+6)*P
    (-u^3+u^2+u)*P + (u+1)*phi(P) = 3*r*P

    bw6_761_g1_mult_by_r_alt(P) is this function with appropriate omega, u0.
    """
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = xP + P
    phi_Q = bw6_phi(Q, omega)
    R = -x3P + x2P + xP + phi_Q
    return R

def bw6_bls12_g1_mult_by_3r_trace_0_mod_r_u_alt(P, omega, u0):
    """Return 3*r*P assuming phi(P) = (omega*x, y) = ((t-2) mod r)*P

    Multiply P by 3*r in order to be able to check if 3*r*P == 0 in E(Fp)
    r = (u-1)^2/3*(u^2-u+1) + u
    assume that (trace mod r) = 0 mod u

    with  t   = -x^5 + 3*x^4 - 3*x^3 + x
    then  t-2 = -x^5 + 3*x^4 - 3*x^3 + x - 2
    and   3*r = x^6 - 2*x^5 + 2*x^3 + x + 1
    then  3*r = (-x-1)*(t-2) - x^3 + x^2 - 1
         -3*r = (x+1)*(t-2) + x^3 - x^2 + 1
    -> -3*r*P = (x+1)*phi(P) + (x^3 - x^2 + 1)*P
    """
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = xP + P
    phi_Q = bw6_phi(Q, omega)
    R = x3P - x2P + P + phi_Q
    return -R

def bw6_bls12_g1_check_membership_trace_u(P, omega, u0):
    R = bw6_bls12_g1_mult_by_3r_trace_0_mod_r_u(P, omega, u0)
    return R == P.curve()(0)

def bw6_bls12_g1_mult_by_3r_trace_3_mod_r_u(P, omega, u0):
    """Return 3*r*P assuming phi(P) = (omega*x, y) = (-(t-1) mod r)*P

    Multiply P by 3*r in order to be able to check if 3*r*P == 0 in E(Fp)
    r = (u-1)^2/3*(u^2-u+1) + u
    assume that (trace mod r) = 3 mod u
    with  t   = x^5 - 3*x^4 + 3*x^3 - x + 3
    and  -t+1 = -x^5 + 3*x^4 - 3*x^3 + x - 2
    and   3*r = x^6 - 2*x^5 + 2*x^3 + x + 1

    then -3*r = x^3 - x^2 + 1 + (x+1)*(-t+1)
          3*r = -x^3 + x^2 - 1 - (x+1)*(-t+1)
    and 3*r*P = (-x^3 + x^2 - 1)*P - (x+1)*phi(P)

    (u+1)*P + (-u^3+u^2+u)*phi(P) = r*(3*u^2-6*u+3)
    (u^3-u^2+1)*P + (u+1)*phi(P) = -3*r

    """
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = xP + P
    phi_Q = bw6_phi(Q, omega)
    R = -x3P + x2P - P - phi_Q
    return R

def bw6_bls12_g1_check_membership_trace_3u(P, omega, u0):
    R = bw6_bls12_g1_mult_by_3r_trace_3_mod_r_u(P, omega, u0)
    return R == P.curve()(0)

def bw6_bls12_g1_mult_by_3r_trace_3_mod_r_u_alt(P, omega, u0):
    """Return 3*r*P assuming phi(P) = (omega*x, y) = ((t-2) mod r)*P

    Multiply P by 3*r in order to be able to check if 3*r*P == 0 in E(Fp)
    r = (u-1)^2/3*(u^2-u+1) + u
    assume that (trace mod r) = 3 mod u
    with  t   = x^5 - 3*x^4 + 3*x^3 - x + 3
    and   t-2 = x^5 - 3*x^4 + 3*x^3 - x + 1
    and   3*r = x^6 - 2*x^5 + 2*x^3 + x + 1

    then  3*r = (x+1)*(t-2) - x^3 + x^2 + x
    and 3*r*P = (-x^3 + x^2 + x)*P + (x+1)*phi(P)

    (-u^3+u^2+u)*P + (u+1)*phi(P) = 3*r
    """
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = xP + P
    phi_Q = bw6_phi(Q, omega)
    R = -x3P + x2P + xP + phi_Q
    return R

def bw6_bls12_g1_mult_by_cofactor_trace_0_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    cx = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t-2) + (u^2 - 3*u + 3)*(u^2 - u + 1) - hy;

    Assume that eigenvalue of endomorphism phi(x,y) = (omega*x,y) is exactly -(t-1).
    Phi with (-omega-1) of eigenvalue (t-2) would give the function
    bw6_bls12_g1_mult_by_cofactor_trace_0_mod_r_u_alt below.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx
    l0 = (ht^2+3*hy^2)/4(u^3-u^2+1) - ht*(u-1)^2 - (ht-3*hy)/2
    l1 = (ht^2+3*hy^2)/4*(u+1) - (ht+3*hy)/2*(u-1)^2 - ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    L0 = d1*(u3P-u2P+P) - ht*(u2P -2*uP + P) - d3*P
    L1 = d1*(uP+P) - d2*(u2P -2*uP + P) - ht*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g1_mult_by_cofactor_trace_0_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    cx = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t-2) + (u^2 - 3*u + 3)*(u^2 - u + 1) - hy;

    Assume that eigenvalue of endomorphism phi(x,y) = (omega'*x,y) is (t-2).
    Phi with the other omega of eigenvalue -(t-1) gives the function
    bw6_bls12_g1_mult_by_cofactor_trace_0_mod_r_u above.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)/4*(u+1) - (ht+3*hy)/2*(u^2-2*u+1) - ht
    l1 = (ht^2+3*hy^2)/4*(u^3-u^2+1) - ht*(u^2-2*u+1) - (ht-3*hy)/2
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^3-u^2-u) - (ht-3*hy)/2*(u^2-2*u+2)  + ht
    l1 = -(ht^2+3*hy^2)/4*(u+1) + (ht+3*hy)/2*(u^2-2*u+1) + ht

    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    L0 = d1*(uP+P) - d2*(u2P -2*uP + P) - ht*P
    L1 = d1*(u3P-u2P+P) - ht*(u2P -2*uP + P) - d3*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g1_mult_by_cofactor_trace_3_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    cx = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t-2) + (u^2 - 3*u + 3)*(u^2 - u + 1) + hy;

    Assume that eigenvalue of endomorphism phi(x,y) = (omega*x,y) is exactly -(t-1).
    Phi with (-omega-1) of eigenvalue (t-2) would give the function
    bw6_bls12_g1_mult_by_cofactor_trace_3_mod_r_u_alt below.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx
    l0 = (ht^2+3*hy^2)/4(u^3-u^2+1) + ht + (ht+3*hy)/2*(u-1)^2
    l1 = (ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^2-2*u+2) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    L0 = d1*(u3P-u2P+P) + ht*P + d2*(u2P-2*uP+P)
    L1 = d1*(uP+P) - d3*(u2P -2*(uP - P)) + ht*P
    #or:
    #l0 = (ht^2+3*hy^2)/4*(-u-1) + (ht-3*hy)/2*(u^2-2*u+2) - ht
    #l1 = (ht^2+3*hy^2)/4*(u^3-u^2-u) + 4*ht*(u-1)^2 + (ht-3*hy)/2
    #L0 = -d1*(uP+P) + d3*(u2P -2*(uP-P)) - ht*P
    #L1 = d1*(u3P-u2P-uP) + ht*(u2P-2*uP+P) + d2*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g1_mult_by_cofactor_trace_3_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    cx = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t-2) + (u^2 - 3*u + 3)*(u^2 - u + 1) + hy

    Assume that eigenvalue of endomorphism phi(x,y) = (omega'*x,y) is (t-2).
    Phi with the other omega of eigenvalue -(t-1) gives the function
    bw6_bls12_g1_mult_by_cofactor_trace_3_mod_r_u above.

    l0 + l1*eigenvalue_phi_mod_cx = 0 mod cx

    l0 = (ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^2-2*u+2) + ht
    l1 = (ht^2+3*hy^2)/4*(u^3-u^2+1) + (ht+3*hy)/2*(u^2-2*u+1) + ht
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^3-u^2-u) + ht*(u^2-2*u+1) + (ht-3*hy)/2
    l1 = -(ht^2+3*hy^2)/4*(u+1) + (ht-3*hy)/2*(u^2-2*u+2) - ht

    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    htP = ht*P
    L0 = d1*(uP+P) - d3*(u2P -2*(uP - P)) + htP
    L1 = d1*(u3P-u2P+P) + d2*(u2P-2*uP+P) + htP
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g2_mult_by_cofactor_trace_0_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    cx = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t-1) + (u^4-4*u^3+7*u^2-6*u+3) + hy
    cx = (ht^2+3*hy^2)/4*r + (ht-hy)/2*t + (u^4-4*u^3+7*u^2-6*u+3) + hy + (-ht+hy)/2
    cx = (ht^2+3*hy^2)/4*r + (ht-hy)/2*t + (u^4-4*u^3+7*u^2-6*u+3) + (-ht+3hy)/2

    The eigenvalue mod r of endomorphism phi(x,y) = (omega*x,y)
    is either -(t0-1) or (t0-2) where t0 is t mod r.
    Note that the trace of E mod r and the trace of E' mod r are the same

    l0 + l1*eigenvalue_phi_mod_c2x = 0 mod c2x

    l0 = (ht^2+3*hy^2)/4*(u+1) + (ht+3*hy)/2*(u^2-2*u+2) -ht
    l1 = (ht^2+3*hy^2)/4*(u^3-u^2+1) -(ht-3*hy)/2*(u^2-2*u+1) -ht
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^3-u^2-u) - ht*(u^2-2*u+1) - (ht+3*hy)/2
    l1 = -(ht^2+3*hy^2)/4*(u+1) -(ht+3*hy)/2*(u^2-2*u+2) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    htP = ht*P
    L0 = d1 * (uP + P) + d2 * (u2P-2*(uP-P)) - htP
    L1 = d1 * (u3P-u2P+P) - d3*(u2P-2*uP+P) - htP
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g2_mult_by_cofactor_trace_0_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    c2x = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t-1) + (u^4-4*u^3+7*u^2-6*u+3) + hy

    The eigenvalue mod r of endomorphism phi(x,y) = (omega*x,y)
    is either -(t0-1) or (t0-2) where t0 is t mod r.
    This function is with the other choice of eigenvalue compared to
    bw6_bls12_g2_mult_by_cofactor_trace_0_mod_r_u
    Note that the trace of E mod r and the trace of E' mod r are the same

    l0 + l1*eigenvalue_phi_mod_c2x = 0 mod c2x

    l0 = -(ht^2+3*hy^2)/4*(u+1) -(ht+3*hy)/2*(u^2-2*u+2) + ht
    l1 = (ht^2+3*hy^2)/4*(u^3-u^2-u) - ht*(u^2-2*u+1) - (ht+3*hy)/2
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^3-u^2+1) - (ht-3*hy)/2*(u^2-2*u+1) - ht
    l1 = (ht^2+3*hy^2)/4*(u+1) + (ht+3*hy)/2*(u^2-2*u+2) - ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    #d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    L0 = -d1 * (uP + P) - d2 * (u2P-2*(uP-P)) + ht*P
    L1 = d1 * (u3P-u2P-uP) - ht*(u2P-2*uP+P) - d2*P
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g2_mult_by_cofactor_trace_3_mod_r_u(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    c2x = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t-1) + (u^4-4*u^3+7*u^2-6*u+3) - hy

    The eigenvalue mod r of endomorphism phi(x,y) = (omega*x,y)
    is either -(t3-1) or (t3-2) where t3 is t mod r.
    Note that the trace of E mod r and the trace of E' mod r are the same

    l0 + l1*eigenvalue_phi_mod_c2x = 0 mod c2x

    l0 = -(ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^2-2*u+1) - ht
    l1 = (ht^2+3*hy^2)/4*(u^3-u^2-u) + (ht+3*hy)/2*(u^2-2*u+2) - ht
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^3-u^2+1) + ht*(u^2-2*u+1) + (ht+3*hy)/2
    l1 = (ht^2+3*hy^2)/4*(u+1) + (ht-3*hy)/2*(u^2-2*u+1) + ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    htP = ht*P
    L0 = -d1*(uP+P) - d3*(u2P-2*uP+P) - htP
    L1 = d1*(u3P-u2P-uP) + d2*(u2P-2*(uP-P)) - htP
    return L0 + bw6_phi(L1, omega)

def bw6_bls12_g2_mult_by_cofactor_trace_3_mod_r_u_alt(P, omega, u0, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fq)
    c2x = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t-1) + (u^4-4*u^3+7*u^2-6*u+3) - hy

    The eigenvalue mod r of endomorphism phi(x,y) = (omega*x,y)
    is either -(t3-1) or (t3-2) where t3 is t mod r.
    This function is with the other choice of eigenvalue compared to
    bw6_bls12_g2_mult_by_cofactor_trace_3_mod_r_u
    Note that the trace of E mod r and the trace of E' mod r are the same

    l0 + l1*eigenvalue_phi_mod_c2x = 0 mod c2x

    l0 = (ht^2+3*hy^2)/4*(u+1) + (ht-3*hy)/2*(u^2-2*u+1) + ht
    l1 = (ht^2+3*hy^2)/4*(u^3-u^2+1) + ht*(u^2-2*u+1) + (ht+3*hy)/2
    equivalently:
    l0 = (ht^2+3*hy^2)/4*(u^3-u^2-u) + (ht+3*hy)/2*(u^2-2*u+2) - ht
    l1 = -(ht^2+3*hy^2)/4*(u+1) - (ht-3*hy)/2*(u^2-2*u+1) - ht
    """
    d1 = (ht**2+3*hy**2)//4
    d2 = (ht+3*hy)//2
    d3 = (ht-3*hy)//2
    uP = u0*P
    u2P = u0*uP
    u3P = u0*u2P
    L0 = d1*(uP+P) + d3*(u2P-2*uP+P) + ht*P
    L1 = d1*(u3P-u2P+P) + ht*(u2P-2*uP+P) + d2*P
    return L0 + bw6_phi(L1, omega)

def final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_subfunction(B, u, ht, hy):
    """
    return B^((ht^2+3*hy^2)/4*r + 3*((ht+hy)/2*t0/3 + (d-1)/3) + 1)

    INPUT:
    - `B`: a finite field element, with B.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t0
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t0^2-4*q = -D*y^2

    RETURN: B^e

    parameters are
    r       = (u-1)^2/3*(u^4-u^2+1) + u
    r      == (u^6 - 2*u^5 + 2*u^3 + u + 1)/3
    d       = u^4 - 4*u^3 + 7*u^2 - 6*u + 3
    (d-1)/3 = (u - 1)/3 * (u-1) * ((u-1)^2 + 1)
    t3      = u^5 - 3*u^4 + 3*u^3 - u + 3
    t3/3   == (u-1)^2/3*(u^3-u^2+1) + (u-1)/3 + 1

    (u^3-u^2+1) == (u+1)*(u^2-2*u+2) - 1

    a sequence to obtain the exponents r, t3/3 and (d-1)/3 is
    C = (u-1)/3
    D = (u-1)*C
    D == (u-1)^2/3
    E = D * ((u-1)^2 + 1)
    E == (d-1)/3
    D = -D
    F = D+1
    G = (u+1) * E + F
    G == (u-1)^2/3*(u^3-u^2+1) + 1
    H = G + C
    H == t3/3
    I = (u+1)*(G+D) - F
    I == r

    cost e((u-1)/3) + 3e(u-1) + 2e(u+1) + e(d1) + e(d2) + 10 M + S + 2cj
    """
    # cost e((u-1)/3) + 3e(u-1) + 2e(u+1) + e(d2) + e(d1) + 10M + S + 2 cj
    C = B**((u-1)//3)
    D = C**(u-1)
    E = (D**(u-1))**(u-1) * D             # B^((d-1)/3)
    D = D.conjugate()
    Fc = D * B
    G = E**(u+1) * Fc
    H = G * C                             # B^(t3/3)
    I = (G * D)**(u+1) * Fc.conjugate()   # B^r
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+hy)//2
    J = H**d1 * E
    K = J**2 * J * B * I**d2
    return K

def final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_alt(m, u, ht, hy):
    """
    Exponentiation to (u^3-u^2+1)*Phi_k(p(u))/r(u)
    e0 = (ht^2+3*hy^2)/4*r*(-u-1) + (ht+hy)/2*tr3*(-u-1) + cc3*(-u-1) - 3
    e1 = (ht^2+3*hy^2)/4*r*(u^3-u^2+1) + (ht+hy)/2*tr3*(u^3-u^2+1) + cc3*(u^3-u^2+1) + 3*(u^2-2*u+2)

    c_bw3 = (ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3 - ht

    (e0 + e1*q)
    = ((-u-1) + (u^3-u^2+1)*q)*((ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3) + 3*(-1 + (u^2-2*u+2)*q)
    = ((-u-1) + (u^3-u^2+1)*q)*(c_bw3 + ht) + 3*(-1 + (u^2-2*u+2)*q)
    = ((u+1)*(-1 + (u^2-2*u+2)*q) - q)*(c_bw3 + ht) + 3*(-1 + (u^2-2*u+2)*q)

    = ((u+1)*((u^2-2*u+2)*q -1) - q)*(c_bw3 + ht) + 3*((u^2-2*u+2)*q -1)
    = ((u+1)*(((u-1)^2+1)*q -1) - q)*(c_bw3 + ht) + 3*(((u-1)^2+1)*q -1)

    cc3 = u^4 - 4*u^3 + 7*u^2 - 6*u + 3
    tr3 = u^5 - 3*u^4 + 3*u^3 - u + 3
    r   = (u^6 - 2*u^5 + 2*u^3 + u + 1)/3

    see final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_subfunction

    cost subfunction e((u-1)/3) + 3e(u-1) + 2e(u+1) + e(d2) + e(d1) + 10M +  S       + 2cj
    cost here                     2e(u-1) +  e(u+1)                 +  5M +  S + frb + 2cj
    cost total       e((u-1)/3) + 5e(u-1) + 3e(u+1) + e(d1) + e(d2) + 15M + 2S + frb + 4cj
    """
    # cost 2e(u-1) + e(u+1) + 5 M + S + frb + 2cj
    mq = m.frobenius().conjugate()
    A = mq**(u-1)
    A = A**(u-1)
    A = (m * mq * A).conjugate()      # A = m**(((u-1)^2+1)q - 1)
    B = A**(u+1) * mq                 # B = m**((u+1) + (u^3-u^2-u)*q)
    A = A**2 * A                      # A = m**(3*(-1 + (u^2-2*u+2)*q))
    C = final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_subfunction(B, u, ht, hy)
    return A * C

def final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u(m, u, ht, hy):
    """
    Exponentiation to (u+1)*Phi_k(p(u))/r(u)

    e0 = (ht^2+3*hy^2)/4*r*(u^3-u^2-u) + (ht+hy)/2*tr3*(u^3-u^2-u) + cc3*(u^3-u^2-u) + 3*(u^2-2*u+1)
    e1 = (ht^2+3*hy^2)/4*r*(u+1) + (ht+hy)/2*tr3*(u+1) + cc3*(u+1) + 3

    ((u^3-u^2-u) + (u+1)*q)*((ht^2+3*hy^2)/4*r + (ht+hy)/2*tr3 + cc3) + 3*((u^2-2*u+1) + q)
    with (u^3-u^2-u) + (u+1)*q == (u+1)(u^2-2*u+1 + q) - 1
    ((u+1)(u^2-2*u+1 + q) - 1)*((ht^2+3*hy^2)/4*r + 3*((ht+hy)/2*tr3/3 + (cc3-1)/3) + 1) + 3*(u^2-2*u+1 + q)

    cost subfunction e((u-1)/3) + 3e(u-1) + 2e(u+1) + e(d2) + e(d1) + 10M +  S       + 2cj
    cost here                     2e(u-1) +  e(u+1)                 +  4M +  S + frb +  cj
    cost total       e((u-1)/3) + 5e(u-1) + 3e(u+1) + e(d1) + e(d2) + 14M + 2S + frb + 3cj
    """
    # cost 2e(u-1) + e(u+1) + 4M + S + frb + cj
    A = m**(u-1)
    A = A**(u-1)
    A = A * m.frobenius()        # A = m^((u-1)^2 + q)
    B = A**(u+1) * m.conjugate() # B = m**((u^3-u^2-u) + (u+1)*q)
    A = A**2 * A                 # A = m**(3*((u^2-2*u+1) + q))
    C = final_exp_hard_bw6_bls12_trace_3_mod_r_mod_u_subfunction(B, u, ht, hy)
    return A * C

def final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_subfunction(B, u, ht, hy):
    """
    return B^((ht^2+3*hy^2)/4*r + 3*((ht+hy)/2*t0/3 + (d-1)/3) + 1)

    INPUT:
    - `B`: a finite field element, with B.conjugate() available (absolute extension of a prime field)
    - `u`: relative integer, the curve seed
    - `ht`: relative integer, curve parameter, lifting cofactor for the trace t0
    - `hy`: relative integer, curve parameter, lifting cofactor for the square part y in t0^2-4*q = -D*y^2

    RETURN: B^e

    parameters are
    cc0 = u^4 - 4*u^3 + 7*u^2 - 6*u + 3
    tr0 = -u^5 + 3*u^4 - 3*u^3 + u
    r = (u^6 - 2*u^5 + 2*u^3 + u + 1)/3

    a sequence to obtain the exponents r, t0/3 and (d-1)/3 is
    a = (u-1)/3
    b = (u-1)*a
    b == (u-1)^2/3
    c = b * ((u-1)^2 + 1)
    c == (cc0-1)/3
    e = -(u+1) * c + b - a
    e == tr0/3
    f = -(u+1)*(e + b) + a + 1
    f == r

    cost
    e((u-1)/3) + 3e(u-1) + 2e(u+1) + e(d1) + e(d2) + 10 M + S + 2cj
    """
    # cost e((u-1)/3) + 3e(u-1) +2e(u+1) + e(d1) + e(d2) + 10 M + S + 2cj
    C = B**((u-1)//3)
    D = C**(u-1)
    E = (D**(u-1))**(u-1) * D                # B^((d-1)/3)
    F = (E**(u+1) * C).conjugate() * D       # B^(tr0/3)
    G = ((F * D)**(u+1)).conjugate() * C * B # B^r
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-hy)//2
    H = F**d1 * E
    H = H**2 * H * B * G**d2
    return H

def final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u(m, u, ht, hy):
    """
    Exponentiation to (u+1)*Phi_k(p(u))/r(u)
    e0 = -((ht^2+3*hy^2)/4*r*(u^3-u^2+1) + (ht-hy)/2*tr0*(u^3-u^2+1) + cc0*(u^3-u^2+1)) + 3*(u^2-2*u+2)
    e1 = (ht^2+3*hy^2)/4*r*(u+1) + (ht-hy)/2*tr0*(u+1) + cc0*(u+1) - 3

    (-(u^3-u^2+1) + (u+1)*q)*((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0) - 3*(-(u^2-2*u+2) + q)
    with (-(u^3-u^2+1) + (u+1)*q) == (-(u^2-2*u+2) + q)*(u+1) + 1
                                  == (-((u-1)^2+1) + q)*(u+1) + 1

    ((q - (u-1)^2-1)*(u+1) + 1)*((ht^2+3*hy^2)/4*r + 3*((ht-hy)/2*tr0/3 + (cc0-1)/3) + 1) - 3*(q - (u-1)^2-1)
    see final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_subfunction

    cost subfunction e((u-1)/3) + 3e(u-1) +2e(u+1) + e(d1) + e(d2) + 10M +  S       + 2cj
    cost here                     2e(u-1) + e(u+1)                 +  5M +  S + frb + 2cj
    cost total       e((u-1)/3) + 5e(u-1) +3e(u+1) + e(d1) + e(d2) + 15M + 2S + frb + 4cj
    """
    # cost: 2*e(u-1) + e(u+1) + 5M + S + frb + 2 cj
    A = m**(u-1)
    A = A**(u-1)
    A = (m * A).conjugate() * m.frobenius()          # A = m^(q - (u-1)^2 - 1)
    B = A**(u+1) * m                                 # B = m**(-(u^3-u^2-u) + (u+1)*q)
    A = A**2 * A
    A = A.conjugate()                                # A = m**(-3*(q - (u-1)^2 - 1))
    C = final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_subfunction(B, u, ht, hy)
    return A * C

def final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_alt(m, u, ht, hy):
    """
    Exponentiation to (u^3-u^2-u)*Phi_k(p(u))/r(u)

    e0 = ((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0)*(u+1) - 3
    e1 = ((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0)*(u^3-u^2-u) - 3*(u^2-2*u+1)

    ((u+1) + (u^3-u^2-u)*q)*((ht^2+3*hy^2)/4*r + (ht-hy)/2*tr0 + cc0) - 3*(1 + (u^2-2*u+1)*q)
    (u^3-u^2-u) == (u+1)*(u^2-2*u+1) - 1
    (u^3-u^2-u) == (u+1)*((u-1)^2) - 1
    ((u+1)(1 + (u-1)^2*q) - q)*((ht^2+3*hy^2)/4*r + 3*((ht-hy)/2*tr0/3 + (cc0-1)/3) + 1) - 3*(1 + (u-1)^2*q)
    see final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_subfunction

    cost subfunction e((u-1)/3) + 3e(u-1) +2e(u+1) + e(d1) + e(d2) + 10M +  S       + 2cj
    cost here                     2e(u-1) + e(u+1)                 +  4M +  S + frb + 2cj
    cost total       e((u-1)/3) + 5e(u-1) +3e(u+1) + e(d1) + e(d2) + 14M + 2S + frb + 4cj
    """
    # cost: 2*e(u-1) + e(u+1) + 4M + S + frb + 2 cj
    mq = m.frobenius()
    A = mq**(u-1)
    A = A**(u-1)
    A = m * A                                # A = m^((u-1)^2*q + 1)
    B = A**(u+1) * mq.conjugate()            # B = m**((u+1) + (u^3-u^2-u)*q)
    A = A**2 * A
    A = A.conjugate()                        # A = m**(-3*((u-1)^2*q + 1))
    C = final_exp_hard_bw6_bls12_trace_0_mod_r_mod_u_subfunction(B, u, ht, hy)
    return A * C

####### BW6_761_BLS12

def bw6_761_phi(P):
    # omegx := (103*x^11 - 482*x^10 + 732*x^9 + 62*x^8 - 1249*x^7 + 1041*x^6 + 214*x^5 - 761*x^4 + 576*x^3 + 11*x^2 - 265*x + 66)/21;
    Fp = P[0].parent()
    omega = Fp(0x531DC16C6ECD27AA846C61024E4CCA6C1F31E53BD9603C2D17BE416C5E4426EE4A737F73B6F952AB5E57926FA701848E0A235A0A398300C65759FC45183151F2F082D4DCB5E37CB6290012D96F8819C547BA8A4000002F962140000000002A)
    return P.curve()([omega*P[0], P[1]])

def bw6_761_g1_mult_by_cofactor(P):
    # multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    # cx = (103*x^6 - 173*x^5 - 96*x^4 + 293*x^3 + 21*x^2 + 52*x + 172)/3
    u0 = 0x8508C00000000001
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = (7*x2P +89*xP + 130*P)
    phi_Q = bw6_761_phi(Q)
    R = (103*x3P -83*x2P -40*xP + 136*P) + phi_Q
    return R

def bw6_761_g1_mult_by_cofactor_alt(P):
    # multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    # cx = (103*x^6 - 173*x^5 - 96*x^4 + 293*x^3 + 21*x^2 + 52*x + 172)/3
    # alternative formula
    u0 = 0x8508C00000000001
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = (103*x3P - 90*x2P - 129*xP + 6*P)
    phi_Q = bw6_761_phi(Q)
    R = -(7*x2P + 89*xP + 130*P) + phi_Q
    return R

def bw6_761_g1_mult_by_r(P):
    """ Multiply P by a multiple of r in order to be able to check if r*P == 0 in E(Fp)
    returns (u+1)*P + (u^3-u^2+1)Phi(P)
    """
    u0 = 0x8508C00000000001
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = (x3P - x2P + P)
    phi_Q = bw6_761_phi(Q)
    R = xP + P + phi_Q
    return R

def bw6_761_g1_mult_by_r_alt(P):
    # multiply P by a multiple of r in order to be able to check if r*P == 0 in E(Fp)
    # alternative formula
    u0 = 0x8508C00000000001
    xP = u0*P
    x2P = u0*xP
    x3P = u0*x2P
    Q = -(xP + P)
    phi_Q = bw6_761_phi(Q)
    R = (x3P - x2P - xP) + phi_Q
    return R

def bw6_761_g1_check_membership(P):
    R = bw6_761_g1_mult_by_r(P)
    return R == P.curve()(0)

def bw6_761_g1_check_membership_alt(P):
    # alternative formula
    R = bw6_761_g1_mult_by_r_alt(P)
    return R == P.curve()(0)

def bw6_761_g2_mult_by_cofactor(Q):
    # multiply point Q by a multiple of the cofactor c' where c'*r = #E'(Fp)
    # cx = (103*x^6 - 173*x^5 - 96*x^4 + 293*x^3 + 21*x^2 + 52*x + 151)/3
    # note that only the constant coefficient of cx differs from cx for G1
    u0 = 0x8508C00000000001
    xQ = u0*Q
    x2Q = u0*xQ
    x3Q = u0*x2Q
    S = (7*x2Q -117*xQ - 109*Q)
    phi_S = bw6_761_phi(S)
    R = (103*x3Q -83*x2Q -143*xQ + 27*Q) + phi_S
    return R

def bw6_761_g2_mult_by_cofactor_alt(Q):
    # multiply point Q by a multiple of the cofactor c' where c'*r = #E'(Fp)
    # cx = (103*x^6 - 173*x^5 - 96*x^4 + 293*x^3 + 21*x^2 + 52*x + 151)/3
    # note that only the constant coefficient of cx differs from cx for G1
    # alternative formula
    u0 = 0x8508C00000000001
    xQ = u0*Q
    x2Q = u0*xQ
    x3Q = u0*x2Q
    S = (103*x3Q - 90*x2Q - 26*xQ + 136*Q)
    phi_S = bw6_761_phi(S)
    R = (-7*x2Q +117*xQ + 109*Q) + phi_S
    return R

def bw6_761_g2_mult_by_r(Q):
    # multiply Q by a multiple of r in order to be able to check if r*Q == 0 in E'(Fp)
    u0 = 0x8508C00000000001
    xQ = u0*Q
    x2Q = u0*xQ
    x3Q = u0*x2Q
    S = (x3Q -x2Q - xQ)
    phi_S = bw6_761_phi(S)
    R = (-xQ - Q) + phi_S
    return R

def bw6_761_g2_mult_by_r_alt(Q):
    # multiply Q by a multiple of r in order to be able to check if r*Q == 0 in E'(Fp)
    # alternative formula
    u0 = 0x8508C00000000001
    xQ = u0*Q
    x2Q = u0*xQ
    x3Q = u0*x2Q
    S = (xQ + Q)
    phi_S = bw6_761_phi(S)
    R = (x3Q - x2Q + Q) + phi_S
    return R

def bw6_761_g2_check_membership(Q):
    # check that r*Q == 0 in E'(Fp)
    R = bw6_761_g2_mult_by_r(Q)
    return R == Q.curve()(0)

def bw6_761_g2_check_membership_alt(Q):
    # check that r*Q == 0 in E'(Fp)
    # alternative formula
    R = bw6_761_g2_mult_by_r_alt(Q)
    return R == Q.curve()(0)

