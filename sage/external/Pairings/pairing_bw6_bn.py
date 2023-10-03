from pairing import *

####### optimal ate pairing
def miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u(Q, P, u):
    """
    Return f_{2u,Q}(P)*frobenius(f_{2u(3u+1)+1,Q}(P))

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 0
    trace = -18*x^3 - 18*x^2 - 9*x + r(x) * ht
    v = 3*u+1
    as f_{2u,Q}(P) * frobenius((f_2u)_{v,[u]Q}(P))
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    """
    if u < 0:
        m_2u, u2Q = miller_function_ate((Q[0], -Q[1]), P, 0, -2*u)
    else:
        m_2u, u2Q = miller_function_ate(Q, P, 0, 2*u)
    Z1 = 1/u2Q[2]
    Z2 = Z1**2
    u2Q = (u2Q[0]*Z2, u2Q[1]*Z1*Z2, 1, 1)
    v = 3*u+1
    if v < 0:
        m_2u_inv = m_2u.conjugate()
        u2Q = (u2Q[0], -u2Q[1], 1, 1)
        m_2uv, u2vQ = miller_function_ate(u2Q, P, 0, -v, m0=m_2u_inv)
    else:
        m_2uv, u2vQ = miller_function_ate(u2Q, P, 0, v, m0=m_2u)
    l, u2v_1Q = add_line_j(u2vQ, (Q[0], Q[1]), (P[0], P[1]))
    m_2uv_1 = m_2uv * l
    return m_2u * m_2uv_1.frobenius()

def miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_2naf(Q, P, u):
    """
    Return f_{2u,Q}(P)*frobenius(f_{2u(3u+1)+1,Q}(P))

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 0
    trace = -18*x^3 - 18*x^2 - 9*x + r(x) * ht
    v = 3*u+1
    as f_{2u,Q}(P) * frobenius((f_2u)_{v,[u]Q}(P))
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    u and v in 2-naf form
    """
    if u < 0:
        m_2u, u2Q = miller_function_ate_2naf((Q[0], -Q[1]), P, 0, -2*u)
    else:
        m_2u, u2Q = miller_function_ate_2naf(Q, P, 0, 2*u)
    Z1 = 1/u2Q[2]
    Z2 = Z1**2
    u2Q = (u2Q[0]*Z2, u2Q[1]*Z1*Z2, 1, 1)
    v = 3*u+1
    m_2u_inv = m_2u.conjugate()
    if v < 0:
        u2Q = (u2Q[0], -u2Q[1], 1, 1)
        m_2uv, u2vQ = miller_function_ate_2naf(u2Q, P, 0, -v, m0=m_2u_inv, m0_inv=m_2u)
    else:
        m_2uv, u2vQ = miller_function_ate_2naf(u2Q, P, 0, v, m0=m_2u, m0_inv=m_2u_inv)
    l, u2v_1Q = add_line_j(u2vQ, (Q[0], Q[1]), (P[0], P[1]))
    m_2uv_1 = m_2uv * l
    return m_2u * m_2uv_1.frobenius()

def miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u(Q, P, u):
    """
    Return frobenius(f_{2u,Q}(P))*f_{2u(3u+1)+1,Q}(P)

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 3
    trace = 18*x^3 + 18*x^2 + 9*x + 3 + r(x) * ht
    v = 3*u+1
    as frobenius(f_{2u,Q}(P)) * (f_2u)_{v,[u]Q}(P)
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    """
    if u < 0:
        m_2u, u2Q = miller_function_ate((Q[0], -Q[1]), P, 0, -2*u)
    else:
        m_2u, u2Q = miller_function_ate(Q, P, 0, 2*u)
    Z1 = 1/u2Q[2]
    Z2 = Z1**2
    u2Q = (u2Q[0]*Z2, u2Q[1]*Z1*Z2, 1, 1)
    v = 3*u+1
    if v < 0:
        m_2u_inv = m_2u.conjugate()
        u2Q = (u2Q[0], -u2Q[1], 1, 1)
        m_2uv, u2vQ = miller_function_ate(u2Q, P, 0, -v, m0=m_2u_inv)
    else:
        m_2uv, u2vQ = miller_function_ate(u2Q, P, 0, v, m0=m_2u)
    l, u2v_1Q = add_line_j(u2vQ, (Q[0], Q[1]), (P[0], P[1]))
    m_2uv_1 = m_2uv * l
    return m_2u.frobenius() * m_2uv_1

def miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_2naf(Q, P, u):
    """
    Return frobenius(f_{2u,Q}(P))*f_{2u(3u+1)+1,Q}(P)

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 3
    trace = 18*x^3 + 18*x^2 + 9*x + 3 + r(x) * ht
    v = 3*u+1
    as frobenius(f_{2u,Q}(P)) * (f_2u)_{v,[u]Q}(P)
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    u, v in 2-naf
    """
    if u < 0:
        m_2u, u2Q = miller_function_ate_2naf((Q[0], -Q[1]), P, 0, -2*u)
    else:
        m_2u, u2Q = miller_function_ate_2naf(Q, P, 0, 2*u)
    Z1 = 1/u2Q[2]
    Z2 = Z1**2
    u2Q = (u2Q[0]*Z2, u2Q[1]*Z1*Z2, 1, 1)
    v = 3*u+1
    m_2u_inv = m_2u.conjugate()
    if v < 0:
        u2Q = (u2Q[0], -u2Q[1], 1, 1)
        m_2uv, u2vQ = miller_function_ate_2naf(u2Q, P, 0, -v, m0=m_2u_inv, m0_inv=m_2u)
    else:
        m_2uv, u2vQ = miller_function_ate_2naf(u2Q, P, 0, v, m0=m_2u, m0_inv=m_2u_inv)
    l, u2v_1Q = add_line_j(u2vQ, (Q[0], Q[1]), (P[0], P[1]))
    m_2uv_1 = m_2uv * l
    return m_2u.frobenius() * m_2uv_1

def miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_aklgl(Q,P,b_t,u,Fp6,D_twist=False,xi=None):
    """
    Return f_{2u,Q}(P)*frobenius(f_{2u(3u+1)+1,Q}(P))

    INPUT:
    -`Q`: point on the twist of degree 6 defined over Fp
    -`P`: point on E(Fp)
    -`b_t`: coefficient of the twist equation in y^2 = x^3 + b_t
    -`u`: seed
    -`Fp6`: degree 6 extension over a prime field
    -`D_twist`: whether the twisted curve is a D-twist of an M-twist
    -`xi`: the residue in Fp so that Fp6 = Fp[x]/(x^6-xi)

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 0
    trace = -18*x^3 - 18*x^2 - 9*x + r(x) * ht
    v = 3*u+1
    as f_{2u,Q}(P) * frobenius((f_2u)_{v,[u]Q}(P))
    with AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    """
    PP = (P[0], P[1])
    if xi is None:
        xi = -Fp6.modulus().constant_coefficient() # works with absolute and towering of extension
    if u < 0:
        m_2u, u2Q = miller_function_ate_aklgl((Q[0], -Q[1], 1), PP, b_t, -2*u, Fp6, D_twist=D_twist)
    else:
        m_2u, u2Q = miller_function_ate_aklgl((Q[0], Q[1], 1), PP, b_t, 2*u, Fp6, D_twist=D_twist)
    Z1 = 1/u2Q[2]
    u2Q = (u2Q[0]*Z1, u2Q[1]*Z1, 1)
    v = 3*u+1
    if v < 0:
        m_2u_inv = m_2u.conjugate()
        u2Q = (u2Q[0], -u2Q[1], 1)
        m_2uv, u2vQ = miller_function_ate_aklgl(u2Q, PP, b_t, -v, Fp6, m0=m_2u_inv, D_twist=D_twist)
    else:
        m_2uv, u2vQ = miller_function_ate_aklgl(u2Q, PP, b_t, v, Fp6, m0=m_2u, D_twist=D_twist)
    l, u2v_1Q = add_line_h_a0_twist6_aklgl(u2vQ, (Q[0], Q[1]), PP, D_twist=D_twist)
    if D_twist:
        m_2uv_1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_2uv, xi, Fp6)
    else:
        m_2uv_1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_2uv, xi, Fp6)
    return m_2u * m_2uv_1.frobenius()

def miller_loop_opt_ate_bw6_bn_trace_0_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fp6, D_twist=False, xi=None):
    """
    Return f_{2u,Q}(P)*frobenius(f_{2u(3u+1)+1,Q}(P))

    INPUT:
    -`Q`: point on the twist of degree 6 defined over Fp
    -`P`: point on E(Fp)
    -`b_t`: coefficient of the twist equation in y^2 = x^3 + b_t
    -`u`: seed
    -`Fp6`: degree 6 extension over a prime field
    -`D_twist`: whether the twisted curve is a D-twist of an M-twist
    -`xi`: the residue in Fp so that Fp6 = Fp[x]/(x^6-xi)

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 0
    trace = -18*x^3 - 18*x^2 - 9*x + r(x) * ht
    v = 3*u+1
    as f_{2u,Q}(P) * frobenius((f_2u)_{v,[u]Q}(P))
    with AKLGL formulas
    u, v in 2-naf
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    """
    PP = (P[0], P[1])
    if xi is None:
        xi = -Fp6.modulus().constant_coefficient() # works with absolute and towering of extension
    if u < 0:
        m_2u, u2Q = miller_function_ate_2naf_aklgl((Q[0], -Q[1], 1), PP, b_t, -2*u, Fp6, D_twist=D_twist)
    else:
        m_2u, u2Q = miller_function_ate_2naf_aklgl((Q[0], Q[1], 1), PP, b_t, 2*u, Fp6, D_twist=D_twist)
    Z1 = 1/u2Q[2]
    u2Q = (u2Q[0]*Z1, u2Q[1]*Z1, 1)
    v = 3*u+1
    m_2u_inv = m_2u.conjugate()
    if v < 0:
        u2Q = (u2Q[0], -u2Q[1], 1)
        m_2uv, u2vQ = miller_function_ate_2naf_aklgl(u2Q, PP, b_t, -v, Fp6, m0=m_2u_inv, m0_inv=m_2u, D_twist=D_twist)
    else:
        m_2uv, u2vQ = miller_function_ate_2naf_aklgl(u2Q, PP, b_t, v, Fp6, m0=m_2u, m0_inv=m_2u_inv, D_twist=D_twist)
    l, u2v_1Q = add_line_h_a0_twist6_aklgl(u2vQ, (Q[0], Q[1]), PP, D_twist=D_twist)
    if D_twist:
        m_2uv_1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_2uv, xi, Fp6)
    else:
        m_2uv_1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_2uv, xi, Fp6)
    return m_2u * m_2uv_1.frobenius()

def miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_aklgl(Q,P,b_t,u,Fp6,D_twist=False,xi=None):
    """
    Return frobenius(f_{2u,Q}(P))*f_{2u(3u+1)+1,Q}(P)

    INPUT:
    -`Q`: point on the twist of degree 6 defined over Fp
    -`P`: point on E(Fp)
    -`b_t`: coefficient of the twist equation in y^2 = x^3 + b_t
    -`u`: seed
    -`Fp6`: degree 6 extension over a prime field
    -`D_twist`: whether the twisted curve is a D-twist of an M-twist
    -`xi`: the residue in Fp so that Fp6 = Fp[x]/(x^6-xi)

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 3
    trace = 18*x^3 + 18*x^2 + 9*x + 3 + r(x) * ht
    v = 3*u+1
    as frobenius(f_{2u,Q}(P)) * (f_2u)_{v,[u]Q}(P)
    with AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    """
    PP = (P[0], P[1])
    if xi is None:
        xi = -Fp6.modulus().constant_coefficient() # works with absolute and towering of extension
    if u < 0:
        m_2u, u2Q = miller_function_ate_aklgl((Q[0], -Q[1], 1), PP, b_t, -2*u, Fp6, D_twist=D_twist)
    else:
        m_2u, u2Q = miller_function_ate_aklgl((Q[0], Q[1], 1), PP, b_t, 2*u, Fp6, D_twist=D_twist)
    Z1 = 1/u2Q[2]
    u2Q = (u2Q[0]*Z1, u2Q[1]*Z1, 1)
    v = 3*u+1
    if v < 0:
        m_2u_inv = m_2u.conjugate()
        u2Q = (u2Q[0], -u2Q[1], 1)
        m_2uv, u2vQ = miller_function_ate_aklgl(u2Q, PP, b_t, -v, Fp6, m0=m_2u_inv, D_twist=D_twist)
    else:
        m_2uv, u2vQ = miller_function_ate_aklgl(u2Q, PP, b_t, v, Fp6, m0=m_2u, D_twist=D_twist)
    l, u2v_1Q = add_line_h_a0_twist6_aklgl(u2vQ, (Q[0], Q[1]), PP, D_twist=D_twist)
    if D_twist:
        m_2uv_1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_2uv, xi, Fp6)
    else:
        m_2uv_1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_2uv, xi, Fp6)
    return m_2u.frobenius() * m_2uv_1

def miller_loop_opt_ate_bw6_bn_trace_3_mod_r_mod_u_aklgl_2naf(Q, P, b_t, u, Fp6, D_twist=False, xi=None):
    """
    Return frobenius(f_{2u,Q}(P))*f_{2u(3u+1)+1,Q}(P)

    INPUT:
    -`Q`: point on the twist of degree 6 defined over Fp
    -`P`: point on E(Fp)
    -`b_t`: coefficient of the twist equation in y^2 = x^3 + b_t
    -`u`: seed
    -`Fp6`: degree 6 extension over a prime field
    -`D_twist`: whether the twisted curve is a D-twist of an M-twist
    -`xi`: the residue in Fp so that Fp6 = Fp[x]/(x^6-xi)

    Optimized optimal ate Miller loop for BW6-BN
    with trace mod r mod u = 3
    trace = 18*x^3 + 18*x^2 + 9*x + 3 + r(x) * ht
    v = 3*u+1
    as frobenius(f_{2u,Q}(P)) * (f_2u)_{v,[u]Q}(P)
    with AKLGL formulas
    u, v in 2-naf
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.

    The same applies for v = 3*u+1
    """
    PP = (P[0], P[1])
    if xi is None:
        xi = -Fp6.modulus().constant_coefficient() # works with absolute and towering of extension
    if u < 0:
        m_2u, u2Q = miller_function_ate_2naf_aklgl((Q[0], -Q[1], 1), PP, b_t, -2*u, Fp6, D_twist=D_twist)
    else:
        m_2u, u2Q = miller_function_ate_2naf_aklgl((Q[0], Q[1], 1), PP, b_t, 2*u, Fp6, D_twist=D_twist)
    Z1 = 1/u2Q[2]
    u2Q = (u2Q[0]*Z1, u2Q[1]*Z1, 1)
    v = 3*u+1
    m_2u_inv = m_2u.conjugate()
    if v < 0:
        u2Q = (u2Q[0], -u2Q[1], 1)
        m_2uv, u2vQ = miller_function_ate_2naf_aklgl(u2Q, PP, b_t, -v, Fp6, m0=m_2u_inv, m0_inv=m_2u, D_twist=D_twist)
    else:
        m_2uv, u2vQ = miller_function_ate_2naf_aklgl(u2Q, PP, b_t, v, Fp6, m0=m_2u, m0_inv=m_2u_inv, D_twist=D_twist)
    l, u2v_1Q = add_line_h_a0_twist6_aklgl(u2vQ, (Q[0], Q[1]), PP, D_twist=D_twist)
    if D_twist:
        m_2uv_1 = sparse_mult_d6_twist(l[0], l[1], l[3], m_2uv, xi, Fp6)
    else:
        m_2uv_1 = sparse_mult_m6_twist(l[0], l[2], l[3], m_2uv, xi, Fp6)
    return m_2u.frobenius() * m_2uv_1

####### final exponentiation
def final_exp_hard_bw6_bn_trace_0_mod_r_mod_u(m, u, ht, hy):
    """
    c0 + c1*q = (2*u) * Phi_k(p)/r
    with cc0 = Phi_6(tr0-1)/(3*r)
    c0 = (ht^2+3*hy^2)/4*r*(6*u^2+2*u+1) + (ht-hy)/2*tr0*(6*u^2+2*u+1) + cc0*(6*u^2+2*u+1) - 3*u - 1
    c1 = (ht^2+3*hy^2)/4*r*(2*u) + 1*(ht-hy)/2*tr0*(2*u) + cc0*(2*u) - 1
    c0 + c1*q = (6*u^2+2*u+1 + (2*u)*q)*((ht^2+3*y^2)/4*r + (ht-hy)/2*tr0 + cc0) - 3*u - 1 - q

    6 exp(u) + exp(d2) + exp(d1) + 6 S + 15 M + 2 Frobenius + 3 conjugate

    r = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1
    t0 = -18*u^3 - 18*u^2 - 9*u
    c0 = 3*u^2 + 3*u + 1
    use the formulas:
    t0 == -3*(u*(2*c0 + 1))
    r == 2*(-u*t0 + c0) - 1
    """
    mu = m**u
    m0 = mu**2
    mc0 = m0 * mu * m # m**(3*u+1)
    mc = (mc0 * m.frobenius()).conjugate()
    m1 = (mc0**2)**u * m
    ma = m0.frobenius() * m1
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-hy)//2
    mac = ma**u
    mac = (mac * ma)**u     # ma^(u^2 + u)
    mac = mac**2 * mac * ma # ma^(3*u^2 + 3*u + 1)
    mat = mac**2 * ma
    mat = mat**u
    mat = mat**2 * mat
    mar = (mat**u * mac)**2 * ma.conjugate()
    mat = mat.conjugate()
    mb = mar**d2 * mat**d1 * mac
    return mb * mc

def final_exp_hard_bw6_bn_trace_0_mod_r_mod_u_alt(m, u, ht, hy):
    """
    c0 + c1*q = (6*u^2+4*u+1) * Phi_k(p)/r
    with cc0 = Phi_6(tr0-1)/(3*r)
    c0 = (ht^2+3*y^2)/4*r*(-2*u) + (ht-hy)/2*tr0*(-2*u) + cc0*(-2*u) + 1
    c1 = (ht^2+3*y^2)/4*r*(6*u^2+4*u+1) + (ht-hy)/2*tr0*(6*u^2+4*u+1) + cc0*(6*u^2+4*u+1) -3*u - 2
    c0 + c1*q = (-2*u + (6*u^2+4*u+1)*q)*((ht^2+3*y^2)/4*r + (ht-hy)/2*tr0 + cc0) + 1 - (3*u+2)*q

    6 exp(u) + exp(d2) + exp(d1) + 7 S + 15 M + 2 Frobenius + 4 conjugate

    r = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1
    t0 = -18*u^3 - 18*u^2 - 9*u
    c0 = 3*u^2 + 3*u + 1
    use the formulas:
    t0 == -3*(u*(2*c0 + 1))
    r == 2*(-u*t0 + c0) - 1
    """
    mu = m**u
    m0 = mu**2
    mc1 = m0 * mu * m**2 # m**(3*u+2)
    mc = (mc1.frobenius()).conjugate() * m
    m1 = (mc1**2)**u * m
    ma = m0.conjugate() * m1.frobenius()
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-hy)//2
    mac = ma**u
    mac = (mac * ma)**u     # ma^(u^2 + u)
    mac = mac**2 * mac * ma # ma^(3*u^2 + 3*u + 1)
    mat = mac**2 * ma
    mat = mat**u
    mat = mat**2 * mat
    mar = (mat**u * mac)**2 * ma.conjugate()
    mat = mat.conjugate()
    mb = mar**d2 * mat**d1 * mac
    return mb * mc

def final_exp_hard_bw6_bn_trace_3_mod_r_mod_u(m, u, ht, hy):
    """
    c0 + c1*q = 2*u * Phi_k(p)/r
    with cc3 = Phi_6(tr3-1)/(3*r)
    c0 = -(ht^2+3*hy^2)/4*r*(6*u^2+4*u+1) - (ht+hy)/2*tr3*(6*u^2+4*u+1) - cc3*(6*u^2+4*u+1) - (3*u+2)
    c1 = (ht^2+3*hy^2)/4*r*(2*u) + (ht+hy)/2*tr3*(2*u) + cc3*(2*u) + 1
    c0 + c1*q = (-(6*u^2+4*u+1) + 2*u*q)*((ht^2+3*y^2)/4*r + (ht+hy)/2*tr3 + cc3) - (3*u+2) + q

    6 exp(u) + exp(d2) + exp(d1) + 7 S + 16 M + 2 Frobenius + 2 conjugate

    r = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1
    t3 = 18*u^3 + 18*u^2 + 9*u + 3
    c3 = 3*u^2 + 3*u + 1
    with the formulas
    t3 == 3*(u*(2*c3 + 1) + 1)
    r == 2*u*(t3 + 3*u) + 1
    """
    mu = m**u
    m0 = mu**2
    mc0 = (m0 * mu * m**2).conjugate()
    mc = mc0 * m.frobenius()
    m1 = (mc0**2)**u * m.conjugate()
    ma = m1 * m0.frobenius()
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+hy)//2
    mac = ma**u
    ma3u = mac**2 * mac
    mac = ma3u**u     # ma^(u^2 + u)
    mac = mac * ma3u * ma # ma^(3*u^2 + 3*u + 1)
    mat = (mac**2 * ma)**u * ma
    mat = mat**2 * mat
    mar = ((mat * ma3u)**u)**2 * ma
    mb = mar**d2 * mat**d1 * mac
    return mb * mc

def final_exp_hard_bw6_bn_trace_3_mod_r_mod_u_alt(m, u, ht, hy):
    """
    c0 + c1*q = (6*u^2+2*u+1) * Phi_k(p)/r
    with cc3 = Phi_6(tr3-1)/(3*r)
    c0 = (ht^2+3*hy^2)/4*r*(2*u) + (ht+hy)/2*tr3*(2*u) + cc3*(2*u) + 1
    c1 = (ht^2+3*hy^2)/4*r*(6*u^2+2*u+1) + (ht+hy)/2*tr3*(6*u^2+2*u+1) + cc3*(6*u^2+2*u+1) + 3*u+1
    c0 + c1*q = (2*u + (6*u^2+2*u+1)*q)*((ht^2+3*y^2)/4*r + (ht+hy)/2*tr3 + cc3) + 1 + (3*u+1)*q

    6 exp(u) + exp(d2) + exp(d1) + 6 S + 16 M + 2 Frobenius

    r = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1
    t3 = 18*u^3 + 18*u^2 + 9*u + 3
    c3 = 3*u^2 + 3*u + 1
    with the formulas
    t3 == 3*(u*(2*c3 + 1) + 1)
    r == 2*u*(t3 + 3*u) + 1
    """
    mu = m**u
    m0 = mu**2
    mc1 = (m0 * mu * m) # m**(3*u+1)
    mc = mc1.frobenius() * m
    m1 = (mc1**2)**u * m
    ma = m0 * m1.frobenius()
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+hy)//2
    mac = ma**u
    ma3u = mac**2 * mac
    mac = ma3u**u     # ma^(u^2 + u)
    mac = mac * ma3u * ma # ma^(3*u^2 + 3*u + 1)
    mat = (mac**2 * ma)**u * ma
    mat = mat**2 * mat
    mar = ((mat * ma3u)**u)**2 * ma
    mb = mar**d2 * mat**d1 * mac
    return mb * mc

####### co-factor multiplication on G1
def bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-2) + Phi_6(t0-1)/(3*r) - hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(2u)        + (ht+3*hy)/2*u     + hy
    c1 = (ht^2+3*hy^2)/4*(6u^2+4u+1) - (ht-3*hy)/2*(u+1) - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t0-1) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = d2*u2P
    Q2 = (3*u+2)*P2 + d2*P
    P1 = d1*uP
    P0 = hy*P
    Q1 = -d0*(uP+P)
    Q0 = -P0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-2) + Phi_6(t0-1)/(3*r) - hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P):
    c0 = (ht^2+3*hy^2)/4*(6u^2+2u+1) - (ht+hy)/2*(2u+1) + hy*u
    c1 = (ht^2+3*hy^2)/4*(-2u) - (ht+3*hy)/2*u - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t0-1) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht+hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = -d2*u2P
    P2 = (-3*u-1)*Q2 + d2*P
    P1 = -d0*(u2P+P)
    P0 = hy*uP
    Q1 = -d1*uP
    Q0 = -hy*P
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u_alt(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-2) + Phi_6(t0-1)/(3*r) - hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6u^2+4u+1) - (ht-3*hy)/2*(u+1) - hy
    c1 = (ht^2+3*hy^2)/4*2u + (ht+3*hy)/2*u + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t0-2) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = d2*u2P
    P2 = (3*u+2)*Q2 + d2*P
    P1 = -d0*(uP+P)
    Q1 = d1*uP
    Q0 = hy*P
    P0 = -Q0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_0_mod_r_u_alt_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-2) + Phi_6(t0-1)/(3*r) - hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(-2*u)        - (ht+3*hy)/2*u       - hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) - (ht+hy)  /2*(2*u+1) + hy*u
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t0-2) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht+hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = -d2*u2P
    Q2 = -(3*u+1)*P2 + d2*P
    P1 = -d1*uP
    P0 = -hy*P
    Q1 = -d0*(u2P+P)
    Q0 = hy*uP
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-2) + Phi_6(t3-1)/(3*r) + hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(-2*u)        - (ht-3*hy)/2* (u+1) - hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) + (ht+3*hy)/2* u     + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t3-1) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = -d2*u2P
    Q2 = -(3*u+1)*P2 + d2*P
    Q1 = d1*uP
    Q0 = hy*P
    P1 = -d0*(uP+P)
    P0 = -Q0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-2) + Phi_6(t3-1)/(3*r) + hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) + (ht+  hy)/2 *(2*u+1) - hy*u
    c1 = (ht^2+3*hy^2)/4*(2*u)         + (ht-3*hy)/2 *(u+1)   + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t3-1) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-3*hy)//2
    d0 = (ht+hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = d2*u2P
    P2 = (3*u+2)*Q2 + d2*P
    Q1 = d1*(uP+P)
    Q0 = hy*P
    P1 = d0*(u2P+P)
    P0 = -hy*uP
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u_alt(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-2) + Phi_6(t3-1)/(3*r) + hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) + (ht+3*hy)/2* u     + hy
    c1 = (ht^2+3*hy^2)/4*(-2*u)        - (ht-3*hy)/2* (u+1) - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t3-2) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = -d2*u2P
    P2 = -(3*u+1)*Q2 + d2*P
    P1 = d1*uP
    P0 = hy*P
    Q1 = -d0*(uP+P)
    Q0 = -P0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_cofactor_trace_3_mod_r_u_alt_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c where c*r = #E(Fp)
    c = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-2) + Phi_6(t3-1)/(3*r) + hy
    fast cofactor multiplication c*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(2*u)         + (ht-3*hy)/2* (u+1)   + hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) + (ht+  hy)/2* (2*u+1) - hy*u
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t3-2) mod r
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-3*hy)//2
    d0 = (ht+hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = d2*u2P
    Q2 = (3*u+2)*P2 + d2*P
    P1 = d1*(uP+P)
    P0 = hy*P
    Q1 = d0*(u2P+P)
    Q0 = -hy*uP
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

####### membership testing on G1
def bw6_bn_g1_mult_by_r_trace_0_mod_r_u(P, omega, u):
    """
    fast multiplication by r
    lambda = -q mod r = 18*u^3+18*u^2+9*u+1
    (6*u^2+4*u+1) + (18*u^3+18*u^2+9*u+1) * (2*u) = r
    """
    Q = 2*u*P
    return (3*u+2)*Q + P + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_r_trace_0_mod_r_u_alt(P, omega, u):
    """
    fast multiplication by r
    lambda = q-1 mod r = -18*u^3-18*u^2-9*u-2
    (6*u^2+2*u+1) + (-18*u^3-18*u^2-9*u-2) * (-2*u) = r
    """
    Q = -2*u*P
    return -(3*u+1)*Q + P + bw6_phi(Q, omega)

def bw6_bn_g1_check_membership_trace_0_mod_r_u(P, omega, u):
    R = bw6_bn_g1_mult_by_r_trace_0_mod_r_u(P, omega, u)
    return R == P.curve()(0)

def bw6_bn_g1_check_membership_trace_0_mod_r_u_alt(P, omega, u):
    R = bw6_bn_g1_mult_by_r_trace_0_mod_r_u_alt(P, omega, u)
    return R == P.curve()(0)

def bw6_bn_g1_mult_by_r_trace_3_mod_r_u(P, omega, u):
    """
    fast multiplication by r
    lambda = -q mod r = -18*u^3-18*u^2-9*u-2
    (6*u^2+2*u+1) + (-18*u^3-18*u^2-9*u-2) * (-2*u) = r
    """
    Q = -2*u*P
    return -(3*u+1)*Q + P + bw6_phi(Q, omega)

def bw6_bn_g1_mult_by_r_trace_3_mod_r_u_alt(P, omega, u):
    """
    fast multiplication by r
    lambda = q-1 mod r = 18*u^3+18*u^2+9*u+1
    (6*u^2+4*u+1) + (18*u^3+18*u^2+9*u+1) * (2*u) = r
    """
    Q = 2*u*P
    return (3*u+2)*Q + P + bw6_phi(Q, omega)

def bw6_bn_g1_check_membership_trace_3_mod_r_u(P, omega, u):
    R = bw6_bn_g1_mult_by_r_trace_3_mod_r_u(P, omega, u)
    return R == P.curve()(0)

def bw6_bn_g1_check_membership_trace_3_mod_r_u_alt(P, omega, u):
    R = bw6_bn_g1_mult_by_r_trace_3_mod_r_u_alt(P, omega, u)
    return R == P.curve()(0)

####### co-factor multiplication on G2
def bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-1) + Phi_6(t0-1)/(3*r) + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t0-1) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) - (ht-3*hy)/2*u      + hy
    c1 = (ht^2+3*hy^2)/4*(-2*u)        + (ht+3*hy)/2*(u+1)  - hy
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-3*hy)//2
    d0 = (ht+3*hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = -d2*u2P
    P2 = -(3*u+1)*Q2 + d2*P
    P1 = -d1*uP
    P0 = hy*P
    Q1 = d0*(uP+P)
    Q0 = -P0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-1) + Phi_6(t0-1)/(3*r) + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t0-1) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(2*u)         - (ht+3*hy)/2*(u+1) + hy
    #c0= (ht^2+3*hy^2)/4*(2*u)         - (ht+  hy)/2 - (ht+3*hy)/2*u
    c1 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) - (ht-  hy)/2       - ht*u
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = d2*u2P
    Q2 = (3*u+2)*P2 + d2*P
    P1 = -d1*(uP+P)
    P0 = hy*P
    Q1 = -d0*P
    Q0 = -ht*uP
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_alt(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-1) + Phi_6(t0-1)/(3*r) + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t0-2) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(-2*u)        + (ht+3*hy)/2*(u+1) - hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) - (ht-3*hy)/2*u     + hy
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = -d2*u2P
    Q2 = -(3*u+1)*P2 + d2*P
    P1 = d1*(uP+P)
    P0 = -hy*P
    Q1 = -d0*uP
    Q0 = -P0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_0_mod_r_u_alt_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht-hy)/2*(t0-1) + Phi_6(t0-1)/(3*r) + hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t0-2) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) - (ht-  hy)/2       - ht*u
    c1 = (ht^2+3*hy^2)/4*(2*u)         - (ht+3*hy)/2*(u+1) + hy
    #c1= (ht^2+3*hy^2)/4*(2*u)         - (ht+  hy)/2       - (ht+3*hy)/2*u
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-hy)//2
    d0 = (ht+3*hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = d2*u2P
    P2 = (3*u+2)*Q2 + d2*P
    P1 = -d1*P
    P0 = -ht*uP
    Q1 = -d0*(uP+P)
    Q0 = hy*P
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-1) + Phi_6(t3-1)/(3*r) - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t3-1) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(-2*u)        + (ht-3*hy)/2*u - hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) + (ht-hy)/2     + ht*u
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-3*hy)//2
    d0 = (ht-hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = -d2*u2P
    Q2 = -(3*u+1)*P2 + d2*P
    P1 = d1*uP
    P0 = -hy*P
    Q1 = d0*P
    Q0 = ht*uP
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-1) + Phi_6(t3-1)/(3*r) - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    -(t3-1) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) + (ht+3*hy)/2*(u+1) - hy
    c1 = (ht^2+3*hy^2)/4*(2*u)         - (ht-3*hy)/2*u     + hy
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht+3*hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = d2*u2P
    P2 = (3*u+2)*Q2 + d2*P
    P1 = d1*(uP+P)
    P0 = -hy*P
    Q1 = -d0*uP
    Q0 = -P0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_alt(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-1) + Phi_6(t3-1)/(3*r) - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t3-2) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(2*u)         - (ht-3*hy)/2*u     + hy
    c1 = (ht^2+3*hy^2)/4*(6*u^2+4*u+1) + (ht+3*hy)/2*(u+1) - hy
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-3*hy)//2
    d0 = (ht+3*hy)//2
    uP = u*P
    u2P = 2*uP
    P2 = d2*u2P
    Q2 = (3*u+2)*P2 + d2*P
    P1 = -d1*uP
    P0 = hy*P
    Q1 = d0*(uP+P)
    Q0 = -P0
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

def bw6_bn_g2_mult_by_cofactor_trace_3_mod_r_u_alt_(P, omega, u, ht, hy):
    """
    multiply point P by a multiple of the cofactor c2 where c2*r = #E'(Fp)
    c2 = (ht^2+3*hy^2)/4*r + (ht+hy)/2*(t3-1) + Phi_6(t3-1)/(3*r) - hy
    assume that the endomorphism phi (x,y) -> (omega*x,y) has eigenvalue
    (t3-2) mod r
    fast cofactor multiplication c2*P = c0*P + c1*phi(P) with
    c0 = (ht^2+3*hy^2)/4*(6*u^2+2*u+1) + (ht-  hy)/2       + ht*u
    c1 = (ht^2+3*hy^2)/4*(-2*u)        + (ht-3*hy)/2*u     - hy
    """
    d2 = (ht**2+3*hy**2)//4
    d1 = (ht-hy)//2
    d0 = (ht-3*hy)//2
    uP = u*P
    u2P = 2*uP
    Q2 = -d2*u2P
    P2 = -(3*u+1)*Q2 + d2*P
    P1 = d1*P
    P0 = ht*uP
    Q1 = d0*uP
    Q0 = -hy*P
    Q = Q0 + Q1 + Q2
    return P2 + P1 + P0 + bw6_phi(Q, omega)

"""
####### membership testing on G2
def bw6_bn_g2_mult_by_r_trace_0_mod_r_u(P, omega, u0):
def bw6_bn_g2_check_membership_trace_0_mod_r_u(P, omega, u0):
def bw6_bn_g2_mult_by_r_trace_3_mod_r_u(P, omega, u0):
def bw6_bn_g2_check_membership_trace_3_mod_r_u(P, omega, u0):
"""
