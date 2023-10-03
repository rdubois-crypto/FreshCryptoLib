from external.Pairings.pairing import *
from external.Pairings.pairing_cocks_pinch import miller_function_tate_aklgl_a0_2_parallel, miller_function_tate_a0_twist6_aklgl_2_multi_scalar

"""
Tate, ate, optimal Tate, optimal ate pairings on BN curves

"""
def miller_loop_ate_bn(Q, P, u):
    """
    Return f_{6*u^2,Q}(P) where t(u)-1 = 6*u^2

    INPUT:
    - `Q`: r-torsion point on G2 not compressed, in E(Fpk) or E(Fqd)
    - `P`: r-torsion point on G1 in E(Fp)
    - `u`: seed for the curve coefficients

    Return a value m in the field of coefficients of the coordinates of Q
    """
    T = 6*u**2 # always positive
    return miller_function_ate(Q, P, 0, T)

def miller_loop_opt_ate_bn(Q, P, u):
    """
    Return f_{6*u+2,Q}(P) l_{[6u+2]Q,pi(Q)}(P) l_{[6u+2]Q+pi(Q), pi^2(-Q)}(P)

    INPUT:
    - `Q`: r-torsion point on G2 not compressed, in E(Fpk) with absolute extension Fpk allowing frobenius
    - `P`: r-torsion point on G1 in E(Fp)
    - `u`: seed for the curve coefficients

    The Frobenius map in Sage can be used only on absolute extensions.
    assume Q has coordinates in Fp12 as an absolute extension over Fp.

    Alternative formula with one more Frobenius on G2 instead to compute pi^3(Q):
    f_{6*u+2,Q}(P) l_{[6u+2]Q,pi(Q)}(P) l_{-pi^2(Q), pi^3(Q)}(P)
    """
    T = 6*u+2
    if T < 0:
        m, S = miller_function_ate((Q[0], -Q[1]), P, 0, -T)
    else:
        m, S = miller_function_ate(Q, P, 0, T)
    pQ = (Q[0].frobenius(), Q[1].frobenius())
    ln1, R1 = add_line_j(S, pQ, (P[0], P[1]))
    p2_Q = (Q[0].frobenius(2), -Q[1].frobenius(2))
    # alternative formula:
    # p3Q = (p2_Q[0].frobenius(), -p2_Q[1].frobenius())
    # ln2, R2 = add_line_affine_j(p2_Q, p3Q, (P[0],P[1]))
    # return m * ln1 * ln2
    ln2, R2 = add_line_j(R1, p2_Q, (P[0],P[1]))
    return m * ln1 * ln2

def miller_loop_opt_ate_bn_2naf(Q, P, u):
    """
    Return f_{6*u+2,Q}(P) l_{[6u+2]Q,pi(Q)}(P) l_{[6u+2]Q+pi(Q), pi^2(-Q)}(P)
    with the Miller loops in 2-naf form

    INPUT:
    - `Q`: r-torsion point on G2 not compressed, in E(Fpk) with absolute extension Fpk allowing frobenius
    - `P`: r-torsion point on G1 in E(Fp)
    - `u`: seed for the curve coefficients

    The Frobenius map in Sage can be used only on absolute extensions.
    assume Q has coordinates in Fp12 as an absolute extension over Fp.

    Alternative formula with one more Frobenius on G2 instead to compute pi^3(Q):
    f_{6*u+2,Q}(P) l_{[6u+2]Q,pi(Q)}(P) l_{-pi^2(Q), pi^3(Q)}(P)
    """
    T = 6*u+2
    if T < 0:
        m, S = miller_function_ate_2naf((Q[0], -Q[1]), P, 0, -T)
    else:
        m, S = miller_function_ate_2naf(Q, P, 0, T)
    pQ = (Q[0].frobenius(), Q[1].frobenius())
    ln1, R1 = add_line_j(S, pQ, (P[0], P[1]))
    p2_Q = (Q[0].frobenius(2), -Q[1].frobenius(2))
    # alternative formula:
    # p3Q = (p2_Q[0].frobenius(), -p2_Q[1].frobenius())
    # ln2, R2 = add_line_affine_j(p2_Q, p3Q, (P[0],P[1]))
    # return m * ln1 * ln2
    ln2, R2 = add_line_j(R1, p2_Q, (P[0], P[1]))
    return m * ln1 * ln2

def miller_loop_opt_tate_bn(P, Q, u):
    """
    Return f_{2*u+1,P}(Q) (f_{6*u^2+2*u, P}(Q))^(p^2)

    INPUT:
    - `P`: r-torsion point on G1 in E(Fp)
    - `Q`: r-torsion point on G2 not compressed, in E(Fpk) with absolute extension Fpk allowing frobenius
    - `u`: seed of the curve parameters

    The Frobenius map in Sage can be used only on absolute extensions.
    assume Q has coordinates in Fp12 as an absolute extension over Fp.

    Assuming that psi: (x,y) -> (omega*x, -y) has eigenvalue (tx-1)^2 mod r = -36*x^3-18*x^2-6*x-1
    it is the same as f_{2*u+1,P}(Q) f_{6*u^2+2*u, psi(P)}(Q)
    Computation can be improved as:
    m = f_{2*u, P}(Q); return m * l_{[2*u]P, P}(Q) ((m)_{3*u+1, [2*u]P}(Q))^(p^2)
    """
    if u < 0:
        m2u, P2u = miller_function_tate((P[0],-P[1]),Q,0,-2*u)
    else:
        m2u, P2u = miller_function_tate(P,Q,0,2*u)
    v = 3*u+1
    # affine coordinates for [2*u]*P
    z = 1/P2u[2]
    z2 = z**2
    z3 = z2 * z
    P2u = (P2u[0]*z2, P2u[1]*z3)
    ln, P2u1 = add_line_affine_j((P2u[0],P2u[1]),(P[0],P[1]),(Q[0],Q[1]))
    if v < 0:
        # invert m2u because m2u^v = m2u^(-|v|) = (1/m2u)^|v|
        mv, Pv = miller_function_tate((P2u[0],-P2u[1]),Q,0,-v, m0=m2u.frobenius(6))
    else:
        mv, Pv = miller_function_tate(P2u,Q,0,v, m0=m2u)
    return m2u * ln * mv.frobenius(2)

def miller_loop_opt_tate_bn_alt(P, Q, u):
    """
    Return f_{6*u^2+4*u+1,P}(Q) (f_{-2*u-1, psi(P)}(Q))^(p^2)

    INPUT:
    - `P`: r-torsion point on G1 in E(Fp)
    - `Q`: r-torsion point on G2 not compressed, in E(Fpk) with absolute extension Fpk allowing frobenius
    - `u`: seed of the curve parameters

    The Frobenius map in Sage can be used only on absolute extensions.
    assume Q has coordinates in Fp12 as an absolute extension over Fp.

    Assuming that psi: (x,y) -> (omega*x, -y) has eigenvalue (tx-1)^2 mod r = -36*x^3-18*x^2-6*x-1
    it is the same as f_{6*u^2+4*u+1,P}(Q) f_{-2*u-1, psi(P)}(Q)
    Computation can be improved as (with inversion as powering to p^6)
    m = f_{2*u, P}(Q); return (m * l_{[2*u]P, P})^(p^8) (m)_{3*u+2, [2*u]P}(Q) l_{[6u^2+4u]P, P}(Q)
    """
    if u < 0:
        m2u, P2u = miller_function_tate((P[0],-P[1]),Q,0,-2*u)
    else:
        m2u, P2u = miller_function_tate(P,Q,0,2*u)
    v = 3*u+2
    # affine coordinates for [2*u]*P
    z = 1/P2u[2]
    z2 = z**2
    z3 = z2 * z
    P2u = (P2u[0]*z2, P2u[1]*z3)
    ln, P2u1 = add_line_affine_j((P2u[0],P2u[1]),(P[0],P[1]),(Q[0],Q[1]))
    if v < 0:
        # invert m2u because m2u^v = m2u^(-|v|) = (1/m2u)^|v|
        mv, Pv = miller_function_tate((P2u[0],-P2u[1]),Q,0,-v, m0=m2u.frobenius(6))
    else:
        mv, Pv = miller_function_tate(P2u,Q,0,v, m0=m2u)
    ln2, Pv1 = add_line_j(Pv,(P[0],P[1]),(Q[0],Q[1]))
    return (m2u * ln).frobenius(8) * (mv * ln2)

def miller_loop_opt_tate_bn_aklgl_a0(P, Q, b, u, Fq6, map_Fq6_Fpk, D_twist=False):
    """
    INPUT:
    - `P`: on E (for G1)
    - `Q`: on E2 over Fp2 (the twist, for G2)
    - `b`: curve coefficient of E (to add and double P)
    - `u`: seed of the curve parameters
    - `Fq6`: degree 6 extension over Fq
    - `D_twist`: True/False (about E2)

    Assuming that psi: (x,y) -> (omega*x, -y) has eigenvalue (tx-1)^2 mod r = -36*x^3-18*x^2-6*x-1
    Return f_{2*u+1,P}(Q) f_{6*u^2+2*u, psi(P)}(Q)
    Alt formula: f_{6*u^2+4*u+1,P}(Q) f_{-2*u-1, psi(P)}(Q)

    can be improved as
    m = f_{2*u, P}(Q)
    Return m * l_{[2*u]P, P}(Q) ((m)_{3*u+1, [2*u]P}(Q))^(p^2)

    alt formula can be improved as (with inversion as powering to p^6)
    m = f_{2*u, P}(Q)
    Return (m * l_{[2*u]P, P})^(p^8) (m)_{3*u+2, [2*u]P}(Q) l_{[6u^2+4u]P, P}(Q)
    remember to swap the value D_twist for Tate (compared to ate)
    """
    if u < 0:
        m2u, P2u = miller_function_tate_aklgl((P[0],-P[1]),Q,b,-2*u,Fq6,D_twist=D_twist)
    else:
        m2u, P2u = miller_function_tate_aklgl(P,Q,b,2*u,Fq6,D_twist=D_twist)
    v = 3*u+1
    # affine coordinates for [2*u]*P from homogeneous coord
    z = 1/P2u[2]
    P2u = (P2u[0]*z, P2u[1]*z)
    ln, P2u1 = add_line_h_a0_twist6_aklgl((P2u[0],P2u[1],1), (P[0],P[1]), (Q[0],Q[1]), D_twist=not D_twist)
    if v < 0:
        # invert m2u because m2u^v = m2u^(-|v|) = (1/m2u)^|v|
        # m2u.frobenius(6) is not available because m2u is in a relative extension, Sage does not know how to do it
        coeffs_m2u = m2u.list() # a list of 6 coefficients, each coefficient is in Fq = Fp2 where Frobenius stands
        for i in range(1, len(coeffs_m2u), 2):
            coeffs_m2u[i] = -coeffs_m2u[i]
        m2u_p6 = Fq6(coeffs_m2u)
        #m2u_inv = 1/m2u
        mv, Pv = miller_function_tate_aklgl((P2u[0],-P2u[1]),Q,b,-v,Fq6,D_twist=D_twist,m0=m2u_p6)
    else:
        mv, Pv = miller_function_tate_aklgl(P2u,Q,b,v,Fq6,D_twist=D_twist,m0=m2u)
    xi = -Fq6.modulus().constant_coefficient()
    if not D_twist:
        m2u_ln = sparse_mult_d6_twist(ln[0],ln[1],ln[3], m2u, xi, Fq6)
    else:
        m2u_ln = sparse_mult_m6_twist(ln[0],ln[2],ln[3], m2u, xi, Fq6)
    return map_Fq6_Fpk(m2u_ln) * (map_Fq6_Fpk(mv)).frobenius(2)

def miller_loop_opt_tate_bn_aklgl_a0_2_parallel(P, Q, u, b, Fq6, map_Fq6_Fpk, D_twist=False, omega=None, xi=None):
    """Miller loop with AKLGL formulas, no optimisation (two independent Miller loops, nothing shared)

    INPUT:
    - `P`: r-torsion point in G1 on E(Fp)
    - `Q`: r-torsion point in compressed G2 on E2(Fq) (twist curve)
    - `b`: E2 curve coefficient (for Q)
    - `omega`: Fp-element s.t. omega^2 + omega + 1 = 0 mod p, so that
               lambda = p^2 mod r is the eigenvalue of psi: (x,y) -> (omega*x, -y) on G1

    RETURN: f_{a0+a1*lambda,P}(Q) = f_{a0, P}(Q) * f_{a1, psi(P)}(Q)
    where [a0+a1*lambda]P = a0 + a1*psi(P) = O 
    a0 = 2*u+1, a1 = 6*u^2+2*u

    multiply the sparse lines or tangent two by two before accumulating
    """
    a0 = 2*u+1
    a1 = 6*u**2+2*u
    return miller_function_tate_aklgl_a0_2_parallel(P, Q, b, a0, a1, Fq6, map_Fq6_Fpk, D_twist=D_twist, omega=omega, xi=xi)

def miller_loop_opt_tate_bn_a0_twist6_aklgl_2_multi_scalar(P, Q, b, u, omega, xi, Fq6, D_twist=False):
    """Miller loop with AKLGL formulas and multi-scalar technique

    INPUT:
    - `P`: point on E(Fp)
    - `Q`: point on E2(Fq) (twist curve)
    - `b`: E2 curve coefficient (for Q)
    - `omega`: Fp-element s.t. omega^2 + omega + 1 = 0 mod p, so that
               lambda = p^4 mod r = -p^2-1 mod r is the eigenvalue of
                psi: (x,y) -> (omega*x, y) on G1

    RETURN: f_{a0+a1*lambda,P}(Q) where [a0+a1*lambda]P = a0 + a1*psi(P) = O
    a0 = 2*u+1, a1 = -(6*u^2+2*u)

    multiply the sparse lines or tangent two by two before accumulating
    """
    a0 = 2*u+1
    a1 = -(6*u**2+2*u)
    return miller_function_tate_a0_twist6_aklgl_2_multi_scalar(P, Q, b, a0, a1, omega, xi, Fq6, D_twist=D_twist)

def frobenius_map_G2(Q, xi_p13, xi_p12):
    """
    map (xQ,yQ) -> (xQ^p * xi^((p-1)/3), yQ^p * xi^((p-1)/2)) = (p % r)*Q
    INPUT:
    -`Q`: point on E'(Fq) (the compressed form of G2) of order r
    For D_twist:
    -`xi_p13`: xi^((p-1)/3)
    -`xi_p12`: xi^((p-1)/2)
    For M-twist:
    -`xi_p13`: 1/xi^((p-1)/3)
    -`xi_p12`: 1/xi^((p-1)/2)

    Return: a point pi(Q) on E'(Fq)
    """
    pQ = ((Q[0]).frobenius() * xi_p13, (Q[1]).frobenius() * xi_p12)
    return Q.curve()(pQ)

def miller_loop_opt_ate_bn_aklgl(Q,P,b_t,u,Fq6,D_twist=False,xi=None):
    """return f_{6u+2,Q}(P) * l_{[6u+2]Q,pi(Q)}(P) * l_{[6u+2]Q+pi(Q), pi^2(-Q)}(P)

    Optimized optimal ate Miller loop for BN curves and AKLGL formulas
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.
    """
    if xi is None:
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    P_ = (P[0], P[1])
    T = 6*u + 2
    if T < 0:
        Q_ = (Q[0], -Q[1], 1)
        m_T, TQ = miller_function_ate_aklgl(Q_, P_, b_t, -T, Fq6, D_twist=D_twist)
    else:
        Q_ = (Q[0], Q[1], 1)
        m_T, TQ = miller_function_ate_aklgl(Q_, P_, b_t, T, Fq6, D_twist=D_twist)
    # pQ with Q in E(Fp2) (before applying the twist)
    # then the output is, with xi^((p-1)/3), xi^((p-1)/2)
    p = Fq6.characteristic()
    xi_p16 = xi**((p-1)//6)
    if not D_twist:
        xi_p16 = 1/xi_p16
    xi_p13 = xi_p16**2
    xi_p12 = xi_p13 * xi_p16
    pQ = ((Q[0]).frobenius() * xi_p13, (Q[1]).frobenius() * xi_p12)
    # xi^((p^2-1)/6) = xi^((p-1)/6)*xi^(p+1)
    xi_p216 = xi_p16.frobenius() * xi_p16
    xi_p213 = xi_p216**2
    xi_p212 = xi_p213 * xi_p216
    p2_Q = (Q[0].frobenius(2) * xi_p213, -Q[1].frobenius(2) * xi_p212)

    l1, TpQ = add_line_h_a0_twist6_aklgl(TQ, pQ, P_, D_twist=D_twist)
    l2, Tp_p2Q = add_line_h_a0_twist6_aklgl(TpQ, p2_Q, P_, D_twist=D_twist)
    if D_twist:
        l1l2 = sparse_sparse_mult_d6_twist(l1[0], l1[1], l1[3], l2[0], l2[1], l2[3], xi, Fq6)
    else:
        l1l2 = sparse_sparse_mult_m6_twist(l1[0], l1[2], l1[3], l2[0], l2[2], l2[3], xi, Fq6)
    return m_T * l1l2

def miller_loop_opt_ate_bn_aklgl_2naf(Q,P,b_t,u,Fq6,D_twist=False,xi=None):
    """return f_{6u+2,Q}(P) * l_{[6u+2]Q,pi(Q)}(P) * l_{[6u+2]Q+pi(Q), pi^2(-Q)}(P)

    Optimized optimal ate Miller loop for BN curves and AKLGL formulas, 2-NAF
    Q, P are r-torsion points in G2, G1, in affine coordinates,
    u is the seed (integer) of the curve coefficients

    If u < 0, then f_{|u|, -Q}(P) is computed thanks to the formula
    f_{ij,Q} = f_{i,Q}^j*f_{j,[i]Q} and with i=-1, j=|u|:
    f_{-|u|,Q} = f_{-1,Q}^u*f_{|u|,-Q} and since f_{-1,Q} is a vectical line,
    it is discarded: f_{-|u|,Q} = f_{|u|,-Q}.
    """
    if xi is None:
        xi = -Fq6.modulus().constant_coefficient() # works with absolute and towering of extension
    P_ = (P[0], P[1])
    T = 6*u + 2
    if T < 0:
        Q_ = (Q[0], -Q[1], 1)
        m_T, TQ = miller_function_ate_2naf_aklgl(Q_, P_, b_t, -T, Fq6, D_twist=D_twist)
    else:
        Q_ = (Q[0], Q[1], 1)
        m_T, TQ = miller_function_ate_2naf_aklgl(Q_, P_, b_t, T, Fq6, D_twist=D_twist)
    # pQ with Q in E(Fp2) (before applying the twist)
    # then the output is, with xi^((p-1)/3), xi^((p-1)/2)
    p = Fq6.characteristic()
    xi_p16 = xi**((p-1)//6)
    if not D_twist:
        xi_p16 = 1/xi_p16
    xi_p13 = xi_p16**2
    xi_p12 = xi_p13 * xi_p16
    pQ = ((Q[0]).frobenius() * xi_p13, (Q[1]).frobenius() * xi_p12)
    # xi^((p^2-1)/6) = xi^((p-1)/6)*xi^(p+1)
    xi_p216 = xi_p16.frobenius() * xi_p16
    xi_p213 = xi_p216**2
    xi_p212 = xi_p213 * xi_p216
    p2_Q = (Q[0].frobenius(2) * xi_p213, -Q[1].frobenius(2) * xi_p212)

    l1, TpQ = add_line_h_a0_twist6_aklgl(TQ, pQ, P_, D_twist=D_twist)
    l2, Tp_p2Q = add_line_h_a0_twist6_aklgl(TpQ, p2_Q, P_, D_twist=D_twist)
    if D_twist:
        l1l2 = sparse_sparse_mult_d6_twist(l1[0], l1[1], l1[3], l2[0], l2[1], l2[3], xi, Fq6)
    else:
        l1l2 = sparse_sparse_mult_m6_twist(l1[0], l1[2], l1[3], l2[0], l2[2], l2[3], xi, Fq6)
    return m_T * l1l2

def final_exp_hard_bn(m, u):
    """
    Final exponentiation, hard part: m^(h*(p^4-p^2+1)//r) where h(x) = 12*x^3 + 6*x^2 + 2*x

    INPUT:
    - `m`: in GF(p^12) as an absolute extension to allow for m.frobenius()
    - `u`: BN curve seed

    Formulas from Fuentes-Castaneda et al, SAC'2011
    l0 = 1+6*x + 12*x**2 + 12*x**3
    l1 = 4*x + 6*x**2 + 12*x**3
    l2 = 6*x + 6*x**2 + 12*x**3
    l3 = -1 + 4*x + 6*x**2 + 12*x**3

    l0 + l1*px + l2*px^2 + l3*px3 == (12*x^3 + 6*x^2 + 2*x) * (px^4 - px^2 + 1)//rx

    cost: 3 exp(u) + 3 S + 10 M + 3 frob
    """
    fx = m**u
    f2x = fx**2
    f4x = f2x**2
    f6x = f4x * f2x
    f6x2 = f6x**u
    f12x2 = f6x2**2
    f12x3 = f12x2**u
    a = f12x3 * f6x2 * f6x
    b = a * f2x.frobenius(6)
    res = (a * f6x2 * m) * b.frobenius() * a.frobenius(2) * (b * m.frobenius(6)).frobenius(3)
    return res
