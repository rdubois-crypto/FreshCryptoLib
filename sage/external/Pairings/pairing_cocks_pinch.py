# pairing computation for Cocks-Pinch curves, for example the 2nd curves in a 2-chain from BLS-12 or BLS-24.
from external.Pairings.pairing import *

def precompute_qj(Q0, Q1, Q2, Q3):
    """
    Precompute all possible sums of points in {Q0, Q1, Q2, Q3}

    INPUT:
    - `Qi`: four points on the same curve

    For example Q0 = Q, Q1 = pi_p(Q), Q2 = pi_p2(Q), Q3 = pi_p3(Q)
    Compute all sums of 2 to 4 distinct points.
    There are 15 such points.
    """
    # precomputations
    Qj = [None]*15     # l3 l2 l1 l0
    Qj[1-1] = Q0       # 0001
    Qj[2-1] = Q1       # 0010
    Q01 = add_cln_b0_affine(Q0, Q1)
    Qj[3-1] = Q01      # 0011
    Qj[4-1] = Q2       # 0100
    Q02 = add_cln_b0_affine(Q0, Q2)
    Qj[5-1] = Q02      # 0101
    Q12 = add_cln_b0_affine(Q1, Q2)
    Qj[6-1] = Q12      # 0110
    Q012 = add_cln_b0_affine(Q0, Q12)
    Qj[7-1] = Q012     # 0111
    Qj[8-1] = Q3       # 1000
    Q03 = add_cln_b0_affine(Q0, Q3)
    Qj[9-1] = Q03      # 1001
    Q13 = add_cln_b0_affine(Q1, Q3)
    Qj[10-1] = Q13     # 1010
    Q013 = add_cln_b0_affine(Q01, Q3)
    Qj[11-1] = Q013    # 1011
    Q23 = add_cln_b0_affine(Q2, Q3)
    Qj[12-1] = Q23     # 1100
    Q023 = add_cln_b0_affine(Q02, Q3)
    Qj[13-1] = Q023    # 1101
    Q123 = add_cln_b0_affine(Q12, Q3)
    Qj[14-1] = Q123    # 1110
    Q0123 = add_cln_b0_affine(Q012, Q3)
    Qj[15-1] = Q0123   # 1111
    return Qj

def precompute_qj_with_lines(Q0, Q1, Q2, Q3, PP):
    """
    Precompute all possible sums of points in {Q0, Q1, Q2, Q3} and product of lines evaluated at PP

    INPUT:
    - `Qi`: four points on the same curve

    For example Q0 = Q, Q1 = pi_p(Q), Q2 = pi_p2(Q), Q3 = pi_p3(Q)
    Compute all sums of 2 to 4 distinct points.
    There are 15 such points.
    Compute the product of the lines involved in the sum, evaluated at PP.

    Needed in miller_function_ate_cln_b0_4_multi_scalar
    """
    # precomputations
    Qj = [None]*15             # l3 l2 l1 l0
    Qj[1-1] = [Q0, 1]          # 0001
    Qj[2-1] = [Q1, 1]          # 0010
    ln01, Q01 = add_line_cln_b0_affine(Q0, Q1, PP)
    Qj[3-1] = [Q01, ln01]      # 0011
    Qj[4-1] = [Q2, 1]          # 0100
    ln02, Q02 = add_line_cln_b0_affine(Q0, Q2, PP)
    Qj[5-1] = [Q02, ln02]      # 0101
    ln12, Q12 = add_line_cln_b0_affine(Q1, Q2, PP)
    Qj[6-1] = [Q12, ln12]      # 0110
    ln012, Q012 = add_line_cln_b0_affine(Q0, Q12, PP)
    Qj[7-1] = [Q012, ln012*ln12]    # 0111
    Qj[8-1] = [Q3, 1]          # 1000
    ln03, Q03 = add_line_cln_b0_affine(Q0, Q3, PP)
    Qj[9-1] = [Q03, ln03]      # 1001
    ln13, Q13 = add_line_cln_b0_affine(Q1, Q3, PP)
    Qj[10-1] = [Q13, ln13]     # 1010
    ln013, Q013 = add_line_cln_b0_affine(Q01, Q3, PP)
    Qj[11-1] = [Q013, ln013*ln01]   # 1011
    ln23, Q23 = add_line_cln_b0_affine(Q2, Q3, PP)
    Qj[12-1] = [Q23, ln23]     # 1100
    ln023, Q023 = add_line_cln_b0_affine(Q02, Q3, PP)
    Qj[13-1] = [Q023, ln023*ln02]   # 1101
    ln123, Q123 = add_line_cln_b0_affine(Q12, Q3, PP)
    Qj[14-1] = [Q123, ln123*ln12]   # 1110
    ln0123, Q0123 = add_line_cln_b0_affine(Q012, Q3, PP)
    Qj[15-1] = [Q0123, ln0123*Qj[7-1][1]] # 1111
    return Qj

def precompute_qj_optimized_frobenius(Q0, Q1, Q2, Q3):
    # precomputations
    Qj = [None]*15     # l3 l2 l1 l0
    Qj[1-1] = Q0       # 0001
    Qj[2-1] = Q1       # 0010
    Q01 = add_cln_b0_affine(Q0, Q1)
    Qj[3-1] = Q01      # 0011
    Qj[4-1] = Q2       # 0100
    Q02 = add_cln_b0_affine(Q0, Q2)
    Qj[5-1] = Q02      # 0101
    #Q12 = add_cln_b0_affine(Q1, Q2)
    Q12 = Q0.curve()([(Q01[0]).frobenius(), (Q01[1]).frobenius()])
    Qj[6-1] = Q12      # 0110
    Q012 = add_cln_b0_affine(Q0, Q12)
    Qj[7-1] = Q012     # 0111
    Qj[8-1] = Q3       # 1000
    #Q23 = add_cln_b0_affine(Q2, Q3)
    Q23 = Q0.curve()([(Q12[0]).frobenius(), (Q12[1]).frobenius()])
    Qj[12-1] = Q23     # 1100
    #Q03 = add_cln_b0_affine(Q0, Q3)
    Q03 = Q0.curve()([(Q23[0]).frobenius(), (Q23[1]).frobenius()])
    Qj[9-1] = Q03      # 1001
    #Q13 = add_cln_b0_affine(Q1, Q3)
    Q13 = Q0.curve()([(Q02[0]).frobenius(), (Q02[1]).frobenius()])
    Qj[10-1] = Q13     # 1010
    #Q123 = add_cln_b0_affine(Q12, Q3)
    Q123 = Q0.curve()([(Q012[0]).frobenius(), (Q012[1]).frobenius()])
    Qj[14-1] = Q123    # 1110
    #Q023 = add_cln_b0_affine(Q02, Q3)
    Q023 = Q0.curve()([(Q123[0]).frobenius(), (Q123[1]).frobenius()])
    Qj[13-1] = Q023    # 1101
    #Q013 = add_cln_b0_affine(Q01, Q3)
    Q013 = Q0.curve()([(Q023[0]).frobenius(), (Q023[1]).frobenius()])
    Qj[11-1] = Q013    # 1011
    Q0123 = add_cln_b0_affine(Q012, Q3)
    Qj[15-1] = Q0123   # 1111
    return Qj

def precompute_qj_with_lines_optimized_frobenius(Q0, Q1, Q2, Q3, PP):
    """
    precomputations for 4-multi-scalar multiplication
    Cost 4 affine point additions in G2,
    7 frobenius maps on G2 (14 frobenius maps on sparse coordinates),
    7 Frobenius powers in Fp8 on 4 sparse lines and 3 dense lines,
    4 sparse-sparse mults in Fp8 and one mult in Fp8

    Works only for all-positive scalars a_i

    """
    Qj = [None]*15             # l3 l2 l1 l0
    Qj[1-1] = [Q0, 1]          # 0001
    Qj[2-1] = [Q1, 1]          # 0010
    ln01, Q01 = add_line_cln_b0_affine(Q0, Q1, PP)
    Qj[3-1] = [Q01, ln01]      # 0011
    Qj[4-1] = [Q2, 1]          # 0100
    ln02, Q02 = add_line_cln_b0_affine(Q0, Q2, PP)
    Qj[5-1] = [Q02, ln02]      # 0101
    #ln12, Q12 = add_line_cln_b0_affine(Q1, Q2, PP)
    Q12 = [(Q01[0]).frobenius(), (Q01[1]).frobenius()]; ln12 = ln01.frobenius()
    Qj[6-1] = [Q12, ln12]      # 0110
    ln012, Q012 = add_line_cln_b0_affine(Q0, Q12, PP)
    Qj[7-1] = [Q012, ln012*ln12]    # 0111
    Qj[8-1] = [Q3, 1]          # 1000
    #ln23, Q23 = add_line_cln_b0_affine(Q2, Q3, PP)
    Q23 = [(Q12[0]).frobenius(), (Q12[1]).frobenius()]; ln23 = ln12.frobenius()
    Qj[12-1] = [Q23, ln23]     # 1100
    #ln03, Q03 = add_line_cln_b0_affine(Q0, Q3, PP)
    Q03 = [(Q23[0]).frobenius(), (Q23[1]).frobenius()]; ln03 = ln23.frobenius()
    Qj[9-1] = [Q03, ln03]      # 1001
    #ln13, Q13 = add_line_cln_b0_affine(Q1, Q3, PP)
    Q13 = [(Q02[0]).frobenius(), (Q02[1]).frobenius()]; ln13 = ln02.frobenius()
    Qj[10-1] = [Q13, ln13]     # 1010
    #ln123, Q123 = add_line_cln_b0_affine(Q12, Q3, PP)
    Q123 = [(Q012[0]).frobenius(), (Q012[1]).frobenius()]; ln123 = ln012.frobenius()
    Qj[14-1] = [Q123, ln123*ln12]   # 1110
    #ln023, Q023 = add_line_cln_b0_affine(Q02, Q3, PP)
    Q023 = [(Q123[0]).frobenius(), (Q123[1]).frobenius()]; ln023 = ln123.frobenius()
    Qj[13-1] = [Q023, ln023*ln02]   # 1101
    #ln013, Q013 = add_line_cln_b0_affine(Q01, Q3, PP)
    Q013 = [(Q023[0]).frobenius(), (Q023[1]).frobenius()]; ln013 = ln023.frobenius()
    Qj[11-1] = [Q013, ln013*ln01]   # 1011
    ln0123, Q0123 = add_line_cln_b0_affine(Q012, Q3, PP)
    Qj[15-1] = [Q0123, ln0123*Qj[7-1][1]] # 1111

def precompute_qj_with_lines_a0_sextic_twist(Q0, Q1, Q2, Q3, PP, Fq6, D_twist):
    """
    Precompute all possible sums of points in {Q0, Q1, Q2, Q3} and product of lines evaluated at PP

    INPUT:
    - `Qi`: four points on the same curve
    - `PP`: point in affine coordinates on Fp

    For example Q0 = Q, Q1 = pi_p(Q), Q2 = pi_p2(Q), Q3 = pi_p3(Q)
    Compute all sums of 2 to 4 distinct points.
    There are 15 such points.
    Compute the product of the lines involved in the sum, evaluated at PP.

    Needed in miller_function_ate_aklgl_a0_4_multi_scalar
    """
    # precomputations
    # 1. precomputations for Zhang Lin formula: s
    PP = (PP[0], PP[1], PP[0]**2)
    #
    if Fq6.base_ring().degree() == 1:
        xi = -Fq6.polynomial().constant_coefficient()
    else:
        xi = -Fq6.modulus().constant_coefficient()

    Qj = [None]*15             # l3 l2 l1 l0
    Qj[1-1] = [Q0, 1]          # 0001
    Qj[2-1] = [Q1, 1]          # 0010

    ln01, Q01 = add_line_zl_a0_affine(Q0, Q1, PP, xi, True, D_twist, False)
    ln01 = Fq6(ln01)
    Q01 = (Q01[0], Q01[1])
    Qj[3-1] = [Q01, ln01]      # 0011
    Qj[4-1] = [Q2, 1]          # 0100
    ln02, Q02 = add_line_zl_a0_affine(Q0, Q2, PP, xi, True, D_twist, False)
    ln02 = Fq6(ln02)
    Q02 = (Q02[0], Q02[1])
    Qj[5-1] = [Q02, ln02]      # 0101
    ln12, Q12 = add_line_zl_a0_affine(Q1, Q2, PP, xi, True, D_twist, False)
    ln12 = Fq6(ln12)
    Q12 = (Q12[0], Q12[1])
    Qj[6-1] = [Q12, ln12]      # 0110
    ln012, Q012 = add_line_zl_a0_affine(Q0, Q12, PP, xi, True, D_twist, False)
    ln012 = Fq6(ln012)
    Q012 = (Q012[0], Q012[1])
    Qj[7-1] = [Q012, ln012*ln12]    # 0111
    Qj[8-1] = [Q3, 1]          # 1000
    ln03, Q03 = add_line_zl_a0_affine(Q0, Q3, PP, xi, True, D_twist, False)
    ln03 = Fq6(ln03)
    Q03 = (Q03[0], Q03[1])
    Qj[9-1] = [Q03, ln03]      # 1001
    ln13, Q13 = add_line_zl_a0_affine(Q1, Q3, PP, xi, True, D_twist, False)
    ln13 = Fq6(ln13)
    Q13 = (Q13[0], Q13[1])
    Qj[10-1] = [Q13, ln13]     # 1010
    ln013, Q013 = add_line_zl_a0_affine(Q01, Q3, PP, xi, True, D_twist, False)
    ln013 = Fq6(ln013)
    Q013 = (Q013[0], Q013[1])
    Qj[11-1] = [Q013, ln013*ln01]   # 1011
    ln23, Q23 = add_line_zl_a0_affine(Q2, Q3, PP, xi, True, D_twist, False)
    ln23 = Fq6(ln23)
    Q23 = (Q23[0], Q23[1])
    Qj[12-1] = [Q23, ln23]     # 1100
    ln023, Q023 = add_line_zl_a0_affine(Q02, Q3, PP, xi, True, D_twist, False)
    ln023 = Fq6(ln023)
    Q023 = (Q023[0], Q023[1])
    Qj[13-1] = [Q023, ln023*ln02]   # 1101
    ln123, Q123 = add_line_zl_a0_affine(Q12, Q3, PP, xi, True, D_twist, False)
    ln123 = Fq6(ln123)
    Q123 = (Q123[0], Q123[1])
    Qj[14-1] = [Q123, ln123*ln12]   # 1110
    ln0123, Q0123 = add_line_zl_a0_affine(Q012, Q3, PP, xi, True, D_twist, False)
    ln0123 = Fq6(ln0123)
    Q0123 = (Q0123[0], Q0123[1])
    Qj[15-1] = [Q0123, ln0123*Qj[7-1][1]] # 1111
    return Qj

def miller_function_tate_cln_b0_2_parallel(P, Q, a, a0, a1, omega=None):
    """Miller loop with Costello-Lauter-Naehrig formulas
    and two Miller functions that can be done in parallel on two cores

    INPUT:
    - `P`: point on E(Fp)
    - `Q`: point on E2(Fq)
    - `a`: E2 curve coefficient (for Q)
    - `scalars`: a0, a1, omega, scalars s.t. a0 + a1*(p**2%r) = 0 mod r
    and omega^2 + 1 = 0 mod p, so that lambda is the eigenvalue of
    psi: (x,y) -> (-x, omega*y) on G1

    f_{a0,P}(Q) * f_{a1,P}(Q)^(p^2)

    No lines needed as a0 + a1*psi(P) = O
    multiply the sparse lines or tangent two by two before accumulating
    """
    f_a0_P_Q, a0P = miller_function_tate_cln_b0(P, Q, a, a0)
    f_a1_P_Q, a1P = miller_function_tate_cln_b0(P, Q, a, a1)
    if omega is not None:
        psi_a1P = (-a1P[0], a1P[1]*omega, a1P[2])
        ln, a0P_psi_a1P = add_line_cln_b0_with_z(a0P, psi_a1P, (Q[0], Q[1])) # it should be a vertical
        return f_a0_P_Q * f_a1_P_Q.frobenius(2), a0P_psi_a1P
    return f_a0_P_Q * f_a1_P_Q.frobenius(2)

def miller_function_tate_cln_b0_2_parallel_alt(P, Q, a, a0, a1, omega):
    """Miller loop with Costello-Lauter-Naehrig formulas
    and two Miller loops that can be done in parallel on two cores

    INPUT:
    - `P`: point on E(Fp)
    - `Q`: point on E2(Fq)
    - `a`: E2 curve coefficient (for Q)
    - `scalars`: a0, a1, omega, scalars s.t. a0 + a1*(p**2%r) = 0 mod r
    and omega^2 + 1 = 0 mod p, so that lambda is the eigenvalue of
    psi: (x,y) -> (-x, omega*y) on G1

    f_{a0,P}(Q) * f_{a1,psi(P)}(Q)

    No lines needed as a0 + a1*psi(P) = O
    multiply the sparse lines or tangent two by two before accumulating
    """
    f_a0_P_Q, a0P = miller_function_tate_cln_b0(P, Q, a, a0)
    psiP = (-P[0], P[1]*omega)
    f_a1_psiP_Q, a1psiP = miller_function_tate_cln_b0(psiP, Q, a, a1)
    ln, a0P_a1psiP = add_line_cln_b0_with_z(a0P, a1psiP, (Q[0], Q[1])) # it should be a vertical
    # the point is useless but it is only for checking that it is actually 0
    return f_a0_P_Q * f_a1_psiP_Q, a0P_a1psiP

def miller_function_tate_cln_b0_2_multi_scalar(P, Q, a, a0, a1, omega):
    """Miller loop with Costello-Lauter-Naehrig formulas
    and multi-scalar technique

    INPUT:
    - `P`: point on E(Fp)
    - `Q`: point on E2(Fq)
    - `a`: E2 curve coefficient (for Q)
    - `scalars`: a0, a1, omega, scalars s.t. a0 + a1*(p**2%r) = 0 mod r
    and omega^2 + 1 = 0 mod p, so that lambda is the eigenvalue of
    psi: (x,y) -> (-x, omega*y) on G1

    f_{a0+a1*lambda,P}(Q)

    where [a0+a1*lambda]P = a0 + a1*psi(P) = O
    multiply the sparse lines or tangent two by two before accumulating
    """
    P0 = (P[0], P[1])
    P1 = (-P[0], P[1]*omega)
    assert P1 in P.curve()
    QQ = (Q[0], Q[1])

    if a0 < 0:
        a0 = -a0
        P0 = (P[0], -P[1])
    if a1 < 0:
        a1 = -a1
        P1 = (P1[0], -P1[1])

    l0 = Integer(a0).digits(2)
    l1 = Integer(a1).digits(2)
    ll0 = len(l0)
    ll1 = len(l1)
    max_l = max(ll0, ll1)
    l0 = l0 + [0] * (max_l - ll0)
    l1 = l1 + [0] * (max_l - ll1)
    # precompute P0 + P1
    ln01, P01 = add_line_cln_b0_affine(P0, P1, QQ)
    Pj = [P0, P1, P01]
    i = max_l - 1
    j = l0[i] + l1[i]*2
    assert j > 0
    S = Pj[j-1]
    S = (S[0], S[1], 1)
    if j == 3:
        m = ln01
    else:
        m = 1
    for i in range(max_l-2, -1, -1):
        # with P0, P1
        m = m**2
        tg, S = double_line_cln_b0(S, QQ, a)
        j = l0[i] + 2*l1[i]
        if j > 0:
            ln, S = add_line_cln_b0(S, Pj[j-1], QQ)
            ln = tg * ln # sparse-sparse
            m = m * ln
            if j == 3:
                m = m * ln01
        else:
            m = m * tg
    return m, S

def miller_function_ate_cln_b0_4_parallel(Q, P, a, scalars):
    """Miller loop with Costello-Lauter-Naehrig formulas
    and 4 Miller function in parallel

    INPUT:
    - `scalars`: [a0,a1,a2,a3]
    - `Q`: point on E2(Fq)
    - `P`: point on E(Fp)
    - `a`: E2 curve coefficient (for Q)

    f_{a0,Q}(P) * f_{a1,Q}(P)^p * f_{a2,Q}(P)^p^2 * f_{a3,Q}(P)^p^3
    * l_{[a0]Q, [a1]pQ}(P) * l_{[a2]p2Q, [a3]p3Q}(P))
    """
    assert len(scalars) == 4
    a0,a1,a2,a3 = scalars
    f_a0, a0Q = miller_function_ate_cln_b0(Q, P, a, a0)
    f_a1, a1Q = miller_function_ate_cln_b0(Q, P, a, a1)
    f_a2, a2Q = miller_function_ate_cln_b0(Q, P, a, a2)
    f_a3, a3Q = miller_function_ate_cln_b0(Q, P, a, a3)

    pi_a1Q = (a1Q[0].frobenius(), a1Q[1].frobenius(), a1Q[2].frobenius())
    ln01, a0Q_pi_a1Q = add_line_cln_b0_with_z(a0Q, pi_a1Q, (P[0], P[1]))
    pi2_a2Q = (a2Q[0].frobenius(2), a2Q[1].frobenius(2), a2Q[2].frobenius(2))
    pi3_a3Q = (a3Q[0].frobenius(3), a3Q[1].frobenius(3), a3Q[2].frobenius(3))
    ln23, pi2_a2Q_pi3_a3Q = add_line_cln_b0_with_z(pi2_a2Q, pi3_a3Q, (P[0], P[1]))

    f = f_a0 * f_a1.frobenius() * f_a2.frobenius(2) * f_a3.frobenius(3) * ln01 * ln23
    # finally
    ln0123, ai_pi_Q = add_line_cln_b0_with_z(a0Q_pi_a1Q, pi2_a2Q_pi3_a3Q, (P[0], P[1]))
    return f, ai_pi_Q

def miller_function_ate_cln_b0_4_parallel_alt(Q, P, a, scalars):
    """Miller loop with Costello-Lauter-Naehrig formulas
    and 4 Miller function in parallel

    INPUT:
    - `scalars`: [a0,a1,a2,a3]
    - `Q`: point on E2(Fq)
    - `P`: point on E(Fp)
    - `a`: E2 curve coefficient (for Q)

    f_{a0,Q}(P) * f_{a1,[p]Q}(P) * f_{a2,[p^2]Q}(P) * f_{a3,[p^3]Q}(P)
    * l_{[a0]Q, [a1]pQ}(P) * l_{[a2]p2Q, [a3]p3Q}(P))
    """
    assert len(scalars) == 4
    a0,a1,a2,a3 = scalars
    Q1 = (Q[0].frobenius(), Q[1].frobenius())
    Q2 = (Q[0].frobenius(2), Q[1].frobenius(2))
    Q3 = (Q[0].frobenius(3), Q[1].frobenius(3))
    f_a0, a0Q = miller_function_ate_cln_b0(Q, P, a, a0)
    f_a1, a1piQ = miller_function_ate_cln_b0(Q1, P, a, a1)
    f_a2, a2pi2Q = miller_function_ate_cln_b0(Q2, P, a, a2)
    f_a3, a3pi3Q = miller_function_ate_cln_b0(Q3, P, a, a3)

    ln01, a0Q_a1piQ = add_line_cln_b0_with_z(a0Q, a1piQ, (P[0], P[1]))
    ln23, a2pi2Q_a3pi3Q = add_line_cln_b0_with_z(a2pi2Q, a3pi3Q, (P[0], P[1]))

    f = f_a0 * f_a1 * f_a2 * f_a3 * ln01 * ln23
    # finally
    ln0123, ai_pi_Q = add_line_cln_b0_with_z(a0Q_a1piQ, a2pi2Q_a3pi3Q, (P[0], P[1]))
    return f, ai_pi_Q

def miller_function_ate_cln_b0_4_multi_scalar(Q, P, a, scalars):
    """Miller loop with Costello-Lauter-Naehrig formulas
    and multi-scalar technique

    INPUT:
    - `scalars`: [a0,a1,a2,a3]
    - `Q`: point on E2(Fq)
    - `P`: point on E(Fp)
    - `a`: E2 curve coefficient (for Q)

    f_{a0 + a1*p + a2*p^2 + a3*p^3,Q}(P)
    accumulate the lines or tangent two by two before multiplying
    """
    assert len(scalars) == 4
    a0,a1,a2,a3 = scalars
    all_positive = a0 >= 0 and a1 >= 0 and a2 >= 0 and a3 >= 0
    all_negative = a0 < 0 and a1 < 0 and a2 < 0 and a3 < 0

    PP = (P[0], P[1])
    Q0 = (Q[0], Q[1])
    if a0 < 0:
        a0 = -a0
        Q0 = (Q0[0], -Q0[1])
    Q1 = (Q[0].frobenius(), Q[1].frobenius())
    assert Q1 in Q.curve()
    if a1 < 0:
        a1 = -a1
        Q1 = (Q1[0], -Q1[1])
    Q2 = (Q[0].frobenius(2), Q[1].frobenius(2))
    assert Q2 in Q.curve()
    if a2 < 0:
        a2 = -a2
        Q2 = (Q2[0], -Q2[1])
    Q3 = (Q[0].frobenius(3), Q[1].frobenius(3))
    assert Q3 in Q.curve()
    if a3 < 0:
        a3 = -a3
        Q3 = (Q3[0], -Q3[1])
    if all_positive or all_negative:
        Qj = precompute_qj_with_lines_optimized_frobenius(Q0, Q1, Q2, Q3, PP)
    else:
        Qj = precompute_qj_with_lines(Q0, Q1, Q2, Q3, PP)
    m = 1
    l0 = Integer(a0).digits(2)
    l1 = Integer(a1).digits(2)
    l2 = Integer(a2).digits(2)
    l3 = Integer(a3).digits(2)
    ll0 = len(l0)
    ll1 = len(l1)
    ll2 = len(l2)
    ll3 = len(l3)
    max_l = max([ll0, ll1, ll2, ll3])
    l0 = l0 + [0]*(max_l - len(l0))
    l1 = l1 + [0]*(max_l - len(l1))
    l2 = l2 + [0]*(max_l - len(l2))
    l3 = l3 + [0]*(max_l - len(l3))
    # init
    i = max_l - 1
    j = l0[i] + l1[i]*2 + l2[i]*4 + l3[i]*8
    assert j > 0
    S = Qj[j-1][0]
    S = (S[0], S[1], 1)
    m = Qj[j-1][1]
    for i in range(max_l - 2, -1, -1):
        m = m**2
        tg, S = double_line_cln_b0(S, PP, a)
        j = l0[i] + l1[i]*2 + l2[i]*4 + l3[i]*8
        if j > 0:
            ln, S = add_line_cln_b0(S, Qj[j-1][0], PP)
            tg = tg * ln # sparse-sparse
            m = m * tg
            ln = Qj[j-1][1]
            if ln != 1:
                m = m * ln
        else:
            m = m * tg
    return m, S

def miller_function_tate_aklgl_a0_2_parallel(P, Q, b, a0, a1, Fq6, map_Fq6_Fpk, D_twist=False, omega=None, xi=None):
    """Miller loop with AKLGL formulas
    and two Miller functions that can be done in parallel on two cores

    INPUT:
    - `P`: point on E(Fp)
    - `Q`: point on E2(Fq)
    - `b`: E2 curve coefficient (for Q)
    - `scalars`: a0, a1, omega, scalars s.t. a0 + a1*(p**(k//3) % r) = 0 mod r
    and omega^2 + omega + 1 = 0 mod p, so that lambda=(p**(k//3) % r) is the eigenvalue of
    psi: (x,y) -> (omega*x, y) on G1

    f_{a0,P}(Q) * f_{a1,P}(Q)^(p^(k//3))

    No lines needed as a0 + a1*psi(P) = O
    multiply the sparse lines or tangent two by two before accumulating
    """
    e = Fq6.base_ring().degree()
    f_a0_P_Q, a0P = miller_function_tate_aklgl(P, Q, b, a0, Fq6, D_twist=D_twist, xi=xi)
    f_a1_P_Q, a1P = miller_function_tate_aklgl(P, Q, b, a1, Fq6, D_twist=D_twist, xi=xi)
    if omega is not None: # apply endomorphism
        psi_a1P = (omega*a1P[0], a1P[1], a1P[2])

        ln, a0P_psi_a1P = add_line_h_a0_twist6_aklgl_with_z(a0P, psi_a1P, (Q[0], Q[1]), D_twist=not D_twist) # it should be a vertical

        return map_Fq6_Fpk(f_a0_P_Q) * map_Fq6_Fpk(f_a1_P_Q).frobenius(e), a0P_psi_a1P
    return map_Fq6_Fpk(f_a0_P_Q) * map_Fq6_Fpk(f_a1_P_Q).frobenius(e), (a0P, a1P)

def miller_function_tate_a0_twist6_aklgl_2_multi_scalar(P, Q, b, a0, a1, omega, xi, Fq6, D_twist=False):
    """Miller loop with AKLGL formulas
    and multi-scalar technique

    INPUT:
    - `P`: point on E(Fp)
    - `Q`: point on E2(Fq) (twist curve)
    - `b`: E2 curve coefficient (for Q)
    - `scalars`: a0, a1, omega, scalars s.t. a0 + a1*(p**2%r) = 0 mod r
    and omega^2 + omega + 1 = 0 mod p, so that lambda = p^4 mod r = -p^2-1 mod r is the eigenvalue of
    psi: (x,y) -> (omega*x, y) on G1

    f_{a0+a1*lambda,P}(Q)

    where [a0+a1*lambda]P = a0 + a1*psi(P) = O
    multiply the sparse lines or tangent two by two before accumulating
    """
    P0 = (P[0], P[1])
    P1 = (omega*P[0], P[1])
    assert P1 in P.curve()
    QQ = (Q[0], Q[1])
    if a0 < 0:
        a0 = -a0
        P0 = (P[0], -P[1])
    if a1 < 0:
        a1 = -a1
        P1 = (P1[0], -P1[1])

    l0 = Integer(a0).digits(2)
    l1 = Integer(a1).digits(2)
    ll0 = len(l0)
    ll1 = len(l1)
    max_l = max(ll0, ll1)
    l0 = l0 + [0] * (max_l - ll0)
    l1 = l1 + [0] * (max_l - ll1)
    # precompute P0 + P1
    ln01, P01 = add_line_h_a0_twist6_aklgl([P0[0], P0[1], 1], P1, QQ, D_twist=not D_twist)
    P01 = [P01[0]/P01[2], P01[1]/P01[2]]
    if D_twist:
        assert ln01[1] == 0
        c = 2
    else:
        assert ln01[2] == 0
        c = 1
    la0, lac, la3 = ln01[0], ln01[c], ln01[3]
    Pj = [P0, P1, P01]
    i = max_l - 1
    j = l0[i] + l1[i]*2
    assert j > 0
    S = Pj[j-1]
    S = (S[0], S[1], 1)
    if j == 3:
        m = Fq6(ln01)
    else:
        m = Fq6(1)
    for i in range(max_l-2, -1, -1):
        # with P0, P1
        m = m**2
        tg, S = double_line_h_a0_twist6_aklgl(S, QQ, b, D_twist=not D_twist)
        j = l0[i] + 2*l1[i]
        if j > 0:
            ln, S = add_line_h_a0_twist6_aklgl(S, Pj[j-1], QQ, D_twist=not D_twist)
            #ln = tg * ln # sparse-sparse
            if D_twist:
                ln = sparse_sparse_mult_m6_twist(ln[0], ln[2], ln[3], tg[0], tg[2], tg[3], xi, Fq6)
            else:
                ln = sparse_sparse_mult_d6_twist(ln[0], ln[1], ln[3], tg[0], tg[1], tg[3], xi, Fq6)
            m = m * ln
            if j == 3:
                #m = m * ln01
                if D_twist:
                    m = sparse_mult_m6_twist(ln01[0],ln01[2],ln01[3], m, xi, Fq6)
                else:
                    m = sparse_mult_d6_twist(ln01[0],ln01[1],ln01[3], m, xi, Fq6)
        else:
            #m = m * tg
            if D_twist:
                m = sparse_mult_m6_twist(tg[0], tg[2], tg[3], m, xi, Fq6)
            else:
                m = sparse_mult_d6_twist(tg[0], tg[1], tg[3], m, xi, Fq6)
    return m, S


def miller_function_ate_a0_twist6_aklgl_4_multi_scalar(Q, P, b, scalars, Fq6, xi_p_1_3, xi_p_1_2, xi_p2_1_3, xi_p2_1_2, D_twist=False):
    """Miller loop with explicit twist of degree 6 and AKLGL formulas,
    and multi-scalar technique

    INPUT:
    - `Q`: point on E2(Fq)
    - `P`: point on E(Fp)
    - `b`: E2 curve coefficient (for Q)
    - `scalars`: [a0,a1,a2,a3]
    - `Fq6`: degree 6 extension of Fq, of generator w, s.t. w^6 = xi
    - `xi_p_1_3`: precomputed value of xi^((p-1)/3)
    - `xi_p_1_2`: precomputed value of xi^((p-1)/2)

    f_{a0 + a1*p + a2*p^2 + a3*p^3,Q}(P)
    accumulate the lines or tangent two by two before multiplying
    """
    if Fq6.base_ring().degree() == 1:
        xi = -Fq6.polynomial().constant_coefficient()
    else:
        xi = -Fq6.modulus().constant_coefficient()

    assert len(scalars) == 4
    a0,a1,a2,a3 = scalars
    all_positive = a0 >= 0 and a1 >= 0 and a2 >= 0 and a3 >= 0
    all_negative = a0 < 0 and a1 < 0 and a2 < 0 and a3 < 0

    PP = (P[0], P[1])
    Q0 = (Q[0], Q[1])
    if a0 < 0:
        a0 = -a0
        Q0 = (Q0[0], -Q0[1])
    # apply the endomorphism:
    # w^p = w^(p-1)*w = xi^((p-1)/6)*w
    # w^(2p) = xi^((p-1)/3)*w^2
    # w^(3p) = xi^((p-1)/2)*w^3
    if D_twist: #(x,y) -> (x*w^2, y*w^3) -> (x^p*xi^((p-1)/3)*w^2, y^p*xi^((p-1)/2)*w^3)
        Q1 = (Q[0].frobenius() * xi_p_1_3, Q[1].frobenius() * xi_p_1_2)
    else:
        Q1 = (Q[0].frobenius() / xi_p_1_3, Q[1].frobenius() / xi_p_1_2)
    assert Q1 in Q.curve()
    if a1 < 0:
        a1 = -a1
        Q1 = (Q1[0], -Q1[1])
    if D_twist: #(x,y) -> (x*w^2, y*w^3) -> (x^(p^2)*xi^((p+1)*(p-1)/3)*w^2, y^(p^2)*xi^((p+1)*(p-1)/2)*w^3)
        Q2 = (Q[0].frobenius(2) * xi_p2_1_3, Q[1].frobenius(2) * xi_p2_1_2)
    else:
        Q2 = (Q[0].frobenius(2) / xi_p2_1_3, Q[1].frobenius(2) / xi_p2_1_2)
    assert Q2 in Q.curve()
    if a2 < 0:
        a2 = -a2
        Q2 = (Q2[0], -Q2[1])
    if D_twist: #(x,y) -> (x*w^2, y*w^3) -> (x^(p^3)*xi^((p^3-1)/3)*w^2, y^(p^3)*xi^((p^3-1)/2)*w^3)
        # and xi^((p^3-1)/3) = xi^((p^2+p+1)*(p-1)/3) = xi^((p-1)/3)^p^2 * xi*((p^2-1)/3)
        Q3 = (Q[0].frobenius(3) * xi_p_1_3 * xi_p2_1_3, Q[1].frobenius(3) * xi_p_1_2 * xi_p2_1_2)
    else:
        Q3 = (Q[0].frobenius(3) / (xi_p_1_3 * xi_p2_1_3), Q[1].frobenius(3) / (xi_p_1_2 * xi_p2_1_2))
    assert Q3 in Q.curve()
    if a3 < 0:
        a3 = -a3
        Q3 = (Q3[0], -Q3[1])
    Qj = precompute_qj_with_lines_a0_sextic_twist(Q0, Q1, Q2, Q3, PP, Fq6, D_twist)
    m = 1
    l0 = Integer(a0).digits(2)
    l1 = Integer(a1).digits(2)
    l2 = Integer(a2).digits(2)
    l3 = Integer(a3).digits(2)
    ll0 = len(l0)
    ll1 = len(l1)
    ll2 = len(l2)
    ll3 = len(l3)
    max_l = max([ll0, ll1, ll2, ll3])
    l0 = l0 + [0]*(max_l - len(l0))
    l1 = l1 + [0]*(max_l - len(l1))
    l2 = l2 + [0]*(max_l - len(l2))
    l3 = l3 + [0]*(max_l - len(l3))
    # init
    i = max_l - 1
    j = l0[i] + l1[i]*2 + l2[i]*4 + l3[i]*8
    assert j > 0
    S = Qj[j-1][0]
    S = (S[0], S[1], 1)
    m = Qj[j-1][1]
    for i in range(max_l - 2, -1, -1):
        m = m**2
        tg, S = double_line_h_a0_twist6_aklgl(S, PP, b, D_twist)
        j = l0[i] + l1[i]*2 + l2[i]*4 + l3[i]*8
        if j > 0:
            ln, S = add_line_h_a0_twist6_aklgl(S, Qj[j-1][0], PP, D_twist)
            if D_twist:
                tg = sparse_sparse_mult_d6_twist(tg[0], tg[1], tg[3], ln[0], ln[1], ln[3], xi, Fq6)
            else:
                tg = sparse_sparse_mult_m6_twist(tg[0], tg[2], tg[3], ln[0], ln[2], ln[3], xi, Fq6)
            m = m * tg
            ln = Qj[j-1][1]
            if ln != 1:
                # sparse or non-sparse? for now, ln is already in Fq6
                m = m * ln
        else:
            if D_twist:
                m = sparse_mult_d6_twist(tg[0], tg[1], tg[3], m, xi, Fq6)
            else:
                m = sparse_mult_m6_twist(tg[0], tg[2], tg[3], m, xi, Fq6)
    return m, S
