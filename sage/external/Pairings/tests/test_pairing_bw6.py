from pairing import final_exp_easy_k6

def test_bw6_g1_mult_by_cofactor(function_name, E, r, c, u, omega, ht, hy):
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        rP = r*P
        assert rP != E0 and c*rP == E0 and c*P != E0
        cP = function_name(P, omega, u, ht, hy)
        crP = function_name(rP, omega, u, ht, hy)
        ok = cP != E0 and r*cP == E0 and crP == E0
        i = i+1
    print("test {}: {}".format(function_name.__name__, ok))
    return ok

def test_bw6_g1_mult_by_r(function_name, E, r, c, u, omega):
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        cP = c*P
        assert cP != E0 and r*cP == E0 and r*P != E0
        rP = function_name(P, omega, u)
        rcP = function_name(cP, omega, u)
        ok = rP != E0 and c*rP == E0 and rcP == E0
        i = i+1
    print("test {}: {}".format(function_name.__name__, ok))
    return ok

def test_bw6_g1_r_subgroup_membersip_testing(function_name, E, r, c, u, omega):
    E0 = E(0)
    ok = True
    i = 0
    while ok and i < 10:
        P = E.random_element()
        cP = c*P
        assert cP != E0 and r*cP == E0 and r*P != E0
        ok = not function_name(P, omega, u) and function_name(cP, omega, u)
        i = i+1
    print("test {}: {}".format(function_name.__name__, ok))
    return ok

def test_final_exp_hard_bw6(function_name, Fpk, u, ht, hy, r, t, expected_exp):
    """
    Exponentiation to some multiple of Phi_k(p(u))/r(u) whose value is given in expected_exp
    """
    i = 0
    ok = True
    p = Fpk.characteristic()
    Fpk_1 = Fpk(1)

    while ok and i < 10:
        f = Fpk.random_element()
        m = final_exp_easy_k6(f)
        g = function_name(m, u, ht, hy)
        h = m**expected_exp
        ok1 = g == h
        ok2 = g**r == Fpk_1
        ok = ok1 and ok2
        i = i+1
    print("test {}: {} ({} {})".format(function_name.__name__, ok, ok1, ok2))
    return ok
