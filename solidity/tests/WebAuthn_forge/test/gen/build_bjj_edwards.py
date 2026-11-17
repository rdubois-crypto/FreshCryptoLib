"""usage : python3 build_bjj_edwards.py <FCL_babyjubjub_edwards.sol>
Generates FCL_babyjubjub_edwards.sol: BabyJubJub EdDSA verification in extended Edwards coordinates on the
a=-1 model (complete formulas add-2008-hwcd-3 / dbl-2008-hwcd), joint 2+2 bit window, batch normalized table,
lazy reduction. Baseline against the Weierstrass isomorphism version."""
import re
import sys

import bjj

r, l, AE, DE = bjj.r, bjj.l, bjj.AE, bjj.DE
SQ = bjj.sqrt((-AE) % r)           # s^2 = -a: (x, y) -> (s.x, y) maps a.x^2+y^2=1+d.x^2y^2 to -x^2+y^2=1+d'.x^2y^2
D1 = (-DE) * pow(AE, -1, r) % r    # d' = -d/a (non square: complete formulas)
K2D = 2 * D1 % r
iso = lambda P: (SQ * P[0] % r, P[1])
B8 = iso(bjj.B8)
B2 = iso(bjj.ed_mul(2, bjj.B8))
B3 = iso(bjj.ed_mul(3, bjj.B8))
H = lambda v: "0x%064x" % v
DMAP = bytes([4, 8, 12, 5, 6, 7, 9, 10, 11, 13, 14, 15]).hex().ljust(64, "0")

# projective doubling without T: (X,Y,Z) -> 2P
DBL = """{
    let A := mulmod(X, X, p)
    let B := mulmod(Y, Y, p)
    let E := add(mulmod(add(X, Y), add(X, Y), p), sub(_2P, add(A, B)))
    let G := add(B, sub(p, A))
    let F := add(G, sub(p, mulmod(add(Z, Z), Z, p)))
    X := mulmod(E, F, p)
    Y := mulmod(G, sub(_2P, add(A, B)), p)
    Z := mulmod(F, G, p)
}"""
# doubling with T output (before an addition)
DBLT = """{
    let A := mulmod(X, X, p)
    let B := mulmod(Y, Y, p)
    let E := add(mulmod(add(X, Y), add(X, Y), p), sub(_2P, add(A, B)))
    let G := add(B, sub(p, A))
    let F := add(G, sub(p, mulmod(add(Z, Z), Z, p)))
    B := sub(_2P, add(A, B))
    X := mulmod(E, F, p)
    Y := mulmod(G, B, p)
    T := mulmod(E, B, p)
    Z := mulmod(F, G, p)
}"""
# complete mixed addition, table entry [y-x, y+x, 2d'.x.y], no T output
MADD = """{
    let A := mulmod(add(Y, sub(p, X)), mload(e), p)
    let B := mulmod(add(Y, X), mload(add(e, 0x20)), p)
    let C := mulmod(T, mload(add(e, 0x40)), p)
    let D := add(Z, Z)
    let E := add(B, sub(p, A))
    B := add(B, A)
    A := add(D, sub(p, C))
    D := add(D, C)
    X := mulmod(E, A, p)
    Y := mulmod(D, B, p)
    Z := mulmod(A, D, p)
}"""

src = f'''//********************************************************************************************/
//  ___           _       ___               _         _    _ _
// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__
// | _| '_/ -_|_-< ' \\  | (__| '_| || | '_ \\  _/ _ \\ | |__| | '_ \\
// |_||_| \\___/__/_||_|  \\___|_|  \\_, | .__/\\__\\___/ |____|_|_.__/
//                                |__/|_|
///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project
///* License: This software is licensed under MIT License
///* This Code may be reused including license and copyright notice.
///* See LICENSE file at the root folder of the project.
///* FILE: FCL_babyjubjub_edwards.sol
///*
///*
///* DESCRIPTION: BabyJubJub (circomlib) EdDSA verification in extended twisted Edwards coordinates,
///*  on the isomorphic a=-1 model (x,y)->(sqrt(-168700).x, y), d'=-168696/168700, as ed25519.
///*  d' is a non square mod r: the formulas of Hisil, Wong, Carter, Dawson 2008 are complete,
///*  no exceptional case, no point at infinity to handle. Mixed addition with [y-x, y+x, 2d'xy]: 6 mulmod.
///*  Joint 2-bit window table (iB8+jQ, Q=8A) normalized with a single inversion, lazy reduction.
///*
//**************************************************************************************/
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

library FCL_BabyJubJub_Edwards {{
    uint256 constant p = {r};
    uint256 constant n = {l};
    uint256 constant ED_A = {AE};
    uint256 constant ED_D = {DE};
    // isomorphic a=-1 model: (x, y) -> (SQ.x, y), -x^2+y^2 = 1 + D1.x^2y^2
    uint256 constant SQ = {H(SQ)};
    uint256 constant D1 = {H(D1)};
    uint256 constant K2D = {H(K2D)};
    uint256 constant minus_2 = {H(r - 2)};
    uint256 constant _2P = {H(2 * r)};
    // image of Base8 in the a=-1 model
    uint256 constant B8X = {H(B8[0])};
    uint256 constant B8Y = {H(B8[1])};
    // table index (i+4j) of the 12 points computed in projective form: Q,2Q,3Q,Q+B8,..,3Q+3B8
    uint256 constant _DMAP = 0x{DMAP};

    function ed_isOnCurve(uint256 x, uint256 y) internal pure returns (bool) {{
        if (x >= p || y >= p) return false;
        uint256 x2 = mulmod(x, x, p);
        uint256 y2 = mulmod(y, y, p);
        return addmod(mulmod(ED_A, x2, p), y2, p) == addmod(1, mulmod(ED_D, mulmod(x2, y2, p), p), p);
    }}

    /**
     * @dev u.B8 + v.(8A) in projective coordinates (X:Y:Z), A on the curve (affine), any 256-bit scalars
     */
    function ed_mulmuladd8(uint256 Ax, uint256 Ay, uint256 scalar_u, uint256 scalar_v)
        internal
        view
        returns (uint256 X, uint256 Y, uint256 Z)
    {{
        assembly ("memory-safe") {{
            //extended (x,y,z,t) addition of two points stored in memory, a=-1, add-2008-hwcd-3, complete
            function addEmem(p1, p2) -> x3, y3, z3, t3 {{
                let x1 := mload(p1)
                let y1 := mload(add(p1, 0x20))
                x3 := mulmod(addmod(y1, sub(p, x1), p), addmod(mload(add(p2, 0x20)), sub(p, mload(p2)), p), p) //A
                y3 := mulmod(add(y1, x1), add(mload(add(p2, 0x20)), mload(p2)), p) //B
                t3 := mulmod(mulmod(mload(add(p1, 0x60)), mload(add(p2, 0x60)), p), K2D, p) //C
                z3 := mulmod(add(mload(add(p1, 0x40)), mload(add(p1, 0x40))), mload(add(p2, 0x40)), p) //D
                let E := addmod(y3, sub(p, x3), p)
                let Hh := addmod(y3, x3, p)
                let F := addmod(z3, sub(p, t3), p)
                let G := addmod(z3, t3, p)
                x3 := mulmod(E, F, p)
                y3 := mulmod(G, Hh, p)
                t3 := mulmod(E, Hh, p)
                z3 := mulmod(F, G, p)
            }}
            function dblE(x1, y1, z1) -> x3, y3, z3, t3 {{
                let A := mulmod(x1, x1, p)
                let B := mulmod(y1, y1, p)
                let E := addmod(mulmod(add(x1, y1), add(x1, y1), p), sub(_2P, add(A, B)), p)
                let G := addmod(B, sub(p, A), p)
                let F := addmod(G, sub(p, mulmod(add(z1, z1), z1, p)), p)
                B := sub(_2P, add(A, B))
                x3 := mulmod(E, F, p)
                y3 := mulmod(G, B, p)
                t3 := mulmod(E, B, p)
                z3 := mulmod(F, G, p)
            }}
            function st4(ptr, x, y) {{
                mstore(ptr, addmod(y, sub(p, x), p))
                mstore(add(ptr, 0x20), addmod(y, x, p))
                mstore(add(ptr, 0x40), mulmod(mulmod(x, y, p), K2D, p))
            }}
            function stP(ptr, x, y, z, t) {{
                mstore(ptr, x)
                mstore(add(ptr, 0x20), y)
                mstore(add(ptr, 0x40), z)
                mstore(add(ptr, 0x60), t)
            }}
            function addSt(src, dst, g) {{
                let a1, b1, c1, d1 := addEmem(src, g)
                stP(dst, a1, b1, c1, d1)
            }}
            function buildQ(T, qx, qy) {{
                // 12 points in extended coordinates: Q, 2Q, 3Q, Q+iB8, 2Q+iB8, 3Q+iB8 ; Q = 8A
                let S0 := add(T, 0x800)
                {{
                    let ta, tb, tc, td := dblE(qx, qy, 1)
                    ta, tb, tc, td := dblE(ta, tb, tc)
                    ta, tb, tc, td := dblE(ta, tb, tc)
                    stP(S0, ta, tb, tc, td)
                    let xa, xb, xc, xd := dblE(ta, tb, tc)
                    stP(add(S0, 0x80), xa, xb, xc, xd)
                }}
                addSt(S0, add(S0, 0x100), add(S0, 0x80))
                // iB8 in extended coordinates
                let GS := add(T, 0x1000)
                stP(GS, B8X, B8Y, 1, mulmod(B8X, B8Y, p))
                stP(add(GS, 0x80), {H(B2[0])}, {H(B2[1])}, 1, {H(B2[0]*B2[1] % r)})
                stP(add(GS, 0x100), {H(B3[0])}, {H(B3[1])}, 1, {H(B3[0]*B3[1] % r)})
                addSt(S0, add(S0, 0x180), GS)
                addSt(S0, add(S0, 0x200), add(GS, 0x80))
                addSt(S0, add(S0, 0x280), add(GS, 0x100))
                addSt(add(S0, 0x80), add(S0, 0x300), GS)
                addSt(add(S0, 0x80), add(S0, 0x380), add(GS, 0x80))
                addSt(add(S0, 0x80), add(S0, 0x400), add(GS, 0x100))
                addSt(add(S0, 0x100), add(S0, 0x480), GS)
                addSt(add(S0, 0x100), add(S0, 0x500), add(GS, 0x80))
                addSt(add(S0, 0x100), add(S0, 0x580), add(GS, 0x100))
                // Montgomery batch inversion of the 12 Z (never 0, complete formulas)
                let PF := add(S0, 0x600)
                let ta := 1
                for {{ let k := 0 }} lt(k, 0x600) {{ k := add(k, 0x80) }} {{
                    ta := mulmod(ta, mload(add(S0, add(k, 0x40))), p)
                    mstore(add(PF, shr(2, k)), ta)
                }}
                {{
                    let M0 := add(PF, 0x180)
                    mstore(M0, 0x20)
                    mstore(add(M0, 0x20), 0x20)
                    mstore(add(M0, 0x40), 0x20)
                    mstore(add(M0, 0x60), ta)
                    mstore(add(M0, 0x80), minus_2)
                    mstore(add(M0, 0xa0), p)
                    if iszero(staticcall(gas(), 0x05, M0, 0xc0, M0, 0x20)) {{ revert(0, 0) }}
                    ta := mload(M0)
                }}
                for {{ let k := 0x580 }} lt(k, 0x600) {{ k := sub(k, 0x80) }} {{
                    let pt := add(S0, k)
                    let tb := ta
                    if k {{ tb := mulmod(ta, mload(add(PF, sub(shr(2, k), 0x20))), p) }}
                    ta := mulmod(ta, mload(add(pt, 0x40)), p)
                    st4(add(T, shl(7, byte(shr(7, k), _DMAP))), mulmod(mload(pt), tb, p), mulmod(mload(add(pt, 0x20)), tb, p))
                }}
            }}

            let T0 := mload(0x40)
            st4(T0, 0, 1)
            st4(add(T0, 0x80), B8X, B8Y)
            st4(add(T0, 0x100), {H(B2[0])}, {H(B2[1])})
            st4(add(T0, 0x180), {H(B3[0])}, {H(B3[1])})
            buildQ(T0, mulmod(SQ, Ax, p), Ay)

            X := 0
            Y := 1
            Z := 1
            let T := 0
            Ax := 252
            if iszero(or(scalar_u, scalar_v)) {{ Ax := not(0) }} // u = v = 0: identity
            for {{}} and(lt(Ax, 256), iszero(or(shr(Ax, scalar_u), shr(Ax, scalar_v)))) {{ Ax := sub(Ax, 4) }} {{}}
            Ay := T0
            for {{}} lt(Ax, 256) {{ Ax := sub(Ax, 4) }} {{
                //@@DBL@@
                //@@DBLT@@
                {{
                    let e := or(and(shl(5, shr(Ax, scalar_u)), 0x180), and(shl(7, shr(Ax, scalar_v)), 0x600))
                    if e {{
                        e := add(Ay, e)
                        //@@MADD@@
                    }}
                }}
                //@@DBL@@
                //@@DBLT@@
                {{
                    let e := or(and(shl(7, shr(Ax, scalar_u)), 0x180), and(shl(9, shr(Ax, scalar_v)), 0x600))
                    if e {{
                        e := add(Ay, e)
                        //@@MADD@@
                    }}
                }}
            }}
        }}
    }}

    /**
     * @dev Raw circomlib EdDSA check S.B8 == R8 + 8.hm.A (hm already hashed), circomlib coordinates,
     *      same input checks as circomlibjs (A and R8 on curve, S < subOrder). No inversion besides the table one.
     */
    function verifyRaw(uint256 hm, uint256 R8x, uint256 R8y, uint256 S, uint256 Ax, uint256 Ay)
        internal
        view
        returns (bool)
    {{
        if (S >= n || !ed_isOnCurve(R8x, R8y) || !ed_isOnCurve(Ax, Ay)) return false;
        // result (X:Y:Z) stored in (Ax, Ay, S)
        (Ax, Ay, S) = ed_mulmuladd8(Ax, Ay, S, addmod(n - hm % n, 0, n));
        return Ax == mulmod(mulmod(SQ, R8x, p), S, p) && Ay == mulmod(R8y, S, p);
    }}
}} //EOF
'''
macros = {"DBL": DBL, "DBLT": DBLT, "MADD": MADD}
out = []
for line in src.split("\n"):
    m = re.match(r"^(\s*)//@@(\w+)@@\s*$", line)
    if m:
        out += [m.group(1) + x for x in macros[m.group(2)].split("\n")]
    else:
        out.append(line)
open(sys.argv[1], "w").write("\n".join(out))
print("ok")
