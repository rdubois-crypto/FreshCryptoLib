//********************************************************************************************/
//  ___           _       ___               _         _    _ _
// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__
// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
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

library FCL_BabyJubJub_Edwards {
    uint256 constant p = 21888242871839275222246405745257275088548364400416034343698204186575808495617;
    uint256 constant n = 2736030358979909402780800718157159386076813972158567259200215660948447373041;
    uint256 constant ED_A = 168700;
    uint256 constant ED_D = 168696;
    // isomorphic a=-1 model: (x, y) -> (SQ.x, y), -x^2+y^2 = 1 + D1.x^2y^2
    uint256 constant SQ = 0x0e0ff35005d6e54b2c7ac0ce84e6117774a5f6ca321ea1176c7be2c1d8b2d768;
    uint256 constant D1 = 0x1aee90f15f2189693df072d799fd11fc039b2959ebb7c867d075ca8cf4d7eb8e;
    uint256 constant K2D = 0x0578d36fdd1172a8c3909ff8b278cb9adf026a6b5db6203e5d099f85f9afd71b;
    uint256 constant minus_2 = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593efffffff;
    uint256 constant _2P = 0x60c89ce5c263405370a08b6d0302b0ba5067d090f372e12287c3eb27e0000002;
    // image of Base8 in the a=-1 model
    uint256 constant B8X = 0x1b024eef745002f42e018e11bfe76f10f0b773ad92c6ce1951e8dc38f1aac062;
    uint256 constant B8Y = 0x25797203f7a0b24925572e1cd16bf9edfce0051fb9e133774b3c257a872d7d8b;
    // table index (i+4j) of the 12 points computed in projective form: Q,2Q,3Q,Q+B8,..,3Q+3B8
    uint256 constant _DMAP = 0x04080c050607090a0b0d0e0f0000000000000000000000000000000000000000;

    function ed_isOnCurve(uint256 x, uint256 y) internal pure returns (bool) {
        if (x >= p || y >= p) return false;
        uint256 x2 = mulmod(x, x, p);
        uint256 y2 = mulmod(y, y, p);
        return addmod(mulmod(ED_A, x2, p), y2, p) == addmod(1, mulmod(ED_D, mulmod(x2, y2, p), p), p);
    }

    /**
     * @dev u.B8 + v.(8A) in projective coordinates (X:Y:Z), A on the curve (affine), any 256-bit scalars
     */
    function ed_mulmuladd8(uint256 Ax, uint256 Ay, uint256 scalar_u, uint256 scalar_v)
        internal
        view
        returns (uint256 X, uint256 Y, uint256 Z)
    {
        assembly ("memory-safe") {
            //extended (x,y,z,t) addition of two points stored in memory, a=-1, add-2008-hwcd-3, complete
            function addEmem(p1, p2) -> x3, y3, z3, t3 {
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
            }
            function dblE(x1, y1, z1) -> x3, y3, z3, t3 {
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
            }
            function st4(ptr, x, y) {
                mstore(ptr, addmod(y, sub(p, x), p))
                mstore(add(ptr, 0x20), addmod(y, x, p))
                mstore(add(ptr, 0x40), mulmod(mulmod(x, y, p), K2D, p))
            }
            function stP(ptr, x, y, z, t) {
                mstore(ptr, x)
                mstore(add(ptr, 0x20), y)
                mstore(add(ptr, 0x40), z)
                mstore(add(ptr, 0x60), t)
            }
            function addSt(src, dst, g) {
                let a1, b1, c1, d1 := addEmem(src, g)
                stP(dst, a1, b1, c1, d1)
            }
            function buildQ(T, qx, qy) {
                // 12 points in extended coordinates: Q, 2Q, 3Q, Q+iB8, 2Q+iB8, 3Q+iB8 ; Q = 8A
                let S0 := add(T, 0x800)
                {
                    let ta, tb, tc, td := dblE(qx, qy, 1)
                    ta, tb, tc, td := dblE(ta, tb, tc)
                    ta, tb, tc, td := dblE(ta, tb, tc)
                    stP(S0, ta, tb, tc, td)
                    let xa, xb, xc, xd := dblE(ta, tb, tc)
                    stP(add(S0, 0x80), xa, xb, xc, xd)
                }
                addSt(S0, add(S0, 0x100), add(S0, 0x80))
                // iB8 in extended coordinates
                let GS := add(T, 0x1000)
                stP(GS, B8X, B8Y, 1, mulmod(B8X, B8Y, p))
                stP(add(GS, 0x80), 0x22b9a68c2ed57742dbc273e1f2e4924057c3d286a1f7b7bcead33e9ed2acb925, 0x01666cafbf0a30da8b9ebeaf848a1da067a892296f1043188e1705402b6d6853, 1, 0x15c43c85a1d12f2b14066926b827babe148fc621faa7ff01a81645b70dd6fd21)
                stP(add(GS, 0x100), 0x1cf6edd694ac30cc90d09aadba508c800dbc7b64072d066291d20368230d976d, 0x21d66f0e2295ae954494f25889f9319cc1b4df71eff3f46ba9e4631b43fd7c95, 1, 0x299cb47eeb747bc6606ec27fa6cfc1eadc4d5f696b149cd0f4e85192f9e11c4f)
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
                for { let k := 0 } lt(k, 0x600) { k := add(k, 0x80) } {
                    ta := mulmod(ta, mload(add(S0, add(k, 0x40))), p)
                    mstore(add(PF, shr(2, k)), ta)
                }
                {
                    let M0 := add(PF, 0x180)
                    mstore(M0, 0x20)
                    mstore(add(M0, 0x20), 0x20)
                    mstore(add(M0, 0x40), 0x20)
                    mstore(add(M0, 0x60), ta)
                    mstore(add(M0, 0x80), minus_2)
                    mstore(add(M0, 0xa0), p)
                    if iszero(staticcall(gas(), 0x05, M0, 0xc0, M0, 0x20)) { revert(0, 0) }
                    ta := mload(M0)
                }
                for { let k := 0x580 } lt(k, 0x600) { k := sub(k, 0x80) } {
                    let pt := add(S0, k)
                    let tb := ta
                    if k { tb := mulmod(ta, mload(add(PF, sub(shr(2, k), 0x20))), p) }
                    ta := mulmod(ta, mload(add(pt, 0x40)), p)
                    st4(add(T, shl(7, byte(shr(7, k), _DMAP))), mulmod(mload(pt), tb, p), mulmod(mload(add(pt, 0x20)), tb, p))
                }
            }

            let T0 := mload(0x40)
            st4(T0, 0, 1)
            st4(add(T0, 0x80), B8X, B8Y)
            st4(add(T0, 0x100), 0x22b9a68c2ed57742dbc273e1f2e4924057c3d286a1f7b7bcead33e9ed2acb925, 0x01666cafbf0a30da8b9ebeaf848a1da067a892296f1043188e1705402b6d6853)
            st4(add(T0, 0x180), 0x1cf6edd694ac30cc90d09aadba508c800dbc7b64072d066291d20368230d976d, 0x21d66f0e2295ae954494f25889f9319cc1b4df71eff3f46ba9e4631b43fd7c95)
            buildQ(T0, mulmod(SQ, Ax, p), Ay)

            X := 0
            Y := 1
            Z := 1
            let T := 0
            Ax := 252
            if iszero(or(scalar_u, scalar_v)) { Ax := not(0) } // u = v = 0: identity
            for {} and(lt(Ax, 256), iszero(or(shr(Ax, scalar_u), shr(Ax, scalar_v)))) { Ax := sub(Ax, 4) } {}
            Ay := T0
            for {} lt(Ax, 256) { Ax := sub(Ax, 4) } {
                {
                    let A := mulmod(X, X, p)
                    let B := mulmod(Y, Y, p)
                    let E := add(mulmod(add(X, Y), add(X, Y), p), sub(_2P, add(A, B)))
                    let G := add(B, sub(p, A))
                    let F := add(G, sub(p, mulmod(add(Z, Z), Z, p)))
                    X := mulmod(E, F, p)
                    Y := mulmod(G, sub(_2P, add(A, B)), p)
                    Z := mulmod(F, G, p)
                }
                {
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
                }
                {
                    let e := or(and(shl(5, shr(Ax, scalar_u)), 0x180), and(shl(7, shr(Ax, scalar_v)), 0x600))
                    if e {
                        e := add(Ay, e)
                        {
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
                        }
                    }
                }
                {
                    let A := mulmod(X, X, p)
                    let B := mulmod(Y, Y, p)
                    let E := add(mulmod(add(X, Y), add(X, Y), p), sub(_2P, add(A, B)))
                    let G := add(B, sub(p, A))
                    let F := add(G, sub(p, mulmod(add(Z, Z), Z, p)))
                    X := mulmod(E, F, p)
                    Y := mulmod(G, sub(_2P, add(A, B)), p)
                    Z := mulmod(F, G, p)
                }
                {
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
                }
                {
                    let e := or(and(shl(7, shr(Ax, scalar_u)), 0x180), and(shl(9, shr(Ax, scalar_v)), 0x600))
                    if e {
                        e := add(Ay, e)
                        {
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
                        }
                    }
                }
            }
        }
    }

    /**
     * @dev Raw circomlib EdDSA check S.B8 == R8 + 8.hm.A (hm already hashed), circomlib coordinates,
     *      same input checks as circomlibjs (A and R8 on curve, S < subOrder). No inversion besides the table one.
     */
    function verifyRaw(uint256 hm, uint256 R8x, uint256 R8y, uint256 S, uint256 Ax, uint256 Ay)
        internal
        view
        returns (bool)
    {
        if (S >= n || !ed_isOnCurve(R8x, R8y) || !ed_isOnCurve(Ax, Ay)) return false;
        // result (X:Y:Z) stored in (Ax, Ay, S)
        (Ax, Ay, S) = ed_mulmuladd8(Ax, Ay, S, addmod(n - hm % n, 0, n));
        return Ax == mulmod(mulmod(SQ, R8x, p), S, p) && Ay == mulmod(R8y, S, p);
    }
} //EOF
