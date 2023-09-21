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
///* FILE: FCL_elliptic.sol
///*
///*
///* DESCRIPTION: Implementation of RFC8032 using FCL Shamir's trick
///*
//**************************************************************************************/
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import "./FCL_ed25519.sol" as Curve;
import {p, gx, gy, n, d, deux_d, MINUS_2, MINUS_1, MODEXP_PRECOMPILE, SqrtMod, pModInv} from "./FCL_ed25519.sol";

library Edwards {
    using {Curve.pModInv, Curve.nModInv} for uint256;
    /**
     * @notice Convert from projective coordinates to affine coordinates
     *
     * @param x The X-coordinate of the point in XYZZ representation
     * @param y The Y-coordinate of the point in XYZZ representation
     * @param z The ZZ value of the point in XYZZ representation
     * @return x1 The X-coordinate of the point in affine representation
     * @return y1 The Y-coordinate of the point in affine representation
     */
    function ed_z2Aff(uint256 x, uint256 y, uint256 z) internal returns (uint256 x1, uint256 y1) {
        // 1/zzz
        uint256 zInv = z.pModInv();

        // Y/zzz -- OUTPUT
        y1 = mulmod(y, zInv, p);
        // X/zz -- OUTPUT
        x1 = mulmod(x, zInv, p);
    }

    /**
     * @notice Extract x coordinates from y
     *
     * @param y The y-coordinate of the point in affine representation
     * @return x The X-coordinate of the point in affine representation
    */
    function ed_decompress(uint256 y, uint256 sign) internal returns (uint256 x)
    {
        uint256 x2;
        uint256 y2=mulmod(y,y,p);
         x2 = mulmod(addmod(y2,MINUS_1,p) , pModInv( addmod(mulmod(d,y2,p),1,p) ) ,p);
        x=SqrtMod(x2);
        if((x&1)!=sign){
            x=p-x;
        }
    }

    function ed_isOnCurve(uint256 x, uint256 y, uint256 z) internal returns (bool b) {
        (x, y) = ed_z2Aff(x, y, z);

        uint256 x2 = mulmod(x, x, p);
        uint256 y2 = mulmod(y, y, p);
        uint256 dy2 = mulmod(d, y2, p); //dy2
        uint256 dy2x2 = mulmod(x2, dy2, p); //dx2y2

        dy2x2 = addmod(x2, addmod(1, dy2x2, p), p); //x2+dx2y2+1 == y2 ?

        return (addmod(p - y2, dy2x2, p) == 0);
    }

    function ed_Add(uint256 x1, uint256 y1, uint256 z1, uint256 t1, uint256 x2, uint256 y2, uint256 z2, uint256 t2)
        internal
        pure
        returns (uint256 x3, uint256 y3, uint256 z3, uint256 t3)
    {
        unchecked {
            assembly {
                x3 := addmod(y1, sub(p, x1), p) //   = (Y1-X1)
                t3 := addmod(y2, sub(p, x2), p) //     (Y2-X2)

                y3 := mulmod(x3, t3, p) // A = (Y1-X1)*(Y2-X2)
                x3 := mulmod(addmod(x1, y1, p), addmod(x2, y2, p), p) // B = (Y1+X1)*(Y2+X2)

                let P3 := mulmod(mulmod(t1, t2, p), deux_d, p) //  C = T1*2*d*T2
                t3 := mulmod(z1, z2, p)
                let P4 := mulmod(t3, 2, p) //   D = Z1*2*Z2
                let P5 := addmod(x3, sub(p, y3), p) //E = B-A
                let P6 := addmod(x3, y3, p) //H=B+A
                t3 := mulmod(P5, P6, p) //T3 = E*H,  not required for Dbl input
                z3 := addmod(P4, sub(p, P3), p) //F = D-C
                x3 := mulmod(P5, z3, p) // X3 = E*F
                P5 := addmod(P4, P3, p) //G = D+C
                y3 := mulmod(P5, P6, p) //Y3 = G*H
                z3 := mulmod(z3, P5, p) //  Z3 = F*G
            }
        }
        return (x3, y3, z3, t3);
    }

    function ed_AddN(uint256 x1, uint256 y1, uint256 z1, uint256 t1, uint256 x2, uint256 y2, uint256 t2)
        internal
        pure
        returns (uint256 x3, uint256 y3, uint256 z3, uint256 t3)
    {
        unchecked {
            assembly {
                x3 := addmod(y1, sub(p, x1), p) //   = (Y1-X1)
                t3 := addmod(y2, sub(p, x2), p) //     (Y2-X2)

                y3 := mulmod(x3, t3, p) // A = (Y1-X1)*(Y2-X2)
                x3 := mulmod(addmod(x1, y1, p), addmod(x2, y2, p), p) // B = (Y1+X1)*(Y2+X2)

                let P3 := mulmod(mulmod(t1, t2, p), deux_d, p) //  C = T1*2*d*T2
                //  t3 := mulmod(z1, z2, p), z2 is normalized
                let P4 := mulmod(z1, 2, p) //   D = Z1*2*Z2
                let P5 := addmod(x3, sub(p, y3), p) //E = B-A
                let P6 := addmod(x3, y3, p) //H=B+A
                t3 := mulmod(P5, P6, p) //T3 = E*H,  not required for Dbl input
                z3 := addmod(P4, sub(p, P3), p) //F = D-C
                x3 := mulmod(P5, z3, p) // X3 = E*F
                P5 := addmod(P4, P3, p) //G = D+C
                y3 := mulmod(P5, P6, p) //Y3 = G*H
                z3 := mulmod(z3, P5, p) //  Z3 = F*G
            }
        }
        return (x3, y3, z3, t3);
    }

    function ed_Dbl(uint256 x1, uint256 y1, uint256 z1)
        internal
        pure
        returns (uint256 x3, uint256 y3, uint256 z3, uint256 t3)
    {
        unchecked {
            assembly {
                x3 := mulmod(x1, x1, p) //A = X1^2
                t3 := mulmod(y1, y1, p) //B = Y1^2
                let P6 := addmod(x3, t3, p) //H=B+A
                let P5 := addmod(x1, y1, p) //   (X1+Y1)
                P5 := mulmod(P5, P5, p)
                P5 := addmod(P6, sub(p, P5), p) //E = H-(X1+Y1)^2

                y3 := mulmod(2, mulmod(z1, z1, p), p) // C = 2*Z1^2
                let P0 := addmod(x3, sub(p, t3), p) //G = A-B
                y3 := addmod(y3, P0, p) //F =C+G, c

                t3 := mulmod(P5, P6, p) //T3 = E*H
                z3 := mulmod(P0, y3, p) //  Z3 = F*G
                x3 := mulmod(y3, P5, p) // X3 = E*F
                y3 := mulmod(P0, P6, p) //Y3 = G*H
            }
        }
        return (x3, y3, z3, t3);
    }

    function ed_AddAff(uint256 x1, uint256 y1, uint256 x2, uint256 y2) internal returns (uint256 x3, uint256 y3) {
        uint256 z3;
        uint256 t3;
        (x3, y3, z3, t3) = ed_Add(x1, y1, 1, mulmod(x1, y1, p), x2, y2, 1, mulmod(x2, y2, p));
        return ed_z2Aff(x3, y3, z3);
    }
    /**
     * @dev Computation of uG+vQ using Strauss-Shamir's trick, G basepoint, Q public key
     * tested KO for now
     */

    function ed_mulmuladd(
        uint256 Q0,
        uint256 Q1, //affine rep for input point Q
        uint256 scalar_u,
        uint256 scalar_v
    ) internal returns (uint256 X) {
        uint256 zz;
        uint256 t;
        uint256 Y;
        uint256 index = 255;
        uint256[6] memory T;
        // uint256 H0;
        // uint256 H1;

        unchecked {
            if (scalar_u == 0 && scalar_v == 0) return 0;

            (T[0], T[1]) = ed_AddAff(gx, gy, Q0, Q1);

            assembly {
                for { let T4 := add(shl(1, and(shr(index, scalar_v), 1)), and(shr(index, scalar_u), 1)) } eq(T4, 0) {
                    index := sub(index, 1)
                    T4 := add(shl(1, and(shr(index, scalar_v), 1)), and(shr(index, scalar_u), 1))
                } {}
                zz := add(shl(1, and(shr(index, scalar_v), 1)), and(shr(index, scalar_u), 1))

                if eq(zz, 1) {
                    X := gx
                    Y := gy
                    t := mulmod(gx, gy, p)
                }
                if eq(zz, 2) {
                    X := Q0
                    Y := Q1
                    t := mulmod(Q0, Q1, p)
                }
                if eq(zz, 3) {
                    X := mload(T)
                    Y := mload(add(T, 0x20))
                    t := mulmod(X, Y, p)
                }

                index := sub(index, 1)
                zz := 1

                for {} gt(MINUS_1, index) { index := sub(index, 1) } {
                    // inlined EcZZ_Dbl
                    let A := mulmod(X, X, p) //A = X1^2
                    let B := mulmod(Y, Y, p) //B = Y1^2
                    let H := addmod(A, B, p) //H=B+A
                    let E := addmod(X, Y, p) //   E = H-(X1+Y1)^2, X, Y available
                    E := mulmod(E, E, p)
                    E := addmod(H, sub(p, E), p) //E
                    Y := mulmod(2, mulmod(zz, zz, p), p) // C = 2*Z1^2
                    let G := addmod(A, sub(p, B), p) //G = A-B
                    Y := addmod(Y, G, p) //F =C+G
                    t := mulmod(E, H, p) //T3 = E*H
                    zz := mulmod(G, Y, p) //  Z3 = F*G
                    X := mulmod(E, Y, p) // X3 = E*F
                    Y := mulmod(G, H, p) //Y3 = G*H

                    {
                        //value of dibit
                        let T4 := add(shl(1, and(shr(index, scalar_v), 1)), and(shr(index, scalar_u), 1))

                        if iszero(T4) { continue } // if T4!=0
                        if eq(T4, 1) {
                            A := gx
                            B := gy
                        }
                        if eq(T4, 2) {
                            A := Q0
                            B := Q1
                        }
                        if eq(T4, 3) {
                            A := mload(T)
                            B := mload(add(T, 0x20))
                        }

                        // inlined EcZZ_AddN
                        H := mulmod(mulmod(t, mulmod(A, B, p), p), deux_d, p) //  C = T1*2*d*T2

                        T4 := addmod(X, sub(p, Y), p) //   = (Y1-X1)
                        E := addmod(A, sub(p, B), p) //     (Y2-X2)

                        t := mulmod(E, T4, p) // A = (Y1-X1)*(Y2-X2)
                        B := mulmod(addmod(X, Y, p), addmod(A, B, p), p) // B = (Y1+X1)*(Y2+X2)

                        //  t3 := mulmod(z1, z2, p), z2 is normalized
                        G := mulmod(zz, 2, p) //   D = Z1*2*Z2
                        E := addmod(B, sub(p, t), p) //E = B-A

                        A := addmod(G, H, p) //G = D+C
                        zz := addmod(G, sub(p, H), p) //F = D-C
                        H := addmod(B, t, p) //H=B+A
                        X := mulmod(E, zz, p) // X3 = E*F

                        Y := mulmod(A, H, p) //Y3 = G*H
                        zz := mulmod(zz, A, p) //  Z3 = F*G
                    }
                } //end loop
                mstore(add(T, 0x60), zz)
                //(X,Y)=ecZZ_SetAff(X,Y,zz, zzz);
                //T[0] = inverseModp_Hard(T[0], p); //1/zzz, inline modular inversion using precompile:
                // Define length of base, exponent and modulus. 0x20 == 32 bytes
                mstore(T, 0x20)
                mstore(add(T, 0x20), 0x20)
                mstore(add(T, 0x40), 0x20)
                // Define variables base, exponent and modulus
                //mstore(add(pointer, 0x60), u)
                mstore(add(T, 0x80), MINUS_2)
                mstore(add(T, 0xa0), p)

                // Call the precompiled contract 0x05 = ModExp
                if iszero(call(not(0), 0x05, 0, T, 0xc0, T, 0x20)) { revert(0, 0) }

                //Y:=mulmod(Y,zzz,p)//Y/zzz
                //zz :=mulmod(zz, mload(T),p) //1/z
                //zz:= mulmod(zz,zz,p) //1/zz
                X := mulmod(X, mload(T), p) //X/zz
            } //end assembly
        } //end unchecked

        return X;
    }
}
