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
///* DESCRIPTION: modified XYZZ system coordinates for EVM elliptic point multiplication
///*  optimization
///*
//**************************************************************************************/

pragma solidity >=0.8.19 <0.9.0;


import "./FCL_ed25519.sol" as Curve;
import { p,  n, d, deux_d, MINUS_2, MINUS_1, MODEXP_PRECOMPILE } from "./FCL_ed25519.sol";

library EDDSA{
     using { Curve.pModInv, Curve.nModInv } for uint256;
  /**
     * @notice Convert from projective coordinates to affine coordinates
     *
     * @param x The X-coordinate of the point in XYZZ representation
     * @param y The Y-coordinate of the point in XYZZ representation
     * @param z The ZZ value of the point in XYZZ representation
     * @return x1 The X-coordinate of the point in affine representation
     * @return y1 The Y-coordinate of the point in affine representation
     */
    function z2Aff(uint256 x, uint256 y, uint256 z, uint256 zzz) internal returns (uint256 x1, uint256 y1) {
        // 1/zzz
        uint256 zInv = z.pModInv();

        // Y/zzz -- OUTPUT
        y1 = mulmod(y, zInv, p);
        // X/zz -- OUTPUT
        x1 = mulmod(x, zInv, p);
    }

    function ecZZ_Add(uint256 x1, uint256 y1, uint256 z1, uint256 t1, uint256 x2, uint256 y2, uint256 z2, uint256 t2)
        internal
        pure
        returns (uint256 x3, uint256 y3, uint256 z3, uint256 t3)
    {
        unchecked {
            assembly {
             x3:= addmod(y1, sub(p, x1),p)//   A = (Y1-X1)
             t3:= addmod(y2, sub(p, x2),p)//   A = (Y2-X2)
             let P2:= mulmod(x3, t3, p) // A = (Y1-X1)*(Y2-X2)
             let P0:= mulmod(addmod(x1,y1,p), addmod(x2,y2,p),p) // B = (Y1+X1)*(Y2+X2)
             let P1:= mulmod(t1,t2,p) 
             let  P3:= mulmod(P1,deux_d,p) //  C = T1*2*d*T2
             P1:= mulmod(z1,z2,p) 
             let P4:= mulmod(P1,2,p) //   D = Z1*2*Z2
             let  P5:=addmod(P0, sub(p,P2),p)  //E = B-A
             let P6:=addmod(P0,P2,p) //H=B+A
             t3:=mulmod(P5,P6,p) //T3 = E*H
             //F = D-C
             x3:=mulmod(P5,addmod(P4, sub(p,P3),p), p) // X3 = E*F
             P5:=addmod(P4, P3,p) //G = D+C
             y3:=mulmod(P5, P6,p) //Y3 = G*H
                           
                 
            }
        }
        return (x3, y3, z3, t3);
    }


    
}