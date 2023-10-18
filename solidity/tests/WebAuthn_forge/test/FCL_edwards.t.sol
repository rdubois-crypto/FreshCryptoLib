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
///* FILE: FCL_Edwards.t.sol
///*
///*
///* DESCRIPTION: test file for Edwards curves
///*
//**************************************************************************************/
// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.20;

import "forge-std/Test.sol";

import "@solidity/FCL_edwards.sol";

contract EdwardsTest is Test {
    // Multiply an elliptic curve point by a scalar
    function multiplyScalar(uint256 x0, uint256 y0, uint256 scalar)
        public
        returns (
            //returns (uint256 x3, uint256 y3, uint256 z3, uint256 t3){
            uint256 x1,
            uint256 y1
        )
    {
        if (scalar == 0) {
            return (0, 1);
        } else if (scalar == 1) {
            return (x0, y0);
        }

        uint256 base2X = x0;
        uint256 base2Y = y0;
        uint256 base2Z = 1;
        uint256 base2t = mulmod(x0, y0, p);
        //uint256 x1; uint256 y1;

        x1 = x0;
        y1 = y0;
        uint256 z1 = 1;
        uint256 t1 = base2t;

        if (scalar % 2 == 0) {
            x1 = 0;
            y1 = 1;
            z1 = 1;
            t1 = 0;
        }

        scalar = scalar >> 1;

        while (scalar > 0) {
            //(base2X, base2Y, base2Z,base2t ) = Edwards.ed_Dbl(base2X, base2Y, base2Z);
            (base2X, base2Y, base2Z, base2t) = Edwards.ed_Dbl(base2X, base2Y, base2Z);

            if (scalar & 1 == 1) {
                (x1, y1, z1, t1) = Edwards.ed_Add(base2X, base2Y, base2Z, base2t, x1, y1, z1, t1);
            }

            scalar = scalar >> 1;
        }

        return Edwards.ed_z2Aff(x1, y1, z1);
        //return (x1, y1, z1,t1);
    }

    function test_z2Aff(uint256 scramble) public {
        vm.assume(scramble != 0);
        vm.assume(scramble < p);

        if (scramble != 0) {
            uint256 x = mulmod(gx, scramble, p);
            uint256 y = mulmod(gy, scramble, p);
            uint256 z = scramble;
            // uint256 t = mulmod(x, y, p);
            (x, y) = Edwards.ed_z2Aff(x, y, z);
            assertEq(x, gx);
            assertEq(y, gy);
        }
    }

    function test_add() public {
        uint256 x;
        uint256 y;
        uint256 z;
        uint256 t;
        uint256 x4;
        uint256 y4;
        uint256 z4;
        uint256 t4;
        uint256 minus_gx = p - gx; //-gy
        uint256 mt;

        (x, y, z, t) = Edwards.ed_Add(gx, gy, 1, mulmod(gx, gy, p), gx, gy, 1, mulmod(gx, gy, p));

        for (uint256 i = 1; i < 100; i++) {
            (x, y, z, t) = Edwards.ed_Add(x, y, z, t, x, y, z, t); //P=2P
            minus_gx = p - x;
            mt = mulmod(minus_gx, y, p);
            (x4, y4, z4, t4) = Edwards.ed_Add(x, y, z, t, x, y, z, t); //2P
            (x4, y4, z4, t4) = Edwards.ed_Add(x, y, z, t, x, y, z, t); //4P

            (x4, y4, z4, t4) = Edwards.ed_Add(x4, y4, z4, t4, minus_gx, y, z, mt); //4G-G=3P

            (x4, y4, z4, t4) = Edwards.ed_Add(x4, y4, z4, t4, minus_gx, y, z, mt); //4G-G=2P

            (x4, y4, z4, t4) = Edwards.ed_Add(x4, y4, z4, t4, minus_gx, y, z, mt); //4G-G=P
        }

        (x, y) = Edwards.ed_z2Aff(x4, y4, z4);
        (x4, y4) = Edwards.ed_z2Aff(x4, y4, z4);
        assertEq(x4, x);
    }

    function test_dbl() public {
        uint256 x;
        uint256 y;
        uint256 z;
        uint256 t;
        uint256 x4;
        uint256 y4;
        uint256 z4;
        uint256 t4;

        (x, y, z, t) = Edwards.ed_Add(gx, gy, 1, mulmod(gx, gy, p), gx, gy, 1, mulmod(gx, gy, p)); //2G
        (x4, y4, z4, t4) = Edwards.ed_Dbl(gx, gy, 1); //2G
        (x4, y4) = Edwards.ed_z2Aff(x4, y4, z4);
        (x, y) = Edwards.ed_z2Aff(x, y, z);

        assertEq(x4, x);

        assertEq(y4, y);
    }

    function test_oncurve() public {
        assertEq(Edwards.ed_isOnCurve(gx, gy, 1), true);
    }

    function test_mul() public {
        uint256 x_res1;
        // uint256 x_res;
        uint256 y_res1;
        // uint256 z_res1;
        // uint256 t_res1;
        //(n-1)G==-G ?
        (x_res1, y_res1) = multiplyScalar(gx, gy, 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ec);

        console.log("---res mul:", x_res1, y_res1);

        assertEq(p - x_res1, gx);
        //(n+1)G==-G ?
        (x_res1, y_res1) = multiplyScalar(gx, gy, 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ee);
        assertEq(x_res1, gx);
    }

    function test_mulmuladd() public {
        uint256 x_res1;
        //uint256 x_res;
        uint256 y_res1;
        // uint256 z_res1;
        //uint256 t_res1;
        //(n-1)G==-G ?
        (x_res1) = Edwards.ed_mulmuladd(gx, gy, 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ec, 0);

        console.log("---res mul:", x_res1, y_res1);

        assertEq(p - x_res1, gx);
    }

    function test_Sqrtmod() public {
        uint256 val = mulmod(gx, gx, p);
        uint256 rac = SqrtMod(val);
        console.log("rac=", rac);
        assertEq(mulmod(rac, rac, p), val);
    }

    function test_ed_decompress() public {
        uint256 x = Edwards.ed_decompress(gy, 1);
        if (Edwards.ed_isOnCurve(x, gy, 1) != true) {
            revert();
        }
        if ((x != gx) && (p - x != gx)) {
            revert();
        }
    }

    function test_mulmuladd_opp() public {
        uint256 x = Edwards.ed_mulmuladd(
            gx,
            p - gy,
            0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3e0,
            0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3e0
        );
        assertEq(x, 0);
    }
}
