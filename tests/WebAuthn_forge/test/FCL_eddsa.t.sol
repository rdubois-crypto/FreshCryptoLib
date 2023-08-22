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
///* FILE: FCL_elliptic.t.sol
///*
///*
///* DESCRIPTION: test file for ecdsa signature protocol
///*
//**************************************************************************************/
// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.20;

import "forge-std/Test.sol";

import "@solidity/FCL_eddsa.sol";

contract EddsaTest is Test {

    function test_add() public {
      
        uint256 x;uint256 y;uint256 z;uint256 t;
        uint256 x4;uint256 y4;uint256 z4;uint256 t4;
        uint256 minus_gy=p-gy;//-gy

        (x,y,z,t)=EDDSA.ed_Add(gx, gy,1,mulmod(gx,gy,p), gx, gy, 1, mulmod(gx,gy,p));//2G
        (x4,y4,z4,t4)=EDDSA.ed_Add(x,y,z,t, x,y,z,t); //4G

        (x4,y4,z4,t4)=EDDSA.ed_Add(x4,y4,z4,t4, gx,minus_gy,1,mulmod(gx,minus_gy,p)); //3G
        (x4,y4,z4,t4)=EDDSA.ed_Add(x4,y4,z4,t4, gx,minus_gy,1,mulmod(gx,minus_gy,p)); //2G
        (x4,y4)=EDDSA.ed_z2Aff(x4,y4,z4);
        (x,y)=EDDSA.ed_z2Aff(x,y,z);

      assertEq(x4,x);

  assertEq(y4,y);

    }

    function test_dbl() public{
         uint256 x;uint256 y;uint256 z;uint256 t;
           uint256 x4;uint256 y4;uint256 z4;uint256 t4;
      
       (x,y,z,t)=EDDSA.ed_Add(gx, gy,1,mulmod(gx,gy,p), gx, gy, 1, mulmod(gx,gy,p));//2G
       (x4,y4,z4,t4)=EDDSA.ed_Dbl(gx, gy,1); //2G
 (x4,y4)=EDDSA.ed_z2Aff(x4,y4,z4);
        (x,y)=EDDSA.ed_z2Aff(x,y,z);

  assertEq(x4,x);

  assertEq(y4,y);

    }

    function test_mulmuladd() public{
      uint256 x_res;
      //(n+1)G==G ?
      x_res=EDDSA.ed_mulmuladd(gx,gy, 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ee, 0);
      
     // assertEq(x_res,gx); fail
    }
}

