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
///* FILE: FCL_eddsa.t.sol
///*
///*
///* DESCRIPTION: test file for eddsa signature protocol
///*
//**************************************************************************************/
// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.20;

import "forge-std/Test.sol";

import "@solidity/FCL_eddsa.sol";

contract EddsaTest is Test {

// Multiply an elliptic curve point by a scalar
    function multiplyScalar(uint256 x0, uint256 y0, uint256 scalar) public  
    //returns (uint256 x3, uint256 y3, uint256 z3, uint256 t3){
    returns (uint256 x1, uint256 y1) {
        if (scalar == 0) {
            return (0,1);
        } else if (scalar == 1) {
            return (x0, y0);
        } 

        uint256 base2X = x0;
        uint256 base2Y = y0;
        uint256 base2Z = 1;
        uint256 base2t = mulmod(x0,y0,p);
        //uint256 x1; uint256 y1;
       
        x1 = x0;
        y1 = y0;
         uint256 z1 = 1;
        uint256 t1 = base2t;
      
        if (scalar % 2 == 0) {
            x1 =0;y1=1;z1=1;t1=0;
        }

        scalar = scalar >> 1;

        while (scalar > 0) {
            //(base2X, base2Y, base2Z,base2t ) = EDDSA.ed_Dbl(base2X, base2Y, base2Z);
            (base2X, base2Y, base2Z,base2t ) = EDDSA.ed_Add(base2X, base2Y, base2Z,base2t,base2X, base2Y, base2Z,base2t);

            if (scalar &1==1) {
                (x1, y1, z1,t1) = EDDSA.ed_Add(base2X, base2Y, base2Z,base2t, x1, y1, z1,t1);
            }
           
            scalar = scalar >> 1;
        }
 
        return EDDSA.ed_z2Aff(x1, y1, z1);
        //return (x1, y1, z1,t1);
    }

    function test_z2Aff(uint256 scramble) public
    {
      if(scramble!=0)
      {
        uint256 x=mulmod(gx,scramble,p);
        uint256 y=mulmod(gy,scramble,p);
        uint256 z=scramble;
        uint256 t=mulmod(x,y,p) ;
        (x,y)=EDDSA.ed_z2Aff(x,y,z);
        assertEq(x,gx);
        assertEq(y,gy);
        
      }

    }

    function test_add() public {
      
        uint256 x;uint256 y;uint256 z;uint256 t;
        uint256 x4;uint256 y4;uint256 z4;uint256 t4;
        uint256 minus_gx=p-gx;//-gy

        (x,y,z,t)=EDDSA.ed_Add(gx, gy,1,mulmod(gx,gy,p), gx, gy, 1, mulmod(gx,gy,p));//2G
       
        console.log("is oncurve post add=",EDDSA.ed_isOnCurve(x,y,z,t));
        console.log("x=",x);
        (x4,y4,z4,t4)=EDDSA.ed_Add(x,y,z,t, x,y,z,t); //4G
 console.log("is oncurve post add 4G=",EDDSA.ed_isOnCurve(x4,y4,z4,t4));
       
        (x4,y4,z4,t4)=EDDSA.ed_Add(x4,y4,z4,t4, minus_gx,gy,1,mulmod(minus_gx,gy,p)); //4G-G=3G
         console.log("is oncurve post add 3G=",EDDSA.ed_isOnCurve(x4,y4,z4,t4));
       
        (x4,y4,z4,t4)=EDDSA.ed_Add(x4,y4,z4,t4, minus_gx,gy,1,mulmod(minus_gx,gy,p)); //3G-G=2G
         console.log("is oncurve post add 2G=",EDDSA.ed_isOnCurve(x4,y4,z4,t4));
       
        (x4,y4,z4,t4)=EDDSA.ed_Add(x4,y4,z4,t4, minus_gx,gy,1,mulmod(minus_gx,gy,p)); //2G-G=G  
         console.log("is oncurve post add G=",EDDSA.ed_isOnCurve(x4,y4,z4,t4));
       
        console.log("P'=",x4,y4);
        console.log("",z4,t4);

        (x4,y4)=EDDSA.ed_z2Aff(x4,y4,z4);
       
        console.log("P=",x4,y4);
      
      assertEq(x4,gx);
   
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

    function test_oncurve() public
    {

      assertEq(EDDSA.ed_isOnCurve(gx,gy,1,mulmod(gx,gy,p)),true);
    }

    function test_mulmuladd() public{
      uint256 x_res1; uint256 x_res;uint256 y_res1;
       uint256 z_res1;  uint256 t_res1;
      //(n-1)G==-G ?
      (x_res1,y_res1)=multiplyScalar(gx,gy,0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ec );
      //x_res=EDDSA.ed_mulmuladd(gx,gy, 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3d, 0);
      
      console.log("---res mul:", x_res1, y_res1);
     
      assertEq(p-x_res1,gx); 
       //(n+1)G==-G ?
      (x_res1,y_res1)=multiplyScalar(gx,gy,0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ee );
      assertEq(x_res1,gx); 
     

    }
}

