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
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import "forge-std/Test.sol";
import "@solidity/FCL_elliptic.sol";
import "@solidity/FCL_ecdsa.sol";

//testing edge case as suggested by Mikhail in commit 5d3c3f77f0d296f095bb071e7df5278a1c0cc1be
contract edgemultTest is Test {
 /* vector from http://point-at-infinity.org/ecc/nisttv
//k = 115792089210356248762697446949407573529996955224135760342422259061068512044367
//x = 7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978
//y = F888AAEE24712FC0D6C26539608BCF244582521AC3167DD661FB4862DD878C2E*/
//edge case for Shamir 
function test_edgeMul() public returns (bool)
{
 uint256[3] memory vec=[
  115792089210356248762697446949407573529996955224135760342422259061068512044367,
  0x7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978,
  0xF888AAEE24712FC0D6C26539608BCF244582521AC3167DD661FB4862DD878C2E
 ];
 uint256 resX;
 uint256 resY;
 uint256[4] memory Q=[uint256(0),0,0,0];

 //(resX, resY)=ec_scalarmulN(vec[0], vec[1], vec[2]);
 resX=FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(Q[0], Q[1], vec[0], 0);
 assertEq(0x7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978, resX);

 return true;
}

}
