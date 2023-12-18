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
import "@solidity/FCL_ecdsa_utils.sol";

//testing edge case as suggested by Mikhail in commit 5d3c3f77f0d296f095bb071e7df5278a1c0cc1be
contract edgemultTest is Test {
    //naive mul for low weight coeff
    function ecAff_naivemul(uint256 a, uint256 x, uint256 y) private view returns (uint256 ax, uint256 ay) {
        while (a > 0) {
            (ax, ay) = FCL_Elliptic_ZZ.ecAff_add(ax, ay, x, y);
            a = a - 1;
        }
    }

    /* vector from http://point-at-infinity.org/ecc/nisttv
    //k = 115792089210356248762697446949407573529996955224135760342422259061068512044367
    //x = 7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978
    //y = F888AAEE24712FC0D6C26539608BCF244582521AC3167DD661FB4862DD878C2E*/
    //edge case for Shamir
    function test_EdgeMul_InlinedDbl() public returns (bool) {
        uint256[3] memory vec = [
            115792089210356248762697446949407573529996955224135760342422259061068512044367,
            0x7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978,
            0xF888AAEE24712FC0D6C26539608BCF244582521AC3167DD661FB4862DD878C2E
        ];
        uint256 resX;
        uint256 resY;
        uint256[4] memory Q = [uint256(0), 0, 0, 0];

        //(resX, resY)=ec_scalarmulN(vec[0], vec[1], vec[2]);
        resX = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(Q[0], Q[1], vec[0], 0);
        assertEq(0x7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978, resX);

        //edge case from niordell
        uint256[4] memory vec2 = [
            102369864249653057322725350723741461599905180004905897298779971437827381725266, //x
            14047598098721058250371778545974983789701612908526165355421494088134814672697, //y
            94632330233094393099906091027057584450760066982961548963789323460936666616340, //u
            23658082558273598274976522756764396112690016745740387240947330865234166656879
        ]; //v

        //expected result using FCL_elliptic.sage, dark side:
        //_G_POINT*94632330233094393099906091027057584450760066982961548963789323460936666616340+_G_CURVE(102369864249653057322725350723741461599905180004905897298779971437827381725266, 14047598098721058250371778545974983789701612908526165355421494088134814672697)*23658082558273598274976522756764396112690016745740387240947330865234166656879
        //(93995665850302450053183256960521438033484268364047930968443817833761593125805 : 60765861213361593633751918097312828188566711467069305801019119884414110226811 : 1)

        resX = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(vec2[0], vec2[1], vec2[2], vec2[3]);
        console.log("resX=%x", resX);

        assertEq(93995665850302450053183256960521438033484268364047930968443817833761593125805, resX);
    }

    //edge case: where Q=n-1, making H precomputed value neutral
    function test_EdgeCase_keqnm1() external {
        uint256 s = FCL_Elliptic_ZZ.n - 1;
        (uint256 Qx, uint256 Qy) = FCL_ecdsa_utils.ecdsa_derivKpub(s);
        console.log("Qx = %d", Qx);
        console.log("Qy = %d", Qy);

        uint256 u = 3;
        uint256 v = 1;

        uint256 r1 = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(Qx, Qy, u, v);
        console.log("r1 = %d", r1);

        (uint256 uGx, uint256 uGy) = ecAff_naivemul(u, FCL_Elliptic_ZZ.gx, FCL_Elliptic_ZZ.gy);
        (uint256 vQx, uint256 vQy) = ecAff_naivemul(v, Qx, Qy);
        (uint256 r2,) = FCL_Elliptic_ZZ.ecAff_add(uGx, uGy, vQx, vQy);
        console.log("r2 = %d", r2);

        (uint256 twoGx,) = ecAff_naivemul(2, FCL_Elliptic_ZZ.gx, FCL_Elliptic_ZZ.gy);
        console.log("2Gx = %d", twoGx);

        assertEq(r1, r2);
    }
}
