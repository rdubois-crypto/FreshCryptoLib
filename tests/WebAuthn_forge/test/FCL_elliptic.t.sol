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
///* DESCRIPTION: test file for elliptic primitives
///*
//**************************************************************************************/
pragma solidity >=0.8.19 <0.9.0;

import "forge-std/Test.sol";
import "../src/FCL_elliptic.sol";
import "../src/FCL_Webauthn.sol";


//external implementation to bench
import "../src/ECops.sol";
import "../src/Secp256r1.sol";
import "../src/Secp256r1_maxrobot.sol";

contract ArithmeticTest is Test {
    //curve prime field modulus
    uint256 constant p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;
    //short weierstrass first coefficient
    uint256 constant a = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC;
    //short weierstrass second coefficient
    uint256 constant b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B;
    //generating point affine coordinates
    uint256 constant gx = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296;
    uint256 constant gy = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5;
    //curve order (number of points)
    uint256 constant n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;
    /* -2 mod p constant, used to speed up inversion and doubling (avoid negation)*/
    uint256 constant minus_2 = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFD;
    /* -2 mod n constant, used to speed up inversion*/
    uint256 constant minus_2modn = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC63254F;

    uint256 constant minus_1 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF;

    uint256 constant _prec_address = 0xcaca;
    uint256 constant _NUM_TEST_ECMULMULADD=100;
    uint256 constant _NUM_TEST_DBL=100;

    function test_Fuzz_InVmodn(uint256 i_u256_a) public {
        vm.assume(i_u256_a < FCL_Elliptic_ZZ.n);
        vm.assume(i_u256_a != 0);

        uint256 res = FCL_Elliptic_ZZ.FCL_nModInv(i_u256_a);

        assertEq(mulmod(res, i_u256_a, FCL_Elliptic_ZZ.n), 1);
    }

    function test_Fuzz_InVmodp(uint256 i_u256_a) public {
        vm.assume(i_u256_a < FCL_Elliptic_ZZ.p);
        vm.assume(i_u256_a != 0);

        uint256 res = FCL_Elliptic_ZZ.FCL_pModInv(i_u256_a);

        assertEq(mulmod(res, i_u256_a, FCL_Elliptic_ZZ.p), 1);
    }
    //ecAff_isOnCurve


    //check consistency of ecmulmuladd and Add
    function test_invariant_FCL_Ecmulmuladd() public {
        uint256 ecpoint_Rx=0;
        uint256 ecpoint_Ry=0;
         uint256 ecpoint_Rzz=1;
        uint256 ecpoint_Rzzz=1;
       uint256  checkpointGasLeft ;
	uint256  checkpointGasLeft2 ;

        //Uncomment the library to test/bench here
        //(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz, ecpoint_Rzzz )= FCL_Elliptic_ZZ.ecZZ_Dbl(gx, gy,1,1);
        (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz )= ECops.twiceProj(gx, gy,1);


        checkpointGasLeft=gasleft() ;

        for(uint256 i=3;i<=_NUM_TEST_ECMULMULADD;i++){
          // (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz, ecpoint_Rzzz )= FCL_Elliptic_ZZ.ecZZ_AddN(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz, ecpoint_Rzzz , gx, gy);

             (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz)= ECops.addProj(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz,gx,gy,1);
        }
        checkpointGasLeft2=gasleft() ;
        console.log("Add number test and gas cost :", _NUM_TEST_ECMULMULADD, checkpointGasLeft - checkpointGasLeft2 - 100);

//        (ecpoint_Rx, ecpoint_Ry)=FCL_Elliptic_ZZ.ecZZ_SetAff(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz, ecpoint_Rzzz);
        (ecpoint_Rx, ecpoint_Ry)=ECops.toAffinePoint(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz);

        //generate a small multiple of G using scalar
        assertEq(ecpoint_Rx, FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(0,0, _NUM_TEST_ECMULMULADD, 0));

          checkpointGasLeft=gasleft() ;
        for(uint256 i=3;i<=_NUM_TEST_DBL;i++){
           //(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz, ecpoint_Rzzz )=FCL_Elliptic_ZZ.ecZZ_Dbl(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz, ecpoint_Rzzz );

              (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz )= ECops.twiceProj(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz);

        }
          checkpointGasLeft2=gasleft() ;
        console.log("Dbl number test and gas cost :", _NUM_TEST_ECMULMULADD, checkpointGasLeft - checkpointGasLeft2 - 100);


    }



}

