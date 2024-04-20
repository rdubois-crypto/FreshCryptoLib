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
pragma solidity ^0.8.20;

import "forge-std/Test.sol";
import "@solidity/FCL_elliptic.sol";
import "@solidity/FCL_ecdsa.sol";
import "@solidity/FCL_ecdsa_utils.sol";
import "@solidity/FCL_Webauthn.sol";

//external implementation to bench
import "@external/ECops.sol";
import "@external/Secp256r1.sol";
import "@external/Secp256r1_maxrobot.sol";

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
    uint256 constant _NUM_TEST_ECMULMULADD = 1000;
    uint256 constant _NUM_TEST_DBL = 100;

    function test_ecAff_isOnCurve_returnsFalse_whenX0() public {
        assertFalse(FCL_Elliptic_ZZ.ecAff_isOnCurve(0, gy));
    }

    function test_ecAff_isOnCurve_returnsFalse_whenY0() public {
        assertFalse(FCL_Elliptic_ZZ.ecAff_isOnCurve(gx, 0));
    }

    function test_ecAff_isOnCurve_returnsFalse_whenXGreaterThanEqualP(uint256 x) public {
        vm.assume(x >= p);
        assertFalse(FCL_Elliptic_ZZ.ecAff_isOnCurve(x, gy));
    }

    function test_ecAff_isOnCurve_returnsFalse_whenYGreaterThanEqualP(uint256 y) public {
        vm.assume(y >= p);
        assertFalse(FCL_Elliptic_ZZ.ecAff_isOnCurve(gx, y));
    }

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

    function test_Fuzz_SqrtMod(uint256 i_u256a) public {
        vm.assume(i_u256a < FCL_Elliptic_ZZ.p);
        uint256 sqrt = FCL_Elliptic_ZZ.SqrtMod(i_u256a);
        bool flag = (mulmod(sqrt, sqrt, FCL_Elliptic_ZZ.p) == i_u256a);
        if (flag == false) {
            i_u256a = mulmod(i_u256a, p - 1, p); //if a is not a square, -a shall be
            sqrt = FCL_Elliptic_ZZ.SqrtMod(i_u256a);
            flag = (mulmod(sqrt, sqrt, FCL_Elliptic_ZZ.p) == i_u256a);
        }
        assertEq(flag, true);
    }

    function test_Fuzz_ecDecompress(uint256 x, uint256 parity) public {
        vm.assume(x < FCL_Elliptic_ZZ.p);
        vm.assume(x > 0);

        uint256 y = FCL_Elliptic_ZZ.ec_Decompress(x, parity & 1);
        if (y != FCL_Elliptic_ZZ._NOTONCURVE) {
            console.log("on curve y=", y);
            assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(x, y), true);
            assertEq(y & 1, parity & 1);
        }
    }

    //testing Coron Shuffling
    function test_Fuzz_Coronize(uint256 x, uint256 alpha, uint256 beta) public {
        uint256 x2;
        uint256 y2;
        uint256 zz2;
        uint256 zzz2;

        vm.assume(x < FCL_Elliptic_ZZ.p);
        vm.assume(x > 0);
        vm.assume(alpha < FCL_Elliptic_ZZ.p);
        vm.assume(alpha > 0);
        vm.assume(beta < FCL_Elliptic_ZZ.p);
        vm.assume(beta > 0);

        vm.assume(FCL_Elliptic_ZZ.ec_Decompress(x, 0) != FCL_Elliptic_ZZ._NOTONCURVE);
        uint256 y = FCL_Elliptic_ZZ.ec_Decompress(x, 0);
        assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(x, y), true);

        (x2, y2, zz2, zzz2) = FCL_Elliptic_ZZ.ecZZ_Coronize(alpha, x, y, 1, 1);
        (x2, y2, zz2, zzz2) = FCL_Elliptic_ZZ.ecZZ_Coronize(beta, x2, y2, zz2, zzz2);

        (x2, y2) = FCL_Elliptic_ZZ.ecZZ_SetAff(x2, y2, zz2, zzz2);

        assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(x2, y2), true);

        assertEq(x, x2);
    }

    //TODO:fuzzing checking compliancy of ecAdd and ecAddN
    function test_Fuzz_ecAdd(uint256 x, uint256 x2, uint256 rand1, uint256 rand2) public {
        vm.assume(x < FCL_Elliptic_ZZ.p);
        vm.assume(x > 0);
        vm.assume(x2 < FCL_Elliptic_ZZ.p);
        vm.assume(x2 > 0);
        vm.assume(x != x2);
        vm.assume(rand2 < FCL_Elliptic_ZZ.p);
        vm.assume(rand2 > 0);
        vm.assume(rand1 < FCL_Elliptic_ZZ.p);
        vm.assume(rand1 > 0);

        vm.assume(FCL_Elliptic_ZZ.ec_Decompress(x, 0) != FCL_Elliptic_ZZ._NOTONCURVE);
        vm.assume(FCL_Elliptic_ZZ.ec_Decompress(x2, 0) != FCL_Elliptic_ZZ._NOTONCURVE);

        uint256 y = FCL_Elliptic_ZZ.ec_Decompress(x, 0);
        uint256 y2 = FCL_Elliptic_ZZ.ec_Decompress(x2, 0);
        uint256 zz;
        uint256 zzz;

        uint256[4] memory P2;
        uint256[4] memory P3;

        uint256 radd;
        uint256 raddN;

        assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(x, y), true);
        assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(x2, y2), true);

        (x, y, zz, zzz) = FCL_Elliptic_ZZ.ecZZ_Coronize(rand1, x, y, 1, 1);
        (radd, raddN) = FCL_Elliptic_ZZ.ecZZ_SetAff(x, y, zz, zzz);
        assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(radd, raddN), true);

        (P3[0], P3[1], P3[2], P3[3]) = FCL_Elliptic_ZZ.ecZZ_AddN(x, y, zz, zzz, x2, y2);
        (radd,) = FCL_Elliptic_ZZ.ecZZ_SetAff(P3[0], P3[1], P3[2], P3[3]);
        // assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(radd,raddN), true);

        (x2, y2, P2[2], P2[3]) = FCL_Elliptic_ZZ.ecZZ_Coronize(rand2, x2, y2, 1, 1);

        (rand1, rand2) = FCL_Elliptic_ZZ.ecZZ_SetAff(x2, y2, P2[2], P2[3]);
        //      assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(radd,raddN), true);

        (P3[0], P3[1], P3[2], P3[3]) = FCL_Elliptic_ZZ.ecZZ_Add(x, y, zz, zzz, x2, y2, P2[2], P2[3]);

        (raddN,) = FCL_Elliptic_ZZ.ecZZ_SetAff(P3[0], P3[1], P3[2], P3[3]);

        assertEq(radd, raddN);
    }

    function test_Invariant_edge() public {
        //choose Q=2P, then verify duplication is ok
        uint256[4] memory Q;
        (Q[0], Q[1], Q[2], Q[3]) = FCL_Elliptic_ZZ.ecZZ_Dbl(gx, gy, 1, 1);
        uint256[4] memory _4P;
        (_4P[0], _4P[1], _4P[2], _4P[3]) = FCL_Elliptic_ZZ.ecZZ_Dbl(Q[0], Q[1], Q[2], Q[3]);
        uint256 _4P_res1;

        (_4P_res1,) = FCL_Elliptic_ZZ.ecZZ_SetAff(_4P[0], _4P[1], _4P[2], _4P[3]);

        uint256 _4P_res2 = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(gx, gy, 4, 0);
        assertEq(_4P_res1, _4P_res2);

        uint256[2] memory nQ;
        (nQ[0], nQ[1]) = FCL_Elliptic_ZZ.ecZZ_SetAff(Q[0], Q[1], Q[2], Q[3]);
        uint256 _4P_res3 = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(nQ[0], nQ[1], 2, 1);

        assertEq(_4P_res1, _4P_res3);
    }

    //check consistency of ecmulmuladd and Add
    function test_invariant_FCL_Ecmulmuladd() public {
        uint256 ecpoint_Rx = 0;
        uint256 ecpoint_Ry = 0;
        uint256 ecpoint_Rzz = 1;
        uint256 ecpoint_Rzzz = 1;
        uint256 checkpointGasLeft;
        uint256 checkpointGasLeft2;

        //Uncomment the library to test/bench here, todo: replace with switch case
        (ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, ecpoint_Rzzz) = FCL_Elliptic_ZZ.ecZZ_Dbl(gx, gy, 1, 1);
        //(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz )= ECops.twiceProj(gx, gy,1);
        // (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz )=Secp256r1._modifiedJacobianDouble(gx, gy,1);
        //(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz) = Secp256r1_maxrobot._jDouble(gx, gy, 1);

        checkpointGasLeft = gasleft();

        for (uint256 i = 3; i <= _NUM_TEST_ECMULMULADD; i++) {
            (ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, ecpoint_Rzzz) =
                FCL_Elliptic_ZZ.ecZZ_AddN(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, ecpoint_Rzzz, gx, gy);

            //(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz)= ECops.addProj(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz,gx,gy,1);
            // (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz )=Secp256r1._jAdd(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz,gx,gy,1);
            // (ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz) =                Secp256r1_maxrobot._jAdd(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, gx, gy, 1);
        }
        checkpointGasLeft2 = gasleft();
        console.log(
            "Add number test and gas cost :", _NUM_TEST_ECMULMULADD, checkpointGasLeft - checkpointGasLeft2 - 100
        );

        (ecpoint_Rx, ecpoint_Ry) = FCL_Elliptic_ZZ.ecZZ_SetAff(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, ecpoint_Rzzz);
        // (ecpoint_Rx, ecpoint_Ry)=ECops.toAffinePoint(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz);
        //(ecpoint_Rx, ecpoint_Ry) = Secp256r1._affineFromJacobian(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz);

        //generate a small multiple of G using scalar
        assertEq(ecpoint_Rx, FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(0, 0, _NUM_TEST_ECMULMULADD, 0));

        checkpointGasLeft = gasleft();
        for (uint256 i = 3; i <= _NUM_TEST_DBL; i++) {
            (ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, ecpoint_Rzzz) =
                FCL_Elliptic_ZZ.ecZZ_Dbl(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz, ecpoint_Rzzz);

            //   (ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz )= Secp256r1._modifiedJacobianDouble(ecpoint_Rx, ecpoint_Ry,ecpoint_Rzz);
            //(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz) = Secp256r1_maxrobot._jDouble(ecpoint_Rx, ecpoint_Ry, ecpoint_Rzz);
        }
        checkpointGasLeft2 = gasleft();
        console.log(
            "Dbl number test and gas cost :", _NUM_TEST_ECMULMULADD, checkpointGasLeft - checkpointGasLeft2 - 100
        );
    }
}
