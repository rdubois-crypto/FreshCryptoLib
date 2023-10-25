//********************************************************************************************/
//  ___           _       ___               _         _    _ _
// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__
// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
//                                |__/|_|
///* Copyright (C) 2023 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project
///* License: This software is licensed under MIT License
///* This Code may be reused including license and copyright notice.
///* See LICENSE file at the root folder of the project.
///* FILE: FCL_ecdsa_utils.t.sol
///*
///*
///* DESCRIPTION: test file for ecdsa signature protocol utilitary functions
///*
//**************************************************************************************/
// SPDX-License-Identifier: MIT

pragma solidity >=0.8.19 <0.9.0;

import "forge-std/Test.sol";
import "@solidity/FCL_elliptic.sol";

import "@solidity/FCL_ecdsa.sol";
import "@solidity/FCL_Webauthn.sol";

//external implementation to bench
import "@external/Secp256r1.sol";
import "@external/Secp256r1_maxrobot.sol";
import "@external/ECops.sol";

contract EcdsaUtilsTest is Test {
    //Fuzz over public key verification using Decompress/Derivation
    function test_Fuzz_SigVerif(uint256 k, uint256 kpriv, uint256 message) public {
        vm.assume(k < FCL_Elliptic_ZZ.n);
        vm.assume(k > 1);
        vm.assume(kpriv < FCL_Elliptic_ZZ.n);
        vm.assume(kpriv > 1);

        vm.assume(message < FCL_Elliptic_ZZ.n);
        vm.assume(message > 1);

        k = FCL_Elliptic_ZZ.n - k; //ensure high hamming weight of fuzzing vectors
        // kpriv=115792089210356248762697446949407573529996955224135760342422259061068512044365;

        uint256 xpub = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(0, 0, kpriv, 0); //deriv public key
        uint256 ypub = FCL_Elliptic_ZZ.ec_Decompress(xpub, 0);
        uint256 r;
        uint256 s;
        assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(xpub, ypub), true);

        (r, s) = FCL_ecdsa_utils.ecdsa_sign(bytes32(message), k, kpriv);

        bool res1 = FCL_ecdsa.ecdsa_verify(bytes32(message), r, s, xpub, ypub);
        bool res2 = FCL_ecdsa.ecdsa_verify(bytes32(message), r, s, xpub, FCL_Elliptic_ZZ.p - ypub);
        bool res = res1 || res2;

        assertEq(res, true);
    }

    function test_derivKpub(uint256 kpriv) public view {
        vm.assume(kpriv < FCL_Elliptic_ZZ.n);
        vm.assume(kpriv > 1);
        uint256 x;
        uint256 y;

        (x, y) = FCL_ecdsa_utils.ecdsa_derivKpub(kpriv);
        if (FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(x, y, kpriv, FCL_Elliptic_ZZ.n - 1) != 0) {
            revert();
        }
        console.log(FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(x, FCL_Elliptic_ZZ.p - y, kpriv, kpriv));
    }

    function test_Fuzz_Precompute() public {
        uint256 kpriv = 3;
        vm.assume(kpriv < FCL_Elliptic_ZZ.n);
        vm.assume(kpriv > 2);
        uint256 x;
        uint256 y;

        (x, y) = FCL_ecdsa_utils.ecdsa_derivKpub(kpriv);

        uint256[2][256] memory Precs = FCL_ecdsa_utils.Precalc_8dim(x, y);

        bytes memory prec = abi.encodePacked(Precs);
        address a_prec;
        assembly {
            a_prec := create2(0, add(prec, 0x20), mload(prec), 0)
        }

        // console.logBytes(a_prec.code);
        vm.etch(address(uint160(a_prec)), prec);

        uint256 x1 = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(x, y, 66, 12);
        uint256 x2 = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S8_extcode(66, 12, a_prec);

        assertEq(x1, x2);
    }
}
