// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.13;

import "forge-std/Test.sol";
import "../src/FCL_elliptic.sol";
import "../src/FCL_Webauthn.sol";

//external implementation to bench
import "../src/ECops.sol";
import "../src/Secp256r1.sol";
import "../src/Secp256r1_maxrobot.sol";

// library elliptic solidity from orbs network
contract wrap_ecdsa_orbs {
    uint256 constant gx = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296;
    uint256 constant gy = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5;
    //curve order (number of points)
    uint256 constant n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;

    function wrap_ecdsa_core(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q) public returns (bool) {
        if (rs[0] == 0 || rs[0] >= n || rs[1] == 0 || rs[1] >= n) {
            return false;
        }

        if (!FCL_Elliptic_ZZ.ecAff_isOnCurve(Q[0], Q[1])) {
            return false;
        }

        uint256 sInv = FCL_Elliptic_ZZ.FCL_nModInv(rs[1]);

        uint256 scalar_u = mulmod(uint256(message), sInv, n);
        uint256 scalar_v = mulmod(rs[0], sInv, n);
        uint256[2] memory P1;
        uint256[2] memory P2;
        (P1[0], P1[1]) = ECops.multiplyScalar(gx, gy, scalar_u);
        (P2[0], P2[1]) = ECops.multiplyScalar(Q[0], Q[1], scalar_v);

        uint256 x1;
        (x1,) = ECops.add(P1[0], P1[1], P2[0], P2[1]);
        assembly {
            x1 := addmod(x1, sub(n, calldataload(rs)), n)
        }
        //return true;
        return x1 == 0;
    }
}

// library from obvioustech
contract wrap_ecdsa_obvious {
    function wrap_ecdsa_core(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q)
        public
        view
        returns (bool)
    {
        PassKeyId memory pass = PassKeyId(Q[0], Q[1], "unused");
        return Secp256r1.Verify(pass, rs[0], rs[1], uint256(message));
    }
}

// library from maxrobot
contract wrap_ecdsa_maxrobot {
    function wrap_ecdsa_core(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q)
        public
        view
        returns (bool)
    {
        return Secp256r1_maxrobot.Verify(Q[0], Q[1], rs, uint256(message));
    }
}

// library FreshCryptoLib without precomputations
contract Wrap_ecdsa_FCL {
    function wrap_ecdsa_core(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q) public returns (bool) {
        return FCL_Elliptic_ZZ.ecdsa_verify(message, rs, Q);
    }

    constructor() {}
}

// library FreshCryptoLib with precomputations
contract Wrap_ecdsa_precal {
    address precomputations;

    function wrap_ecdsa_core(bytes32 message, uint256[2] calldata rs) public returns (bool) {
        return FCL_Elliptic_ZZ.ecdsa_precomputed_verify(message, rs, precomputations);
    }

    constructor(address bytecode) {
        precomputations = bytecode;
    }
}

contract EcdsaTest is Test {
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

    /*//////////////////////////////////////////////////////////////
                                 TESTS
    //////////////////////////////////////////////////////////////*/

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

    function test_Invariant_edge() public {
        console.logString(vm.projectRoot());
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

    //testing Wychproof vectors
    function test_Invariant_ecZZ_mulmuladd_S_asm() public {
        string memory deployData;
        uint256[2] memory key;
        uint256 numtests;
        (key, deployData, numtests) = wychproof_keyload();
        uint256 checkpointGasLeft;
        uint256 checkpointGasLeft2;

        bool res = FCL_Elliptic_ZZ.ecAff_isOnCurve(key[0], key[1]);
        assertEq(res, true);

        uint256[2] memory rs;
        string memory title;
        string memory snum = "1";
        for (uint256 i = 1; i <= numtests; i++) {
            snum = vm.toString(i);
            uint256 message;
            (rs, message, title) = wychproof_vecload(deployData, snum);

            vm.prank(vm.addr(5));

            checkpointGasLeft = gasleft();
            wrap_ecdsa_maxrobot wrap = new wrap_ecdsa_maxrobot();

            checkpointGasLeft2 = gasleft();
            console.log("deployment no prec:", checkpointGasLeft - checkpointGasLeft2 - 100);

            checkpointGasLeft = gasleft();
            Wrap_ecdsa_precal wrap2 = new Wrap_ecdsa_precal(address(uint160(_prec_address)));
            checkpointGasLeft2 = gasleft();
            console.log("deployment with prec, (no table cost):", checkpointGasLeft - checkpointGasLeft2 - 100);

            console.log("message:", message);
            console.log("sigx:", rs[0]);
            console.log("sigy:", rs[1]);

            checkpointGasLeft = gasleft();
            res = wrap.wrap_ecdsa_core(bytes32(message), rs, key);
            checkpointGasLeft2 = gasleft();
            console.log("signature verif no prec:", checkpointGasLeft - checkpointGasLeft2 - 100);

            assertEq(res, true);
            checkpointGasLeft = gasleft();
            res = wrap2.wrap_ecdsa_core(bytes32(message), rs);
            checkpointGasLeft2 = gasleft();
            console.log("signature verif with prec:", checkpointGasLeft - checkpointGasLeft2 - 100);

            // ensure both implementations return the same result
            assertEq(res, true);
        }
    }

    /*//////////////////////////////////////////////////////////////
                           INTERNAL FUNCTIONS
    //////////////////////////////////////////////////////////////*/

    function wychproof_keyload() public returns (uint256[2] memory key, string memory deployData, uint256 numtests) {
        deployData = vm.readFile("test/vectors_wychproof/vec_sec256r1_valid.json");

        uint256 wx = vm.parseJsonUint(deployData, ".NumberOfTests");
        console.log("NumberOfTests:", wx);
        key[0] = vm.parseJsonUint(deployData, ".keyx");
        console.log("key_x:", key[0]);
        key[1] = vm.parseJsonUint(deployData, ".keyy");
        console.log("key_y:", key[1]);
        bool res = FCL_Elliptic_ZZ.ecAff_isOnCurve(key[0], key[1]);
        assertEq(res, true);

        bytes memory precompute = precompute_shamir_table(key[0], key[1]);
        verify_precompute(precompute);

        console.log("Is key on curve:", res);

        return (key, deployData, wx);
    }

    //load a single test vector
    function wychproof_vecload(string memory deployData, string memory snum)
        public
        returns (uint256[2] memory rs, uint256 message, string memory title)
    {
        title = string(vm.parseJson(deployData, string.concat(".test_", snum)));

        console.log("\n test:", snum, title);
        console.log("\n ||");

        rs[0] = vm.parseJsonUint(deployData, string.concat(".sigx_", snum));
        rs[1] = vm.parseJsonUint(deployData, string.concat(".sigy_", snum));
        message = vm.parseJsonUint(deployData, string.concat(".msg_", snum));
    }

    function precompute_shamir_table(uint256 C0, uint256 C1) private returns (bytes memory precompute) {
        vm.setEnv("C0", vm.toString(C0));
        vm.setEnv("C1", vm.toString(C1));

        // Precompute a 8 dimensional table for Shamir's trick from C0 and C1
        // and return the table as a bytes
        string[] memory inputs = new string[](2);
        inputs[0] = "sage";
        inputs[1] = "test/sage/shamir-precomputation.sage";
        precompute = vm.ffi(inputs);
    }

    function verify_precompute(bytes memory prec) private returns (bool) {
        // address of the precomputations bytecode contract
        address a_prec = address(uint160(_prec_address));
        // set the precomputed points as the bytecode of the contract
        vm.etch(a_prec, prec);

        //pointer to an elliptic point
        uint256[2] memory px;

        // check the precomputations are correct, all point are on curve P256
        for (uint256 i = 1; i < 256; i++) {
            uint256 offset = 64 * i;

            assembly {
                extcodecopy(a_prec, px, offset, 64)
            }

            assertEq(FCL_Elliptic_ZZ.ecAff_isOnCurve(px[0], px[1]), true);
        }

        return true;
    }
}
