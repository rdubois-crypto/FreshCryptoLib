// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.13;

import "forge-std/Test.sol";
import "../src/FCL_elliptic.sol";

import "../src/FCL_Webauthn.sol";

contract Wrap_ecdsa {
    function wrap_ecdsa_core(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q) public returns (bool) {
        return FCL_Elliptic_ZZ.ecdsa_verify(message, rs, Q);
    }

    constructor() {}
}

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

    function uintToBytes(uint256 v) public pure returns (bytes32 ret) {
        if (v == 0) {
            ret = "0";
        } else {
            while (v > 0) {
                ret = bytes32(uint256(ret) / (2 ** 8));
                ret |= bytes32(((v % 10) + 48) * 2 ** (8 * 31));
                v /= 10;
            }
        }
        return ret;
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

    function write_precalcsage(uint256 C0, uint256 C1) public {
        string memory line = "cd ../../sage/FCL_ecdsa_precompute;rm -f *json;sage -c 'C0=";
        line = string.concat(line, vm.toString(C0));
        line = string.concat(line, ";C1=");
        line = string.concat(line, vm.toString(C1));
        line = string.concat(
            line,
            ";load(\"FCL_ecdsa_precompute.sage\")' ;cp fcl_ecdsa_precbytecode.json ../../tests/WebAuthn_forge/test/vectors_sage/fcl_ecdsa_precbytecode.json;"
        );

        vm.writeLine("scriptz.sh", line);

        string[] memory inputs = new string[](2);

        inputs[0] = "bash";
        inputs[1] = "scriptz.sh";
        vm.ffi(inputs);
        vm.removeFile("scriptz.sh");
        console.log("precalc done");
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
        write_precalcsage(key[0], key[1]);
        load_precalc();

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

    //testing Wychproof vectors
    function test_Invariant_ecZZ_mulmuladd_S_asm() public {
        string memory deployData;
        uint256[2] memory key;
        uint256 numtests;
        (key, deployData, numtests) = wychproof_keyload();

        bool res = FCL_Elliptic_ZZ.ecAff_isOnCurve(key[0], key[1]);
        bool res2;
        assertEq(res, true);

        uint256[2] memory rs;
        string memory title;
        string memory snum = "1";
        for (uint256 i = 1; i <= numtests; i++) {
            snum = vm.toString(i);
            uint256 message;
            (rs, message, title) = wychproof_vecload(deployData, snum);

            vm.prank(vm.addr(5));
            Wrap_ecdsa wrap = new Wrap_ecdsa();
            Wrap_ecdsa_precal wrap2 = new Wrap_ecdsa_precal(address(uint160(_prec_address)));

            console.log("message:", message);
            console.log("sigx:", rs[0]);
            console.log("sigy:", rs[1]);

            res = wrap.wrap_ecdsa_core(bytes32(message), rs, key);
            res2 = wrap2.wrap_ecdsa_core(bytes32(message), rs);
            // ensure both implementations return the same result
            assertEq(res, res2);
        }
    }

    function load_precalc() public returns (bool) {
        string memory deployData = vm.readFile("test/vectors_sage/fcl_ecdsa_precbytecode.json");
        bytes memory prec = abi.decode(vm.parseJson(deployData, ".Bytecode"), (bytes));
        address a_prec; //address of the precomputations bytecode contract
        a_prec = address(uint160(_prec_address));
        vm.etch(a_prec, prec);

        uint256[2] memory px; //pointer to an elliptic point

        //check precomputations are correct, all point on curve P256
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
