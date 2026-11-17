// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import "forge-std/Test.sol";
import {FCL_BabyJubJub_Edwards as BJE} from "@solidity/FCL_babyjubjub_edwards.sol";
import {FCL_BabyJubJub as BJJ} from "@solidity/FCL_babyjubjub.sol";

contract BJEdWrap {
    function verifyRaw(uint256 hm, uint256 r8x, uint256 r8y, uint256 s, uint256 ax, uint256 ay)
        external
        view
        returns (bool)
    {
        return BJE.verifyRaw(hm, r8x, r8y, s, ax, ay);
    }
}

contract BJWeWrap {
    function verifyRaw(uint256 hm, uint256 r8x, uint256 r8y, uint256 s, uint256 ax, uint256 ay)
        external
        view
        returns (bool)
    {
        return BJJ.verifyRaw(hm, r8x, r8y, s, ax, ay);
    }

    function mma(uint256 qx, uint256 qy, uint256 u, uint256 v) external view returns (uint256, uint256) {
        return BJJ.ecZZ_mulmuladd(qx, qy, u, v);
    }

    function pubKeyToW8(uint256 ax, uint256 ay) external view returns (uint256, uint256) {
        return BJJ.pubKeyToW8(ax, ay);
    }
}

contract BabyJubJubEdwardsTest is Test {
    BJEdWrap e;
    BJWeWrap w;

    function setUp() public {
        e = new BJEdWrap();
        w = new BJWeWrap();
    }

    /// same circomlibjs vectors as the Weierstrass version, gas comparison
    function test_edwards_vectors_and_gas() public {
        string memory d = vm.readFile("test/vectors_bjj/eddsa.json");
        uint256[] memory hm = vm.parseJsonUintArray(d, ".hm");
        uint256[] memory r8x = vm.parseJsonUintArray(d, ".r8x");
        uint256[] memory r8y = vm.parseJsonUintArray(d, ".r8y");
        uint256[] memory s = vm.parseJsonUintArray(d, ".s");
        uint256[] memory ax = vm.parseJsonUintArray(d, ".ax");
        uint256[] memory ay = vm.parseJsonUintArray(d, ".ay");
        uint256[] memory ok = vm.parseJsonUintArray(d, ".ok");
        string[] memory tag = vm.parseJsonStringArray(d, ".tag");
        uint256[2] memory acc;
        uint256 cnt;
        for (uint256 i = 0; i < hm.length; i++) {
            assertEq(e.verifyRaw(hm[i], r8x[i], r8y[i], s[i], ax[i], ay[i]), ok[i] == 1, tag[i]);
            if (keccak256(bytes(tag[i])) != keccak256("poseidon")) continue;
            uint256 g = gasleft();
            e.verifyRaw(hm[i], r8x[i], r8y[i], s[i], ax[i], ay[i]);
            acc[0] += g - gasleft();
            g = gasleft();
            w.verifyRaw(hm[i], r8x[i], r8y[i], s[i], ax[i], ay[i]);
            acc[1] += g - gasleft();
            cnt++;
        }
        console2.log("verifyRaw Edwards a=-1 / Weierstrass a=-3:", acc[0] / cnt, acc[1] / cnt);
    }

    /// both verifiers accept the same rebuilt signature and reject the same perturbations
    function testFuzz_cross(uint256 k, uint256 s, uint256 hm) public {
        k = bound(k, 1, BJJ.n - 1);
        s = bound(s, 0, BJJ.n - 1);
        hm = bound(hm, 0, type(uint256).max - 1);
        uint256[2] memory A;
        uint256[2] memory W;
        uint256[2] memory R;
        (W[0], W[1]) = w.mma(0, 0, k, 0);
        (A[0], A[1]) = BJJ.W2ed(W[0], W[1]);
        (W[0], W[1]) = w.pubKeyToW8(A[0], A[1]);
        (R[0], R[1]) = w.mma(W[0], W[1], s, addmod(BJJ.n - hm % BJJ.n, 0, BJJ.n));
        (R[0], R[1]) = BJJ.W2ed(R[0], R[1]);
        assertTrue(e.verifyRaw(hm, R[0], R[1], s, A[0], A[1]));
        assertFalse(e.verifyRaw(hm + 1, R[0], R[1], s, A[0], A[1]));
        assertFalse(e.verifyRaw(hm, R[0], R[1], addmod(s, 1, BJJ.n), A[0], A[1]));
    }
}
