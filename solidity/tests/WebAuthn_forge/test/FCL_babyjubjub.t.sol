// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import "forge-std/Test.sol";
import {FCL_BabyJubJub as BJJ} from "@solidity/FCL_babyjubjub.sol";

/// external wrapper: gas measurement and memory isolation, as a calling contract would do
contract BJJWrap {
    function verifyRaw(uint256 hm, uint256 r8x, uint256 r8y, uint256 s, uint256 ax, uint256 ay)
        external
        view
        returns (bool)
    {
        return BJJ.verifyRaw(hm, r8x, r8y, s, ax, ay);
    }

    function verifyRawW(uint256 hm, uint256 r8x, uint256 r8y, uint256 s, uint256 wx, uint256 wy)
        external
        view
        returns (bool)
    {
        return BJJ.verifyRawW(hm, r8x, r8y, s, wx, wy);
    }

    function pubKeyToW8(uint256 ax, uint256 ay) external view returns (uint256, uint256) {
        return BJJ.pubKeyToW8(ax, ay);
    }

    function mma(uint256 qx, uint256 qy, uint256 u, uint256 v) external view returns (uint256, uint256) {
        return BJJ.ecZZ_mulmuladd(qx, qy, u, v);
    }
}

struct Vec {
    uint256[] hm;
    uint256[] r8x;
    uint256[] r8y;
    uint256[] s;
    uint256[] ax;
    uint256[] ay;
    uint256[] ok;
    string[] tag;
}

contract BabyJubJubTest is Test {
    BJJWrap w;
    string eddsa;
    string weier;

    function setUp() public {
        w = new BJJWrap();
        eddsa = vm.readFile("test/vectors_bjj/eddsa.json");
        weier = vm.readFile("test/vectors_bjj/weier.json");
    }

    function arr(string memory d, string memory k) internal returns (uint256[] memory) {
        return vm.parseJsonUintArray(d, k);
    }

    /// the isomorphism maps Base8 to the engine generator, both curve equations hold
    function test_isomorphism() public {
        (uint256 X, uint256 Y) = BJJ.ed2W(BJJ.B8X, BJJ.B8Y);
        assertEq(X, BJJ.gx);
        assertEq(Y, BJJ.gy);
        assertTrue(BJJ.ecAff_isOnCurve(X, Y));
        assertTrue(BJJ.ed_isOnCurve(BJJ.B8X, BJJ.B8Y));
        (uint256 x, uint256 y) = BJJ.W2ed(X, Y);
        assertEq(x, BJJ.B8X);
        assertEq(y, BJJ.B8Y);
        // n.G = O
        (X, Y) = w.mma(0, 0, BJJ.n, 0);
        assertEq(X, 0);
        assertEq(Y, 0);
    }

    /// pubKeyToW8 == phi(8A) computed in python (test/gen/bjj.py)
    function test_pubKeyToW8() public {
        uint256[] memory ax = arr(weier, ".w8.ax");
        uint256[] memory ay = arr(weier, ".w8.ay");
        uint256[] memory wx = arr(weier, ".w8.wx");
        uint256[] memory wy = arr(weier, ".w8.wy");
        for (uint256 i = 0; i < ax.length; i++) {
            (uint256 X, uint256 Y) = w.pubKeyToW8(ax[i], ay[i]);
            assertEq(X, wx[i]);
            assertEq(Y, wy[i]);
            (uint256 ex, uint256 ey) = BJJ.W2ed(X, Y);
            assertTrue(BJJ.ed_isOnCurve(ex, ey));
        }
    }

    /// uG+vQ on the isomorphic curve: random, degenerate (Q = +/-c.G, point at infinity) and rare loop paths
    function test_mulmuladd() public {
        uint256[] memory qx = arr(weier, ".mma.qx");
        uint256[] memory qy = arr(weier, ".mma.qy");
        uint256[] memory u = arr(weier, ".mma.u");
        uint256[] memory v = arr(weier, ".mma.v");
        uint256[] memory x = arr(weier, ".mma.x");
        uint256[] memory y = arr(weier, ".mma.y");
        for (uint256 i = 0; i < qx.length; i++) {
            (uint256 X, uint256 Y) = w.mma(qx[i], qy[i], u[i], v[i]);
            assertEq(X, x[i], vm.toString(i));
            assertEq(Y, y[i], vm.toString(i));
        }
    }

    /// circomlibjs signatures (Poseidon, MiMC7, MiMCSponge, Pedersen), invalid and small order cases
    function load() internal returns (Vec memory V) {
        V.hm = arr(eddsa, ".hm");
        V.r8x = arr(eddsa, ".r8x");
        V.r8y = arr(eddsa, ".r8y");
        V.s = arr(eddsa, ".s");
        V.ax = arr(eddsa, ".ax");
        V.ay = arr(eddsa, ".ay");
        V.ok = arr(eddsa, ".ok");
        V.tag = vm.parseJsonStringArray(eddsa, ".tag");
    }

    function check(Vec memory V, uint256 i) internal {
        bool expected = V.ok[i] == 1;
        assertEq(w.verifyRaw(V.hm[i], V.r8x[i], V.r8y[i], V.s[i], V.ax[i], V.ay[i]), expected, V.tag[i]);
        if (BJJ.ed_isOnCurve(V.ax[i], V.ay[i])) {
            (uint256 wx, uint256 wy) = w.pubKeyToW8(V.ax[i], V.ay[i]);
            assertEq(w.verifyRawW(V.hm[i], V.r8x[i], V.r8y[i], V.s[i], wx, wy), expected, V.tag[i]);
        }
    }

    function test_circomlib_vectors() public {
        Vec memory V = load();
        uint256 valid;
        for (uint256 i = 0; i < V.hm.length; i++) {
            check(V, i);
            if (V.ok[i] == 1) valid++;
        }
        console2.log("circomlib vectors:", V.hm.length, "valid:", valid);
    }

    /// average gas over the valid Poseidon signatures
    function test_gas() public {
        Vec memory V = load();
        uint256[3] memory acc;
        uint256 cnt;
        for (uint256 i = 0; i < V.hm.length; i++) {
            if (keccak256(bytes(V.tag[i])) != keccak256("poseidon")) continue;
            (uint256 wx, uint256 wy) = w.pubKeyToW8(V.ax[i], V.ay[i]);
            uint256 g = gasleft();
            w.verifyRaw(V.hm[i], V.r8x[i], V.r8y[i], V.s[i], V.ax[i], V.ay[i]);
            acc[0] += g - gasleft();
            g = gasleft();
            w.verifyRawW(V.hm[i], V.r8x[i], V.r8y[i], V.s[i], wx, wy);
            acc[1] += g - gasleft();
            g = gasleft();
            w.pubKeyToW8(V.ax[i], V.ay[i]);
            acc[2] += g - gasleft();
            cnt++;
        }
        console2.log("verifyRaw  (Edwards key)      :", acc[0] / cnt);
        console2.log("verifyRawW (precomputed W8)  :", acc[1] / cnt);
        console2.log("pubKeyToW8                   :", acc[2] / cnt);
    }

    /// random consistency: R8 = S.B8 - 8hm.A must verify, hm+1 must not
    function testFuzz_consistency(uint256 k, uint256 s, uint256 hm) public {
        k = bound(k, 1, BJJ.n - 1);
        s = bound(s, 0, BJJ.n - 1);
        hm = bound(hm, 0, type(uint256).max - 1);
        uint256[2] memory A;
        uint256[2] memory W;
        uint256[2] memory R;
        // A = k.B8 (Edwards), W8 = 8A
        (W[0], W[1]) = w.mma(0, 0, k, 0);
        (A[0], A[1]) = BJJ.W2ed(W[0], W[1]);
        (W[0], W[1]) = w.pubKeyToW8(A[0], A[1]);
        // R8 = s.B8 - hm.W8 computed by the engine, mapped back to Edwards
        (R[0], R[1]) = w.mma(W[0], W[1], s, addmod(BJJ.n - hm % BJJ.n, 0, BJJ.n));
        (R[0], R[1]) = BJJ.W2ed(R[0], R[1]);
        assertTrue(w.verifyRaw(hm, R[0], R[1], s, A[0], A[1]));
        assertFalse(w.verifyRaw(hm + 1, R[0], R[1], s, A[0], A[1]));
    }
}
