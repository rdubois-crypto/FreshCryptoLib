// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

// optimized FCL_elliptic: differential against the original library (test/ref), python oracle (test/gen),
// Wycheproof, 8-dim precomputed table variants, gas measurements.

import "forge-std/Test.sol";
import {FCL_Elliptic_ZZ_REF} from "./ref/FCL_elliptic_ref.sol";
import {FCL_Elliptic_ZZ} from "@solidity/FCL_elliptic.sol";
import {FCL_ecdsa} from "@solidity/FCL_ecdsa.sol";

contract RefW_opt {
    uint256 constant n = FCL_Elliptic_ZZ_REF.n;

    function smma(uint256 qx, uint256 qy, uint256 u, uint256 v) external view returns (uint256) {
        return FCL_Elliptic_ZZ_REF.ecZZ_mulmuladd_S_asm(qx, qy, u, v);
    }

    function mma(uint256 qx, uint256 qy, uint256 u, uint256 v) external view returns (uint256, uint256) {
        return FCL_Elliptic_ZZ_REF.ecZZ_mulmuladd(qx, qy, u, v);
    }

    // verbatim FCL_ecdsa.ecdsa_verify wired to the reference library
    function verify(bytes32 m, uint256 r, uint256 s, uint256 qx, uint256 qy) external view returns (bool) {
        if (r == 0 || r >= n || s == 0 || s >= n) return false;
        if (!FCL_Elliptic_ZZ_REF.ecAff_isOnCurve(qx, qy)) return false;
        uint256 sInv = FCL_Elliptic_ZZ_REF.FCL_nModInv(s);
        uint256 x1 = FCL_Elliptic_ZZ_REF.ecZZ_mulmuladd_S_asm(qx, qy, mulmod(uint256(m), sInv, n), mulmod(r, sInv, n));
        return addmod(x1, n - r, n) == 0;
    }
}

contract OptW_opt {
    function smma(uint256 qx, uint256 qy, uint256 u, uint256 v) external view returns (uint256) {
        return FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(qx, qy, u, v);
    }

    function mma(uint256 qx, uint256 qy, uint256 u, uint256 v) external view returns (uint256, uint256) {
        return FCL_Elliptic_ZZ.ecZZ_mulmuladd(qx, qy, u, v);
    }

    function verify(bytes32 m, uint256 r, uint256 s, uint256 qx, uint256 qy) external view returns (bool) {
        return FCL_ecdsa.ecdsa_verify(m, r, s, qx, qy);
    }
}

contract FCL_elliptic_opt_Test is Test {
    RefW_opt refw;
    OptW_opt optw;
    string vec;

    function setUp() public {
        refw = new RefW_opt();
        optw = new OptW_opt();
        vec = vm.readFile("test/vectors_opt/p256.json");
    }

    function arr(string memory k) internal returns (uint256[] memory) {
        return vm.parseJsonUintArray(vec, string.concat(".", k));
    }

    /// optimized mulmuladd == python oracle on every vector, edge cases and rare loop paths included
    function test_mulmuladd_vectors() public {
        uint256[] memory qx = arr("mma_qx");
        uint256[] memory qy = arr("mma_qy");
        uint256[] memory u = arr("mma_u");
        uint256[] memory v = arr("mma_v");
        uint256[] memory x = arr("mma_x");
        uint256[] memory y = arr("mma_y");
        uint256 refBad;
        for (uint256 i = 0; i < qx.length; i++) {
            assertEq(optw.smma(qx[i], qy[i], u[i], v[i]), x[i], "smma");
            (uint256 X, uint256 Y) = optw.mma(qx[i], qy[i], u[i], v[i]);
            assertEq(X, x[i], "mma.x");
            assertEq(Y, y[i], "mma.y");
            (uint256 Xr, uint256 Yr) = refw.mma(qx[i], qy[i], u[i], v[i]);
            if (refw.smma(qx[i], qy[i], u[i], v[i]) != x[i] || Xr != x[i] || Yr != y[i]) refBad++;
        }
        console2.log("mulmuladd vectors:", qx.length, "wrong in the original library:", refBad);
    }

    function test_ecdsa_vectors() public {
        uint256[] memory h = arr("ecdsa_h");
        uint256[] memory r = arr("ecdsa_r");
        uint256[] memory s = arr("ecdsa_s");
        uint256[] memory qx = arr("ecdsa_qx");
        uint256[] memory qy = arr("ecdsa_qy");
        uint256[] memory ok = arr("ecdsa_ok");
        for (uint256 i = 0; i < h.length; i++) {
            assertEq(optw.verify(bytes32(h[i]), r[i], s[i], qx[i], qy[i]), ok[i] == 1);
            assertEq(refw.verify(bytes32(h[i]), r[i], s[i], qx[i], qy[i]), ok[i] == 1);
        }
    }

    function wycheproof(string memory file, bool expected) internal returns (uint256 cnt) {
        string memory d = vm.readFile(file);
        uint256 qx = vm.parseJsonUint(d, ".keyx");
        uint256 qy = vm.parseJsonUint(d, ".keyy");
        uint256 nb = vm.parseJsonUint(d, ".NumberOfTests");
        for (uint256 i = 1; i <= nb; i++) {
            string memory k = vm.toString(i);
            uint256 m = vm.parseJsonUint(d, string.concat(".msg_", k));
            uint256 r = vm.parseJsonUint(d, string.concat(".sigx_", k));
            uint256 s = vm.parseJsonUint(d, string.concat(".sigy_", k));
            assertEq(optw.verify(bytes32(m), r, s, qx, qy), expected, string.concat(file, " #", k));
            assertEq(refw.verify(bytes32(m), r, s, qx, qy), expected, string.concat("ref ", k));
            cnt++;
        }
    }

    function test_wycheproof() public {
        uint256 a = wycheproof("test/vectors_wychproof/vec_sec256r1_valid.json", true);
        uint256 b = wycheproof("test/vectors_wychproof/vec_sec256r1_invalid.json", false);
        console2.log("wycheproof valid / invalid:", a, b);
    }

    /// random differential against the original library (generic keys, where it is correct)
    function testFuzz_diff(uint256 q, uint256 u, uint256 v) public {
        q = bound(q, 1, FCL_Elliptic_ZZ.n - 1);
        u = bound(u, 0, FCL_Elliptic_ZZ.n - 1);
        v = bound(v, 0, FCL_Elliptic_ZZ.n - 1);
        (uint256 qx, uint256 qy) = refw.mma(0, 0, q, 0);
        vm.assume(qx != 0);
        // the original ecZZ_mulmuladd is wrong for Q = -G (and small multiples), covered by the python oracle
        vm.assume(FCL_Elliptic_ZZ.ecZZ_smallMultipleOfG(qx, qy) == 0);
        (uint256 X, uint256 Y) = optw.mma(qx, qy, u, v);
        (uint256 Xr, uint256 Yr) = refw.mma(qx, qy, u, v);
        assertEq(X, Xr);
        assertEq(Y, Yr);
        assertEq(optw.smma(qx, qy, u, v), refw.smma(qx, qy, u, v));
    }

    /// average gas over random valid signatures
    function test_gas() public {
        uint256[] memory h = arr("ecdsa_h");
        uint256[] memory r = arr("ecdsa_r");
        uint256[] memory s = arr("ecdsa_s");
        uint256[] memory qx = arr("ecdsa_qx");
        uint256[] memory qy = arr("ecdsa_qy");
        uint256[6] memory acc;
        uint256 g;
        uint256 cnt = h.length;
        for (uint256 i = 0; i < cnt; i++) {
            g = gasleft();
            refw.verify(bytes32(h[i]), r[i], s[i], qx[i], qy[i]);
            acc[0] += g - gasleft();
            g = gasleft();
            optw.verify(bytes32(h[i]), r[i], s[i], qx[i], qy[i]);
            acc[1] += g - gasleft();
            g = gasleft();
            refw.smma(qx[i], qy[i], h[i] % FCL_Elliptic_ZZ.n, r[i]);
            acc[2] += g - gasleft();
            g = gasleft();
            optw.smma(qx[i], qy[i], h[i] % FCL_Elliptic_ZZ.n, r[i]);
            acc[3] += g - gasleft();
            g = gasleft();
            refw.mma(qx[i], qy[i], h[i] % FCL_Elliptic_ZZ.n, r[i]);
            acc[4] += g - gasleft();
            g = gasleft();
            optw.mma(qx[i], qy[i], h[i] % FCL_Elliptic_ZZ.n, r[i]);
            acc[5] += g - gasleft();
        }
        console2.log("ecdsa_verify          ref / opt:", acc[0] / cnt, acc[1] / cnt);
        console2.log("ecZZ_mulmuladd_S_asm  ref / opt:", acc[2] / cnt, acc[3] / cnt);
        console2.log("ecZZ_mulmuladd        ref / opt:", acc[4] / cnt, acc[5] / cnt);
    }
}

contract S8Ref {
    function ext(uint256 u, uint256 v, address t) external view returns (uint256) {
        return FCL_Elliptic_ZZ_REF.ecZZ_mulmuladd_S8_extcode(u, v, t);
    }

    function hack(uint256 u, uint256 v, uint256 ptr) external view returns (uint256) {
        return FCL_Elliptic_ZZ_REF.ecZZ_mulmuladd_S8_hackmem(u, v, ptr);
    }
}

contract S8Opt {
    function ext(uint256 u, uint256 v, address t) external view returns (uint256) {
        return FCL_Elliptic_ZZ.ecZZ_mulmuladd_S8_extcode(u, v, t);
    }

    function hack(uint256 u, uint256 v, uint256 ptr) external view returns (uint256) {
        return FCL_Elliptic_ZZ.ecZZ_mulmuladd_S8_hackmem(u, v, ptr);
    }
}

contract FCL_elliptic_opt_S8_Test is Test {
    address constant TABLE = address(0x7AB1E);
    S8Ref ref;
    S8Opt opt;
    S8Ref refH;
    S8Opt optH;
    uint256 ptrRef;
    uint256 ptrOpt;
    string vec;

    function setUp() public {
        vec = vm.readFile("test/vectors_opt/p256_s8.json");
        bytes memory tab = vm.parseJsonBytes(vec, ".table");
        vm.etch(TABLE, tab);
        ref = new S8Ref();
        opt = new S8Opt();
        // hackmem: the table is appended to the code of the calling contract
        ptrRef = address(ref).code.length;
        ptrOpt = address(opt).code.length;
        refH = S8Ref(address(0xAAA1));
        optH = S8Opt(address(0xAAA2));
        vm.etch(address(refH), bytes.concat(address(ref).code, tab));
        vm.etch(address(optH), bytes.concat(address(opt).code, tab));
    }

    function test_s8_2G() public {
        // table of Q = 2G: rare paths (accumulator = +/- table entry)
        string memory d = vm.readFile("test/vectors_opt/p256_s8_2G.json");
        bytes memory tab = vm.parseJsonBytes(d, ".table");
        vm.etch(TABLE, tab);
        vm.etch(address(refH), bytes.concat(address(ref).code, tab));
        vm.etch(address(optH), bytes.concat(address(opt).code, tab));
        uint256[] memory u = vm.parseJsonUintArray(d, ".u");
        uint256[] memory v = vm.parseJsonUintArray(d, ".v");
        uint256[] memory x = vm.parseJsonUintArray(d, ".x");
        for (uint256 i = 0; i < u.length; i++) {
            assertEq(opt.ext(u[i], v[i], TABLE), x[i], "opt ext 2G");
            assertEq(optH.hack(u[i], v[i], ptrOpt), x[i], "opt hack 2G");
        }
    }

    function test_s8() public {
        uint256[] memory u = vm.parseJsonUintArray(vec, ".u");
        uint256[] memory v = vm.parseJsonUintArray(vec, ".v");
        uint256[] memory x = vm.parseJsonUintArray(vec, ".x");
        uint256[4] memory acc;
        uint256 cnt;
        for (uint256 i = 0; i < u.length; i++) {
            assertEq(ref.ext(u[i], v[i], TABLE), x[i], "ref ext");
            assertEq(opt.ext(u[i], v[i], TABLE), x[i], "opt ext");
            if (refH.hack(u[i], v[i], ptrRef) != x[i]) console2.log("original S8_hackmem wrong on case", i);
            assertEq(optH.hack(u[i], v[i], ptrOpt), x[i], "opt hack");
            if (i < 16) {
                uint256 g = gasleft();
                ref.ext(u[i], v[i], TABLE);
                acc[0] += g - gasleft();
                g = gasleft();
                opt.ext(u[i], v[i], TABLE);
                acc[1] += g - gasleft();
                g = gasleft();
                refH.hack(u[i], v[i], ptrRef);
                acc[2] += g - gasleft();
                g = gasleft();
                optH.hack(u[i], v[i], ptrOpt);
                acc[3] += g - gasleft();
                cnt++;
            }
        }
        console2.log("S8_extcode ref / opt:", acc[0] / cnt, acc[1] / cnt);
        console2.log("S8_hackmem ref / opt:", acc[2] / cnt, acc[3] / cnt);
    }
}
