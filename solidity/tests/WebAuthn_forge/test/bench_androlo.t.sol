pragma solidity ^0.8.20;

import "forge-std/Test.sol";

import "@external/ECCMath.sol";
import "@external/Secp256k1.sol";

contract test_Androlo is Test {
    // Field size
    uint256 constant pp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;

    // Base point (generator) G
    uint256 constant Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798;
    uint256 constant Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8;

    // Order of G
    uint256 constant nn = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;

    // Maximum value of s
    uint256 constant lowSmax = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B20A0;

    uint256 constant _NUM_TEST_ECMULMULADD = 1000;
    uint256 constant _NUM_TEST_DBL = 100;

    function test_DblnAdd() public view returns (bool) {
        uint256[3] memory ecpoint_R = [uint256(0), uint256(0), uint256(1)];
        uint256[2] memory ecpoint_G = [Gx, Gy];

        uint256 checkpointGasLeft;
        uint256 checkpointGasLeft2;

        checkpointGasLeft = gasleft();
        for (uint256 i = 3; i <= _NUM_TEST_ECMULMULADD; i++) {
            (ecpoint_R) = Secp256k1._addMixed(ecpoint_R, ecpoint_G);
        }

        checkpointGasLeft2 = gasleft();
        console.log(
            "Androlo:Add number test and gas cost :",
            _NUM_TEST_ECMULMULADD,
            checkpointGasLeft - checkpointGasLeft2 - 100
        );
        checkpointGasLeft = gasleft();

        for (uint256 i = 3; i <= _NUM_TEST_ECMULMULADD; i++) {
            (ecpoint_R) = Secp256k1._double(ecpoint_R);
        }
        checkpointGasLeft2 = gasleft();
        console.log(
            "Androlo:Dbls number test and gas cost :", _NUM_TEST_DBL, checkpointGasLeft - checkpointGasLeft2 - 100
        );

        return true;
    }

    /// @dev See Curve.validateSignature
    function bench_Signature(uint256 e, uint256[2] memory rs, uint256[2] memory Q) internal pure returns (bool) {
        uint256 n = nn;
        uint256 p = pp;
        if (rs[0] == 0 || rs[0] >= n || rs[1] == 0 || rs[1] > lowSmax) {
            return false;
        }
        if (!Secp256k1.isPubKey(Q)) {
            return false;
        }

        uint256 sInv = ECCMath.invmod(rs[1], n);

        uint256[3] memory u1G = Secp256k1._mul(mulmod(e, sInv, n), [Gx, Gy]);

        uint256[3] memory u2Q = Secp256k1._mul(mulmod(rs[0], sInv, n), Q);
        uint256[3] memory P = Secp256k1._add(u1G, u2Q);

        if (P[2] == 0) {
            return false;
        }

        uint256 Px = ECCMath.invmod(P[2], p); // need Px/Pz^2
        Px = mulmod(P[0], mulmod(Px, Px, p), p);
        //return Px % n == rs[0];
        return true; //comparizon set to true: the aim is only to bench gas cost here
    }

    function test_sig() public view returns (bool) {
        uint256 checkpointGasLeft;
        uint256 checkpointGasLeft2;

        uint256[2] memory dummy_rs = [
            0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364140,
            0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B2070
        ];
        uint256 dummy_e = 0xEBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364140;
        bool res;
        uint256[2] memory Q = [
            0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5,
            0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a
        ];
        checkpointGasLeft = gasleft();
        res = bench_Signature(dummy_e, dummy_rs, Q);
        checkpointGasLeft2 = gasleft();
        console.log("Androlo Sec256k1 ecdsa:", checkpointGasLeft - checkpointGasLeft2 - 100);

        return res;
    }
}
