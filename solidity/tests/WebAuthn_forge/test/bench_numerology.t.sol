pragma solidity ^0.8.20;

import "forge-std/Test.sol";
import "@external/Numerology.sol";

contract test_Numerology is Test {
    uint256 constant _NUM_TEST_ECMULMULADD = 1000;
    uint256 constant _NUM_TEST_DBL = 100;
    //only for bench purposes

    function test_DblnAdd() public view returns (bool) {
        uint256[3] memory ecpoint_R = [
            0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
            0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8,
            uint256(1)
        ];
        uint256[3] memory ecpoint_G = [
            0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5,
            0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a,
            uint256(1)
        ];

        uint256 checkpointGasLeft;
        uint256 checkpointGasLeft2;

        checkpointGasLeft = gasleft();
        for (uint256 i = 3; i <= _NUM_TEST_ECMULMULADD; i++) {
            (ecpoint_R) = Numerology.addJac(ecpoint_R, ecpoint_G);
        }

        checkpointGasLeft2 = gasleft();
        console.log(
            "Numerology:Add number test and gas cost :",
            _NUM_TEST_ECMULMULADD,
            checkpointGasLeft - checkpointGasLeft2 - 100
        );

        checkpointGasLeft = gasleft();

        for (uint256 i = 3; i <= _NUM_TEST_DBL; i++) {
            (ecpoint_R) = Numerology.doubleJacobian(ecpoint_R);
        }
        checkpointGasLeft2 = gasleft();
        console.log(
            "Numerology:Dbls number test and gas cost :", _NUM_TEST_DBL, checkpointGasLeft - checkpointGasLeft2 - 100
        );

        return true;
    }

    function test_proof_verification() public view returns (bool) {
        // Verifier.deployed().then(function(inst) { return inst.test_proof_verification.estimateGas(); })

        int256[4] memory k_l = [
            int256(-89243190524605339210527649141408088119),
            int256(-53877858828609620138203152946894934485),
            int256(-185204247857117235934281322466442848518),
            int256(-7585701889390054782280085152653861472)
        ];

        uint256[4] memory P_Q = [
            0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
            0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8,
            0xc6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5,
            0x1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a
        ];
        uint256 checkpointGasLeft;
        uint256 checkpointGasLeft2;

        uint256[4] memory wnaf;
        uint256 max_count = 0;
        uint256 count = 0;
        checkpointGasLeft = gasleft();

        for (uint256 j = 0; j < 4; j++) {
            (wnaf[j], count) = Numerology._wnaf(k_l[j]);
            if (count > max_count) {
                max_count = count;
            }
        }

        uint256[3] memory kP_lQ = Numerology._sim_mul_wnaf(wnaf, max_count, P_Q);
        checkpointGasLeft2 = gasleft();
        console.log("Numerology Sec256k1 ecmulmuladd:", checkpointGasLeft - checkpointGasLeft2 - 100);

        uint256[3] memory expected = [
            0x7635e27fba8e1f779dcfdde1b1eacbe0571fbe39ecf6056d29ba4bd3ef5e22f2,
            0x197888e5cec769ac2f1eb65dbcbc0e49c00a8cdf01f8030d8286b68c1933fb18,
            1
        ];

        return Numerology.eqJacobian(kP_lQ, expected);
    }

    function test_proof_verification_hack() public pure returns (bool) {
        // Verifier.deployed().then(function(inst) { return inst.test_proof_verification_hack.estimateGas(); })

        uint256[2] memory k_l = [
            0x1787f38d854231dfec2b27a0f621414d10bfa95970b3e576aed29e1e8287e51e,
            0x60d4afe7ba1967e921f22bb9fbd277f10cfe04810451d7c549ae9099c4a8b677
        ];

        uint256[4] memory P_Q = [
            0xa6ecb3f599964fe04c72e486a8f90172493c21f4185f1ab9a7fe05659480c548,
            0xdf67fd3f4255826c234a5262adc70e14a6d42f13ee55b65e885e666e1dd5d3f5,
            0x3f75f99c8df97f477b700407dd7a956cee2e4915f2df4fa9b11c403935c05dc5,
            0x497495d2e486fd4cff4929265eac91711014b9d6ad1c2220d709e20d466627fd
        ];

        uint256[6] memory expected = [
            0xaddcb45773b26a2f8ac2143627d54f47a12aab533dc1b41b4e791985e9eca496, // kP_x
            0x72da5adb3a30a2cf147d309b0cf58c76b322c82a5edae164e13dbeed6429c41d, // kP_y
            0xf07716879380e987f8b5551a1d989068d0003061088a869a33ceb9848771c6fd, // lQ_x
            0x2447ed4564b75b0f9ff84013aaa37c2ab67a2c621b0edc91a06895f19a93aebb, // lQ_y
            0x9ca8f6ff6a2eb6f62787f70b9f7c4939d1a3890ec87343e4f6716f9f6867eb86, // Rx
            0x290c40f22995dc8b956d2c63ec060d332d082124d638ed618891171db8bc206f // Ry
        ];

        bool sum_is_correct = Numerology.eqJacobian(
            Numerology.addJac([expected[0], expected[1], 1], [expected[2], expected[3], 1]),
            [expected[4], expected[5], 1]
        );
        bool kP_is_correct = Numerology.ecmulVerify(P_Q[0], P_Q[1], k_l[0], expected[0], expected[1]);
        bool lQ_is_correct = Numerology.ecmulVerify(P_Q[2], P_Q[3], k_l[1], expected[2], expected[3]);

        return sum_is_correct && kP_is_correct && lQ_is_correct;
    }

    function test_add_eq_jac() public pure returns (bool) {
        uint256 e0 = 0xaddcb45773b26a2f8ac2143627d54f47a12aab533dc1b41b4e791985e9eca496; // kP_x
        uint256 e1 = 0x72da5adb3a30a2cf147d309b0cf58c76b322c82a5edae164e13dbeed6429c41d; // kP_y
        uint256 e2 = 0xf07716879380e987f8b5551a1d989068d0003061088a869a33ceb9848771c6fd; // lQ_x
        uint256 e3 = 0x2447ed4564b75b0f9ff84013aaa37c2ab67a2c621b0edc91a06895f19a93aebb; // lQ_y
        uint256 e4 = 0x9ca8f6ff6a2eb6f62787f70b9f7c4939d1a3890ec87343e4f6716f9f6867eb86; // Rx
        uint256 e5 = 0x290c40f22995dc8b956d2c63ec060d332d082124d638ed618891171db8bc206f; // Ry

        return Numerology.eqAffineJacobian([e4, e5], Numerology.addAffineJacobian([e0, e1], [e2, e3]));
    }
}
