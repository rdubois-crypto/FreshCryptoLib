// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import {BaseScript} from "./BaseScript.sol";
import {FCL_Elliptic_ZZ} from "../../../src/FCL_elliptic.sol";
import {FCL_ecdsa} from "../../../src/FCL_elliptic.sol";

/// @notice Wrap the FCL_Elliptic library in a contract to be able to deploy it
contract LibraryWrapper {
    function ecdsa_verify(bytes32 message, uint256 r, uint256 s, uint256 Qx, uint256 Qy) external view returns (bool) {
        return FCL_ecdsa.ecdsa_verify(message, r, s, Qx, Qy);
    }

    function ecdsa_verify(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q)
        external
        view
        returns (bool)
    {
        return FCL_ecdsa.ecdsa_verify(message, rs, Q);
    }

    function ecdsa_precomputed_verify(bytes32 message, uint256[2] calldata rs, address Shamir8)
        external
        view
        returns (bool)
    {
        return FCL_Elliptic_ZZ.ecdsa_precomputed_verify(message, rs, Shamir8);
    }

    function ecdsa_sign(bytes32 message, uint256 k, uint256 kpriv) external view returns (uint256 r, uint256 s) {
        return FCL_ecdsa.ecdsa_sign(message, k, kpriv);
    }
}

/// @notice This script deploys the FCL_Elliptic library and the wrapper contract
contract MyScript is BaseScript {
    function run() external broadcast returns (address addressOfLibrary) {
        // deploy the library contract and return the address
        addressOfLibrary = address(new LibraryWrapper{salt:0}());
    }
}

/*
    In the tests/WebAuthn_forge/script directory, run the following command to deploy the library:

    ℹ️ RUN THIS SCRIPT USING A LEDGER:
    forge script DeployElliptic.s.sol:MyScript --rpc-url <RPC_URL> --ledger --sender <ACCOUNT_ADDRESS> \
    [--broadcast]

    ℹ️ RUN THIS SCRIPT WITH AN ARBITRARY PRIVATE KEY (NOT RECOMMENDED):
    PRIVATE_KEY=<PRIVATE_KEY> forge script DeployElliptic.s.sol:MyScript --rpc-url <RPC_URL> [--broadcast]

    ℹ️ RUN THIS SCRIPT ON ANVIL IN DEFAULT MODE:
    forge script DeployElliptic.s.sol:MyScript --rpc-url http://127.0.0.1:8545 --broadcast --sender \
    0xf39Fd6e51aad88F6F4ce6aB8827279cffFb92266 --mnemonics "test test test test test test test test test test test junk"

    ℹ️ CALL THE LIBRARY ONCE DEPLOYED:
    cast call <CONTRACT_ADDRESS> verify(bytes32,uint256,uint256,uint256,uint256)" <MESSAGE> <R> <S> <QX> <QY>

    example:
        cast call 0xe7f1725E7734CE288F8367e1Bb143E90bb3F0512 \
        "ecdsa_verify(bytes32,uint256,uint256,uint256,uint256)" \
        0xbb5a52f42f9c9261ed4361f59422a1e30036e7c32b270c8807a419feca605023 \
        19738613187745101558623338726804762177711919211234071563652772152683725073944 \
        34753961278895633991577816754222591531863837041401341770838584739693604822390 \
        18614955573315897657680976650685450080931919913269223958732452353593824192568 \
        90223116347859880166570198725387569567414254547569925327988539833150573990206
*/
