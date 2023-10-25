// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import {BaseScript} from "./BaseScript.sol";
import {FCL_ecdsa} from "@solidity/FCL_ecdsa.sol";
import {FCL_ecdsa_utils} from "@solidity/FCL_ecdsa_utils.sol";
import {FCL_Elliptic_ZZ} from "@solidity/FCL_elliptic.sol";

/// @notice Wrap the FCL_ecdsa library in a contract to be able to deploy it
contract FCL_ecdsa_wrapper {
    function ecdsa_verify(bytes32 message, uint256 r, uint256 s, uint256 Qx, uint256 Qy) external view returns (bool) {
        return FCL_ecdsa.ecdsa_verify(message, r, s, Qx, Qy);
    }
}

contract FCL_all_wrapper {
    function ecdsa_verify(bytes32 message, uint256 r, uint256 s, uint256 Qx, uint256 Qy) external view returns (bool) {
        return FCL_ecdsa.ecdsa_verify(message, r, s, Qx, Qy);
    }

    function ecdsa_verify(bytes32 message, uint256[2] calldata rs, uint256[2] calldata Q)
        external
        view
        returns (bool)
    {
        return FCL_ecdsa_utils.ecdsa_verify(message, rs, Q);
    }

    function ecdsa_precomputed_verify(bytes32 message, uint256 r, uint256 s, address Shamir8)
        external
        view
        returns (bool)
    {
        return FCL_ecdsa.ecdsa_precomputed_verify(message, r, s, Shamir8);
    }

    function ecdsa_sign(bytes32 message, uint256 k, uint256 kpriv) external view returns (uint256 r, uint256 s) {
        return FCL_ecdsa_utils.ecdsa_sign(message, k, kpriv);
    }

    function ecdsa_DerivKpub(uint256 kpriv) external view returns (uint256 x, uint256 y) {
        return FCL_ecdsa_utils.ecdsa_derivKpub(kpriv);
    }

    function ecdsa_GenKeyPair() external view returns (uint256 kpriv, uint256 x, uint256 y) {
        kpriv = block.prevrandao;
        (x, y) = FCL_ecdsa_utils.ecdsa_derivKpub(kpriv);
    }

    function ecdsa_precalc_8dim(uint256 Qx, uint256 Qy) external view returns (uint256[2][256] memory Prec) {
        return FCL_ecdsa_utils.Precalc_8dim(Qx, Qy);
    }
}

/// @notice This script deploys the FCL_Elliptic library and the wrapper contract
contract MyScript is BaseScript {
    function run() external broadcast returns (address addressOfLibrary) {
        // deploy the library contract and return the address
        addressOfLibrary = address(new FCL_ecdsa_wrapper{salt:0}());
    }
}

contract Script_Deploy_FCL_all is BaseScript {
    function run() external broadcast returns (address addressOfLibrary) {
        // deploy the library contract and return the address
        addressOfLibrary = address(new FCL_all_wrapper{salt:0}());
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

  
*/
