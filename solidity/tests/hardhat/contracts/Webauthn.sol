// SPDX-License-Identifier: Apache-2.0
pragma solidity ^0.8.0;

import {Base64URL} from "./Base64URL.sol";
import {FCL_Elliptic_ZZ} from "./FCL_elliptic.sol";
import {FCL_WebAuthn} from "./FCL_Webauthn.sol";
import "hardhat/console.sol";



error InvalidAuthenticatorData();
error InvalidClientData();
error InvalidSignature();

contract Webauthn {
    uint256 public counter;

    
    function ecdsa_verif( bytes32 hash,  uint[2] calldata rs,
        uint[2] calldata Q)  public  returns (bool)
    {
        // bytes32 message = sha256(verifyData);
        console.log("hash=", uint(hash));
        console.log("rs0=", rs[0]);
        uint256 gasleft1 = gasleft();
        bool result=FCL_Elliptic_ZZ.ecdsa_verify(bytes32(hash), rs, Q);
        uint256 gasleft2 = gasleft();
        uint256 gasused = gasleft1 - gasleft2;
        if(result){
            console.log("gasused=%s", gasused);
        }
        console.log("result= %s", result);
    }
    



    function validate(
        bytes calldata authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes calldata clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] calldata rs,
        uint[2] calldata Q
    ) public {
        if (
            !FCL_WebAuthn.checkSignature(
                authenticatorData,
                authenticatorDataFlagMask,
                clientData,
                clientChallenge,
                clientChallengeDataOffset,
                rs,
                Q
            )
        ) {
            revert InvalidSignature();
        }
        counter++;
    }

}
