// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.20;

import "forge-std/Test.sol";
//FreshCryptoLib implementation
import "@solidity/FCL_elliptic.sol";
import "@solidity/FCL_Webauthn.sol";

//orbs_network implementation
import "@external/ECops.sol";

contract WebAuthn_bench {
    uint256 constant _FCL_ID = 0;
    uint256 constant _ORBS_ID = 1;
    uint256 constant _ALEMBICH_ID = 2;

    uint256 constant _MAXID = 0;

    uint256 corelib_ID;
    uint256 dataPointer;

    function set_dataPointer(uint256 i_dataPointer) public {
        dataPointer = i_dataPointer;
    }

    constructor(uint256 libID) {
        if (libID > _MAXID) {
            revert();
        }

        corelib_ID = libID;
    }
}
