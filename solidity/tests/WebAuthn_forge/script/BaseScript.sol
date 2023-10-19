// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

import {Script} from "../lib/forge-std/src/Script.sol";

abstract contract BaseScript is Script {
    /// @notice this modifier can be used as a generic broadcast solution. It will automatically either:
    ///        - use the private key provided as an environment variable to broadcast
    ///        - or starts the hardware wallet flow if the correct flags are provided and the env variable is not set
    modifier broadcast() {
        uint256 privateKey = vm.envOr("PRIVATE_KEY", uint256(0));
        privateKey != 0 ? vm.startBroadcast(privateKey) : vm.startBroadcast();

        _;

        vm.stopBroadcast();
    }
}
