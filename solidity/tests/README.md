# Fresh Crypto Lib (FCL) : tests

A set of examples of integration of the FCL.

## Requirements

In order to run the tests in this directory [sage](https://www.sagemath.org/) is required.

On macos you can install it using brew by following command:

```bash
brew install --cask sage
```

On linux you can install it using apt by following command:

```bash
apt install sagemath
```

## Content

Here is a comprehensive listing of the subdirectories:

### Precompiles:

This directory contain the logic to generate the precompiles.

`test_ecrecover.sol`: generate vectors to assess sage implementation of ecrecover


### cast:

Contains automation script to test deployed contracts using cast commands:
* cast_call.sh : contains only 'cast call' commands (testing transactions are off chain)
* cast_send.sh : contains 'cast send' commands, sending true transactions on chain


### hardhat:

**!beware that hardhat is not maintained, forge is now the preferred environment.**

Implementation of the ECDSA P256 for WebAuthn using [hardhat](https://hardhat.org/).
WebAuthn part forked from https://github.com/btchip/Webauthn.sol

Optimization of ecdsa verification using FCL_elliptic.sol library.

Initial gas cost per verification: 1.15M

Final gas cost verification:

-208K without precomputations,

-139K with precomputations.

#### Specific requirements

- In order to run the tests in this directory [hardhat](https://hardhat.org/) is required. Follow this link for the proper installation: https://hardhat.org/getting-started/#installation

### WebAuthn_forge:

Implementation of the ECDSA P256 for WebAuthn using [forge by foundry](https://book.getfoundry.sh/).
WebAuthn part forked from https://github.com/btchip/Webauthn.sol
The framework includes the following libraries to bench and compare execution time:
- orbs-network
- Numerology
- obvioustech

Final gas cost verification:

-202K without precomputations,
-69K with precomputations.


#### Specific requirements

- In order to run the tests in this directory [forge by foundry](https://book.getfoundry.sh/) is required. Follow this link for the proper installation: https://book.getfoundry.sh/getting-started/installation

- Install the dependencies by running `forge install` in the `WebAuthn_forge` directory.

- Execute the remapping command by running `forge remappings`. This command get the automatically inferred remappings for the project.
