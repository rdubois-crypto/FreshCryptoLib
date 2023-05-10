# Fresh Crypto Lib (FCL) : tests

A set of examples of integration of the FCL.

## Content
### Precompiles: 
test_ecrecover.sol: generate vectors to assess sage implementation of ecrecover

### WebAuthn_hardhat: 

Implementation of the ECDSA P256 for WebAuthn using hardhat.
WebAuthn part forked from https://github.com/btchip/Webauthn.sol

Optimization of ecdsa verification using FCL_elliptic.sol library.

Initial gas cost per verification: 1.15M 

Final gas cost verification: 

-257K without precomputations, 

-139K with precomputations.


### WebAuthn_hardhat: 
Implementation of the ECDSA P256 for WebAuthn using forge.
WebAuthn part forked from https://github.com/btchip/Webauthn.sol

requirements:
- install forge
execution:
>forge test -vv --ffi
