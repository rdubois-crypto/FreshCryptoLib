# Fresh Crypto Lib (FCL) : Examples


## Content
### WebAuthn: 
Implementation of the ECDSA P256
Forked from https://github.com/btchip/Webauthn.sol

Optimization of ecdsa verification using FCL_elliptic.sol library.

Initial gas cost per verification: 1.15M

Final gas cost verification: 

-288K without precomputations, 

-152K with precomputations.

