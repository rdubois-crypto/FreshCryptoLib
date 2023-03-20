# Fresh Crypto Lib (FCL) : solidity


## Content

### FCL_elliptic.sol: 
Implementation of the ECDSA P256 using XYZZ coordinates.


Gas cost of ecdsa verification:

- 272.6K without precomputations,

- 152K with precomputations (16kb, using FCL_ecdsa_precompute.sage).

