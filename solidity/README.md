# Fresh Crypto Lib (FCL) : solidity


## Content

### FCL_elliptic.sol: 
Optimized Implementation of the ECDSA P256 using XYZZ coordinates.

Gas cost of ecdsa verification (with 100000 runs in configuration files):
- 200K without precomputation table
- 75K with precomputations table of 16kb
### FCL_Webauthn:
implementation of WebAuthn authentication mechanism on top of P256/sec256r1 ecdsa

 
gas cost of WebAuthn verification (to be check again in forge):

- 257K without precomputations,
- 137K with precomputations (16kb, using FCL_ecdsa_precompute.sage).
