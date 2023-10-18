# Fresh Crypto Lib (FCL) : solidity


## Content

### FCL_elliptic.sol: 
Optimized Implementation of the ECDSA P256 using XYZZ coordinates.

Gas cost of ecdsa verification (with 100000 runs in configuration files):
- 200K without precomputation table
- 75K with precomputations table of 16kb

Full test with WychProof vectors and comparizon to existing libraries (orbs-network, obvioustech) is available in tests/WebAuthn_forge.

### FCL_sha512.sol: 

A sha512 implementation of a one step SHA512 (processing Init/Update/Final at once) for message of size lesser than 61 bytes.
Intended to be used by EdDSA only.

### FCL_eddsa.sol: 

Optimized implementation of EDDSA over ed25519 as specified by RFC 8032 https://datatracker.ietf.org/doc/html/rfc8032.
API are not mature, only one block message can be handled

Gas cost: 270K, (compared to best of our knowledge here at 1.2M: https://github.com/javgh/ed25519-solidity



### FCL_Webauthn:
implementation of WebAuthn authentication mechanism on top of P256/sec256r1 ecdsa

 
gas cost of WebAuthn verification (to be check again in forge):

- 257K without precomputations,
- 137K with precomputations (16kb, using FCL_ecdsa_precompute.sage).


### Testing:
See directory tests/Webauthn_forge. (hardhat environment is deprecated.) 
