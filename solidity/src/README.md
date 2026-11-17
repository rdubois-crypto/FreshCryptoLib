# Fresh Crypto Lib (FCL) : solidity


## Content

### FCL_elliptic.sol: 
Optimized Implementation of the ECDSA P256 using XYZZ coordinates.

uG+vQ processes 2 bits of u and 2 bits of v per step with the 15 points iG+jQ. The 12 points depending on Q are
normalized with a single modular inversion, the doubling costs 11 mulmod, the mixed addition needs no `sub`.
The variants using a 16kb precomputed table transpose the scalars into nibbles once.

Gas cost (forge, solc 0.8.30, 1000000 optimizer runs, external call included). Osaka prices the modexp
precompile three times higher than Prague (EIP-7883), hence two columns:

| | Osaka | Prague |
|---|---:|---:|
| ecdsa_verify, no precomputation | 157.7K (was 212.5K) | 149.5K (was 204.3K) |
| ecdsa_precomputed_hackmem, 16kb table | 62.0K (was 74.5K) | 56.5K (was 69.0K) |

The runtime code of a contract exposing only ecdsa_verify grows from 5.7kb to 10.4kb.

Full test with WychProof vectors and comparizon to existing libraries (orbs-network, obvioustech) is available in tests/WebAuthn_forge.
The optimized engine is also checked against a python oracle, including the degenerate keys Q = -G and Q = c.G for
small c, and against the previous version kept in tests/WebAuthn_forge/test/ref (FCL_elliptic_opt.t.sol).

### FCL_sha512.sol: 

A sha512 implementation of a one step SHA512 (processing Init/Update/Final at once) for message of size lesser than 61 bytes.
Intended to be used by EdDSA only.

### FCL_eddsa.sol: 

Optimized implementation of EDDSA over ed25519 as specified by RFC 8032 https://datatracker.ietf.org/doc/html/rfc8032.
API are not mature, only one block message can be handled

Gas cost: 270K, (compared to best of our knowledge here at 1.2M: https://github.com/javgh/ed25519-solidity



### FCL_Webauthn:
implementation of WebAuthn authentication mechanism on top of P256/sec256r1 ecdsa

 
gas cost of WebAuthn verification without precomputations (test_webauthn_Base64URL_checkSignature):

- 169.6K under Osaka, 161.4K under Prague (was 227.2K under Osaka).

### FCL_babyjubjub.sol:

EdDSA verification over BabyJubJub as produced by circomlib. The twisted Edwards curve
168700x^2+y^2=1+168696x^2y^2 is mapped through its Montgomery form to a short Weierstrass curve, then rescaled
to a=-3 since -3/a is a fourth power mod r: the FCL_elliptic.sol engine applies unchanged, with lazy reduction
as the field has 254 bits.

- verifyRaw(hm, R8x, R8y, S, Ax, Ay): checks S.B8 = R8 + 8.hm.A, hm being the hashed challenge
  (Poseidon, MiMC7, MiMCSponge or Pedersen in circomlib), with the input checks of circomlibjs.
- verifyRawW(hm, R8x, R8y, S, Wx, Wy) with (Wx, Wy) = pubKeyToW8(Ax, Ay): same check with a key precomputed once.

| | Osaka | Prague |
|---|---:|---:|
| verifyRaw | 138.0K | 132.6K |
| verifyRawW | 132.0K | 129.3K |
| pubKeyToW8 | 6.7K | 4.0K |

### FCL_babyjubjub_edwards.sol:

verifyRaw on the a=-1 model x -> sqrt(-168700).x, where the extended formulas are complete: 120.9K under Osaka,
118.2K under Prague, about 12 % below the Weierstrass version.


### Testing:
See directory tests/Webauthn_forge. (hardhat environment is deprecated.) 
