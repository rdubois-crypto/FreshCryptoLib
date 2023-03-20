# Fresh Crypto Lib (FCL)

The Fresh Crypto Lib is a set of functions for blockchain systems such as Wallet, SmarContract.


## Content
The implemented content is:
- sec256r1 optimizations for EVM chains and Starknet
- an ecdaa implementation using blockchain primitives
- a Starknet dedicated musig2 implementation


### Directory solidity:
FCL_elliptic.sol: an EVM optimized implementation of sec256r1, using language hacks and specificities.

FCL_ecdaa.sol: an EVM version of the ECDAA anonymous attestation for anonymous airdrops

### Directory cairo:
FCL_ec_mulmuladd.cairo: an implementation of the operation aP+bQ (addition of the results of two distincts point multiplication by scalar a and b). It uses the Shamir's trick with the windowing method.
signature_opt.cairo : optimisation of ECDSA verification using ec_mulmuladd_W function

FCL_cairo_secp : optimization of the ECDSA function over sec256r1, using starkware implementation with ec_mulmuladdW_sec256k1 (original implementation from Starkware commons here:https://github.com/https://github.com/starkware-libs/cairo-lang/tree/master/src/starkware/cairo/common/cairo_secp)

FCL_cairo_secp256k1 : optimization of the ECDSA function over sec256k1 using ec_mulmuladdW_sec256r1 (original implementation from Cartridge here:https://github.com/cartridge-gg/cairo-secp256r1).


FCL_musig2: Original implementation of the Schnorr verification algorithm. Please note that it is a custom implementation (cryptographically equivalent, but not identical to BlockStream implementation).
Namely arbitrary domain separator, choice of hash, byte ordering and annoying little choices are not compatible with Musig2 BIP proposal.

### Directory sage:

### Directory examples:
Examples of integration of FCL sources.


## License 
License: This software is licensed under MIT License (see LICENSE FILE at root directory of project).


