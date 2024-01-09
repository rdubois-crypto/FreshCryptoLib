# Fresh Crypto Lib (FCL)

The Fresh Crypto Lib is a set of functions for blockchain systems such as Wallet, SmartContracts.


## Content

The implemented content is:
- sec256r1 optimizations for EVM chains and Starknet
- an ecdaa implementation using blockchain primitives
- a Starknet dedicated musig2 implementation


### Directory solidity:

* FCL_ecdsa.sol: an EVM optimized implementation of ecdsa over sec256r1(P256), using language hacks and specificities, as described in paper https://eprint.iacr.org/2023/939.pdf.
* FCL_eddsa.sol   : an EVM optimized implementation of ed25519, using same paper tricks.
* FCL_sha512.sol : implementation of the SHA512 primitive (single bloc implementation)
* FLC_Webauthn.sol: implementation of the WebAuthn2/FIDO2 authentication over ECDSA with P256
<!--- FCL_ecdaa.sol: an EVM version of the ECDAA anonymous attestation for anonymous airdrops -->


| PR # | Create2 | Mainnets | Testnets |
|--------:|---------|:--:|:----|
||         |  |         |
|[46](https://github.com/rdubois-crypto/FreshCryptoLib/pull/46)| 0xE9399D1183a5cf9E14B120875A616b6E2bcB840a    | [Polygon](https://polygonscan.com/address/0xe9399d1183a5cf9e14b120875a616b6e2bcb840a#code)  | [Optimism](https://goerli-optimism.etherscan.io/address/0xe9399d1183a5cf9e14b120875a616b6e2bcb840a#code), [Sepolia](https://sepolia.etherscan.io/address/0xe9399d1183a5cf9e14b120875a616b6e2bcb840a#code), [Linea](https://explorer.goerli.linea.build/address/0xE9399D1183a5cf9E14B120875A616b6E2bcB840a/contracts#address-tabs)  |  

(code is verified except for linea, one could check that bytecode is identical).

### Directory cairo0.9:

* FCL_ec_mulmuladd.cairo: an implementation of the operation aP+bQ (addition of the results of two distincts point multiplication by scalar a and b). It uses the Shamir's trick with the windowing method.
signature_opt.cairo : optimisation of ECDSA verification using ec_mulmuladd_W function

<!---* FCL_cairo_secp : optimization of the ECDSA function over sec256r1, using starkware implementation with ec_mulmuladdW_sec256k1 (original implementation from Starkware commons here:https://github.com/https://github.com/starkware-libs/cairo-lang/tree/master/src/starkware/cairo/common/cairo_secp)-->

<!---* FCL_cairo_secp256k1 : optimization of the ECDSA function over sec256k1 using ec_mulmuladdW_sec256r1 (original implementation from Cartridge here:https://github.com/cartridge-gg/cairo-secp256r1) -->


<!---*FCL_musig2: Original implementation of the Schnorr verification algorithm. Please note that it is a custom implementation (cryptographically equivalent, but not identical to BlockStream implementation).
Namely arbitrary domain separator, choice of hash, byte ordering and annoying little choices are not compatible with Musig2 BIP proposal.-->

Note : The language is now deprecated since its transition from python-like to rust-like language.

### Directory sage:

* FCL_ecdsa_precompute.sage : precompute bytecode contract to speed up ecdsa verification for a given key.
* FCL_ecdaa : sage reference for a blockchain implementation of ECDAA
* FCL_pairings : sage implementation of curve and pairing computation over BN254 (aka altbn128) and BLS12381 using INRIA sources.


### Acknowledments:


#### Building Blocks

The following repos are used as building blocks in the FCL:
* Aurore Guillevic's Gitlab at INRIA: https://gitlab.inria.fr/tnfs-alpha/alpha/-/tree/190b87732901750ed1438a8cf340571531d32230/sage/tnfs for its generic sagemath BN and BLS curves and pairing implementation.
* Paul Miller **Noble** javascript library for its G1 implementation of BN254 and BLS12, and keccak256. https://paulmillr.com/noble/

#### Benchmark

The following repos have been used in benchmarks:
* Alembic/cometh:https://github.com/alembic-tech/P256-verify-signature/blob/main/contracts/EllipticCurve.sol
* MaxRobot : https://github.com/maxrobot/elliptic-solidity
* Numerology : https://github.com/nucypher/numerology
* Obvious : https://github.com/itsobvioustech/aa-passkeys-wallet


#### Testing 

* Invariant test vectors are extracted from the wycheproof project : https://github.com/google/wycheproof

### FCL in the wild

* Academic paper: https://eprint.iacr.org/2023/939 for EthCC2023
* Alembic : https://github.com/alembic-tech/p256-signer/blob/main/contracts/FCL/FCL_elliptic.sol
* Braavos https://github.com/myBraavos/efficient-secp256r1/blob/develop/src/secp256r1/ec_mulmuladd.cairo
* Cartridge https://github.com/cartridge-gg/cairo-secp256r1/pull/3
* EIP665 PR#7515 :https://github.com/ethereum/EIPs/pull/7515
* Presentation made at EthCC 2023 in Paris: https://www.youtube.com/live/Rlq21oA_FA8
* Forum DAO,  :https://github.com/forumdaos/forum-contracts/tree/main/src/libraries
* Daimo, ethereum payments : https://github.com/daimo-eth/p256-verifier/blob/master/src/P256Verifier.sol
* Safe : https://github.com/safe-global/safe-modules/tree/master/4337/contracts/test/FCL
* Wallet Abstraction, EthGlobal NY hackathon finalist : https://github.com/qd-qd/wallet-abstraction
* Wax : Wallet Account Experiments, PSE team (EF funded) https://github.com/getwax/wax/commit/fd246685f67ebb56f5ebb82cb48df56088e28c7b  

Do not hesitate to reach us to be listed.
## License 
License: This software is licensed under MIT License (see LICENSE FILE at root directory of project).


