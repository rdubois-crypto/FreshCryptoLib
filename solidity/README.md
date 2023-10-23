# Fresh Crypto Lib (FCL) : SOLIDITY

The Fresh Crypto Lib is a set of functions for blockchain systems such as Wallet, SmartContracts.
This directory contains implementation in solidity.


## Content

* FCL_ecdsa.sol: an EVM optimized implementation of ecdsa over sec256r1(P256), using language hacks and specificities, as described in paper https://eprint.iacr.org/2023/939.pdf.
* FCL_eddsa.sol   : an EVM optimized implementation of ed25519, using same paper tricks.
* FCL_sha512.sol : implementation of the SHA512 primitive (single bloc implementation)
* FLC_Webauthn.sol: implementation of the WebAuthn2/FIDO2 authentication over ECDSA with P256
<!--- FCL_ecdaa.sol: an EVM version of the ECDAA anonymous attestation for anonymous airdrops -->

## Compiling/Testing

To compile the sources, go to tests/WebAuthn_forge, then run

>make test

CI's are implemented, PR will be successfull if linting and all tests are OK. There are usefull commands in the makefile, such as 

>make lint-write 

Fix linting problem

## Deploying

>make deploy

Deploy the ecdsa_verify contract and verify it on chain. You need to provide a private key and an etherscan API key to the makefile. More RPC's, ChainID's and faucets are available in utils/rpc.md.


## Current Deployments

| Commit # | Create2 | Mainnets | Testnets | Functions |
|--------:|---------|:--:|:----|:----|
||         |  |         | |
|[37](https://github.com/rdubois-crypto/FreshCryptoLib/commit/29f60f19d3a07ec501ce36429f9688d9be372368)| 0xebcaae4af6844b5f24a4730c5f58130977e62a2b    | [Polygon](https://polygonscan.com/address/0xebcaae4af6844b5f24a4730c5f58130977e62a2b#code)  | [Optimism](https://goerli-optimism.etherscan.io/address/0xebcaae4af6844b5f24a4730c5f58130977e62a2b#code), [Sepolia](https://sepolia.etherscan.io/address/0xebcaae4af6844b5f24a4730c5f58130977e62a2b#code), [Linea](https://explorer.goerli.linea.build/address/0xEBCaaE4Af6844B5F24A4730C5f58130977E62A2B/contracts#address-tabs)  |  ecdsa_sign, ecdsa_verify |
|[43](https://github.com/rdubois-crypto/FreshCryptoLib/pull/43)| 0xbba700416337e99deb4b823a752e35db1d329c0d    | N/A |  [Sepolia](https://sepolia.etherscan.io/address/0xbba700416337e99deb4b823a752e35db1d329c0d#code) |  ecdsa_sign, ecdsa_verify |

## On chain commands

### Available commands
The following functions are available on-chain. 

Minimal verification contract:
* ecdsa_verify(bytes32 message, uint256 r, uint256 s, uint256 Qx, uint256 Qy) : verify message signature (r,s) with public key (Qx, Qy)

All utils (Signature, Keygen, derivation, Precomputations) 
* ecdsa_sign(bytes32 message, uint256 k, uint256 kpriv) : sign message with private key kpriv and nonce k

(Last function is provided for testing purpose only).


### Full example

Here is an example of a successfull ecdsa verification using forge cast over Sepolia :

#### derivate the public key from private key sk=3 (just for testing, insecure)
** cast call 0xbba700416337e99deb4b823a752e35db1d329c0d "ecdsa_DerivKpub(uint256)" 3 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public **
 
return value 0x5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032

Qx=0x5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c

Qcy0x8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032

#### Sign a message (just for testing, insecure)
** cast call 0xbba700416337e99deb4b823a752e35db1d329c0d "ecdsa_sign(bytes32,uint256,uint256)" 0xbb5a52f42f9c9261ed4361f59422a1e30036e7c32b270c8807a419feca605023 12 3 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public **

return 0x741dd5bda817d95e4626537320e5d55179983028b2f82c99d500c5ee8624e3c4974efc58adfdad357aa487b13f3c58272d20327820a078e930c5f2ccc63a8f2b

r=0x741dd5bda817d95e4626537320e5d55179983028b2f82c99d500c5ee8624e3c4

s=0x974efc58adfdad357aa487b13f3c58272d20327820a078e930c5f2ccc63a8f2b

#### Verify a message
**cast call 0xbba700416337e99deb4b823a752e35db1d329c0d "ecdsa_verify(bytes32,uint256,uint256,uint256,uint256)" 0xbb5a52f42f9c9261ed4361f59422a1e30036e7c32b270c8807a419feca605023 0x741dd5bda817d95e4626537320e5d55179983028b2f82c99d500c5ee8624e3c4 0x974efc58adfdad357aa487b13f3c58272d20327820a078e930c5f2ccc63a8f2b 0x5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c  0x8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public **

return true

