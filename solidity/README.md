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

Deploy the ecdsa_verify contract and verify it on chain. You need to provide a private key and an etherscan API key to the makefile.


## Current Deployments

| Commit # | Create2 | Mainnets | Testnets |
|--------:|---------|:--:|:----|
||         |  |         |
|[37](https://github.com/rdubois-crypto/FreshCryptoLib/commit/29f60f19d3a07ec501ce36429f9688d9be372368)| 0xebcaae4af6844b5f24a4730c5f58130977e62a2b    | [Polygon](https://polygonscan.com/address/0xebcaae4af6844b5f24a4730c5f58130977e62a2b#code)  | [Optimism](https://goerli-optimism.etherscan.io/address/0xebcaae4af6844b5f24a4730c5f58130977e62a2b#code), [Sepolia](https://sepolia.etherscan.io/address/0xebcaae4af6844b5f24a4730c5f58130977e62a2b#code), [Linea](https://explorer.goerli.linea.build/address/0xEBCaaE4Af6844B5F24A4730C5f58130977E62A2B/contracts#address-tabs)  |  


## On chain commands

The following functions are available on-chain. 

* ecdsa_verify(bytes32 message, uint256 r, uint256 s, uint256 Qx, uint256 Qy) : verify message signature (r,s) with public key (Qx, Qy)
* ecdsa_sign(bytes32 message, uint256 k, uint256 kpriv) : sign message with private key kpriv and nonce k

(Last function is provided for testing purpose only).



Here is an example of a successfull ecdsa verification using forge cast over Sepolia :

cast call 0xebcaae4af6844b5f24a4730c5f58130977e62a2b         "ecdsa_verify(bytes32,uint256,uint256,uint256,uint256)"         0xbb5a52f42f9c9261ed4361f59422a1e30036e7c32b270c8807a419feca605023         19738613187745101558623338726804762177711919211234071563652772152683725073944         34753961278895633991577816754222591531863837041401341770838584739693604822390         18614955573315897657680976650685450080931919913269223958732452353593824192568         90223116347859880166570198725387569567414254547569925327988539833150573990205 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public

