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

will fix any linting problem

## Deploying

>make deploy 
Will deploy the ecdsa_verify contract and verify it on chain. You need to provide a private key and an etherscan API key to the makefile.



## Current Deployments

Current create2 common to all networks address is :0xfdfbd703a269a0cb8e04304a14fbb616de68c424.

Code deployment can be check here:

[Sepolia](https://sepolia.etherscan.io/address/0xfdfbd703a269a0cb8e04304a14fbb616de68c424#code)
