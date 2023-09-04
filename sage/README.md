	# Fresh Crypto Lib (FCL): sage

Some sagemath utilities for the FCL. You will need to have sagemath installed.


## Requirements

* solidity sha3 compatibility : Some functions require the old keccak_256 as implemented in Ethereum to function. You might need to install pysha3 manually from this repo 
https://github.com/tiran/pysha3/tree/coverity_scan. (git clone, then make install).
* add sage directory to python path: type >export PYTHONPATH=$PYTHONPATH:$(pwd) in sage directory


## Content



### Directory FCL_common:
Little routines common to all projects.


### Directory FCL_ecdaa:
Implementation of ecdaa over BN254 and BLS12 for blockchain implementation reference


### Directory FCL_ecdsa:
Signing, Verifying with key recovery ECDSA (canonical and non canonical versions).


### Directory FCL_ethereum:
FCL_ethereum.sage: implementations of Ethereum precompiles:
- ecrecover

FCL_ethereumhack.sage: implementation of hacky mul and  mulmuladd.

### Directory FCL_ecdsa_precompute:
Contains the precomputation of contract bytecode to be used by FCL_elliptic.sol.

Type make gentab to generate precompiled bytecode contract for the public key of coordinates (C0,C1).
Public key  may be changed by modifying the makefile constants C0 and C1.


### Directory FCL_starknet:

A Starknet implementation of Pedersen Hash, compatible with its cairo counterpart, but faster. This implementation was used by our team to win the DNA challenge of starknet Lisbon:
https://blog.ledger.com/starknet-ctf/


### Directory FCL_schnorr:
A x-only version of the schnorr algorithm, compatible with BIP340 and BIP Musig2.



## License 
License: This software is licensed under MIT License (see LICENSE FILE at root directory of project).

## external ressources:

* Pairing implementation used in protocols (ecdaa, ethereum pairing precompiles) comes from the impressive work of Aurore Guillevic's Gitlab at INRIA: https://gitlab.inria.fr/tnfs-alpha/alpha/-/tree/190b87732901750ed1438a8cf340571531d32230/sage/tnfs
