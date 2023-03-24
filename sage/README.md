# Fresh Crypto Lib (FCL): sage

Some sagemath utilities for the FCL. You will need to have sagemath installed.

## Content



### Directory FCL_common:
Little routines common to all projects.

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

### Directory FCL_schnorr:
A x-only version of the schnorr algorithm, compatible with BIP340 and BIP Musig2.



## License 
License: This software is licensed under MIT License (see LICENSE FILE at root directory of project).


