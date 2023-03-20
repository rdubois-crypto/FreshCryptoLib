# Webauthn.sol
Let's use that extra world computer power ...

## installing/compiling
You will need to have hardhat installed. 

To obtain gas comparisons commands is:
>export REPORT_GAS=1;npx hardhat test

## testing
There are two function implementing the ecdsa verification:
- ecdsa_verify takes as input the message hash, the r,s signature parts and the Qx, Qy public key (uncompressed) coordinates at output the boolean as result
- ecdsa_precomputed_verify takes as input only the message hash and r,s signature, expanded public key for shamir's trick being stored in precomputed.js

To modify the precomputed public key, use the "FCL_ecdsa_precompute.sage" with constant (C0,C1)=Qx, Qy to regenerate precomputed.js.
the test script push the table as a contract to later read  coefficients using extcodecopy opcode.

## Results
Initial gas cost per verification of forked repo : 1.15M 

Final gas cost verification with FCL_elliptic.sol: 

288K without precomputations, 

152K with precomputations.

The project illustrates the different attempts of optimization before reaching the current FCL_elliptic.sol.

Validate_prec is commented, as this prime attempt (no assembly, lots of access to uint[] memory variables) is to expensive (1.5M gas, 30M of deployment)

Validate_prec2 makes use of the sstore2 library

Validate_prec3 makes a direct use of the extcodecopy command to access precomputation table stored as a contract.


