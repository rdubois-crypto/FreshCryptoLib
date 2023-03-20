# Webauthn.sol
Let's use that extra world computer power ...

## installing/compiling
You will need to have hardhat installed. 

To obtain gas comparisons commands is:
>export REPORT_GAS=1;npx hardhat test


## Results
Initial gas cost per verification of forked repo : 1.15M 
Final gas cost verification with FCL_elliptic.sol: 288K without precomputations, 152K with precomputations.

The project illustrates the different attempts of optimization before reaching the current FCL_elliptic.sol.

Validate_prec is commented, as this prime attempt (no assembly, lots of access to uint[] memory variables) is to expensive (1.5M gas, 30M of deployment)
Validate_prec2 makes use of the sstore2 library
Validate_prec3 makes a direct use of the extcodecopy command to access precomputation table stored as a contract.


