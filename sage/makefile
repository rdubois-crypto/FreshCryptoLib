# if your imports fail, try this:
export:
	export PYTHONPATH=$$PYTHONPATH:$$(pwd)
#export PYTHONPATH=$PYTHONPATH:$(pwd)

## This is a simple makefile to fasten Musig2 compilation and testing


#compile Cairo and Sage and Compare Results
test_bls12_381:
	sage FCL_pairings/arithmetic/curves/bls12_381.sage
test_altbn128:
	sage FCL_pairings/arithmetic/curves/atlbn128.py

test_ecdaa:
	sage FCL_ecdaa/tests/test_ECDAA.py

rfc_ecdaa:
	sage FCL_ecdaa/rfc/gen_rfc_ECDAA.py 

