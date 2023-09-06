
# if your imports fail, try this:

#export PYTHONPATH=$PYTHONPATH:$(pwd)


#compile Cairo and Sage and Compare Results
test_bls12_381:
	sage FCL_pairings/arithmetic/curves/bls12_381.sage
test_altbn128:
#	export PYTHONPATH=$$PYTHONPATH:$$(pwd)
	sage FCL_pairings/arithmetic/curves/atlbn128.py

test_ecdaa:
	sage FCL_ecdaa/tests/test_ECDAA.py

rfc_ecdaa:
	sage FCL_ecdaa/rfc/gen_rfc_ECDAA.py 
	
test_starknet:
	sage FCL_starknet/tests/FCL_starknet_test.py

test_musig2:
	#sage -c '_MU=2;nb_users=4; size_message=1;seed=3;load("FCL_musig2/tests/FCL_musig2_test.py");'   
	sage FCL_musig2/tests/FCL_musig2_test.py
all:
	make test_starknet;make test_altbn128;make test_ecdaa;make test_musig2;

