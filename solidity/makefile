.PHONY: tests

default: tests

#test using forge	
tests:
	cd tests/WebAuthn_forge; make clean; make lint_write; make test;

#deploy the FCL on chain
deploy:
	./utils/deploy.sh
#once deployed, test using cast
cast_test:
	 ./solidity/tests/cast/cast_call.sh 	
cast_send:


