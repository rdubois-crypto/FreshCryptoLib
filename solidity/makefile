.PHONY: tests

default:
	tests
tests:
	cd tests/WebAuthn_forge; make clean; make lint_write; make test;


#this ugly copy is related to unability of forge to deal with alias and verify contracts	
deploy_ecdsa_sepolia:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<API_KEY> PRIVATE_KEY=<PRIVATE_KEY> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 11155111 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --sender 0xA099BA0dAe1f54f1Fd7238bd9d6885Af2427A28C;\
cd ..;rm -Rf deploy/

deploy_ecdsa_holesky:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<API_KEY> PRIVATE_KEY=<PRIVATE_KEY> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 17000 --rpc-url https://ethereum-holesky.publicnode.com --sender 0xA099BA0dAe1f54f1Fd7238bd9d6885Af2427A28C;\
cd ..;rm -Rf deploy/

