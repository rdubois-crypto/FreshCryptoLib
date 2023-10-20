.PHONY: tests

default: unsafe_deploy_ecdsa_sepolia;
	
tests:
	cd tests/WebAuthn_forge; make clean; make lint_write; make test;


#this ugly copy is related to unability of forge to deal with alias and verify contracts	

#safe way, change the sender by the Ethereum corresponding address of Nano, you must enable bling signing in the settings
#it is recommended to use a separate Nano only for deployments (do not mix with other funds)
# for now the verification seems to fail with the --ledger option (see issue)
deploy_ecdsa_secp256r1_ledger_sepolia:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
	ETHERSCAN_API_KEY=<> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 11155111 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --ledger --sender 0x936632cC3B9BC47ad23D41dC7cc200015c447f71 ;\
	

#unsafe way
unsafe_deploy_ecdsa_sepolia:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<> PRIVATE_KEY=<> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 11155111 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --sender 


unsafe_deploy_ecdsa_holesky:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<API_KEY> PRIVATE_KEY=<PRIVATE_KEY> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 17000 --rpc-url https://ethereum-holesky.publicnode.com --sender ;\
cd ..;rm -Rf deploy/



unsafe_deploy_ecdsa_goerli:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
	ETHERSCAN_API_KEY=<> PRIVATE_KEY=<> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 5 --rpc-url https://ethereum-goerli.publicnode.com --sender 



#https://explorer.goerli.linea.build/address/0x024c647BB9fb8740E9D3BC07991c88f78DeA72Cd/contracts#address-tabs

RPC: https://mainnet.optimism.io
