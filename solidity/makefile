.PHONY: tests

default: 
	deploy_FCL_all_goerlibase
	
tests:
	cd tests/WebAuthn_forge; make clean; make lint_write; make test;


#this ugly copy is related to unability of forge to deal with alias and verify contracts	

#safe way, change the sender by the Ethereum corresponding address of Nano, you must enable bling signing in the settings
#it is recommended to use a separate Nano only for deployments (do not mix with other funds)
deploy_FCL_all_linea:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=0ba7c1c10fb849dbac880aca39a0fb32 forge script script/DeployElliptic.s.sol:Script_Deploy_FCL_all --ledger --broadcast --verify --chain-id 59140 --rpc-url https://rpc.goerli.linea.build --sender 0x936632cC3B9BC47ad23D41dC7cc200015c447f71

deploy_FCL_all_sepolia:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
	ETHERSCAN_API_KEY=HURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5XHURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5X  forge script script/DeployElliptic.s.sol:Script_Deploy_FCL_all --ledger --broadcast --verify --chain-id 11155111 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public  --sender 0x936632cC3B9BC47ad23D41dC7cc200015c447f71 \

	
deploy_FCL_all_goerlibase:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
	ETHERSCAN_API_KEY=0ba7c1c10fb849dbac880aca39a0fb32  forge script script/DeployElliptic.s.sol:Script_Deploy_FCL_all --ledger --broadcast --verify --chain-id 84531 --rpc-url https://rpc.notadegen.com/base/goerli  --sender 0x936632cC3B9BC47ad23D41dC7cc200015c447f71 


#unsafe way, use it only with test accounts you could loose
unsafe_deploy_ecdsa_verify_sepolia:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<> PRIVATE_KEY=<> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 11155111 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --sender 


unsafe_deploy_ecdsa_verify_holesky:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<API_KEY> PRIVATE_KEY=<PRIVATE_KEY> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 17000 --rpc-url https://ethereum-holesky.publicnode.com --sender ;\
cd ..;rm -Rf deploy/



unsafe_deploy_ecdsa_verify_goerli:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
	ETHERSCAN_API_KEY=<> PRIVATE_KEY=<> forge script script/DeployElliptic.s.sol:MyScript --broadcast --verify --chain-id 5 --rpc-url https://ethereum-goerli.publicnode.com --sender 

unsafe_deploy_FCL_ALL_sepolia:
	mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\
	cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
ETHERSCAN_API_KEY=<> PRIVATE_KEY=<> forge script script/DeployElliptic.s.sol:Script_Deploy_FCL_all --broadcast --verify --chain-id 11155111 --rpc-url https://ethereum-sepolia.blockpi.network/v1/rpc/public --sender 


#https://explorer.goerli.linea.build/address/0x024c647BB9fb8740E9D3BC07991c88f78DeA72Cd/contracts#address-tabs

RPC: https://mainnet.optimism.io
