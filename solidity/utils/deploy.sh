#this script is meant to be run from solidity/

#configuration

#the way to sign the transactions to send, it is recommended to use a ledger dedicated to test and deployment 
SIGNER="--ledger" #safe version, open Ethereum application on your ledger, allow blind signing
#SIGNER="PRIVATE_KEY=<>", unsafe version, using a private key of EOA
#PRIVATE_KEY
#opmainnet --chain-id 10 --rpc-url https://mainnet.optimism.io
#chain-id 137 --rpc-url https://polygon.llamarpc.com --sender 0x936632cC3B9BC47ad23D41dC7cc200015c447f71


#the public key related to your ledger/private key
SENDER=0x936632cC3B9BC47ad23D41dC7cc200015c447f71
#the script path to deploy
SCRIPT_PATH=script/DeployElliptic.s.sol
SCRIPT_FUNCTION=:Script_Deploy_FCL_all

#MAINNET
POLYGON_CHAINID=137
POLYGON_RPC=https://polygon.llamarpc.com

#TESTNETS
SEPOLIA_CHAINID=11155111
SEPOLIA_RPC=https://ethereum-sepolia.blockpi.network/v1/rpc/public

LINEA_CHAINID=59140
LINEA_TESTNET_RPC=https://rpc.goerli.linea.build

BASE_TESTNET_CHAINID=84531
BASE_TESTNET_RPC=https://goerli.base.org 

OP_TESNET_CHAINID=420
OP_TESTNET_RPC=https://optimism-goerli.public.blastapi.io


ALL_RPC=($SEPOLIA_RPC $LINEA_TESTNET_RPC $BASE_TESTNET_RPC $OP_TESTNET_RPC $POLYGON_RPC)

ALL_TESTNETWORKS=("SEPOLIA" "LINEA TESTNET" "BASE TESTNET" "OP TESTNET" "POLYGON MAINNET")
#ALL_NETWORKS=("POLYGON MAINNET")
ALL_CHAINID=($SEPOLIA_CHAINID $LINEA_CHAINID $BASE_TESTNET_CHAINID $OP_TESNET_CHAINID $POLYGON_CHAINID)


#the api key for block explorer verification 
SEPOLIA_API_KEY=HURV4UYJZCCUTXEYM73M6J6CIJE1KN1W5X
INFURA_API_KEY=0ba7c1c10fb849dbac880aca39a0fb32
LINEA_API_KEY=93H7P72W4KM16VPVTJECGCSQRXYE4MR211
BASE_API_KEY=$INFURA_API_KEY
OP_API_KEY=FV931ZWRMJCQWPHSJQ3KMPHI3CH48AFA7R
POLYGON_API_KEY=7SXPV7TTIZJXZ226UZDRYHXBV5R2AGTJYY

ALL_API_KEY=($SEPOLIA_API_KEY $LINEA_API_KEY $BASE_API_KEY $OP_API_KEY $POLYGON_API_KEY)

LAST_PR=$(git log --merges --oneline |head -n 1)
#copy FCL sources and script
mkdir deploy; cd deploy; forge init --no-git ;cp ../utils/foundry.toml_deploy foundry.toml;\cp  -R ../src/ .; cp -R ../tests/WebAuthn_forge/script/ .; \
echo "/*related PR:"$LAST_PR"*/" >> $SCRIPT_PATH

echo "******************** BEGIN DEPLOYMENT OF PR:"$LAST_PR

#deploy and verify library on all networks, polygon need --legacy
for i in ${!ALL_TESTNETWORKS[@]}; do
  echo "Chain $i is ${ALL_NETWORKS[$i]} "
  echo "     ChainID: ${ALL_CHAINID[$i]}"
  CHAIN_ID="${ALL_CHAINID[$i]}"
  
  echo "     RPC: ${ALL_RPC[$i]}" 
  RPC="${ALL_RPC[$i]}"
  API_KEY=${ALL_API_KEY[$i]}
  
  ETHERSCAN_API_KEY=$API_KEY  forge script $SCRIPT_PATH$SCRIPT_FUNCTION  --broadcast --verify --chain-id $CHAIN_ID $SIGNER --rpc-url $RPC  --sender $SENDER 
  

done

#special deployment for polygon
CHAIN_ID=$POLYGON_CHAINID
RPC=$POLYGON_RPC

ETHERSCAN_API_KEY=$API_KEY  forge script $SCRIPT_PATH$SCRIPT_FUNCTION  --broadcast --verify --legacy --chain-id $CHAIN_ID $SIGNER --rpc-url $RPC  --sender $SENDER 
  
 	
#cleaning deployment\
cd ..;rm -Rf deploy/ 


