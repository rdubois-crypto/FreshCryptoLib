#the way to sign the transactions to send 
SIGNER="--ledger" #safe version, open Ethereum application on your ledger, allow blind signing
#$SIGNER="PRIVATE_KEY=<>", unsafe version, using a private key of EOA

#list of RPCs
LINEA_TESTNET=https://rpc.goerli.linea.build
BASESCAN_TESTNET_RPC=https://goerli.base.org 

FCL_ADDRESS=0xE9399D1183a5cf9E14B120875A616b6E2bcB840a
DAIMO_ADDRESS=0xc2b78104907F722DABAc4C69f826a522B2754De4 #Daimo is only testnet is basescan as of 25/10/23

# pick a network and RPC
NETWORK="BASE"
RPC=$BASESCAN_TESTNET_RPC

#pick a seed
SEED=$( echo $RANDOM)
#pick random values for input (testing purpose)
KPRIV=$(cast keccak "1"$SEED)
NONCE=$(cast keccak "2"$SEED)
MESSAGE=$(cast keccak "3"$SEED)

#derivate the public key value from private
KPUB=$(cast call $FCL_ADDRESS "ecdsa_DerivKpub(uint256)" $KPRIV --rpc-url $RPC)
Qx=${KPUB:2:64}
Qy=${KPUB:66:64}

echo "Kpub:"$Qx" "$Qy

########################################################################################################
#send precomputations over network

echo "network:"$NETWORK
echo "seed:"$SEED


PREC=$(cast call $FCL_ADDRESS "ecdsa_precalc_8dim(uint256, uint256)" $Qx  $Qy --rpc-url $RPC)
DEPLOY_CONSTANT="0x600B5981380380925939F3" #prefix constant to deploy bytecode
BYTECODE=$DEPLOY_CONSTANT${PREC:2:${#PREC}}

echo $(cast send $SIGNER --rpc-url $RPC  --create $BYTECODE)  >> res.dat
ADDRESS_PREC=$(awk -v col=6 '{print $col}' res.dat)
echo "Deployed precomputation address:"$ADDRESS_PREC


for i in 1 5
do



########################################################################################################

echo "message:"$MESSAGE

echo "*********************************** VERIFICATIONS, classic API"
echo ""

########################################################################################################
echo "######################"
#compute signature (r,s) of message with nonce and private key
SIG=$(cast call $FCL_ADDRESS "ecdsa_sign(bytes32, uint256, uint256)" $MESSAGE $NONCE $KPRIV --rpc-url $RPC)
r=${SIG:2:64}
s=${SIG:66:64}

echo "signature:"$SIG
echo "r:"$r
echo "n s:"$s

#verify signature of message with public key


echo "######################"
VERIF=$(cast call $FCL_ADDRESS "ecdsa_verify(bytes32, uint256, uint256, uint256, uint256)" $MESSAGE $r $s $Qx $Qy --rpc-url $RPC)

echo "verif:"$VERIF
echo "######################"
cast send $SIGNER $FCL_ADDRESS "ecdsa_verify(bytes32, uint256, uint256, uint256, uint256)" $MESSAGE $r $s $Qx $Qy --rpc-url $RPC


########################################################################################################


#verify using precomputations at address $ADDRESS_PREC
VERIF2=$(cast call $FCL_ADDRESS "ecdsa_precomputed_verify(bytes32,uint256,uint256,address)" $MESSAGE $r $s $ADDRESS_PREC --rpc-url $RPC)
cast send $SIGNER $FCL_ADDRESS "ecdsa_precomputed_verify(bytes32,uint256,uint256,address)" $MESSAGE $r $s $ADDRESS_PREC --rpc-url $RPC
echo "verif with prec:"$VERIF
echo "######################"



EIP7212_INPUT=$MESSAGE$r$s$Qx$Qy
echo "EIP7212 input:"$EIP7212_INPUT
VERIF2=$(cast call $FCL_ADDRESS $EIP7212_INPUT --rpc-url $RPC)
echo "*** verif with FCL (precompile):"$VERIF2" with RPC="$RPC

echo "*********************************** VERIFICATIONS, precompile API"
echo ""

EIP7212_INPUT=$MESSAGE$r$s$Qx$Qy
echo "EIP7212 input:"$EIP7212_INPUT
VERIF2=$(cast call $FCL_ADDRESS $EIP7212_INPUT --rpc-url $RPC)
echo "*** FCL verif result:"$VERIF2" with RPC="$RPC
cast send $SIGNER $FCL_ADDRESS $EIP7212_INPUT --rpc-url $RPC


VERIF2=$(cast call $DAIMO_ADDRESS $EIP7212_INPUT --rpc-url $RPC)
echo "*** Daimo verif result:"$VERIF2" with RPC="$RPC
cast send $SIGNER $DAIMO_ADDRESS $EIP7212_INPUT --rpc-url $RPC


#randomize input (except key) for next iteration
SEED=$( echo $RANDOM)
NONCE=$(cast keccak "2"$SEED)
MESSAGE=$(cast keccak "3"$SEED)

done
rm res.dat



