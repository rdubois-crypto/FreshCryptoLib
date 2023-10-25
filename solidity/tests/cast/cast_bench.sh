#the way to sign the transactions to send 
SIGNER="--ledger" #safe version, open Ethereum application on your ledger, allow blind signing

#list of RPCs
LINEA_TESTNET_RPC=https://rpc.goerli.linea.build
SEPOLIA_TESTNET_RPC=https://ethereum-sepolia.blockpi.network/v1/rpc/public
BASESCAN_TESTNET=https://goerli.base.org 


FCL_ADDRESS=0xe55ccF9bf490a42AE6E31ab1A429c475c571a05d   #This FCL available on Sepolia, basescan and linea testnets
DAIMO_ADDRESS=0xc2b78104907F722DABAc4C69f826a522B2754De4 #Daimo is only testnet is basescan as of 25/10/23
VYPER_ADDRESS=0xD99D0f622506C2521cceb80B78CAeBE1798C7Ed5 #Vyper is available over sepolia
 
 
 
# pick a network and RPC
for iteration in 1 2
do
echo
RPC=$BASESCAN_TESTNET
#pick a seed
SEED=$( echo $RANDOM)
#pick random values for input (testing purpose)
KPRIV=$(cast keccak "1"$SEED)
NONCE=$(cast keccak "2"$SEED)
MESSAGE=$(cast keccak "3"$SEED)

echo "seed:"$SEED

#derivate the public key value from private
KPUB=$(cast call $SIGNER $FCL_ADDRESS "ecdsa_DerivKpub(uint256)" $KPRIV --rpc-url $RPC)
Qx=${KPUB:2:64}
Qy=${KPUB:66:64}

echo "Kpub:"$Qx" "$Qy

#compute signature (r,s) of message with nonce and private key
SIG=$(cast call $FCL_ADDRESS "ecdsa_sign(bytes32, uint256, uint256)" $MESSAGE $NONCE $KPRIV --rpc-url $RPC)
r=${SIG:2:64}
s=${SIG:66:64}

echo "signature:"$SIG

#verify signature of message with public key
VERIF=$(cast send $SIGNER $FCL_ADDRESS "ecdsa_verify(bytes32, uint256, uint256, uint256, uint256)" $MESSAGE $r $s $Qx $Qy --rpc-url $RPC)

echo "*** verif with FCL:"$VERIF" with RPC="$RPC

#testing Vyper
RPC=$SEPOLIA_TESTNET_RPC
EIP7212_INPUT=$MESSAGE$r$s$Qx$Qy

echo "EIP7212 input:"$EIP7212_INPUT
VERIF2=$(cast send $SIGNER $VYPER_ADDRESS $EIP7212_INPUT --rpc-url $RPC)

echo "*** verif vyper:"$VERIF2" with RPC="$RPC

#testing Daimo
RPC=$BASESCAN_TESTNET
VERIF3=$(cast send $SIGNER $DAIMO_ADDRESS $EIP7212_INPUT --rpc-url $RPC)
echo "verif with daimo:"$VERIF3" with RPC="$RPC


# chaîne de caractères
done

