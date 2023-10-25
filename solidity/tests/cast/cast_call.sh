#list of RPCs
LINEA_TESTNET=https://rpc.goerli.linea.build

FCL_ADDRESS=0xe55ccF9bf490a42AE6E31ab1A429c475c571a05d
 
# pick a network and RPC
NETWORK="LINEA_TEST"
RPC=$LINEA_TESTNET

#pick a seed
SEED=$( echo $RANDOM)
#pick random values for input (testing purpose)
KPRIV=$(cast keccak "1"$SEED)
NONCE=$(cast keccak "2"$SEED)
MESSAGE=$(cast keccak "3"$SEED)

echo "network:"$NETWORK
echo "seed:"$SEED

#derivate the public key value from private
KPUB=$(cast call $FCL_ADDRESS "ecdsa_DerivKpub(uint256)" $KPRIV --rpc-url $RPC)
Qx=${KPUB:2:64}
Qy=${KPUB:66:64}

echo "Kpub:"$Qx" "$Qy

#compute signature (r,s) of message with nonce and private key
SIG=$(cast call $FCL_ADDRESS "ecdsa_sign(bytes32, uint256, uint256)" $MESSAGE $NONCE $KPRIV --rpc-url $RPC)
r=${SIG:2:64}
s=${SIG:66:64}

echo "signature:"$SIG

#verify signature of message with public key
VERIF=$(cast call $FCL_ADDRESS "ecdsa_verify(bytes32, uint256, uint256, uint256, uint256)" $MESSAGE $r $s $Qx $Qy --rpc-url $RPC)

echo "verif:"$VERIF

# chaîne de caractères
if [ $VERIF = 0x0000000000000000000000000000000000000000000000000000000000000001 ] 
then 
   echo "OK"
else
   echo "KO"   
fi


