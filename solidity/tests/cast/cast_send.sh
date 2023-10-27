#the way to sign the transactions to send 
SIGNER="--ledger" #safe version, open Ethereum application on your ledger, allow blind signing
#$SIGNER="PRIVATE_KEY=<>", unsafe version, using a private key of EOA

#list of RPCs
LINEA_TESTNET=https://rpc.goerli.linea.build
FCL_ADDRESS=0xEd0D252a3A26FB0269333BD6Cc720a8a68a68fcb
 
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

echo "message:"$MESSAGE

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
echo "r:"$r
echo "n s:"$s

#verify signature of message with public key
VERIF=$(cast call $FCL_ADDRESS "ecdsa_verify(bytes32, uint256, uint256, uint256, uint256)" $MESSAGE $r $s $Qx $Qy --rpc-url $RPC)

echo "verif:"$VERIF

PREC=$(cast call $FCL_ADDRESS "ecdsa_precalc_8dim(uint256, uint256)" $Qx  $Qy --rpc-url $RPC)
DEPLOY_CONSTANT="0x600B5981380380925939F3" #prefix constant to deploy bytecode
BYTECODE=$DEPLOY_CONSTANT${PREC:2:${#PREC}}

#send precomputations over network
echo $(cast send $SIGNER --rpc-url https://rpc.goerli.linea.build  --create $BYTECODE)  >> res.dat
ADDRESS_PREC=$(awk -v col=6 '{print $col}' res.dat)
echo "Deployed precomputation address:"$ADDRESS_PREC


#verify using precomputations at address $ADDRESS_PREC
VERIF2=$(cast call $FCL_ADDRESS "ecdsa_precomputed_verify(bytes32,uint256,uint256,address)" $MESSAGE $r $s $ADDRESS_PREC --rpc-url $RPC)

echo "verif with prec:"$VERIF

# chaîne de caractères
if [ $VERIF = 0x0000000000000000000000000000000000000000000000000000000000000001 ] 
then 
   echo "OK"
else
   echo "KO"   
fi


