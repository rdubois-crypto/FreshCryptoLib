#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project 
#///* License: This software is licensed under MIT License 	 
#///* See LICENSE file at the root folder of the project.				 
#///* FILE: FCL_ethereum.sage						         
#///* 											 
#///* 											 
#///* DESCRIPTION: emulation of Ethereum Precompiled contracts
#///* 
#//**************************************************************************************/

import hashlib, secrets


load('../FCL_common/FCL_elliptic.sage');
load('../FCL_common/FCL_bn_io.sage');
    
from sha3 import keccak_256

    
#//initialize a curve with ethereum sec256K1 parameters
_CURVE_256K1, _G_256K1 = FCL_ec_Init_Curve(sec256k_p, sec256k_a, sec256k_b, sec256k_gx, sec256k_gy, sec256k_n);
_256K1_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
_SEC256R1_BYTESIZE=32;
_MASK20B =2^160-1;


def FCL_ecRecoverPoint_from_hash(hash, v,r,s):
 Fq=GF(_256K1_ORDER);
 if r>_256K1_ORDER:
  return False;
 if s>_256K1_ORDER: 
  return False;
 
 #in ethereum 
 y=FCL_ec_decompress(_CURVE_256K1, r,v-27);
 R=_CURVE_256K1([r,y]);
 rm1=(Fq(r)^-1);
 u1=int(hash*rm1*-1);
 u2=int(s*rm1);
 
 Q=_G_256K1*u1+R*u2;
 
 return Q;


#//the actual precompile 0x05
def FCL_ecRecover(hash, v,r,s):
  Q=FCL_ecRecoverPoint_from_hash(hash, v, r, s);
  return FCL_ethereum_PubkeyToAddress(Q); #the 20 LSB bytes of keccak hash used as address
 
 

#https://en.wikipedia.org/wiki/Elliptic_Curve_Digital_Signature_Algorithm:public key recovery section
#in EVM, note that parity is encoded by constants v being equal to 28 (even y) or 27 (odd y)
def FCL_ecRecoverPoint(v,r,s,msg):
 h= int('0x'+keccak_256(msg).hexdigest(),16);
 
 return FCL_ecRecoverPoint_from_hash(h,v,r,s);


#//Convert Public Key to Address
def FCL_ethereum_PubkeyToAddress(pk):

  pre_h=FCL_BN_to_bytes(int(pk[0]), _SEC256R1_BYTESIZE)+FCL_BN_to_bytes(int(pk[1]), _SEC256R1_BYTESIZE);
  
  h= int('0x'+keccak_256(pre_h).hexdigest(),16);
  
  return h&_MASK20B;


def FCL_ecRecover_from_message(v,r,s,msg):

  Q=FCL_ecRecoverPoint(v,r,s,msg);
  
  return FCL_ethereum_PubkeyToAddress(Q); #the 20 LSB bytes of keccak hash used as address
 

 
def examples():
 print("\n***** FCL_ethereum_PubkeyToAddress() \n Example extracted from https://www.npmjs.com/package/ethereum-public-key-to-address:");
 
 #https://www.npmjs.com/package/ethereum-public-key-to-address
 Qx=0xe68acfc0253a10620dff706b0a1b1f1f5833ea3beb3bde2250d5f271f3563606;
 Qy=0x672ebc45e0b7ea2e816ecb70ca03137b1c9476eec63d4632e990020b7b6fba39;
 Qt=_CURVE_256K1([Qx,Qy]);
 print(" PubKey raw:",Qt);

 hQt=FCL_ethereum_PubkeyToAddress(Qt);
 print(" PubKey address:",hex(hQt));

 print("\n ***** FCL_ecRecover() \n Example extracted from solidity call, see tests/Precompiles/test_ecrecover.sol");

 # example from EVM using solidity code
 expected=0x2B891c3102159D688631E55a445eac38Ec9edE21;
 
 r=0xe6aa80563d9931d611917eb184059d1eaa287d627406dab4460a319c59c071a2;
 s=0x11aca1b3ecce3b4329139715e2b86dac82e596f6619f976ce186af9c4bf940d8;
 h= 0xfeeaf9a319cce6371b255b4db39156d8fbc67fb4b4cd5debc30da0abc69e420b;
 v= 28;
 Q=FCL_ecRecoverPoint_from_hash(h, v,r,s);
 print(" recovery Point:",Q);
 hash=FCL_ecRecover(h,v,r,s);
 print(" recovery address:",hex(hash));
 
 r=0xe6aa80563d9931d611917eb184059d1eaa287d627406dab4460a319c59c071a2; 	
 s=0xee535e4c1331c4bcd6ec68ea1d47925237c945f04da908cede4baef0843d0069 
 h= 0xfeeaf9a319cce6371b255b4db39156d8fbc67fb4b4cd5debc30da0abc69e420b;
 v= 28;

 hash=FCL_ecRecover(h,v,r,s);
 print(" recovery address with -s:",hex(hash));
 
 
 return 0;






