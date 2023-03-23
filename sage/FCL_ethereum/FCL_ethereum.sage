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

#https://en.wikipedia.org/wiki/Elliptic_Curve_Digital_Signature_Algorithm:public key recovery section
#in EVM, note that parity is encoded by constants v being equal to 28 (even y) or 27 (odd y)
def FCL_ecRecoverPoint(v,r,s,msg):
 Fq=GF(_256K1_ORDER);
 if r>_256K1_ORDER:
  return False;
 if s>_256K1_ORDER: 
  return False;
 
 #in ethereum 
 y=FCL_ec_decompress(_CURVE_256K1, r, 28-v);
 R=_CURVE_256K1([r,y]);
 e= int('0x'+keccak_256(msg).hexdigest(),16);
 rm1=(Fq(r)^-1)
 u1=int(e*rm1);
 u2=int(s*rm1);
 
 Q=_G_256K1*u1+R*u2;
 
 return Q;

#//Convert Public Key to Address
def FCL_ethereum_PubkeyToAddress(pk):

  pre_h=FCL_BN_to_bytes(int(pk[0]), _SEC256R1_BYTESIZE)+FCL_BN_to_bytes(int(pk[1]), _SEC256R1_BYTESIZE);
  
  h= int('0x'+keccak_256(pre_h).hexdigest(),16);
  
  return h&_MASK20B;


#//the actual precompile 0x05
def FCL_ecRecover(v,r,s,msg):
  Q=FCL_ecRecoverPoint(v,r,s,msg);
  
  return FCL_ethereum_PubkeyToAddress(Q); #the 20 LSB bytes of keccak hash used as address
 
def examples(): 
#example extracted from https://cryptobook.nakov.com/digital-signatures/ecdsa-sign-verify-examples
 msg="Message for signing".encode('iso-8859-1');
 print("m:\n",msg);
 sk=0x68abc765746a33272e47b0a96a0b4184048f70354221e04219fbc223bfe79794
 pk=sk*_G_256K1;
 print("pub:",hex(pk[0]), hex(pk[1])); 	
 r = 0x4cddf146c578d20a31fa6128e5d9afe6ac666e5ef5899f2767cacb56a42703cc;
 s = 0x3847036857aa3f077a2e142eee707e5af2653baa99b9d10764a0be3d61595dbb;
 v = 0x0;

 Q=FCL_ecRecoverPoint(28, r,s,msg);
 print("pre recovery:",hex(pk[0]), hex(pk[1])); 	
 h=FCL_ecRecover(28, r,s,msg);
 print("recovery hash:",hex(h)); 	

 #https://www.npmjs.com/package/ethereum-public-key-to-address
 Qx=0xe68acfc0253a10620dff706b0a1b1f1f5833ea3beb3bde2250d5f271f3563606;
 Qy=0x672ebc45e0b7ea2e816ecb70ca03137b1c9476eec63d4632e990020b7b6fba39;
 Qt=_CURVE_256K1([Qx,Qy]);
 hQt=FCL_ethereum_PubkeyToAddress(Qt);
 print("PubKey address;",hex(hQt));

 return 0;


