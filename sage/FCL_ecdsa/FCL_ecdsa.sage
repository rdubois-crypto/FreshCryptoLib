#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project 
#///* License: This software is licensed under MIT License 	 
#///* See LICENSE file at the root folder of the project.				 
#///* FILE: FCL_schnorr.sage						         
#///* 											 
#///* 											 
#///* DESCRIPTION: A x-only version of schnorr signature algorithme over elliptic curve
#///* This version is compatible with BIP340 and BIP Musig2
#///* It is an additive version of https://en.wikipedia.org/wiki/Schnorr_signature,
#///* (elliptic curve replaces prime field), with a tiny difference (substraction is performed
#///* by verifier)
#//**************************************************************************************/

from hashlib import *
#from sha3 import keccak_256

load('../FCL_common/FCL_bn_io.sage');
load('../FCL_common/FCL_elliptic.sage');

    
_BITCOIN_TAG ='BIP0340/challenge'
_ETHER_TAG =  'EIPXXXX/challenge'
#_ETHER_TAG2=keccak_256("\x19Ethereum Signed Message:\n32" + keccak_256(message))
_STARKNET_TAG='SIPXXXX/challenge'

_BITCOIN_HASH=sha256
#_ETHER_HASH=keccak_256
#_STARKNET_HASH=poseidon

#set global variables to bitcoin settings
_G_MODULUS=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
_G_BYTESIZE =32;
_G_HASH = _BITCOIN_HASH
_G_CURVE, _G_POINT = FCL_ec_Init_Curve(sec256k_p, sec256k_a, sec256k_b, sec256k_gx, sec256k_gy, sec256k_n);
_G_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;



#set global variables to stark settings
_G_CURVE, _G_POINT = FCL_ec_Init_Curve(stark_p, stark_a, stark_b, stark_gx, stark_gy, stark_n);
_G_ORDER=stark_n;

#set global variables to ethereum settings
#_G_HASH = _ETHER_HASH
_G_CURVE, _G_POINT = FCL_ec_Init_Curve(sec256k_p, sec256k_a, sec256k_b, sec256k_gx, sec256k_gy, sec256k_n);
_G_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
_G_ALPHA= 3; #3 is a non square for this prime field
Fq=GF(_G_ORDER);
Fp=_G_POINT[0].parent();


#set global variables to P256 settings
_G_CURVE, _G_POINT = FCL_ec_Init_Curve(sec256p_p, sec256p_a, sec256p_b, sec256p_gx, sec256p_gy, sec256p_n);
_G_ORDER=sec256p_n;


def FCL_Hash2sec256k1(h):
 Fp=GF(sec256k_p);
 y2=Fp(h^3+sec256k_b);
 if is_square(y2==false):
 	h=3*(h+sec256k_b);
 	print("swap");	
 	y2=h^3+sec256k_b;
 	print("square:",is_square(y2));
 return (h, y2);

def FCL_ecdsa_keygen(random_k):
 pk=random_k*_G_POINT;
 sk=random_k %_G_ORDER;
 
 return pk ,sk;


#//ecdsa signature with external randomness
def FCL_ecdsa_sign_hash(curve, G, hash, seckey, k):
 #e=H(msg)
 e=hash;
 kG=k*G;

 r=int(kG[0])% _G_ORDER;
 v=int(kG[1])&1;
 s=(e+r*seckey)*Fq(k)^-1   
 return int(r),int(s), v;
  
#//ecdsa signature with external randomness
def FCL_ecdsa_sign_core(curve, G, msg, seckey, k):
 #e=H(msg)
 e=int('0x'+_G_HASH(msg).hexdigest(),16);

 r,s,v=FCL_ecdsa_sign_hash(curve, G, e, seckey, k);
 
 return int(r),int(s), v;

#canonization bip, BIP62
def FCL_ecdsa_sign_canonic(curve, G, msg, seckey, k):
 r,s, v=FCL_ecdsa_sign_core(curve, G, msg, seckey, k);
 if s>(_G_ORDER>>1):
  s=_G_ORDER-s;
  v=1-v;
 return r,s, v;


def FCL_ecdsa_verify_core(curve, G, e, Q, r,s):

 sm1=Fq(s)^-1;
 u1=int(Fq(e)*sm1);
 u2=int(Fq(r)*sm1);
 kG=u1*G+u2*Q;
 print("\n u1",u1,"\n u2",u2);
 
 print("\n u1",hex(u1),"\n u2",hex(u2));
 rc=int(kG[0]);
 
 return (Fq(rc)==Fq(r));


def FCL_ecdsa_verify(curve, G, msg, Q, r,s):
 e=int('0x'+_G_HASH(msg).hexdigest(),16);
 
 return FCL_ecdsa_verify_core(curve, G, e, Q, r,s);
 

def FCL_ecdsa_recovery(curve, G, hash, parity,r,s):
 Fq=GF(_G_ORDER);
 if r>_G_ORDER:
  return False;
 if s>_G_ORDER: 
  return False;
 
 #in ethereum 
 y=FCL_ec_decompress(curve, r, parity);
 R=curve([r,y]);
 #print("rec R=",R);
 rm1=(Fq(r)^-1)
 u1=int(hash*(rm1*-1));
 u2=int(s*rm1);
 
 Q=G*u1+R*u2;
 	
 return Q;

_NTESTS=10;
	
def test_consistency():
 for i in [0.._NTESTS]: 
   pubkey, seckey=FCL_ecdsa_keygen(randint(0,_G_ORDER-1)); 
   k=randint(0,_G_ORDER-1);
   m=randint(0,_G_ORDER-1);
   msg=FCL_BN_to_bytes(m,_G_BYTESIZE);
   h=int('0x'+_G_HASH(msg).hexdigest(),16);
   r,s,v=FCL_ecdsa_sign_canonic(_G_CURVE, _G_POINT, msg, seckey, k);
   flag=FCL_ecdsa_verify(_G_CURVE, _G_POINT, msg, pubkey, r,s);
   Q=FCL_ecdsa_recovery(_G_CURVE, _G_POINT,h, v ,r,s  );
   if(flag!=True):
    print("Verification failed !!"); 
    return False;
   if(Q!=pubkey):
    print("Recovery failed !!"); 
    return False;
 
   print("\n Verification and Recovery OK");
 return True;

#https://starkscan.co/contract/0x053a2e69119c26977102dae51ba3e87e01e2c43161615aa5af73dd4483dbd73c
#https://starkscan.co/contract/0x053a2e69119c26977102dae51ba3e87e01e2c43161615aa5af73dd4483dbd73c#read-write-contract0x31d4839cf06868be8d891e486af2765f7e67acd8babaa087bc2d3b8ed9cc046
#FCL_ecdsa_verify(_G_CURVE, _G_POINT, 0x05f32d2947ac403194b1b788a5828f05b5ef89a577f72f71c33171c75900b8de, 0x31d4839cf06868be8d891e486af2765f7e67acd8babaa087bc2d3b8ed9cc046, 



