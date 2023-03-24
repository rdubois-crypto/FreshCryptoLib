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
from sha3 import keccak_256

load('../FCL_common/FCL_bn_io.sage');
load('../FCL_common/FCL_elliptic.sage');
    
_BITCOIN_TAG ='BIP0340/challenge'
_ETHER_TAG =  'EIPXXXX/challenge'
_STARKNET_TAG='SIPXXXX/challenge'

_BITCOIN_HASH=sha256
_ETHER_HASH=keccak_256
#_STARKNET_HASH=poseidon

#set global variables to bitcoin settings
_G_BYTESIZE =32;
_G_HASH = _BITCOIN_HASH
_G_CURVE, _G_POINT = FCL_ec_Init_Curve(sec256k_p, sec256k_a, sec256k_b, sec256k_gx, sec256k_gy, sec256k_n);
_G_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;
#set global variables to ethereum settings
_G_HASH = _ETHER_HASH
_G_CURVE, _G_POINT = FCL_ec_Init_Curve(sec256k_p, sec256k_a, sec256k_b, sec256k_gx, sec256k_gy, sec256k_n);
_G_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;


Fq=GF(_G_ORDER);

#The function hashtag(x) where tag is a UTF-8 encoded tag name and x is a byte array 
# returns the 32-byte hash SHA256(SHA256(tag) || SHA256(tag) || x).

#bip340 like hashing, note that the double concatenation of tag is 
def hashtag(tag, hf, msg):
    tag_hash = hf(tag.encode()).digest()
    return int('0x'+hf(tag_hash + tag_hash + msg).hexdigest(),16);
    
#def hashtag_skeccak(tag: str, msg: bytes) -> bytes:
#def hashtag_poseidon


	
#print( _G_HASH(FCL_BN_to_bytes(0x616263,3)).hexdigest());

def FCL_schnorr_keygen(random_k):
 pk=random_k*_G_POINT;
 sk=random_k %_G_ORDER;
 if( int(pk[1])%2 !=0):
   pk=-pk;
   sk=_G_ORDER -sk;
 return pk[0],sk;
 
#//schnorr signature with external randomness
def FCL_schnorr_sign(curve, G, msg, seckey, k):
  r=k*G;
  if( int(r[1])%2!=0):
   r=-r;
   
  if(int(r[1])%2):
    r=curve(r[0], int(-r[1]));
 
  pubkey_x= (seckey*G)[0];
  m=FCL_BN_to_bytes(int(r[0]), _G_BYTESIZE)+FCL_BN_to_bytes(int(pubkey_x), _G_BYTESIZE)+msg;
  e=   hashtag(_BITCOIN_TAG, _BITCOIN_HASH,  m);
  s= (k+seckey*e)  %_G_ORDER;
  
  return int(r[0]),s,e;



#// take as input directly the correct hash (to mimic contract)
def FCL_schnorr_verify_fromhash(curve, G,  msg, pubkey_x, s,e):
    y= FCL_ec_decompress_even(curve, pubkey_x);
    Pub=curve([pubkey_x,y]);
    r=s*G-e*Pub;
        
    m=FCL_BN_to_bytes(int(r[0]), _G_BYTESIZE)+FCL_BN_to_bytes(int(pubkey_x), _G_BYTESIZE)+msg;
  
    e_r=   hashtag(_BITCOIN_TAG, _BITCOIN_HASH, m);
    if( e_r != e):
      return False;
    return True;      


def FCL_schnorr_verify_bip340(curve, G, msg, pubkey_x, s,rx):
    y= FCL_ec_decompress_even(curve, pubkey_x);
    Pub=curve([pubkey_x,y]);
        
    m=FCL_BN_to_bytes(int(rx), _G_BYTESIZE)+FCL_BN_to_bytes(int(pubkey_x), _G_BYTESIZE)+msg;
    
    y= FCL_ec_decompress_even(curve, pubkey_x);
    Pub=curve([pubkey_x,y]);
   
    e_r=   hashtag(_BITCOIN_TAG, _BITCOIN_HASH, m);
    
    r=s*G-e_r*Pub;
    
    if( e_r != e):
      return False;
   
    return True;      




pubkey_x, seckey=FCL_schnorr_keygen(randint(0,_G_ORDER-1));
k=randint(0,_G_ORDER-1);
msg=FCL_BN_to_bytes(randint(0,_G_ORDER-1),_G_BYTESIZE);
r,s,e=FCL_schnorr_sign(_G_CURVE, _G_POINT, msg, seckey, k);
flag=FCL_schnorr_verify_fromhash(_G_CURVE, _G_POINT, msg, pubkey_x, s,e);
print("verif from hash=",flag);
flag=FCL_schnorr_verify_bip340(_G_CURVE, _G_POINT, msg, pubkey_x, s,r);
print("verif from point",flag);



quit;



# Schnorr with short hash, as proven secure in https://eprint.iacr.org/2019/1105.pdf
def FCL_shortSchnorr_verify():
    return False;


