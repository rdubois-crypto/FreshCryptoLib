#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project 
#///* License: This software is licensed under MIT License 	 
#///* See LICENSE file at the root folder of the project.				 
#///* FILE: FCL_test_ecdsa.sage						         
#///* 											 
#///* 											 
#///* DESCRIPTION: A x-only version of schnorr signature algorithme over elliptic curve
#///* This version is compatible with BIP340 and BIP Musig2
#///* It is an additive version of https://en.wikipedia.org/wiki/Schnorr_signature,
#///* (elliptic curve replaces prime field), with a tiny difference (substraction is performed
#///* by verifier)
#//**************************************************************************************/
import sys
from hashlib import *

load('../FCL_common/FCL_bn_io.sage');
load('../FCL_common/FCL_elliptic.sage');
load('../FCL_ecdsa/FCL_ecdsa.sage');

#vectors from wychproof
key=_G_CURVE( [ 18614955573315897657680976650685450080931919913269223958732452353593824192568 , 90223116347859880166570198725387569567414254547569925327988539833150573990206]);
 
test1=[" 0_signature_malleability ",  0xbb5a52f42f9c9261ed4361f59422a1e30036e7c32b270c8807a419feca605023 , 19738613187745101558623338726804762177711919211234071563652772152683725073944 ,  34753961278895633991577816754222591531863837041401341770838584739693604822390];

test2=[" 2_valid ", 0xbb5a52f42f9c9261ed4361f59422a1e30036e7c32b270c8807a419feca605023 , 19738613187745101558623338726804762177711919211234071563652772152683725073944 ,  81038127931460614771119630195184981998133118182734418571583674321374907221979];
#test_consistency();

#print(FCL_ecdsa_verify_core(_G_CURVE, _G_POINT, test1[1], key, test1[2], test1[3] ));

#print(FCL_ecdsa_verify_core(_G_CURVE, _G_POINT, test1[1], key, test1[2], _G_ORDER-test1[3] ));

#print(FCL_ecdsa_verify_core(_G_CURVE, _G_POINT, test2[1], key, test2[2], test2[3] ));

def gen_vec():
   title="sage random";
   _NVECS=57;
   pubkey, seckey=FCL_ecdsa_keygen(randint(0,_G_ORDER-1)); 
   print("{\n \"keyx\":", int(pubkey[0]), ",\"keyy\":", int(pubkey[1]), ",")

   for i in [1.._NVECS]: 
    
  
    k=randint(0,_G_ORDER-1);
    m=randint(0,_G_ORDER-1);
    msg=FCL_BN_to_bytes(m,_G_BYTESIZE);
    h=int('0x'+_G_HASH(msg).hexdigest(),16);
    r,s,v=FCL_ecdsa_sign_canonic(_G_CURVE, _G_POINT, msg, seckey, k);
    print(" \"test_"+str(i)+"\":\"", title, "\", \"msg_"+str(i)+"\": \""+hex(h)+"\"",", \"sigx_"+str(i)+"\":", r, ", \"sigy_"+str(i)+"\":", s,",");
   
   
   print("\"NumberOfTests:\"",i,"\n}");
             
gen_vec();
