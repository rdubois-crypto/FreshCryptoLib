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
#///* DESCRIPTION: generating additional edge cases random test vectors for ecdsa
#//**************************************************************************************/
import sys
from hashlib import *


load('../FCL_common/FCL_bn_io.sage');
load('../FCL_common/FCL_elliptic.sage');
load('../FCL_ecdsa/FCL_ecdsa.sage');


def gen_vec_edge():
   title="ecmulmul edge case random";
   _NVECS=57;
   pubkey, seckey=FCL_ecdsa_keygen(randint(0,_G_ORDER-1)); 
   print("{\n \"keyx\":", int(pubkey[0]), ",\"keyy\":", int(pubkey[1]), ",")
   Fq=GF(_G_ORDER);

   for i in [1.._NVECS]: 
    
    #msg=FCL_BN_to_bytes(m,_G_BYTESIZE);
    u2=randint(0,_G_ORDER-1)&(2^(8*_G_BYTESIZE)-3);#u2Q=u2.sec.P
    
    #/*choose u1P+u2Q==2P, so that last addition is a doubling: u1=(2-u2.sec) [n]*/ 
    
    u1=Fq(2)-Fq(u2*seckey);
    print("\n i",i,"u1=",hex(u1), "u2=",hex(u2));
    
    u1G=_G_POINT*int(u1);
    u2Q=pubkey*u2;
    
    print(" sub:",u1G+u2Q-2*_G_POINT);
    
    #h=int('0x'+_G_HASH(msg).hexdigest(),16);
    #r,s,v=FCL_ecdsa_sign_canonic(_G_CURVE, _G_POINT, msg, seckey, k);
    #print(" \"test_"+str(i)+"\":\"", title, "\", \"msg_"+str(i)+"\": \""+hex(h)+"\"",", \"sigx_"+str(i)+"\":", r, ", \"sigy_"+str(i)+"\":", s,",");
   
   
   print("\"NumberOfTests:\"",i,"\n}");
             
gen_vec_edge();
