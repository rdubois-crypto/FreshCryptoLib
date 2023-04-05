#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project 
#///* License: This software is licensed under MIT License 	 
#///* See LICENSE file at the root folder of the project.				 
#///* FILE: FCL_ethereum_hack.sage						         
#///* 											 
#///* 											 
#///* DESCRIPTION: implement ethereum hacky use of ecrecover
#///* https://ethresear.ch/t/you-can-kinda-abuse-ecrecover-to-do-ecmul-in-secp256k1-today/2384
#//**************************************************************************************/

import hashlib, secrets


load('../FCL_common/FCL_elliptic.sage');
load('../FCL_common/FCL_bn_io.sage');
load('FCL_ethereum.sage');
    
from sha3 import keccak_256

# simulating tools
def FCL_hackmul_Point(k, G):
  v=27+(int(G[1])%2);
  r=int(G[0]);
  s=(r*k)%_256K1_ORDER ; 
  return FCL_ecRecoverPoint_from_hash(0, v,r,s);

 
# compute sG+eQ, r=Qx,  hash=-rs, v=v, r=r, s=-er 
def FCL_hackmulmuladd_Point(s,  e, Q):
  v=27+(int(Q[1])%2);
  r=int(Q[0]);
  hash=_256K1_ORDER-((r*s)%_256K1_ORDER);
  s=_256K1_ORDER-((e*r)%_256K1_ORDER);
  
  return FCL_ecRecoverPoint_from_hash(hash, v,r,s);
 
 
#unfortunately precompiles return hash 
def FCL_hackmul(k, G):
  v=27+int(G[1])%2;
  r=G[0];
  s=(r*k)%_256K1_ORDER ; 
  return FCL_ecRecover(0, v,r,s);
  	
# compute sG+eQ, r=Qx,  hash=-rs, v=v, r=r, s=-er 
def FCL_hackmulmuladd(s, G, e, Q):
  v=27+(int(G[1])%2);
  r=G[0];
  hash=(-r*s)%_256K1_ORDER;
  s=(-e*s)%_256K1_ORDER;
  
  return FCL_ecRecover(hash, v,r,s);
  
def test_hacks():
 flag=True;
 
 k=randint(0,_256K1_ORDER-1);
 P1=FCL_hackmul_Point(k, _G_256K1);
 P2=k*_G_256K1;
 
 print("flag:",P1==P2);
 Q=P1;
 e=randint(0,_256K1_ORDER-1);

 P1=FCL_hackmulmuladd_Point(k,e,Q);
 P2=k*_G_256K1-e*Q;
 print("flag:",P1==P2);

 
 return flag;  
  
  
  
