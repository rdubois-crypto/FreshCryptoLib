#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project 
#///* License: This software is licensed under MIT License 	 
#///* See LICENSE file at the root folder of the project.
##/* 											 */
##/* 											 */
##/* DESCRIPTION: Pedersen Hash modified Stark version as coded in */
##/* https://github.com/starkware-libs/cairo-lang/blob/master/src/starkware/crypto/signature/fast_pedersen_hash.py
##/* spec is here:https://docs.starkware.co/starkex/pedersen-hash-function.html
##/source code 		 */
##/**************************************************************************************/

	
from FCL_starknet.FCL_pedersen_hash import *




#https://github.com/xJonathanLEI/starknet-rs/blob/master/starknet-crypto/src/pedersen_hash.rs
#a="03d937c035c878245caf64531a5756109c53068da139362728feb561405371cb",
#  b=          "0208a0a10250e382e1e4bbe2880906c2791bf6275695e02fbbc6aeff9cd8b31a",
#expected:030e480bed5fe53fa909cc0f8c4d99b8f9f2c016be4c41e13a4848797979c662
def test_vec():
 
 a= 0x03d937c035c878245caf64531a5756109c53068da139362728feb561405371cb;
 b= 0x0208a0a10250e382e1e4bbe2880906c2791bf6275695e02fbbc6aeff9cd8b31a;
 
 expected = 0x030e480bed5fe53fa909cc0f8c4d99b8f9f2c016be4c41e13a4848797979c662;
 res=false;
 
 print("\n *********************************** ");
 print("\n Test FCL sage starknet module: ");
 print("\n *********************************** ");
 print("   Unitary starknet-rs test from https://github.com/xJonathanLEI/starknet-rs/blob/master/starknet-crypto/src/pedersen_hash.rs: ");
 
# print("ped:",hex(pedersen(a,b)));
 
 if pedersen(a,b)== expected:
   print("OK");
   res=true;
 else:
  print("KO");
 
 return res;

def test_vec_long():
 tv=[0x049d36570d4e46f48e99674bd3fcc84644ddd6b96f7c741b1562b82f9e004dc7,0x02e269d930f6d7ab92b15ce8ff9f5e63709391617e3465fff79ba6baf278ce60,0x03d937c035c878245caf64531a5756109c53068da139362728feb561405371cb,    0x0208a0a10250e382e1e4bbe2880906c2791bf6275695e02fbbc6aeff9cd8b31a,0x030e480bed5fe53fa909cc0f8c4d99b8f9f2c016be4c41e13a4848797979c662];
 expected=0x01a152d943c2b63cc61a485bc292b1f1be01d863e660c3e5f8ff440e3b6cfb9f;
 
 res= pedersen_hash(tv,Integer(5));
 if res== expected:
   print("OK");
   return true;
 else:
  print("KO");
  
  return false;    

if __name__ == "__main__":

 arithmetic(False)
 test_vec();   
 test_vec_long();   


