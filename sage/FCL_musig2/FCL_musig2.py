#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project 
#///* License: This software is licensed under MIT License 	 
#///* See LICENSE file at the root folder of the project.
##/* FILE: musig2.sage							             	  */
##/* 											  */
##/* 											  */
##/* DESCRIPTION: 2 round_multisignature Musig2 signatures verification sage simulation*/
##/* This is a high level simulation for validation purpose				  */
##/* https:##eprint.iacr.org/2020/1261.pdf             				  */
##/* note that some constant aggregating values could be precomputed			 
##***********************************************************/


#import hash function and curve parameters
from FCL_starknet.FCL_pedersen_hash import *


_MU=2;nb_users=4; size_message=1;seed=3;


# note: as much as possible, loop variable i is used to adress user in [0..n-1], 
# while j adress the position in a vector (of dimension n) of values .
_VERBOSE=false;

##***********************IO functions**********************************************/
def StringToHex(String):
	return String.encode('utf-8').hex();

def CairoDeclare_from_Point(comment, Point):
	if(_VERBOSE):
	  print("local ", comment,":(felt, felt)=(", hex(Point[0]), ",",hex(Point[1]),");");
	return 0;

def CairoDeclare_from_sacalar(comment, Scalar):
	if(_VERBOSE):
	  print("local ", comment,":felt=", hex(Scalar), ";");
	return 0;
#seed=3;
#set a fixed/variable seed for debug/generation purposes
set_random_seed(seed)
#print("set seed with", seed)
	
##***********************hash functions*******************************************/
# The following constants are used for domain separation of the H_agg, H_nonce and H_sig functions
# (No oracle from one can provide oracle for the other)
_SEPARATION_AGG=3
_SEPARATION_NON=2
_SEPARATION_SIG=1

		
### H_aggregate = pedersen_hash_state(L||X)
# expected public keys format is (int, int)
def H_agg(Coefficients, n_users, Xx, Xy, curve_order):
	Fq=GF(curve_order);
	#concatenate L with X, each public key is represented as a couple of felt
	Input=[_SEPARATION_AGG];			#separating the domains
	Input+=Coefficients+[Xx, Xy];			#Append X to L
	Hagg=pedersen_hash(Input, 2*(n_users+1)+1);	#H(L||X)
	return int(Fq(Hagg));

### H_Nonce = pedersen_hash_state(Xtilde||Ris||m)
#Xtilde is the public key aggregation computed at first round
def H_non(KeyAgg, nb_users, vec_R, message, size_message, curve_order, _MU):
	#concatenate L with X, each public key is represented as a couple of felt
	Input_hash=[_SEPARATION_NON];#separating the domains
	Input_hash+=[int(KeyAgg[0]), int(KeyAgg[1])];			#Append X, Ris and m
	Fq=GF(curve_order);
	
	for j in range(0,_MU):
		Input_hash=Input_hash+[int(vec_R[j][0]), int(vec_R[j][1])];
		
	for i in range(0,size_message):
		Input_hash=Input_hash+[message[i]];
	
	Hnon=pedersen_hash(Input_hash, 2*(_MU+1)+size_message);	#H(L||X)
	
	return int(Fq(Hnon));

def H_sig(KeyAgg, R, m, size_message, curve_order):
	Fq=GF(curve_order);
	Input=[_SEPARATION_SIG];#separating the domains
	Input+=[int(KeyAgg[0])];
	Input+=[int(R[0])];
	for cpt_i in range(0,size_message):
		Input=Input+[m[cpt_i]];
		
	#print("Input to pedersen=",Input, "size=", len(Input));	
	Hsig=pedersen_hash(Input, len(Input));
	
	return int(Fq(Hsig));


##***********************Aggregator functions****************************************/

###Computation of sum Pki^ai, ai=H(L, X_i)
# expected public keys format is (int, int)
def Musig2_KeyAgg(Curve, curve_generator, L, nb_users, curve_order):
	sum=0*curve_generator;				# initiate sum at infty_point
	for i in range(0,nb_users):
		ai=H_agg(L, nb_users, L[2*i], L[2*i+1], curve_order);	#i-th ai
		sum=sum+ai*Curve(L[2*i], L[2*i+1]);	#Xi^ai, where Xi is the ith public key
	return sum;

#input is n_users vector of size v
def Musig2_Sig1Agg(vec_Ri, nb_users, _MU):
	infty_point=0*curve_Generator;
	Aggregated_Ri=[0 for i in range(0,_MU)];
	
	for j in range(0,_MU):
		Aggregated_Ri[j]=infty_point;
	
	for i in range(0,nb_users):
		for j in range(0,_MU):#sum the contribution to previous ones
			Aggregated_Ri[j]+=vec_Ri[i][j];	
	return Aggregated_Ri;

def Musig2_Sig2Agg(vec_s, curve_order, nb_users):
	Fq=GF(curve_order);
	s=Fq(0);
	for i in range(0,nb_users):
		s+=vec_s[i];
	
	return int(s);

##***********************Single user functions***************************************/

### Key generation
def Musig2_KeyGen(Curve, curve_generator, curve_order):
	Fq=GF(curve_order);
	privatekey=int(Fq.random_element());
	publickey=curve_generator*privatekey;
	#/* compatibility BIP340*/
	if( (int(publickey[1])&1)==1):#BIP340 compatibility
		privatekey=curve_order-privatekey;
		publickey=curve_generator*privatekey;		
			
	return [int(privatekey), publickey];
	
	
### Round1 from single signer view
# expected public key as a point
# in a practical version, random element shall be replaced by rfc6979 adaptation
def Musig2_Sign_Round1(curve_order, n_users,curve_Generator, _MU):
	Fq=GF(curve_order);	
	
	Ri= [0 for i in range(0,_MU)] ;
	
	nonces=[0 for i in range(0,n_users)];
	
	for j in range(0,_MU):
		nonces[j]=int(Fq.random_element());
		Ri[j]=nonces[j]*curve_Generator;     #Rij=rij*G;
	
	return [nonces, Ri];	
	
### Round2 from single signer view
def Musig2_Sign_Round2_all( 
		Curve, curve_generator,curve_order, nb_users, KeyAgg,#common parameters
		ai, privatekey_xi,					#user's data
		vec_R, vec_nonces,					#First round output
		message,  message_feltlength, _MU):			#input message
	Fq=GF(curve_order);
	b=H_non(KeyAgg, nb_users, vec_R, message,  message_feltlength, curve_order, _MU);
	
	R=0*curve_generator;
	for j in range(0,_MU):
		R=R+((b**(j))*vec_R[j]) ;
	c=(H_sig(KeyAgg, R, message, message_feltlength, curve_order));
	s=(c*ai*privatekey_xi);
	for j in range(0,_MU):
		s=s+(vec_nonces[j])*(Fq(b)**j);
	return [R,s, c];


##***********************Verifier functions******************************************/
def Musig_Verif_Core(Curve, curve_Generator,R,s,KeyAgg, c):
	#print("\s=", s, "\n c=", c);
	
	Gpows=s*curve_Generator;
	Xpowc=c*KeyAgg;
	#print("\nGpows=", Gpows, "\nXpowc=", Xpowc, "\nR+Xpowc=",R+Xpowc);
	
	return (Gpows==R+Xpowc);
	
#stark curve parameters from https://docs.starkware.co/starkex/stark-curve.html






