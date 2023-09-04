##*************************************************************************************/
##/* Copyright (C) 2022 - Renaud Dubois - This file is part of cairo_musig2 project	 */
##/* License: This software is licensed under a dual BSD and GPL v2 license. 	 */
##/* See LICENSE file at the root folder of the project.				 */
##/* FILE: musig2.sage							             	  */
##/* 											  */
##/* 											  */
##/* DESCRIPTION: 2 round_multisignature Musig2 signatures verification sage simulation*/
##/* This is a high level simulation for validation purpose				  */
##/* https:##eprint.iacr.org/2020/1261.pdf             				  */
##/* note that some constant aggregating values could be precomputed			  */
##**************************************************************************************/
## Stark curve parameters extracted from
## https:github.com/starkware-libs/cairo-lang/blob/master/src/starkware/cairo/common/ec.cairo

load('../FCL_starknet/FCL_pedersen_hash.sage');

# note: as much as possible, loop variable i is used to adress user in [0..n-1], 
# while j adress the position in a vector (of dimension n) of values .

##***********************IO functions**********************************************/
def StringToHex(String):
	return String.encode('utf-8').hex();

def CairoDeclare_from_Point(comment, Point):
	print("local ", comment,":(felt, felt)=(", hex(Point[0]), ",",hex(Point[1]),");");
	return 0;

def CairoDeclare_from_sacalar(comment, Scalar):
	print("local ", comment,":felt=", hex(Scalar), ";");
	return 0;

#set a fixed/variable seed for debug/generation purposes
set_random_seed(seed)
print("set seed with", seed)
	
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
def H_non(KeyAgg, nb_users, vec_R, message, size_message, curve_order):
	#concatenate L with X, each public key is represented as a couple of felt
	Input_hash=[_SEPARATION_NON];#separating the domains
	Input_hash+=[int(KeyAgg[0]), int(KeyAgg[1])];			#Append X, Ris and m
	Fq=GF(curve_order);
	
	for j in [0.._MU-1]:
		Input_hash=Input_hash+[int(vec_R[j][0]), int(vec_R[j][1])];
		
	for i in [0..size_message-1]:
		Input_hash=Input_hash+[message[i]];
	
	Hnon=pedersen_hash(Input_hash, 2*(_MU+1)+size_message);	#H(L||X)
	
	return int(Fq(Hnon));

def H_sig(KeyAgg, R, m, size_message, curve_order):
	Fq=GF(curve_order);
	Input=[_SEPARATION_SIG];#separating the domains
	Input+=[int(KeyAgg[0])];
	Input+=[int(R[0])];
	for cpt_i in [0..size_message-1]:
		Input=Input+[m[cpt_i]];
		
	print("Input to pedersen=",Input, "size=", len(Input));	
	Hsig=pedersen_hash(Input, len(Input));
	
	return int(Fq(Hsig));


##***********************Aggregator functions****************************************/

###Computation of sum Pki^ai, ai=H(L, X_i)
# expected public keys format is (int, int)
def Musig2_KeyAgg(Curve, curve_generator, L, n_users, curve_order):
	sum=0*curve_generator;				# initiate sum at infty_point
	for i in [0..n_users-1]:
		ai=H_agg(L, nb_users, L[2*i], L[2*i+1], curve_order);	#i-th ai
		sum=sum+ai*Curve(L[2*i], L[2*i+1]);	#Xi^ai, where Xi is the ith public key
	return sum;

#input is n_users vector of size v
def Musig2_Sig1Agg(vec_Ri):
	infty_point=0*curve_Generator;
	Aggregated_Ri=[0.._MU-1];
	
	for j in [0.._MU-1]:
		Aggregated_Ri[j]=infty_point;
	
	for i in [0..nb_users-1]:
		for j in [0.._MU-1]:#sum the contribution to previous ones
			Aggregated_Ri[j]+=vec_Ri[i][j];	
	return Aggregated_Ri;
	

def Musig2_Sig2Agg(vec_s):
	s=Fq(0);
	for i in [0..nb_users-1]:
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
def Musig2_Sign_Round1(curve_order, n_users,curve_Generator):
	Fq=GF(curve_order);	
	Ri=[0.._MU-1];
	nonces=[0..n_users-1];
	
	for j in [0.._MU-1]:
		nonces[j]=int(Fq.random_element());
		Ri[j]=nonces[j]*curve_Generator;     #Rij=rij*G;
	
	return [nonces, Ri];	
	
### Round2 from single signer view
def Musig2_Sign_Round2_all( 
		Curve, curve_generator,curve_order, nb_users, KeyAgg,#common parameters
		ai, privatekey_xi,					#user's data
		vec_R, vec_nonces,					#First round output
		message,  message_feltlength):			#input message
	Fq=GF(curve_order);
	b=H_non(KeyAgg, nb_users, vec_R, message,  message_feltlength, curve_order);
	
	R=0*curve_generator;
	for j in [0.._MU-1]:
		R=R+(b^(j)*vec_R[j]) ;
	c=(H_sig(KeyAgg, R, message, message_feltlength, curve_order));
	s=(c*ai*privatekey_xi);
	for j in [0.._MU-1]:
		s=s+(vec_nonces[j])*Fq(b)^j;
	return [R,s, c];


##***********************Verifier functions******************************************/
def Musig_Verif_Core(Curve, curve_Generator,R,s,KeyAgg, c):
	Gpows=s*curve_Generator;
	Xpowc=c*KeyAgg;
	#("\nGpows=", Gpows, "\nXpowc=", Xpowc, "\nR+Xpowc=",R+Xpowc);
	
	return (Gpows==R+Xpowc);
	
#stark curve parameters from https://docs.starkware.co/starkex/stark-curve.html

curve_characteristic=2^251+17*2^192+1     
is_prime(curve_characteristic); #should return true
beta = 0x6f21413efbe40de150e596d72f7a8c5609ad26c15c915c1f4cdfcb99cee9e89
Stark_order=0x800000000000010ffffffffffffffffb781126dcae7b2321e66a241adc64d2f
GEN_X = 0x1ef15c18599971b7beced415a40f0c7deacfd9b0d1819e03d723d8bc943cfca;
GEN_Y = 0x5668060aa49730b7be4801df46ec62de53ecd11abe43a32873000c36e8dc1f;	
[Curve,curve_Generator, P0, P1, P2, P3, Shift]=Init_Stark(curve_characteristic,1, beta,GEN_X, GEN_Y,Stark_order) ;






