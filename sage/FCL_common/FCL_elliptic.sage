#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project */
#///* License: This software is licensed under MIT License 	 */
#///* See LICENSE file at the root folder of the project.				 */
#///* FILE: FCL_elliptic.sage						         */
#///* 											 */
#///* 											 */
#///* DESCRIPTION: precompute a 8 dimensional table for Shamir's trick from a public key
#///* 
#//**************************************************************************************/


def Init_Curve(curve_characteristic,curve_a, curve_b,Gx, Gy, curve_Order):    
	Fp=GF(curve_characteristic); 				#Initialize Prime field of Point
	Fq=GF(curve_Order);					#Initialize Prime field of scalars
	Curve=EllipticCurve(Fp, [curve_a, curve_b]);		#Initialize Elliptic curve
	curve_Generator=Curve([Gx, Gy]);
	
	return [Curve,curve_Generator];


#//initialize elliptic curve , short weierstrass form	
def FCL_ec_Init_Curve(curve_characteristic,curve_a, curve_b,Gx, Gy, curve_Order):    
	Fp=GF(curve_characteristic); 				#Initialize Prime field of Point
	Fq=GF(curve_Order);					#Initialize Prime field of scalars
	Curve=EllipticCurve(Fp, [curve_a, curve_b]);		#Initialize Elliptic curve
	curve_Generator=Curve([Gx, Gy]);
	
	return [Curve,curve_Generator];



#//decompress even y from x value to point (x,y)    
def FCL_ec_decompress_even(curve, pubkey_x):    
  y2=pubkey_x**3+ curve.a4()*pubkey_x+curve.a6();
  y=sqrt(y2);
  if (int(y)%2):
   return -y;
  return y;	
  

#//decompress even y from x value to point (x,y)    
def FCL_ec_decompress(curve, pubkey_x, parity):    
  y2=pubkey_x**3+ curve.a4()*pubkey_x+curve.a6();
  y=sqrt(y2);
 
  if ((int(y)%2)!=parity):
   return -y;
  return y;	
    

#//Curve secp256r1, aka p256	
#//curve prime field modulus
sec256p_p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;
#//short weierstrass first coefficient (=a4 in sage)
sec256p_a =0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC;
#//short weierstrass second coefficient    
sec256p_b =0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B;
#//generating point affine coordinates  (=a6 in sage)  
sec256p_gx =0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296;
sec256p_gy =0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5;
#//curve order (number of points)
sec256p_n =0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;    	



#//Curve secp256k1, aka bitcoin curve, aka ethereum curve
sec256k_p=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;
sec256k_a=0;
sec256k_b=7;
sec256k_gx=0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798 ;
sec256k_gy=0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8;
sec256k_n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;

#//Init_Curve(sec256p_p, sec256p_a, sec256p_b, sec256p_gx, sec256p_gy, sec256p_n);

stark_p=2^251+17*2^192+1     
stark_a=1;
stark_b=0x6f21413efbe40de150e596d72f7a8c5609ad26c15c915c1f4cdfcb99cee9e89;
stark_q=0x800000000000010ffffffffffffffffb781126dcae7b2321e66a241adc64d2f;
stark_gx = 0x1ef15c18599971b7beced415a40f0c7deacfd9b0d1819e03d723d8bc943cfca;
stark_gy = 0x5668060aa49730b7be4801df46ec62de53ecd11abe43a32873000c36e8dc1f;	
stark_n=0x800000000000010ffffffffffffffffb781126dcae7b2321e66a241adc64d2f;

#//https://github.com/bellesmarta/baby_jubjub is compliant with https://github.com/iden3/circomlibjs/blob/4f094c5be05c1f0210924a3ab204d8fd8da69f49/src/babyjub.js in non reduced Ted form
#//https://github.com/iden3/circomlib/blob/master/test/babyjub.js
#//https://github.com/iden3/circomlibjs/blob/4f094c5be05c1f0210924a3ab204d8fd8da69f49/test/eddsa.js
#//here it is a twisted edwards curve:https://hyperelliptic.org/EFD/g1p/auto-twisted.html
#//generate poseidon:https://github.com/iden3/circomlibjs/blob/main/src/poseidon_gencontract.js
babyjj_p=21888242871839275222246405745257275088548364400416034343698204186575808495617;
babyjj_n=21888242871839275222246405745257275088614511777268538073601725287587578984328;
babyjj_A=168700;
babyjj_D=168696;
#//https://github.com/bellesmarta/baby_jubjub, unreduced
babyjj_gx=995203441582195749578291179787384436505546430278305826713579947235728471134;
babyjj_gy=5472060717959818805561601436314318772137091100104008585924551046643952123905;


#nes_p=next_prime(sec256k_p);
#Fp=GF(nes_p); 	
#Curve=EllipticCurve(Fp, [sec256k_a, sec256k_b]);	
#order=Curve.order();
#while(is_prime(order)==false):
# nes_p=next_prime(nes_p);
# Fp=GF(nes_p); 	
# Curve=EllipticCurve(Fp, [sec256k_a, sec256k_b]);	
# order=Curve.order();
#print ("nes_p=",nes_p);

#nes_p= 115792089237316195423570985008687907853269984665640564039457584007908834744347;
#Fp=GF(nes_p); 	
#Curve=EllipticCurve(Fp, [sec256k_a, sec256k_b]);	
#q=Curve.order();
#Fq=GF(q); 	
#Curveq=EllipticCurve(Fq, [sec256k_a, sec256k_b]);	
#orderq=Curveq.order();








