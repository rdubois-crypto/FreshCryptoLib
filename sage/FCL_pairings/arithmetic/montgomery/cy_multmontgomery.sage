#/**********************************************************************************/
#/* Copyright (C) 2022 - Renaud Dubois
#/* This file is part of cy_lib project
#/* License: This software is licensed under a dual BSD and GPL v2 license.
##/**********************************************************************************/

#/*  Object : 	CYLIB_SAGE 		*/
#/*  Langage : 	sage	 		*/
#/*  Description :The CYLIB_SAGE is a simulation of the CYLIB library, designed for simulation and testing of the CYLIB.
#/*  FILE : cy_multmontgomery.sage is a simulation of a generic size Montgomery multiplier

	
load('../common/io_conversions.sage');
	
 #/***************** Montgomery CIOS multiplication *******/
 
def MontGomery_Mult(a,b,p,n0, r, size):
	tu64_carry=0;
	l_carry=0;
	tn=0
	u=0
	
	mask = (2^g_size)-1
	ze=tn&mask
	r=copy(matrix(1, size+1)[0])
	
	for i in [0..size-1]:
 		tu64_carry=0;
 		l_carry=0;
 		
 		#tn=r[j]+a[j]*b[i]+carry
 		print("\n *********************** i=",i);
 		for j in [0..size-1]:
 			print("a=",hex(a[j]), "b", hex(b[i]), "carry", hex(l_carry), "r[j]",hex(r[j]))
 			tn=MUL_MWORD(tn, b[i], a[j]); 	
 			tn=ADDTO_DWORD(tn, tu64_carry, r[j]);
 			print("\n i",i,"j=",j, "tn=",hex(tn));	
 			tn=ADDTO_DWORD(tn, tu64_carry, l_carry);
 			print("tn=",hex(tn))

 			print("\n tn0:", hex(tn&mask));
 			print("\n tn1:", hex(tn>>g_size));
 			
 			r[j]=COPY_MWORD(r[j],tn&mask);
 			l_carry=COPY_MWORD(l_carry,tn>>g_size);
 			print("\n lcarry=",hex(l_carry))
 		r[size]=ADDCTO_MWORD(r[size], tu64_carry, l_carry);# r[size]+=carry;//devil carry here
 		print("\n r loop1: ");
 		PRINT_HEXATAB(r, size+1);
 		temp=tu64_carry;
 		print("\n temp=",hex(temp))
 		u=WORD_MUL_LOW(u, r[0], n0, mask);
 		print(" (u,r0,n0)=",hex(u), hex(r[0]), hex(n0))
 		tn=MUL_MWORD(tn, u, p[0]); 
 		tn=ADDTO_DWORD(tn, tu64_carry, r[0]);
 		l_carry=COPY_MWORD(l_carry,tn>>g_size);
        
 		
 		print("\n begin loop2 with (u,tn)=",hex(u), hex(tn));
        #/*R=R+u[i]*Modulus*/
 		for j in [1..size-1]:	
        	#/*(tn[0],tn[1])=r[j]+u*Modulus[j+1]+carry */
 			tn=MUL_MWORD(tn, u, p[j]);
 			print("tn=",hex(tn), "lcarry=", hex(l_carry), "r[j]=",hex(r[j]))
 			tn=ADDTO_DWORD(tn, tu64_carry, r[j] );
 			print("tn=",hex(tn))
 			tn=ADDTO_DWORD(tn, tu64_carry, l_carry);
 			print("tn=",hex(tn))
 
 			print("\n i= j=",i,j);
 			print("\n tn0:", hex(tn&mask));
 			print("\n tn1:", hex(tn>>g_size));
 			
           	#/*Allocating and shifting */
 			r[j-1]=COPY_MWORD(r[j-1],tn&mask);
 			l_carry=COPY_MWORD(l_carry,tn>>g_size);
        #endfor
 		tn=COPY_MWORD(tn, r[size]);
 		tn=ADDTO_DWORD(tn,tu64_carry,  l_carry );
 		r[size-1]=COPY_MWORD((r)[size-1],tn&mask);#//(r)[size-1]=tn[0];
 		print("\n r loop2: ", r);

 		r[size]=ZERO_MWORD((r)[size]);
 		r[size]=ADDTO_MWORD((r)[size], tu64_carry, temp);
 		r[size]=ADDTO_MWORD((r)[size], tu64_carry, tn>>g_size);
 		print("\n r[size]:", r[size]);
 		PRINT_HEXATAB(r, size+1);
	return r;

####################################		
def test_modulus(n, size, a, b):
	Fp=FiniteField(n)
	tn=ceil(ceil(log(n)/log(2))/size)
	print("tn=",tn)
	R=(2^size)^tn;
	n0=mod(-n, 2^size)^-1
	Rm1 = Fp(2^size)^-1
	
	print("a=",(Fp(a)))
	print("aR=",hex(Fp(a*R)))
	
	print("b=",(Fp(b)))
	print("bR=",hex(Fp(b*R)))
	
	print("n0=", n0, hex(n0))
	print("a*b*R mod n=", Fp(a*b*R), hex(Fp(a*b*R)))
	print("a*b mod n=", Fp(a*b), hex(Fp(a*b)))
	
	return 0;		
		
#################################################################################################################################
#/***************** Tests *******/
#################################################################################################################################


#vector tests for nsize 128, wordsize128

#print("###########################################  vector tests for nsize 128, wordsize128)

g_size = 128

#test_modulus(n,128,a, b)
n0_expected=0xf31c335364fbed41065a0cd88ff7a41

#mod_aR=Conv_word(aR,128)
#mod_bR=Conv_word(bR, 128)
#mod_n=Conv_word(n, 128)
#r=MontGomery_Mult(mod_aR, mod_bR,mod_n,n0_expected, r, 2)
#PRINT_HEXATAB(r,3)
#mod_aR=MontGomery_Mult(mod_bR, r,mod_n,n0_expected, mod_aR, 2)
#PRINT_HEXATAB(mod_aR,3)
#mod_bR=MontGomery_Mult(mod_aR, r,mod_n,n0_expected, mod_bR, 2)
#PRINT_HEXATAB(mod_bR,3)



#vector test from    https://www.researchgate.net/publication/4107322_Montgomery_modular_multiplication_architecture_for_public_key_cryptosystems/link/58778eda08ae6eb871d152b5/download 
#n =       0xeee74404d129949520704c5bf5814703
#aR = 	  0x223520375bd184b2bac64c9d1a6c55fa
#bR = 	  0xd857ed0720d590d61f05c150e1e40917
#expected= 0xe3bd635debc8021ea0208d75df078ea6


#sec256k1 modulus
n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F

aR = 0x329efd7b3f30d8b2c2476c804506484dfcd150436a0cbf871ad7a1c78e2ac345
bR = 0x49e541bf948579eb4e92ebf5e8f0f80f184c6db04df8670013dde86d17cdf296
expected =0xf3d1ecf7de8a92956ee05242b9e6342d8a4be8f1965c3693686e0ba06791184d
a=0xac373f3fbdf8edff670873cee395a935329efd7b3f30d8b2c2476c7f98cf067d
b=0xc1554db541bedcf0648f35904d527cb649e541bf948579eb4e92ebf5279ba778

size=floor(log(2^128)/log(2))  

Fp=FiniteField(n)

R=2^g_size;
n0=mod(-n, 2^g_size)^-1
Rm1 = Fp(2^g_size)^-1

Mr=Fp(aR*bR*Rm1)

vec_a=Conv_word_MSB(aR, 8)
vec_b=Conv_word_MSB(bR, g_size)
vec_p=Conv_word_MSB(n, g_size)
r=O

hex(a*Rm1*b*Rm1*R);


#MontGomery_Mult(vec_a,vec_b,vec_p , n0, r, 2)

#import subprocess

#a=subprocess.check_output(["./test"])

Conv_word_MSB(aR, 8)
Conv_word_MSB(bR, 8)
Conv_word_MSB(expected, 8)






