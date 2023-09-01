#/**********************************************************************************/
#/* Copyright (C) 2022 - This file is part of cy_lib project				 	   */
#/* License: This software is licensed under a dual BSD and GPL v2 license. 	   */
#/* See LICENSE file at the root folder of the project.							   */
#/* FILE: precalc_montg.sage										   		       */
#/* 																			   */
#/* 																			   */
#/* DESCRIPTION: computation of Montgomery constants							   */
#/**********************************************************************************/

load('../common/io_conversions.sage');

####################################		
def precompute_montgomery_constants(name, n, sizeword, sizecurve):
	print("Montgomery constants for",name)
	Fp=FiniteField(n)
	tn=ceil(sizecurve/sizeword)
	print("tn=",tn)
	R=(2^sizecurve);
	n0=mod(-n, 2^sizeword)^-1
	Rmodp = (Fp(2)^sizecurve)
	
	print("n0=", n0, hex(n0))
	print("2^size mod p =", Rmodp, hex(Rmodp))
	#note that this constant is used in bolos, we can deduce that a montgomery reduction instead of an interleaved multiplication is used
	print("2^(2*size mod p) =", Rmodp^2, hex(Rmodp^2))
		
	
	return [lift(n0), lift(Rmodp), lift(Rmodp^2)];		
	

#/*********** Write as machine words ************/ 
def Conv_word(A, size_word):
	sizeA=ceil(log(A)/log(2))  
	sizeM=ceil(sizeA/size_word)
	M=copy(matrix(1, sizeM)[0])
	
	mask = 2^size_word-1
	print("{");
	for i in [0..sizeM-1]:
		M[i]= A 	& mask;
		print(" ",hex(M[i]),",");
		A = A >> size_word
	print("};");
	return M;	
	
	
#sec256k1 modulus
n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
name="p_256k1"
#precompute for a 128 bits emulation by 64 bit words
res=precompute_montgomery_constants(name, n, 128, 256);
#Conv_word_MSB(res[0], 8);
#Conv_word_MSB(res[1], 8);

a=0xac373f3fbdf8edff670873cee395a935329efd7b3f30d8b2c2476c7f98cf067d;
b=0xc1554db541bedcf0648f35904d527cb649e541bf948579eb4e92ebf5279ba778;
Fp=FiniteField(n)

print("aR=", hex(Fp(a*res[1])));


#sec256r1 modulus
#n=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
#name="p_256r1"
#precompute for a 128 bits emulation by 64 bit words
#res=precompute_montgomery_constants(name, n, 128, 256);
#Conv_word_MSB(res[0], 8);
#Conv_word_MSB(res[1], 8);

#Fp=FiniteField(n)

#a=0xac373f3fbdf8edff670873cee395a935329efd7b3f30d8b2c2476c7f98cf067d;
#b=0xc1554db541bedcf0648f35904d527cb649e541bf948579eb4e92ebf5279ba778;

#print("aR=", Fp(a*res[1]));


#sec384 modulus
#n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
#name="p_384r1"

#res=precompute_montgomery_constants(name, n, 128, 384);

#Fp:=GF(StringToInteger("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",16));
#res:=(Fp!(2^256))^-1;
#IntegerToString(Integers()!res,16);


#Fp:=GF(StringToInteger("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF",16);








