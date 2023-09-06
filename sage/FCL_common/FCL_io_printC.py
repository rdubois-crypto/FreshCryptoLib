#/**********************************************************************************/
#/* Copyright (C) 2022 - Renaud Dubois
#/* This file is part of cy_lib project
#/* License: This software is licensed under a dual BSD and GPL v2 license.
#/* DESCRIPTION: Conversions from integer to C representation */
#/**********************************************************************************/

#/*  Object : 	CYLIB_SAGE 		*/
#/*  Langage : 	sage	 		*/
#/*  Description :The CYLIB_SAGE is a simulation of the CYLIB library, designed for simulation and testing of the CYLIB.
#/*  FILE : cy_multmontgomery.sage is a simulation of a generic size Montgomery multiplier
#/*  Date:
#/*  Version:
#/*  Last modification:

WORD64_BITS = 64
WORD64_MASK = 0xffffffffffffffff
 
WORD128_BITS = 128
WORD128_MASK = 0xffffffffffffffffffffffffffffffff
 
import math
#/*********** Global variables ********************/ 

#/* target architecture machine words size */
g_size = -1

#/*********** Simulating machine words ************/ 
def Conv_word_LSB(A, size_word):
	sizeA=ceil(log(A)/log(2))  
	print("sizeA=",sizeA)
	sizeM=ceil(sizeA/size_word)
	print("sizeM=",sizeM)
	#M=copy(matrix(1, sizeM)[0])
	M= [0 for i in range(sizeM)] ;
	mask = 2**size_word-1
	for i in range(0,sizeM):
		M[i]= A 	& mask;
		print("i=",hex(M[i]));
		A = A >> size_word
	return M;	

#display a variable in a C like style
def Conv_word_MSB(A, size_word):
	sizeA=math.ceil(math.log(A)/math.log(2))  
	sizeM=math.ceil(sizeA/size_word)
	#M=copy(matrix(1, sizeM)[0])
	M= [0 for i in range(sizeM)] ;
	
	mask = 2**size_word-1
	#print("{", end='');
	for i in range(0,sizeM-1):
		M[sizeM-1-i]= A 	& mask;
#		/*print(" ",hex(M[i]),",", end='');*/
		A = A >> size_word
	M[0]= A 	& mask;
#	/*print(" ",hex(M[i]),"};");*/
	
	return M;	
		
#display a variable in a C like style with prefix and suffix comments
def Print_C_MSB(prefix, name, A, size_word, suffix):
	print(prefix, end='');
	sizeA=math.ceil(math.log(A)/math.log(2))  
	sizeM=math.ceil(sizeA/size_word)
	
	print("\n uint8_t ", name, "[", sizeM,"]=", end='');
	Conv_word_MSB(A, size_word);
	print(suffix);

#write a variable in a C like style with prefix and suffix comments into a file
def fprint_c_MSB(f, prefix, name, A, size_word, suffix):
	
	concat_str=prefix;
	sizeA=math.ceil(math.log(A)/math.log(2))  
	sizeM=math.ceil(sizeA/size_word)
	
	concat_str=concat_str+"\n uint8_t "+ name+ "["+str( sizeM)+"]={";
	M=Conv_word_MSB(A, size_word);
	for i in range(0,sizeM-1):
		concat_str=concat_str+" "+hex(M[i])+",";
	concat_str=concat_str+" "+hex(M[sizeM-1])+"};";
	concat_str=concat_str+suffix;
	
	f.write(concat_str);
	
	return concat_str;
		
def Conv_Num(M, size_word):
	A=0;
	for i in range(0..len(u)):
		A=A+M[i]<<(size_word*i);
	return A;

def ZERO_MWORD(a) :
	a=0
	return a;
	
def MUL_MWORD(tn, b, a):
	tn=a*b
	return tn;

def ADDTO_DWORD(tn, tu64_carry, l_carry):
	tn=tn+l_carry;
	return tn;

def  ADDTO_MWORD(r, tu64_carry, l_carry):
	r+=l_carry
	
	return r;
	
def ADDCTO_MWORD(r, tu64_carry, l_carry):
	r+=l_carry
	tu64_carry= r>>g_size;
	return r;

def WORD_MUL_LOW(u, r, n0, mask):
	u=(r*n0)&mask;
	return u;

def COPY_MWORD(l_carry,tn):
	l_carry=tn;
	return l_carry;

def PRINT_HEXATAB(r, size):
	for i in range(0..size):
		print(hex(r[i]))
	return 0;	
	
	
	
	
	
	
