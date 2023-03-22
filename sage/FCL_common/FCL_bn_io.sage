#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _    
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__ 
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|                        
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project */
#///* License: This software is licensed under MIT License 	 */
#///* See LICENSE file at the root folder of the project.				 */
#///* FILE: FCL_bn_io.sage						         */
#///* 											 */
#///* 											 */
#///* DESCRIPTION: precompute a 8 dimensional table for Shamir's trick from a public key
#///* 
#//**************************************************************************************/

from sage.crypto.util import bin_to_ascii

#Convert a number to a chain of words of bitsize size_word, ordered from MSB to LSB
def FCL_int_to_words_MSB(A, size_word):
	sizeA=ceil(log(A)/log(2))  
	sizeM=ceil(sizeA/size_word)
	M=copy(matrix(1, sizeM)[0])
	
	mask = 2^size_word-1
	for i in [0..sizeM-2]:
		M[sizeM-1-i]= A 	& mask;
		A = A >> size_word
	M[0]= A 	& mask;
	
	return M;	
	

#Convert an int to the input of hashlib function
#inVal: integer reprensenting the number to convert
#size: size in bytes of the output 
def FCL_BN_to_bytes(inVal, size):
 bin_a=(bin(inVal)[2:]);
 bin_a=bin_a.zfill(size*8);
 
 return bin_to_ascii(bin_a).encode('iso-8859-1');
 	
