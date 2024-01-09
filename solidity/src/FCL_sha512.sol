//********************************************************************************************/
//  ___           _       ___               _         _    _ _
// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__
// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
//                                |__/|_|
///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project
///* License: This software is licensed under MIT License
///* This Code may be reused including license and copyright notice.
///* See LICENSE file at the root folder of the project.
///* FILE: FCL_sha512.sol
///*
///*
///* DESCRIPTION: dummy SHA512 translation from https://opensource.apple.com/source/network_cmds/network_cmds-511/unbound/compat/sha512.c.auto.html
///* This work is still WIP, not working yet, for gas efficiency the function only process messages of size multiple of 64
//**************************************************************************************/
//Initialize hash values:
//(first 32 bits of the fractional parts of the square roots of the first 8 primes 2..19):
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

library sha512 {
    uint256 constant SHA512_BLOCK_LENGTH8 = 128;
    uint256 constant SHA512_BLOCK_LENGTH64 = SHA512_BLOCK_LENGTH8 >> 3;
    uint256 constant SHA512_SHORT_BLOCK_LENGTH64 = SHA512_BLOCK_LENGTH64 - 4;
    uint256 constant SHA512_DIGEST_LENGTH = 64;

    struct SHA512_CTX {
        uint64[8] state;
        uint256 usedspace64;
        uint256 bitcount;
        uint64[SHA512_BLOCK_LENGTH64] buffer;
    }

   
    function Swap512(uint256[2] memory w) internal pure returns(uint256[2] memory r)
    {
        r[0]=Swap256(w[1]);
        r[1]=Swap256(w[0]);

        return r;
    }

    function Swap256(uint256 w) internal pure returns (uint256 x)
    {
        return uint256(Swap128(uint128(w>>128)))^(uint256(Swap128(uint128(w&0xffffffffffffffffffffffffffffffff)))<<128);
    }

    function Swap128(uint128 w) internal pure returns (uint128 x)
    {
        return uint128(Swap64(uint64(w>>64)))^(uint128(Swap64(uint64(w&0xffffffffffffffff)))<<64);

    }

    function Swap64(uint64 w) internal pure returns (uint64 x){
     uint64 tmp= (w >> 32) | (w << 32);
	 tmp = ((tmp & 0xff00ff00ff00ff00) >> 8) |    ((tmp & 0x00ff00ff00ff00ff) << 8); 
	 x = ((tmp & 0xffff0000ffff0000) >> 16) |   ((tmp & 0x0000ffff0000ffff) << 16); 
    }

    function Sigma1_512(uint64 h) internal pure returns (uint64 x){

        return (((h) >> (14)) | ((h) << (64 - (14)))) ^ (((h) >> (18)) | ((h) << (64 - (18))))^ (((h) >> (41)) | ((h) << (64 - (41))));
    }

    function Sigma0_512(uint64 h) internal pure returns (uint64 x){
        return (((h) >> (28)) | ((h) << (64 - (28)))) ^ (((h) >> (34)) | ((h) << (64 - (34))))^ (((h) >> (39)) | ((h) << (64 - (39))));
    }

    function sigma1_512(uint64 h) internal pure returns (uint64 x){
          return (((h) >> (19)) | ((h) << (64 - (19)))) ^ (((h) >> (61)) | ((h) << (64 - (61))))^ (((h) >> (6)) );

    }

    function sigma0_512(uint64 h) internal pure returns (uint64 x){
          return (((h) >> (1)) | ((h) << (64 - (1)))) ^ (((h) >> (8)) | ((h) << (64 - (8))))^ (((h) >> (7)) );

    }

    function Ch(uint64 x,uint64 y,uint64 z)	internal pure returns(uint64 r){
        return (((x) & (y)) ^ ((~(x)) & (z)));
    }

    function Maj(uint64 x,uint64 y,uint64 z) internal pure  returns(uint64 r ){
	return (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)));
    }


    function Sha_Init() public pure returns (SHA512_CTX memory context){
         context.state[0]=   0x6a09e667f3bcc908;
         context.state[1]=   0xbb67ae8584caa73b;
         context.state[2]=   0x3c6ef372fe94f82b;
         context.state[3]=   0xa54ff53a5f1d36f1;
         context.state[4]=   0x510e527fade682d1;
         context.state[5]=   0x9b05688c2b3e6c1f;
         context.state[6]=   0x1f83d9abfb41bd6b;
         context.state[7]=   0x5be0cd19137e2179;
        
      
       context.usedspace64=0;
       for(uint i=0;i<SHA512_BLOCK_LENGTH64;i++) {
        context.buffer[i]=0;
       }
    }

    
    function SHA512_Transform(SHA512_CTX memory i_context, uint64[SHA512_BLOCK_LENGTH64] memory data) internal view returns(SHA512_CTX memory context)  {
        unchecked {
            
        context=i_context;

        uint64 a = context.state[0];
        uint64 b = context.state[1];
        uint64 c = context.state[2];
        uint64 d = context.state[3];
        uint64 e = context.state[4];
        uint64 f = context.state[5];
        uint64 g = context.state[6];
        uint64 h = context.state[7];
        uint64 j = 0;
    
        do {
            context.buffer[j] = Swap64(data[j]);   
           
            uint64 T1 = h + (((((e) >> (14)) | ((e) << (64 - (14)))) ^ (((e) >> (18)) | ((e) << (64 - (18))))^ (((e) >> (41)) | ((e) << (64 - (41)))))) + Ch(e, f, g) + uint64(k512(j)) + context.buffer[j];           
            uint64 T2 = Sigma0_512(a) + Maj(a, b, c);
            h = g;
            g = f;
            f = e;
            e = d + T1;
            d = c;
            c = b;
            b = a;
            a = T1 + T2;
            j++;
        } while (j < 16);

        do {
           
            /* Part of the message block expansion: */
            uint64 T1 = context.buffer[(j + 1) & 0x0f];
            T1 = sigma0_512(T1);
            uint64 T2 = context.buffer[(j+14)&0x0f];
		    T2 =  sigma1_512(T2);
          
            /* Apply the SHA-512 compression function to update a..h */
            T1 = h + Sigma1_512(e) + Ch(e, f, g) + uint64(k512(j)) +
		     (context.buffer[j&0x0f] += T2 + context.buffer[(j+9)&0x0f] + T1);
		      T2 = Sigma0_512(a) + Maj(a, b, c);

            h = g;
            g = f;
            f = e;
            e = d + T1;
            d = c;
            c = b;
            b = a;
            a = T1 + T2;

            j++;
        } while (j < 80);
        
        /* Compute the current intermediate hash value */
        context.state[0] += a;
        context.state[1] += b;
        context.state[2] += c;
        context.state[3] += d;
        context.state[4] += e;
        context.state[5] += f;
        context.state[6] += g;
        context.state[7] += h;
        }
        return context;
    }


    function k512(uint j) internal view returns (uint64 r)
    {
        uint256[1] memory T;
         assembly{
        extcodecopy(0xcaca, T, mul(j,8), 8)
    }
        r=uint64(T[0]>>192);

        return r;
    }
    //update hash with a 64-multiple bit block 

    function Sha_Update64(SHA512_CTX memory context, uint64[SHA512_BLOCK_LENGTH64] calldata datain, uint256 len64)
        public view
    {
        uint256 freespace64;

        if (len64 == 0) {
            return;
        }

        if (context.usedspace64 > 0) {
            /* Calculate how much free space is available in the buffer */
            freespace64 = SHA512_BLOCK_LENGTH64 - context.usedspace64;

            if (len64 >= freespace64) {
                /* Fill the buffer completely and process it */
            for(uint i=0; i<freespace64;i++){
             //   ADDINC128(context.bitcount, freespace << 3);
                len64 -= freespace64;
                context.buffer[context.usedspace64]=datain[i];
                SHA512_Transform(context, context.buffer);
            }
            }
            else {
                /* The buffer is not yet full */
                for(uint i=0; i<len64;i++){
                    context.buffer[context.usedspace64]=datain[i];
                    

                }
                //ADDINC128(context->bitcount, len << 3);
                /* Clean up: */
           //     usedspace = freespace = 0;
                return;
            }
        
        }
    }

    function Sha_Last(SHA512_CTX memory context) internal view{
      if (context.usedspace64 > 0) {
		/* Begin padding with a 1 bit: */
		context.buffer[context.usedspace64++] = 0x80;

		if (context.usedspace64 <= SHA512_SHORT_BLOCK_LENGTH64) {
			/* Set-up for the last transform: */
            for(uint i=0;i<SHA512_SHORT_BLOCK_LENGTH64 - context.usedspace64;i++){
			 context.buffer[context.usedspace64]=0;

            }
		} else {
			if (context.usedspace64 < SHA512_BLOCK_LENGTH64) {
				for(uint i=0;i<SHA512_BLOCK_LENGTH64 - context.usedspace64;i++){
			      context.buffer[context.usedspace64]=0;  
			    }
            }
			/* Do second-to-last transform: */
			//SHA512_Transform(context, context.buffer);

			/* And set-up for the last transform: */
			//MEMSET_BZERO(context->buffer, SHA512_BLOCK_LENGTH - 2);
		}

	}
    else {
		/* Prepare for final transform: */
		//MEMSET_BZERO(context->buffer, SHA512_SHORT_BLOCK_LENGTH);

		/* Begin padding with a 1 bit: */
		context.buffer[0] = 0x80;
	}
	/* Store the length of input data (in bits): */
	//cast_var.theChars = context->buffer;
	//cast_var.theLongs[SHA512_SHORT_BLOCK_LENGTH / 8] = context->bitcount[1];
	//cast_var.theLongs[SHA512_SHORT_BLOCK_LENGTH / 8 + 1] = context->bitcount[0];

	/* final transform: */
	SHA512_Transform(context, context.buffer);
    
    }


// a single step single block (padding included) of sha512 processing
function SHA512(uint64[16] memory data) internal view returns(uint256 low, uint256 high){
  uint64[16] memory buffer;
  uint256[1] memory T;

  uint64 a = 0x6a09e667f3bcc908;
  uint64 b = 0xbb67ae8584caa73b;
  uint64 c = 0x3c6ef372fe94f82b;
  uint64 d = 0xa54ff53a5f1d36f1;
  uint64 e = 0x510e527fade682d1;
  uint64 f = 0x9b05688c2b3e6c1f;
  uint64 g = 0x1f83d9abfb41bd6b;
  uint64 h = 0x5be0cd19137e2179;
  uint64 j = 0;
 unchecked{       
 do {
          
assembly{
            let T1:= mload(add(data,mul(32,j))) // buffer[j] =T1= (data[j]);   
            
            mstore(add(buffer, mul(32,j)), T1)    


            extcodecopy(0xcaca, T, mul(j,8), 8)
            T1:=    add(T1, shr(192, mload(T)))
            T1:=    add(add(h,T1),xor( and(e,f), and(not(e), g)))            


            T1:= and(0xffffffffffffffff,add(T1, xor(xor( or(shr(14,e), shl(50,e)) , or(shr(18,e), shl(46,e))), or(shr(41,e), shl(23,e)))      ))
           
            let T2:= xor(xor( or(shr(28,a), shl(36,a)) , or(shr(34,a), shl(30,a))), or(shr(39,a), shl(25,a)))   
            T2:= and(0xffffffffffffffff,add(T2, xor(xor(and(a,b), and(b,c)), and(a,c)))) //MAJ
            h := g
            g := f
            f := e
            e := and(0xffffffffffffffff,add(d , T1))
            d := c
            c := b
            b := a
            a := and(0xffffffffffffffff,add(T1 , T2))
            j :=add(j,1)
            }
                      

        } while (j < 16);

        do {
           
            /* Part of the message block expansion: */
            //uint64 T1 = buffer[(j + 1) & 0x0f];
            uint64 T1;uint64 T2;
               assembly{
                T1:= mload(add(buffer, mul(32,and(0x0f, add(j,1)))))   
                  T1:=  xor(xor( or(shr(1,T1), shl(63,T1)) , or(shr(8,T1), shl(56,T1))), shr(7,T1)     ) 
                T2:=mload(add(buffer, mul(32,and(0x0f, add(j,14)))))   
            
                T2:=  xor(xor( or(shr(19,T2), shl(45,T2)) , or(shr(61,T2), shl(3,T2))), shr(6,T2)     ) 
                T1:=add(T1,mload(add(buffer, mul(32,and(0x0f, add(j,9))))))   
                T2:=add(T2,T1)

                let addr:=add(buffer, mul(32,and(0x0f, j)))
                mstore(addr, and(0xffffffffffffffff,add(mload(addr), T2) ))
                  T1 :=mload(addr) 
          
            T1:=    add(add(h,T1),xor( and(e,f), and(not(e), g)))            
             T1:= and(0xffffffffffffffff,add(T1, xor(xor( or(shr(14,e), shl(50,e)) , or(shr(18,e), shl(46,e))), or(shr(41,e), shl(23,e)))      ))
           
            extcodecopy(0xcaca, T, mul(j,8), 8)
            T1:=    add(T1, shr(192, mload(T)))

            let T3:= xor(xor( or(shr(28,a), shl(36,a)) , or(shr(34,a), shl(30,a))), or(shr(39,a), shl(25,a)))   
            T3:= and(0xffffffffffffffff,add(T3, xor(xor(and(a,b), and(b,c)), and(a,c)))) //MAJ
            h := g
            g := f
            f := e
            e := and(0xffffffffffffffff,add(d , T1))
            d := c
            c := b
            b := a
            a := and(0xffffffffffffffff,add(T1 , T3))
            j :=add(j,1)
            }
        } while (j < 80);
 
         a+=0x6a09e667f3bcc908;
           b += 0xbb67ae8584caa73b;
   c += 0x3c6ef372fe94f82b;
   d += 0xa54ff53a5f1d36f1;
   e += 0x510e527fade682d1;
   f += 0x9b05688c2b3e6c1f;
   g += 0x1f83d9abfb41bd6b;
   h += 0x5be0cd19137e2179;
}
         
        low=(uint256(a)<<192)+(uint256(b)<<128)+(uint256(c)<<64)+d;
        high=(uint256(e)<<192)+(uint256(f)<<128)+(uint256(g)<<64)+h;

        return (low, high);
}

}
