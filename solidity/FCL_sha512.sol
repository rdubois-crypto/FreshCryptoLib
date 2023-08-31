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
///* This work is still WIP, not working yet
//**************************************************************************************/
//Initialize hash values:
//(first 32 bits of the fractional parts of the square roots of the first 8 primes 2..19):
contract sha512 {
    int256 constant SHA512_BLOCK_LENGTH8 = 128;
    int256 constant SHA512_BLOCK_LENGTH64 = SHA512_BLOCK_LENGTH8 >> 3;

    uint256 constant SHA512_DIGEST_LENGTH = 64;

    struct SHA512_CTX {
        uint64[8] state;
        uint256[2] bitcount;
        bytes buffer;
    }

    function Sha_Init(SHA512_CTX memory context) public {
        uint64[8] memory InitialK = [
            0x6a09e667f3bcc908,
            0xbb67ae8584caa73b,
            0x3c6ef372fe94f82b,
            0xa54ff53a5f1d36f1,
            0x510e527fade682d1,
            0x9b05688c2b3e6c1f,
            0x1f83d9abfb41bd6b,
            0x5be0cd19137e2179
        ];
    }

    function SHA512_Transform(SHA512_CTX memory context, uint64[SHA512_BLOCK_LENGTH64] memory data) internal {
        uint64[80] memory K512 = [
            0x428a2f98d728ae22,
            0x7137449123ef65cd,
            0xb5c0fbcfec4d3b2f,
            0xe9b5dba58189dbbc,
            0x3956c25bf348b538,
            0x59f111f1b605d019,
            0x923f82a4af194f9b,
            0xab1c5ed5da6d8118,
            0xd807aa98a3030242,
            0x12835b0145706fbe,
            0x243185be4ee4b28c,
            0x550c7dc3d5ffb4e2,
            0x72be5d74f27b896f,
            0x80deb1fe3b1696b1,
            0x9bdc06a725c71235,
            0xc19bf174cf692694,
            0xe49b69c19ef14ad2,
            0xefbe4786384f25e3,
            0x0fc19dc68b8cd5b5,
            0x240ca1cc77ac9c65,
            0x2de92c6f592b0275,
            0x4a7484aa6ea6e483,
            0x5cb0a9dcbd41fbd4,
            0x76f988da831153b5,
            0x983e5152ee66dfab,
            0xa831c66d2db43210,
            0xb00327c898fb213f,
            0xbf597fc7beef0ee4,
            0xc6e00bf33da88fc2,
            0xd5a79147930aa725,
            0x06ca6351e003826f,
            0x142929670a0e6e70,
            0x27b70a8546d22ffc,
            0x2e1b21385c26c926,
            0x4d2c6dfc5ac42aed,
            0x53380d139d95b3df,
            0x650a73548baf63de,
            0x766a0abb3c77b2a8,
            0x81c2c92e47edaee6,
            0x92722c851482353b,
            0xa2bfe8a14cf10364,
            0xa81a664bbc423001,
            0xc24b8b70d0f89791,
            0xc76c51a30654be30,
            0xd192e819d6ef5218,
            0xd69906245565a910,
            0xf40e35855771202a,
            0x106aa07032bbd1b8,
            0x19a4c116b8d2d0c8,
            0x1e376c085141ab53,
            0x2748774cdf8eeb99,
            0x34b0bcb5e19b48a8,
            0x391c0cb3c5c95a63,
            0x4ed8aa4ae3418acb,
            0x5b9cca4f7763e373,
            0x682e6ff3d6b2b8a3,
            0x748f82ee5defb2fc,
            0x78a5636f43172f60,
            0x84c87814a1f0ab72,
            0x8cc702081a6439ec,
            0x90befffa23631e28,
            0xa4506cebde82bde9,
            0xbef9a3f7b2c67915,
            0xc67178f2e372532b,
            0xca273eceea26619c,
            0xd186b8c721c0c207,
            0xeada7dd6cde0eb1e,
            0xf57d4f7fee6ed178,
            0x06f067aa72176fba,
            0x0a637dc5a2c898a6,
            0x113f9804bef90dae,
            0x1b710b35131c471b,
            0x28db77f523047d84,
            0x32caab7b40c72493,
            0x3c9ebe0a15c9bebc,
            0x431d67c49c100d4c,
            0x4cc5d4becb3e42b6,
            0x597f299cfc657e2a,
            0x5fcb6fab3ad6faec,
            0x6c44198c4a475817
        ];

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
            //T1 = h + Sigma1_512(e) + Ch(e, f, g) + K512[j] + (W512[j] = *data++);
            uint64 T1 = h + (((h) >> (14)) | ((h) << (64 - (14)))) ^ (((h) >> (18)) | ((h) << (64 - (18))))
                ^ (((h) >> (41)) | ((h) << (64 - (41))));
            context.buffer[j] = data[j];

            T1 += (((e) & (f)) ^ ((~(e)) & (g))) + K512[j] + context.buffer[j];
            //T2 = Sigma0_512(a) + Maj(a, b, c);
            uint64 T2 = (((a) >> (14)) | ((a) << (64 - (1)))) ^ (((a) >> (18)) | ((a) << (64 - (8))))
                ^ ((a) >> (7)) + (((a) & (b)) ^ ((a) & (c)) ^ ((b) & (c)));
            h = g;
            g = f;
            f = e;
            e = d + T1;
            d = c;
            c = b;
            b = a;
            a = T1 + T2;
        } while (j < 16);

        do {
            /* Part of the message block expansion: */
            s0 = W512[(j + 1) & 0x0f];
            //s0 = sigma0_512(s0);
            s0 = (((s0) >> (14)) | ((s0) << (64 - (1)))) ^ (((s0) >> (18)) | ((s0) << (64 - (8)))) ^ ((s0) >> (7));
            s1 = W512[(j + 14) & 0x0f];
            s1 = (((s1) >> (14)) | ((s1) << (64 - (14)))) ^ (((s1) >> (18)) | ((s1) << (64 - (18))))
                ^ (((s1) >> (41)) | ((s1) << (64 - (41)))); //s1 =  sigma1_512(s1);

            /* Apply the SHA-512 compression function to update a..h */
            //T1 = h + Sigma1_512(e) + Ch(e, f, g) + K512[j] +		     (W512[j&0x0f] += s1 + W512[(j+9)&0x0f] + s0);
            T1 = h + (((e) >> (14)) | ((e) << (64 - (14)))) ^ (((e) >> (18)) | ((e) << (64 - (18))))
                ^ (((e) >> (41)) | ((e) << (64 - (41))));
            context.buffer[j & 0x0f] += s1 + context.buffer[(j + 9) & 0x0f] + s0;
            T1 += (((e) & (f)) ^ ((~(e)) & (g))) + K512[j] + context.buffer[j & 0x0f];
            //T2 = Sigma0_512(a) + Maj(a, b, c);
            T2 = T2 = (((a) >> (14)) | ((a) << (64 - (1)))) ^ (((a) >> (18)) | ((a) << (64 - (8))))
                ^ ((a) >> (7)) + (((a) & (b)) ^ ((a) & (c)) ^ ((b) & (c)));
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

    //update hash with a 64 bit block multiple
    function Sha_Update64(SHA512_CTX memory context, uint64[SHA512_BLOCK_LENGTH64] calldata datain, uint256 len)
        public
    {
        uint256 usedspace64;
        uint256 freespace64;

        if (len == 0) {
            return;
        }
        usedspace64 = (context.bitcount[0] >> 6) % SHA512_BLOCK_LENGTH;
        if (usedspace > 0) {
            /* Calculate how much free space is available in the buffer */
            freespace = SHA512_BLOCK_LENGTH - usedspace;

            if (len >= freespace) {
                /* Fill the buffer completely and process it */
                MEMCPY_BCOPY(context.buffer[usedspace], data, freespace);
                ADDINC128(context.bitcount, freespace << 3);
                len -= freespace;
                data += freespace;
                //SHA512_Transform(context, (sha2_word64*)context->buffer);
            } else {
                /* The buffer is not yet full */
                //MEMCPY_BCOPY(&context->buffer[usedspace], data, len);
                //ADDINC128(context->bitcount, len << 3);
                /* Clean up: */
                usedspace = freespace = 0;
                return;
            }
        }
    }

    function Sha_Last(uint256 a) internal {}
}
