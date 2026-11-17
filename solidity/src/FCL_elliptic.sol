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
///* FILE: FCL_elliptic.sol
///*
///*
///* DESCRIPTION: modified XYZZ system coordinates for EVM elliptic point multiplication
///*  optimization
///*  uG+vQ: joint 2-bit window (15 points iG+jQ, batch normalized with one inversion),
///*  halved doubling (11 mulmod), mixed addition on a sign alternating accumulator (no sub),
///*  8-dim Shamir variants: scalars transposed in nibbles once, one table offset per round.
///*
//**************************************************************************************/
//* WARNING: this code SHALL not be used for non prime order curves for security reasons.
// Code is optimized for a=-3 only curves with prime order, constant like -1, -2 shall be replaced
// if ever used for other curve than sec256R1
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

library FCL_Elliptic_ZZ {
    // Set parameters for curve sec256r1.

    // address of the ModExp precompiled contract (Arbitrary-precision exponentiation under modulo)
    address constant MODEXP_PRECOMPILE = 0x0000000000000000000000000000000000000005;
    //curve prime field modulus
    uint256 constant p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;
    //short weierstrass first coefficient
    uint256 constant a = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC;
    //short weierstrass second coefficient
    uint256 constant b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B;
    //generating point affine coordinates
    uint256 constant gx = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296;
    uint256 constant gy = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5;
    //curve order (number of points)
    uint256 constant n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;
    /* -2 mod p constant, used to speed up inversion and doubling (avoid negation)*/
    uint256 constant minus_2 = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFD;
    /* -2 mod n constant, used to speed up inversion*/
    uint256 constant minus_2modn = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC63254F;

    uint256 constant minus_1 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF;
    //P+1 div 4
    uint256 constant pp1div4=0x3fffffffc0000000400000000000000000000000400000000000000000000000;
    //arbitrary constant to express no quadratic residuosity
    uint256 constant _NOTSQUARE=0xFFFFFFFF00000002000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;
    uint256 constant _NOTONCURVE=0xFFFFFFFF00000003000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;

    // 3/2 mod p, used by the halved doubling (U=Y instead of U=2Y)
    uint256 constant _C3H = 0x7fffffff80000000800000000000000000000000800000000000000000000001;
    // 2G and 3G, constant part of the 2-bit joint window table
    uint256 constant _G2X = 0x7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978;
    uint256 constant _G2Y = 0x07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1;
    uint256 constant _G3X = 0x5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c;
    uint256 constant _G3Y = 0x8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032;
    // table index (i+4j) of the 11 points computed in XYZZ: 2Q,3Q,Q+G,Q+2G,Q+3G,2Q+G,..,3Q+3G
    uint256 constant _DMAP = 0x080c050607090a0b0d0e0f000000000000000000000000000000000000000000;

    /**
     * /* inversion mod n via a^(n-2), use of precompiled using little Fermat theorem
     */
    function FCL_nModInv(uint256 u) internal view returns (uint256 result) {
        assembly {
            let pointer := mload(0x40)
            // Define length of base, exponent and modulus. 0x20 == 32 bytes
            mstore(pointer, 0x20)
            mstore(add(pointer, 0x20), 0x20)
            mstore(add(pointer, 0x40), 0x20)
            // Define variables base, exponent and modulus
            mstore(add(pointer, 0x60), u)
            mstore(add(pointer, 0x80), minus_2modn)
            mstore(add(pointer, 0xa0), n)

            // Call the precompiled contract 0x05 = ModExp
            if iszero(staticcall(not(0), 0x05, pointer, 0xc0, pointer, 0x20)) { revert(0, 0) }
            result := mload(pointer)
        }
    }
    /**
     * /* @dev inversion mod nusing little Fermat theorem via a^(n-2), use of precompiled
     */

    function FCL_pModInv(uint256 u) internal view returns (uint256 result) {
        assembly {
            let pointer := mload(0x40)
            // Define length of base, exponent and modulus. 0x20 == 32 bytes
            mstore(pointer, 0x20)
            mstore(add(pointer, 0x20), 0x20)
            mstore(add(pointer, 0x40), 0x20)
            // Define variables base, exponent and modulus
            mstore(add(pointer, 0x60), u)
            mstore(add(pointer, 0x80), minus_2)
            mstore(add(pointer, 0xa0), p)

            // Call the precompiled contract 0x05 = ModExp
            if iszero(staticcall(not(0), 0x05, pointer, 0xc0, pointer, 0x20)) { revert(0, 0) }
            result := mload(pointer)
        }
    }

    //Coron projective shuffling, take as input alpha as blinding factor
   function ecZZ_Coronize(uint256 alpha, uint256 x, uint256 y,  uint256 zz, uint256 zzz) internal pure  returns (uint256 x3, uint256 y3, uint256 zz3, uint256 zzz3)
   {
       
        uint256 alpha2=mulmod(alpha,alpha,p);
       
        x3=mulmod(alpha2, x,p); //alpha^-2.x
        y3=mulmod(mulmod(alpha, alpha2,p), y,p);

        zz3=mulmod(zz,alpha2,p);//alpha^2 zz
        zzz3=mulmod(zzz,mulmod(alpha, alpha2,p),p);//alpha^3 zzz
        
        return (x3, y3, zz3, zzz3);
   }


 function ecZZ_Add(uint256 x1, uint256 y1, uint256 zz1, uint256 zzz1, uint256 x2, uint256 y2, uint256 zz2, uint256 zzz2) internal pure  returns (uint256 x3, uint256 y3, uint256 zz3, uint256 zzz3)
  {
    uint256 u1=mulmod(x1,zz2,p); // U1 = X1*ZZ2
    uint256 u2=mulmod(x2, zz1,p);               //  U2 = X2*ZZ1
    u2=addmod(u2, p-u1, p);//  P = U2-U1
    x1=mulmod(u2, u2, p);//PP
    x2=mulmod(x1, u2, p);//PPP
    
    zz3=mulmod(x1, mulmod(zz1, zz2, p),p);//ZZ3 = ZZ1*ZZ2*PP  
    zzz3=mulmod(zzz1, mulmod(zzz2, x2, p),p);//ZZZ3 = ZZZ1*ZZZ2*PPP

    zz1=mulmod(y1, zzz2,p);  // S1 = Y1*ZZZ2
    zz2=mulmod(y2, zzz1, p);    // S2 = Y2*ZZZ1 
    zz2=addmod(zz2, p-zz1, p);//R = S2-S1
    zzz1=mulmod(u1, x1,p); //Q = U1*PP
    x3= addmod(addmod(mulmod(zz2, zz2, p), p-x2,p), mulmod(minus_2, zzz1,p),p); //X3 = R2-PPP-2*Q
    y3=addmod( mulmod(zz2, addmod(zzz1, p-x3, p),p), p-mulmod(zz1, x2, p),p);//R*(Q-X3)-S1*PPP

    return (x3, y3, zz3, zzz3);
  }

/// @notice Calculate one modular square root of a given integer. Assume that p=3 mod 4.
/// @dev Uses the ModExp precompiled contract at address 0x05 for fast computation using little Fermat theorem
/// @param self The integer of which to find the modular inverse
/// @return result The modular inverse of the input integer. If the modular inverse doesn't exist, it revert the tx

function SqrtMod(uint256 self) internal view returns (uint256 result){
 assembly ("memory-safe") {
        // load the free memory pointer value
        let pointer := mload(0x40)

        // Define length of base (Bsize)
        mstore(pointer, 0x20)
        // Define the exponent size (Esize)
        mstore(add(pointer, 0x20), 0x20)
        // Define the modulus size (Msize)
        mstore(add(pointer, 0x40), 0x20)
        // Define variables base (B)
        mstore(add(pointer, 0x60), self)
        // Define the exponent (E)
        mstore(add(pointer, 0x80), pp1div4)
        // We save the point of the last argument, it will be override by the result
        // of the precompile call in order to avoid paying for the memory expansion properly
        let _result := add(pointer, 0xa0)
        // Define the modulus (M)
        mstore(_result, p)

        // Call the precompiled ModExp (0x05) https://www.evm.codes/precompiled#0x05
        if iszero(
            staticcall(
                not(0), // amount of gas to send
                MODEXP_PRECOMPILE, // target
                pointer, // argsOffset
                0xc0, // argsSize (6 * 32 bytes)
                _result, // retOffset (we override M to avoid paying for the memory expansion)
                0x20 // retSize (32 bytes)
            )
        ) { revert(0, 0) }

  result := mload(_result)
//  result :=addmod(result,0,p)
 }
   if(mulmod(result,result,p)!=self){
     result=_NOTSQUARE;
   }
  
   return result;
}
    /**
     * /* @dev Convert from affine rep to XYZZ rep
     */
    function ecAff_SetZZ(uint256 x0, uint256 y0) internal pure returns (uint256[4] memory P) {
        unchecked {
            P[2] = 1; //ZZ
            P[3] = 1; //ZZZ
            P[0] = x0;
            P[1] = y0;
        }
    }

    function ec_Decompress(uint256 x, uint256 parity) internal view returns(uint256 y){ 

        uint256 y2=mulmod(x,mulmod(x,x,p),p);//x3
        y2=addmod(b,addmod(y2,mulmod(x,a,p),p),p);//x3+ax+b

        y=SqrtMod(y2);
        if(y==_NOTSQUARE){
           return _NOTONCURVE;
        }
        if((y&1)!=(parity&1)){
            y=p-y;
        }
    }

    /**
     * /* @dev Convert from XYZZ rep to affine rep
     */
    /*    https://hyperelliptic.org/EFD/g1p/auto-shortw-xyzz-3.html#addition-add-2008-s*/
    function ecZZ_SetAff(uint256 x, uint256 y, uint256 zz, uint256 zzz) internal view returns (uint256 x1, uint256 y1) {
        uint256 zzzInv = FCL_pModInv(zzz); //1/zzz
        y1 = mulmod(y, zzzInv, p); //Y/zzz
        uint256 _b = mulmod(zz, zzzInv, p); //1/z
        zzzInv = mulmod(_b, _b, p); //1/zz
        x1 = mulmod(x, zzzInv, p); //X/zz
    }

    /**
     * /* @dev Sutherland2008 doubling
     */
    /* The "dbl-2008-s-1" doubling formulas */

    function ecZZ_Dbl(uint256 x, uint256 y, uint256 zz, uint256 zzz)
        internal
        pure
        returns (uint256 P0, uint256 P1, uint256 P2, uint256 P3)
    {
        unchecked {
            assembly {
                P0 := mulmod(2, y, p) //U = 2*Y1
                P2 := mulmod(P0, P0, p) // V=U^2
                P3 := mulmod(x, P2, p) // S = X1*V
                P1 := mulmod(P0, P2, p) // W=UV
                P2 := mulmod(P2, zz, p) //zz3=V*ZZ1
                zz := mulmod(3, mulmod(addmod(x, sub(p, zz), p), addmod(x, zz, p), p), p) //M=3*(X1-ZZ1)*(X1+ZZ1)
                P0 := addmod(mulmod(zz, zz, p), mulmod(minus_2, P3, p), p) //X3=M^2-2S
                x := mulmod(zz, addmod(P3, sub(p, P0), p), p) //M(S-X3)
                P3 := mulmod(P1, zzz, p) //zzz3=W*zzz1
                P1 := addmod(x, sub(p, mulmod(P1, y, p)), p) //Y3= M(S-X3)-W*Y1
            }
        }
        return (P0, P1, P2, P3);
    }

    /**
     * @dev Sutherland2008 add a ZZ point with a normalized point and greedy formulae
     * warning: assume that P1(x1,y1)!=P2(x2,y2), true in multiplication loop with prime order (cofactor 1)
     */

    function ecZZ_AddN(uint256 x1, uint256 y1, uint256 zz1, uint256 zzz1, uint256 x2, uint256 y2)
        internal
        pure
        returns (uint256 P0, uint256 P1, uint256 P2, uint256 P3)
    {
        unchecked {
            if (y1 == 0) {
                return (x2, y2, 1, 1);
            }

            assembly {
                y1 := sub(p, y1)
                y2 := addmod(mulmod(y2, zzz1, p), y1, p)
                x2 := addmod(mulmod(x2, zz1, p), sub(p, x1), p)
                P0 := mulmod(x2, x2, p) //PP = P^2
                P1 := mulmod(P0, x2, p) //PPP = P*PP
                P2 := mulmod(zz1, P0, p) ////ZZ3 = ZZ1*PP
                P3 := mulmod(zzz1, P1, p) ////ZZZ3 = ZZZ1*PPP
                zz1 := mulmod(x1, P0, p) //Q = X1*PP
                P0 := addmod(addmod(mulmod(y2, y2, p), sub(p, P1), p), mulmod(minus_2, zz1, p), p) //R^2-PPP-2*Q
                P1 := addmod(mulmod(addmod(zz1, sub(p, P0), p), y2, p), mulmod(y1, P1, p), p) //R*(Q-X3)
            }
            //end assembly
        } //end unchecked
        return (P0, P1, P2, P3);
    }

    /**
     * @dev Return the zero curve in XYZZ coordinates.
     */
    function ecZZ_SetZero() internal pure returns (uint256 x, uint256 y, uint256 zz, uint256 zzz) {
        return (0, 0, 0, 0);
    }
    /**
     * @dev Check if point is the neutral of the curve
     */

    // uint256 x0, uint256 y0, uint256 zz0, uint256 zzz0
    function ecZZ_IsZero(uint256, uint256 y0, uint256, uint256) internal pure returns (bool) {
        return y0 == 0;
    }
    /**
     * @dev Return the zero curve in affine coordinates. Compatible with the double formulae (no special case)
     */

    function ecAff_SetZero() internal pure returns (uint256 x, uint256 y) {
        return (0, 0);
    }

    /**
     * @dev Check if the curve is the zero curve in affine rep.
     */
    // uint256 x, uint256 y)
    function ecAff_IsZero(uint256, uint256 y) internal pure returns (bool flag) {
        return (y == 0);
    }

    /**
     * @dev Check if a point in affine coordinates is on the curve (reject Neutral that is indeed on the curve).
     */
    function ecAff_isOnCurve(uint256 x, uint256 y) internal pure returns (bool) {
        if (x >= p || y >= p || ((x == 0) && (y == 0))) {
            return false;
        }
        unchecked {
            uint256 LHS = mulmod(y, y, p); // y^2
            uint256 RHS = addmod(mulmod(mulmod(x, x, p), x, p), mulmod(x, a, p), p); // x^3+ax
            RHS = addmod(RHS, b, p); // x^3 + a*x + b

            return LHS == RHS;
        }
    }

    /**
     * @dev Add two elliptic curve points in affine coordinates. Deal with P=Q
     */

    function ecAff_add(uint256 x0, uint256 y0, uint256 x1, uint256 y1) internal view returns (uint256, uint256) {
        uint256 zz0;
        uint256 zzz0;

        if (ecAff_IsZero(x0, y0)) return (x1, y1);
        if (ecAff_IsZero(x1, y1)) return (x0, y0);
        if((x0==x1)&&(y0==y1)) {
            (x0, y0, zz0, zzz0) = ecZZ_Dbl(x0, y0,1,1);
        }
        else{
            (x0, y0, zz0, zzz0) = ecZZ_AddN(x0, y0, 1, 1, x1, y1);
        }

        return ecZZ_SetAff(x0, y0, zz0, zzz0);
    }

    /**
     * @dev Detects keys Q = +/-c.G with c in {1, 2, 3, 1/2, 3/2, 1/3, 2/3} (mod n).
     *      For those keys some entries iG+jQ of the 2-bit joint window table are the point at infinity
     *      or require a doubling, so the multiplication is rewritten as (u + c.v).G.
     *      Returns the multiplier c (mod n), or 0 when Q is not such a key.
     */
    function ecZZ_smallMultipleOfG(uint256 Q0, uint256 Q1) internal pure returns (uint256 c) {
        assembly ("memory-safe") {
            let y
            switch Q0
            case 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296 {
                c := 1
                y := 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
            }
            case 0x7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978 {
                c := 2
                y := 0x07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1
            }
            case 0x5ecbe4d1a6330a44c8f7ef951d4bf165e6c6b721efada985fb41661bc6e7fd6c {
                c := 3
                y := 0x8734640c4998ff7e374b06ce1a64a2ecd82ab036384fb83d9a79b127a27d5032
            }
            case 0x2afa386b3f2bdcdb83f4d83f8fa3874d7b74dcb454bd644fdd6bf3d1f2da8db6 {
                c := 0x7fffffff800000007fffffffffffffffde737d56d38bcf4279dce5617e3192a9
                y := 0x8de7b41d3557a9cc9d4ac90ef7ad299a51759b030e1472b36b952a7686908d63
            }
            case 0xd897744db0ca3e5a506e1b02b731d784541db35865c2debc1dad0e1a757d18fc {
                c := 0x7fffffff800000007fffffffffffffffde737d56d38bcf4279dce5617e3192aa
                y := 0x31dbef292eeb2b81d165130d97b1cb644673439fad03103c722441787749e444
            }
            case 0x517d3f033d9b7d1994d200de245f8952bf5ac043d4014ca9af9ec20fee5119c8 {
                c := 0xaaaaaaaa00000000aaaaaaaaaaaaaaaa7def51c91a0fbf034d26872ca84218e1
                y := 0x7f60b6f5d7769fe2872c0c1f7f63193bcef9fc2a853f6518257d015e3baf805c
            }
            case 0x69b9b213061d41b3d9ac285f203e898e6257c079ba2dd20137ad032b155c6c12 {
                c := 0x555555550000000055555555555555553ef7a8e48d07df81a693439654210c71
                y := 0x0d6678ca01e385621ae907c5d7bd47b5f4be12a8ec26037836d5fd46e8d4d5de
            }
            if c { if iszero(eq(Q1, y)) { c := sub(n, c) } }
        }
    }

    /**
     * @dev Core of uG+vQ: joint 2-bit window (Shamir/Strauss), table of the 15 points iG+jQ,
     *      normalized to affine with a single modular inversion (Montgomery batch trick).
     *      Returns the result in XYZZ coordinates (ZZ==0 iff point at infinity, then X==0).
     *      Table entry k (k = i + 4j) is stored at T + 128k as [-x, y, -y, x].
     *      Q = (0,0) is read as the point at infinity, any 256-bit scalars are accepted.
     */
    function ecZZ_mulmuladd_core(uint256 Q0, uint256 Q1, uint256 scalar_u, uint256 scalar_v)
        internal
        view
        returns (uint256 X, uint256 Y, uint256 ZZ, uint256 ZZZ)
    {
        unchecked {
            if (Q0 == 0 && Q1 == 0) scalar_v = 0;
            if (scalar_v != 0) {
                uint256 c = ecZZ_smallMultipleOfG(Q0, Q1);
                if (c != 0) {
                    scalar_u = addmod(scalar_u, mulmod(c, scalar_v, n), n);
                    scalar_v = 0;
                }
            }
            if (scalar_u == 0 && scalar_v == 0) return (0, 0, 0, 0);
        }
        assembly ("memory-safe") {
            /* I. Precomputations */
            function st4(ptr, x, y) {
                mstore(ptr, sub(p, x))
                mstore(add(ptr, 0x20), y)
                mstore(add(ptr, 0x40), sub(p, y))
                mstore(add(ptr, 0x60), x)
            }
            function stZ(ptr, x, y, zz, zzz) {
                mstore(ptr, x)
                mstore(add(ptr, 0x20), y)
                mstore(add(ptr, 0x40), zz)
                mstore(add(ptr, 0x60), zzz)
            }
            //generic doubling, output 2P (not -2P)
            function dblStd(x, y, zz, zzz) -> X3, Y3, ZZ3, ZZZ3 {
                let U := addmod(y, y, p)
                let V := mulmod(U, U, p)
                let W := mulmod(U, V, p)
                let S := mulmod(x, V, p)
                let M := mulmod(3, mulmod(addmod(x, sub(p, zz), p), addmod(x, zz, p), p), p)
                X3 := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                Y3 := addmod(mulmod(M, addmod(S, sub(p, X3), p), p), sub(p, mulmod(W, y, p)), p)
                ZZ3 := mulmod(V, zz, p)
                ZZZ3 := mulmod(W, zzz, p)
            }
            //P1 (XYZZ) + P2 (affine), P1 != +/-P2, P1 != O guaranteed by the small multiple filter
            function addStd(x1, y1, zz1, zzz1, x2, y2) -> X3, Y3, ZZ3, ZZZ3 {
                let P := addmod(mulmod(x2, zz1, p), sub(p, x1), p)
                let R := addmod(mulmod(y2, zzz1, p), sub(p, y1), p)
                let PP := mulmod(P, P, p)
                let PPP := mulmod(P, PP, p)
                let Qq := mulmod(x1, PP, p)
                X3 := addmod(addmod(mulmod(R, R, p), sub(p, PPP), p), mulmod(minus_2, Qq, p), p)
                Y3 := addmod(mulmod(R, addmod(Qq, sub(p, X3), p), p), sub(p, mulmod(y1, PPP, p)), p)
                ZZ3 := mulmod(zz1, PP, p)
                ZZZ3 := mulmod(zzz1, PPP, p)
            }

            //halved doubling, output -2P (FCL negated-Y convention), rare path only
            function dblNeg(x, y, zz, zzz) -> X3, Y3, ZZ3, ZZZ3 {
                let M := mulmod(mulmod(addmod(x, sub(p, zz), p), addmod(x, zz, p), p), _C3H, p)
                let V := mulmod(y, y, p)
                ZZ3 := mulmod(V, zz, p)
                let S := mulmod(x, V, p)
                V := mulmod(y, V, p)
                ZZZ3 := mulmod(V, zzz, p)
                X3 := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                Y3 := addmod(mulmod(V, y, p), mulmod(M, addmod(X3, sub(p, S), p), p), p)
            }
            function addSt(src, dst, x2, y2) {
                let ta, tb, tc, td :=
                    addStd(mload(src), mload(add(src, 0x20)), mload(add(src, 0x40)), mload(add(src, 0x60)), x2, y2)
                stZ(dst, ta, tb, tc, td)
            }

            //table entries depending on Q: Q, 2Q, 3Q, jQ+iG, computed in XYZZ then batch normalized
            function buildQ(T, qx, qy) {
                st4(add(T, 0x200), qx, qy)
                // 11 points in XYZZ: 2Q, 3Q, Q+iG, 2Q+iG, 3Q+iG (i=1..3), table index in _DMAP
                let S0 := add(T, 0x800)
                {
                    let ta, tb, tc, td := dblStd(qx, qy, 1, 1)
                    stZ(S0, ta, tb, tc, td)
                    ta, tb, tc, td := addStd(qx, qy, 1, 1, gx, gy)
                    stZ(add(S0, 0x100), ta, tb, tc, td)
                    ta, tb, tc, td := addStd(qx, qy, 1, 1, _G2X, _G2Y)
                    stZ(add(S0, 0x180), ta, tb, tc, td)
                    ta, tb, tc, td := addStd(qx, qy, 1, 1, _G3X, _G3Y)
                    stZ(add(S0, 0x200), ta, tb, tc, td)
                }
                addSt(S0, add(S0, 0x80), qx, qy)
                addSt(S0, add(S0, 0x280), gx, gy)
                addSt(S0, add(S0, 0x300), _G2X, _G2Y)
                addSt(S0, add(S0, 0x380), _G3X, _G3Y)
                addSt(add(S0, 0x80), add(S0, 0x400), gx, gy)
                addSt(add(S0, 0x80), add(S0, 0x480), _G2X, _G2Y)
                addSt(add(S0, 0x80), add(S0, 0x500), _G3X, _G3Y)

                // Montgomery batch inversion of the 11 ZZZ
                let PF := add(S0, 0x580)
                let ta := 1
                for { let k := 0 } lt(k, 0x580) { k := add(k, 0x80) } {
                    ta := mulmod(ta, mload(add(S0, add(k, 0x60))), p)
                    mstore(add(PF, shr(2, k)), ta)
                }
                {
                    let M0 := add(PF, 0x160)
                    mstore(M0, 0x20)
                    mstore(add(M0, 0x20), 0x20)
                    mstore(add(M0, 0x40), 0x20)
                    mstore(add(M0, 0x60), ta)
                    mstore(add(M0, 0x80), minus_2)
                    mstore(add(M0, 0xa0), p)
                    if iszero(staticcall(gas(), 0x05, M0, 0xc0, M0, 0x20)) { revert(0, 0) }
                    ta := mload(M0) //1/prod(ZZZ)
                }
                for { let k := 0x500 } lt(k, 0x580) { k := sub(k, 0x80) } {
                    let pt := add(S0, k)
                    let tb := ta //1/ZZZ_k
                    if k { tb := mulmod(ta, mload(add(PF, sub(shr(2, k), 0x20))), p) }
                    ta := mulmod(ta, mload(add(pt, 0x60)), p)
                    let tc := mulmod(mload(add(pt, 0x40)), tb, p) //ZZ/ZZZ = 1/Z
                    st4(
                        add(T, shl(7, byte(shr(7, k), _DMAP))),
                        mulmod(mload(pt), mulmod(tc, tc, p), p),
                        mulmod(mload(add(pt, 0x20)), tb, p)
                    )
                }
            }

            let T := mload(0x40)
            st4(add(T, 0x80), gx, gy)
            st4(add(T, 0x100), _G2X, _G2Y)
            st4(add(T, 0x180), _G3X, _G3Y)

            if scalar_v { buildQ(T, Q0, Q1) }

            /* II. Main loop, 4 bits (2 of u, 2 of v) per iteration, accumulator sign alternates:
                   first half adds -T[d] to -(4A), second half adds T[d] to (4A) */
            Q1 := 252
            for {} iszero(or(shr(Q1, scalar_u), shr(Q1, scalar_v))) { Q1 := sub(Q1, 4) } {}
            Q0 := T
            X := 0
            Y := 0
            ZZ := 0
            ZZZ := 0
            for {} lt(Q1, 256) { Q1 := sub(Q1, 4) } {
                {
                    let M := mulmod(mulmod(addmod(X, sub(p, ZZ), p), addmod(X, ZZ, p), p), _C3H, p)
                    let V := mulmod(Y, Y, p)
                    ZZZ := mulmod(mulmod(Y, V, p), ZZZ, p)
                    ZZ := mulmod(V, ZZ, p)
                    let S := mulmod(X, V, p)
                    X := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                    Y := addmod(mulmod(V, V, p), mulmod(M, addmod(X, sub(p, S), p), p), p)
                }
                {
                    let M := mulmod(mulmod(addmod(X, sub(p, ZZ), p), addmod(X, ZZ, p), p), _C3H, p)
                    let V := mulmod(Y, Y, p)
                    ZZZ := mulmod(mulmod(Y, V, p), ZZZ, p)
                    ZZ := mulmod(V, ZZ, p)
                    let S := mulmod(X, V, p)
                    X := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                    Y := addmod(mulmod(V, V, p), mulmod(M, addmod(X, sub(p, S), p), p), p)
                }
                {
                    let e :=
                        or(and(shl(5, shr(Q1, scalar_u)), 0x180), and(shl(7, shr(Q1, scalar_v)), 0x600))
                    switch e
                    case 0 { Y := sub(p, Y) }
                    default {
                        e := add(Q0, e)
                        let Pn := addmod(mulmod(mload(e), ZZ, p), X, p)
                        let R := addmod(mulmod(mload(add(e, 0x40)), ZZZ, p), Y, p)
                        switch Pn
                        case 0 {
                            switch ZZ
                            case 0 {
                                X := mload(add(e, 0x60))
                                Y := mload(add(e, 0x40))
                                ZZ := 1
                                ZZZ := 1
                            }
                            default {
                                switch R
                                case 0 { X, Y, ZZ, ZZZ := dblNeg(X, Y, ZZ, ZZZ) }
                                default {
                                    X := 0
                                    ZZ := 0
                                    ZZZ := 0
                                }
                            }
                        }
                        default {
                            let PP := mulmod(Pn, Pn, p)
                            Pn := mulmod(Pn, PP, p)
                            ZZ := mulmod(ZZ, PP, p)
                            ZZZ := mulmod(ZZZ, Pn, p)
                            PP := mulmod(X, PP, p)
                            X := addmod(addmod(mulmod(R, R, p), Pn, p), mulmod(minus_2, PP, p), p)
                            Y := addmod(mulmod(R, addmod(X, sub(p, PP), p), p), mulmod(Y, Pn, p), p)
                        }
                    }
                }
                {
                    let M := mulmod(mulmod(addmod(X, sub(p, ZZ), p), addmod(X, ZZ, p), p), _C3H, p)
                    let V := mulmod(Y, Y, p)
                    ZZZ := mulmod(mulmod(Y, V, p), ZZZ, p)
                    ZZ := mulmod(V, ZZ, p)
                    let S := mulmod(X, V, p)
                    X := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                    Y := addmod(mulmod(V, V, p), mulmod(M, addmod(X, sub(p, S), p), p), p)
                }
                {
                    let M := mulmod(mulmod(addmod(X, sub(p, ZZ), p), addmod(X, ZZ, p), p), _C3H, p)
                    let V := mulmod(Y, Y, p)
                    ZZZ := mulmod(mulmod(Y, V, p), ZZZ, p)
                    ZZ := mulmod(V, ZZ, p)
                    let S := mulmod(X, V, p)
                    X := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                    Y := addmod(mulmod(V, V, p), mulmod(M, addmod(X, sub(p, S), p), p), p)
                }
                {
                    let e := or(and(shl(7, shr(Q1, scalar_u)), 0x180), and(shl(9, shr(Q1, scalar_v)), 0x600))
                    switch e
                    case 0 { Y := sub(p, Y) }
                    default {
                        e := add(Q0, e)
                        let Pn := addmod(mulmod(mload(e), ZZ, p), X, p)
                        let R := addmod(mulmod(mload(add(e, 0x20)), ZZZ, p), Y, p)
                        switch Pn
                        case 0 {
                            switch ZZ
                            case 0 {
                                X := mload(add(e, 0x60))
                                Y := mload(add(e, 0x20))
                                ZZ := 1
                                ZZZ := 1
                            }
                            default {
                                switch R
                                case 0 { X, Y, ZZ, ZZZ := dblNeg(X, Y, ZZ, ZZZ) }
                                default {
                                    X := 0
                                    ZZ := 0
                                    ZZZ := 0
                                }
                            }
                        }
                        default {
                            let PP := mulmod(Pn, Pn, p)
                            Pn := mulmod(Pn, PP, p)
                            ZZ := mulmod(ZZ, PP, p)
                            ZZZ := mulmod(ZZZ, Pn, p)
                            PP := mulmod(X, PP, p)
                            X := addmod(addmod(mulmod(R, R, p), Pn, p), mulmod(minus_2, PP, p), p)
                            Y := addmod(mulmod(R, addmod(X, sub(p, PP), p), p), mulmod(Y, Pn, p), p)
                        }
                    }
                }
            }
        }
    }

    /**
     * @dev Computation of uG+vQ using Strauss-Shamir's trick, G basepoint, Q public key
     *       Returns only x for ECDSA use
     *      */
    function ecZZ_mulmuladd_S_asm(
        uint256 Q0,
        uint256 Q1, //affine rep for input point Q
        uint256 scalar_u,
        uint256 scalar_v
    ) internal view returns (uint256 X) {
        uint256 ZZ;
        (X,, ZZ,) = ecZZ_mulmuladd_core(Q0, Q1, scalar_u, scalar_v);
        assembly ("memory-safe") {
            let T := mload(0x40)
            mstore(T, 0x20)
            mstore(add(T, 0x20), 0x20)
            mstore(add(T, 0x40), 0x20)
            mstore(add(T, 0x60), ZZ)
            mstore(add(T, 0x80), minus_2)
            mstore(add(T, 0xa0), p)
            if iszero(staticcall(gas(), 0x05, T, 0xc0, T, 0x20)) { revert(0, 0) }
            X := mulmod(X, mload(T), p) //X/zz
        }
    }

    /**
     * @dev Computation of uG+vQ using Strauss-Shamir's trick, G basepoint, Q public key
     *       Returns affine representation of point (normalized), (0,0) for the point at infinity
     *      */
    function ecZZ_mulmuladd(
        uint256 Q0,
        uint256 Q1, //affine rep for input point Q
        uint256 scalar_u,
        uint256 scalar_v
    ) internal view returns (uint256 X, uint256 Y) {
        uint256 ZZ;
        uint256 ZZZ;
        (X, Y, ZZ, ZZZ) = ecZZ_mulmuladd_core(Q0, Q1, scalar_u, scalar_v);
        assembly ("memory-safe") {
            let T := mload(0x40)
            mstore(T, 0x20)
            mstore(add(T, 0x20), 0x20)
            mstore(add(T, 0x40), 0x20)
            mstore(add(T, 0x60), ZZZ)
            mstore(add(T, 0x80), minus_2)
            mstore(add(T, 0xa0), p)
            if iszero(staticcall(gas(), 0x05, T, 0xc0, T, 0x20)) { revert(0, 0) }
            ZZZ := mload(T)
            Y := mulmod(Y, ZZZ, p) //Y/zzz
            ZZ := mulmod(ZZ, ZZZ, p) //1/z
            X := mulmod(X, mulmod(ZZ, ZZ, p), p) //X/zz
        }
    }

    /**
     * @dev Transposes the 4 limbs of 64 bits of a scalar into 64 nibbles:
     *      nibble j = (bit 192+j, bit 128+j, bit 64+j, bit j), i.e. the 4 bits used together by the 8-dim Shamir trick
     */
    function ecZZ_transpose4x64(uint256 s) internal pure returns (uint256 t) {
        assembly ("memory-safe") {
            function spread(x) -> r {
                r := and(or(x, shl(96, x)), 0xffffffff000000000000000000000000ffffffff)
                r := and(or(r, shl(48, r)), 0xffff000000000000ffff000000000000ffff000000000000ffff)
                r := and(or(r, shl(24, r)), 0xff000000ff000000ff000000ff000000ff000000ff000000ff000000ff)
                r := and(or(r, shl(12, r)), 0xf000f000f000f000f000f000f000f000f000f000f000f000f000f000f000f)
                r := and(or(r, shl(6, r)), 0x303030303030303030303030303030303030303030303030303030303030303)
                r := and(or(r, shl(3, r)), 0x1111111111111111111111111111111111111111111111111111111111111111)
            }
            let m := 0xffffffffffffffff
            t :=
                or(
                    or(spread(and(s, m)), shl(1, spread(and(shr(64, s), m)))),
                    or(shl(2, spread(and(shr(128, s), m))), shl(3, spread(shr(192, s))))
                )
        }
    }

    //8 dimensions Shamir's trick, using precomputations stored in Shamir8,  stored as Bytecode of an external
    //contract at given address dataPointer
    //(thx to Lakhdar https://github.com/Kelvyne for EVM storage explanations and tricks)
    // the external tool to generate tables from public key is in the /sage directory
    function ecZZ_mulmuladd_S8_extcode(uint256 scalar_u, uint256 scalar_v, address dataPointer)
        internal
        view
        returns (uint256 X /*, uint Y*/ )
    {
        scalar_u = ecZZ_transpose4x64(scalar_u);
        scalar_v = ecZZ_transpose4x64(scalar_v);
        assembly ("memory-safe") {
            //halved doubling, output -2P (FCL negated-Y convention), rare path only
            function dblNeg(x, y, zz, zzz) -> X3, Y3, ZZ3, ZZZ3 {
                let M := mulmod(mulmod(addmod(x, sub(p, zz), p), addmod(x, zz, p), p), _C3H, p)
                let V := mulmod(y, y, p)
                ZZZ3 := mulmod(mulmod(y, V, p), zzz, p)
                ZZ3 := mulmod(V, zz, p)
                let S := mulmod(x, V, p)
                X3 := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                Y3 := addmod(mulmod(V, V, p), mulmod(M, addmod(X3, sub(p, S), p), p), p)
            }
            if or(scalar_u, scalar_v) {
                let k := 252
                for {} iszero(or(shr(k, scalar_u), shr(k, scalar_v))) { k := sub(k, 4) } {}
                extcodecopy(
                    dataPointer, 0, or(and(shl(6, shr(k, scalar_u)), 0x3c0), and(shl(10, shr(k, scalar_v)), 0x3c00)), 64
                )
                X := mload(0)
                let Y := mload(0x20)
                let ZZ := 1
                let ZZZ := 1
                for { k := sub(k, 4) } lt(k, 256) { k := sub(k, 4) } {
                    {
                        let M := mulmod(mulmod(addmod(X, sub(p, ZZ), p), addmod(X, ZZ, p), p), _C3H, p)
                        let V := mulmod(Y, Y, p)
                        ZZZ := mulmod(mulmod(Y, V, p), ZZZ, p)
                        ZZ := mulmod(V, ZZ, p)
                        let S := mulmod(X, V, p)
                        X := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                        Y := addmod(mulmod(V, V, p), mulmod(M, addmod(X, sub(p, S), p), p), p)
                    }
                    let e := or(and(shl(6, shr(k, scalar_u)), 0x3c0), and(shl(10, shr(k, scalar_v)), 0x3c00))
                    switch e
                    case 0 { Y := sub(p, Y) }
                    default {
                        extcodecopy(dataPointer, 0, e, 64)
                        let Pn := addmod(X, sub(p, mulmod(mload(0), ZZ, p)), p)
                        let R := addmod(mulmod(mload(0x20), ZZZ, p), Y, p)
                        switch Pn
                        case 0 {
                            switch ZZ
                            case 0 {
                                X := mload(0)
                                Y := mload(0x20)
                                ZZ := 1
                                ZZZ := 1
                            }
                            default {
                                switch R
                                case 0 { X, Y, ZZ, ZZZ := dblNeg(X, Y, ZZ, ZZZ) }
                                default {
                                    X := 0
                                    ZZ := 0
                                    ZZZ := 0
                                }
                            }
                        }
                        default {
                            let PP := mulmod(Pn, Pn, p)
                            Pn := mulmod(Pn, PP, p)
                            ZZ := mulmod(ZZ, PP, p)
                            ZZZ := mulmod(ZZZ, Pn, p)
                            PP := mulmod(X, PP, p)
                            X := addmod(addmod(mulmod(R, R, p), Pn, p), mulmod(minus_2, PP, p), p)
                            Y := addmod(mulmod(R, addmod(X, sub(p, PP), p), p), mulmod(Y, Pn, p), p)
                        }
                    }
                }
                let T := mload(0x40)
                mstore(T, 0x20)
                mstore(add(T, 0x20), 0x20)
                mstore(add(T, 0x40), 0x20)
                mstore(add(T, 0x60), ZZ)
                mstore(add(T, 0x80), minus_2)
                mstore(add(T, 0xa0), p)
                if iszero(staticcall(gas(), 0x05, T, 0xc0, T, 0x20)) { revert(0, 0) }
                X := mulmod(X, mload(T), p) //X/zz
            }
        }
    }

    // improving the extcodecopy trick : append array at end of contract
    function ecZZ_mulmuladd_S8_hackmem(uint256 scalar_u, uint256 scalar_v, uint256 dataPointer)
        internal
        view
        returns (uint256 X /*, uint Y*/ )
    {
        scalar_u = ecZZ_transpose4x64(scalar_u);
        scalar_v = ecZZ_transpose4x64(scalar_v);
        assembly ("memory-safe") {
            //halved doubling, output -2P (FCL negated-Y convention), rare path only
            function dblNeg(x, y, zz, zzz) -> X3, Y3, ZZ3, ZZZ3 {
                let M := mulmod(mulmod(addmod(x, sub(p, zz), p), addmod(x, zz, p), p), _C3H, p)
                let V := mulmod(y, y, p)
                ZZZ3 := mulmod(mulmod(y, V, p), zzz, p)
                ZZ3 := mulmod(V, zz, p)
                let S := mulmod(x, V, p)
                X3 := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                Y3 := addmod(mulmod(V, V, p), mulmod(M, addmod(X3, sub(p, S), p), p), p)
            }
            if or(scalar_u, scalar_v) {
                let k := 252
                for {} iszero(or(shr(k, scalar_u), shr(k, scalar_v))) { k := sub(k, 4) } {}
                codecopy(
                    0,
                    add(dataPointer, or(and(shl(6, shr(k, scalar_u)), 0x3c0), and(shl(10, shr(k, scalar_v)), 0x3c00))),
                    64
                )
                X := mload(0)
                let Y := mload(0x20)
                let ZZ := 1
                let ZZZ := 1
                for { k := sub(k, 4) } lt(k, 256) { k := sub(k, 4) } {
                    {
                        let M := mulmod(mulmod(addmod(X, sub(p, ZZ), p), addmod(X, ZZ, p), p), _C3H, p)
                        let V := mulmod(Y, Y, p)
                        ZZZ := mulmod(mulmod(Y, V, p), ZZZ, p)
                        ZZ := mulmod(V, ZZ, p)
                        let S := mulmod(X, V, p)
                        X := addmod(mulmod(M, M, p), mulmod(minus_2, S, p), p)
                        Y := addmod(mulmod(V, V, p), mulmod(M, addmod(X, sub(p, S), p), p), p)
                    }
                    let e := or(and(shl(6, shr(k, scalar_u)), 0x3c0), and(shl(10, shr(k, scalar_v)), 0x3c00))
                    switch e
                    case 0 { Y := sub(p, Y) }
                    default {
                        codecopy(0, add(dataPointer, e), 64)
                        let Pn := addmod(X, sub(p, mulmod(mload(0), ZZ, p)), p)
                        let R := addmod(mulmod(mload(0x20), ZZZ, p), Y, p)
                        switch Pn
                        case 0 {
                            switch ZZ
                            case 0 {
                                X := mload(0)
                                Y := mload(0x20)
                                ZZ := 1
                                ZZZ := 1
                            }
                            default {
                                switch R
                                case 0 { X, Y, ZZ, ZZZ := dblNeg(X, Y, ZZ, ZZZ) }
                                default {
                                    X := 0
                                    ZZ := 0
                                    ZZZ := 0
                                }
                            }
                        }
                        default {
                            let PP := mulmod(Pn, Pn, p)
                            Pn := mulmod(Pn, PP, p)
                            ZZ := mulmod(ZZ, PP, p)
                            ZZZ := mulmod(ZZZ, Pn, p)
                            PP := mulmod(X, PP, p)
                            X := addmod(addmod(mulmod(R, R, p), Pn, p), mulmod(minus_2, PP, p), p)
                            Y := addmod(mulmod(R, addmod(X, sub(p, PP), p), p), mulmod(Y, Pn, p), p)
                        }
                    }
                }
                let T := mload(0x40)
                mstore(T, 0x20)
                mstore(add(T, 0x20), 0x20)
                mstore(add(T, 0x40), 0x20)
                mstore(add(T, 0x60), ZZ)
                mstore(add(T, 0x80), minus_2)
                mstore(add(T, 0xa0), p)
                if iszero(staticcall(gas(), 0x05, T, 0xc0, T, 0x20)) { revert(0, 0) }
                X := mulmod(X, mload(T), p) //X/zz
            }
        }
    }

    /**
     * @dev ECDSA verification using a precomputed table of multiples of P and Q stored in contract at address Shamir8
     *     generation of contract bytecode for precomputations is done using sagemath code
     *     (see sage directory, WebAuthn_precompute.sage)
     */

    /**
     * @dev ECDSA verification using a precomputed table of multiples of P and Q appended at end of contract at address endcontract
     *     generation of contract bytecode for precomputations is done using sagemath code
     *     (see sage directory, WebAuthn_precompute.sage)
     */

    function ecdsa_precomputed_hackmem(bytes32 message, uint256[2] calldata rs, uint256 endcontract)
        internal view
        returns (bool)
    {
        uint256 r = rs[0];
        uint256 s = rs[1];
        if (r == 0 || r >= n || s == 0 || s >= n) {
            return false;
        }
        /* Q is pushed via bytecode assumed to be correct
        if (!isOnCurve(Q[0], Q[1])) {
            return false;
        }*/

        uint256 sInv = FCL_nModInv(s);
        uint256 X;

        //Shamir 8 dimensions
        X = ecZZ_mulmuladd_S8_hackmem(mulmod(uint256(message), sInv, n), mulmod(r, sInv, n), endcontract);

        assembly {
            X := addmod(X, sub(n, r), n)
        }
        return X == 0;
    } //end  ecdsa_precomputed_verify()



} //EOF
