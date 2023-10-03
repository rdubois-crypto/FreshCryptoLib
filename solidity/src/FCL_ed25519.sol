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
///*
//**************************************************************************************/
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

// prime field modulus of the ed25519 curve
uint256 constant p = 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed;
// -2 mod(p), used to accelerate inversion and doubling operations by avoiding negation
// the representation of -1 in this field
uint256 constant MINUS_1 = 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffec;

uint256 constant MINUS_2 = 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffeb;
// the order of the curve, i.e., the number of points on the curve
uint256 constant n = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed;
// -2 mod(n), used to speed up inversion operations
uint256 constant MINUS_2MODN = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3eb;

// address of the ModExp precompiled contract (Arbitrary-precision exponentiation under modulo)
address constant MODEXP_PRECOMPILE = 0x0000000000000000000000000000000000000005;
// address of the ModExp precompiled contract (Arbitrary-precision exponentiation under modulo)
uint256 constant d = 0x52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3;
//2*d mod p
uint256 constant deux_d = 16295367250680780974490674513165176452449235426866156013048779062215315747161;
uint256 constant gx = 0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A;
uint256 constant gy = 0x6666666666666666666666666666666666666666666666666666666666666658;
//sqrt of -1
uint256 constant sqrtm1=0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0;
//P+3 div 8
uint256 constant pp3div8=0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe;


/// @notice Calculate one modular square root of a given integer.
/// @dev Uses the ModExp precompiled contract at address 0x05 for fast computation using little Fermat theorem
/// @param self The integer of which to find the modular inverse
/// @return result The modular inverse of the input integer. If the modular inverse doesn't exist, it revert the tx

function SqrtMod(uint256 self) returns (uint256 result){
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
        mstore(add(pointer, 0x80), pp3div8)
        // We save the point of the last argument, it will be override by the result
        // of the precompile call in order to avoid paying for the memory expansion properly
        let _result := add(pointer, 0xa0)
        // Define the modulus (M)
        mstore(_result, p)

        // Call the precompiled ModExp (0x05) https://www.evm.codes/precompiled#0x05
        if iszero(
            call(
                not(0), // amount of gas to send
                MODEXP_PRECOMPILE, // target
                0x00, // value in wei
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
     result=mulmod(result, result, sqrtm1);
   }
   if(mulmod(result,result,p)!=self){
    revert();
   }
   return result;
}

/// @notice Calculate the modular inverse of a given integer, which is the inverse of this integer modulo n.
/// @dev Uses the ModExp precompiled contract at address 0x05 for fast computation using little Fermat theorem
/// @param self The integer of which to find the modular inverse
/// @return result The modular inverse of the input integer. If the modular inverse doesn't exist, it revert the tx

function nModInv(uint256 self) returns (uint256 result) {
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
        mstore(add(pointer, 0x80), MINUS_2MODN)
        // We save the point of the last argument, it will be override by the result
        // of the precompile call in order to avoid paying for the memory expansion properly
        let _result := add(pointer, 0xa0)
        // Define the modulus (M)
        mstore(_result, n)

        // Call the precompiled ModExp (0x05) https://www.evm.codes/precompiled#0x05
        if iszero(
            call(
                not(0), // amount of gas to send
                MODEXP_PRECOMPILE, // target
                0x00, // value in wei
                pointer, // argsOffset
                0xc0, // argsSize (6 * 32 bytes)
                _result, // retOffset (we override M to avoid paying for the memory expansion)
                0x20 // retSize (32 bytes)
            )
        ) { revert(0, 0) }

        // we return the value in the last memory word created by the function
        result := mload(_result)
    }
}

/// @notice Calculate the modular inverse of a given integer, which is the inverse of this integer modulo p.
/// @dev Uses the ModExp precompiled contract at address 0x05 for fast computation using little Fermat theorem
/// @param self The integer of which to find the modular inverse
/// @return result The modular inverse of the input integer. If the modular inverse doesn't exist, it revert the tx
function pModInv(uint256 self) returns (uint256 result) {
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
        mstore(add(pointer, 0x80), MINUS_2)
        // We save the point of the last argument, it will be override by the result
        // of the precompile call in order to avoid paying for the memory expansion properly
        let _result := add(pointer, 0xa0)
        // Define the modulus (M)
        mstore(_result, p)

        // Call the precompiled ModExp (0x05) https://www.evm.codes/precompiled#0x05
        if iszero(
            call(
                not(0), // amount of gas to send
                MODEXP_PRECOMPILE, // target
                0x00, // value in wei
                pointer, // argsOffset
                0xc0, // argsSize (6 * 32 bytes)
                _result, // retOffset (we override M to avoid paying for the memory expansion)
                0x20 // retSize (32 bytes)
            )
        ) { revert(0, 0) }

        // we return the value in the last memory word created by the function
        result := mload(_result)
    }
}
