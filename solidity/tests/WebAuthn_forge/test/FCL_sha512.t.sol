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
///* FILE: FCL_sha512.t.sol
///*
///*
///* DESCRIPTION: test file for sha512 hash function
///*
//**************************************************************************************/

pragma solidity ^0.8.20;

import "forge-std/Test.sol";
import "@solidity/FCL_sha512.sol";

contract Sha512Test is Test {
    constructor() {
        SHA_setUp();
    }

    sha512.SHA512_CTX ctx;

    address constant KK512 = address(uint160(0xcaca));

    uint256 constant SHA512_BLOCK_LENGTH8 = 128;
    uint256 constant SHA512_BLOCK_LENGTH64 = SHA512_BLOCK_LENGTH8 >> 3;
    uint256 constant SHA512_SHORT_BLOCK_LENGTH64 = SHA512_BLOCK_LENGTH64 - 4;
    uint256 constant SHA512_DIGEST_LENGTH = 64;

    function Sha_Init() internal pure returns (sha512.SHA512_CTX memory context) {
        context.state[0] = 0x6a09e667f3bcc908;
        context.state[1] = 0xbb67ae8584caa73b;
        context.state[2] = 0x3c6ef372fe94f82b;
        context.state[3] = 0xa54ff53a5f1d36f1;
        context.state[4] = 0x510e527fade682d1;
        context.state[5] = 0x9b05688c2b3e6c1f;
        context.state[6] = 0x1f83d9abfb41bd6b;
        context.state[7] = 0x5be0cd19137e2179;

        context.usedspace64 = 0;
        for (uint256 i = 0; i < SHA512_BLOCK_LENGTH64; i++) {
            context.buffer[i] = 0;
        }
    }

    function SHA_setUp() internal {
        uint256 _prec_address = 0xcaca;
        // uint256[4] memory T;
        //uint64 val;

        string memory deployData = vm.readFile("test/sha512_const.json");
        bytes memory K512 = abi.decode(vm.parseJson(deployData, ".K512"), (bytes));
        address a_K512 = address(uint160(_prec_address));
        vm.etch(a_K512, K512); //todo : replace with create

        /*assembly{
        extcodecopy(a_K512, T, 0, 8)
    }
    val=uint64(T[0]>>192);

    console.log("lu: ", string(K512));

    console.log("lu: %x %x", T[0], val);*/
    }

    function displayCtx(string memory comment, sha512.SHA512_CTX memory context) internal view {
        console.log(comment);
        console.log("State:%x %x", context.state[0], context.state[1]);
        console.log("%x %x", context.state[2], context.state[3]);
        console.log("%x %x", context.state[4], context.state[5]);
        console.log("%x %x", context.state[6], context.state[7]);

        console.log("Buffer:");
        for (uint256 i = 0; i < 16; i++) {
            console.log("%x", context.buffer[i]);
        }
    }

    function test_Transform() public {
        ctx = Sha_Init();
        //  displayCtx("***After Init",ctx);
        ctx.buffer[0] = 0x80636261; //"message abc";
        ctx.buffer[15] = 0x1800000000000000; //"padding"

        ctx = sha512.SHA512_Transform(ctx, ctx.buffer);
        // displayCtx("***After Transform", ctx);
        assertEq(ctx.state[0], 0xddaf35a193617aba);
    }

    //reference norm sha512(abc)
    function test_SHA512_abc() public {
        uint64[16] memory buffer = [
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0),
            uint64(0)
        ];

        buffer[0] = 0x6162638000000000; //"message abc";
        buffer[15] = 0x18; //"padding"

        uint256 lowh;
        uint256 highh;
        (lowh, highh) = sha512.SHA512(buffer);

        assertEq(lowh, 0xddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a);
    }

    //intermediate values obtained through simultations with python
    //input to h: 0x0x6291d657deec24024827e69c3abe01a30ce548a284743a445e3680d7db5ac3acfc51cd8e6218a1a38da47ed00230f0580816ed13ba3303ac5deb911548908025af82
    //66 bytes=0x210 bits
    //h=0xbf62c3fb850acebf2d240df6fe5f136359ab6728da6056e3c6ddabb4ae5748549ec08df799a1bc959b0558f8675832c0648b4a939956f62e8ff39319ffb4bf09
    //h mod q= 0x60ab51a60e3f1ceb60549479b152ae2f4a41d9dd8da0f6c3ef2892d51118e95

    //reduce a 512 bits number modulo N
    function red512Modq(uint256[2] memory val) public pure returns (uint256) {
        //2²⁵²+27742317777372353535851937790883648493
        //0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
        //0xffffffffffffffffffffffffffffffec6ef5bf4737dcf70d6ec31748d98951d is 2^256%N
        return addmod(
            mulmod(
                val[0],
                0xffffffffffffffffffffffffffffffec6ef5bf4737dcf70d6ec31748d98951d,
                0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
            ),
            val[1],
            0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
        );
    }

    function test_red512modq() public view {
        uint256[2] memory tstv = [
            0xbf62c3fb850acebf2d240df6fe5f136359ab6728da6056e3c6ddabb4ae574854,
            0x9ec08df799a1bc959b0558f8675832c0648b4a939956f62e8ff39319ffb4bf09
        ];

        uint256 res = red512Modq(tstv);

        console.log("reduction %x", res);
    }

    function test_SHA512_ed25519_3() public {
        uint256[2] memory res;
        //Rs, A, msg
        ctx.buffer[0] = 0x6291d657deec2402; //Rs
        ctx.buffer[1] = 0x4827e69c3abe01a3;
        ctx.buffer[2] = 0x0ce548a284743a44; //
        ctx.buffer[3] = 0x5e3680d7db5ac3ac;
        ctx.buffer[4] = 0xfc51cd8e6218a1a3; //pubY
        ctx.buffer[5] = 0x8da47ed00230f058;
        ctx.buffer[6] = 0x0816ed13ba3303ac;
        ctx.buffer[7] = 0x5deb911548908025;
        ctx.buffer[8] = 0xaf82800000000000;
        ctx.buffer[15] = 0x210; //padding, 66bytes=0x210 bits

        (res[0], res[1]) = sha512.SHA512(ctx.buffer);

        //console.log("hash=%x %x", res[0], res[1]);
        res = sha512.Swap512(res); //endianness curse

        assertEq(res[0], 0xbf62c3fb850acebf2d240df6fe5f136359ab6728da6056e3c6ddabb4ae574854);
        assertEq(res[1], 0x9ec08df799a1bc959b0558f8675832c0648b4a939956f62e8ff39319ffb4bf09);
    }
}
