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
///* FILE: FCL_Edwards.t.sol
///*
///*
///* DESCRIPTION: test file for Edwards curves
///*
//**************************************************************************************/
// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.20;

import "forge-std/Test.sol";

import "@solidity/FCL_eddsa.sol";

contract EDDSATest is Test {
    //RFC8032, test vector 3
    function test_RFC8032_3() public pure {
        uint64[16] memory buffer;
        //public key:43933056957747458452560886832567536073542840507013052263144963060608791330050,
        //16962727616734173323702303146057009569815335830970791807500022961899349823996
        /*uint256[2] memory kpub = [
            43933056957747458452560886832567536073542840507013052263144963060608791330050,
            16962727616734173323702303146057009569815335830970791807500022961899349823996
        ];*/

        //sig value:0x6291d657deec24024827e69c3abe01a30ce548a284743a445e3680d7db5ac3ac18ff9b538d16f290ae67f760984dc6594a7c15e9716ed28dc027beceea1ec40a
        //msg:af82
        buffer[0] = 0x6291d657deec2402; //Rs
        buffer[1] = 0x4827e69c3abe01a3;
        buffer[2] = 0x0ce548a284743a44;
        buffer[3] = 0x5e3680d7db5ac3ac;
        buffer[4] = 0xfc51cd8e6218a1a3; //public y value, swapped
        buffer[5] = 0x8da47ed00230f058;
        buffer[6] = 0x0816ed13ba3303ac;
        buffer[7] = 0x5deb911548908025;
        buffer[8] = 0xaf82800000000000; //msg+padd
        buffer[15] = 0x210; //end of padding, 66bytes=0x210 bits
    }
}
