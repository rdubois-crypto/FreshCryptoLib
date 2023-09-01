#//********************************************************************************************/
#//  ___           _       ___               _         _    _ _
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|
#///* Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project */
#///* License: This software is licensed under MIT License 	 */
#///* See LICENSE file at the root folder of the project.				 */
#///* FILE: FCL_test_ecrecover.sol						         */
#///* 											 */
#///* 											 */
#///* DESCRIPTION: generate a single vector to assess sage FCL_ethereum.sage implementation
#///* of ecrecover
#///* You can directly copy paste it in Remix:https://remix.ethereum.org/,
#///* click compile, deploy then test_recover will modify and display debug value
#//**************************************************************************************/

// SPDX-License-Identifier: MIT
pragma solidity ^0.8.17;

contract exampleRecover {

 //86607385674704723724606695130422278402341833367494746190489019919529834957660
 bytes32 r=0xe6aa80563d9931d611917eb184059d1eaa287d627406dab4460a319c59c071a2;
 bytes32 s=0x11aca1b3ecce3b4329139715e2b86dac82e596f6619f976ce186af9c4bf940d8;
 bytes32 h= 0xfeeaf9a319cce6371b255b4db39156d8fbc67fb4b4cd5debc30da0abc69e420b;
 uint8 v=28;

 address debug;

 function test_recover() public view returns (address b){
   address a=ecrecover(h,v,r,s);

   return a;
 }
}
