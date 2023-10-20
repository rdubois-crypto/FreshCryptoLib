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
///* FILE: FCL_eddsa.sol
///*
///*
///* DESCRIPTION: Implementation of RFC8032 using FCL Shamir's trick
///*
//**************************************************************************************/
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;

//import "./FCL_ed25519.sol" as Curve;
import "@solidity/FCL_sha512.sol";
import "@solidity/FCL_edwards.sol";


import {p, gx, gy, n, d, deux_d, MINUS_2, MINUS_1, MODEXP_PRECOMPILE} from "./FCL_ed25519.sol";

library EDDSA {



function SHA512modq(uint64[16] memory Data) internal view returns (uint256 h)
{

 uint256[2] memory val;
 (val[0], val[1])=sha512.SHA512(Data);

  return addmod(mulmod(val[0],0xffffffffffffffffffffffffffffffec6ef5bf4737dcf70d6ec31748d98951d, 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed)
                ,val[1],0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed);


}

//to reduce verification cost, the front shall present the concatenation
// R || public ||Msg ||Padding || S  
function Verify( uint256[2] calldata Q, uint64[16] memory RsPubMsgPad, uint256 s) internal returns(bool res)
  {

    //todo check decompression of R
    uint256 h=SHA512modq(RsPubMsgPad);
    
    uint256 x=Edwards.ed_mulmuladd(Q[0], Q[1], s,p-h);
    return(x==0);
  }




}
