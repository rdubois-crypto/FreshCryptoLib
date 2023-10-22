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
///* FILE: FCL_ecdsa.sol
///*
///*
///* DESCRIPTION: ecdsa verification implementation
///*
//**************************************************************************************/
//* WARNING: this code SHALL not be used for non prime order curves for security reasons.
// Code is optimized for a=-3 only curves with prime order, constant like -1, -2 shall be replaced
// if ever used for other curve than sec256R1
// SPDX-License-Identifier: MIT
pragma solidity >=0.8.19 <0.9.0;


import {FCL_Elliptic_ZZ} from "./FCL_elliptic.sol";



library FCL_ecdsa {
    // Set parameters for curve sec256r1.public
      //curve order (number of points)
    uint256 constant n = FCL_Elliptic_ZZ.n;
  
    /**
     * @dev ECDSA verification, given , signature, and public key.
     */

    /**
     * @dev ECDSA verification, given , signature, and public key, no calldata version
     */
    function ecdsa_verify(bytes32 message, uint256 r, uint256 s, uint256 Qx, uint256 Qy)  internal view returns (bool){

        if (r == 0 || r >= FCL_Elliptic_ZZ.n || s == 0 || s >= FCL_Elliptic_ZZ.n) {
            return false;
        }
        
        if (!FCL_Elliptic_ZZ.ecAff_isOnCurve(Qx, Qy)) {
            return false;
        }

        uint256 sInv = FCL_Elliptic_ZZ.FCL_nModInv(s);

        uint256 scalar_u = mulmod(uint256(message), sInv, FCL_Elliptic_ZZ.n);
        uint256 scalar_v = mulmod(r, sInv, FCL_Elliptic_ZZ.n);
        uint256 x1;

        x1 = FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(Qx, Qy, scalar_u, scalar_v);

        x1= addmod(x1, n-r,n );
        
      
       
        return x1 == 0;

    }

    function ec_recover_r1(uint256 h, uint256 v, uint256 r, uint256 s) internal view returns (address)
    {
         if (r == 0 || r >= FCL_Elliptic_ZZ.n || s == 0 || s >= FCL_Elliptic_ZZ.n) {
            return address(0);
        }
        uint256 y=FCL_Elliptic_ZZ.ec_Decompress(r, v-27);
        uint256 rinv=FCL_Elliptic_ZZ.FCL_nModInv(r);
        uint256 u1=mulmod(FCL_Elliptic_ZZ.n-addmod(0,h,FCL_Elliptic_ZZ.n), rinv,FCL_Elliptic_ZZ.n);//-hr^-1
        uint256 u2=mulmod(s, rinv,FCL_Elliptic_ZZ.n);//sr^-1

        uint256 Qx;
        uint256 Qy;
        (Qx,Qy)=FCL_Elliptic_ZZ.ecZZ_mulmuladd(r,y, u1, u2);

        return address(uint160(uint256(keccak256(abi.encodePacked(Qx, Qy)))));
    }

/*
    //ecdsa signature for test purpose only (who would like to have a private key onchain anyway ?)
    //K is nonce, kpriv is private key
    function ecdsa_sign(bytes32 message, uint256 k , uint256 kpriv) internal view returns(uint256 r, uint256 s)
    {
        r=FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(0,0, k, 0) ;//Calculate the curve point k.G (abuse ecmulmul add with v=0)
        r=addmod(0,r, FCL_Elliptic_ZZ.n); 
        s=mulmod(FCL_Elliptic_ZZ.FCL_nModInv(k), addmod(uint256(message), mulmod(r, kpriv, FCL_Elliptic_ZZ.n),FCL_Elliptic_ZZ.n),FCL_Elliptic_ZZ.n);//s=k^-1.(h+r.kpriv)

        
        if(r==0||s==0){
            revert();
        }


    }

    //ecdsa key derivation
    //kpriv is private key return (x,y) coordinates of associated Pubkey
    function ecdsa_derivKpub(uint256 kpriv) internal view returns(uint256 x, uint256 y)
    {
        
        x=FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(0,0, kpriv, 0) ;//Calculate the curve point k.G (abuse ecmulmul add with v=0)
        y=FCL_Elliptic_ZZ.ec_Decompress(x, 1);
       
        if (FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(x, y, kpriv, FCL_Elliptic_ZZ.n - 1) != 0) //extract correct y value
        {
            y=FCL_Elliptic_ZZ.p-y;
        }        

    }*/
 
}
